#include "common/zhangIntegerDecisionProof.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool condition, const char* reason)
{
    if (!condition) throw std::runtime_error(reason);
}

int main()
{
    ZhangDecisionProofDictionary dictionary;
    std::string reason;
    const std::string payload("row=2*x\nKEY=a\0b\n", 17);
    require(dictionary.intern("decision-1", payload, reason), "first intern failed");
    require(dictionary.intern("decision-1", payload, reason) && dictionary.records().size() == 1,
        "identical content was duplicated");
    require(!dictionary.intern("decision-1", payload + "different", reason) &&
        reason == "INTERNED_PROOF_ID_CONTENT_CONFLICT", "same ID different content accepted");
    require(!dictionary.intern("", payload, reason), "empty ID accepted");
    require(dictionary.intern("decision-2", "multiline\nproof detail\n", reason), "second intern failed");
    std::stringstream wire;
    dictionary.write(wire);
    require(wire.str().rfind("INTERNED_V1\n2\n", 0) == 0, "wrong dictionary version/count");
    ZhangDecisionProofDictionary restored;
    require(restored.read(wire, reason) && restored.records() == dictionary.records(),
        "lossless dictionary round trip failed");
    std::vector<std::string> resolved;
    require(restored.resolve({"decision-2", "decision-1", "decision-1"}, resolved, reason) &&
        resolved.size() == 3 && resolved[1] == payload && resolved[2] == payload,
        "proof references lost order or duplicate uses");
    const auto before = resolved;
    require(!restored.resolve({"missing"}, resolved, reason) && resolved == before,
        "unknown proof ID resolved or partially replaced output");
    const auto recordsBefore = restored.records();
    std::stringstream truncated("INTERNED_V1\n1\n3\nabc\n9999999999999999999999999999999999999\n");
    require(!restored.read(truncated, reason) && restored.records() == recordsBefore,
        "invalid length altered dictionary");
    std::stringstream oversized;
    dictionary.write(oversized);
    require(!restored.read(oversized, reason, 100, 3) && restored.records() == recordsBefore,
        "oversize proof bypassed reader limits");
    std::stringstream conflict("INTERNED_V1\n2\n1\na\n1\nx\n1\na\n1\ny\n");
    require(!restored.read(conflict, reason) &&
        reason == "INTERNED_PROOF_ID_CONTENT_CONFLICT" && restored.records() == recordsBefore,
        "wire ID-content conflict accepted");
    std::stringstream oldFormat("proof_detail=legacy\n");
    require(!restored.read(oldFormat, reason), "legacy text guessed as interned dictionary");
    require(!zhangInternHandoffProofsEnabled(), "new wire format enabled without adapted reader");

    auto parent = std::make_shared<ZhangIntegerDecisionProof>();
    parent->id = "parent";
    parent->originalStatement = "a=1";
    parent->observationProvenance = "epoch";
    parent->conditionalFailureBound = 0.0001;
    auto child = std::make_shared<ZhangIntegerDecisionProof>();
    child->id = "child";
    child->originalStatement = "b=2\ncomplete witness";
    child->observationProvenance = "epoch";
    child->conditionalFailureBound = 0.0002;
    child->parents = {parent};
    const auto content = zhangDecisionProofAtomContent(*child);
    require(dictionary.intern(child->id, content, reason), "canonical atom intern failed");
    auto sameIdDifferentContent = std::make_shared<ZhangIntegerDecisionProof>(*child);
    sameIdDifferentContent->originalStatement += "changed";
    require(!dictionary.intern(child->id, zhangDecisionProofAtomContent(*sameIdDifferentContent), reason),
        "canonical atom content conflict missed");
    // Text interning does not relax the decision authority identity contract.
    require(!zhangDecisionRiskClosure({child, sameIdDifferentContent}).valid,
        "publication dictionary relaxed risk ID collision check");
    require(zhangDecisionRiskClosure({child, child}).valid,
        "same immutable proof unexpectedly rejected");

    // Reconstruct per-product legacy detail vectors from the real new wire
    // references; quoting/newlines are not searched as substring markers.
    const std::string detail1 = "id=\"decision-1\" source=\"epoch\" statement=\"a=1\"";
    const std::string detail2 = "id=\"decision-2\" source=\"epoch\" statement=\"b=2\n"
        "handoff_proof_dictionary_begin\ncomplete witness\"";
    ZhangDecisionProofDictionary handoff;
    require(handoff.intern("decision-1", detail1, reason) &&
        handoff.intern("decision-2", detail2, reason), "handoff intern failed");
    auto handoffWire = [&](const std::string& rows)
    {
        std::ostringstream stream;
        stream << "epoch=2019-01-01 02:00:00\n"
            << "handoff_proof_format=INTERNED_V1\n" << rows
            << "handoff_proof_dictionary_begin\n";
        handoff.write(stream);
        stream << "handoff_proof_dictionary_end\n";
        return stream.str();
    };
    const std::string row1 = "handoff_satellite=G02 observable=L1 proof_count=2 "
        "proof=\"decision-1\" proof=\"decision-2\" proof_detail_ref=\"decision-2\" "
        "proof_detail_ref=\"decision-1\"\n";
    const std::string row2 = "handoff_satellite=G03 observable=L1 proof_count=1 "
        "proof=\"decision-1\" proof_detail_ref=\"decision-1\"\n";
    std::istringstream certificateWire(handoffWire(row1 + row2));
    ZhangHandoffProofCertificate certificate;
    require(zhangReadHandoffProofCertificate(certificateWire, certificate, reason) &&
        certificate.interned && certificate.rows.size() == 2,
        "handoff full certificate reference round trip failed");
    std::vector<std::string> reconstructed;
    require(certificate.dictionary.resolve(certificate.rows[0].detailReferences,
        reconstructed, reason) && reconstructed == std::vector<std::string>{detail2, detail1},
        "legacy detail bytes/order were not reconstructed");
    auto rejectCertificate = [&](const std::string& text, const char* why)
    {
        std::istringstream input(text);
        const auto originalRecords = certificate.dictionary.records();
        require(!zhangReadHandoffProofCertificate(input, certificate, reason) &&
            certificate.dictionary.records() == originalRecords, why);
    };
    rejectCertificate(handoffWire("handoff_satellite=G02 proof_count=1 proof=\"unknown\" "
        "proof_detail_ref=\"unknown\"\n"), "unknown dictionary reference accepted");
    rejectCertificate(handoffWire("handoff_satellite=G02 proof_count=1 proof=\"decision-1\" "
        "proof_detail_ref=\"decision-2\"\n"), "proof/ref identity mismatch accepted");
    rejectCertificate(handoffWire("handoff_satellite=G02 proof_count=2 proof=\"decision-1\" "
        "proof_detail_ref=\"decision-1\"\n"), "proof reference count mismatch accepted");
    rejectCertificate("handoff_proof_format=INTERNED_V1\n" + row1,
        "missing dictionary accepted as legacy");
    rejectCertificate(handoffWire(row1) + "handoff_satellite=G02 proof=\"unterminated",
        "truncated quoted new-format field accepted");
    rejectCertificate("handoff_proof_format=INTERNED_V2\n", "unsupported new format accepted");
    std::ostringstream legacy;
    legacy << "handoff_satellite=G02 proof_count=1 proof=\"decision-2\" proof_detail="
        << std::quoted(detail2) << '\n';
    std::istringstream legacyWire(legacy.str());
    require(zhangReadHandoffProofCertificate(legacyWire, certificate, reason) &&
        !certificate.interned, "legacy inline details no longer readable");
    std::cout << "PASS proof_dictionary lossless_roundtrip_and_conflicts\n";
}
