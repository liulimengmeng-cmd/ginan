#pragma once

#include <cmath>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <istream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Identity belongs to the accepted decision, NEVER to a transported row hash.
// Bounds are conditional on all parents being correct; union accounting visits
// every original decision once. This is a model-bound contract, not empirical
// calibration of the GNSS receiver or of an adaptive PAR procedure.
struct ZhangIntegerDecisionProof;
using ZhangDecisionProofPtr = std::shared_ptr<const ZhangIntegerDecisionProof>;
using ZhangDecisionProofs = std::vector<ZhangDecisionProofPtr>;
struct ZhangIntegerDecisionProof
{
    std::string id;
    std::string originalStatement;
    std::string observationProvenance;
    double conditionalFailureBound = 1;
    ZhangDecisionProofs parents;
};

/** Lossless publication dictionary, scoped by its owner to one epoch/commit.
 * It stores the writer's exact payload bytes, never a truncated summary or a
 * hash in place of the proof. Admission/risk accounting remains independent
 * of this optional wire representation. */
class ZhangDecisionProofDictionary
{
public:
    const std::map<std::string, std::string>& records() const { return records_; }

    bool intern(const std::string& id, const std::string& content,
        std::string& reason)
    {
        if (id.empty() || content.empty())
        { reason = "EMPTY_INTERNED_PROOF"; return false; }
        const auto prior = records_.find(id);
        if (prior != records_.end())
        {
            if (prior->second != content)
            { reason = "INTERNED_PROOF_ID_CONTENT_CONFLICT"; return false; }
            reason = "INTERNED_PROOF_REUSED";
            return true;
        }
        records_.emplace(id, content);
        reason = "INTERNED_PROOF_ADDED";
        return true;
    }

    const std::string* resolve(const std::string& id) const
    {
        const auto found = records_.find(id);
        return found == records_.end() ? nullptr : &found->second;
    }

    bool resolve(const std::vector<std::string>& ids,
        std::vector<std::string>& contents, std::string& reason) const
    {
        std::vector<std::string> resolved;
        resolved.reserve(ids.size());
        for (const auto& id : ids)
        {
            const auto content = resolve(id);
            if (!content)
            { reason = "INTERNED_PROOF_REFERENCE_MISSING"; return false; }
            resolved.push_back(*content);
        }
        contents = std::move(resolved);
        reason = "INTERNED_PROOF_REFERENCES_RESOLVED";
        return true;
    }

    void write(std::ostream& stream) const
    {
        stream << "INTERNED_V1\n" << records_.size() << '\n';
        for (const auto& [id, content] : records_)
        {
            writeField(stream, id);
            writeField(stream, content);
        }
    }

    // Explicit size limits prevent malformed length fields allocating an
    // unbounded amount of memory. A failed read leaves the prior dictionary
    // unchanged, including when a repeated ID has conflicting content.
    bool read(std::istream& stream, std::string& reason,
        std::size_t maximumRecords = 1000000,
        std::size_t maximumBytes = 1024ULL * 1024ULL * 1024ULL)
    {
        std::string version;
        if (!std::getline(stream, version) || version != "INTERNED_V1")
        { reason = "INTERNED_PROOF_FORMAT_UNSUPPORTED"; return false; }
        std::size_t count = 0;
        if (!readLength(stream, count) || count > maximumRecords)
        { reason = "INTERNED_PROOF_RECORD_LIMIT"; return false; }
        ZhangDecisionProofDictionary proposed;
        std::size_t remaining = maximumBytes;
        for (std::size_t i = 0; i < count; ++i)
        {
            std::string id, content;
            if (!readField(stream, id, remaining) ||
                !readField(stream, content, remaining))
            { reason = "INTERNED_PROOF_TRUNCATED_OR_OVERSIZE"; return false; }
            if (!proposed.intern(id, content, reason)) return false;
        }
        records_ = std::move(proposed.records_);
        reason = "INTERNED_PROOF_DICTIONARY_READ";
        return true;
    }

private:
    static void writeField(std::ostream& stream, const std::string& value)
    {
        stream << value.size() << '\n';
        stream.write(value.data(), static_cast<std::streamsize>(value.size()));
        stream.put('\n');
    }
    static bool readLength(std::istream& stream, std::size_t& size)
    {
        std::string line;
        if (!std::getline(stream, line) || line.empty()) return false;
        std::size_t value = 0;
        for (const unsigned char digit : line)
        {
            if (digit < '0' || digit > '9' ||
                value > (std::numeric_limits<std::size_t>::max() - (digit - '0')) / 10)
                return false;
            value = value * 10 + (digit - '0');
        }
        size = value;
        return true;
    }
    static bool readField(std::istream& stream, std::string& value,
        std::size_t& remaining)
    {
        std::size_t length = 0;
        if (!readLength(stream, length) || length > remaining ||
            length > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()))
            return false;
        std::string field(length, '\0');
        if (!stream.read(field.data(), static_cast<std::streamsize>(length)) ||
            stream.get() != '\n') return false;
        remaining -= length;
        value = std::move(field);
        return true;
    }
    std::map<std::string, std::string> records_;
};

// Enabling a new writer format without an adapted consumer is deliberately
// impossible through this helper. The default capability argument is false.
inline bool zhangInternHandoffProofsEnabled(bool readerSupportsInternedV1 = false)
{
    const auto setting = std::getenv("ZHANG_R51_INTERN_HANDOFF_PROOFS");
    return readerSupportsInternedV1 && setting && std::string(setting) == "1";
}

/** Canonical atom content for dictionaries. Parent IDs retain their original
 * order and the floating bound is encoded by its complete IEEE bits. */
inline std::string zhangDecisionProofAtomContent(const ZhangIntegerDecisionProof& proof)
{
    std::ostringstream stream;
    auto field = [&](const std::string& text)
    {
        stream << text.size() << '\n';
        stream.write(text.data(), static_cast<std::streamsize>(text.size()));
        stream.put('\n');
    };
    field(proof.id);
    field(proof.originalStatement);
    field(proof.observationProvenance);
    stream << std::bit_cast<std::uint64_t>(proof.conditionalFailureBound) << '\n'
        << proof.parents.size() << '\n';
    for (const auto& parent : proof.parents) field(parent ? parent->id : std::string{});
    return stream.str();
}

/** Legacy detail payloads remain byte-identical inside INTERNED_V1. The ID
 * embedded in each original statement is the dictionary key, not a hash or a
 * satellite-specific alias. This prevents an ID being attached to a different
 * original decision when reconstructing legacy per-product details. */
inline bool zhangHandoffProofStatementId(const std::string& statement,
    std::string& id)
{
    std::istringstream input(statement);
    std::string prefix(3, '\0');
    if (!input.read(prefix.data(), 3) || prefix != "id=" || input.peek() != '"' ||
        !(input >> std::quoted(id)) || id.empty()) return false;
    return true;
}

struct ZhangHandoffProofReferenceRow
{
    std::vector<std::string> proofIds;
    std::vector<std::string> detailReferences;
};

struct ZhangHandoffProofCertificate
{
    bool interned = false;
    ZhangDecisionProofDictionary dictionary;
    std::vector<ZhangHandoffProofReferenceRow> rows;
};

/** Parse/validate the optional dictionary portion of a full bundle certificate.
 * Old inline proof_detail certificates remain accepted unchanged. Quoted
 * values may contain spaces, escaped quotes, or newlines; markers inside them
 * are data. A malformed new format is never silently treated as legacy. */
inline bool zhangReadHandoffProofCertificate(std::istream& input,
    ZhangHandoffProofCertificate& certificate, std::string& reason)
{
    ZhangHandoffProofCertificate proposed;
    bool formatSeen = false, dictionarySeen = false, newFieldsSeen = false;
    bool handoffRowsBeforeFormat = false;
    bool malformedLogicalLine = false;
    auto logicalLine = [&](std::string& line)
    {
        line.clear();
        bool quoted = false, escaped = false;
        char c;
        while (input.get(c))
        {
            if (c == '\n' && !quoted) return true;
            line.push_back(c);
            if (quoted && escaped) escaped = false;
            else if (quoted && c == '\\') escaped = true;
            else if (c == '"') quoted = !quoted;
        }
        malformedLogicalLine = quoted || escaped;
        return !line.empty() && !malformedLogicalLine;
    };
    std::string line;
    while (logicalLine(line))
    {
        if (line.rfind("handoff_proof_format=", 0) == 0)
        {
            if (formatSeen || handoffRowsBeforeFormat || line != "handoff_proof_format=INTERNED_V1")
            { reason = "HANDOFF_PROOF_FORMAT_UNSUPPORTED_OR_DUPLICATE"; return false; }
            formatSeen = true;
            proposed.interned = true;
            continue;
        }
        if (line == "handoff_proof_dictionary_begin")
        {
            if (!formatSeen || dictionarySeen)
            { reason = "HANDOFF_PROOF_DICTIONARY_UNEXPECTED"; return false; }
            if (!proposed.dictionary.read(input, reason)) return false;
            if (!std::getline(input, line) || line != "handoff_proof_dictionary_end")
            { reason = "HANDOFF_PROOF_DICTIONARY_END_MISSING"; return false; }
            dictionarySeen = true;
            continue;
        }
        if (line.rfind("handoff_satellite=", 0) != 0) continue;
        if (!formatSeen) handoffRowsBeforeFormat = true;
        // Read field values with the same quoting rules as the legacy writer.
        // Unrelated metadata fields are retained by the certificate itself.
        std::istringstream fields(line);
        ZhangHandoffProofReferenceRow references;
        std::size_t declaredCount = 0;
        bool countSeen = false, inlineDetails = false;
        while (fields >> std::ws && fields.peek() != std::char_traits<char>::eof())
        {
            std::string key, value;
            char c;
            while (fields.get(c) && c != '=' && c != ' ' && c != '\t') key.push_back(c);
            if (!fields || c != '=')
            { reason = "HANDOFF_PROOF_FIELD_INVALID"; return false; }
            if (fields.peek() == '"')
            {
                if (!(fields >> std::quoted(value)))
                { reason = "HANDOFF_PROOF_QUOTED_VALUE_INVALID"; return false; }
            }
            else if (!(fields >> value))
            { reason = "HANDOFF_PROOF_FIELD_VALUE_MISSING"; return false; }
            if (key == "proof_count")
            {
                if (countSeen || value.empty())
                { reason = "HANDOFF_PROOF_COUNT_INVALID"; return false; }
                for (const unsigned char digit : value)
                {
                    if (digit < '0' || digit > '9' || declaredCount >
                        (std::numeric_limits<std::size_t>::max() - (digit - '0')) / 10)
                    { reason = "HANDOFF_PROOF_COUNT_INVALID"; return false; }
                    declaredCount = declaredCount * 10 + (digit - '0');
                }
                countSeen = true;
            }
            else if (key == "proof") references.proofIds.push_back(std::move(value));
            else if (key == "proof_detail_ref")
            {
                newFieldsSeen = true;
                references.detailReferences.push_back(std::move(value));
            }
            else if (key == "proof_detail") inlineDetails = true;
        }
        if (formatSeen)
        {
            if (inlineDetails || (!countSeen && (!references.proofIds.empty() ||
                    !references.detailReferences.empty())) ||
                declaredCount != references.proofIds.size() ||
                declaredCount != references.detailReferences.size())
            { reason = "HANDOFF_PROOF_REFERENCE_COUNT_MISMATCH"; return false; }
            std::multiset<std::string> ids(references.proofIds.begin(), references.proofIds.end());
            const std::multiset<std::string> refs(references.detailReferences.begin(),
                references.detailReferences.end());
            if (ids != refs)
            { reason = "HANDOFF_PROOF_REFERENCE_ID_MISMATCH"; return false; }
            proposed.rows.push_back(std::move(references));
        }
    }
    if (input.bad() || (formatSeen && (!input.eof() || malformedLogicalLine)))
    { reason = "HANDOFF_PROOF_CERTIFICATE_READ_FAILED"; return false; }
    if (newFieldsSeen && !formatSeen)
    { reason = "HANDOFF_PROOF_FORMAT_MISSING"; return false; }
    if (formatSeen && !dictionarySeen)
    { reason = "HANDOFF_PROOF_DICTIONARY_MISSING"; return false; }
    if (formatSeen)
    {
        for (const auto& [id, content] : proposed.dictionary.records())
        {
            std::string embeddedId;
            if (!zhangHandoffProofStatementId(content, embeddedId) || embeddedId != id)
            { reason = "HANDOFF_PROOF_DICTIONARY_PAYLOAD_ID_MISMATCH"; return false; }
        }
        for (const auto& references : proposed.rows)
            for (const auto& id : references.detailReferences)
                if (!proposed.dictionary.resolve(id))
                { reason = "INTERNED_PROOF_REFERENCE_MISSING"; return false; }
    }
    certificate = std::move(proposed);
    reason = certificate.interned ? "HANDOFF_INTERNED_PROOFS_VERIFIED" : "HANDOFF_LEGACY_INLINE_PROOFS";
    return true;
}

struct ZhangDecisionRiskClosure
{
    bool valid = false;
    double bound = 0;
    std::string reason = "MISSING_DECISION_PROOF";
    std::map<std::string, ZhangDecisionProofPtr> atoms;
};

struct ZhangIncrementalDecisionRisk
{
    bool valid = false;
    int baselineAtoms = 0;
    int unionAtoms = 0;
    int addedAtoms = 0;
    double bound = 0;
    std::string reason = "NOT_EVALUATED";
};

inline ZhangDecisionRiskClosure zhangDecisionRiskClosure(const ZhangDecisionProofs& roots)
{
    ZhangDecisionRiskClosure result;
    std::set<std::string> visiting;
    auto visit = [&](auto&& self, const ZhangDecisionProofPtr& proof) -> bool
    {
        if (!proof || proof->id.empty() || proof->originalStatement.empty() ||
            proof->observationProvenance.empty() ||
            !std::isfinite(proof->conditionalFailureBound) ||
            proof->conditionalFailureBound < 0 || proof->conditionalFailureBound > 1)
            return false;
        if (visiting.count(proof->id))
        { result.reason = "DECISION_DEPENDENCY_CYCLE"; return false; }
        auto prior = result.atoms.find(proof->id);
        if (prior != result.atoms.end())
        {
            // Two different immutable objects claiming an ID are not silently
            // merged, even if their current numerical integer values agree.
            if (prior->second.get() != proof.get())
            { result.reason = "DECISION_ID_COLLISION"; return false; }
            return true;
        }
        visiting.insert(proof->id);
        for (const auto& parent : proof->parents)
            if (!self(self, parent)) return false;
        visiting.erase(proof->id);
        result.atoms.emplace(proof->id, proof);
        result.bound += proof->conditionalFailureBound;
        return true;
    };
    for (const auto& root : roots) if (!visit(visit, root)) return result;
    result.valid = true;
    result.reason = "DECISION_DEPENDENCIES_CLOSED";
    return result;
}

inline ZhangDecisionProofs zhangMergeDecisionProofs(
    const ZhangDecisionProofs& first, const ZhangDecisionProofs& second)
{
    auto result = first;
    for (const auto& proof : second)
    {
        bool present = false;
        for (const auto& old : result) present |= old.get() == proof.get();
        if (!present) result.push_back(proof);
    }
    return result;
}

/** Risk contributed by candidate decisions that are not already ancestors of
 * the baseline certificate.  Comparing two accumulated floating-point sums is
 * unsafe here: a shared proof can be rounded differently after a view or HNF
 * change.  Immutable proof IDs define the statistical decision exactly once;
 * an ID collision remains fail-closed through zhangDecisionRiskClosure. */
inline ZhangIncrementalDecisionRisk zhangIncrementalDecisionRisk(
    const ZhangDecisionProofs& baseline,
    const ZhangDecisionProofs& candidate)
{
    ZhangIncrementalDecisionRisk result;
    const auto baselineClosure = zhangDecisionRiskClosure(baseline);
    if (!baselineClosure.valid)
    {
        result.reason = baselineClosure.reason;
        return result;
    }
    const auto unionClosure = zhangDecisionRiskClosure(
        zhangMergeDecisionProofs(baseline, candidate));
    if (!unionClosure.valid)
    {
        result.reason = unionClosure.reason;
        return result;
    }
    result.baselineAtoms = static_cast<int>(baselineClosure.atoms.size());
    result.unionAtoms = static_cast<int>(unionClosure.atoms.size());
    for (const auto& [id, proof] : unionClosure.atoms)
    {
        if (baselineClosure.atoms.contains(id)) continue;
        result.bound = std::min(1.0,
            result.bound + proof->conditionalFailureBound);
        result.addedAtoms++;
    }
    result.valid = true;
    result.reason = "NONE";
    return result;
}
