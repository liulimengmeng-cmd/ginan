#pragma once

#include "common/zhangProductDatumTransition.hpp"
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>

// Experiment-only matching of a PREVALIDATED common transition. This module
// has no KF/reference/product/QC access and cannot change the common plan.
// Runtime integration must supply actual endpoint/commit/state evidence and
// atomically commit the returned receipt with the existing release operation.
namespace zhangSingleEvent
{
using Key = std::pair<std::string, std::string>;

inline void require(bool value, const char* reason)
{ if (!value) throw std::runtime_error(reason); }

inline bool sha256(const std::string& value)
{
    return value.size() == 64 && std::all_of(value.begin(), value.end(),
        [](char c) { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'); });
}

inline Key key(const ZhangProductDatumTransition& value)
{ return {value.target.satellite, value.target.observable}; }

inline ZhangProductDatumTransition certificate(const std::string& wire)
{
    ZhangProductDatumTransition out;
    require(zhangReadDatumWireProof(wire, out) &&
        zhangDatumTransitionPureIntegerValid(out) &&
        zhangDatumWireProof(out) == wire, "SINGLE_EVENT_NONCANONICAL_CERTIFICATE");
    return out;
}

struct Plan
{
    std::string receiver, manifestSha, experimentPlanSha, streamSha;
    std::string parentCheckpointSha, availableCommitSha;
    long long eventEpoch = 0, availableEpoch = 0, applyEpoch = 0;
    // Complete transaction: all nonzero directions, both frozen frequencies.
    // Full wire bytes bind source/target, cycle payload and dual-frequency proof.
    std::vector<std::string> certificates;
};

inline std::map<Key, ZhangProductDatumTransition> validate(const Plan& plan)
{
    require(!plan.receiver.empty() && sha256(plan.manifestSha) && sha256(plan.experimentPlanSha) &&
        sha256(plan.streamSha) && sha256(plan.parentCheckpointSha) && sha256(plan.availableCommitSha),
        "SINGLE_EVENT_PLAN_BINDING_INVALID");
    // V1 admits normal first consumption only, not delayed historical events.
    require(plan.eventEpoch > 0 && plan.eventEpoch % 30 == 0 &&
        plan.availableEpoch == plan.eventEpoch && plan.applyEpoch == plan.eventEpoch,
        "SINGLE_EVENT_PLAN_TIME_OUTSIDE_V1");
    require(!plan.certificates.empty() && std::is_sorted(plan.certificates.begin(), plan.certificates.end()),
        "SINGLE_EVENT_PLAN_CERTIFICATE_ORDER");
    std::map<Key, ZhangProductDatumTransition> targets;
    std::map<std::string, std::vector<ZhangProductDatumTransition>> pairs;
    std::set<char> systems;
    for (const auto& wire : plan.certificates)
    {
        const auto value = certificate(wire);
        require(value.correctionShiftCycles != 0 && value.target.epoch == plan.eventEpoch &&
            value.source.epoch == plan.eventEpoch - 30,
            "SINGLE_EVENT_PLAN_NOT_DIRECT_NONZERO_EVENT");
        require(targets.emplace(key(value), value).second, "SINGLE_EVENT_DUPLICATE_PLAN_ENDPOINT");
        pairs[value.target.satellite].push_back(value);
        systems.insert(value.target.satellite.front());
    }
    require(systems.size() == 1 && std::string("GEC").find(*systems.begin()) != std::string::npos,
        "SINGLE_EVENT_PLAN_REQUIRES_ONE_GEC_SYSTEM");
    const std::map<char, std::set<std::string>> signals{
        {'G', {"L1C", "L2W"}}, {'E', {"L1C", "L5Q"}}, {'C', {"L1P", "L5P"}}};
    for (const auto& [satellite, pair] : pairs)
    {
        require(pair.size() == 2 &&
            std::set<std::string>{pair[0].target.observable, pair[1].target.observable} == signals.at(satellite.front()) &&
            zhangDatumDualPairValid(pair[0], pair[1]), "SINGLE_EVENT_INCOMPLETE_DUAL_TRANSACTION");
    }
    return targets;
}

struct HeldEndpoint
{
    bool existedBefore = false, existsAfterCommon = false;
    bool authorizedBefore = false, authorizedAfterCommon = false;
    std::string sourceEndpoint, targetEndpoint;
    std::string referenceBefore, referenceAfterCommon;
    // Stable scientific ownership/lineage, not an integer index into the KF.
    std::string stateOwner, physicalGenerationBefore, physicalGenerationAfter;
    std::string statisticalGenerationBefore, statisticalGenerationAfter;
    std::vector<std::string> consumedCertificates;
};

inline constexpr const char* stage =
    "AFTER_ALL_COMMON_REFERENCE_AND_LEGAL_TRANSPORT_BEFORE_CURRENT_OBSERVATION_UPDATE";

struct CommonEvidence
{
    std::string receiver, streamSha, availableCommitSha, parentCheckpointSha;
    std::string processingStage, commonPlanSha, commonStateSha;
    long long epoch = 0;
    std::map<Key, HeldEndpoint> held;
    // Classified by the ORIGINAL common transition logic, before consulting
    // this selector. Only noncommon transactions in the selected system.
    std::vector<std::string> noncommonCertificates;
};

enum class Status { PENDING, ELIGIBLE, INELIGIBLE, MISSED, ALREADY_TERMINAL };
enum class Action { KEEP, RESET_SELECTED_EVENT };

struct Decision
{
    Status status = Status::PENDING;
    std::string reason, manifestSha, commonPlanSha, commonStateSha;
    long long epoch = 0;
    std::vector<Key> targets;
    std::vector<std::string> stateOwners;
};

struct Receipt
{
    std::string manifestSha, experimentPlanSha, parentCheckpointSha;
    std::string commonPlanSha, commonStateSha, committedStateSha;
    long long epoch = 0;
    Status status = Status::PENDING;
    Action action = Action::KEEP;
    std::vector<Key> targets;
    std::vector<std::string> stateOwners;
    std::size_t releasedRank = 0;
    std::string reason;
};

inline Decision evaluate(const Plan& plan, const CommonEvidence& evidence, const Receipt* prior = nullptr)
{
    const auto targets = validate(plan);
    Decision result;
    result.manifestSha = plan.manifestSha;
    result.epoch = evidence.epoch;
    if (prior)
    {
        require(prior->manifestSha == plan.manifestSha && prior->experimentPlanSha == plan.experimentPlanSha &&
            prior->parentCheckpointSha == plan.parentCheckpointSha,
            "SINGLE_EVENT_RECEIPT_BINDING_CONFLICT");
        require(prior->status == Status::ELIGIBLE || prior->status == Status::INELIGIBLE || prior->status == Status::MISSED,
            "SINGLE_EVENT_RECEIPT_NOT_TERMINAL");
        require(prior->epoch <= evidence.epoch && prior->epoch >= plan.applyEpoch,
            "SINGLE_EVENT_RECEIPT_TIME_CONFLICT");
        require(sha256(prior->committedStateSha), "SINGLE_EVENT_RECEIPT_STATE_DIGEST_INVALID");
        if (prior->status == Status::ELIGIBLE)
        {
            std::vector<Key> expected;
            for (const auto& [target, ignored] : targets) expected.push_back(target);
            require(prior->epoch == plan.applyEpoch && prior->targets == expected &&
                prior->stateOwners.size() == expected.size() && sha256(prior->commonPlanSha) &&
                sha256(prior->commonStateSha), "SINGLE_EVENT_RECEIPT_TARGETS_INVALID");
            require(prior->action == Action::KEEP ? prior->releasedRank == 0 :
                prior->action == Action::RESET_SELECTED_EVENT && prior->releasedRank == expected.size(),
                "SINGLE_EVENT_RECEIPT_RELEASE_RANK_INVALID");
        }
        else require(prior->releasedRank == 0 && prior->targets.empty() && prior->stateOwners.empty() &&
            (prior->status == Status::MISSED ? prior->epoch > plan.applyEpoch : prior->epoch == plan.applyEpoch),
            "SINGLE_EVENT_INELIGIBLE_RECEIPT_INVALID");
        result.status = Status::ALREADY_TERMINAL;
        result.reason = "ALREADY_PROCESSED_NO_SECOND_INTERVENTION";
        return result;
    }
    if (evidence.epoch < plan.applyEpoch)
    { result.reason = "WAITING_PREDECLARED_APPLICATION"; return result; }
    if (evidence.epoch > plan.applyEpoch)
    { result.status = Status::MISSED; result.reason = "MISSED_EXPECTED_APPLICATION"; return result; }
    // Identity/config errors invalidate implementation acceptance. They are
    // not scientific ineligibility that can be skipped to pick another event.
    require(evidence.receiver == plan.receiver && evidence.streamSha == plan.streamSha &&
        evidence.parentCheckpointSha == plan.parentCheckpointSha &&
        evidence.availableCommitSha == plan.availableCommitSha && evidence.processingStage == stage &&
        sha256(evidence.commonPlanSha) && sha256(evidence.commonStateSha),
        "SINGLE_EVENT_COMMON_EVIDENCE_BINDING_CONFLICT");
    result.commonPlanSha = evidence.commonPlanSha;
    result.commonStateSha = evidence.commonStateSha;
    auto reject = [&](const char* reason) {
        result.status = Status::INELIGIBLE; result.reason = reason;
        result.targets.clear(); result.stateOwners.clear(); return result;
    };
    std::map<std::tuple<std::string, std::string, double>, std::string> unique;
    std::set<std::string> actual;
    for (const auto& wire : evidence.noncommonCertificates)
    {
        const auto value = certificate(wire);
        require(value.target.satellite.front() == targets.begin()->first.first.front(),
            "SINGLE_EVENT_COMMON_EVIDENCE_SYSTEM_CONFLICT");
        const auto identity = std::make_tuple(value.target.satellite, value.target.observable, value.target.epoch);
        const auto [where, inserted] = unique.emplace(identity, wire);
        require(inserted || where->second == wire, "SINGLE_EVENT_CONFLICTING_CERTIFICATE_DUPLICATE");
        actual.insert(wire); // Exact duplicate delivery is idempotent, not a second event.
    }
    const std::set<std::string> declared(plan.certificates.begin(), plan.certificates.end());
    if (actual != declared) return reject("MIXED_OR_NONMATCHING_NONCOMMON_TRANSACTION");
    for (const auto& [target, value] : targets)
    {
        const auto found = evidence.held.find(target);
        if (found == evidence.held.end() || !found->second.existedBefore || !found->second.existsAfterCommon)
            return reject("NO_STATE_NO_INTERVENTION");
        const auto& held = found->second;
        if (held.referenceBefore == target.first || held.referenceAfterCommon == target.first)
            return reject("REFERENCE_EVENT_OUTSIDE_V1");
        if (!held.authorizedBefore || !held.authorizedAfterCommon)
            return reject("TARGET_INTEGER_AUTHORIZATION_NOT_RETAINED");
        if (held.physicalGenerationBefore.empty() || held.statisticalGenerationBefore.empty() ||
            held.physicalGenerationBefore != held.physicalGenerationAfter ||
            held.statisticalGenerationBefore != held.statisticalGenerationAfter)
            return reject("TARGET_ORDINARY_LIFECYCLE_CHANGED");
        require(!held.stateOwner.empty(), "SINGLE_EVENT_STATE_OWNER_MISSING");
        if (held.sourceEndpoint != zhangDatumEndpointIdentity(value.source) ||
            held.targetEndpoint != zhangDatumEndpointIdentity(value.target))
            return reject("TARGET_APPLIED_ENDPOINT_DOES_NOT_MATCH");
        const auto wire = zhangDatumWireProof(value);
        if (std::count(held.consumedCertificates.begin(), held.consumedCertificates.end(), wire) != 1)
            return reject("TARGET_NOT_CONSUMED_EXACTLY_ONCE");
        result.targets.push_back(target); result.stateOwners.push_back(held.stateOwner);
    }
    require(std::set<std::string>(result.stateOwners.begin(), result.stateOwners.end()).size() == result.targets.size(),
        "SINGLE_EVENT_NONINDEPENDENT_STATE_OWNERS");
    result.status = Status::ELIGIBLE;
    result.reason = "MATCHED_COMPLETE_DIRECT_DUAL_NONREFERENCE_TRANSACTION";
    return result;
}

inline Receipt receipt(const Plan& plan, const Decision& decision, Action action,
    const std::string& committedStateSha, std::size_t releasedRank)
{
    validate(plan);
    require(decision.manifestSha == plan.manifestSha && sha256(committedStateSha) &&
        (decision.status == Status::ELIGIBLE || decision.status == Status::INELIGIBLE || decision.status == Status::MISSED),
        "SINGLE_EVENT_INVALID_COMMIT_RECEIPT");
    const bool release = decision.status == Status::ELIGIBLE && action == Action::RESET_SELECTED_EVENT;
    require(action == Action::KEEP || action == Action::RESET_SELECTED_EVENT, "SINGLE_EVENT_ACTION_INVALID");
    require(release ? releasedRank == decision.targets.size() && releasedRank > 0 : releasedRank == 0,
        "SINGLE_EVENT_PARTIAL_OR_UNEXPECTED_RELEASE");
    Receipt out;
    out.manifestSha = plan.manifestSha; out.experimentPlanSha = plan.experimentPlanSha;
    out.parentCheckpointSha = plan.parentCheckpointSha;
    out.commonPlanSha = decision.commonPlanSha; out.commonStateSha = decision.commonStateSha;
    out.committedStateSha = committedStateSha; out.epoch = decision.epoch;
    out.status = decision.status; out.action = action;
    out.targets = decision.targets; out.stateOwners = decision.stateOwners;
    out.releasedRank = releasedRank; out.reason = decision.reason;
    return out;
}
} // namespace zhangSingleEvent
