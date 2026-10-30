#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <utility>

/** A target is a request on one COMPLETE current affine conditioner domain.
 * Version is local to the immutable posterior/chart/policy scope. The owner
 * advances it only after an exact, feasible, rank-growing affine union commits.
 * A fingerprint, when supplied, must describe all current H AND b, not merely
 * the last accepted row, the target residual or a matrix-rank diagnostic.
 *
 * The scheduler never decides that two integer images are equivalent: its
 * canonicalTarget must come from the existing exact target canonicalizer on
 * THIS domain. In particular, it must not divide out integer factors/gcds.
 */
struct ZhangR49BridgeRequestIdentity
{
    std::string posteriorId, chartId, policyId;
    std::uint64_t domainVersion = 0;
    std::string canonicalTarget;
    std::string conditionerFingerprint;

    bool valid() const
    {
        return !posteriorId.empty() && !chartId.empty() && !policyId.empty() &&
            !canonicalTarget.empty();
    }

    std::string key() const
    {
        if (!valid()) return {};
        // Length framing prevents names containing '|', ':' or brackets from
        // colliding. Version zero is a real initial domain, not an unset value.
        std::string result = "R49_VERSIONED_TARGET_V1|";
        const auto append = [&](const std::string& value)
        { result += std::to_string(value.size()) + ":" + value; };
        append(posteriorId);
        append(chartId);
        append(policyId);
        append(std::to_string(domainVersion));
        append(conditionerFingerprint);
        append(canonicalTarget);
        return result;
    }
};

template<class Domain>
ZhangR49BridgeRequestIdentity zhangR49VersionedBridgeRequest(
    const Domain& domain, std::string canonicalTarget,
    std::string completeConditionerFingerprint = {})
{
    return {domain.posteriorId, domain.chartId, domain.policyId, domain.version,
        std::move(canonicalTarget), std::move(completeConditionerFingerprint)};
}

struct ZhangR49BridgeScheduler {
    // These are monotone search-accounting fields, deliberately NOT part of
    // the route trial below. Rejecting/rolling back a candidate cannot refund
    // a real ILS call or its consumed family allocation.
    int attempts=0;
    double spentRisk=0;
    std::set<std::string> testedRequests;

    bool claimRequest(const ZhangR49BridgeRequestIdentity& request)
    {
        const auto identity = request.key();
        return !identity.empty() && testedRequests.insert(identity).second;
    }

    // Claiming and consuming a call are separate: exact target compilation
    // may reject a request before the ILS callback is reached. Charge a real
    // call immediately BEFORE invoking that callback, including failed ILS.
    bool consumeSearchCall(int maximumAttempts)
    {
        if (maximumAttempts <= 0 || attempts < 0 || attempts >= maximumAttempts)
            return false;
        ++attempts;
        return true;
    }

    // Account only the search receipt's nonnegative spent allocation. The
    // caller still owns its existing overall family-risk/budget gate. Neither
    // a failed ratio nor a rejected trial may pass a negative refund here.
    bool chargeRisk(double spent)
    {
        if (!std::isfinite(spent) || spent < 0 || !std::isfinite(spentRisk) ||
            spentRisk < 0 || !std::isfinite(spentRisk + spent)) return false;
        spentRisk += spent;
        return true;
    }
};

template<class Domain, class Validator>
bool zhangR49CommitTrial(Domain& safe, Domain trial, Validator validate) {
    if(!validate(trial)) return false;
    safe=std::move(trial);
    return true;
}

/** Commit only an exact domain extension which gained independent rank.
 * old/new rank must be computed by the caller from the COMPLETE affine union
 * with the existing exact backend; accepted proposal row count is not rank.
 * validate still owns integer feasibility, RHS consistency, proof closure,
 * numerical validity and unchanged-posterior identity. A same-rank result may
 * be a duplicate/entailed row, so it must not create a new retry opportunity.
 */
template<class Route, class Validator>
bool zhangR49CommitGrowingTrial(Route& safe, Route trial,
    int oldExactRank, int newExactRank, Validator validate)
{
    if (oldExactRank < 0 || newExactRank <= oldExactRank ||
        safe.domain.posteriorId.empty() || safe.domain.chartId.empty() ||
        safe.domain.policyId.empty() ||
        safe.domain.posteriorId != trial.domain.posteriorId ||
        safe.domain.chartId != trial.domain.chartId ||
        safe.domain.policyId != trial.domain.policyId ||
        safe.domain.version != trial.domain.version ||
        safe.domain.version == std::numeric_limits<std::uint64_t>::max() ||
        !validate(trial)) return false;
    trial.domain.version = safe.domain.version + 1;
    safe = std::move(trial);
    return true;
}
