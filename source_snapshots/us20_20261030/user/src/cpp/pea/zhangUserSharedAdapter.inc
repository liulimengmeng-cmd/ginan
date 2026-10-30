// Included by ppp_ambres.cpp after the candidate policy helpers. The adapter
// owns no integer decisions and never changes FLOAT observations or moments.
#include "common/zhangUserArcGeneration.hpp"
#include "common/zhangUserArPolicy.hpp"
#include "common/zhangUserGonsisComponents.hpp"
#include "common/zhangUserProductChart.hpp"

struct ZhangUserSharedEndpoint
{
    SatSys satellite;
    E_ObsCode firstCode = E_ObsCode::NONE;
    E_ObsCode secondCode = E_ObsCode::NONE;
    string identity;
    string firstIntegerIdentity;
    string secondIntegerIdentity;
    string physicalIdentity;
    string stateIdentity;
    string productIdentity;
    // Non-empty only for GONSIS pair-scoped search. Equality authorises
    // within-component SD rows, not individual coordinates relative to the
    // KF chart reference. The default legacy adapter keeps its old authority.
    string componentDomain;
    VectorXd firstRow;
    VectorXd secondRow;
    bool firstReferenceZero = false;
    bool secondReferenceZero = false;
    bool eligible = false;
    std::size_t effectiveSamples = 0;
    double elevation = std::numeric_limits<double>::quiet_NaN();
};

static bool zhangUserSharedCodes(E_Sys system, E_ObsCode& first, E_ObsCode& second)
{
    if (acsConfig.zhangPppAr.user_adapter)
    {
        const auto configured = acsConfig.zhangPppAr.baseline_observables.find(system);
        if (configured == acsConfig.zhangPppAr.baseline_observables.end() ||
            configured->second.size() != 2) return false;
        first = configured->second[0];
        second = configured->second[1];
    }
    else
    {
        // Pick the first configured physical code on each of two frequencies,
        // not two alternatives on the same frequency. GPS fixture: L1C/L2W.
        const auto configured = acsConfig.code_priorities.find(system);
        if (configured == acsConfig.code_priorities.end()) return false;
        const auto frequencies = code2Freq.find(system);
        if (frequencies == code2Freq.end()) return false;
        E_FType firstFrequency = NONE;
        for (const auto code : configured->second)
        {
            const auto frequency = frequencies->second.find(code);
            if (frequency == frequencies->second.end()) continue;
            if (first == E_ObsCode::NONE)
            {
                first = code;
                firstFrequency = frequency->second;
            }
            else if (frequency->second != firstFrequency)
            {
                second = code;
                break;
            }
        }
    }
    if (first == E_ObsCode::NONE || second == E_ObsCode::NONE) return false;
    const auto frequencies = code2Freq.find(system);
    return frequencies != code2Freq.end() && frequencies->second.count(first) &&
        frequencies->second.count(second) &&
        frequencies->second.at(first) != frequencies->second.at(second) &&
        zhangE27Wavelength(system, first) > 0 && zhangE27Wavelength(system, second) > 0;
}

// Require a real, current, direct phase OSB. The legacy measurement lookup can
// interpolate/extrapolate and overwrites the returned tini with the query time;
// that returned value cannot identify an immutable continuity interval. Read
// the exact stored interval instead. This is an external-product eligibility
// check, NOT an internal server-lattice certificate or a covariance rewrite.
static string zhangUserSharedExternalPhaseIdentity(
    const SatSys& satellite, E_ObsCode physicalCode, GTime time)
{
    E_ObsCode lookupCode = physicalCode;
    if (acsConfig.ssrInOpts.one_freq_phase_bias)
    {
        const auto frequencies = code2Freq.find(satellite.sys);
        if (frequencies == code2Freq.end() || !frequencies->second.count(physicalCode)) return {};
        lookupCode = freq2CodeHax(satellite.sys, frequencies->second.at(physicalCode));
    }
    const string id = satellite.id() + ":" + string(1, satellite.sysChar());
    const auto device = biasMaps[PHAS].find(id);
    if (device == biasMaps[PHAS].end()) return {};
    const auto first = device->second.find(lookupCode);
    if (first == device->second.end()) return {};
    const auto zero = first->second.find(E_ObsCode::NONE);
    if (zero == first->second.end()) return {};
    const auto selected = zero->second.lower_bound(time);
    if (selected == zero->second.end()) return {};
    const auto& entry = selected->second;
    if (entry.source != "bsx" || entry.measType != PHAS || entry.Sat != satellite ||
        !entry.name.empty() || entry.cod1 != lookupCode || entry.cod2 != E_ObsCode::NONE ||
        entry.tini == GTime::noTime() || entry.tfin == GTime::noTime() ||
        time < entry.tini || !(time < entry.tfin) ||
        !std::isfinite(entry.bias) || !std::isfinite(entry.var) || entry.var < 0 ||
        !std::isfinite(entry.slop) || !std::isfinite(entry.slpv) || entry.slpv < 0)
        return {};
    std::ostringstream token;
    token << std::setprecision(17) << "EXTERNAL_PHASE_OSB:" << id << ":"
          << static_cast<int>(physicalCode) << ":" << static_cast<int>(lookupCode) << ":"
          << entry.tini.bigTime << ":" << entry.tfin.bigTime << ":"
          << entry.refTime.bigTime << ":" << entry.bias << ":" << entry.slop;
    return token.str();
}

static string zhangUserSharedInternalPhaseIdentity(
    const KFState& state, const string& receiver, const SatSys& satellite,
    E_ObsCode code, GTime time, bool pairScoped = false,
    ZhangUserGonsisComponentCertificate* component = nullptr)
{
    if(zhangUserDatumPolicyRequested())
    {
        ZhangInternalProduct current;
        if(!queryZhangInternalProduct(time,satellite,code,current) || !current.pppar_usable ||
            !current.dual_frequency_ar_valid || !zhangFormalPppArProductSolution(current.solution) ||
            !zhangUserChartHistoryAuthorized(state,receiver,satellite,static_cast<int>(code)))return {};
    }
    const auto retained=zhangUserChartAnchor(state,receiver,satellite,static_cast<int>(code));
    if (zhangUserGonsisComponentCertificateValid(retained))
    {
        if (!pairScoped && !zhangPppArUserAmbiguityIntegerValid(state,receiver,satellite,code)) return {};
        if (component) *component=retained;
        return zhangUserProductHistoryIdentity(state,receiver,satellite,static_cast<int>(code),retained);
    }
    if (zhangUserChartIsolated(state,receiver,satellite.sys,static_cast<int>(code)) ||
        state.metaDataMap.count(zhangUserChartKey(receiver,satellite,static_cast<int>(code))+":UNANCHORED_MAP")) return {};
    ZhangInternalProduct product;
    if (!queryZhangInternalProduct(time, satellite, code, product) ||
        !product.pppar_usable || !product.dual_frequency_ar_valid) return {};
    const ZhangUserGonsisComponentCertificate certificate{
        product.pppar_usable, product.dual_frequency_ar_valid,
        product.integer_component_id, product.integer_component_version,
        product.integer_alignment_generation};
    if (pairScoped)
    {
        if (!zhangFormalPppArProductSolution(product.solution) ||
            !zhangUserGonsisComponentCertificateValid(certificate)) return {};
        if (component) *component = certificate;
        return zhangUserProductHistoryIdentity(state,receiver,satellite,static_cast<int>(code),certificate);
    }
    else if (!zhangPppArUserAmbiguityIntegerValid(state, receiver, satellite, code)) return {};
    return "INTERNAL_CERTIFIED_PHASE:" + satellite.id() + ":" +
        std::to_string(static_cast<int>(code)) + ":" +
        zhangUserArEncodedField(product.phase_product_segment_id) + ":" +
        std::to_string(product.discontinuity_counter) + ":" +
        std::to_string(product.datum_version) + ":" +
        zhangUserArEncodedField(product.integer_component_id) + ":" +
        std::to_string(product.integer_component_version) + ":" +
        std::to_string(product.integer_alignment_generation) + ":PUBLIC:" +
        zhangUserArEncodedField(product.public_continuity_metadata
            ? product.public_phase_segment_id : "LEGACY") + ":IA:" +
        std::to_string(zhangUserArcGeneration(state.metaDataMap,
            zhangUserIntegerAuthorizationGenerationKey(receiver, satellite.id(),
                static_cast<int>(code))));
}

static bool zhangUserSharedEndpointRow(
    const KFState& state, Receiver* receiver, const SatSys& satellite,
    E_ObsCode code, VectorXd& row, bool& referenceZero)
{
    referenceZero = false;
    row = VectorXd::Zero(state.x.size());
    KFKey key;
    key.type = KF::AMBIGUITY;
    key.str = receiver->id;
    key.Sat = satellite;
    key.num = static_cast<int>(code);
    const auto found = state.kfIndexMap.find(key);
    if (found != state.kfIndexMap.end())
    {
        if (!state.stateTransitionMap.count(key) || found->second < 0 ||
            found->second >= state.x.size() || !std::isfinite(state.x(found->second)) ||
            !std::isfinite(state.P(found->second, found->second)) ||
            state.P(found->second, found->second) < 0) return false;
        row(found->second) = 1;
        return !acsConfig.zhangPppAr.user_adapter ||
            zhangUserCanonicalAmbiguityRow(state,receiver->id,satellite,static_cast<int>(code),row);
    }
    // Zero is legal only for the explicitly certified omitted chart reference.
    // A missing ordinary state or an unknown physical variable is never zero.
    const bool certifiedReference = acsConfig.zhangPppAr.user_adapter &&
        zhangPppArUserReferenceAmbiguity(state, receiver->id, satellite, code) &&
        zhangPppArUserAmbiguityIntegerValid(state, receiver->id, satellite, code);
    referenceZero = certifiedReference;
    return certifiedReference;
}

static auto& zhangUserSeparatedArcTrackers()
{
    static map<pair<double,std::size_t>, ZhangUserArSeparatedArcTracker> trackers;
    return trackers;
}

static vector<ZhangUserSharedEndpoint> buildZhangUserSharedEndpoints(
    Trace& trace, const KFState& state, Receiver* receiver, E_Sys system, GTime time,
    bool pairScopedComponents = false)
{
    vector<ZhangUserSharedEndpoint> result;
    if (!receiver || time == GTime::noTime() || time != state.time) return result;
    if (pairScopedComponents && !zhangUserGonsisSearchRequested())
        throw std::runtime_error("PAIR_SCOPED_COMPONENTS_REQUIRE_GONSIS_PROFILE");
    E_ObsCode first = E_ObsCode::NONE;
    E_ObsCode second = E_ObsCode::NONE;
    if (!zhangUserSharedCodes(system, first, second)) return result;
    const double geometryTolerance = acsConfig.assign_closest_epoch
        ? acsConfig.epoch_interval / 2 : acsConfig.epoch_tolerance;
    const string runtime = zhangAmbresRuntimeId(state);
    const auto reliability = zhangUserArReliabilityPolicy();
    ZhangUserArExperimentPolicy reliabilityArcPolicy;
    reliabilityArcPolicy.minimumArcSeconds = reliability.minimumArcSeconds;
    reliabilityArcPolicy.minimumArcSamples = reliability.minimumArcSamples;
    auto& separatedTrackers=zhangUserSeparatedArcTrackers();
    auto& separatedArcs = separatedTrackers.try_emplace(
        std::make_pair(reliability.minimumArcSeconds, reliability.minimumArcSamples),
        reliabilityArcPolicy).first->second;
    set<SatSys> seen;
    for (const auto& obs : only<GObs>(receiver->obsList))
    {
        if (obs.Sat.sys != system || !seen.insert(obs.Sat).second) continue;
        ZhangUserSharedEndpoint endpoint;
        endpoint.satellite = obs.Sat;
        endpoint.firstCode = first;
        endpoint.secondCode = second;
        endpoint.elevation = obs.satStat_ptr ? obs.satStat_ptr->el
            : std::numeric_limits<double>::quiet_NaN();
        const double observationAge = (obs.time - time).to_double();
        const auto elevation = zhangUserArElevationDecision(
            zhangUserArElevationEnabled(), acsConfig.ambrOpts.elevation_mask_deg,
            endpoint.elevation, obs.satStat_ptr != nullptr, observationAge, geometryTolerance);
        const bool baseValid = !obs.exclude && obs.satStat_ptr && obs.time != GTime::noTime() &&
            std::isfinite(observationAge) && std::abs(observationAge) <= geometryTolerance;
        bool valid = baseValid;
        bool physicallyObserved = baseValid;
        bool stateAvailable = true;
        bool productAvailable = true;
        string firstIdentity;
        string secondIdentity;
        string firstPhysical, secondPhysical, firstState, secondState, firstProduct, secondProduct;
        ZhangUserGonsisComponentCertificate firstComponent, secondComponent;
        for (const auto code : {first, second})
        {
            bool observed = false;
            bool physicalSignalPresent = false;
            unsigned long long trackingGeneration = 0;
            if (baseValid)
                for (const auto& [frequency, signal] : obs.sigs)
                {
                    if (signal.code != code || signal.invalid || !signal.L ||
                        !std::isfinite(signal.L)) continue;
                    const auto tracking = obs.satStat_ptr->sigStatMap.find(ft2string(frequency));
                    const auto phase = obs.satStat_ptr->sigStatMap.find(enum_to_string(code));
                    if (tracking != obs.satStat_ptr->sigStatMap.end())
                    {
                        physicalSignalPresent = true;
                        trackingGeneration = tracking->second.tracking.eventSequence;
                    }
                    if (!signal.P || !std::isfinite(signal.P) || (signal.LLI & 0x03)) continue;
                    if (tracking == obs.satStat_ptr->sigStatMap.end() ||
                        tracking->second.slip.any || tracking->second.savedSlip.any ||
                        tracking->second.phaseQuarantined || tracking->second.codeQuarantined ||
                        (phase != obs.satStat_ptr->sigStatMap.end() &&
                         (phase->second.phaseRejectCount || phase->second.phaseQuarantined))) continue;
                    observed = true;
                    trackingGeneration = tracking->second.tracking.eventSequence;
                    break;
                }
            string productIdentity = acsConfig.zhangPppAr.user_adapter
                ? zhangUserSharedInternalPhaseIdentity(state, receiver->id, obs.Sat, code, time,
                    pairScopedComponents, code == first ? &firstComponent : &secondComponent)
                : zhangUserSharedExternalPhaseIdentity(obs.Sat, code, time);
            auto& row = code == first ? endpoint.firstRow : endpoint.secondRow;
            auto& referenceZero = code == first ? endpoint.firstReferenceZero : endpoint.secondReferenceZero;
            const bool rowValid = zhangUserSharedEndpointRow(state, receiver, obs.Sat, code, row, referenceZero);
            valid = valid && observed && !productIdentity.empty() && rowValid;
            physicallyObserved = physicallyObserved && physicalSignalPresent;
            stateAvailable = stateAvailable && rowValid;
            productAvailable = productAvailable && !productIdentity.empty();
            const string generationKey = zhangUserArcGenerationKey(
                acsConfig.zhangPppAr.user_adapter, receiver->id, obs.Sat.id(), static_cast<int>(code));
            const string trackingIdentity = receiver->id + ":" + obs.Sat.id() + ":" +
                std::to_string(static_cast<int>(code)) + ":T" + std::to_string(trackingGeneration);
            string stateIdentity = receiver->id + ":" + obs.Sat.id() + ":" +
                std::to_string(static_cast<int>(code)) + ":G" +
                std::to_string(zhangUserArcGeneration(state.metaDataMap, generationKey)) + ":T" +
                std::to_string(trackingGeneration);
            // Internal product-coordinate retirement is not a physical break.
            // Consume BOTH independently maintained generations: otherwise a
            // same-product reset/re-add could inherit old state-age/FIX locks.
            // External WUM keeps its existing state-generation identity.
            if (acsConfig.zhangPppAr.user_adapter)
                stateIdentity += ":SG" + std::to_string(zhangUserArcGeneration(state.metaDataMap,
                    zhangUserArcGenerationKey(false,receiver->id,obs.Sat.id(),static_cast<int>(code))));
            auto& identity = code == first ? firstIdentity : secondIdentity;
            identity = zhangUserArEncodedField(stateIdentity) + zhangUserArEncodedField(productIdentity);
            (code == first ? firstPhysical : secondPhysical) = trackingIdentity;
            (code == first ? firstState : secondState) = stateIdentity;
            (code == first ? firstProduct : secondProduct) = productIdentity;
        }
        if (pairScopedComponents)
        {
            endpoint.componentDomain = acsConfig.zhangPppAr.user_adapter
                ? zhangUserGonsisDualComponentDomain(firstComponent, secondComponent)
                : "EXTERNAL_OSB_CURRENT_CALIBRATION";
            valid = valid && !endpoint.componentDomain.empty();
        }
        const string traceComponentDomain = endpoint.componentDomain;
        const string rawIdentity = valid
            ? zhangUserArEncodedField(firstIdentity) + zhangUserArEncodedField(secondIdentity) : string{};
        const auto legacyArc = zhangUserExperimentArcs().update(
            {runtime, receiver->id, obs.Sat.id()}, double(time.bigTime), rawIdentity, valid);
        const string physicalIdentity = zhangUserArEncodedField(firstPhysical) +
            zhangUserArEncodedField(secondPhysical);
        const string stateIdentity = zhangUserArEncodedField(firstState) +
            zhangUserArEncodedField(secondState);
        const string productIdentity = zhangUserArEncodedField(firstProduct) +
            zhangUserArEncodedField(secondProduct);
        const auto ages = separatedArcs.update({runtime, receiver->id, obs.Sat.id()},
            double(time.bigTime), physicalIdentity, physicallyObserved,
            productIdentity, productAvailable, valid, stateIdentity, stateAvailable);
        const bool separated = reliability.safe || reliability.gonsis;
        const auto& arc = separated ? ages.effective : legacyArc;
        if (arc.valid)
        {
            endpoint.identity = separated ? rawIdentity :
                rawIdentity + "|arc_start=" + std::to_string(arc.first);
            endpoint.firstIntegerIdentity = firstIdentity;
            endpoint.secondIntegerIdentity = secondIdentity;
            endpoint.physicalIdentity = physicalIdentity;
            endpoint.stateIdentity = stateIdentity;
            endpoint.productIdentity = productIdentity;
            endpoint.eligible = elevation.eligible && arc.eligible;
            endpoint.effectiveSamples = ages.effective.count;
            result.push_back(std::move(endpoint));
        }
        trace << "\nZHANG_USER_SHARED_ENDPOINT time=" << time.to_string(0)
              << " receiver=" << receiver->id << " satellite=" << obs.Sat.id()
              << " provider=" << (acsConfig.zhangPppAr.user_adapter ? "INTERNAL" : "EXTERNAL_OSB")
              << " first_code=" << enum_to_string(first) << " second_code=" << enum_to_string(second)
              << " first_certified_reference_zero=" << endpoint.firstReferenceZero
              << " second_certified_reference_zero=" << endpoint.secondReferenceZero
              << " reference_zero_authority=ACTUAL_USER_REFERENCE_AND_INTEGER_VALID_GETTERS"
              << " integer_authority=" << (pairScopedComponents ? "WITHIN_IDENTIFIED_INTEGER_CHART_PAIR" : "KF_CHART_REFERENCE_PAIR")
              << " component_domain=" << (traceComponentDomain.empty() ? "NONE" : traceComponentDomain)
              << " physical_valid=" << arc.valid << " search_eligible="
              << (arc.valid && elevation.eligible && arc.eligible)
              << " elevation_deg=" << zhangDiagnosticDouble(obs.satStat_ptr ? obs.satStat_ptr->el * R2D : 0)
              << " arc_age_seconds=" << zhangDiagnosticDouble(arc.age) << " arc_samples=" << arc.count
              << " physical_age_seconds=" << zhangDiagnosticDouble(ages.physical.age)
              << " state_age_seconds=" << zhangDiagnosticDouble(ages.state.age)
              << " product_age_seconds=" << zhangDiagnosticDouble(ages.product.age)
              << " effective_age_seconds=" << zhangDiagnosticDouble(ages.effective.age)
              << " effective_samples=" << ages.effective.count
              << " physical_reason=" << ages.physical.reason
              << " state_reason=" << ages.state.reason
              << " product_reason=" << ages.product.reason
              << " reliability_profile=" << reliability.name
              << " arc_reason=" << arc.reason << " elevation_reason=" << elevation.reason
              << " gate_scope=CANDIDATE_ONLY float_observations_mutated=0";
    }
    return result;
}
