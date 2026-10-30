// User-only GONSIS UCDF strategy, included after the existing real product
// adapters. The uploaded GONSIS decisions are independently reimplemented;
// this is not its PPP engine, IF/MW path, stochastic model or observation QC.
#include "common/zhangUserGonsis.hpp"
#include "common/zhangUserMultiGnss.hpp"
#include "common/zhangUserGonsisOutput.hpp"
#include "common/zhangUserGonsisEpoch.hpp"
#include "pea/zhangUserRegional.hpp"

struct ZhangUserGonsisFloatObservations
{
    GTime time;
    std::uint64_t commitSequence = 0;
    map<KFKey,int> columns;
    VectorXd floatMean;
    MatrixXd floatCovariance;
    MatrixXd design;
    VectorXd residuals;
    VectorXd variances;
    MatrixXd measurementCovariance;
    vector<KFKey> measurements;
    vector<ZhangUserGonsisMeasurementKind> kinds;
    vector<int> sourceMeasurementRows;
    vector<string> measurementUnits;
    vector<string> measurementIdentities;
    Eigen::Index regionalIonoRows = 0;
    Eigen::Index expectedRegionalIonoRows = -1;
    Eigen::Index regionalTropRows = 0, expectedRegionalTropRows = -1;
    bool regionalMode = false;
    string regionalCaptureStatus = "NOT_CAPTURED";
    string captureFailureReason = "NOT_CAPTURED";
    bool valid = false;
};

static string zhangUserGonsisResidualRowIdentity(const KFKey& key, int sourceRow)
{
    // KFKey's display string can truncate long comments. Keep complete typed
    // fields length-framed, with the original final-measurement row ordinal.
    std::ostringstream out;
    auto field=[&](const string& value){out<<value.size()<<':'<<value;};
    field(std::to_string(sourceRow)); field(std::to_string(static_cast<int>(key.type)));
    field(key.str); field(std::to_string(static_cast<int>(key.Sat.sys)));
    field(std::to_string(key.Sat.prn)); field(std::to_string(key.num)); field(key.comment);
    return out.str();
}

// One epoch per authoritative runtime. No historic observation matrices or
// integers are retained. A capture replaces the previous epoch atomically.
static map<string,ZhangUserGonsisFloatObservations> zhangUserGonsisObservationSnapshots;
struct ZhangUserGonsisPublishedReference
{
    GTime time;
    ZhangUserGonsisReferenceHint hint;
};
using ZhangUserGonsisReferenceKey = tuple<string,string,E_Sys,string>;
static map<ZhangUserGonsisReferenceKey,ZhangUserGonsisPublishedReference> zhangUserGonsisReferences;
struct ZhangUserGonsisFixedLock
{
    string identity;
    int count = 0;
    GTime lastIncrement;
};
// Source lock increases after a published FIX, not on every tracking sample.
// One record per satellite/current identity; this stores no integer value.
static map<ZhangUserGonsisReferenceKey,map<string,ZhangUserGonsisFixedLock>> zhangUserGonsisFixedLocks;

void captureZhangUserGonsisFloatObservations(const KFState& state, const KFMeas& measurement)
{
    if (!zhangUserArReliabilityPolicy().gonsis) return;
    if (acsConfig.zhangFullRank.enable || acsConfig.zhangPppAr.output_products)
        throw std::runtime_error("GONSIS_UDUC_OBSERVATION_CAPTURE_IS_USER_ONLY");
    ZhangUserGonsisFloatObservations snapshot;
    snapshot.time = state.time;
    snapshot.commitSequence = state.factorCommitSequence;
    snapshot.columns = state.kfIndexMap;
    snapshot.floatMean = state.x;
    snapshot.floatCovariance = state.P;
    snapshot.regionalMode = zhangUserRegionalEnabled();
    const int count = measurement.obsKeys.size();
    const bool aligned = state.lastPppTransactionCommitted && state.lastPppQcConverged &&
        measurement.time == state.time && state.time != GTime::noTime() &&
        measurement.H.cols() == state.x.size() && measurement.H.rows() == count &&
        measurement.VV.size() == count && measurement.R.rows() == count && measurement.R.cols() == count;
    vector<int> selected;
    if (aligned)
    {
        snapshot.expectedRegionalIonoRows = 0; snapshot.expectedRegionalTropRows=0;
        if (snapshot.regionalMode && measurement.metaDataMaps.size()!=static_cast<size_t>(count))
            throw std::runtime_error("GONSIS_REGIONAL_FINAL_MEASUREMENT_METADATA_DIMENSION_MISMATCH");
        for (int i=0;i<count;++i)
        {
            const auto& key=measurement.obsKeys[i];
            const bool marker=i<measurement.metaDataMaps.size() &&
                measurement.metaDataMaps[i].contains(ZHANG_USER_REGIONAL_SD_RESIDUAL_MARKER) &&
                measurement.metaDataMaps[i].at(ZHANG_USER_REGIONAL_SD_RESIDUAL_MARKER)!=nullptr;
            const bool labelled=key.comment.starts_with(ZHANG_USER_REGIONAL_SD_METRE_LABEL);
            if (marker) ++snapshot.expectedRegionalIonoRows;
            bool regionalRow=false; string rowReason;
            if (!zhangUserGonsisRegionalResidualRowContract(snapshot.regionalMode,
                    key.type==KF::IONOSPHERIC,marker,labelled,regionalRow,rowReason))
                throw std::runtime_error(rowReason+":source_row="+std::to_string(i));
            const bool tropMarker=i<measurement.metaDataMaps.size() &&
                measurement.metaDataMaps[i].contains(ZHANG_USER_REGIONAL_TROP_MARKER) &&
                measurement.metaDataMaps[i].at(ZHANG_USER_REGIONAL_TROP_MARKER)!=nullptr;
            const bool tropLabel=key.comment.starts_with(ZHANG_USER_REGIONAL_TROP_LABEL);
            if (tropMarker) ++snapshot.expectedRegionalTropRows;
            if ((tropMarker || tropLabel) && (!snapshot.regionalMode || !tropMarker || !tropLabel || key.type!=KF::TROP))
                throw std::runtime_error("GONSIS_REGIONAL_TROP_ROW_IDENTITY_OR_UNIT_MISMATCH");
            if (key.type==KF::CODE_MEAS || key.type==KF::PHAS_MEAS || regionalRow || tropMarker)
                selected.push_back(i);
        }
    }
    snapshot.design.resize(selected.size(),state.x.size());
    snapshot.residuals.resize(selected.size());
    snapshot.variances.resize(selected.size());
    snapshot.measurementCovariance.resize(selected.size(),selected.size());
    for (int row=0;row<selected.size();++row)
    {
        const int old = selected[row];
        snapshot.design.row(row)=measurement.H.row(old);
        snapshot.residuals(row)=measurement.VV(old);
        snapshot.variances(row)=measurement.R(old,old);
        for (int column=0;column<selected.size();++column)
            snapshot.measurementCovariance(row,column)=measurement.R(old,selected[column]);
        snapshot.measurements.push_back(measurement.obsKeys[old]);
        const auto kind=measurement.obsKeys[old].type==KF::TROP ?
            ZhangUserGonsisMeasurementKind::REGIONAL_TROP_ZTD :
            measurement.obsKeys[old].type==KF::IONOSPHERIC ?
            ZhangUserGonsisMeasurementKind::REGIONAL_IONO_SD :
            measurement.obsKeys[old].type==KF::PHAS_MEAS ? ZhangUserGonsisMeasurementKind::PHASE :
            ZhangUserGonsisMeasurementKind::CODE;
        snapshot.kinds.push_back(kind);
        if (kind==ZhangUserGonsisMeasurementKind::REGIONAL_IONO_SD) ++snapshot.regionalIonoRows;
        if (kind==ZhangUserGonsisMeasurementKind::REGIONAL_TROP_ZTD) ++snapshot.regionalTropRows;
        snapshot.sourceMeasurementRows.push_back(old);
        snapshot.measurementUnits.push_back("m");
        snapshot.measurementIdentities.push_back(zhangUserGonsisResidualRowIdentity(measurement.obsKeys[old],old));
    }
    if (aligned && snapshot.regionalTropRows!=snapshot.expectedRegionalTropRows)
        throw std::runtime_error("GONSIS_REGIONAL_TROP_FINAL_CAPTURE_COUNT_MISMATCH");
    if (aligned && snapshot.regionalIonoRows!=snapshot.expectedRegionalIonoRows)
        throw std::runtime_error("GONSIS_REGIONAL_IONO_FINAL_CAPTURE_COUNT_MISMATCH");
    snapshot.regionalCaptureStatus = !aligned ? "FINAL_MEASUREMENT_NOT_ALIGNED_NO_COUNT_ASSERTION" :
        !snapshot.regionalMode ? "REGIONAL_MODE_OFF" :
        snapshot.regionalIonoRows ? "CURRENT_FINAL_REGIONAL_SD_ROWS_BOUND" :
        "CURRENT_FINAL_MEASUREMENT_CONFIRMS_NO_REGIONAL_PSEUDO";
    const auto covarianceFailure=selected.empty() ? string("EMPTY_MEASUREMENT_SET") :
        zhangUserGonsisDetail::covarianceFailure(snapshot.measurementCovariance);
    snapshot.valid = aligned && !selected.empty() && snapshot.floatMean.allFinite() && snapshot.floatCovariance.allFinite() &&
        snapshot.design.allFinite() && snapshot.residuals.allFinite() && snapshot.variances.allFinite() &&
        (snapshot.variances.array()>0).all() && covarianceFailure.empty();
    snapshot.captureFailureReason = snapshot.valid ? "NONE" : !aligned ? "FINAL_FLOAT_MEASUREMENT_NOT_ALIGNED" :
        !covarianceFailure.empty() ? "MEASUREMENT_"+covarianceFailure : "CAPTURE_MOMENTS_INVALID";
    zhangUserGonsisObservationSnapshots[zhangAmbresRuntimeId(state)] = std::move(snapshot);
}

static std::array<Eigen::Index,3> zhangUserGonsisPositionIndices(const KFState& state, const string& receiver)
{
    std::array<Eigen::Index,3> indices{-1,-1,-1};
    for (const auto& [key,index] : state.kfIndexMap)
        if (key.type == KF::REC_POS && key.str == receiver && key.num>=0 && key.num<3)
            indices[key.num]=index;
    return indices;
}

static void resolveGonsisUserWideLaneL1Par(Trace& trace, KFState& state)
{
    if (!zhangUserGonsisSearchRequested() || !zhangRatioOnly() ||
        acsConfig.ambrOpts.mode != E_ARmode::LAMBDA_ALT)
        throw std::runtime_error("GONSIS_UDUC_REQUIRES_USER_PROFILE_LAMBDA_ALT_RATIO_ONLY");
    const auto policy = zhangUserArReliabilityPolicy();
    static auto replayControls = ZhangUserGonsisReplayControls::fromEnvironment();
    static const bool exactIlsReuse = zhangUserGonsisReplaySwitch("ZHANG_USER_GONSIS_ILS_CACHE", true);
    static const bool targetFailedL1 = zhangUserGonsisReplaySwitch("ZHANG_USER_GONSIS_TARGETED_PAR", true);
    string replayEpoch = state.time.to_string(0);
    if (replayEpoch.size() >= 11) replayEpoch[10] = 'T';
    const string provider = acsConfig.zhangPppAr.user_adapter ? "internal" : "wum";
    const string runtime = zhangAmbresRuntimeId(state);
    const KFState* owner = state.factorCallbackOwnerOverride;
    if (!owner) throw std::runtime_error("GONSIS_UDUC_REQUIRES_AUTHORITATIVE_FLOAT_QUERY_OWNER");
    const string sourceRuntime = zhangAmbresRuntimeId(*owner);
    const auto found = zhangUserGonsisObservationSnapshots.find(sourceRuntime);
    const bool snapshotReady = found != zhangUserGonsisObservationSnapshots.end() && found->second.valid &&
        found->second.time == state.time && found->second.commitSequence == state.factorCommitSequence &&
        found->second.columns == state.kfIndexMap && found->second.floatMean.size() == state.x.size() &&
        found->second.floatCovariance.rows() == state.P.rows() && found->second.floatCovariance.cols() == state.P.cols() &&
        (found->second.floatMean.array() == state.x.array()).all() &&
        (found->second.floatCovariance.array() == state.P.array()).all();
    trace << "\nZHANG_USER_GONSIS_POLICY time=" << state.time.to_string(0)
          << " profile=GONSIS_UDUC provider=" << provider
          << " ratio_threshold=" << acsConfig.ambrOpts.ratioThres << " ratio_comparison=STRICT_GREATER"
          << " source_example_ratio=3 wl_residual_max_cycles=0.25 minimum_wl_pairs=" << policy.minimumParRank
          << " bootstrap_veto=0 conditional_sigma_veto=0 confirmation_epochs=1"
          << " float_position_average_sigma_max_m=0.5 fixed_position_sigma_norm_max_m=0.15"
          << " fixed_residual_sigma_max=4 wl_constraint_variance_cycles2=1e-8 if_constraint_variance_m2=9e-6"
          << " implemented_scope=GEC_DUAL_BAND_UC_CORE new_satellite_historical_ratio_fallback=NOT_PORTED"
          << " reference_promotion_counter=PUBLISHED_FIX_EPOCHS integer_counter_feedback=0"
          << " product_error_policy=" << zhangUserProductErrorPolicyName(zhangUserProductErrorPolicy())
          << " search_scope=CURRENT_PRODUCT_COMPONENT_PAIRS component_ils=INDEPENDENT"
          << " output_conditioning=ONE_JOINT_FULL_COVARIANCE_TRANSACTION"
          << " ils_reuse=" << (exactIlsReuse ? "BIT_IDENTICAL_MOMENTS_AND_TARGETS_CURRENT_DOMAIN_ONLY" : "DISABLED")
          << " par_removal=" << (targetFailedL1 ? "ACTUAL_FAILED_L1_ENDPOINTS" : "LEGACY_ALL_ACTIVE")
          << " replay_enabled=" << replayControls.enabled << " replay_selected_receivers_only=1"
          << " observation_updates_per_epoch=1 same_epoch_measurement_reuse=0 integer_history=0"
          << " risk_claim=UNCALIBRATED_MODEL_CONDITIONAL_NOT_FFRT nis_evaluations=0";
    using Group = pair<string,E_Sys>;
    map<Group,Receiver*> groups;
    for (const auto& [key,index] : state.kfIndexMap)
        if (key.type == KF::AMBIGUITY && key.rec_ptr && acsConfig.solve_amb_for[key.Sat.sys])
            groups[{key.str,key.Sat.sys}]=key.rec_ptr;
    struct Receipt
    {
        string receiver;
        E_Sys system;
        int eligibleDimension = 0;
        int pairs = 0;
        string reference = "NONE";
        string status;
        string componentDomain = "NONE";
        ZhangUserGonsisReferenceHint hint;
        vector<pair<string,string>> fixedIdentities;
    };
    vector<Receipt> receipts;
    vector<MatrixXd> allRows;
    vector<VectorXd> allValues, allVariances;
    bool replayJointOutput = false;
    for (const auto& [group,rec] : groups)
    {
        const auto& [receiver,system] = group;
        Receipt receipt{receiver,system};
        if (!state.lastPppTransactionCommitted || !state.lastPppQcConverged || !snapshotReady)
        {
            receipt.status = !state.lastPppTransactionCommitted || !state.lastPppQcConverged
                ? "FLOAT_QC_GATE" : "SAME_EPOCH_OBSERVATION_EVIDENCE_UNAVAILABLE";
            receipts.push_back(std::move(receipt));
            continue;
        }
        if (!zhangUserGecSystem(system))
        {
            receipt.status="GONSIS_UDUC_UNSUPPORTED_CONSTELLATION";
            receipts.push_back(std::move(receipt));
            continue;
        }
        // Track valid physical/state/product ages even while the FLOAT
        // position gate postpones ILS. This is metadata, not old integers.
        const auto allEndpoints=buildZhangUserSharedEndpoints(trace,state,rec,system,state.time,true);
        const auto position = zhangUserGonsisPositionIndices(state,receiver);
        const auto startGate = zhangValidateUserGonsisPosition(state.P,state.P,position);
        if (!std::isfinite(startGate.floatAverageSigma) || startGate.floatAverageSigma>0.5 ||
            position[0]<0 || position[1]<0 || position[2]<0)
        {
            receipt.status="FLOAT_POSITION_PRECISION_GATE";
            receipts.push_back(std::move(receipt));
            continue;
        }
        if (allEndpoints.empty())
        {
            receipt.status="INSUFFICIENT_LEGAL_DUAL_BAND_SATELLITES";
            receipts.push_back(std::move(receipt));
            continue;
        }
        const bool replayThisEpoch = replayControls.select(sourceRuntime,receiver,replayEpoch);
        replayJointOutput |= replayThisEpoch;
        vector<string> domainNames;
        for (const auto& endpoint : allEndpoints) domainNames.push_back(endpoint.componentDomain);
        const auto domains=zhangUserGonsisPartitionDomains(domainNames);
        // Reference/lock metadata contains no integer history. Retire inactive
        // domain versions rather than growing it for the entire product history.
        auto retireInactive = [&](auto& records)
        {
            for (auto entry=records.begin();entry!=records.end();)
            {
                const auto& key=entry->first;
                if (std::get<0>(key)==runtime && std::get<1>(key)==receiver &&
                    std::get<2>(key)==system && !domains.count(std::get<3>(key)))
                    entry=records.erase(entry);
                else ++entry;
            }
        };
        retireInactive(zhangUserGonsisReferences);
        retireInactive(zhangUserGonsisFixedLocks);
        for (const auto& [componentDomain,indices] : domains)
        {
            Receipt receipt{receiver,system};
            receipt.componentDomain=componentDomain;
            vector<ZhangUserSharedEndpoint> endpoints;
            for (const auto index : indices) endpoints.push_back(allEndpoints[index]);
            for (const auto& endpoint : endpoints)
                if (!zhangUserGecPair(system, endpoint.firstCode, endpoint.secondCode))
                    throw std::runtime_error("GONSIS_GEC_PRODUCT_SIGNAL_CONVENTION_MISMATCH");
            const double lambda1=zhangE27Wavelength(system,endpoints.front().firstCode);
            const double lambda2=zhangE27Wavelength(system,endpoints.front().secondCode);
            const auto referenceKey=std::make_tuple(runtime,receiver,system,componentDomain);
            auto& fixedLocks=zhangUserGonsisFixedLocks[referenceKey];
            MatrixXd physicalTransform(2*endpoints.size(),state.x.size());
            vector<ZhangUserGonsisSatellite> satellites;
            int eligible = 0;
            for (int i=0;i<endpoints.size();++i)
            {
                const auto& endpoint=endpoints[i];
                physicalTransform.row(2*i)=endpoint.firstRow.transpose();
                physicalTransform.row(2*i+1)=endpoint.secondRow.transpose();
                const double firstVariance=(endpoint.firstRow.transpose()*state.P*endpoint.firstRow)(0,0);
                const double secondVariance=(endpoint.secondRow.transpose()*state.P*endpoint.secondRow)(0,0);
                auto& lock=fixedLocks[endpoint.satellite.id()];
                if (lock.identity!=endpoint.identity ||
                    (lock.lastIncrement!=GTime::noTime() && state.time<lock.lastIncrement))
                    lock={endpoint.identity,0,GTime::noTime()};
                satellites.push_back({endpoint.satellite.id(),endpoint.identity,endpoint.elevation,
                    lock.count,lambda1*lambda1*firstVariance+lambda2*lambda2*secondVariance,
                    endpoint.eligible});
                eligible += endpoint.eligible;
            }
            receipt.eligibleDimension=2*std::max(0,eligible-1);
            ZhangUserGonsisOptions options;
            options.ratioThreshold=acsConfig.ambrOpts.ratioThres;
            options.minimumWideLanePairs=policy.minimumParRank;
            options.maximumTrials=policy.maximumParTrials;
            options.wavelength1Metres=lambda1;
            options.wavelength2Metres=lambda2;
            options.reuseIdenticalIls=exactIlsReuse;
            options.targetFailedL1SubsetForPar=targetFailedL1;
            const string replayContext="time="+replayEpoch+" receiver="+receiver+" provider="+provider+
                " component_domain="+componentDomain;
            int replayTrialCounter=0;
            if (replayThisEpoch)
            {
                trace << "\nZHANG_USER_GONSIS_REPLAY_FRAME " << replayContext
                      << " schema=1 integer_history=0 coordinate_order=PHYSICAL_INTERLEAVED_N1_N2"
                      << " input_state_dimension=" << state.x.size() << " physical_dimension=" << physicalTransform.rows();
                zhangUserGonsisReplayMatrix(trace,replayContext+" trial=0","full_FLOAT_mean",state.x);
                zhangUserGonsisReplayMatrix(trace,replayContext+" trial=0","full_FLOAT_Q",state.P);
                zhangUserGonsisReplayMatrix(trace,replayContext+" trial=0","physical_transform",physicalTransform);
                const auto columnPrecision=trace.precision();
                trace << std::setprecision(std::numeric_limits<double>::max_digits10);
                for (int index=0;index<endpoints.size();++index)
                    trace << "\nZHANG_USER_GONSIS_REPLAY_PHYSICAL_COLUMN " << replayContext
                          << " satellite=" << endpoints[index].satellite.id()
                          << " first_index=" << 2*index << " second_index=" << 2*index+1
                          << " first_code=" << enum_to_string(endpoints[index].firstCode)
                          << " second_code=" << enum_to_string(endpoints[index].secondCode)
                          << " eligible=" << endpoints[index].eligible
                          << " elevation_radians=" << satellites[index].elevationRadians
                          << " published_fix_samples=" << satellites[index].effectiveSamples
                          << " metric_variance_m2=" << satellites[index].metricVarianceMetres2
                          << " identity_hex=" << zhangUserGonsisReplayHex(endpoints[index].identity);
                trace << std::setprecision(columnPrecision);
                for (const auto& [key,index] : state.kfIndexMap)
                    trace << "\nZHANG_USER_GONSIS_REPLAY_STATE_COLUMN " << replayContext
                          << " index=" << index << " type=" << enum_to_string(key.type)
                          << " satellite=" << key.Sat.id() << " receiver_hex=" << zhangUserGonsisReplayHex(key.str)
                          << " number=" << key.num;
                options.replayTrial = [&](const ZhangUserGonsisTrial& detail,const ZhangUserGonsisReplayView& view)
                {
                    // Context and all matrices are synchronous borrowed views;
                    // no historic matrix captures or delayed callbacks exist.
                    const string context=replayContext+" trial="+std::to_string(++replayTrialCounter);
                    const auto precision=trace.precision();
                    trace << std::setprecision(std::numeric_limits<double>::max_digits10)
                          << "\nZHANG_USER_GONSIS_REPLAY_BEGIN " << context << " schema=1"
                          << " reference=" << detail.referenceSatellite << " reference_index=" << detail.referenceIndex
                          << " reference_identity_hex=" << zhangUserGonsisReplayHex(detail.referenceIdentity)
                          << " reason=" << detail.reason << " solver_reason="
                          << (detail.solverReason.empty() ? "NOT_RUN" : detail.solverReason)
                          << " complete=" << detail.ilsComplete << " ratio=" << detail.ratio
                          << " first_squared_cost=" << view.solution.bestSquaredCost
                          << " second_squared_cost=" << view.solution.secondSquaredCost
                          << " lambda1_m=" << lambda1 << " lambda2_m=" << lambda2
                          << " WL_variance_cycles2=" << options.wideLaneConstraintVarianceCycles2
                          << " IF_variance_m2=" << options.ifConstraintVarianceMetres2
                          << " ratio_threshold=" << options.ratioThreshold
                          << " WL_residual_max_cycles=" << options.wideLaneResidualCycles
                          << " minimum_WL_pairs=" << options.minimumWideLanePairs
                          << " maximum_trials=" << options.maximumTrials
                          << " targeted_PAR=" << options.targetFailedL1SubsetForPar
                          << " removed_index=" << detail.removedSatelliteIndex
                          << " removal_pool=" << detail.removalPool << " cache_reused=" << detail.ilsCacheReused
                          << " top_two_available=" << (view.solution.replayTopTwo.firstReduced.size()>0)
                          << " input_coordinates_valid=" << view.solution.replayTopTwo.inputCoordinatesValid;
                    trace << " active_indices=";
                    for (int index=0;index<detail.satelliteIndices.size();++index)
                        trace << (index ? "," : "") << detail.satelliteIndices[index];
                    trace << " WL_passed_indices=";
                    if (detail.wideLaneSatelliteIndices.empty()) trace << "EMPTY";
                    for (int index=0;index<detail.wideLaneSatelliteIndices.size();++index)
                        trace << (index ? "," : "") << detail.wideLaneSatelliteIndices[index];
                    trace << std::setprecision(precision);
                    zhangUserGonsisReplayMatrix(trace,context,"physical_frozen_mean",view.rawFrozenMean);
                    zhangUserGonsisReplayMatrix(trace,context,"physical_frozen_Q_raw",view.rawFrozenCovariance);
                    zhangUserGonsisReplayMatrix(trace,context,"physical_frozen_Q_used",view.usedFrozenCovariance);
                    MatrixXd candidates=MatrixXd::Zero(detail.wideLaneCandidateSatelliteIndices.size(),view.rawFrozenMean.size());
                    for (int row=0;row<detail.wideLaneCandidateSatelliteIndices.size();++row)
                    {
                        const int index=detail.wideLaneCandidateSatelliteIndices[row];
                        candidates(row,2*index)=1; candidates(row,2*index+1)=-1;
                        candidates(row,2*detail.referenceIndex)=-1; candidates(row,2*detail.referenceIndex+1)=1;
                    }
                    const MatrixXd candidateQ=candidates*view.usedFrozenCovariance*candidates.transpose();
                    const Eigen::Map<const VectorXd> candidateMeans(detail.wideLaneCandidateMeansCycles.data(),detail.wideLaneCandidateMeansCycles.size());
                    const Eigen::Map<const VectorXd> candidateIntegers(detail.wideLaneCandidateIntegers.data(),detail.wideLaneCandidateIntegers.size());
                    const Eigen::Map<const VectorXd> candidateResiduals(detail.wideLaneCandidateResidualsCycles.data(),detail.wideLaneCandidateResidualsCycles.size());
                    zhangUserGonsisReplayMatrix(trace,context,"WL_candidate_rows",candidates);
                    zhangUserGonsisReplayMatrix(trace,context,"WL_candidate_mean_actual_dot",candidateMeans);
                    zhangUserGonsisReplayMatrix(trace,context,"WL_candidate_Q_derived_no_sigma_gate",candidateQ);
                    zhangUserGonsisReplayMatrix(trace,context,"WL_candidate_integers_actual_round",candidateIntegers);
                    zhangUserGonsisReplayMatrix(trace,context,"WL_candidate_residual_actual",candidateResiduals);
                    zhangUserGonsisReplayMatrix(trace,context,"WL_passed_rows",view.wideLaneRows);
                    zhangUserGonsisReplayMatrix(trace,context,"WL_passed_integers",view.wideLaneIntegers);
                    zhangUserGonsisReplayMatrix(trace,context,"WL_conditioned_mean",view.wlConditionedMean);
                    zhangUserGonsisReplayMatrix(trace,context,"WL_conditioned_Q",view.wlConditionedCovariance);
                    zhangUserGonsisReplayMatrix(trace,context,"L1_target_rows",view.l1Rows);
                    zhangUserGonsisReplayMatrix(trace,context,"IF_target_rows",view.ifRows);
                    zhangUserGonsisReplayMatrix(trace,context,"L1_search_mean",view.l1Mean);
                    zhangUserGonsisReplayMatrix(trace,context,"L1_search_Q",view.l1Covariance);
                    zhangUserGonsisReplayMatrix(trace,context,"L1_best_integer_input",view.solution.integers);
                    const auto& topTwo=view.solution.replayTopTwo;
                    zhangUserGonsisReplayMatrix(trace,context,"LAMBDA_integer_basis",topTwo.integerBasis);
                    zhangUserGonsisReplayMatrix(trace,context,"LAMBDA_best_reduced",topTwo.firstReduced);
                    zhangUserGonsisReplayMatrix(trace,context,"LAMBDA_second_reduced",topTwo.secondReduced);
                    zhangUserGonsisReplayMatrix(trace,context,"LAMBDA_best_input_exact",topTwo.firstInput);
                    zhangUserGonsisReplayMatrix(trace,context,"LAMBDA_second_input_exact",topTwo.secondInput);
                    trace << "\nZHANG_USER_GONSIS_REPLAY_END " << context << " matrix_log=COMPLETE"
                          << " candidate_coordinate_contract=" << (!topTwo.inputCoordinatesValid ? "UNAVAILABLE" :
                              topTwo.firstInput.size()==view.solution.integers.size() &&
                              (topTwo.firstInput.array()==view.solution.integers.array()).all() ? "PASS" : "FAILED")
                          << " covariance_saturated_or_jittered=0 actual_ILS_options=LAMBDA_ALT_FULL_RANK_TOP2";
                };
            }
            auto solver = [&](const VectorXd& mean,const MatrixXd& covariance)
            {
                ZhangUserGonsisIlsSolution solution;
                GinAR_mtx matrix;
                matrix.aflt=mean; matrix.Paflt=covariance;
                GinAR_opt searchOptions;
                searchOptions.mode=E_ARmode::LAMBDA_ALT;
                searchOptions.ratthr=options.ratioThreshold;
                searchOptions.sucthr=acsConfig.ambrOpts.succsThres;
                searchOptions.nset=2;
                searchOptions.nitr=acsConfig.ambrOpts.AR_max_itr;
                searchOptions.min_lambda_fix_count=mean.size();
                searchOptions.max_lambda_fix_count=mean.size();
                searchOptions.force_bootstrap_success=false;
                searchOptions.maximum_lambda_conditional_sigma=0;
                searchOptions.lambda_candidate_nis_alpha=0;
                if (replayThisEpoch)
                    searchOptions.lambdaReplayTopTwo = [&](const ZhangUserGonsisLambdaTopTwo& detail)
                    { solution.replayTopTwo=detail; };
                GNSS_AR(trace,matrix,searchOptions);
                if (matrix.searchDiagnostic.localNisExecuted || matrix.lambda_candidate_nis_valid)
                    throw std::runtime_error("GONSIS_UDUC_UNEXPECTED_NIS_EXECUTION");
                if ((matrix.aflt.array()!=mean.array()).any() || (matrix.Paflt.array()!=covariance.array()).any())
                    throw std::runtime_error("GONSIS_UDUC_FROZEN_SEARCH_MOMENTS_MUTATED");
                solution.reason=matrix.searchDiagnostic.reason;
                solution.bestSquaredCost=matrix.searchDiagnostic.bestSquaredDistance;
                solution.secondSquaredCost=matrix.searchDiagnostic.secondSquaredDistance;
                solution.complete=matrix.searchDiagnostic.ilsComplete && matrix.searchDiagnostic.ratioExecuted &&
                    matrix.Ztrs.rows()==mean.size() && matrix.Ztrs.cols()==mean.size() && matrix.zfix.size()==mean.size();
                if (solution.complete)
                {
                    ZhangExactMatrix rows; ZhangExactVector values;
                    if (!zhangUserSearchExact(matrix.Ztrs,matrix.zfix,rows,values))
                        throw std::runtime_error("GONSIS_UDUC_LAMBDA_NONINTEGER_CHART");
                    const auto recovered=zhangExactAffineIntegerQuotient(rows,values,mean.size());
                    if (!recovered.valid || recovered.quotientRank ||
                        !zhangUserExactToDouble(recovered.particularSolution,solution.integers))
                        throw std::runtime_error("GONSIS_UDUC_LAMBDA_FULL_INTEGER_RECOVERY_FAILED");
                }
                const auto searchPrecision=trace.precision();
                trace << std::setprecision(std::numeric_limits<double>::max_digits10)
                      << "\nZHANG_USER_GONSIS_L1_SEARCH time=" << state.time.to_string(0)
                      << " receiver=" << receiver << " dimension=" << mean.size()
                      << " component_domain=" << componentDomain
                      << " complete=" << solution.complete << " reason=" << solution.reason
                      << " first_squared_distance=" << solution.bestSquaredCost
                      << " second_squared_distance=" << solution.secondSquaredCost
                      << " bootstrap_full_diagnostic=" << matrix.searchDiagnostic.bootstrapFull
                      << " bootstrap_best_1d_diagnostic=" << matrix.searchDiagnostic.bootstrapBest1d
                      << " bootstrap_veto=0 same_full_covariance=1 nis_evaluations=0";
                trace << std::setprecision(searchPrecision);
                return solution;
            };
            ZhangUserGonsisReferenceHint referenceHint;
            const auto oldReference=zhangUserGonsisReferences.find(referenceKey);
            if (oldReference!=zhangUserGonsisReferences.end())
            {
                const double age=(state.time-oldReference->second.time).to_double();
                if (age>0 && age<=std::max(45.0,1.5*acsConfig.epoch_interval))
                    referenceHint=oldReference->second.hint;
            }
            const auto result=zhangResolveUserGonsis(physicalTransform*state.x,
                physicalTransform*state.P*physicalTransform.transpose(),satellites,solver,options,referenceHint);
            const auto trialPrecision=trace.precision();
            for (int trial=0;trial<result.trials.size();++trial)
            {
                const auto& detail=result.trials[trial];
                trace << "\nZHANG_USER_GONSIS_TRIAL time=" << state.time.to_string(0)
                      << " receiver=" << receiver << " system=" << enum_to_string(system)
                      << " component_domain=" << componentDomain
                      << " trial=" << trial+1 << " reference=" << detail.referenceSatellite
                      << " satellites=" << detail.satelliteIndices.size()
                      << " passed_wl_pairs=" << detail.wideLaneSatelliteIndices.size()
                      << " l1_dimension=" << detail.requestedL1Rank << " complete=" << detail.ilsComplete
                      << " ratio=" << std::setprecision(17) << detail.ratio << std::setprecision(trialPrecision)
                      << " reason=" << detail.reason
                      << " next_removed_satellite=" << (detail.removedSatelliteIndex<0 ? "NONE"
                          : satellites[detail.removedSatelliteIndex].satellite)
                      << " removal_pool=" << detail.removalPool << " ils_cache_reused=" << detail.ilsCacheReused
                      << " retries_start_from_same_float=1 nis_evaluations=0";
                for (int i=0;i<detail.wideLaneCandidateSatelliteIndices.size();++i)
                    trace << "\nZHANG_USER_GONSIS_WL time=" << state.time.to_string(0)
                          << " receiver=" << receiver << " trial=" << trial+1
                          << " component_domain=" << componentDomain
                          << " reference=" << detail.referenceSatellite
                          << " satellite=" << satellites[detail.wideLaneCandidateSatelliteIndices[i]].satellite
                          << " residual_cycles=" << zhangDiagnosticDouble(detail.wideLaneCandidateResidualsCycles[i])
                          << " threshold_cycles=0.25 bootstrap_veto=0";
            }
            receipt.status=result.reason;
            trace << "\nZHANG_USER_GONSIS_SEARCH_REUSE time=" << state.time.to_string(0)
                  << " receiver=" << receiver << " component_domain=" << componentDomain
                  << " actual_ils_calls=" << result.ilsCalls << " exact_cache_hits=" << result.ilsCacheHits
                  << " cache_scope=CURRENT_FROZEN_DOMAIN_FULL_NUMERIC_KEY";
            if (result.accepted)
            {
                const int pairs=result.wideLaneIntegers.size();
                MatrixXd integerRows(2*pairs,physicalTransform.rows());
                VectorXd integerValues(2*pairs);
                integerRows.topRows(pairs)=result.wideLaneRows;
                integerRows.bottomRows(pairs)=result.l1Rows;
                integerValues.head(pairs)=result.wideLaneIntegers;
                integerValues.tail(pairs)=result.l1Integers;
                ZhangExactMatrix exactRows; ZhangExactVector exactValues;
                if (!zhangUserSearchExact(integerRows,integerValues,exactRows,exactValues) ||
                    !zhangAuditUserIntegerUnion(exactRows,exactValues,physicalTransform.rows()).valid)
                    throw std::runtime_error("GONSIS_UDUC_FINAL_INTEGER_UNION_INVALID");
                allRows.push_back(result.constraintRows*physicalTransform);
                allValues.push_back(result.constraintValues);
                allVariances.push_back(result.constraintVariances);
                receipt.pairs=pairs;
                receipt.reference=result.referenceSatellite;
                receipt.hint={result.referenceSatellite,result.referenceIdentity,true};
                for (const auto index : result.acceptedSatelliteIndices)
                    receipt.fixedIdentities.push_back({satellites[index].satellite,satellites[index].identity});
            }
            receipts.push_back(std::move(receipt));
        }
    }
    int rowCount=0;
    for (const auto& rows : allRows) rowCount+=rows.rows();
    bool committed=true;
    string outputFailure;
    if (rowCount)
    {
        MatrixXd rows(rowCount,state.x.size());
        VectorXd values(rowCount),variances(rowCount);
        int offset=0;
        for (int i=0;i<allRows.size();++i)
        {
            rows.middleRows(offset,allRows[i].rows())=allRows[i];
            values.segment(offset,allValues[i].size())=allValues[i];
            variances.segment(offset,allVariances[i].size())=allVariances[i];
            offset+=allRows[i].rows();
        }
        const string jointReplayContext="time="+replayEpoch+" receiver=JOINT provider="+provider+
            " component_domain=JOINT_OUTPUT trial=0";
        if (replayJointOutput)
        {
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_FLOAT_mean",state.x);
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_FLOAT_Q",state.P);
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_constraint_rows",rows);
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_constraint_values",values);
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_constraint_variances",variances);
        }
        const auto output=zhangConditionUserGonsisOutput(state.x,state.P,rows,values,variances);
        if (replayJointOutput)
        {
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_output_mean",output.mean);
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_output_Q",output.covariance);
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_constraint_innovation",output.innovation);
            zhangUserGonsisReplayMatrix(trace,jointReplayContext,"joint_constraint_residual",output.constraintResidual);
            trace << "\nZHANG_USER_GONSIS_REPLAY_OUTPUT " << jointReplayContext
                  << " conditioning_valid=" << output.valid << " conditioning_reason=" << output.reason;
        }
        committed=output.valid;
        if (!committed) outputFailure=output.reason;
        if (committed)
            for (const auto& receipt : receipts)
            {
                if (!receipt.pairs) continue;
                const auto position=zhangValidateUserGonsisPosition(state.P,output.covariance,
                    zhangUserGonsisPositionIndices(state,receipt.receiver));
                trace << "\nZHANG_USER_GONSIS_POSITION time=" << state.time.to_string(0)
                      << " receiver=" << receipt.receiver << " float_average_sigma_m=" << position.floatAverageSigma
                      << " fixed_sigma_norm_m=" << position.fixedSigmaNorm << " status=" << position.reason;
                if (!position.valid)
                { if (committed) outputFailure=position.reason; committed=false; }
            }
        // Record the same candidate's residual margins even if its position
        // precision gate rejected it. This is read-only evaluation, NOT a
        // second conditioning or an alternate FLOAT root.
        if (output.valid)
        {
            const auto& snapshot=found->second;
            if (replayJointOutput)
            {
                zhangUserGonsisReplayMatrix(trace,jointReplayContext,"retained_design",snapshot.design);
                zhangUserGonsisReplayMatrix(trace,jointReplayContext,"FLOAT_postfit_residual",snapshot.residuals);
                zhangUserGonsisReplayMatrix(trace,jointReplayContext,"retained_measurement_variances",snapshot.variances);
                zhangUserGonsisReplayMatrix(trace,jointReplayContext,"retained_full_measurement_R",snapshot.measurementCovariance);
            }
            const auto residualGate=zhangValidateUserGonsisResiduals(
                state.x,output.mean,snapshot.design,snapshot.residuals,snapshot.variances,snapshot.kinds,
                &snapshot.measurementCovariance,snapshot.expectedRegionalIonoRows,snapshot.expectedRegionalTropRows);
            const string residualContext="time="+state.time.to_string(0)+" receiver=JOINT attempt=EPOCH_ONLY";
            zhangTraceUserGonsisResidualGate(trace,residualContext,residualGate);
            trace << "\nZHANG_USER_GONSIS_REGIONAL_CAPTURE " << residualContext
                  << " status=" << snapshot.regionalCaptureStatus
                  << " regional_mode=" << snapshot.regionalMode
                  << " actual_trop_rows=" << snapshot.regionalTropRows << " expected_trop_rows=" << snapshot.expectedRegionalTropRows
                  << " actual_rows=" << snapshot.regionalIonoRows
                  << " expected_rows=" << snapshot.expectedRegionalIonoRows
                  << " marker_matched=" << (snapshot.regionalIonoRows==snapshot.expectedRegionalIonoRows)
                  << " same_final_measurement=1 units=m full_CP_PHASE_IONO_R=1 nis_evaluations=0";
            if (replayJointOutput)
            {
                zhangUserGonsisReplayMatrix(trace,jointReplayContext,"fixed_residual_all_rows_m",residualGate.fixedResidual);
                zhangUserGonsisReplayMatrix(trace,jointReplayContext,"measurement_sigmas_m",residualGate.measurementSigmas);
                zhangUserGonsisReplayMatrix(trace,jointReplayContext,"four_sigma_margins_m",residualGate.fourSigmaMarginsMetres);
            }
            if (!residualGate.valid)
            { if (committed) outputFailure=residualGate.reason; committed=false; }
        }
        if (committed)
        {
            state.dx=output.mean-state.x;
            state.x=output.mean; state.P=output.covariance;
        }
    }
    for (const auto& receipt : receipts)
    {
        const auto key=std::make_tuple(runtime,receipt.receiver,receipt.system,receipt.componentDomain);
        if (committed && receipt.pairs)
        {
            zhangUserGonsisReferences[key]={state.time,receipt.hint};
            for (const auto& [satellite,identity] : receipt.fixedIdentities)
            {
                auto& lock=zhangUserGonsisFixedLocks[key][satellite];
                if (lock.identity!=identity) lock={identity,0,GTime::noTime()};
                if (lock.lastIncrement!=state.time)
                {
                    lock.count=std::min(1000000,lock.count+1);
                    lock.lastIncrement=state.time;
                }
            }
        }
        else zhangUserGonsisReferences.erase(key);
        trace << "\nZHANG_USER_GONSIS_COMPONENT_SUMMARY time=" << state.time.to_string(0)
              << " receiver=" << receipt.receiver << " provider=" << provider
              << " system=" << enum_to_string(receipt.system) << " eligible_dimension=" << receipt.eligibleDimension
              << " component_domain=" << receipt.componentDomain
              << " carried_rank=0 newly_fixed_rank=" << (committed ? 2*receipt.pairs : 0)
              << " final_rank=" << (committed ? 2*receipt.pairs : 0)
              << " wl_rank=" << (committed ? receipt.pairs : 0) << " l1_rank=" << (committed ? receipt.pairs : 0)
              << " reference=" << receipt.reference
              << " status=" << (committed ? receipt.status : "OUTPUT_GATE_REJECTED:"+outputFailure)
              << " feedback=EPOCH_ONLY soft_integer_pseudomeasurements=1 integer_truth_certified=0"
              << " nis_evaluations=0 jitter=0 float_feedback=0";
    }
    // Preserve one public summary per receiver/system. All accepted components
    // were jointly conditioned ABOVE on the original full state covariance;
    // there is no sequential observation reuse or block-diagonal shortcut.
    map<Group,vector<const Receipt*>> groupedReceipts;
    for (const auto& receipt : receipts) groupedReceipts[{receipt.receiver,receipt.system}].push_back(&receipt);
    for (const auto& [group,details] : groupedReceipts)
    {
        int eligibleDimension=0,pairs=0,acceptedComponents=0;
        string reference="NONE";
        for (const auto* detail : details)
        {
            eligibleDimension+=detail->eligibleDimension;
            pairs+=detail->pairs;
            if (detail->pairs)
            {
                ++acceptedComponents;
                reference=acceptedComponents==1 ? detail->reference : "MULTI_COMPONENT";
            }
        }
        const string status=!committed ? "OUTPUT_GATE_REJECTED:"+outputFailure
            : details.size()==1 ? details.front()->status
            : pairs ? (acceptedComponents==details.size() ? "ALL_COMPONENTS_ACCEPTED" : "PARTIAL_COMPONENTS_ACCEPTED")
            : "ALL_COMPONENTS_REJECTED";
        trace << "\nZHANG_USER_GONSIS_AR_SUMMARY time=" << state.time.to_string(0)
              << " receiver=" << group.first << " provider=" << provider
              << " system=" << enum_to_string(group.second) << " eligible_dimension=" << eligibleDimension
              << " carried_rank=0 newly_fixed_rank=" << (committed ? 2*pairs : 0)
              << " final_rank=" << (committed ? 2*pairs : 0)
              << " wl_rank=" << (committed ? pairs : 0) << " l1_rank=" << (committed ? pairs : 0)
              << " reference=" << reference << " status=" << status
              << " component_count=" << details.size() << " accepted_components=" << (committed ? acceptedComponents : 0)
              << " feedback=EPOCH_ONLY soft_integer_pseudomeasurements=1 integer_truth_certified=0"
              << " nis_evaluations=0 jitter=0 float_feedback=0";
    }
    trace << "\nZHANG_USER_GONSIS_OUTPUT_TRANSACTION time=" << state.time.to_string(0)
          << " status=" << (committed ? "COMMITTED" : "ROLLED_BACK")
          << " committed_rank=" << (committed ? rowCount : 0) << " reason=" << (committed ? "PASS" : outputFailure)
          << " float_authority_unchanged=1 observation_update=0";
}
