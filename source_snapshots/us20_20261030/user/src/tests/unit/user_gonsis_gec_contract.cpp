#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/biases.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangUserGonsisEpoch.hpp"
#include "pea/zhangPppAr.hpp"
#include "common/zhangUserMultiGnss.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
bool sameMoments(const KFState& first, const KFState& second)
{
    return first.x.size() == second.x.size() && first.P.rows() == second.P.rows() &&
        first.P.cols() == second.P.cols() && first.dx.size() == second.dx.size() &&
        std::memcmp(first.x.data(), second.x.data(), first.x.size() * sizeof(double)) == 0 &&
        std::memcmp(first.P.data(), second.P.data(), first.P.size() * sizeof(double)) == 0 &&
        std::memcmp(first.dx.data(), second.dx.data(), first.dx.size() * sizeof(double)) == 0;
}

// These are synthetic committed FLOAT moments, not another implementation of
// the AR selector, WL/IF construction, ILS or output update. The actual public
// production entry point is linked below, without including its .inc source.
struct Fixture
{
    Receiver receiver;
    KFState state;
    KFMeas measurement;
    unsigned callbackCalls = 0;
    std::vector<KFKey> ambiguityKeys;

    explicit Fixture(unsigned id)
    {
        receiver.id = "GONSIS_UNIT_" + std::to_string(id);
        state.time.bigTime = 100000 + 30 * id;
        if (!bindZhangCheckpointRuntimeId(state, "gonsis-integration-" + std::to_string(id)))
            throw std::runtime_error("fixture runtime binding failed");
        const int dimension = 1 + 3 + 2 * 18;
        state.x = VectorXd::Zero(dimension);
        state.x(0) = 1;
        state.dx = VectorXd::Zero(dimension);
        state.P = MatrixXd::Identity(dimension, dimension) * 0.0001;
        state.P(0, 0) = 0;
        state.kfIndexMap.clear();
        state.kfIndexMap[state.oneKey] = 0;
        state.lastPppTransactionCommitted = true;
        state.lastPppQcConverged = true;
        state.factorCommitSequence = 77 + id;
        for (int axis = 0; axis < 3; ++axis)
        {
            KFKey key;
            key.type = KF::REC_POS;
            key.str = receiver.id;
            key.num = axis;
            key.rec_ptr = &receiver;
            state.kfIndexMap[key] = 1 + axis;
            state.stateTransitionMap[key][key][0] = 1;
            state.x(1 + axis) = axis == 0 ? 6378137 : 0;
        }
        KFMeasEntryList entries;
        for (int satelliteIndex = 0; satelliteIndex < 18; ++satelliteIndex)
        {
            const E_Sys system=satelliteIndex<6 ? E_Sys::GPS : satelliteIndex<12 ? E_Sys::GAL : E_Sys::BDS;
            const int local=satelliteIndex%6;
            const SatSys satellite(system, local+(system==E_Sys::BDS ? 19 : 1));
            const auto codes=zhangUserGecCodes(system);
            const E_FType second=system==E_Sys::GPS ? F2 : F5;
            auto& status = receiver.satStatMap[satellite];
            status.el = (70 - local * 5) * D2R;
            for (const auto frequency : {F1, second}) status.sigStatMap[ft2string(frequency)] = {};
            for (const auto code : {codes.first, codes.second})
                status.sigStatMap[enum_to_string(code)] = {};
            GObs observation;
            observation.Sat = satellite;
            observation.time = state.time;
            observation.satStat_ptr = &status;
            observation.sigs[F1].code = codes.first;
            observation.sigs[second].code = codes.second;
            observation.sigs[F1].L = 1e8;
            observation.sigs[second].L = 8e7;
            observation.sigs[F1].P = observation.sigs[second].P = 2e7;
            receiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(observation));
            for (int band = 0; band < 2; ++band)
            {
                const auto code = band ? codes.second : codes.first;
                KFKey ambiguity;
                ambiguity.type = KF::AMBIGUITY;
                ambiguity.str = receiver.id;
                ambiguity.Sat = satellite;
                ambiguity.num = static_cast<int>(code);
                ambiguity.rec_ptr = &receiver;
                const int column = 4 + 2 * satelliteIndex + band;
                state.kfIndexMap[ambiguity] = column;
                state.stateTransitionMap[ambiguity][ambiguity][0] = 1;
                state.x(column) = (band ? 85 + satelliteIndex : 100 + 2 * satelliteIndex) +
                    (band ? .0005 : .001) * satelliteIndex;
                ambiguityKeys.push_back(ambiguity);
                BiasEntry bias;
                bias.Sat = satellite;
                bias.cod1 = code;
                bias.cod2 = E_ObsCode::NONE;
                bias.measType = PHAS;
                bias.source = "bsx";
                bias.tini.bigTime = 99900;
                bias.tfin.bigTime = 200000;
                bias.refTime.bigTime = 150000;
                bias.bias = .123;
                bias.var = .004;
                biasMaps[PHAS][satellite.id() + ":" + satellite.sysChar()][code][E_ObsCode::NONE][bias.tini] = bias;
                for (const auto type : {KF::CODE_MEAS, KF::PHAS_MEAS})
                {
                    KFMeasEntry entry(&state);
                    entry.obsKey.type = type;
                    entry.obsKey.str = receiver.id;
                    entry.obsKey.Sat = satellite;
                    entry.obsKey.num = static_cast<int>(code);
                    entry.obsKey.rec_ptr = &receiver;
                    for (const auto& [key, index] : state.kfIndexMap)
                        if (key.type == KF::REC_POS)
                            entry.designEntryMap[key] = .2 + .1 * key.num;
                    if (type == KF::PHAS_MEAS)
                        entry.designEntryMap[ambiguity] = genericWavelength.at(band ? second : F1);
                    entry.addNoiseEntry(entry.obsKey, 1, type == KF::PHAS_MEAS ? .003 * .003 : .3 * .3);
                    entry.setValue(0);
                    entries.push_back(entry);
                }
            }
        }
        // Known nonzero cross-system and position/ambiguity covariance.
        VectorXd common=VectorXd::Zero(dimension);
        for(int i=1;i<dimension;++i) common(i)=0.0004*(1+i%7);
        state.P+=common*common.transpose();
        measurement = KFMeas(state, entries, state.time);
        measurement.V.setZero();
        measurement.VV.setZero();
        state.acceptedMeasurementFactorCallback = [&](auto&&...) { ++callbackCalls; return true; };
        state.stateTransitionFactorCallback = [&](auto&&...) { ++callbackCalls; return true; };
        state.exactStateTransformCallback = [&](auto&&...) { ++callbackCalls; return true; };
        state.authoritativeFactorTransactionCallback = [&](auto&&...) { ++callbackCalls; return true; };
    }
};
}

extern "C" int __wrap_main()
{
    unsigned checks = 0;
    auto check = [&](bool value, const char* label)
    {
        ++checks;
        if (!value) throw std::runtime_error(std::string("FAIL ") + label);
    };
    try
    {
        for (const char* name : {"ZHANG_USER_AR_MIN_ARC_SECONDS", "ZHANG_USER_AR_MIN_ARC_SAMPLES",
            "ZHANG_USER_AR_CONFIRM_EPOCHS", "ZHANG_USER_AR_MAX_CONDITIONAL_SIGMA_CYCLES",
            "ZHANG_USER_AR_PAR_MIN_RANK", "ZHANG_USER_AR_PAR_MAX_TRIALS",
            "ZHANG_USER_ALLOW_LEGACY_HOLD_RESEARCH"}) unsetenv(name);
        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE", "GONSIS_UDUC", 1);
        setenv("ZHANG_USER_AR_SEARCH", "GONSIS_UDUC", 1);
        setenv("ZHANG_USER_QC", "CAUSAL_TRACKING", 1);
        setenv("ZHANG_USER_AR_FEEDBACK", "EPOCH_ONLY", 1);
        setenv("ZHANG_USER_ENFORCE_AR_ELEVATION", "1", 1);
        setenv("ZHANG_R51_RATIO_ONLY", "1", 1);
        setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS", "0", 1);
        acsConfig.zhangPppAr.user_adapter = false;
        acsConfig.zhangPppAr.output_products = false;
        acsConfig.zhangFullRank.enable = false;
        acsConfig.phaseClockOsb.enable = false;
        acsConfig.ssrInOpts.one_freq_phase_bias = false;
        for (E_Sys system : {E_Sys::GPS,E_Sys::GAL,E_Sys::BDS})
        {
            const auto codes=zhangUserGecCodes(system);
            acsConfig.solve_amb_for[system]=true;
            acsConfig.code_priorities[system]={codes.first,codes.second};
        }
        acsConfig.ambrOpts.mode = E_ARmode::LAMBDA_ALT;
        acsConfig.ambrOpts.elevation_mask_deg = 15;
        acsConfig.ambrOpts.ratioThres = 2;
        acsConfig.ambrOpts.succsThres = .999;
        acsConfig.ambrOpts.lambda_set = 2;
        acsConfig.ambrOpts.AR_max_itr = 100;
        acsConfig.assign_closest_epoch = false;
        acsConfig.epoch_tolerance = .5;
        auto verifyIsolation = [&](const Fixture& fixture, const KFState& before, const KFState& output)
        {
            check(sameMoments(fixture.state, before), "authoritative FLOAT moments bitwise unchanged");
            check(fixture.state.kfIndexMap == before.kfIndexMap &&
                fixture.state.factorCommitSequence == before.factorCommitSequence &&
                fixture.state.metaDataMap == before.metaDataMap,
                "FLOAT layout sequence metadata unchanged");
            check(fixture.callbackCalls == 0, "private AR invoked no authoritative factor callbacks");
            check(output.factorCommitSequence == before.factorCommitSequence &&
                output.metaDataMap == before.metaDataMap,
                "output publishes no integer history or runtime bookkeeping changes");
            check(bool(output.acceptedMeasurementFactorCallback) &&
                bool(output.stateTransitionFactorCallback) && bool(output.exactStateTransformCallback) &&
                bool(output.authoritativeFactorTransactionCallback) && output.factorCallbackOwnerOverride == nullptr,
                "private bookkeeping callbacks restored on return");
        };
        Fixture success(1);
        const KFState before(success.state);
        captureZhangUserGonsisFloatObservations(success.state, success.measurement);
        KFState output(success.state);
        std::ostringstream successTrace;
        fixEpochOnlyUserAmbiguities(successTrace, output, success.state);
        check((output.x - before.x).norm() > 1e-5 && output.P.trace() < before.P.trace(),
            "production GONSIS changes only private ambiguity moments with soft covariance reduction");
        verifyIsolation(success, before, output);
        check(successTrace.str().find("GONSIS") != std::string::npos &&
            successTrace.str().find("nis_evaluations=0") != std::string::npos,
            "production strategy and zero NIS explicit in receipt");
        check(successTrace.str().find("final_rank=10 wl_rank=5 l1_rank=5") != std::string::npos &&
            successTrace.str().find("status=FULL_SUBSET_ACCEPTED") != std::string::npos &&
            successTrace.str().find("committed_rank=30 reason=PASS") != std::string::npos,
            "complete five-pair WL and L1 search passes the actual production output transaction");


        for (const char* system : {"GPS","GAL","BDS"})
            check(successTrace.str().find(std::string("system=")+system+" eligible_dimension=10")!=std::string::npos,
                "all three actual production integer groups searched");
        // Independent dense full-covariance conditioning oracle: known truth
        // integers and physical IF coefficients, no implementation helpers.
        MatrixXd h=MatrixXd::Zero(30,before.x.size());
        VectorXd truth=before.x;
        for(int i=4;i<truth.size();++i) truth(i)=std::round(truth(i));
        VectorXd variance(30);
        for(int sys=0;sys<3;++sys)
        {
            const double l1=299792458.0/1575.42e6;
            const double l2=299792458.0/(sys==0 ? 1227.60e6 : 1176.45e6);
            const double a=l2*l2/(l2*l2-l1*l1),b=-l1*l1/(l2*l2-l1*l1);
            const int reference=4+12*sys;
            for(int sat=1;sat<6;++sat)
            {
                const int row=10*sys+sat-1,column=reference+2*sat;
                h(row,column)=1;h(row,column+1)=-1;h(row,reference)=-1;h(row,reference+1)=1;
                h(row+5,column)=a*l1;h(row+5,column+1)=b*l2;
                h(row+5,reference)=-a*l1;h(row+5,reference+1)=-b*l2;
                variance(row)=1e-8;variance(row+5)=9e-6;
            }
        }
        MatrixXd s=h*before.P*h.transpose();s.diagonal()+=variance;
        const MatrixXd k=before.P*h.transpose()*s.inverse();
        const VectorXd expected=before.x+k*h*(truth-before.x);
        const MatrixXd expectedP=before.P-k*h*before.P;
        check((output.x-expected).cwiseAbs().maxCoeff()<1e-8,"joint GEC mean agrees with independent dense oracle");
        check((output.P-expectedP).cwiseAbs().maxCoeff()<1e-11,"joint GEC full covariance agrees with oracle");
        check((output.P.block(4,16,12,12)).norm()>1e-10,"cross-system covariance retained");
        check((output.x.head(4)-before.x.head(4)).norm()>1e-10,"joint ambiguity constraints propagate into position");
        std::cout<<"USER_GONSIS_GEC_CONTRACT PASS checks="<<checks
            <<" actual_production_ILS=1 independent_joint_oracle=1 separate_system_datums=1 full_cross_Q=1 FLOAT_bitwise_preserved=1\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return EXIT_FAILURE;}
}
