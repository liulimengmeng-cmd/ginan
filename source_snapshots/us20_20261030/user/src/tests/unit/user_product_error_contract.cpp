#include "common/zhangUserProductError.hpp"
#include "common/receiver.hpp"
#include <Eigen/Cholesky>
#include <cstdlib>
#include <iostream>
#include <sstream>

void chunkFilter(Trace&, KFState&, KFMeas&, ReceiverMap&,
    std::map<std::string, FilterChunk>&, std::map<std::string, PooledTraceFile>&);

// Linked to the real production PEA/KF objects with --wrap=main. The oracle is
// an independently formed joint information-system inverse, not our helper or
// KF's sequential update. Compare the ENTIRE mean and covariance, including
// coordinate/ambiguity/product-error cross blocks. No NIS and no jitter.
extern "C" int __wrap_main()
{
    unsigned checks = 0;
    auto check = [&](bool okay, const char* text)
    {
        ++checks;
        if (!okay) throw std::runtime_error(std::string("FAIL ") + text);
    };
    auto at = [](long double seconds) { GTime time; time.bigTime = seconds; return time; };
    auto record = []()
    {
        ZhangUserConstantProductErrorRecord value;
        value.supported = true;
        value.reason = "INDEPENDENT_FIXTURE_RECORD";
        value.owner = "ZHANG_USER_PRODUCT_ERROR:G06:G:1:1";
        value.identity = "REC_A";
        value.satellite = SatSys("G06");
        value.code = 1;
        value.validFrom = 1000000000;
        value.validTo = 1000001000;
        value.variance = 4;
        return value;
    };
    try
    {
        const bool oldReceiverChunking = acsConfig.pppOpts.receiver_chunking;
        acsConfig.pppOpts.receiver_chunking = true;
        ReceiverMap fixtureReceivers;
        std::map<std::string, PooledTraceFile> fixtureTraceFiles;
        check(zhangParseUserProductErrorPolicy("") == ZhangUserProductErrorPolicy::WHITE, "default white");
        check(zhangParseUserProductErrorPolicy("SEGMENT_CONSTANT") == ZhangUserProductErrorPolicy::SEGMENT_CONSTANT, "explicit segment policy");
        const auto conditionalPolicy = ZhangUserProductErrorPolicy::EXTERNAL_FIXED_CONDITIONAL;
        check(zhangParseUserProductErrorPolicy("EXTERNAL_FIXED_CONDITIONAL") == conditionalPolicy,
            "external conditional mode explicit parser");
        check(std::string(zhangUserProductErrorPolicyName(conditionalPolicy)) == "EXTERNAL_FIXED_CONDITIONAL",
            "conditional mode has distinct receipt name");
        check(zhangParseUserProductErrorUnsupported("FAIL") == ZhangUserProductErrorUnsupported::FAIL, "explicit strict unsupported policy");
        bool invalid = false;
        try { zhangParseUserProductErrorPolicy("INTERPOLATED_CONSTANT"); }
        catch (const std::runtime_error&) { invalid = true; }
        check(invalid, "unsupported policy rejected");
        invalid = false;
        try { zhangParseUserProductErrorPolicy("external_fixed_conditional"); }
        catch (const std::runtime_error&) { invalid = true; }
        check(invalid, "conditional parser does not accept approximate spelling");

        const double conditionalBias = 0.0123;
        const auto fixedVariance = zhangUserExternalBiasVariance(
            conditionalPolicy, false, false, true, true, conditionalBias, 0.04);
        check(fixedVariance.fixedConditional && fixedVariance.rawInputVariance == 0.04 &&
            fixedVariance.modeledVariance == 0 && conditionalBias == 0.0123,
            "nominal unchanged raw BIA variance retained while R omits conditional OSB variance");
        for (const auto oldPolicy : {ZhangUserProductErrorPolicy::WHITE,
                                    ZhangUserProductErrorPolicy::SEGMENT_CONSTANT})
        {
            const auto unchanged = zhangUserExternalBiasVariance(
                oldPolicy, false, false, true, true, conditionalBias, 0.04);
            check(!unchanged.fixedConditional && unchanged.modeledVariance == 0.04 &&
                unchanged.rawInputVariance == 0.04, "old variance policies unchanged by new helper");
        }
        const auto estimated = zhangUserExternalBiasVariance(
            conditionalPolicy, false, true, false, false, conditionalBias, -1);
        check(!estimated.fixedConditional && estimated.modeledVariance == -1 &&
            estimated.rawInputVariance == -1, "explicit estimated bias retains owned noise sentinel");
        for (int fault = 0; fault < 6; ++fault)
        {
            invalid = false;
            try
            {
                zhangUserExternalBiasVariance(conditionalPolicy, fault == 0, false,
                    fault != 1, fault != 2,
                    fault == 3 ? std::numeric_limits<double>::quiet_NaN() : conditionalBias,
                    fault == 4 ? -1 : (fault == 5 ? std::numeric_limits<double>::infinity() : 0.04));
            }
            catch (const std::runtime_error&) { invalid = true; }
            check(invalid, "conditional zero cannot hide internal/missing/invalid nominal products");
        }
        KFState conditionalState;
        zhangUserValidateConditionalProductErrorState(conditionalPolicy, conditionalState);
        const auto conditionalTransitionCount = conditionalState.stateTransitionMap.size();
        KFMeasEntry conditionalMeasurement(&conditionalState);
        conditionalMeasurement.obsKey.type = KF::PHAS_MEAS;
        conditionalMeasurement.addNoiseEntry(conditionalMeasurement.obsKey, 1, 0.0025);
        KFKey conditionalBiasNoise;
        conditionalBiasNoise.type = KF::PHASE_BIAS;
        conditionalMeasurement.addNoiseEntry(conditionalBiasNoise, 1, fixedVariance.modeledVariance);
        check(!conditionalMeasurement.noiseElementMap.count(conditionalBiasNoise) &&
            conditionalMeasurement.noiseElementMap.at(conditionalMeasurement.obsKey) == 0.0025 &&
            conditionalState.stateTransitionMap.size() == conditionalTransitionCount &&
            !zhangUserRequiresJointProductErrorFilter(conditionalState),
            "conditional mode removes only OSB noise and creates no latent state or nominal correction");
        conditionalState.metaDataMap["ZHANG_USER_PRODUCT_ERROR_JOINT_FILTER_REQUIRED"] = "1";
        invalid = false;
        try { zhangUserValidateConditionalProductErrorState(conditionalPolicy, conditionalState); }
        catch (const std::runtime_error&) { invalid = true; }
        check(invalid, "cannot relabel a learned segment-error posterior even after latent retirement");
        zhangUserValidateConditionalProductErrorState(ZhangUserProductErrorPolicy::WHITE, conditionalState);

        KFKey coordinate, ambiguity;
        coordinate.type = KF::REC_POS;
        coordinate.str = "UNIT";
        ambiguity.type = KF::AMBIGUITY;
        ambiguity.str = "UNIT";
        ambiguity.Sat = SatSys("G06");
        ambiguity.num = 1;
        const auto constant = record();
        const KFKey errorKey = zhangUserProductErrorStateKey(constant);
        check(errorKey.type != KF::AMBIGUITY && errorKey.type == KF::CALC, "product error continuous only");
        check(zhangUserProductErrorKey(errorKey), "product state has owned namespace");
        VectorXd priorMean(3); priorMean << 2, -0.4, 0;
        MatrixXd priorP(3, 3); priorP << 9, 0.7, 0, 0.7, 1.5, 0, 0, 0, 4;
        MatrixXd h(3, 3); h << 1, 0.3, 1, 1, -0.5, 1, 1, 0.7, 1;
        VectorXd values(3); values << 5, 4, 6;
        VectorXd observationVariance(3); observationVariance << 0.25, 0.36, 0.49;
        std::ostringstream trace;

        KFState constantState;
        constantState.prefitOpts.sigma_check = false;
        constantState.prefitOpts.omega_test = false;
        constantState.postfitOpts.sigma_check = false;
        constantState.postfitOpts.omega_test = false;
        constantState.postfitOpts.max_iterations = 1;
        constantState.chiSquareTest.enable = false;
        constantState.assume_linearity = true;
        constantState.joseph_stabilisation = true;
        for (int epoch = 0; epoch < 3; ++epoch)
        {
            const auto time = at(constant.validFrom + epoch * 30);
            KFMeasEntry entry(&constantState);
            entry.obsKey.type = KF::PHAS_MEAS;
            entry.obsKey.str = "UNIT";
            entry.obsKey.Sat = SatSys("G06");
            entry.obsKey.num = epoch;
            InitialState posInit, ambInit;
            posInit.estimate = ambInit.estimate = true;
            posInit.x = 2; posInit.P = 9;
            ambInit.x = -0.4; ambInit.P = 1.5;
            if (epoch)
            {
                constantState.getKFValue(coordinate, posInit.x, nullptr, nullptr, false);
                constantState.getKFValue(ambiguity, ambInit.x, nullptr, nullptr, false);
            }
            entry.addDsgnEntry(coordinate, h(epoch, 0), posInit);
            entry.addDsgnEntry(ambiguity, h(epoch, 1), ambInit);
            const auto attached = zhangUserAttachConstantProductError(constantState, entry, constant, time, 1, trace);
            check(attached.priorInjected == (epoch == 0), "prior injected only once");
            check(constantState.initNoiseMap.count(errorKey) == (epoch == 0 ? 1u : 0u), "existing epoch no initial noise");
            entry.addNoiseEntry(entry.obsKey, 1, observationVariance(epoch));
            entry.setValue(values(epoch));
            check(constantState.stateTransition(trace, time), "real state transition succeeded");
            const int pos = constantState.kfIndexMap.at(coordinate);
            const int amb = constantState.kfIndexMap.at(ambiguity);
            if (!epoch) constantState.P(pos, amb) = constantState.P(amb, pos) = 0.7;
            KFMeasEntryList entries;
            entries.push_back(entry);
            KFMeas measurement(constantState, entries, time);
            check(std::abs(measurement.R(0, 0) - observationVariance(epoch)) < 1e-15, "bias std moved to prior not silently lost");
            check(measurement.H(0, constantState.kfIndexMap.at(errorKey)) == 1, "latent error has measurement design");
            std::map<std::string, FilterChunk> chunks;
            chunkFilter(trace, constantState, measurement, fixtureReceivers, chunks, fixtureTraceFiles);
            check(chunks.empty() || (chunks.size() == 1 && chunks.begin()->second.begX == 0 &&
                chunks.begin()->second.numX == constantState.x.size() && chunks.begin()->second.begH == 0 &&
                chunks.begin()->second.numH == measurement.H.rows()), "product uncertainty requires complete measurement chunk");
            check(constantState.filterKalman(trace, measurement, "/PRODUCT_ERROR_UNIT", false, &chunks) == KFFilterResult::COMMITTED, "actual KF update committed");
            const MatrixXd subH = h.topRows(epoch + 1);
            const MatrixXd inverseR = observationVariance.head(epoch + 1).cwiseInverse().asDiagonal();
            const MatrixXd information = priorP.inverse() + subH.transpose() * inverseR * subH;
            const MatrixXd expectedP = information.inverse();
            const VectorXd expectedMean = expectedP *
                (priorP.inverse() * priorMean + subH.transpose() * inverseR * values.head(epoch + 1));
            const std::array<int, 3> indices = {pos, amb, constantState.kfIndexMap.at(errorKey)};
            for (int i = 0; i < 3; ++i)
            {
                check(std::abs(constantState.x(indices[i]) - expectedMean(i)) < 2e-12, "full mean matches independent batch");
                for (int j = 0; j < 3; ++j)
                    check(std::abs(constantState.P(indices[i], indices[j]) - expectedP(i, j)) < 2e-12, "full covariance matches independent batch");
            }
            check(std::abs(constantState.P(pos, indices[2])) > 0.1, "coordinate-error cross covariance retained");
            check(std::abs(constantState.P(amb, indices[2])) > 0.001, "ambiguity-error cross covariance retained");
        }

        // Direct H support is NOT the criterion for a state's update. A
        // temporarily unobserved product error still changes through its P
        // cross block. Real chunkFilter must select the joint update even when
        // the latent error's H column is entirely zero for two further epochs.
        KFState unobserved(constantState);
        const std::array<KFKey, 3> modelKeys = {coordinate, ambiguity, errorKey};
        VectorXd unobservedPriorMean(3);
        MatrixXd unobservedPriorP(3, 3);
        for (int i = 0; i < 3; ++i)
        {
            unobservedPriorMean(i) = unobserved.x(unobserved.kfIndexMap.at(modelKeys[i]));
            for (int j = 0; j < 3; ++j)
                unobservedPriorP(i, j) = unobserved.P(unobserved.kfIndexMap.at(modelKeys[i]),
                    unobserved.kfIndexMap.at(modelKeys[j]));
        }
        MatrixXd unobservedH(2, 3); unobservedH << 1, 0.2, 0, 1, -0.4, 0;
        Vector2d unobservedY; unobservedY << 4.5, 5.2;
        Vector2d unobservedR; unobservedR << 0.16, 0.25;
        for (int epoch = 0; epoch < 2; ++epoch)
        {
            const auto time = at(constant.validFrom + 90 + epoch * 30);
            KFMeasEntry entry(&unobserved);
            entry.obsKey.type = KF::CODE_MEAS; entry.obsKey.str = "UNIT";
            entry.obsKey.Sat = SatSys("G06"); entry.obsKey.num = epoch;
            InitialState posInit, ambInit;
            posInit.estimate = ambInit.estimate = true;
            unobserved.getKFValue(coordinate, posInit.x, nullptr, nullptr, false);
            unobserved.getKFValue(ambiguity, ambInit.x, nullptr, nullptr, false);
            entry.addDsgnEntry(coordinate, unobservedH(epoch, 0), posInit);
            entry.addDsgnEntry(ambiguity, unobservedH(epoch, 1), ambInit);
            entry.addNoiseEntry(entry.obsKey, 1, unobservedR(epoch));
            entry.setValue(unobservedY(epoch));
            check(unobserved.stateTransition(trace, time), "unobserved product cross-P survives chunked prediction");
            KFMeasEntryList entries;
            entries.push_back(entry);
            KFMeas measurement(unobserved, entries, time);
            check(measurement.H.col(unobserved.kfIndexMap.at(errorKey)).norm() == 0,
                "test includes entirely unobserved H-zero product column");
            std::map<std::string, FilterChunk> chunks;
            chunkFilter(trace, unobserved, measurement, fixtureReceivers, chunks, fixtureTraceFiles);
            check(chunks.empty() || (chunks.size() == 1 && chunks.begin()->second.begX == 0 &&
                chunks.begin()->second.numX == unobserved.x.size() && chunks.begin()->second.begH == 0 &&
                chunks.begin()->second.numH == measurement.H.rows()), "H-zero product still requires complete joint filter chunk");
            check(unobserved.filterKalman(trace, measurement, "/PRODUCT_ERROR_H_ZERO_UNIT", false, &chunks)
                == KFFilterResult::COMMITTED, "H-zero product actual joint KF commit");
            const MatrixXd subset = unobservedH.topRows(epoch + 1);
            const MatrixXd precision = unobservedR.head(epoch + 1).cwiseInverse().asDiagonal();
            const MatrixXd expectedP = (unobservedPriorP.inverse() + subset.transpose() * precision * subset).inverse();
            const VectorXd expectedMean = expectedP * (unobservedPriorP.inverse() * unobservedPriorMean +
                subset.transpose() * precision * unobservedY.head(epoch + 1));
            for (int i = 0; i < 3; ++i)
            {
                const int index = unobserved.kfIndexMap.at(modelKeys[i]);
                check(std::abs(unobserved.x(index) - expectedMean(i)) < 2e-12,
                    "H-zero full mean matches independent two-epoch batch");
                for (int j = 0; j < 3; ++j)
                    check(std::abs(unobserved.P(index, unobserved.kfIndexMap.at(modelKeys[j])) - expectedP(i, j)) < 2e-12,
                        "H-zero full covariance matches independent two-epoch batch");
            }
        }
        check(std::abs(unobserved.x(unobserved.kfIndexMap.at(errorKey)) - unobservedPriorMean(2)) > 1e-3,
            "unobserved product-error mean actually updates through cross covariance");
        const int oldPosIndex = constantState.kfIndexMap.at(coordinate);
        const int oldAmbIndex = constantState.kfIndexMap.at(ambiguity);
        const Vector2d oldMean(constantState.x(oldPosIndex), constantState.x(oldAmbIndex));
        Matrix2d oldCovariance;
        oldCovariance << constantState.P(oldPosIndex, oldPosIndex), constantState.P(oldPosIndex, oldAmbIndex),
            constantState.P(oldAmbIndex, oldPosIndex), constantState.P(oldAmbIndex, oldAmbIndex);
        auto replacement = constant;
        replacement.identity = "REC_B";
        replacement.validFrom = constant.validFrom + 90;
        replacement.variance = 16;
        KFMeasEntry boundary(&constantState);
        const auto changed = zhangUserAttachConstantProductError(constantState, boundary, replacement,
            at(replacement.validFrom), 1, trace);
        check(changed.priorInjected && changed.retired == 1, "new record independent prior and old record retired");
        check(!constantState.stateTransitionMap.count(errorKey), "old owned error removed from transition");
        check(constantState.stateTransition(trace, at(replacement.validFrom)), "boundary transition marginalizes old error");
        const int posAfter = constantState.kfIndexMap.at(coordinate);
        const int ambAfter = constantState.kfIndexMap.at(ambiguity);
        check(std::abs(constantState.x(posAfter) - oldMean(0)) < 1e-14 &&
            std::abs(constantState.x(ambAfter) - oldMean(1)) < 1e-14, "retirement does not reset or fix physical states");
        for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j)
        {
            const int indices[2] = {posAfter, ambAfter};
            check(std::abs(constantState.P(indices[i], indices[j]) - oldCovariance(i, j)) < 1e-14,
                "retirement preserves full retained marginal covariance");
        }
        const int newError = constantState.kfIndexMap.at(changed.key);
        check(constantState.x(newError) == 0 && constantState.P(newError, newError) == 16,
            "record replacement does not inherit learned old error");
        check(constantState.P(posAfter, newError) == 0 && constantState.P(ambAfter, newError) == 0,
            "new independent record prior cross blocks not inherited");
        check(zhangUserRetireExpiredProductErrors(constantState, at(replacement.validTo), trace) == 1,
            "expired record retired at half-open boundary");
        check(constantState.stateTransition(trace, at(replacement.validTo)), "expiry marginalization transition succeeded");
        check(!constantState.kfIndexMap.count(changed.key), "no expired latent state accumulation");
        check(constantState.kfIndexMap.count(ambiguity), "product retirement never retires integer physical arc");
        check(constantState.x.size() == 3, "remaining ONE position ambiguity only");
        check(zhangUserRequiresJointProductErrorFilter(constantState),
            "joint-error lifetime marker remains sticky after all owned errors retire");

        // Independent WHITE model gives a different posterior from a constant
        // record: sigma does not average away in the shared-error model.
        MatrixXd whiteH = h.leftCols(2);
        const Matrix2d whitePrior = priorP.topLeftCorner<2, 2>();
        const MatrixXd whiteInverseR = (observationVariance.array() + 4).matrix().cwiseInverse().asDiagonal();
        const Matrix2d whiteP = (whitePrior.inverse() + whiteH.transpose() * whiteInverseR * whiteH).inverse();
        check((whiteP - oldCovariance).norm() > 0.1, "white and segment models are deliberately different");
        KFState whiteState;
        whiteState.prefitOpts.sigma_check = whiteState.prefitOpts.omega_test = false;
        whiteState.postfitOpts.sigma_check = whiteState.postfitOpts.omega_test = false;
        whiteState.postfitOpts.max_iterations = 1;
        whiteState.chiSquareTest.enable = false;
        whiteState.joseph_stabilisation = true;
        for (int epoch = 0; epoch < 3; ++epoch)
        {
            const auto time = at(constant.validFrom + epoch * 30);
            KFMeasEntry entry(&whiteState);
            entry.obsKey.type = KF::PHAS_MEAS; entry.obsKey.str = "UNIT";
            entry.obsKey.Sat = SatSys("G06"); entry.obsKey.num = epoch;
            InitialState posInit, ambInit;
            posInit.estimate = ambInit.estimate = true;
            posInit.x = 2; posInit.P = 9; ambInit.x = -0.4; ambInit.P = 1.5;
            entry.addDsgnEntry(coordinate, h(epoch, 0), posInit);
            entry.addDsgnEntry(ambiguity, h(epoch, 1), ambInit);
            entry.addNoiseEntry(entry.obsKey, 1, observationVariance(epoch));
            KFKey oldWhiteBias; oldWhiteBias.type = KF::PHASE_BIAS; oldWhiteBias.Sat = SatSys("G06");
            entry.addNoiseEntry(oldWhiteBias, 1, constant.variance);
            entry.setValue(values(epoch));
            check(whiteState.stateTransition(trace, time), "white actual state transition");
            const int pos = whiteState.kfIndexMap.at(coordinate);
            const int amb = whiteState.kfIndexMap.at(ambiguity);
            if (!epoch) whiteState.P(pos, amb) = whiteState.P(amb, pos) = 0.7;
            KFMeasEntryList entries;
            entries.push_back(entry);
            KFMeas measurement(whiteState, entries, time);
            check(measurement.R(0, 0) == observationVariance(epoch) + constant.variance, "white sigma remains in R each epoch");
            check(whiteState.filterKalman(trace, measurement, "/PRODUCT_ERROR_WHITE_UNIT") == KFFilterResult::COMMITTED, "white actual KF commit");
            const MatrixXd subset = whiteH.topRows(epoch + 1);
            const MatrixXd precision = (observationVariance.head(epoch + 1).array() + 4).matrix().cwiseInverse().asDiagonal();
            const Matrix2d expectedP = (whitePrior.inverse() + subset.transpose() * precision * subset).inverse();
            const Vector2d expectedMean = expectedP * (whitePrior.inverse() * priorMean.head(2) +
                subset.transpose() * precision * values.head(epoch + 1));
            const int indices[2] = {pos, amb};
            for (int i = 0; i < 2; ++i)
            {
                check(std::abs(whiteState.x(indices[i]) - expectedMean(i)) < 2e-12, "white full mean independent batch");
                for (int j = 0; j < 2; ++j)
                    check(std::abs(whiteState.P(indices[i], indices[j]) - expectedP(i, j)) < 2e-12, "white full covariance independent batch");
            }
            check(!whiteState.kfIndexMap.count(errorKey), "white compatibility route adds no latent state");
        }
        check(!zhangUserProductUncertaintyComplete(constantState), "no evidence not treated complete");
        zhangUserMarkProductUncertaintyIncomplete(constantState, "EXTERNAL_CROSS_SATELLITE_FREQUENCY_CLOCK_COVARIANCE");
        check(!zhangUserProductUncertaintyComplete(constantState), "marginal-only source cannot certify complete model");
        const KFState copied(constantState);
        check(!zhangUserProductUncertaintyComplete(copied), "KF branch preserves uncertainty limitation");

        // Direct record routing rejects unknown std, interpolated aggregates,
        // nonconstant slope and expired intervals. Do not silently fix sigma.
        const auto savedBiasMaps = biasMaps;
        BiasEntry input;
        input.tini = at(1000000000);
        input.tfin = at(1000001000);
        input.refTime = at(1000000500);
        input.Sat = SatSys("G06");
        input.cod1 = E_ObsCode::L1C;
        input.cod2 = E_ObsCode::NONE;
        input.measType = PHAS;
        input.source = "bsx";
        input.bias = 0.0123;
        input.var = 0.04;
        const std::string device = "G06:G";
        auto& records = biasMaps[PHAS][device][input.cod1][input.cod2];
        records.clear(); records[input.tini] = input;
        auto queried = zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            at(1000000030), input.bias, input.var, false);
        check(queried.supported && queried.variance == input.var, "valid immutable record source preserves BIA std");
        check(!zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            input.tfin, input.bias, input.var, false).supported, "interval end not inherited");
        records[input.tini].var = std::numeric_limits<double>::quiet_NaN();
        check(!zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            at(1000000030), input.bias, 0, false).supported, "missing variance not authorized zero constant");
        records[input.tini].var = 0;
        check(!zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            at(1000000030), input.bias, 0, false).supported, "legacy zero or missing std explicitly unsupported");
        const auto conditionalNominal = zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            at(1000000030), input.bias, 0, false, false);
        check(conditionalNominal.supported && conditionalNominal.variance == 0 &&
            conditionalNominal.reason.find("NOT_AUTHORIZED_AS_TRUTH") != std::string::npos,
            "conditional nominal may have unknown std without asserting known exact zero");
        check(!zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            input.tfin, input.bias, 0, false, false).supported,
            "conditional nominal also enforces exclusive validity end");
        records[input.tini].bias = std::numeric_limits<double>::quiet_NaN();
        check(!zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            at(1000000030), input.bias, 0, false, false).supported,
            "conditional omitted covariance cannot hide nonfinite nominal value");
        records[input.tini] = input;
        check(!zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            at(1000000030), input.bias + 0.1, input.var, false, false).supported,
            "conditional route cannot silently replace nominal applied by current model");
        BiasEntry inputCode = input;
        inputCode.cod1 = E_ObsCode::L1C; // GINAN's signal enum is shared by CODE and PHAS.
        inputCode.measType = CODE;
        biasMaps[CODE][device][inputCode.cod1][E_ObsCode::NONE][inputCode.tini] = inputCode;
        const auto conditionalCode = zhangUserDirectConstantBiasRecord(inputCode.Sat,
            inputCode.cod1, CODE, at(1000000030), inputCode.bias, inputCode.var, false, false);
        check(conditionalCode.supported && zhangUserExternalBiasVariance(conditionalPolicy,
            false, false, true, conditionalCode.supported, inputCode.bias, inputCode.var).modeledVariance == 0,
            "same conditional rule applies to real direct nominal CODE record");
        records[input.tini] = input; records[input.tini].slop = 1e-6;
        check(!zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            at(1000000030), input.bias, input.var, false).supported, "slope not frozen as constant");
        records[input.tini] = input;
        BiasEntry next = input; next.tini = input.tfin;
        next.tfin = at(1000002000); next.refTime = at(1000001500);
        records[next.tini] = next;
        check(!zhangUserDirectConstantBiasRecord(input.Sat, input.cod1, PHAS,
            at(1000000030), input.bias, input.var, false).supported, "interpolation constituent model required");
        biasMaps = savedBiasMaps;
        acsConfig.pppOpts.receiver_chunking = oldReceiverChunking;
        std::cout << "PASS user_product_error_contract checks=" << checks
                  << " full_mean_covariance_independent_batch=1 repeated_prior=0 jitter=0 NIS=0\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
