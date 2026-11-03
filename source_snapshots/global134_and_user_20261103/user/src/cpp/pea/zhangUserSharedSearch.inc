// Included inside ppp_ambres.cpp after the provider adapters. The two product
// providers retain different observation/covariance models, but share every
// integer search, partial-fix, exact-union and HOLD operation below.
#include "common/zhangUserIntegerContinuityAudit.hpp"
#include "common/zhangUserProductError.hpp"
struct ZhangUserSharedLedger
{
    vector<string> columns;
    ZhangExactMatrix rows;
    ZhangExactVector values;
    int discardedEqualityRank = 0;
    bool fullIntegerProjectionCertified = true;
};

static map<tuple<string, string, E_Sys>, ZhangUserSharedLedger> zhangUserSharedLedgers;
static ZhangUserIntegerContinuityAudit zhangUserSharedContinuityAudit;

static bool zhangUserSearchExact(const MatrixXd& matrix, const VectorXd& values,
    ZhangExactMatrix& rows, ZhangExactVector& rhs)
{
    if (matrix.rows() != values.size() || !matrix.allFinite() || !values.allFinite()) return false;
    constexpr double limit = 9007199254740992.0;
    rows.clear(); rhs.clear();
    for (int r = 0; r < matrix.rows(); ++r)
    {
        ZhangExactVector row(matrix.cols());
        for (int c = 0; c < matrix.cols(); ++c)
        {
            const double v = matrix(r,c);
            if (std::abs(v) > limit || std::abs(v - std::round(v)) > 1e-8) return false;
            row[c] = static_cast<long long>(std::llround(v));
        }
        const double v = values(r);
        if (std::abs(v) > limit || std::abs(v - std::round(v)) > 1e-8) return false;
        rows.push_back(std::move(row)); rhs.push_back(static_cast<long long>(std::llround(v)));
    }
    return true;
}

// Retiring physical arcs is deliberately conservative: keep exact equality
// consequences involving only surviving arcs. This is NOT advertised as a
// complete existential integer image (latent congruences can be retired).
// Candidate-mask changes instead use a temporary restriction, while the full
// persistent physical ledger remains intact.
static ZhangUserSharedLedger zhangUserReindexLedger(
    const ZhangUserSharedLedger& previous, const vector<string>& columns, int& retired)
{
    ZhangUserSharedLedger out;
    out.columns = columns;
    map<string, int> current;
    for (int c = 0; c < columns.size(); ++c) current[columns[c]] = c;
    vector<int> retained;
    for (int c = 0; c < previous.columns.size(); ++c)
        if (current.count(previous.columns[c])) retained.push_back(c);
    retired = previous.columns.size() - retained.size();
    if (previous.rows.empty()) return out;
    const auto projected = zhangRestrictUserIntegerLedger(
        previous.rows, previous.values, retained, previous.columns.size());
    if (!projected.valid) throw std::runtime_error("USER_HOLD_LEDGER_RETIREMENT_FAILED:" + projected.reason);
    out.discardedEqualityRank = projected.discardedRank;
    out.fullIntegerProjectionCertified = projected.knownFullProjection;
    for (const auto& row : projected.rows)
    {
        ZhangExactVector reordered(columns.size());
        for (int c = 0; c < retained.size(); ++c)
            reordered[current.at(previous.columns[retained[c]])] = row[c];
        out.rows.push_back(std::move(reordered));
    }
    out.values = projected.values;
    return out;
}

static void resolveSharedUserWideLaneL1Par(Trace& trace, KFState& state)
{
    if (!zhangRatioOnly())
        throw std::runtime_error("SHARED_USER_AR_REQUIRES_EXPLICIT_RATIO_ONLY_POLICY");
    if (acsConfig.ambrOpts.mode != E_ARmode::LAMBDA_ALT)
        throw std::runtime_error("SHARED_USER_AR_REQUIRES_LAMBDA_ALT_CONFIGURATION");
    const bool hold = zhangUserQuotientHoldRequested();
    const auto reliability = zhangUserArReliabilityPolicy();
    const string provider = acsConfig.zhangPppAr.user_adapter ? "internal" : "wum";
    const string runtime = zhangAmbresRuntimeId(state);
    trace << "\nZHANG_USER_SHARED_RELIABILITY_POLICY time=" << state.time.to_string(0)
          << " profile=" << reliability.name << " minimum_arc_seconds=" << reliability.minimumArcSeconds
          << " minimum_arc_samples=" << reliability.minimumArcSamples
          << " confirmation_epochs=" << reliability.confirmationEpochs
          << " confirmation_interpretation=CORRELATED_METADATA_NOT_INDEPENDENT_PROBABILITIES"
          << " bootstrap_success_gate=" << reliability.enforceBootstrapSuccess
          << " conditional_sigma_max_cycles=" << reliability.maximumConditionalSigmaCycles
          << " minimum_par_rank=" << reliability.minimumParRank
          << " maximum_par_trials=" << reliability.maximumParTrials
          << " minimum_support_satellites=" << reliability.minimumSupportSatellites
          << " product_uncertainty_complete=" << zhangUserProductUncertaintyComplete(state)
          << " risk_claim=UNCALIBRATED_MODEL_CONDITIONAL_NOT_FFRT nis_evaluations=0";
    map<tuple<string,string,E_Sys>, ZhangUserSharedLedger> pending;
    vector<MatrixXd> outputRows;
    vector<VectorXd> outputValues;
    vector<std::function<void(bool)>> receipts;
    using Group = pair<string,E_Sys>;
    map<Group,Receiver*> groups;
    for (const auto& [key,index] : state.kfIndexMap)
        if (key.type == KF::AMBIGUITY && key.rec_ptr && acsConfig.solve_amb_for[key.Sat.sys])
            groups[{key.str,key.Sat.sys}] = key.rec_ptr;
    auto summary = [&](const string& receiver, E_Sys system, int dimension, int carried,
                       int finalRank, int retired, const string& reference, const string& status)
    {
        receipts.push_back([&,receiver,system,dimension,carried,finalRank,retired,reference,status](bool committed)
        {
        trace << "\nZHANG_USER_SHARED_AR_SUMMARY time=" << state.time.to_string(0)
              << " receiver=" << receiver << " provider=" << provider
              << " system=" << enum_to_string(system) << " eligible_dimension=" << dimension
              << " carried_rank=" << (committed ? carried : 0)
              << " newly_fixed_rank=" << (committed ? finalRank-carried : 0)
              << " final_rank=" << (committed ? finalRank : 0)
              << " quotient_dimension=" << (committed ? dimension-carried : dimension)
              << " retired_columns=" << retired << " reference=" << reference
              << " feedback=" << (hold ? "QUOTIENT_HOLD" : "EPOCH_ONLY")
              << " status=" << (committed ? status : "OUTPUT_TRANSACTION_ROLLED_BACK")
              << " nis_evaluations=0 jitter=0 float_feedback=0";
        });
    };
    if (!state.lastPppTransactionCommitted || !state.lastPppQcConverged)
    {
        // No uncommitted physical state may refresh or consume held evidence.
        if (groups.empty()) summary("NONE", E_Sys::GPS,0,0,0,0,"NONE","FLOAT_QC_GATE");
        for (const auto& [group,rec] : groups)
            summary(group.first,group.second,0,0,0,0,"NONE","FLOAT_QC_GATE");
        for (const auto& receipt:receipts) receipt(true);
        return;
    }
    if (groups.empty()) summary("NONE",E_Sys::GPS,0,0,0,0,"NONE","NO_AMBIGUITY_STATES");
    for (const auto& [group,rec] : groups)
    {
        const auto& [receiver,system] = group;
        const auto endpoints = buildZhangUserSharedEndpoints(trace,state,rec,system,state.time);
        const auto key = std::make_tuple(runtime,receiver,system);
        vector<string> catalogue;
        vector<string> columnOwners;
        for (const auto& endpoint : endpoints)
        {
            catalogue.push_back(endpoint.firstIntegerIdentity + "/L1");
            catalogue.push_back(endpoint.secondIntegerIdentity + "/L2");
            const string stableOwner = zhangUserArEncodedField(receiver) +
                zhangUserArEncodedField(endpoint.satellite.id());
            columnOwners.push_back(stableOwner + ":" + std::to_string(static_cast<int>(endpoint.firstCode)));
            columnOwners.push_back(stableOwner + ":" + std::to_string(static_cast<int>(endpoint.secondCode)));
        }
        int retired = 0;
        ZhangUserSharedLedger ledger;
        ledger.columns = catalogue;
        if (hold)
        {
            const auto old = zhangUserSharedLedgers.find(key);
            if (old != zhangUserSharedLedgers.end())
                ledger = zhangUserReindexLedger(old->second,catalogue,retired);
        }
        trace << "\nZHANG_USER_SHARED_RETIREMENT time=" << state.time.to_string(0)
              << " receiver=" << receiver << " role=PHYSICAL_RETIREMENT"
              << " retired_columns=" << retired << " retained_equality_rank=" << ledger.rows.size()
              << " discarded_equality_rank=" << ledger.discardedEqualityRank
              << " full_integer_projection_certified=" << ledger.fullIntegerProjectionCertified;
        vector<int> candidates;
        for (int i = 0; i < endpoints.size(); ++i)
            if (endpoints[i].eligible) candidates.push_back(i);
        if (candidates.size() < 2)
        {
            // Missing/young candidate epochs interrupt confirmation. Audit
            // still retires identities; no hidden same-arc integer feedback.
            zhangUserSharedContinuityAudit.observe(runtime + ":" + receiver + ":" + enum_to_string(system),
                double(state.time.bigTime), catalogue, {}, {}, reliability.confirmationEpochs,
                reliability.contradictionIsolation, columnOwners);
            if (hold) pending[key] = ledger;
            summary(receiver,system,0,0,0,retired,"NONE","INSUFFICIENT_LEGAL_DUAL_BAND_SATELLITES");
            continue;
        }
        const int reference = *std::min_element(candidates.begin(),candidates.end(),
            [&](int a, int b)
            {
                const double va = (endpoints[a].firstRow.transpose()*state.P*endpoints[a].firstRow)(0,0)
                    + (endpoints[a].secondRow.transpose()*state.P*endpoints[a].secondRow)(0,0);
                const double vb = (endpoints[b].firstRow.transpose()*state.P*endpoints[b].firstRow)(0,0)
                    + (endpoints[b].secondRow.transpose()*state.P*endpoints[b].secondRow)(0,0);
                return va != vb ? va < vb : endpoints[a].satellite < endpoints[b].satellite;
            });
        // The exact domain includes EVERY still-valid physical endpoint. AR
        // elevation is a target mask, not permission to eliminate latent
        // integers: 2*m+l=1 still requires l to be odd when m is not searched.
        vector<int> targets;
        map<int,int> targetIndex;
        for (int endpoint=0; endpoint<endpoints.size(); ++endpoint)
            if (endpoint != reference)
            { targetIndex[endpoint]=targets.size(); targets.push_back(endpoint); }
        vector<int> candidateTargets;
        for (const int candidate : candidates)
            if (candidate != reference) candidateTargets.push_back(targetIndex.at(candidate));
        const int d = targets.size(), ambient = 2*d;
        const int candidateDimension = 2*candidateTargets.size();
        // y=(between-satellite WL, between-satellite L1). This is a primitive
        // integer chart; the common reference removes both receiver datums.
        ZhangExactMatrix physical(ambient,ZhangExactVector(catalogue.size()));
        MatrixXd transform = MatrixXd::Zero(ambient,state.x.size());
        for (int i = 0; i < d; ++i)
        {
            const int s = targets[i];
            physical[i][2*s]=1; physical[i][2*reference]=-1;
            physical[i][2*s+1]=-1; physical[i][2*reference+1]=1;
            physical[d+i][2*s]=1; physical[d+i][2*reference]=-1;
            transform.row(d+i) = (endpoints[s].firstRow-endpoints[reference].firstRow).transpose();
            transform.row(i) = transform.row(d+i)
                - (endpoints[s].secondRow-endpoints[reference].secondRow).transpose();
        }
        const int masked = 2*(endpoints.size()-candidates.size());
        trace << "\nZHANG_USER_SHARED_RETIREMENT time=" << state.time.to_string(0)
              << " receiver=" << receiver << " role=TEMPORARY_CANDIDATE_MASK"
              << " masked_columns=" << masked
              << " discarded_equality_rank=0 full_integer_projection_certified=1"
              << " persistent_ledger_discarded=0 latent_columns_retained=1";
        ZhangExactMatrix held;
        ZhangExactVector heldValues;
        for (int r=0;r<ledger.rows.size();++r)
        {
            const auto member = zhangIntegerRowLatticeContains(physical,ledger.rows[r]);
            if (!member.contained)
                throw std::runtime_error("USER_HOLD_ROW_OUTSIDE_LEGAL_SD_LATTICE");
            held.push_back(member.combination); heldValues.push_back(ledger.values[r]);
        }
        const auto baseline = zhangAuditUserIntegerUnion(held,heldValues,ambient);
        if (!baseline.valid) throw std::runtime_error("USER_HOLD_INFEASIBLE:"+baseline.reason);
        held=baseline.rows; heldValues=baseline.values;
        ZhangExactMatrix candidateSelectors(candidateDimension,ZhangExactVector(ambient));
        for (int i=0;i<candidateTargets.size();++i)
        {
            candidateSelectors[i][candidateTargets[i]]=1;
            candidateSelectors[candidateTargets.size()+i][d+candidateTargets[i]]=1;
        }
        auto candidateImage = zhangProjectUserCandidateIntegerImage(
            zhangExactAffineIntegerQuotient(held,heldValues,ambient),candidateSelectors);
        if (!candidateImage.valid)
            throw std::runtime_error("USER_CANDIDATE_DOMAIN_INVALID:"+candidateImage.reason);
        const int carried = candidateImage.fixedRank;
        const VectorXd mean = transform*state.x;
        const MatrixXd covariance = transform*state.P*transform.transpose();
        // Never take the next epoch's float search from a previous fixed KF.
        auto searchStage = [&](const string& stageName,const ZhangExactMatrix& selectors)
        {
            if (selectors.empty()) return;
            const auto residual = zhangBuildUserResidualQuotient(mean,covariance,held,heldValues);
            if (!residual.valid)
            {
                trace << "\nZHANG_USER_SHARED_STAGE time=" << state.time.to_string(0)
                      << " receiver=" << receiver << " stage=" << stageName
                      << " status=REJECTED reason=" << residual.reason;
                return;
            }
            const auto image = zhangR47CompileProductSearchFrame(selectors,held,heldValues,ambient);
            if (!image.valid) throw std::runtime_error("USER_STAGE_INTEGER_IMAGE_INVALID:"+image.reason);
            if (!image.searchRank)
            {
                trace << "\nZHANG_USER_SHARED_STAGE time=" << state.time.to_string(0)
                      << " receiver=" << receiver << " stage=" << stageName
                      << " search_dimension=0 status=ALREADY_DETERMINED";
                return;
            }
            MatrixXd projector; VectorXd offset;
            if (!zhangUserExactToDouble(image.projector,ambient,projector) ||
                !zhangUserExactToDouble(image.offsets,offset))
                throw std::runtime_error("USER_STAGE_IMAGE_COEFFICIENT_OVERFLOW");
            GinAR_mtx search;
            search.aflt = projector*residual.conditionalMean+offset;
            search.Paflt = projector*residual.conditionalCovariance*projector.transpose();
            search.Paflt = (0.5*(search.Paflt+search.Paflt.transpose())).eval();
            Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(search.Paflt);
            const bool positive = eigen.info()==Eigen::Success && eigen.eigenvalues().allFinite()
                && eigen.eigenvalues().minCoeff()>std::max(1e-14,
                    std::max(1.0,eigen.eigenvalues().cwiseAbs().maxCoeff())*
                    std::numeric_limits<double>::epsilon()*image.searchRank*64);
            if (!positive)
            {
                trace << "\nZHANG_USER_SHARED_STAGE time=" << state.time.to_string(0)
                      << " receiver=" << receiver << " stage=" << stageName
                      << " status=REJECTED reason=RESIDUAL_IMAGE_NOT_PD jitter=0";
                return;
            }
            GinAR_opt options;
            options.mode=E_ARmode::LAMBDA_ALT;
            options.sucthr=acsConfig.ambrOpts.succsThres;
            options.ratthr=acsConfig.ambrOpts.ratioThres;
            options.nset=acsConfig.ambrOpts.lambda_set;
            options.nitr=acsConfig.ambrOpts.AR_max_itr;
            options.min_lambda_fix_count=1;
            options.lambda_candidate_nis_alpha=0;
            MatrixXd physicalRows;
            if (!zhangUserExactToDouble(physical, catalogue.size(), physicalRows))
                throw std::runtime_error("USER_PHYSICAL_SUPPORT_COEFFICIENT_OVERFLOW");
            const MatrixXd physicalSupportTransform = projector * physicalRows;
            const int count=zhangResolveUserSharedPar(trace,search,options,nullptr,&physicalSupportTransform);
            trace << "\nZHANG_USER_SHARED_STAGE time=" << state.time.to_string(0)
                  << " receiver=" << receiver << " stage=" << stageName
                  << " search_dimension=" << image.searchRank << " accepted_rows=" << count
                  << " ratio_threshold=" << options.ratthr << " full_covariance=1 nis_evaluations=0";
            if (count<=0) return;
            ZhangExactMatrix accepted; ZhangExactVector acceptedValues;
            if (!zhangUserSearchExact(search.Ztrs,search.zfix,accepted,acceptedValues))
                throw std::runtime_error("USER_SEARCH_NONINTEGER_DECISION");
            const auto newRows=zhangExactMultiply(accepted,image.projector);
            const auto shifts=zhangExactMatrixTimesColumn(accepted,image.offsets);
            for (int i=0;i<acceptedValues.size();++i) acceptedValues[i]-=shifts[i];
            auto unionRows=held; auto unionValues=heldValues;
            unionRows.insert(unionRows.end(),newRows.begin(),newRows.end());
            unionValues.insert(unionValues.end(),acceptedValues.begin(),acceptedValues.end());
            const auto united=zhangAuditUserIntegerUnion(unionRows,unionValues,ambient);
            if (!united.valid) throw std::runtime_error("USER_SEARCH_AFFINE_UNION_REJECTED:"+united.reason);
            held=united.rows; heldValues=united.values;
        };
        ZhangExactMatrix wl(candidateTargets.size(),ZhangExactVector(ambient));
        for(int i=0;i<candidateTargets.size();++i) wl[i][candidateTargets[i]]=1;
        searchStage("WL",wl);
        // This is a user-search prerequisite, NOT a product certificate.
        // In a confirmed feasible integer domain, 2*WL=6 determines WL=3
        // even if no integer row-combination witness equals the unit WL row.
        const auto wlImage=zhangProjectUserCandidateIntegerImage(
            zhangExactAffineIntegerQuotient(held,heldValues,ambient),wl);
        if(!wlImage.valid) throw std::runtime_error("USER_WL_TARGET_IMAGE_INVALID:"+wlImage.reason);
        ZhangExactMatrix first;
        for(int i=0;i<candidateTargets.size();++i)
        {
            const auto& generators=wlImage.imageGenerators[i];
            const bool determined=std::all_of(generators.begin(),generators.end(),
                [](const ZhangExactInteger& value){return value==0;});
            if(determined)
            { ZhangExactVector row(ambient); row[d+candidateTargets[i]]=1; first.push_back(std::move(row)); }
        }
        searchStage("L1",first);
        auto final = zhangAuditUserIntegerUnion(held,heldValues,ambient);
        if (!final.valid) throw std::runtime_error("USER_FINAL_AFFINE_UNION_INVALID");
        const auto audit = zhangUserSharedContinuityAudit.observe(
            runtime + ":" + receiver + ":" + enum_to_string(system), double(state.time.bigTime),
            catalogue, zhangExactMultiply(final.rows, physical), final.values,
            reliability.confirmationEpochs, reliability.contradictionIsolation, columnOwners);
        if (!audit.valid) throw std::runtime_error("USER_INTEGER_CONTINUITY_AUDIT_INVALID:" + audit.reason);
        trace << "\nZHANG_USER_SHARED_INTEGER_CONTINUITY time=" << state.time.to_string(0)
              << " receiver=" << receiver << " provider=" << provider
              << " physical_domain_contradiction=" << audit.contradiction
              << " isolated_columns=" << audit.isolatedColumns
              << " isolated_owner_records=" << audit.isolatedOwners
              << " retained_rank=" << audit.retainedRank
              << " audit_evidence_rank=" << audit.auditEvidenceRank
              << " audit_evidence_columns=" << audit.auditEvidenceColumns
              << " confirmation_count=" << audit.confirmationCount
              << " confirmation_required=" << reliability.confirmationEpochs
              << " confirmation_eligible=" << audit.confirmationEligible
              << " status=" << audit.reason
              << " old_integer_truth_assumed=0 reference_transport=EXACT_PHYSICAL_REINDEX"
              << " audit_authority=SEARCH_CANDIDATE_METADATA_NOT_PUBLISHED_TRUTH"
              << " estimator_feedback=0 risk_probability_multiplication=0";
        string publicationStatus;
        if (reliability.safe || reliability.confirmationEpochs > 1)
        {
            ZhangExactMatrix auditedRows;
            if (audit.confirmationEligible)
                for (const auto& row : audit.output.rows)
                {
                    const auto member = zhangIntegerRowLatticeContains(physical, row);
                    if (!member.contained)
                        throw std::runtime_error("USER_AUDITED_ROW_OUTSIDE_CURRENT_SD_LATTICE");
                    auditedRows.push_back(member.combination);
                }
            final = zhangAuditUserIntegerUnion(auditedRows,
                audit.confirmationEligible ? audit.output.values : ZhangExactVector{}, ambient);
            if (!final.valid) throw std::runtime_error("USER_AUDITED_FINAL_DOMAIN_INVALID");
            if (!audit.confirmationEligible)
                publicationStatus = audit.isolatedColumns ? "INTEGER_CONTRADICTORY_COLUMNS_ISOLATED_FLOAT" :
                    (audit.retainedRank ? "RELIABILITY_CONFIRMATION_PENDING" : "FLOAT_ONLY");
            else if (audit.isolatedColumns) publicationStatus = "INTEGER_PARTIAL_CONTRADICTORY_COLUMNS_ISOLATED";
        }
        candidateImage=zhangProjectUserCandidateIntegerImage(
            zhangExactAffineIntegerQuotient(final.rows,final.values,ambient),candidateSelectors);
        if (!candidateImage.valid || candidateImage.fixedRank<carried)
            throw std::runtime_error("USER_FINAL_CANDIDATE_DOMAIN_INVALID:"+candidateImage.reason);
        if(final.rank)
        {
            MatrixXd integerRows; VectorXd integers;
            if(!zhangUserExactToDouble(final.rows,ambient,integerRows) ||
               !zhangUserExactToDouble(final.values,integers))
                throw std::runtime_error("USER_OUTPUT_EXACT_COEFFICIENT_OVERFLOW");
            outputRows.push_back(integerRows*transform); outputValues.push_back(integers);
            if(hold)
            {
                auto rows=ledger.rows; auto rhs=ledger.values;
                const auto physicalFixed=zhangExactMultiply(final.rows,physical);
                rows.insert(rows.end(),physicalFixed.begin(),physicalFixed.end());
                rhs.insert(rhs.end(),final.values.begin(),final.values.end());
                const auto united=zhangAuditUserIntegerUnion(rows,rhs,catalogue.size());
                if(!united.valid) throw std::runtime_error("USER_PHYSICAL_LEDGER_UNION_INVALID");
                ledger.rows=united.rows; ledger.values=united.values;
            }
        }
        if(hold) pending[key]=ledger;
        summary(receiver,system,candidateDimension,carried,candidateImage.fixedRank,retired,
            endpoints[reference].satellite.id(),!publicationStatus.empty() ? publicationStatus :
                (candidateImage.fixedRank?"INTEGER_CERTIFIED":"FLOAT_ONLY"));
        trace << "\nZHANG_USER_SHARED_LEDGER time=" << state.time.to_string(0)
              << " receiver=" << receiver << " persistent_physical_columns=" << catalogue.size()
              << " temporarily_masked_columns=" << masked
              << " full_domain_dimension=" << ambient << " full_domain_fixed_rank=" << final.rank
              << " candidate_fixed_rank=" << candidateImage.fixedRank
              << " candidate_free_image_rank=" << candidateImage.freeImageRank
              << " candidate_image_preserves_congruences=1"
              << " retirement_policy=CONSERVATIVE_EQUALITY_CONSEQUENCES"
              << " chart_reference_changes_do_not_reset=1";
    }
    int rowCount=0;
    for(const auto& rows:outputRows) rowCount+=rows.rows();
    if(rowCount)
    {
        MatrixXd rows(rowCount,state.x.size()); VectorXd integers(rowCount);
        int start=0;
        for(int i=0;i<outputRows.size();++i)
        { rows.middleRows(start,outputRows[i].rows())=outputRows[i];
          integers.segment(start,outputValues[i].size())=outputValues[i]; start+=outputRows[i].rows(); }
        const auto conditioned=zhangConditionPosteriorEffectiveIntegers(state.x,state.P,rows,integers);
        ZhangUserOutputValidation outputValidation;
        if(conditioned.valid)
            outputValidation=zhangValidateUserIntegerOutput(conditioned.mean,conditioned.covariance,rows,integers);
        trace << "\nZHANG_USER_SHARED_OUTPUT_NUMERICS time=" << state.time.to_string(0)
              << " status=" << (outputValidation.valid?"PASS":"FAIL")
              << " mean_residual=" << outputValidation.maximumMeanResidual
              << " mean_tolerance=" << outputValidation.maximumMeanTolerance
              << " covariance_residual=" << outputValidation.maximumCovarianceResidual
              << " covariance_tolerance=" << outputValidation.maximumCovarianceTolerance
              << " minimum_eigenvalue=" << outputValidation.minimumEigenvalue
              << " psd_tolerance=" << outputValidation.psdTolerance
              << " jitter=0 nis_evaluations=0";
        if(!conditioned.valid || !outputValidation.valid)
        {
            zhangTransactionalConditioningFailed=true;
            const string failure=conditioned.valid?outputValidation.reason:conditioned.failureReason;
            zhangTransactionalConditioningReason="SHARED_USER_OUTPUT_CONDITION:"+failure;
            trace << "\nZHANG_USER_SHARED_OUTPUT_TRANSACTION time=" << state.time.to_string(0)
                  << " status=ROLLED_BACK committed_rank=0 reason=" << failure;
            for (const auto& receipt:receipts) receipt(false);
            return;
        }
        state.dx=conditioned.mean-state.x;
        state.x=conditioned.mean; state.P=conditioned.covariance;
    }
    if(hold) for(auto& [key,ledger]:pending) zhangUserSharedLedgers[key]=std::move(ledger);
    else for(auto it=zhangUserSharedLedgers.begin();it!=zhangUserSharedLedgers.end();)
        if(std::get<0>(it->first)==runtime) it=zhangUserSharedLedgers.erase(it); else ++it;
    trace << "\nZHANG_USER_SHARED_OUTPUT_TRANSACTION time=" << state.time.to_string(0)
          << " status=COMMITTED committed_rank=" << rowCount
          << " held_ledger_committed=" << hold << " float_authority_unchanged=1";
    for (const auto& receipt:receipts) receipt(true);
}
