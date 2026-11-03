// Included after GONSIS core/adapter. Conditional integer evidence is owned
// separately from FLOAT. Process restart naturally clears this registry.
#include "common/zhangUserConditionalHold.hpp"
static map<tuple<string,string,E_Sys>,ZhangUserConditionalHold> zhangUserConditionalLedgers;

static void resolveConditionalHoldUserAmbiguities(Trace& trace,KFState& state)
{
    if (!zhangUserConditionalHoldRequested() || !zhangUserGonsisSearchRequested() || !zhangRatioOnly() ||
        acsConfig.ambrOpts.mode!=E_ARmode::LAMBDA_ALT)
        throw std::runtime_error("CONDITIONAL_HOLD_REQUIRES_GONSIS_LAMBDA_RATIO_ONLY");
    const auto* owner=state.factorCallbackOwnerOverride;
    if (!owner) throw std::runtime_error("CONDITIONAL_HOLD_REQUIRES_FLOAT_OWNER");
    const string runtime=zhangAmbresRuntimeId(*owner);
    const string branchRuntime=zhangAmbresRuntimeId(state);
    trace << "\nZHANG_USER_GONSIS_POLICY time=" << state.time.to_string(0)
          << " profile=GONSIS_UDUC feedback=CONDITIONAL_HOLD ratio_threshold=" << acsConfig.ambrOpts.ratioThres
          << " ratio_comparison=STRICT_GREATER integer_history=1 float_feedback=0 nis_evaluations=0"
          << " output_conditioning=ONE_JOINT_FULL_COVARIANCE_TRANSACTION";
    const auto found=zhangUserGonsisObservationSnapshots.find(runtime);
    if (found==zhangUserGonsisObservationSnapshots.end() || !found->second.valid ||
        found->second.time!=state.time || found->second.commitSequence!=state.factorCommitSequence ||
        found->second.columns!=state.kfIndexMap || !state.lastPppTransactionCommitted || !state.lastPppQcConverged ||
        found->second.floatMean.size()!=state.x.size() ||
        (found->second.floatMean.array()!=state.x.array()).any() ||
        (found->second.floatCovariance.array()!=state.P.array()).any())
        throw std::runtime_error("CONDITIONAL_HOLD_REQUIRES_SAME_EPOCH_UNMODIFIED_FLOAT");
    using Key=tuple<string,string,E_Sys>;
    struct Frame
    {
        Key key;
        vector<ZhangUserSharedEndpoint> endpoints;
        map<string,string> rawIdentities;
        map<string,ZhangUserGonsisReferenceHint> searchHints;
        vector<ZhangUserHoldNode> active;
        ZhangUserConditionalHold proposed;
        ZhangUserHoldProjection carried,combined;
        int eligible=0,retired=0,searchCalls=0,conflicts=0;
        double lambda1=0,lambda2=0;
    };
    map<pair<string,E_Sys>,Receiver*> groups;
    for (const auto& [key,index]:state.kfIndexMap)
        if (key.type==KF::AMBIGUITY && key.rec_ptr && acsConfig.solve_amb_for[key.Sat.sys])
            groups[{key.str,key.Sat.sys}]=key.rec_ptr;
    vector<Frame> frames;
    auto node=[](const ZhangUserSharedEndpoint& e)
    {return ZhangUserHoldNode{e.satellite.id(),e.identity,e.componentDomain};};
    for (const auto& [group,receiver]:groups)
    {
        if (group.second!=E_Sys::GPS) throw std::runtime_error("CONDITIONAL_HOLD_VALIDATED_GPS_ONLY");
        Frame f; f.key={runtime,group.first,group.second};
        f.endpoints=buildZhangUserSharedEndpoints(trace,state,receiver,group.second,state.time,true);
        for (const auto& e:f.endpoints) f.rawIdentities[e.satellite.id()]=e.identity;
        f.endpoints.erase(std::remove_if(f.endpoints.begin(),f.endpoints.end(),[&](auto& e)
        {
            string authorizations;
            for (const auto code:{e.firstCode,e.secondCode})
            {
                ZhangInternalProduct product;
                if (!queryZhangInternalProduct(state.time,e.satellite,code,product) ||
                    !product.pppar_usable || !product.dual_frequency_ar_valid ||
                    !zhangFormalPppArProductSolution(product.solution)) return true;
                authorizations+=":"+std::to_string(zhangUserArcGeneration(state.metaDataMap,
                    zhangUserIntegerAuthorizationGenerationKey(group.first,e.satellite.id(),static_cast<int>(code))));
            }
            // OFFSET describes the current expression of a canonical integer,
            // not its identity. A certified affine/KF-reference transport must
            // not retire every integer just because those offsets changed.
            // The existing FLOAT lifecycle increments IA on a missing/invalid
            // integer map, even if a REAL FLOAT transport is still authorised.
            e.identity=zhangUserArEncodedField(e.physicalIdentity)+zhangUserArEncodedField(e.stateIdentity)+
                zhangUserArEncodedField(e.componentDomain)+":IA"+authorizations;
            return false;
        }),f.endpoints.end());
        vector<ZhangUserHoldNode> observed;
        for (const auto& e:f.endpoints) observed.push_back(node(e));
        auto& ledger=zhangUserConditionalLedgers[f.key];
        f.retired=ledger.retireReplaced(observed); // invalidation is NOT rolled back by an output failure
        vector<int> order;
        for (int i=0;i<f.endpoints.size();++i) if (f.endpoints[i].eligible) order.push_back(i);
        std::stable_sort(order.begin(),order.end(),[&](int a,int b){return f.endpoints[a].elevation>f.endpoints[b].elevation;});
        for (const int i:order) f.active.push_back(node(f.endpoints[i]));
        f.carried=ledger.project(f.active); f.proposed=ledger;
        if (!f.endpoints.empty())
        {
            f.lambda1=zhangE27Wavelength(group.second,f.endpoints.front().firstCode);
            f.lambda2=zhangE27Wavelength(group.second,f.endpoints.front().secondCode);
        }
        map<string,vector<int>> domains;
        for (int i=0;i<f.endpoints.size();++i) domains[f.endpoints[i].componentDomain].push_back(i);
        auto retireInactive=[&](auto& records)
        {
            for (auto it=records.begin();it!=records.end();)
                if (std::get<0>(it->first)==branchRuntime && std::get<1>(it->first)==group.first &&
                    std::get<2>(it->first)==group.second && !domains.count(std::get<3>(it->first))) it=records.erase(it);
                else ++it;
        };
        retireInactive(zhangUserGonsisReferences); retireInactive(zhangUserGonsisFixedLocks);
        const auto start=zhangValidateUserGonsisPosition(state.P,state.P,zhangUserGonsisPositionIndices(state,group.first));
        for (const auto& [domain,indices]:domains)
        {
            vector<ZhangUserHoldNode> eligible;
            for (const auto& n:f.active) if (n.domain==domain) eligible.push_back(n);
            const int dimension=2*std::max(0,int(eligible.size())-1); f.eligible+=dimension;
            const int carried=2*ledger.project(eligible).edges.size();
            trace << "\nZHANG_USER_CONDITIONAL_HOLD_DOMAIN time=" << state.time.to_string(0)
                  << " receiver=" << group.first << " component_domain=" << domain
                  << " eligible_dimension=" << dimension << " carried_rank=" << carried
                  << " reference=" << (eligible.empty()?"NONE":eligible.front().name)
                  << " continuity_test=EXACT_PHYSICAL_STATE_PRODUCT_CHART_IDENTITIES"
                  << " missing_observation_action=SUSPEND known_identity_change_action=RETIRE_AFFECTED_NODE"
                  << " full_constraint_noise_transport=1";
            if (!dimension || carried==dimension || !std::isfinite(start.floatAverageSigma) || start.floatAverageSigma>0.5) continue;
            MatrixXd transform(2*indices.size(),state.x.size());
            vector<ZhangUserGonsisSatellite> satellites;
            const auto referenceKey=std::make_tuple(branchRuntime,group.first,group.second,domain);
            auto& fixedLocks=zhangUserGonsisFixedLocks[referenceKey];
            for (int i=0;i<indices.size();++i)
            {
                const auto& e=f.endpoints[indices[i]];
                transform.row(2*i)=e.firstRow.transpose(); transform.row(2*i+1)=e.secondRow.transpose();
                const double variance=f.lambda1*f.lambda1*(e.firstRow.transpose()*state.P*e.firstRow)(0,0)+
                    f.lambda2*f.lambda2*(e.secondRow.transpose()*state.P*e.secondRow)(0,0);
                const auto& rawIdentity=f.rawIdentities.at(e.satellite.id());
                auto& lock=fixedLocks[e.satellite.id()];
                if (lock.identity!=rawIdentity || (lock.lastIncrement!=GTime::noTime() && state.time<lock.lastIncrement))
                    lock={rawIdentity,0,GTime::noTime()};
                satellites.push_back({e.satellite.id(),rawIdentity,e.elevation,lock.count,variance,e.eligible});
            }
            ZhangUserGonsisOptions options;
            options.ratioThreshold=acsConfig.ambrOpts.ratioThres;
            options.minimumWideLanePairs=zhangUserArReliabilityPolicy().minimumParRank;
            options.maximumTrials=zhangUserArReliabilityPolicy().maximumParTrials;
            options.wavelength1Metres=f.lambda1; options.wavelength2Metres=f.lambda2;
            auto solver=[&](const VectorXd& mean,const MatrixXd& covariance)
            {
                ZhangUserGonsisIlsSolution out; GinAR_mtx matrix; matrix.aflt=mean; matrix.Paflt=covariance;
                GinAR_opt opt; opt.mode=E_ARmode::LAMBDA_ALT; opt.ratthr=options.ratioThreshold;
                opt.sucthr=acsConfig.ambrOpts.succsThres; opt.nset=2; opt.nitr=acsConfig.ambrOpts.AR_max_itr;
                opt.min_lambda_fix_count=mean.size(); opt.max_lambda_fix_count=mean.size();
                opt.force_bootstrap_success=false; opt.maximum_lambda_conditional_sigma=0; opt.lambda_candidate_nis_alpha=0;
                GNSS_AR(trace,matrix,opt);
                if (matrix.searchDiagnostic.localNisExecuted || matrix.lambda_candidate_nis_valid)
                    throw std::runtime_error("CONDITIONAL_HOLD_UNEXPECTED_NIS");
                out.reason=matrix.searchDiagnostic.reason;
                out.bestSquaredCost=matrix.searchDiagnostic.bestSquaredDistance;
                out.secondSquaredCost=matrix.searchDiagnostic.secondSquaredDistance;
                out.complete=matrix.searchDiagnostic.ilsComplete && matrix.searchDiagnostic.ratioExecuted &&
                    matrix.Ztrs.rows()==mean.size() && matrix.Ztrs.cols()==mean.size() && matrix.zfix.size()==mean.size();
                if (out.complete)
                {
                    ZhangExactMatrix h; ZhangExactVector b;
                    if (!zhangUserSearchExact(matrix.Ztrs,matrix.zfix,h,b))
                        throw std::runtime_error("CONDITIONAL_HOLD_NONINTEGER_SEARCH_CHART");
                    const auto exact=zhangExactAffineIntegerQuotient(h,b,mean.size());
                    if (!exact.valid || exact.quotientRank || !zhangUserExactToDouble(exact.particularSolution,out.integers))
                        throw std::runtime_error("CONDITIONAL_HOLD_INTEGER_RECOVERY_FAILED");
                }
                return out;
            };
            ZhangUserGonsisReferenceHint referenceHint;
            const auto oldReference=zhangUserGonsisReferences.find(referenceKey);
            if (oldReference!=zhangUserGonsisReferences.end())
            {
                const double age=(state.time-oldReference->second.time).to_double();
                if (age>0 && age<=std::max(45.0,1.5*acsConfig.epoch_interval)) referenceHint=oldReference->second.hint;
            }
            const auto result=zhangResolveUserGonsis(transform*state.x,transform*state.P*transform.transpose(),satellites,solver,options,referenceHint);
            f.searchCalls+=result.ilsCalls;
            for (const auto& t:result.trials)
                trace << "\nZHANG_USER_CONDITIONAL_HOLD_SEARCH time=" << state.time.to_string(0)
                      << " receiver=" << group.first << " component_domain=" << domain
                      << " ratio=" << std::setprecision(17) << t.ratio << " threshold=" << options.ratioThreshold
                      << " reason=" << t.reason << " nis_evaluations=0";
            if (result.accepted)
            {
                vector<ZhangUserHoldEdge> candidates;
                for (int i=1;i<result.acceptedSatelliteIndices.size();++i)
                    candidates.push_back({node(f.endpoints[indices[result.acceptedSatelliteIndices[i]]]),
                        node(f.endpoints[indices[result.referenceIndex]]),
                        zhangUserHoldInteger(result.wideLaneIntegers(i-1)),zhangUserHoldInteger(result.l1Integers(i-1))});
                int added=0;
                if (!f.proposed.merge(candidates,options.wideLaneConstraintVarianceCycles2,options.ifConstraintVarianceMetres2,added))
                    ++f.conflicts;
                else f.searchHints[domain]={result.referenceSatellite,result.referenceIdentity,true};
            }
        }
        f.combined=f.proposed.project(f.active); frames.push_back(std::move(f));
    }
    auto condition=[&](bool historicalOnly,ZhangUserGonsisOutput& out,string& reason)
    {
        int count=0;
        for (const auto& f:frames) count+=2*(historicalOnly?f.carried:f.combined).edges.size();
        if (!count) {reason="NO_APPLICABLE_INTEGER_RELATIONS";return false;}
        MatrixXd h=MatrixXd::Zero(count,state.x.size()),r=MatrixXd::Zero(count,count); VectorXd z(count);
        int offset=0;
        for (const auto& f:frames)
        {
            const auto& p=historicalOnly?f.carried:f.combined;
            const double nl=f.lambda1*f.lambda2/(f.lambda1+f.lambda2);
            const double wc=f.lambda1*f.lambda1*f.lambda2/(f.lambda2*f.lambda2-f.lambda1*f.lambda1);
            for (const auto& edge:p.edges)
            {
                auto endpoint=[&](const auto& n)->const ZhangUserSharedEndpoint&
                {for (const auto& e:f.endpoints) if (node(e)==n) return e;throw std::runtime_error("HOLD_MISSING_CURRENT_ENDPOINT");};
                const auto& a=endpoint(edge.a); const auto& b=endpoint(edge.b);
                const VectorXd l1=a.firstRow-b.firstRow,wl=l1-a.secondRow+b.secondRow;
                h.row(offset)=wl.transpose(); h.row(offset+1)=(nl*l1+wc*wl).transpose();
                if (std::abs((long double)edge.wl)>9007199254740991.0L || std::abs((long double)edge.l1)>9007199254740991.0L)
                    throw std::runtime_error("HOLD_PROJECTED_INTEGER_NOT_EXACT_DOUBLE");
                z(offset)=edge.wl; z(offset+1)=nl*edge.l1+wc*edge.wl; offset+=2;
            }
            const int n=p.noise.rows(); r.block(offset-n,offset-n,n,n)=p.noise;
        }
        out=zhangConditionUserGonsisOutput(state.x,state.P,h,z,r.diagonal(),&r);
        if (!out.valid) {reason=out.reason;return false;}
        for (const auto& f:frames)
        {
            if ((historicalOnly?f.carried:f.combined).edges.empty()) continue;
            const auto gate=zhangValidateUserGonsisPosition(state.P,out.covariance,zhangUserGonsisPositionIndices(state,std::get<1>(f.key)));
            if (!gate.valid) {reason=gate.reason;return false;}
        }
        const auto& s=found->second;
        const auto residual=zhangValidateUserGonsisResiduals(state.x,out.mean,s.design,s.residuals,s.variances,s.kinds);
        reason=residual.reason; return residual.valid;
    };
    ZhangUserGonsisOutput output; string reason;
    bool historicalOnly=false;
    bool accepted=condition(false,output,reason);
    if (!accepted)
    {
        const string unionReason=reason; historicalOnly=true; accepted=condition(true,output,reason);
        trace << "\nZHANG_USER_CONDITIONAL_HOLD_RETRY time=" << state.time.to_string(0)
              << " combined_failure=" << unionReason << " held_only_accepted=" << accepted
              << " reason=" << reason << " same_unconditioned_FLOAT=1";
    }
    if (accepted) {state.dx=output.mean-state.x;state.x=output.mean;state.P=output.covariance;}
    for (auto& f:frames)
    {
        const int carried=2*f.carried.edges.size();
        const int finalRank=accepted?2*(historicalOnly?f.carried:f.combined).edges.size():0;
        if (accepted && !historicalOnly) zhangUserConditionalLedgers[f.key]=std::move(f.proposed);
        // A failed held-only numerical/residual gate cannot keep asserting the
        // suspect integer set. Release integer evidence, NOT the physical FLOAT.
        if (!accepted && carried) zhangUserConditionalLedgers[f.key]=ZhangUserConditionalHold{};
        const auto& projection=historicalOnly?f.carried:f.combined;
        // Use the same source-style reference retention/promotion metadata as
        // epoch-only GONSIS. Holding must not silently replace the supplement
        // search's reference policy by an always-highest-elevation policy.
        map<string,vector<ZhangUserHoldNode>> publishedNodes;
        if (accepted)
            for (const auto& e:projection.edges) for (const auto& n:{e.a,e.b})
            {
                auto& nodes=publishedNodes[n.domain];
                if (std::find(nodes.begin(),nodes.end(),n)==nodes.end()) nodes.push_back(n);
            }
        set<string> currentDomains;
        for (const auto& e:f.endpoints) currentDomains.insert(e.componentDomain);
        for (const auto& domain:currentDomains)
        {
            const auto referenceKey=std::make_tuple(branchRuntime,std::get<1>(f.key),std::get<2>(f.key),domain);
            const auto nodes=publishedNodes.find(domain);
            if (nodes==publishedNodes.end()) {zhangUserGonsisReferences.erase(referenceKey);continue;}
            auto& locks=zhangUserGonsisFixedLocks[referenceKey];
            vector<ZhangUserGonsisSatellite> fixed; vector<int> active;
            for (const auto& n:nodes->second)
                for (const auto& e:f.endpoints) if (node(e)==n)
                {
                    const auto& raw=f.rawIdentities.at(n.name); auto& lock=locks[n.name];
                    if (lock.identity!=raw) lock={raw,0,GTime::noTime()};
                    active.push_back(fixed.size());fixed.push_back({n.name,raw,e.elevation,lock.count,0,true});
                }
            ZhangUserGonsisReferenceHint hint;
            const auto old=zhangUserGonsisReferences.find(referenceKey);
            if (old!=zhangUserGonsisReferences.end())
            {
                const double age=(state.time-old->second.time).to_double();
                if (age>0 && age<=std::max(45.0,1.5*acsConfig.epoch_interval)) hint=old->second.hint;
            }
            if (!historicalOnly && f.searchHints.count(domain)) hint=f.searchHints.at(domain);
            else
            {
                const int selected=zhangUserGonsisDetail::reference(active,fixed,hint);
                hint={fixed[selected].satellite,fixed[selected].identity,true};
            }
            zhangUserGonsisReferences[referenceKey]={state.time,hint};
            for (const auto& satellite:fixed)
            {
                auto& lock=locks[satellite.satellite];
                if (lock.lastIncrement!=state.time)
                {lock.count=std::min(1000000,lock.count+1);lock.lastIncrement=state.time;}
            }
        }
        trace << "\nZHANG_USER_GONSIS_AR_SUMMARY time=" << state.time.to_string(0)
              << " receiver=" << std::get<1>(f.key) << " provider=internal system=GPS"
              << " eligible_dimension=" << f.eligible << " carried_rank=" << (accepted?carried:0)
              << " newly_fixed_rank=" << (accepted?std::max(0,finalRank-carried):0) << " final_rank=" << finalRank
              << " wl_rank=" << finalRank/2 << " l1_rank=" << finalRank/2
              << " reference=" << (projection.edges.empty()?"NONE":projection.edges.front().b.name)
              << " status=" << (accepted?(historicalOnly?"HELD_ONLY_ACCEPTED":"CONDITIONAL_UNION_ACCEPTED"):reason)
              << " feedback=CONDITIONAL_HOLD float_feedback=0 nis_evaluations=0"
              << " full_constraint_noise_transport=1 integer_truth_certified=0";
        trace << "\nZHANG_USER_CONDITIONAL_HOLD_AUDIT time=" << state.time.to_string(0)
              << " receiver=" << std::get<1>(f.key) << " retired_nodes=" << f.retired
              << " conflicting_current_sets=" << f.conflicts << " actual_ils_calls=" << f.searchCalls
              << " stored_rank=" << 2*zhangUserConditionalLedgers[f.key].edges.size()
              << " held_gate_failure_release=" << (!accepted && carried>0)
              << " output_gate=" << reason << " observation_update=0 float_authority_unchanged=1";
    }
}
