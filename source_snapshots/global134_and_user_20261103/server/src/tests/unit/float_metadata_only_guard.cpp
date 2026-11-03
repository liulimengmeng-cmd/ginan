#include "common/zhangCheckpoint.hpp"
#include "common/zhangFloatMetadataCommit.hpp"
#include <functional>
#include <iostream>
#include <stdexcept>

static void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main()
{
    KFState state;KFKey a;a.type=KF::PHASE_BIAS;a.Sat=SatSys("G01");
    state.time.bigTime=1000;state.kfIndexMap={{KFState::oneKey,0},{a,1}};
    state.x.resize(2);state.x<<1,2;state.dx=VectorXd::Zero(2);
    state.P=MatrixXd::Zero(2,2);state.P(1,1)=2;
    state.prefitRatios=state.postfitRatios=VectorXd::Zero(2);
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="METADATA_ONLY_GUARD";
    state.factorCommitSequence=17;
    const auto digest=zhangKfNonMetadataMutationDigest(state);
    check(digest.size()==64,"nonmetadata guard could not hash a valid source");
    auto metadataOnly=state;metadataOnly.metaDataMap["UNRELATED_METADATA"]="guard excludes only this entire field";
    check(zhangKfNonMetadataMutationDigest(metadataOnly)==digest,"metadata guard changed when only metaDataMap changed");
    const std::vector<std::function<void(KFState&)>> mutations={
        [](auto&s){s.time.bigTime+=1;},[](auto&s){s.x(1)+=.25;},
        [](auto&s){s.P(0,1)=-0.;},[](auto&s){s.dx(1)+=.25;},
        [](auto&s){s.prefitRatios(1)+=.25;},[](auto&s){s.postfitRatios(1)+=.25;},
        [a](auto&s){s.kfIndexMap.at(a)=0;},[a](auto&s){s.stateTransitionMap[a][a][0]=2;},
        [a](auto&s){s.gaussMarkovTauMap[a]=2;},[a](auto&s){s.gaussMarkovMuMap[a]=2;},
        [a](auto&s){s.procNoiseMap[a]=2;},[a](auto&s){s.initNoiseMap[a]=2;},
        [a](auto&s){s.sigmaMaxMap[a]=2;},[a](auto&s){s.outageLimitMap[a]=2;},
        [a](auto&s){s.exponentialNoiseMap[a]={2,3};},[a](auto&s){s.pseudoStateMap[a][a]=2;},
        [a](auto&s){s.pseudoParentMap[a]=KFState::oneKey;},[a](auto&s){s.errorCountMap[a]=2;},
        [](auto&s){s.filterChunkMap["NEW"].begX=2;},[](auto&s){s.lsqRequired=!s.lsqRequired;},
        [](auto&s){s.sigmaPass=!s.sigmaPass;},[](auto&s){s.chiQCPass=!s.chiQCPass;},
        [](auto&s){s.chi2=2;},[](auto&s){s.dof=2;},[](auto&s){s.chi2PerDof=2;},
        [](auto&s){s.qc=2;},[](auto&s){s.id+="changed";},[](auto&s){s.rts_basename+="changed";},
        [](auto&s){s.output_residuals=!s.output_residuals;},
        [](auto&s){s.outputMongoMeasurements=!s.outputMongoMeasurements;},
        [](auto&s){s.statisticsMap["changed"]=2;},[](auto&s){s.statisticsMapSum["changed"]=2;},
        [](auto&s){++s.factorCommitSequence;}
    };
    check(mutations.size()==33,"test must cover32 nonmetadata fields and the true runtime counter");
    for(const auto& mutate:mutations)
    {auto changed=state;mutate(changed);check(zhangKfNonMetadataMutationDigest(changed)!=digest,"a nonmetadata field mutation escaped its exact raw-byte guard");}
    const auto before=state.metaDataMap;
    auto after=before;const std::string target="ZHANG_FLOAT_TARGET_V2/explicit-target";
    after[target]="verified draft bytes";
    std::string reason;
    check(zhangValidateFloatPublicationMetadataDelta(before,after,{target},{{target,after.at(target)}},reason),
        "explicit exact delta rejected");
    check(!zhangValidateFloatPublicationMetadataDelta(before,after,{},{{target,after.at(target)}},reason),"missing key authority accepted");
    check(!zhangValidateFloatPublicationMetadataDelta(before,after,{target},{{target,"wrong bytes"}},reason),"wrong expected value accepted");
    auto erased=after;erased.erase(ZHANG_CHECKPOINT_RUNTIME_ID_METADATA);
    check(!zhangValidateFloatPublicationMetadataDelta(before,erased,{target},{{target,after.at(target)}},reason),"existing identity deletion accepted");
    auto forged=before;forged[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="FORGED";
    check(!zhangValidateFloatPublicationMetadataDelta(before,forged,{ZHANG_CHECKPOINT_RUNTIME_ID_METADATA},
        {{ZHANG_CHECKPOINT_RUNTIME_ID_METADATA,"FORGED"}},reason),"even an explicit key list permitted runtime-owner replacement");
    std::cout<<"PASS float_metadata_only_guard core_nonmetadata_fields=32 runtime_counter=1 "
        "full_IEEE_including_signed_zero=1 explicit_keys_values_only=1 core_identity_never_authorised=1\n";
}
