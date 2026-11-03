// Actual candidate KF transaction. QR is printed only as a numeric
// diagnostic; the independent literal witness (phase:1, ambiguity:1) alone is
// checked exactly against the ACTUAL production T callback.
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include "pea/zhangReference.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void check(bool value,const char* why)
{if(!value)throw std::runtime_error(why);}

static void dump(const KFState& state,const char* stage,const char* model)
{
    std::vector<ZhangProductFloatTargetRequest> targets;std::string reason;
    check(zhangReadFloatTargets(state.metaDataMap,targets,reason) && targets.size()==1,"probe target directory invalid");
    const auto& t=targets.front();
    std::cout<<"SHEAR_TARGET stage="<<stage<<" model="<<model<<" bound="<<t.numericBindingComplete
        <<" active="<<t.active<<" consumer="<<std::quoted(t.consumer)
        <<" reason="<<std::quoted(t.bindingReason)<<" owner="<<std::quoted(t.scope.runtimeId)
        <<" order_sha256="<<zhangCheckpointSha256(t.scope.stateOrderId)
        <<" generation="<<t.scope.coordinateGeneration<<" actual_factor_sequence="<<state.factorCommitSequence
        <<" symbolic_space="<<int(t.symbolicSpace)<<" offset="<<std::hexfloat<<t.numericOffset<<std::defaultfloat
        <<" chart="<<std::quoted(t.scope.chartId)<<" row=";
    for(const auto&[key,value]:t.numericRow)
        std::cout<<int(key.type)<<':'<<key.Sat.id()<<':'<<key.str<<':'<<key.num<<'='<<std::hexfloat<<value<<std::defaultfloat<<';';
    std::cout<<'\n';
}

static void probe(bool identityModels)
{
    const char* model=identityModels?"ALL_IDENTITY_MODELS":"PRODUCTION_HELPER_DEFAULT_MODELS";
    KFState state;
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=std::string("ACTUAL_FLOAT_TARGET_SHEAR_")+model;
    KFKey phase;phase.type=KF::PHASE_BIAS;phase.Sat=SatSys("G02");phase.num=int(E_ObsCode::L1C);
    KFKey ambiguity;ambiguity.type=KF::AMBIGUITY;ambiguity.Sat=phase.Sat;ambiguity.str="A";ambiguity.num=phase.num;
    KFKey clock;clock.type=KF::SAT_CLOCK;clock.Sat=phase.Sat;
    state.kfIndexMap={{KFState::oneKey,0},{phase,1},{ambiguity,2},{clock,3}};
    state.x.resize(4);state.x<<1,2,3,4;state.dx=VectorXd::Zero(4);
    state.P=MatrixXd::Zero(4,4);state.P(1,1)=4;state.P(2,2)=9;state.P(3,3)=16;
    state.P(1,2)=state.P(2,1)=1;state.P(1,3)=state.P(3,1)=2;state.P(2,3)=state.P(3,2)=3;
    if(identityModels)for(const auto&[key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
    ZhangProductFloatTargetRequest target;
    target.targetId="PENDING_PUBLICATION_DEPENDENCY|ACTUAL_SHEAR_LITERAL";
    target.consumer="PENDING_PUBLICATION_DEPENDENCY";target.lifecycleId="LITERAL_RESERVED_SOURCE_ROW";target.unit="m";
    target.symbolicSpace=ZhangProductFloatTargetSpace::SOURCE_KF_COORDINATE;
    target.scope={zhangProductFloatOwnerId(state),"ACTUAL_SOURCE_KF_CHART",zhangFloatCarrierStateOrder(state.kfIndexMap),state.factorCommitSequence};
    target.numericRow={{phase,1},{ambiguity,2}};target.numericBindingComplete=true;
    std::string reason;std::ostringstream trace;
    check(zhangRegisterGraphProductHandoffTarget(state,target,reason),"actual production registration failed");
    dump(state,"BEFORE",model);
    bool callback=false,literalWitness=false;
    state.exactStateTransformCallback=[&](const KFState& source,GTime,const std::map<KFKey,int>& sourceOrder,
        const std::map<KFKey,int>& destinationOrder,const SparseMatrix<double>& transform,const std::string&,
        const VectorXd&,const MatrixXd&,std::uint64_t before,std::uint64_t after)
    {
        callback=true;
        const std::map<KFKey,double> literalNew{{phase,1},{ambiguity,1}};
        ZhangFloatCarrierSparseRow pulled,old;
        for(const auto&[key,value]:target.numericRow)old[sourceOrder.at(key)]=*zhangFloatCarrierExactDouble(value);
        for(const auto&[key,value]:literalNew)
            for(int column=0;column<transform.cols();++column)
            {
                const double entry=transform.coeff(destinationOrder.at(key),column);
                if(entry!=0)pulled[column]+=*zhangFloatCarrierExactDouble(entry)* *zhangFloatCarrierExactDouble(value);
            }
        for(auto i=pulled.begin();i!=pulled.end();)if(i->second==0)i=pulled.erase(i);else ++i;
        literalWitness=pulled==old;
        VectorXd f=VectorXd::Zero(source.x.size());f(sourceOrder.at(phase))=1;f(sourceOrder.at(ambiguity))=2;
        SparseMatrix<double> transpose=transform.transpose();transpose.makeCompressed();
        Eigen::SparseQR<SparseMatrix<double>,Eigen::COLAMDOrdering<int>> qr;qr.compute(transpose);
        check(qr.info()==Eigen::Success,"independent QR observation failed");
        const VectorXd candidate=qr.solve(f);
        std::cout<<"SHEAR_ACTUAL_CALLBACK model="<<model<<" before_sequence="<<before<<" after_sequence="<<after
            <<" literal_w_phase_1_ambiguity_1_exact_pullback="<<literalWitness
            <<" qr_phase="<<std::hexfloat<<candidate(destinationOrder.at(phase))
            <<" qr_ambiguity="<<candidate(destinationOrder.at(ambiguity))<<std::defaultfloat
            <<" qr_numeric_residual="<<(transpose*candidate-f).norm()<<'\n';
        return true;
    };
    const auto initialSequence=state.factorCommitSequence;
    check(state.applyStateTransform(trace,{{KFState::oneKey,{{KFState::oneKey,1}}},
        {phase,{{phase,1},{ambiguity,1}}},{ambiguity,{{ambiguity,1}}},{clock,{{clock,1}}}},
        "ACTUAL_CANDIDATE_KF_SHEAR_RECOVERY"),"actual production KF shear did not commit");
    check(callback && literalWitness && state.factorCommitSequence==initialSequence+1,"actual callback/exact witness/sequence absent");
    dump(state,"AFTER",model);
    std::vector<ZhangProductFloatTargetRequest> after;
    check(zhangReadFloatTargets(state.metaDataMap,after,reason),"after directory read failed");
    const auto& recovered=after.front();
    check(recovered.numericBindingComplete && recovered.numericRow==std::map<KFKey,double>{{phase,1},{ambiguity,1}} &&
        recovered.scope.runtimeId==target.scope.runtimeId &&
        recovered.scope.stateOrderId==zhangFloatCarrierStateOrder(state.kfIndexMap) &&
        recovered.scope.coordinateGeneration==state.factorCommitSequence &&
        recovered.symbolicSpace==target.symbolicSpace && recovered.numericOffset==target.numericOffset,
        "actual exact fallback did not preserve complete source meaning and current binding");
    const double literalOldMean=2+2*3;
    const double literalNewMean=state.x(state.kfIndexMap.at(phase))+state.x(state.kfIndexMap.at(ambiguity));
    check(literalOldMean==literalNewMean,"literal exact phase+ambiguity witness changed actual target mean");
    VectorXd w=VectorXd::Zero(state.x.size());w(state.kfIndexMap.at(phase))=1;w(state.kfIndexMap.at(ambiguity))=1;
    check((w.transpose()*state.P*w)(0,0)==44 && (w.transpose()*state.P)(state.kfIndexMap.at(clock))==8 &&
        state.P(state.kfIndexMap.at(clock),state.kfIndexMap.at(clock))==16,
        "literal independent full phase/clock covariance oracle differs");
    check(zhangProductFloatTargetMetadataValid(state,reason),"recovered binding invalid on actual posterior");
    std::cout<<"SHEAR_PROBE_OBSERVATION model="<<model<<" actual_after_bound="<<after.front().numericBindingComplete
        <<" independent_literal_exact_witness=1 identical_literal_target_mean=1 new_mean="<<literalNewMean
        <<" literal_variance=44 literal_phase_clock_cross=8 original_objects_unchanged=1 integer_authority_granted=0 scientific_acceptance=0\n";
}

int main()
{
    setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
    probe(false);probe(true);
    return 0;
}
