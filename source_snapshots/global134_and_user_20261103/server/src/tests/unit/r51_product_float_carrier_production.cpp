// Link the actual Reference/PppAr/algebra production objects. Literal small
// Gaussian moments below are independent of the carrier planner/resolver.
#include "common/acsConfig.hpp"
#include "common/constants.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangIntegerProductGainFrontier.hpp"
#include "common/zhangR51PhysicalImage.hpp"
#include "common/zhangWholeFloatSupport.hpp"
#include "pea/zhangReference.hpp"
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool ok,const char* why)
{if(!ok)throw std::runtime_error(why);}

struct Fixture
{
    KFState state;
    KFKey a,b,clock;
    ZhangGraphEdge ea{"A",SatSys("G01")},eb{"B",SatSys("G02")};
    Fixture(const std::string& owner)
    {
        state.time.bigTime=1405216800;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=owner;
        a=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,ea,0);
        b=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,eb,0);
        clock.type=KF::SAT_CLOCK;clock.Sat=SatSys("G01");
        state.kfIndexMap={{KFState::oneKey,0},{a,1},{b,2},{clock,3}};
        state.x.resize(4);state.x<<1,2,5,10;state.dx=VectorXd::Zero(4);
        state.P.resize(4,4);state.P<<0,0,0,0, 0,4,1,2, 0,1,9,3, 0,2,3,16;
        for(const auto& [key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
        require(zhangHandoffFloatRegisterRecord(state,a,ea,int(E_ObsCode::L1C),0),"first RAW registry setup failed");
        require(zhangHandoffFloatRegisterRecord(state,b,eb,int(E_ObsCode::L1C),0),"second RAW registry setup failed");
    }
    ZhangProductFloatTargetRequest reservation(bool complete=true) const
    {
        ZhangProductFloatTargetRequest target;
        target.targetId="LITERAL_PENDING_GAUSSIAN_RESERVATION";
        target.consumer="PENDING_PUBLICATION_DEPENDENCY";
        target.lifecycleId="LITERAL_CURRENT_FLOAT_PHASE";target.unit="m";
        target.scope={zhangProductFloatOwnerId(state),"LITERAL_SOURCE_CHART",
            zhangFloatCarrierStateOrder(state.kfIndexMap),state.factorCommitSequence};
        target.numericRow={{a,1},{b,1}};target.numericBindingComplete=complete;
        target.bindingReason=complete?"ACTUAL_LITERAL_JACOBIAN":"MISSING_OLD_ANCHOR";
        return target;
    }
};

static void rawEqual(const KFState& actual,const KFState& before)
{
    require(actual.kfIndexMap==before.kfIndexMap && actual.metaDataMap==before.metaDataMap &&
        actual.stateTransitionMap==before.stateTransitionMap &&
        actual.factorCommitSequence==before.factorCommitSequence,"raw fallback changed lifecycle fields");
    require(actual.x.size()==before.x.size() && actual.P.size()==before.P.size() &&
        !std::memcmp(actual.x.data(),before.x.data(),actual.x.size()*sizeof(double)) &&
        !std::memcmp(actual.dx.data(),before.dx.data(),actual.dx.size()*sizeof(double)) &&
        !std::memcmp(actual.P.data(),before.P.data(),actual.P.size()*sizeof(double)),
        "raw fallback changed full moment bytes");
}

static void activeDependenciesAreNotCacheEvictions()
{
    for(const bool complete:{true,false})
    {
        Fixture f(complete?"V2_BOUND_RAW_LIFETIME":"V2_UNKNOWN_RAW_LIFETIME");
        std::string reason;std::ostringstream trace;
        auto target=f.reservation(complete);
        require(zhangRegisterGraphProductHandoffTarget(f.state,target,reason),"whole RAW reservation registration failed");
        f.state.time=f.state.time+1000; // older than cache TTL, no legacy DEP tags
        const auto before=f.state;
        require(retireZhangGraphProductEventAuxiliaries(trace,f.state,"V2_ACTIVE_LIFETIME_CONTRACT"),
            "active V2 consumer was treated as an expired inactive RAW cache");
        rawEqual(f.state,before);
        require(f.state.kfIndexMap.contains(f.a) && f.state.kfIndexMap.contains(f.b),"active RAW information was evicted");
        if(!complete)
        {
            require(zhangCompressGraphProductFloatTargets(trace,f.state,reason),"legitimate missing anchor failed FLOAT");
            rawEqual(f.state,before);
        }
    }
}

static void realCompressionAndWholeResolution()
{
    Fixture f("V2_CARRIER_REAL_PRODUCTION");
    std::string reason;std::ostringstream trace;
    auto target=f.reservation();
    require(zhangRegisterGraphProductHandoffTarget(f.state,target,reason),"whole phase registration failed");
    const auto currentDomains=zhangCurrentEpochEvidenceRegistry().size();
    require(zhangCompressGraphProductFloatTargets(trace,f.state,reason),"production compression transaction failed");
    require(reason=="GLOBAL_WHOLE_FLOAT_TARGET_SPACE_COMPRESSED","real target-space compression never committed");
    require(!f.state.kfIndexMap.contains(f.a) && !f.state.kfIndexMap.contains(f.b),"raw auxiliaries were merely retained/renamed");
    std::vector<KFKey> carriers;
    for(const auto& [key,index]:f.state.kfIndexMap)if(zhangIsProductFloatCarrierKey(key))carriers.push_back(key);
    require(carriers.size()==1 && f.state.x.size()==3,"whole target rank did not reduce two RAW directions to one");
    const int carrier=f.state.kfIndexMap.at(carriers.front()),clock=f.state.kfIndexMap.at(f.clock);
    require(f.state.x(carrier)==7 && f.state.P(carrier,carrier)==15 && f.state.P(carrier,clock)==5 &&
        f.state.P(clock,clock)==16,"literal full phase/clock Gaussian oracle differs");
    require(zhangCurrentEpochEvidenceRegistry().size()==currentDomains,"FLOAT compression injected CURRENT integer RHS");
    require(zhangGraphProductHandoffStateValid(f.state,reason),"typed carrier/whole row registry invalid after commit");
    KFKey individual;
    require(!zhangHandoffFloatFind(f.state,{f.ea,int(E_ObsCode::L1C),0},individual),"carrier was forged as a RAW physical identity");
    ZhangProductFloatTargetRequest physical;
    physical.targetId="WHOLE_PHASE_SUM_QUERY";physical.consumer="FLOAT_QUERY";physical.unit="m";
    physical.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
    physical.physicalCoefficients[{f.ea,E_ObsCode::L1C,0}]=1;
    physical.physicalCoefficients[{f.eb,E_ObsCode::L1C,0}]=1;
    std::map<KFKey,double> row;double offset=0;
    require(zhangResolveGraphProductHandoffTarget(f.state,physical,row,offset,reason),
        "known whole phase sum was rejected because its constituent RAW arcs are gone");
    require(row==std::map<KFKey,double>{{carriers.front(),1}} && offset==0,"whole sum did not resolve to the actual carrier");
    physical.physicalCoefficients.erase({f.eb,E_ObsCode::L1C,0});
    require(!zhangResolveGraphProductHandoffTarget(f.state,physical,row,offset,reason),
        "unknown individual direction was silently filled or recovered from a sum");
    std::vector<ZhangProductFloatTargetRequest> bindings;
    require(zhangReadFloatTargets(f.state.metaDataMap,bindings,reason),"new target directory invalid");
    std::size_t symbolic=0;
    for(const auto& binding:bindings)if(binding.consumer=="FLOAT_CARRIER_SYMBOLIC_BINDING")
    {
        require(!binding.active && binding.numericBindingComplete &&
            binding.symbolicSpace==ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE,
            "FLOAT symbolic carrier was published as an active integer constraint");
        ++symbolic;
    }
    require(symbolic==1,"actual carrier phase definition was not bound for whole-target queries");
    const auto before=f.state;
    require(zhangCompressGraphProductFloatTargets(trace,f.state,reason),"no-op second compression failed");
    rawEqual(f.state,before);
}

static void v2ActiveRowsBeyondInactive128()
{
    KFState state;state.time.bigTime=1405216800;
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="V2_129_ACTIVE_RAW_TARGET";
    state.kfIndexMap.clear();state.kfIndexMap[KFState::oneKey]=0;
    std::vector<std::pair<KFKey,ZhangGraphEdge>> raw;
    for(int i=0;i<129;++i)
    {
        const ZhangGraphEdge edge{"R"+std::to_string(i),SatSys("G01")};
        const auto key=zhangGraphOldArcAuxiliaryKey(state,E_ObsCode::L1C,edge,0);
        state.kfIndexMap[key]=i+1;raw.emplace_back(key,edge);
        const auto identity=zhangHandoffFloatIdentityText({edge,int(E_ObsCode::L1C),0});
        state.metaDataMap[std::string(ZHANG_HANDOFF_FLOAT_DEPENDENCY_PREFIX)+identity]=identity;
    }
    state.x=VectorXd::Ones(130);state.dx=VectorXd::Zero(130);state.P=MatrixXd::Identity(130,130);state.P(0,0)=0;
    for(const auto& [key,index]:state.kfIndexMap)state.stateTransitionMap[key][key][0]=1;
    for(const auto& [key,edge]:raw)
        require(zhangHandoffFloatRegisterRecord(state,key,edge,int(E_ObsCode::L1C),0),"129 RAW initialization failed");
    ZhangProductFloatTargetRequest target;
    target.targetId="V2_SUM_OF_129";target.consumer="PENDING_PUBLICATION_DEPENDENCY";
    target.lifecycleId="LITERAL_ACTIVE_129";target.unit="m";
    target.scope={zhangProductFloatOwnerId(state),"LITERAL_129_SOURCE",zhangFloatCarrierStateOrder(state.kfIndexMap),0};
    for(const auto& [key,edge]:raw)target.numericRow[key]=1;
    target.numericBindingComplete=true;
    std::string reason;std::ostringstream trace;
    require(zhangRegisterGraphProductHandoffTarget(state,target,reason),"129 complete reservation failed");
    zhangClearGraphProductHandoffDependencies(state); // V2 is the ONLY remaining lifetime authority
    state.time=state.time+1000;
    const auto before=state;
    require(retireZhangGraphProductEventAuxiliaries(trace,state,"V2_ONLY_ACTIVE_OVER_128"),"V2-only RAW lifetime sweep failed");
    rawEqual(state,before);
    require(zhangHandoffFloatDependencies(state).size()==129,"active V2 directions were charged to the inactive128 budget");
}

static void callbackFailurePreservesEverything()
{
    Fixture f("V2_CARRIER_CALLBACK_ROLLBACK");std::string reason;std::ostringstream trace;
    require(zhangRegisterGraphProductHandoffTarget(f.state,f.reservation(),reason),"rollback reservation registration failed");
    int calls=0;
    f.state.exactStateTransformCallback=[&](const KFState&,GTime,const std::map<KFKey,int>&,
        const std::map<KFKey,int>&,const SparseMatrix<double>&,const std::string& label,
        const VectorXd&,const MatrixXd&,std::uint64_t,std::uint64_t)
        {++calls;require(label==ZHANG_FLOAT_TARGET_COMPRESSION_LABEL,"compression used permissive retirement label");return false;};
    f.state.lastFactorTransactionFailureReason="UNRELATED_PREEXISTING_DIAGNOSTIC";
    const auto before=f.state;
    require(zhangCompressGraphProductFloatTargets(trace,f.state,reason),"optional callback failure poisoned FLOAT");
    rawEqual(f.state,before);
    require(calls==1 && f.state.lastFactorTransactionFailureReason==before.lastFactorTransactionFailureReason,
        "optional failure changed callback ownership or authoritative FLOAT error");
}

static void productionWholeSpanKeepsRealNuisanceAndCurrentRhs()
{
    Fixture f("V2_PRODUCTION_WHOLE_SPAN");
    // Actual current FLOAT coordinates for the two crossing links. The old
    // diagonal links below are compressed by the real production transaction.
    KFKey c,d;c.type=KF::PHASE_BIAS;c.Sat=SatSys("G02");c.str="A_CURRENT_PHI";
    d.type=KF::PHASE_BIAS;d.Sat=SatSys("G01");d.str="B_CURRENT_PHI";
    f.state.kfIndexMap[c]=4;f.state.kfIndexMap[d]=5;
    f.state.x.conservativeResize(6);f.state.x(4)=11;f.state.x(5)=13;
    f.state.dx=VectorXd::Zero(6);
    MatrixXd q=MatrixXd::Zero(6,6);q.topLeftCorner(4,4)=f.state.P;q(4,4)=2;q(5,5)=3;f.state.P=q;
    f.state.stateTransitionMap[c][c][0]=1;f.state.stateTransitionMap[d][d][0]=1;
    std::string reason;std::ostringstream trace;
    require(zhangRegisterGraphProductHandoffTarget(f.state,f.reservation(),reason),"whole-span original reservation failed");
    require(zhangCompressGraphProductFloatTargets(trace,f.state,reason) &&
        reason=="GLOBAL_WHOLE_FLOAT_TARGET_SPACE_COMPRESSED","whole-span real carrier compression failed");
    KFKey carrier;bool haveCarrier=false;
    for(const auto& [key,index]:f.state.kfIndexMap)if(zhangIsProductFloatCarrierKey(key)){carrier=key;haveCarrier=true;}
    require(haveCarrier,"whole-span test has no actual typed carrier");
    // Keep the archived carrier definition, but deliberately expose only the
    // honest actual row S=(Phi_a+Phi_b)/2+g. Its explicit g is a real clock
    // coordinate, not an invented integer or a discarded nuisance column.
    for(auto i=f.state.metaDataMap.begin();i!=f.state.metaDataMap.end();)
    {
        ZhangProductFloatTargetRequest old;
        if(i->first.starts_with(ZHANG_FLOAT_TARGET_V2_PREFIX) && zhangDecodeFloatTarget(i->second,old) &&
            old.consumer=="FLOAT_CARRIER_SYMBOLIC_BINDING")i=f.state.metaDataMap.erase(i);
        else ++i;
    }
    ZhangProductFloatTargetRequest s;
    s.targetId="ACTUAL_HALF_CARRIER_PLUS_G";s.consumer="PENDING_PUBLICATION_DEPENDENCY";
    s.lifecycleId="LITERAL_NUMERICAL_SUPPORT";s.unit="m";
    s.symbolicSpace=ZhangProductFloatTargetSpace::VERSIONED_PHASE_EDGE;
    s.physicalCoefficients[{{"A",SatSys("G01")},E_ObsCode::L1C,0}]=1;
    s.physicalCoefficients[{{"B",SatSys("G02")},E_ObsCode::L1C,0}]=1;
    s.physicalScale=ZhangBridgeSpanRational(1)/2;s.nuisanceCoefficients["g"]=1;
    s.numericRow={{carrier,0.5},{f.clock,1}};s.numericBindingComplete=true;
    s.scope={zhangProductFloatOwnerId(f.state),"ACTUAL_CARRIER_SOURCE",
        zhangFloatCarrierStateOrder(f.state.kfIndexMap),0};
    std::vector<ZhangProductFloatTargetRequest> beforeBindings;
    require(zhangReadFloatTargets(f.state.metaDataMap,beforeBindings,reason),"whole-span binding directory unreadable");
    require(!beforeBindings.empty(),"whole-span fixture lost all bound targets");
    s.scope.coordinateGeneration=beforeBindings.front().scope.coordinateGeneration;
    require(zhangRegisterGraphProductHandoffTarget(f.state,s,reason),"actual half-carrier plus g registration failed");

    const std::vector<std::string> ids={"L1C|A|G01|V0","L1C|B|G02|V0","L1C|A|G02|V0","L1C|B|G01|V0"};
    const std::array<E_ObsCode,2> codes={E_ObsCode::L1C,E_ObsCode::L2W};
    const double lambda=CLIGHT/FREQ1;
    auto arcs=[&](const std::string& id,std::map<KFKey,double>& row,std::string& why) {
        if(id==ids[2] || id==ids[3]){row={{id==ids[2]?c:d,1/lambda}};why="CURRENT_VERSIONED_ARC_ROW";return true;}
        why="RETIRED_ARC_AUXILIARY_STATE_MISSING";return false;
    };
    int verified=0;auto current=[&](std::string&){++verified;return true;};
    const ZhangBridgeSpanMatrix cycle={{1,1,-1,-1}};
    // Independent literal oracle: eliminating the two known crossing Phi rows
    // leaves [lambda/2,lambda/2, g] for S but [1,1,0g] for t. No coefficient
    // can satisfy both old N columns and the g column.
    auto negative=zhangResolveWholeFloatSupportProduction(f.state,E_Sys::GPS,codes,ids,{}, {},cycle,{0},arcs,current);
    require(!negative.valid,"production missing-M shortcut discarded the real explicit g nuisance");
    auto gauge=s;gauge.targetId="ACTUAL_G_ROW";gauge.physicalCoefficients.clear();
    gauge.physicalScale=1;gauge.numericRow={{f.clock,1}};
    require(zhangRegisterGraphProductHandoffTarget(f.state,gauge,reason),"actual gauge cancellation row registration failed");
    auto positive=zhangResolveWholeFloatSupportProduction(f.state,E_Sys::GPS,codes,ids,{}, {},cycle,{0},arcs,current);
    require(positive.valid && positive.oldAux && positive.rows.size()==1,"actual complete nuisance cancellation was rejected");
    const auto& row=positive.rows.front();
    auto coefficient=[&](const KFKey& key){const auto i=row.find(key);return i==row.end()?0:i->second;};
    require(std::abs(coefficient(carrier)-1/lambda)<1e-12 && std::abs(coefficient(c)+1/lambda)<1e-12 &&
        std::abs(coefficient(d)+1/lambda)<1e-12 && coefficient(f.clock)==0 && positive.offsets[0]==0,
        "production support differs from literal C/lambda-Phi_c/lambda-Phi_d/lambda oracle");
    require(!zhangResolveWholeFloatSupportProduction(f.state,E_Sys::GPS,codes,ids,{}, {},{{1,0,0,0}},{0},arcs,current).valid,
        "whole carrier recovered an unknown individual old ambiguity");
    const ZhangExactMatrix h={{1,1,-1,-1}};
    const auto domainImage=zhangR51PhysicalImage({{2,2,-2,-2}}, {},{},{{2,2,-2,-2}},4);
    require(domainImage.valid && domainImage.projector.size()==1 && domainImage.generators==ZhangExactMatrix{{2}},
        "literal nonsaturated 2Z image construction failed");
    const auto imageBefore=domainImage;
    require(zhangResolveWholeFloatSupportProduction(f.state,E_Sys::GPS,codes,ids,{}, {},
        domainImage.projector,domainImage.offsets,arcs,current).valid,"actual PAR image rows failed whole-support assembly");
    auto allMissing=[](const std::string&,std::map<KFKey,double>&,std::string& why){why="RETIRED_ARC_AUXILIARY_STATE_MISSING";return false;};
    auto rhs3=zhangResolveWholeFloatSupportProduction(f.state,E_Sys::GPS,codes,ids,h,{3},cycle,{4},allMissing,current);
    auto rhs5=zhangResolveWholeFloatSupportProduction(f.state,E_Sys::GPS,codes,ids,h,{5},cycle,{4},allMissing,current);
    require(rhs3.valid && rhs5.valid && rhs3.rows[0].empty() && rhs5.rows[0].empty() &&
        rhs3.offsets[0]==7 && rhs5.offsets[0]==9 && verified==2,"support did not use the actual CURRENT b and exact offset");
    require(domainImage.projector==imageBefore.projector && domainImage.offsets==imageBefore.offsets &&
        domainImage.generators==imageBefore.generators && domainImage.particularTarget==imageBefore.particularTarget &&
        domainImage.integerImageSource==imageBefore.integerImageSource,
        "numerical whole support changed integer image/PAR coordinates");
}

int main()
{
    setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF","1",1);
    setenv("ZHANG_R51_PRODUCT_REPAIR","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
    setenv("ZHANG_R51_RATIO_THRESHOLD","2",1);setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1);
    acsConfig.zhangFullRank.enable=true;acsConfig.zhangPppAr.user_adapter=false;
    activeDependenciesAreNotCacheEvictions();v2ActiveRowsBeyondInactive128();
    realCompressionAndWholeResolution();callbackFailurePreservesEverything();
    productionWholeSpanKeepsRealNuisanceAndCurrentRhs();
    std::cout<<"PASS actual_float_carrier_production real_rank_compression=1 whole_sum_without_raw_arcs=1 "
        "independent_unknown_direction_rejected=1 literal_phase_clock_fullQ=1 V2_active_RAW_pinned=1 "
        "UNKNOWN_preserves_RAW=1 V2_only_active_over_128=1 strict_callback_rollback=1 current_integer_rhs_unchanged=1 "
        "repeated_noop_does_not_grow=1 production_whole_span_actual_typed_C_scale=1 hidden_nuisance_rejected=1 "
        "actual_gauge_row_cancellation=1 CURRENT_b_only=1 original_integer_image_unchanged=1\n";
}
