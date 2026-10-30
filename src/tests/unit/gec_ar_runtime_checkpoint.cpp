// Production-linked archive contract. The independent expectations below are
// fixed physical rows/proof IDs and parsed wire fields, not a second importer.
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangMultiGnssInteger.hpp"
#include "common/zhangR51Integrity.hpp"
#include "common/zhangRatioOnly.hpp"
// Relative include also permits the explicitly labelled V6 RED test to use
// this pure size-contract helper while linking the unchanged old importer.
#include "../../cpp/common/zhangGecArCheckpointLimits.hpp"
#include "ambres/GNSSambres.hpp"
#include "pea/zhangPppAr.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/json.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <cstring>
#include <cstdlib>
#include <array>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <sstream>
using std::map; using std::set; using std::string; using std::vector;
#include "gec_private_ambres_wire_schema.hpp"
// The production refresh shared by both epoch-only AR call sites. Keeping
// this declaration here avoids widening the public checkpoint API for tests.
void refreshZhangR51CurrentPhysicalArchives(const KFState& owner);

namespace j = boost::json;
namespace
{
int assertions = 0;
#ifdef ZHANG_TEST_LEGACY_OBJECTS_NO_DISPATCH
constexpr bool productionDispatchContractCompiled=false;
#else
constexpr bool productionDispatchContractCompiled=true;
#endif
const std::string runtime = "GEC_AR_ARCHIVE_INDEPENDENT_CONTRACT";
void check(bool condition, const std::string& why)
{ ++assertions; if (!condition) throw std::runtime_error(why); }
std::string text(const j::value& value)
{ const auto& s=value.as_string(); return {s.data(),s.size()}; }

struct System { E_Sys system; const char* satellite; E_ObsCode first,second; };
const System systems[] = {{E_Sys::GPS,"G01",E_ObsCode::L1C,E_ObsCode::L2W},
    {E_Sys::GAL,"E11",E_ObsCode::L1C,E_ObsCode::L5Q},
    {E_Sys::BDS,"C19",E_ObsCode::L1P,E_ObsCode::L5P}};

void mode(bool ar)
{
    setenv("ZHANG_R51_ENABLE","1",1); setenv("ZHANG_R51_RATIO_ONLY","1",1);
    setenv("ZHANG_R51_RATIO_THRESHOLD","2",1); setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1);
    if (ar) { setenv("ZHANG_MULTI_GNSS_AR_CHECKPOINT","1",1); unsetenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT"); }
    else { unsetenv("ZHANG_MULTI_GNSS_AR_CHECKPOINT"); setenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT","1",1); }
    acsConfig.zhangFullRank.enable = true;
    for (auto& [system,enabled] : acsConfig.process_sys) enabled = false;
    acsConfig.zhangFullRank.sysOpts.clear();
    for (const auto& spec : systems)
    {
        acsConfig.process_sys[spec.system] = true;
        acsConfig.zhangFullRank.sysOpts[spec.system].baseline_observables = {spec.first,spec.second};
    }
}

KFState state()
{
    KFState result; result.time.bigTime = 1405290000;
    result.x = VectorXd::LinSpaced(6,1.25,6.25);
    result.P = MatrixXd::Constant(6,6,.125); result.P.diagonal().setConstant(2.);
    result.dx = VectorXd::Zero(6); result.kfIndexMap.clear();
    int index = 0;
    for (const auto& spec : systems) for (auto code : {spec.first,spec.second})
    {
        KFKey key; key.type=KF::AMBIGUITY; key.Sat=SatSys(spec.satellite); key.str="A";
        key.num=int(code); key.comment="FULL_KEY_METADATA"; key.estimatedTime.bigTime=result.time.bigTime-30;
        result.kfIndexMap.emplace(key,index++);
    }
    std::string reason;
    check(bindZhangCheckpointRuntimeId(result,runtime,&reason),reason);
    return result;
}

ZhangCurrentEpochPhysicalRow physical(E_ObsCode code,const System& spec,int scale=1)
{
    const auto label=enum_to_string(code); const auto sat=SatSys(spec.satellite);
    const auto other=SatSys(spec.system,sat.prn+1).id();
    return {{label+"|A|"+sat.id()+"|V2",scale},{label+"|B|"+other+"|V2",scale},
        {label+"|A|"+other+"|V2",-scale},{label+"|B|"+sat.id()+"|V2",-scale}};
}

void seed(KFState& owner)
{
    zhangCurrentEpochEvidenceRegistry().clear(); zhangMultiGnssCurrentEvidenceAuthorities().clear();
    const auto root=zhangR51NumericRoot(owner.x,owner.P);
    const auto order=zhangCurrentEvidenceStateOrderId(owner.kfIndexMap);
    ZhangMultiGnssCurrentEvidenceAuthority authority;
    authority.epochGpstSeconds=static_cast<long int>(owner.time.bigTime);
    authority.floatRoot=root; authority.floatStateOrderId=order;
    int column=0;
    for (const auto& spec : systems)
    {
        ZhangCurrentEpochEvidence e; e.runtimeId=runtime; e.system=spec.system;
        e.epochGpstSeconds=authority.epochGpstSeconds; e.posteriorId="POSTERIOR_"+enum_to_string(spec.system);
        e.chartId="CHART_"+enum_to_string(spec.system); e.floatStateRoot=root;
        e.floatPosteriorPsdValidated=true; e.floatPosteriorPsdRootId=root;
        e.floatPosteriorPsdStateOrderId=order; e.valid=true; e.chart.columns=6;
        ZhangDecisionProofPtr firstProof;
        KFKey firstKey;
        for (auto code : {spec.first,spec.second})
        {
            const auto row=physical(code,spec);
            const auto chord=enum_to_string(code)+"|A|"+SatSys(spec.satellite).id()+"|V2";
            e.chart.chordIdentities[column]=chord; e.chart.expansions[column++]=row;
            e.chart.cataloguedSignals.insert(enum_to_string(code));
            for (const auto& [id,x] : row) e.chart.activePhysicalIdentities.insert(id);
            auto parent=std::make_shared<ZhangIntegerDecisionProof>();
            parent->id="PARENT_"+enum_to_string(spec.system)+"_"+enum_to_string(code);
            parent->originalStatement="independent controlled parent";
            parent->observationProvenance=runtime+"|PPP_COMMITTED|E"+std::to_string(e.epochGpstSeconds);
            parent->conditionalFailureBound=.00001;
            auto proof=std::make_shared<ZhangIntegerDecisionProof>();
            proof->id="DECISION_"+enum_to_string(spec.system)+"_"+enum_to_string(code);
            proof->originalStatement="independent controlled integer two";
            proof->observationProvenance=parent->observationProvenance;
            proof->conditionalFailureBound=.00002; proof->parents={parent};
            KFKey identity; identity.type=KF::AMBIGUITY; identity.Sat=SatSys(spec.satellite);
            identity.str="A"; identity.num=int(code);
            const auto key=owner.kfIndexMap.find(identity)->first;
            check(e.add(row,2,{proof}),"controlled physical evidence admission");
            check(e.addNumericConstraint({{key,1}},2,row),"controlled numeric binding");
            if (code==spec.first) { firstProof=proof; firstKey=key; }
        }
        const auto overlap=physical(spec.first,spec,2);
        check(e.addEventOverlap(overlap,4,{firstProof},"INDEPENDENT_CURRENT_EVENT_WITNESS"),"event overlap");
        check(e.addNumericConstraint({{firstKey,2}},4,overlap),"event numeric binding");
        e.chart.unexpandedColumns.emplace((column+1)%6,"CONTROLLED_UNKNOWN_COLUMN_NOT_ZERO");
        authority.systemIdentities[spec.system]={e.posteriorId,e.chartId};
        zhangCurrentEpochEvidenceRegistry().emplace(std::make_pair(runtime,spec.system),std::move(e));
    }
    zhangMultiGnssCurrentEvidenceAuthorities().emplace(runtime,std::move(authority));
}

std::string exported(const KFState& owner)
{
    std::string wire,reason;
    check(exportZhangGecArRuntimeSection(owner,runtime,wire,reason),reason);
    return wire;
}

#ifndef ZHANG_TEST_LEGACY_OBJECTS_NO_DISPATCH
void perEpochArDispatchContract(KFState& owner)
{
    KFState augmented=owner;
    acsConfig.reference_clock="NO_REFERENCE";
    acsConfig.reference_bias="NO_REFERENCE";
    acsConfig.process_minimum_constraints=false;
    acsConfig.minconOpts.once_per_epoch=false;
    for(bool fullRank : {false,true})
    for(bool transactional : {false,true})
    {
        acsConfig.zhangFullRank.enable=fullRank;
        acsConfig.zhangPppAr.transactional_integer_fixing=transactional;
        for(bool hold : {false,true})
        {
            auto& selected=selectZhangPerEpochArState(owner,augmented,hold,false);
            const bool usesOwner=hold || (fullRank && transactional);
            check(&selected==(usesOwner?&owner:&augmented),
                "production per-epoch dispatch chose wrong KFState object");
        }
    }
    auto rejected=[&](const std::string& expected,bool inRts=false)
    {
        bool threw=false;
        try { (void)selectZhangPerEpochArState(owner,augmented,false,inRts); }
        catch(const std::runtime_error& error)
        { threw=std::string(error.what())==expected; }
        check(threw,"unsupported transformed-output AR route did not fail closed: "+expected);
    };
    rejected("ZHANG_TRANSACTIONAL_AR_RTS_OWNER_UNSUPPORTED",true);
    acsConfig.reference_clock="PIVOT";
    rejected("ZHANG_TRANSACTIONAL_AR_PREPIVOTED_OWNER_UNSUPPORTED");
    acsConfig.reference_clock="NO_REFERENCE";
    acsConfig.reference_bias="PIVOT";
    rejected("ZHANG_TRANSACTIONAL_AR_PREPIVOTED_OWNER_UNSUPPORTED");
    acsConfig.reference_bias="NO_REFERENCE";
    acsConfig.process_minimum_constraints=true;
    acsConfig.minconOpts.once_per_epoch=true;
    rejected("ZHANG_TRANSACTIONAL_AR_MINCON_OWNER_UNSUPPORTED");
    acsConfig.process_minimum_constraints=false;
    acsConfig.minconOpts.once_per_epoch=false;
    owner.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA]=runtime+"/branch/forged-owner";
    rejected("ZHANG_TRANSACTIONAL_AR_OWNER_IS_DISPOSABLE_BRANCH");
    owner.metaDataMap.erase(ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA);
}
#endif

// This is deliberately the production lifecycle, not a copied resolver. The
// archive-only fixture used to miss the interaction between a restored owner,
// disabled raw-factor capture, and disposable AR KFState copies. A copy keeps
// the root metadata but must never acquire the root object's service lease.
void ownerObjectLifecycleContract(KFState& owner)
{
    acsConfig.zhangPppAr.fixed_lag_factor_capture_shadow=false;
    acsConfig.zhangPppAr.targeted_besd_capture_shadow=false;
    acsConfig.zhangPppAr.user_adapter=false;
    std::string canonical;
    auto result=exportZhangPppArCheckpointSection(owner,runtime,canonical);
    check(result.valid && !canonical.empty(),"owner initial PPP export: "+result.failureReason);

    // Exact old production trigger: outputs.cpp creates augmentedKF by plain
    // copy, without derived branch metadata, while fix_and_hold=false. Turning
    // raw-factor capture off must not release the actual FLOAT service owner.
    configureZhangE18FactorCapture(owner);
    KFState augmentedOutputCopy=owner;
    check(!augmentedOutputCopy.metaDataMap.count(ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA),
        "plain output-copy fixture accidentally has a derived runtime ID");
    (void)selectZhangE18CanonicalProductRelations(augmentedOutputCopy,E_Sys::GPS,{}, {},0);
    std::string afterPlainCopy;
    result=exportZhangPppArCheckpointSection(owner,runtime,afterPlainCopy);
    check(result.valid,"COPIED_OUTPUT_STOLE_OWNER_AFTER_CAPTURE_OFF:"+result.failureReason);
    check(afterPlainCopy==canonical,"plain output copy altered authoritative PPP runtime");
    check(importZhangPppArCheckpointSection(owner,runtime,canonical,false,6).valid,
        "same-owner PPP restore with factor capture disabled");

    const char* roles[]={"fixed-transaction","wide-lane-snapshot","multi-product-fixed"};
    for(const auto role : roles)
    {
        KFState disposable=owner;
        disposable.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA]=
            runtime+"/branch/"+role;
        (void)selectZhangE18CanonicalProductRelations(disposable,E_Sys::GPS,{}, {},0);
        std::string after;
        result=exportZhangPppArCheckpointSection(owner,runtime,after);
        check(result.valid,"RESTORED_OWNER_STOLEN_BY_AR_COPY:"+result.failureReason+
            ":role="+role);
        check(after==canonical,"disposable AR copy altered the authoritative PPP partition");
    }
    for(int round=0;round<2;++round)
    {
        check(importZhangPppArCheckpointSection(owner,runtime,canonical,true,6).valid,
            "same-owner PPP preflight rejected");
        check(importZhangPppArCheckpointSection(owner,runtime,canonical,false,6).valid,
            "same-owner PPP commit rejected");
        configureZhangE18FactorCapture(owner);
        KFState conflict=owner;
        std::string rejected;
        result=exportZhangPppArCheckpointSection(conflict,runtime,rejected);
        check(!result.valid && rejected.empty() &&
            result.failureReason.find("RUNTIME_OBJECT_CONFLICT")!=std::string::npos,
            "distinct live KFState claimed the authoritative runtime ID");
        for(bool preflight : {true,false})
        {
            result=importZhangPppArCheckpointSection(conflict,runtime,canonical,preflight,6);
            check(!result.valid && result.failureReason.find("RUNTIME_OBJECT_CONFLICT")!=std::string::npos,
                "conflicting live KFState stole owner through checkpoint import");
        }
        // Even a rejected copy's lifecycle entry must not deregister the live
        // owner. The owner remains exportable after temporary copies die.
        configureZhangE18FactorCapture(conflict);
        std::string after;
        result=exportZhangPppArCheckpointSection(owner,runtime,after);
        check(result.valid && after==canonical,
            "rejected copy or disabled capture retired the authoritative owner lease");
    }
}

void inspectIndependentEvidence(const j::value& item)
{
    const auto& e=item.as_object();
    const auto system=std::string(e.at("system").as_string());
    const System* expected=nullptr;
    for (const auto& spec : systems)
        if (system==enum_to_string(spec.system)) expected=&spec;
    check(expected!=nullptr,"unknown system in independent wire inventory");
    const std::string satellite=expected->satellite;
    const std::string other=system=="GPS" ? "G02" : system=="GAL" ? "E12" : "C20";
    check(text(e.at("owner"))==runtime && text(e.at("runtime"))==runtime && text(e.at("key_system"))==system,
        "wire owner/system metadata changed");
    check(text(e.at("posterior"))=="POSTERIOR_"+system && text(e.at("chart_id"))=="CHART_"+system,
        "posterior/chart identity changed");
    check(e.at("valid").as_bool() && e.at("psd_validated").as_bool() &&
        e.at("epoch").as_int64()==1405290000,"current/PSD scope lost");
    check((e.at("generation").is_uint64() ? e.at("generation").as_uint64() :
        static_cast<std::uint64_t>(e.at("generation").as_int64()))==3,"constraint-generation changed");
    const auto& chart=e.at("chart").as_object();
    check(chart.at("columns").as_int64()==6 && chart.at("chords").as_array().size()==2 &&
        chart.at("expansions").as_array().size()==2 && chart.at("active").as_array().size()==8 &&
        chart.at("signals").as_array().size()==2,"physical chart inventory changed");
    check(chart.at("unknown").as_array().size()==1 &&
        text(chart.at("unknown").as_array()[0].as_object().at("reason"))=="CONTROLLED_UNKNOWN_COLUMN_NOT_ZERO",
        "unavailable-column diagnostic dropped or silently made zero");
    const E_ObsCode codes[]={expected->first,expected->second};
    for (int r=0;r<2;++r)
    {
        const auto signal=enum_to_string(codes[r]);
        const j::object physicalExpected{{signal+"|A|"+satellite+"|V2","1"},
            {signal+"|B|"+other+"|V2","1"},{signal+"|A|"+other+"|V2","-1"},
            {signal+"|B|"+satellite+"|V2","-1"}};
        check(e.at("rows").as_array()[r].as_object()==physicalExpected && text(e.at("values").as_array()[r])=="2",
            "physical integer equation differs from independent four-edge cycle");
        check(e.at("bindings").as_array()[r].as_object().at("physical").as_object()==physicalExpected &&
            text(e.at("bindings").as_array()[r].as_object().at("value"))=="2",
            "physical/numerical binding equation changed");
        const auto& numeric=e.at("numeric").as_array()[r].as_object();
        check(text(numeric.at("value"))=="2" && numeric.at("row").as_array().size()==1,
            "numeric conditioning integer/term count changed");
        const auto& key=numeric.at("row").as_array()[0].as_object();
        check(key.at("type").as_int64()==int(KF::AMBIGUITY) &&
            key.at("number").as_int64()==int(codes[r]) && text(key.at("satellite"))==satellite &&
            text(key.at("receiver"))=="A" && text(key.at("comment"))=="FULL_KEY_METADATA" &&
            text(key.at("receiver_pointer"))=="" && key.at("coefficient").as_double()==1.,
            "full KF key/finite numeric coefficient metadata changed");
        const auto hex=std::string(key.at("estimated_time").as_string());
        check(hex.size()==sizeof(long double)*2,"long-double key time truncated");
        std::array<unsigned char,sizeof(long double)> raw{};
        for (std::size_t i=0;i<raw.size();++i)
            raw[i]=static_cast<unsigned char>(std::stoul(hex.substr(2*i,2),nullptr,16));
        long double epoch=0; std::memcpy(&epoch,raw.data(),raw.size());
        check(epoch==1405289970.L,"estimated-time value lost");
        const auto& proofGroup=e.at("proof_groups").as_array()[r].as_array();
        check(proofGroup.size()==1 && text(proofGroup[0])=="DECISION_"+system+"_"+signal,
            "proof root identity changed");
    }
    check(text(e.at("event_values").as_array()[0])=="4" &&
        text(e.at("witnesses").as_array()[0])=="INDEPENDENT_CURRENT_EVENT_WITNESS" &&
        e.at("numeric").as_array()[2].as_object().at("row").as_array()[0].as_object()
            .at("coefficient").as_double()==2.,"controlled event-overlap value/witness changed");
    std::set<std::string> atomIds;
    for (const auto& atomValue : e.at("proof_atoms").as_array())
    {
        const auto& atom=atomValue.as_object();
        const auto id=std::string(atom.at("id").as_string());
        atomIds.insert(id);
        check(text(atom.at("provenance"))==runtime+"|PPP_COMMITTED|E1405290000",
            "proof observation provenance changed");
        const bool parent=id.rfind("PARENT_",0)==0;
        check(atom.at("risk").as_double()==(parent?.00001:.00002),"proof risk metadata changed");
        const auto& parents=atom.at("parents").as_array();
        check(parent ? parents.empty() : parents.size()==1 && text(parents[0])=="PARENT_"+id.substr(9),
            "proof parent topology changed");
    }
    const std::set<std::string> expectedIds={"PARENT_"+system+"_"+enum_to_string(expected->first),
        "PARENT_"+system+"_"+enum_to_string(expected->second),
        "DECISION_"+system+"_"+enum_to_string(expected->first),
        "DECISION_"+system+"_"+enum_to_string(expected->second)};
    check(atomIds==expectedIds,"independent proof DAG inventory changed");
}

void sharedParentContract(KFState& owner)
{
    const auto ordinary=exported(owner);
    auto parent=std::make_shared<ZhangIntegerDecisionProof>();
    parent->id="CONTROLLED_BUNDLE_SHARED_PARENT";
    parent->originalStatement="one physical observation decision shared by G/E/C";
    parent->observationProvenance=runtime+"|PPP_COMMITTED|E1405290000";
    parent->conditionalFailureBound=.000003;
    for (auto& [key,e] : zhangCurrentEpochEvidenceRegistry())
    {
        std::map<std::string,ZhangDecisionProofPtr> replacements;
        for (auto* groups : {&e.rowProofs,&e.eventOverlapProofs})
            for (auto& group : *groups) for (auto& proof : group)
            {
                auto found=replacements.find(proof->id);
                if (found==replacements.end())
                {
                    auto copy=std::make_shared<ZhangIntegerDecisionProof>(*proof);
                    copy->parents.push_back(parent);
                    found=replacements.emplace(copy->id,std::move(copy)).first;
                }
                proof=found->second;
            }
        e.invalidateCurrentDomain();
    }
    const auto shared=exported(owner);
    std::string reason;
    for (int round=0;round<2;++round)
    {
        check(importZhangGecArRuntimeSection(owner,runtime,shared,false,6,reason),reason);
        ZhangDecisionProofs jointRoots;
        for (const auto& [key,e] : zhangCurrentEpochEvidenceRegistry())
            for (const auto& group : e.rowProofs)
                jointRoots.insert(jointRoots.end(),group.begin(),group.end());
        const auto joint=zhangDecisionRiskClosure(jointRoots);
        check(joint.valid && joint.atoms.size()==13 &&
            joint.atoms.count("CONTROLLED_BUNDLE_SHARED_PARENT")==1,
            "joint G/E/C restore split one shared ancestor into colliding IDs");
        check(exported(owner)==shared,"shared-parent G/E/C roundtrip changed wire evidence");
    }
    auto bad=j::parse(shared).as_object();
    bool foundParent=false;
    for (auto& value : bad["evidence"].as_array()[1].as_object()["proof_atoms"].as_array())
        if (text(value.as_object().at("id"))=="CONTROLLED_BUNDLE_SHARED_PARENT")
        {
            value.as_object()["statement"]="FORGED_SHARED_PARENT_CONTENT";
            foundParent=true;
        }
    check(foundParent,"shared-parent negative fixture is empty");
    for (bool preflight : {true,false})
    {
        check(!importZhangGecArRuntimeSection(owner,runtime,j::serialize(bad),preflight,6,reason),
            "cross-partition shared proof content conflict accepted");
        check(exported(owner)==shared,"failed shared-parent import mutated live authority");
    }
    // Same content does not license two distinct live immutable objects to
    // claim one decision ID; export must preserve that original fail-closed rule.
    auto duplicateParent=std::make_shared<ZhangIntegerDecisionProof>(*parent);
    auto& gal=zhangCurrentEpochEvidenceRegistry().at({runtime,E_Sys::GAL});
    std::map<std::string,ZhangDecisionProofPtr> replacements;
    for (auto* groups : {&gal.rowProofs,&gal.eventOverlapProofs})
        for (auto& group : *groups) for (auto& proof : group)
        {
            auto found=replacements.find(proof->id);
            if (found==replacements.end())
            {
                auto copy=std::make_shared<ZhangIntegerDecisionProof>(*proof);
                for (auto& p : copy->parents)
                    if (p->id==duplicateParent->id) p=duplicateParent;
                found=replacements.emplace(copy->id,std::move(copy)).first;
            }
            proof=found->second;
        }
    std::string rejected;
    check(!exportZhangGecArRuntimeSection(owner,runtime,rejected,reason) && rejected.empty(),
        "export merged distinct cross-partition objects claiming one proof ID");
    check(importZhangGecArRuntimeSection(owner,runtime,ordinary,false,6,reason),reason);
    check(exported(owner)==ordinary,"shared-parent test did not restore ordinary fixture");
}

gec_test_ambres_wire::AmbresCheckpointEnvelope decodedAmbres(const std::string& wire)
{
    std::istringstream stream(wire,std::ios::binary);
    boost::archive::binary_iarchive archive(stream,boost::archive::no_header);
    gec_test_ambres_wire::AmbresCheckpointEnvelope result; archive>>result; return result;
}
std::string encodedAmbres(const gec_test_ambres_wire::AmbresCheckpointEnvelope& value)
{
    std::ostringstream stream(std::ios::binary);
    boost::archive::binary_oarchive archive(stream,boost::archive::no_header); archive<<value; return stream.str();
}

void ambresSystemSignalContract(KFState& owner)
{
    std::string wire;
    auto result=exportZhangAmbresCheckpointSection(owner,runtime,wire);
    check(result.valid,"empty AMBRES export: "+result.failureReason);
    auto snapshot=decodedAmbres(wire);
    for (const auto& spec : systems)
    {
        gec_test_ambres_wire::AmbresHeldLatticeCheckpoint lattice;
        lattice.ownerRuntimeId=runtime; lattice.system=spec.system;
        gec_test_ambres_wire::AmbresHeldRowCheckpoint row;
        row.coefficients.push_back({spec.first,"A",SatSys(spec.satellite),2,"1"}); row.value="2";
        lattice.rows.push_back(row); snapshot.heldLattices.push_back(lattice);
    }
    auto incoming=encodedAmbres(snapshot);
    result=importZhangAmbresCheckpointSection(owner,runtime,incoming,true);
    check(result.valid,"G/E/C held lattice preflight: "+result.failureReason);
    result=importZhangAmbresCheckpointSection(owner,runtime,incoming,false);
    check(result.valid,"G/E/C held lattice commit: "+result.failureReason);
    result=exportZhangAmbresCheckpointSection(owner,runtime,wire);
    check(result.valid && decodedAmbres(wire).heldLattices.size()==3,"G/E/C held lattice lost");
    const auto canonical=wire;
    result=importZhangAmbresCheckpointSection(owner,runtime,wire,false);
    check(result.valid,"second AMBRES commit: "+result.failureReason);
    result=exportZhangAmbresCheckpointSection(owner,runtime,wire);
    check(result.valid && wire==canonical,"nonempty G/E/C AMBRES not byte roundtrip");
    for (int badCase=0;badCase<4;++badCase)
    {
        auto bad=decodedAmbres(canonical);
        if(badCase==0)bad.heldLattices[1].rows[0].coefficients[0].code=E_ObsCode::L2W;
        if(badCase==1)bad.heldLattices[2].rows[0].coefficients[0].code=E_ObsCode::L5Q;
        if(badCase==2)bad.heldLattices[0].rows[0].coefficients[0].satellite=SatSys("E11");
        if(badCase==3)bad.heldLattices[2].system=E_Sys::GAL;
        check(!importZhangAmbresCheckpointSection(owner,runtime,encodedAmbres(bad),true).valid,
            "AMBRES accepted cross-system/tracking-code injection");
        check(exportZhangAmbresCheckpointSection(owner,runtime,wire).valid && wire==canonical,
            "failed AMBRES preflight mutated live state");
    }
    // The explicit AR mode is the sole extension. Disabling both GEC flags
    // exercises the legacy GPS decoder, not the FLOAT empty-state guard.
    mode(false); unsetenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT");
    check(!importZhangAmbresCheckpointSection(owner,runtime,canonical,true).valid,
        "legacy GPS decoder silently broadened to G/E/C");
    mode(true);
    check(exportZhangAmbresCheckpointSection(owner,runtime,wire).valid && wire==canonical,
        "legacy rejection changed GEC live state");
}

void allEnumPhysicalArchiveContract(KFState& owner)
{
    // Match the actual config parser: it materialises sysOpts for every E_Sys,
    // including NONE and disabled physical systems. Only G/E/C have current
    // evidence. No registry is cleared to make this fixture exportable.
    for(const auto system : magic_enum::enum_values<E_Sys>())
    {
        acsConfig.zhangFullRank.sysOpts.try_emplace(system);
        if(system!=E_Sys::GPS && system!=E_Sys::GAL && system!=E_Sys::BDS)
            acsConfig.process_sys[system]=false;
    }
    const auto privateBefore=exported(owner);
    refreshZhangR51CurrentPhysicalArchives(owner);
    std::string wire;
    auto result=exportZhangAmbresCheckpointSection(owner,runtime,wire);
    check(result.valid,"ALL_ENUM_CURRENT_R51_REFRESH_EXPORT_FAILED:"+result.failureReason);
    const auto canonical=wire;
    const auto decoded=decodedAmbres(canonical);
    check(decoded.r51Archives.size()==3,"refresh emitted missing/inactive R51 archive partitions");
    std::set<E_Sys> actualSystems;
    for(const auto& archive : decoded.r51Archives)
    {
        actualSystems.insert(archive.system);
        const System* spec=nullptr;
        for(const auto& current : systems) if(current.system==archive.system) spec=&current;
        check(spec && archive.owner==runtime && archive.valid,"R51 archive owner/system/validity changed");
        check(archive.rows.size()==2 && archive.values==std::vector<std::string>{"2","2"},
            "R51 refresh dropped exact current integer rows/values");
        std::set<std::string> expectedRoots;
        for(int row=0;row<2;++row)
        {
            const auto code=row?spec->second:spec->first;
            std::map<std::string,std::string> expected;
            for(const auto& [id,value] : physical(code,*spec)) expected[id]=value.convert_to<std::string>();
            check(archive.rows[row]==expected,"R51 exact physical four-edge row changed");
            expectedRoots.insert("DECISION_"+enum_to_string(spec->system)+"_"+enum_to_string(code));
        }
        check(std::set<std::string>(archive.roots.begin(),archive.roots.end())==expectedRoots &&
            archive.proofs.size()==4,"R51 current decision roots/proof DAG missing");
    }
    check(actualSystems==std::set<E_Sys>{E_Sys::GPS,E_Sys::GAL,E_Sys::BDS},
        "all-enum refresh broadened the active constellation set");
    for(int round=0;round<2;++round)
    {
        result=importZhangAmbresCheckpointSection(owner,runtime,canonical,true);
        check(result.valid,"R51 nonempty preflight failed:"+result.failureReason);
        result=importZhangAmbresCheckpointSection(owner,runtime,canonical,false);
        check(result.valid,"R51 nonempty commit failed:"+result.failureReason);
        result=exportZhangAmbresCheckpointSection(owner,runtime,wire);
        check(result.valid && wire==canonical,"R51 G/E/C physical rows/roots/DAG not byte roundtrip");
        check(exported(owner)==privateBefore,"R51 archive restore changed independent current authority");
    }
    const std::vector<std::pair<std::string,std::function<void(
        gec_test_ambres_wire::AmbresCheckpointEnvelope&)>>> mutations={
        {"inactive empty NONE",[](auto& d){auto x=d.r51Archives.front();x.system=E_Sys::NONE;
            x.rows.clear();x.values.clear();x.roots.clear();x.proofs.clear();d.r51Archives.push_back(x);}},
        {"inactive empty GLO",[](auto& d){auto x=d.r51Archives.front();x.system=E_Sys::GLO;
            x.rows.clear();x.values.clear();x.roots.clear();x.proofs.clear();d.r51Archives.push_back(x);}},
        {"inactive real rows",[](auto& d){auto x=d.r51Archives.front();x.system=E_Sys::GLO;d.r51Archives.push_back(x);}},
        {"foreign owner",[](auto& d){d.r51Archives.front().owner="FOREIGN";}},
        {"duplicate partition",[](auto& d){d.r51Archives.push_back(d.r51Archives.front());}},
        {"row RHS mismatch",[](auto& d){d.r51Archives.front().values.pop_back();}},
        {"cross-system physical arc",[](auto& d){d.r51Archives.front().rows.front()["L1C|A|E11|V2"]="1";}},
        {"wrong tracking code",[](auto& d){for(auto& a:d.r51Archives)if(a.system==E_Sys::GAL)
            a.rows.front()["L2W|A|E11|V2"]="1";}},
        {"zero coefficient",[](auto& d){d.r51Archives.front().rows.front().begin()->second="0";}},
        {"noninteger coefficient",[](auto& d){d.r51Archives.front().rows.front().begin()->second="1/2";}},
        {"noninteger RHS",[](auto& d){d.r51Archives.front().values.front()="1/2";}},
        {"missing proof root",[](auto& d){d.r51Archives.front().roots.push_back("MISSING");}},
        {"missing proof parent",[](auto& d){d.r51Archives.front().proofs.front().parents.push_back("MISSING");}},
        {"duplicate proof ID",[](auto& d){auto& a=d.r51Archives.front();a.proofs.push_back(a.proofs.front());}}
    };
    for(const auto& [name,mutate] : mutations)
    {
        auto bad=decoded;mutate(bad);
        const auto payload=encodedAmbres(bad);
        for(bool preflight : {true,false})
        {
            result=importZhangAmbresCheckpointSection(owner,runtime,payload,preflight);
            check(!result.valid,"R51 accepted invalid archive: "+name);
            refreshZhangR51CurrentPhysicalArchives(owner);
            result=exportZhangAmbresCheckpointSection(owner,runtime,wire);
            check(result.valid && wire==canonical,"failed R51 import/refresh mutated live archives: "+name);
            check(exported(owner)==privateBefore,"failed R51 import changed current authority: "+name);
        }
    }
    // Populate an invalid *live* archive through the explicitly legacy mode,
    // whose historical decoder is outside the GEC capability. The new refresh
    // must not erase this evidence to make GEC export pass.
    auto invalidLive=decoded;
    invalidLive.heldLattices.clear(); // legacy held-lattice decoder remains GPS-only
    auto inactive=decoded.r51Archives.front();inactive.system=E_Sys::NONE;
    inactive.rows.clear();inactive.values.clear();inactive.roots.clear();inactive.proofs.clear();
    invalidLive.r51Archives.push_back(inactive);
    unsetenv("ZHANG_MULTI_GNSS_AR_CHECKPOINT");unsetenv("ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT");
    result=importZhangAmbresCheckpointSection(owner,runtime,encodedAmbres(invalidLive),false);
    check(result.valid,"legacy invalid-live fixture setup failed:"+result.failureReason);
    setenv("ZHANG_MULTI_GNSS_AR_CHECKPOINT","1",1);
    refreshZhangR51CurrentPhysicalArchives(owner);
    result=exportZhangAmbresCheckpointSection(owner,runtime,wire);
    check(!result.valid && wire.empty() &&
        result.failureReason.find("AMBRES_GEC_R51_ARCHIVE_SCOPE_INVALID")==0,
        "refresh erased/accepted an existing inactive live R51 archive");
    // Restore the tested valid envelope explicitly; do not use a global clear.
    check(importZhangAmbresCheckpointSection(owner,runtime,canonical,false).valid,
        "R51 fixture valid-envelope restore failed");
    check(exported(owner)==privateBefore,"physical archive test changed current evidence authority");
}

void gecPayloadCapacityContract(KFState& owner)
{
    const std::size_t cap=std::size_t{1024}*1024*1024;
    check(ZHANG_GEC_AR_RUNTIME_MAX_PAYLOAD_BYTES==cap,"shared GEC wire cap is not the reviewed 1 GiB");
    for(const auto bytes : {std::size_t{0},std::size_t{1},std::size_t{256}*1024*1024,
        std::size_t{256}*1024*1024+1,cap-1,cap,cap+1,std::numeric_limits<std::size_t>::max()})
        check(zhangGecArRuntimePayloadSizeValid(bytes)==(bytes!=0 && bytes<=cap),
            "shared GEC wire boundary checker differs at bytes="+std::to_string(bytes));
    // Valid JSON remains valid with trailing whitespace. Cross the old hard
    // limit without manufacturing millions of proofs or allocating the 1 GiB
    // cap. The old native importer rejects this before parsing; the fixed
    // importer must perform its usual scope/proof/core checks and preserve the
    // canonical private evidence. This is not the real 490 MB replay gate.
    const auto canonical=exported(owner);
    std::string padded=canonical;
    padded.resize(std::size_t{256}*1024*1024+1,' ');
    std::string reason;
    check(importZhangGecArRuntimeSection(owner,runtime,padded,true,6,reason),
        "VALID_GEC_JSON_OVER_256M_PREFLIGHT_FAILED:"+reason);
    check(exported(owner)==canonical,"large-wire preflight modified current evidence");
    check(importZhangGecArRuntimeSection(owner,runtime,padded,false,6,reason),
        "VALID_GEC_JSON_OVER_256M_COMMIT_FAILED:"+reason);
    check(exported(owner)==canonical,"large-wire commit changed canonical rows/proof DAG/core scope");
}

void losslessBinary64JsonContract(KFState& owner)
{
    // Oracles are fixed IEEE-754 encodings, NOT values parsed by Boost.JSON.
    // The first encoding is the actual e135 C36/L5P coefficient whose old
    // default parser changed c0278b980bc64528 into c0278b980bc64529 (one ULP).
    struct Case { std::uint64_t coefficient, risk; const char* label; };
    const Case cases[] = {
        {0xc0278b980bc64528ULL,0x3fc6b099e4f4edceULL,"REAL_E135_NEGATIVE_17_DIGIT"},
        {0x40278b980bc64528ULL,0x3fc6b099e4f4edceULL,"POSITIVE_17_DIGIT"},
        {0x3ff0000000000001ULL,0x0000000000000001ULL,"NEXTAFTER_ONE_MIN_SUBNORMAL_RISK"},
        {0xbff0000000000001ULL,0x8000000000000000ULL,"NEGATIVE_NEXTAFTER_ONE_SIGNED_ZERO_RISK"},
        {0x0c0c490bd79fb61fULL,0x3fc6b099e4f4edceULL,"TINY_FINITE_17_DIGIT"},
        {0xf3db96d386d6c8b2ULL,0x3fc6b099e4f4edceULL,"LARGE_FINITE_17_DIGIT"},
        {0x0000000000000001ULL,0x0000000000000001ULL,"MIN_SUBNORMAL"},
        {0x0010000000000000ULL,0x0010000000000000ULL,"MIN_NORMAL"},
        {0x7fefffffffffffffULL,0x3fc6b099e4f4edceULL,"MAX_FINITE"}};
    const auto ordinary=exported(owner);
    const auto originalIndex=owner.kfIndexMap;
    const VectorXd originalX=owner.x, originalDx=owner.dx;
    const MatrixXd originalP=owner.P;
    Receiver sourceReceiver, reboundReceiver;
    sourceReceiver.id=reboundReceiver.id="A";
    auto receiverKeys=[&](Receiver* receiver)
    {
        std::map<KFKey,int> index;
        for(const auto& [key,column]:originalIndex)
        { auto copy=key;copy.rec_ptr=receiver;index.emplace(std::move(copy),column); }
        owner.kfIndexMap=std::move(index);
    };
    std::string reason;
    for(const auto& item:cases)
    {
        owner.kfIndexMap=originalIndex;
        check(importZhangGecArRuntimeSection(owner,runtime,ordinary,false,6,reason),
            "binary64 fixture initial restore:"+reason);
        receiverKeys(&sourceReceiver);
        // Seed the real current registries and modify BOTH copies of each
        // numerical relation, preserving the physical/numerical link exactly.
        seed(owner);
        map<string,ZhangDecisionProofPtr> proofs;
        auto cloneProof=[&](auto&& self,const ZhangDecisionProofPtr& input)->ZhangDecisionProofPtr
        {
            const auto prior=proofs.find(input->id);
            if(prior!=proofs.end())return prior->second;
            auto copy=std::make_shared<ZhangIntegerDecisionProof>(*input);
            copy->parents.clear();
            for(const auto& parent:input->parents)copy->parents.push_back(self(self,parent));
            copy->conditionalFailureBound=std::bit_cast<double>(item.risk);
            proofs.emplace(copy->id,copy);return copy;
        };
        for(auto& [key,evidence]:zhangCurrentEpochEvidenceRegistry())
        {
            std::vector<KFKey> keys;
            for(const auto& [numeric,column]:owner.kfIndexMap)
                if(numeric.Sat.sys==key.second)keys.push_back(numeric);
            check(keys.size()==2,"binary64 fixture needs two real owner keys per system");
            ZhangCoordinateRow replacement;
            // Reverse insertion is deliberate; canonical map ordering must
            // remain stable, while both nonempty receiver pointers rebind.
            replacement.emplace(keys[1],0.125);
            replacement.emplace(keys[0],std::bit_cast<double>(item.coefficient));
            evidence.numericConstraints.front().rowCycles=replacement;
            bool bound=false;
            for(auto& binding:evidence.physicalNumericBindings)
                if(binding.physical==evidence.rows.front())
                { binding.rowCycles=replacement;bound=true; }
            check(bound,"binary64 fixture lacks independent physical binding");
            for(auto* groups:{&evidence.rowProofs,&evidence.eventOverlapProofs})
                for(auto& group:*groups)for(auto& proof:group)proof=cloneProof(cloneProof,proof);
            check(evidence.validatePhysicalNumericBindings(reason),
                "binary64 fixture invalid before production export:"+reason);
        }
        const auto canonical=exported(owner);
        const auto verify=[&]()
        {
            for(const auto& [key,evidence]:zhangCurrentEpochEvidenceRegistry())
            {
                const auto& row=evidence.numericConstraints.front().rowCycles;
                check(row.size()==2,"binary64 import dropped multi-key map term");
                check(std::bit_cast<std::uint64_t>(row.begin()->second)==item.coefficient,
                    string("GEC_AR_BINARY64_IMPORT_CHANGED coefficient:")+item.label);
                check(std::bit_cast<std::uint64_t>(std::next(row.begin())->second)==0x3fc0000000000000ULL,
                    "binary64 import changed independently exact second map coefficient");
                for(const auto& [numeric,coefficient]:row)
                    check(numeric.rec_ptr==&reboundReceiver && numeric.rec_ptr!=&sourceReceiver,
                        "binary64 import did not bind current owner receiver pointer");
                bool bound=false;
                for(const auto& binding:evidence.physicalNumericBindings)
                    if(binding.physical==evidence.rows.front())
                    {
                        check(binding.rowCycles==row,"binary64 binding and numerical map differ");
                        check(std::bit_cast<std::uint64_t>(binding.rowCycles.begin()->second)==item.coefficient,
                            string("GEC_AR_BINARY64_IMPORT_CHANGED binding:")+item.label);
                        bound=true;
                    }
                check(bound,"binary64 import dropped physical binding");
                ZhangDecisionProofs roots;
                for(const auto& group:evidence.rowProofs)roots.insert(roots.end(),group.begin(),group.end());
                const auto closure=zhangDecisionRiskClosure(roots);
                check(closure.valid,"binary64 import broke proof references/pointer interning");
                for(const auto& [id,proof]:closure.atoms)
                    check(std::bit_cast<std::uint64_t>(proof->conditionalFailureBound)==item.risk,
                        string("GEC_AR_BINARY64_IMPORT_CHANGED proof risk:")+item.label);
            }
        };
        receiverKeys(&reboundReceiver);
        for(int round=0;round<2;++round)
        {
            check(importZhangGecArRuntimeSection(owner,runtime,canonical,true,6,reason),
                "binary64 native preflight failed:"+reason);
            check(exported(owner)==canonical,"binary64 preflight mutated complete current wire");
            check(importZhangGecArRuntimeSection(owner,runtime,canonical,false,6,reason),
                "binary64 native commit failed:"+reason);
            verify();
            check(exported(owner)==canonical,"GEC_AR_BINARY64_COMPLETE_WIRE_CHANGED");
        }
        check(std::memcmp(owner.x.data(),originalX.data(),sizeof(double)*originalX.size())==0 &&
            std::memcmp(owner.dx.data(),originalDx.data(),sizeof(double)*originalDx.size())==0 &&
            std::memcmp(owner.P.data(),originalP.data(),sizeof(double)*originalP.size())==0,
            "binary64 archive changed FLOAT x/dx/fullP");
    }
    owner.kfIndexMap=originalIndex;
    check(importZhangGecArRuntimeSection(owner,runtime,ordinary,false,6,reason),
        "binary64 fixture final restore:"+reason);
    check(exported(owner)==ordinary,"binary64 fixture did not restore complete original authority");
}
}

int main()
try
{
    mode(true); auto owner=state();
#ifndef ZHANG_TEST_LEGACY_OBJECTS_NO_DISPATCH
    perEpochArDispatchContract(owner);
#endif
    ownerObjectLifecycleContract(owner); seed(owner);
    losslessBinary64JsonContract(owner);
    const VectorXd x=owner.x; const MatrixXd p=owner.P;
    const auto original=exported(owner); auto doc=j::parse(original);
    // Independent expected wire content, including all systems and graph rows.
    const auto& evidence=doc.as_object().at("evidence").as_array();
    check(evidence.size()==3,"three nonempty current-evidence partitions missing");
    check(doc.as_object().at("authorities").as_array().size()==1,"joint authority missing");
    for (const auto& item : evidence)
    {
        const auto& e=item.as_object();
        check(e.at("rows").as_array().size()==2,"two exact current rows lost");
        check(e.at("event_rows").as_array().size()==1,"event evidence lost");
        check(e.at("bindings").as_array().size()==3,"physical/numerical links lost");
        check(e.at("proof_atoms").as_array().size()==4,"proof DAG lost/shared parent duplicated");
        check(e.at("numeric").as_array().size()==3,"numeric conditioning recipe lost");
        inspectIndependentEvidence(item);
    }
    std::string reason;
    check(importZhangGecArRuntimeSection(owner,runtime,original,true,6,reason),reason);
    check(exported(owner)==original,"preflight modified live evidence");
    zhangCurrentEpochEvidenceRegistry().clear(); zhangMultiGnssCurrentEvidenceAuthorities().clear();
    check(zhangGecArRuntimeRegistriesEmpty(),"empty-registry classifier false");
    check(importZhangGecArRuntimeSection(owner,runtime,original,false,6,reason),reason);
    check(!zhangGecArRuntimeRegistriesEmpty(),"nonempty registry silently dropped");
    check(exported(owner)==original,"nonempty archive not byte-identical after first roundtrip");
    for (const auto& [key,value] : zhangCurrentEpochEvidenceRegistry())
        check(value.domainCache && !value.domainCache->domain && value.domainCache->builds==0 &&
            value.domainCache->hits==0,"checkpoint imported compiled exact cache as authority");
    check(importZhangGecArRuntimeSection(owner,runtime,original,false,6,reason),reason);
    check(exported(owner)==original,"nonempty archive not byte-identical after second roundtrip");
    const std::vector<std::pair<std::string,std::function<void(j::object&)>>> mutations={
        {"schema",[](auto& d){d["schema"]=2;}},
        {"extra field",[](auto& d){d["unknown_authority"]=true;}},
        {"wrong runtime",[](auto& d){d["runtime"]="FOREIGN";}},
        {"wrong epoch",[](auto& d){d["epoch"]=1405290030;}},
        {"wrong dimension",[](auto& d){d["dimension"]=7;}},
        {"wrong float root",[](auto& d){d["float_root"]="FORGED";}},
        {"wrong state order",[](auto& d){d["state_order"]="FORGED";}},
        {"duplicate evidence",[](auto& d){auto& a=d["evidence"].as_array();auto copy=a[0];a.emplace_back(copy);}},
        {"missing numeric binding",[](auto& d){d["evidence"].as_array()[0].as_object()["bindings"].as_array().clear();}},
        {"wrong RHS",[](auto& d){d["evidence"].as_array()[0].as_object()["bindings"].as_array()[0].as_object()["value"]="9";}},
        {"proof orphan",[](auto& d){auto& a=d["evidence"].as_array()[0].as_object()["proof_atoms"].as_array();auto copy=a[0];copy.as_object()["id"]="ORPHAN";a.emplace_back(copy);}},
        {"proof cycle",[](auto& d){auto& p=d["evidence"].as_array()[0].as_object()["proof_atoms"].as_array()[0].as_object();p["parents"]=j::array{p["id"]};}},
        {"missing proof root",[](auto& d){d["evidence"].as_array()[0].as_object()["proof_groups"].as_array()[0]=j::array{"ABSENT"};}},
        {"wrong proof epoch",[](auto& d){d["evidence"].as_array()[0].as_object()["proof_atoms"].as_array()[0].as_object()["provenance"]="OLD_EPOCH";}},
        {"foreign numerical signal",[](auto& d){auto& e=d["evidence"].as_array()[1].as_object();e["numeric"].as_array()[0].as_object()["row"].as_array()[0].as_object()["number"]=int(E_ObsCode::L2W);}},
        {"unknown column made expanded",[](auto& d){auto& e=d["evidence"].as_array()[0].as_object();e["chart"].as_object()["unknown"].as_array()[0].as_object()["column"]=0;}},
        {"wrong GEC physical signal",[](auto& d){auto& e=d["evidence"].as_array()[1].as_object();e["rows"].as_array()[0].as_object()["L2W|A|E11|V2"]="1";}},
        {"authority false posterior",[](auto& d){d["authorities"].as_array()[0].as_object()["systems"].as_object()["GAL"]=j::array{"FALSE","FALSE"};}},
        {"missing authority evidence",[](auto& d){d["evidence"].as_array().erase(d["evidence"].as_array().begin());}}
    };
    for (const auto& [name,mutate] : mutations)
    {
        auto bad=j::parse(original).as_object(); mutate(bad); const auto wire=j::serialize(bad);
        for (bool preflight : {true,false})
        {
            check(!importZhangGecArRuntimeSection(owner,runtime,wire,preflight,6,reason),"accepted corruption: "+name);
            check(exported(owner)==original,"failed import mutated live state: "+name);
        }
    }
    for (const auto& bad : {original+"TRAILING",original.substr(0,original.size()-1),std::string{}})
        check(!importZhangGecArRuntimeSection(owner,runtime,bad,false,6,reason),"malformed wire accepted");
    sharedParentContract(owner);
    mode(false);
    std::string rejected;
    check(!exportZhangGecArRuntimeSection(owner,runtime,rejected,reason),"FLOAT exported GEC integer authority");
    check(!importZhangGecArRuntimeSection(owner,runtime,original,false,6,reason),"FLOAT imported GEC integer authority");
    std::string empty;
    check(!exportZhangAmbresCheckpointSection(owner,runtime,empty).valid,"FLOAT silently ignored live current evidence");
    mode(true);
    ambresSystemSignalContract(owner);
    allEnumPhysicalArchiveContract(owner);
    gecPayloadCapacityContract(owner);
    struct NisResult { bool valid=false; double nis=std::numeric_limits<double>::quiet_NaN();std::string status; };
    bool computed=false;
    auto nis=zhangRatioOnlySkipNis([&]()->NisResult{computed=true;throw std::runtime_error("NIS MUST NOT RUN");});
    check(!computed && !nis.valid && std::isnan(nis.nis),"NIS evaluated or fabricated statistical pass");
    check(zhangR51EffectiveRatioThreshold(3)==2,"ratio2 not preserved");
    check(std::memcmp(owner.x.data(),x.data(),sizeof(double)*x.size())==0 &&
        std::memcmp(owner.P.data(),p.data(),sizeof(double)*p.size())==0,"archive touched FLOAT x/fullP");
    std::cout<<"GEC_AR_ARCHIVE_RESULT "<<j::serialize(j::object{{"status","PASS"},{"assertions",assertions},
        {"nonempty_systems",3},{"complete_private_evidence_roundtrips",2},
        {"shared_parent_joint_GEC_roundtrips",2},
        {"owner_object_lifecycle_after_restore",true},
        {"production_per_epoch_dispatch_contract",productionDispatchContractCompiled},
        {"disposable_AR_copy_roles_tested",3},{"conflicting_live_object_rejected",true},
        {"all_enum_current_R51_archive_refresh",true},{"R51_nonempty_systems",3},
        {"R51_archive_roundtrips",2},{"R51_invalid_wire_cases",14},
        {"existing_invalid_R51_archive_not_erased",true},
        {"GEC_shared_payload_size_contract",true},{"GEC_size_boundary_cases",8},
        {"GEC_over_old_256MiB_native_preflight_and_commit",true},
        {"GEC_large_wire_bytes",std::size_t{256}*1024*1024+1},
        {"GEC_1GiB_cap_boundary_checked_without_cap_allocation",true},
        {"GEC_binary64_fixed_pattern_cases",9},{"GEC_binary64_production_roundtrips_per_case",2},
        {"GEC_actual_e135_coefficient_bits","c0278b980bc64528"},
        {"GEC_nonempty_receiver_pointer_rebinding",true},{"GEC_multiterm_numeric_binding_maps",true},
        {"GEC_hard_double_proof_risk_bits_exact",true},
        {"full_float_x_P_byte_exact",true},{"NIS_computations",0},{"synthetic_archive_contract_only",true}})<<'\n';
    return 0;
}
catch(const std::exception& error)
{std::cerr<<"GEC_AR_ARCHIVE_FAIL "<<error.what()<<'\n';return 1;}
