#include "common/algebra.hpp"
#include "common/zhangUserArExperiment.hpp"
#include "common/zhangUserProductChart.hpp"
#include "common/zhangUs20ArCheckpoint.hpp"
#include "common/zhangUserArPolicy.hpp"
#include "common/zhangUserGonsis.hpp"
#include <boost/json.hpp>
#include <iostream>

extern "C" int __wrap_main(int,char**)
{
    unsigned checks=0;
    auto require=[&](bool ok,const char* reason){++checks;if(!ok)throw std::runtime_error(reason);};
    auto rejected=[&](auto op){bool bad=false;try{op();}catch(...){bad=true;}require(bad,"invalid checkpoint accepted");};
    try
    {
        setenv("ZHANG_USER_DATUM_POLICY","TRANSPORT_KNOWN",1);
        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE","GONSIS_UDUC",1);
        setenv("ZHANG_USER_QC","CAUSAL_TRACKING",1);setenv("ZHANG_USER_AR_SEARCH","GONSIS_UDUC",1);
        setenv("ZHANG_USER_AR_FEEDBACK","EPOCH_ONLY",1);
        for(const auto system:{E_Sys::GPS,E_Sys::GAL,E_Sys::BDS})
        {
            const SatSys sat(system,system==E_Sys::BDS?20:7);const int code=1;
            KFState state;ZhangUserGonsisComponentCertificate cert{true,true,"UNCHANGED_ANCHOR",3,2};
            zhangUserRememberChartAnchor(state,"USER",sat,code,cert);
            const auto id=zhangUserProductHistoryIdentity(state,"USER",sat,code,cert);
            require(!id.empty(),"authorized canonical history missing");
            ZhangUserGonsisReferenceHint hint{sat.id(),id,true};
            std::vector<ZhangUserGonsisSatellite> choices(2);
            choices[0].satellite=sat.id();choices[0].identity=id;choices[0].elevationRadians=0.8;
            choices[1].satellite=SatSys(system,25).id();choices[1].identity="OTHER_VALID_ARC";choices[1].elevationRadians=1.1;
            const std::vector<int> active{0,1};
            require(zhangUserGonsisDetail::reference(active,choices,hint)==0,"valid previous reference ignored");
            for(double offset:{20.,13.,-7.,0.})
            {
                zhangUserSetChartOffset(state,"USER",sat,code,offset);
                require(zhangUserProductHistoryIdentity(state,"USER",sat,code,cert)==id,"coordinate gauge reset canonical identity");
                choices[0].identity=zhangUserProductHistoryIdentity(state,"USER",sat,code,cert);
                require(zhangUserGonsisDetail::reference(active,choices,hint)==0,"gauge broke actual reference hint selection");
            }
            choices[1].effectiveSamples=200;
            require(zhangUserGonsisDetail::reference(active,choices,hint)==0,"published fix count threshold changed");
            choices[1].effectiveSamples=201;
            require(zhangUserGonsisDetail::reference(active,choices,hint)==1,"published fix count failed to affect reference choice");
            choices[1].effectiveSamples=0;
            const auto saved=state.metaDataMap;const auto key=zhangUserChartKey("USER",sat,code);
            state.metaDataMap.erase(key+":OFFSET");require(!zhangUserChartHistoryAuthorized(state,"USER",sat,code),"missing offset accepted");state.metaDataMap=saved;
            state.metaDataMap[key+":UNANCHORED_MAP"]="1";require(!zhangUserChartHistoryAuthorized(state,"USER",sat,code),"unanchored accepted");state.metaDataMap=saved;
            zhangIncrementUserArcGeneration(state.metaDataMap,zhangUserIntegerAuthorizationGenerationKey("USER",sat.id(),code));
            require(zhangUserProductHistoryIdentity(state,"USER",sat,code,cert).empty(),"revoked authorization inherited");
            choices[0].identity=zhangUserProductHistoryIdentity(state,"USER",sat,code,cert);
            require(zhangUserGonsisDetail::reference(active,choices,hint)==1,"revoked reference hint still consumed");
            zhangUserRememberChartAnchor(state,"USER",sat,code,cert);
            require(!zhangUserChartHistoryAuthorized(state,"USER",sat,code),"current certificate silently reauthorized old anchor");
            zhangUserSetChartOffset(state,"USER",sat,code,0);
            require(!zhangUserChartHistoryAuthorized(state,"USER",sat,code),"zero net shift resurrected revoked history");
            zhangUserRetireChartAnchor(state,"USER",sat,code);
            zhangIncrementUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(false,"USER",sat.id(),code));
            zhangUserRememberChartAnchor(state,"USER",sat,code,cert);
            require(zhangUserChartHistoryAuthorized(state,"USER",sat,code) && zhangUserProductHistoryIdentity(state,"USER",sat,code,cert)!=id,"new state reused old authorization history");
            zhangIncrementUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(true,"USER",sat.id(),code));
            require(!zhangUserChartHistoryAuthorized(state,"USER",sat,code),"physical break inherited anchor");
            state.metaDataMap=saved;unsetenv("ZHANG_USER_DATUM_POLICY");
            const auto legacy=zhangUserProductHistoryIdentity(state,"USER",sat,code,cert);
            zhangUserSetChartOffset(state,"USER",sat,code,20);
            require(zhangUserProductHistoryIdentity(state,"USER",sat,code,cert)!=legacy,"legacy identity changed implicitly");
            setenv("ZHANG_USER_DATUM_POLICY","TRANSPORT_KNOWN",1);
        }
        ZhangUserArExperimentPolicy p;p.minimumArcSeconds=0;p.minimumArcSamples=1;
        ZhangUserArSeparatedArcTracker continuous(p),restored(p);
        const ZhangUserArExperimentKey key{"TEST/branch/user-epoch-only","USER","G07"};
        for(double t:{100.,130.,160.})continuous.update(key,t,"P0",true,"ANCHOR0",true,true,"S0",true);
        const auto snapshot=continuous.checkpoint("TEST");restored.restoreCheckpoint(snapshot,"TEST",160);
        auto a=continuous.update(key,190,"P0",true,"ANCHOR0",true,true,"S0",true);
        auto b=restored.update(key,190,"P0",true,"ANCHOR0",true,true,"S0",true);
        require(a.physical.count==4 && b.physical.count==4 && b.product.age==90 && b.effective.age==90,"restart lost causal ages");
        b=restored.update(key,190,"P0",true,"ANCHOR0",true,true,"S0",true);require(b.effective.count==4,"duplicate increments history");
        b=restored.update(key,220,"P0",true,"ANCHOR0",true,true,"S1",true);
        require(b.physical.count==5 && b.product.count==5 && b.state.count==1 && b.effective.count==1,"local state release changed physical or retained state age");
        auto bad=snapshot;bad[0].arcs.begin()->second.last=170;
        rejected([&]{restored.restoreCheckpoint(bad,"TEST",160);});
        require(restored.checkpoint("TEST")[0].arcs.begin()->second.count==5,"failed restore changed live history");
        bad=snapshot;bad[3].policy.minimumArcSamples=2;rejected([&]{restored.restoreCheckpoint(bad,"TEST",160);});
        rejected([&]{restored.restoreCheckpoint(snapshot,"OTHER",160);});

        const string root="AR_CHECKPOINT_TEST";const double time=160;
        const auto original=exportZhangUs20UserArCheckpoint(root,time);
        auto doc=boost::json::parse(original).as_object();
        auto arc=boost::json::object{{"runtime",root+"/branch/user-epoch-only"},{"receiver","USER"},{"object","G07"},
            {"valid",true},{"eligible",true},{"first",100.},{"last",160.},{"age",60.},{"count",3},
            {"identity","PHYSICAL_STATE_CANONICAL_IA0"},{"reason","CONTINUOUS_SAMPLE"}};
        auto tracker=doc.at("legacy_arcs").as_object();tracker["records"]=boost::json::array{arc};doc["legacy_arcs"]=tracker;
        auto separatedTracker=tracker;separatedTracker.at("policy").as_object()["confirmation_epochs"]=3;
        boost::json::array four;for(int i=0;i<4;++i)four.push_back(separatedTracker);
        doc["separated_arcs"]=boost::json::array{boost::json::object{{"minimum_seconds",0.},{"minimum_samples",1},{"trackers",four}}};
        auto group=boost::json::object{{"runtime",root+"/branch/user-epoch-only"},{"receiver","USER"},{"system",int(E_Sys::GPS)},{"domain","SAME_LATTICE"}};
        auto ref=group;ref["time"]="160";ref["satellite"]="G07";ref["identity"]="PHYSICAL_STATE_CANONICAL_IA0";ref["previous_published_fix"]=true;
        doc["references"]=boost::json::array{ref};auto lock=group;
        lock["records"]=boost::json::array{boost::json::object{{"satellite","G07"},{"identity","PHYSICAL_STATE_CANONICAL_IA0"},{"count",3},{"last","160"}}};
        doc["published_fix_locks"]=boost::json::array{lock};const auto payload=boost::json::serialize(doc);
        importZhangUs20UserArCheckpoint(payload,root,time,true);
        require(exportZhangUs20UserArCheckpoint(root,time)==original,"preflight changed AR runtime");
        importZhangUs20UserArCheckpoint(payload,root,time,false);
        require(exportZhangUs20UserArCheckpoint(root,time)==payload,"complete AR checkpoint not exact");
        auto wrong=doc;wrong["runtime"]="OTHER";rejected([&]{importZhangUs20UserArCheckpoint(boost::json::serialize(wrong),root,time,false);});
        wrong=doc;wrong.at("references").as_array()[0].as_object()["time"]="190";rejected([&]{importZhangUs20UserArCheckpoint(boost::json::serialize(wrong),root,time,false);});
        wrong=doc;wrong.at("published_fix_locks").as_array()[0].as_object().at("records").as_array()[0].as_object()["count"]=-1;
        rejected([&]{importZhangUs20UserArCheckpoint(boost::json::serialize(wrong),root,time,false);});
        wrong=doc;wrong.at("references").as_array()[0].as_object()["satellite"]="E07";rejected([&]{importZhangUs20UserArCheckpoint(boost::json::serialize(wrong),root,time,false);});
        require(exportZhangUs20UserArCheckpoint(root,time)==payload,"failed preflight damaged live AR counters");
        wrong=doc;wrong.at("references").as_array()[0].as_object()["time"]="160.0";
        wrong.at("published_fix_locks").as_array()[0].as_object().at("records").as_array()[0].as_object()["count"]=4;
        rejected([&]{importZhangUs20UserArCheckpoint(boost::json::serialize(wrong),root,time,false);});
        require(exportZhangUs20UserArCheckpoint(root,time)==payload,"noncanonical payload partially committed AR counters");
        importZhangUs20UserArCheckpoint(original,root,time,false);
        require(exportZhangUs20UserArCheckpoint(root,time)==original,"empty causal partition restore not exact");
        std::cout<<"USER_AR_CHECKPOINT_AND_CANONICAL_HISTORY_PASS checks="<<checks<<'\n';return 0;
    }
    catch(const std::exception& e){std::cerr<<"FAIL checks="<<checks<<" "<<e.what()<<'\n';return 1;}
}
