#include "common/algebra.hpp"
#include "common/zhangUserResetCheckpoint.hpp"
#include "common/receiver.hpp"
#include "pea/zhangReceiverCheckpoint.hpp"
#include <iostream>

extern "C" int __wrap_main(int,char**)
{
    unsigned checks=0;
    auto require=[&](bool value,const char* reason){++checks;if(!value)throw std::runtime_error(reason);};
    auto rejected=[&](auto operation){bool failed=false;try{operation();}catch(...){failed=true;}require(failed,"negative accepted");};
    try
    {
        namespace fs=std::filesystem;const auto root=fs::current_path()/"user_reset_checkpoint_fixture";
        require(fs::create_directory(root),"fresh fixture required");
        setenv("ZHANG_US20_RESET_CHECKPOINT","1",1);setenv("ZHANG_USER_DATUM_POLICY","TRANSPORT_KNOWN",1);
        setenv("ZHANG_US20_FORK_BRANCH","BASE",1);unsetenv("ZHANG_US20_DUMP_PROVENANCE");
        const std::string cfg="EXACT_SCIENTIFIC_CONFIGURATION",inputs="EXACT_GEC_PRODUCTS",binary="/proc/self/exe";
        KFState state;state.time.bigTime=100000;const auto x=state.x;const auto P=state.P;
        boost::json::object branch{{"canonical_config_sha256",zhangCheckpointSha256(cfg)},
            {"input_manifest_sha256",zhangCheckpointSha256(inputs)},{"policy","TRANSPORT_KNOWN"},
            {"restore_from","BASE"},{"restore_epoch",state.time.to_string(0)}};
        boost::json::object plan{{"schema","US20_GEC_RESET_FORK_PLAN_V1"},
            {"binary_sha256",zhangCheckpointFileSha256(binary)},
            {"branches",boost::json::object{{"BASE",branch},{"RESUME",branch}}}};
        const auto text=boost::json::serialize(plan);const auto file=root/"plan.json";
        std::ofstream(file)<<text;setenv("ZHANG_US20_FORK_PLAN",file.c_str(),1);
        setenv("ZHANG_US20_FORK_PLAN_SHA256",zhangCheckpointSha256(text).c_str(),1);
        auto c=cfg,i=inputs;require(zhangUs20AuthoriseForkProvenance(c,i,binary),"exact approved provenance");
        require(c.find(text)!=std::string::npos && i.find(text)!=std::string::npos,"common immutable fork identity");
        rejected([&]{auto bad=cfg+"changed threshold",in=inputs;zhangUs20AuthoriseForkProvenance(bad,in,binary);});
        rejected([&]{auto config=cfg,bad=inputs+"changed product";zhangUs20AuthoriseForkProvenance(config,bad,binary);});
        setenv("ZHANG_USER_DATUM_POLICY","RESET_KNOWN_LOCAL",1);
        rejected([&]{auto config=cfg,in=inputs;zhangUs20AuthoriseForkProvenance(config,in,binary);});
        setenv("ZHANG_USER_DATUM_POLICY","TRANSPORT_KNOWN",1);
        const auto source=root/"source_ledger",destination=root/"restored_ledger";
        const std::string claim=zhangCheckpointSha256("MTA1:100000"),data="BOUND_FULL_Q_NOISE_ID\nMTA1\n100000\n";
        fs::create_directories(source/claim);std::ofstream(source/claim/"source.txt")<<data;
        setenv("ZHANG_USER_TROP_CONSUMPTION_LEDGER",source.c_str(),1);zhangUs20CaptureLedger(state);
        require(state.x==x && state.P==P,"ledger capture cannot change stochastic state");
        KFState untouched;untouched.time=state.time;
        const auto metadataBefore=untouched.metaDataMap;
        const auto ledgerPayload=zhangUs20CaptureLedgerPayload(untouched);
        require(untouched.metaDataMap==metadataBefore,"supplement capture cannot modify live metadata");
        setenv("ZHANG_USER_TROP_CONSUMPTION_LEDGER",destination.c_str(),1);setenv("ZHANG_US20_FORK_BRANCH","RESUME",1);
        zhangUs20RestoreLedger(state);require(zhangUs20ReadFile(destination/claim/"source.txt")==data,"exact raw consumption receipt");
        zhangUs20RestoreLedgerPayload(ledgerPayload,state.time,true);
        zhangUs20RestoreLedger(state);require(zhangUs20LedgerEntries(destination,100000).size()==1,"repeat restore idempotent");
        auto bad=state;bad.metaDataMap["US20_RESET_CHECKPOINT_SOURCE_BRANCH"]="UNAPPROVED";rejected([&]{zhangUs20RestoreLedger(bad);});
        bad=state;bad.time.bigTime+=30;rejected([&]{zhangUs20RestoreLedger(bad);});
        bad=state;bad.metaDataMap["US20_RESET_CHECKPOINT_LEDGER_V1"]="{\"../escape\":\"X\\nMTA1\\n100000\\n\"}";
        rejected([&]{zhangUs20RestoreLedger(bad);});require(!fs::exists(root/"escape"),"invalid ledger path rejected before writing");
        const auto extra=zhangCheckpointSha256("MTA1:99900");fs::create_directory(destination/extra);
        std::ofstream(destination/extra/"source.txt")<<"X\nMTA1\n99900\n";
        rejected([&]{zhangUs20RestoreLedger(state);});
        rejected([&]{zhangUs20ValidateLedgerItem(zhangCheckpointSha256("MTA1:100030"),"X\nMTA1\n100030\n",100000);});
        ReceiverMap receivers;
        auto& receiver=receivers["PKRD"];receiver.id="PKRD";
        for(const char* id:{"G01","E01","C01"})
        {
            auto& q=receiver.satStatMap[SatSys(id)];
            q.userAcceptedGf.accepted={{99900,0.02,0.003},{99930,0.021,0.003},{99960,0.022,0.003}};
            q.userAcceptedGf.consecutiveSuspects=23;q.userAcceptedGf.rawSeen=true;
            q.userAcceptedGf.lastRawTime=99990;q.userAcceptedGf.firstSuspectThreshold=0.05;
            q.userAcceptedMw=q.userAcceptedCodeGf=q.userAcceptedGf;
            q.userGfDecision={false,false,false,0.023,1.0,0.05,"PERSISTENT_UNRESOLVED_INNOVATION"};
            q.userGfAvailable=true;q.userQcFrequencyA=1;q.userQcFrequencyB=2;
            q.sigStatMap["L1"].userQcState=ZhangUserQcState::SUSPECT;
            q.sigStatMap["L2"].userPhysicalBreakThisEpoch=true;
        }
        const auto qc=exportZhangUs20QcCheckpoint(receivers,"US20_TEST",100000);
        auto restored=receivers;
        for(auto& [sat,q]:restored.at("PKRD").satStatMap)
        {
            q.userAcceptedGf={};q.userAcceptedMw={};q.userAcceptedCodeGf={};q.userGfDecision={};
            q.userGfAvailable=false;q.userQcFrequencyA=q.userQcFrequencyB=0;
            for(auto& [signal,s]:q.sigStatMap){s.userQcState=ZhangUserQcState::ACCEPTED;s.userPhysicalBreakThisEpoch=false;}
        }
        importZhangUs20QcCheckpoint(qc,restored,"US20_TEST",100000);
        require(exportZhangUs20QcCheckpoint(restored,"US20_TEST",100000)==qc,"full GEC causal QC exact restore");
        for(const char* id:{"G01","E01","C01"})
        {
            auto original=receivers.at("PKRD").satStatMap.at(SatSys(id)).userAcceptedGf;
            auto resumed=restored.at("PKRD").satStatMap.at(SatSys(id)).userAcceptedGf;
            const auto a=zhangInspectUserTrend(original,100020,1.0,0.003,0.05,5,0.000001,120);
            const auto b=zhangInspectUserTrend(resumed,100020,1.0,0.003,0.05,5,0.000001,120);
            require(!a.accepted && !b.accepted && !b.warmup && std::string(a.reason)==b.reason,
                "restoration must not readmit unresolved data through warmup");
        }
        rejected([&]{importZhangUs20QcCheckpoint(qc,restored,"OTHER_RUNTIME",100000);});
        rejected([&]{importZhangUs20QcCheckpoint(qc,restored,"US20_TEST",99900);});
        auto incomplete=restored;incomplete.at("PKRD").satStatMap.erase(SatSys("E01"));
        rejected([&]{importZhangUs20QcCheckpoint(qc,incomplete,"US20_TEST",100000);});
        auto invalid=receivers;invalid.at("PKRD").satStatMap.at(SatSys("G01")).userAcceptedGf.accepted[0].time=100030;
        rejected([&]{exportZhangUs20QcCheckpoint(invalid,"US20_TEST",100000);});
        std::cout<<"PASS US20 reset checkpoint exact branch provenance and TROP ledger checks="<<checks<<'\n';return 0;
    }
    catch(const std::exception& error){std::cerr<<"FAIL "<<error.what()<<'\n';return 1;}
}
