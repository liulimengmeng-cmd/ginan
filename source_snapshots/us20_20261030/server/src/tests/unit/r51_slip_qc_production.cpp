// Linked against the actual PEA objects, including acsQC.cpp and checkpoint
// serialization. No replacement implementation of the production adapter.
#include "common/acsConfig.hpp"
#include "common/acsQC.hpp"
#include "common/navigation.hpp"
#include "common/receiver.hpp"
#include "common/observations.hpp"
#include "common/satStat.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangRatioOnly.hpp"
#include "common/rinex.hpp"
#include "common/ephPrecise.hpp"
#include "common/ephemeris.hpp"
#include "common/sinex.hpp"
#include "orbprop/coordinates.hpp"
#include "pea/zhangReceiverCheckpoint.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool condition,const std::string& text) {if(!condition)throw std::runtime_error(text);}
static const std::array<double,3> wavelength{299792458./1575420000,299792458./1227600000,299792458./1176450000};
struct Row {double time,range,el,l1,l2,l5,p1,p2,phaseVar,codeVar;int satellite,lli1=0,lli2=0,lli5=0;};
static int rawReplay(const std::string& root,int n1,int n2,int n3)
{
    std::string sp3=root+"/products/WUM0MGXFIN_20241990000_01D_05M_ORB.SP3";
    readSp3ToNav(sp3,nav,0);
    std::ifstream clock(root+"/products/WUM0MGXFIN_20241990000_01D_30S_CLK.CLK");
    require(bool(clock),"precise clock file missing");
    char clockType=0;double clockVersion=0;E_Sys clockSys=E_Sys::NONE;E_TimeSys clockTime=E_TimeSys::GPST;
    std::map<E_Sys,std::map<int,CodeType>> clockCodes;RinexStation clockStation;ObsList ignored;
    readRnx(clock,clockType,ignored,nav,clockStation,clockVersion,clockSys,clockTime,clockCodes);
    std::ifstream rinex(root+"/data/MAL200KEN_R_20241990000_01D_30S_MO.rnx");
    require(bool(rinex),"MAL2 raw observation missing");
    char type=0;double version=0;E_Sys system=E_Sys::NONE;E_TimeSys timeSystem=E_TimeSys::GPST;
    std::map<E_Sys,std::map<int,CodeType>> codes;RinexStation station;
    Receiver receiver;receiver.id="MAL2";
    receiver.source=root+"/data/MAL200KEN_R_20241990000_01D_30S_MO.rnx";
    setenv("ZHANG_R51_SLIP_RAW_WARMUP","1",1);
    int epochs=0,quarantines=0,breaks=0;double detected=-1;bool seeded=false;
    while(rinex.good())
    {
        ObsList parsed;
        readRnx(rinex,type,parsed,nav,station,version,system,timeSystem,codes);
        if(parsed.empty())continue;
        const double seconds=parsed.front()->time.bigTime-1405209600;
        if(seconds<6600)continue;
        if(seconds>7800)break;
        receiver.aprioriPos=station.pos;
        require(receiver.aprioriPos.norm()>6e6,"RINEX receiver position not read");
        ObsList observations;std::ostringstream trace;
        for(auto& obs:only<GObs>(parsed))
        {
            if(obs.Sat.sys!=E_Sys::GPS)continue;
            obs.mount="MAL2";obs.rec_ptr=&receiver;obs.satStat_ptr=&receiver.satStatMap[obs.Sat];
            obs.satNav_ptr=&nav.satNavMap[obs.Sat];updateLamMap(obs.time,obs);
            for(const auto& [frequency,signals]:obs.sigsLists)
                for(const auto& signal:signals)
                    if(signal.code==E_ObsCode::L1C || signal.code==E_ObsCode::L2W)
                        obs.sigs[frequency]=signal;
            if(!obs.sigs.count(F1) || !obs.sigs.count(F2))continue;
            if(obs.Sat.prn==25 && seconds>=7260)
            {
                obs.sigs[F1].L+=n1;obs.sigs[F2].L+=n2;
                for(auto& signal:obs.sigsLists[F1])signal.L+=n1;
                for(auto& signal:obs.sigsLists[F2])signal.L+=n2;
                for(auto& signal:obs.sigsLists[F5])signal.L+=n3;
            }
            // Force the production adapter's independent PRECISE path, even
            // though the main network normally prefers a KALMAN clock.
            obs.clkSource=E_Source::KALMAN;
            GObs external=obs;
            require(satPosClk(trace,obs.time,external,nav,{E_Source::PRECISE},{E_Source::PRECISE},
                nullptr,nullptr,E_OffsetType::COM),"production precise ephemeris failed");
            obs.rSatApc=external.rSatApc;obs.ephPosValid=obs.ephClkValid=true;obs.posSource=E_Source::PRECISE;
            VectorEcef los;geodist(obs.rSatApc,receiver.aprioriPos,los);
            satazel(ecef2pos(receiver.aprioriPos),los,*obs.satStat_ptr);
            if(obs.satStat_ptr->el<10*D2R)continue;
            const double factor=std::sin(obs.satStat_ptr->el);
            for(auto& [frequency,signal]:obs.sigs)
            {signal.phasVar=std::pow(.002/factor,2);signal.codeVar=std::pow(.4/factor,2);}
            for(auto& [frequency,signals]:obs.sigsLists)
                for(auto& signal:signals)
                {signal.phasVar=std::pow(.002/factor,2);signal.codeVar=std::pow(.4/factor,2);}
            // Exercise migration from a previous raw LC without new detector
            // state at exactly the resumed 02:00 epoch.
            if(seconds==7200)obs.satStat_ptr->slipEvidence={};
            observations.push_back(static_cast<std::shared_ptr<GObs>>(obs));
        }
        clearSlips(observations);
        for(auto& obs:only<GObs>(observations))
        {obs.satStat_ptr->lc_pre=obs.satStat_ptr->lc_new;obs.satStat_ptr->lc_new={};}
        obs2lcs(trace,observations);
        detectslips(trace,observations);
        const auto text=trace.str();seeded |= text.find("ZHANG_SLIP_RAW_WARMUP rec=MAL2 epochs=")!=std::string::npos;
        std::cout<<text;
        for(auto& obs:only<GObs>(observations))
        {
            if(obs.Sat.prn!=25)continue;
            auto& first=obs.satStat_ptr->sigStatMap[ft2string(F1)];
            auto& second=obs.satStat_ptr->sigStatMap[ft2string(F2)];
            const bool quarantined=first.phaseQuarantined || second.phaseQuarantined;
            const bool broken=first.slip.any || second.slip.any;
            if(!n1 && !n2 && !n3)require(!quarantined && !broken,"real RINEX/precise MAL2 G25 false alarm");
            if(seconds>=7260)
            {
                quarantines+=quarantined;breaks+=broken;
                if(detected<0 && (quarantined || broken))detected=seconds;
            }
            // This harness ends at QC; consuming pending flags here is not a
            // replacement for a graph/KF transaction or a product certificate.
        }
        epochs++;
    }
    require(epochs==41 && seeded,"raw replay incomplete or checkpoint seed path not exercised");
    if(n1 || n2)require(detected>=0,"raw production replay missed primary slip injection");
    std::cout<<"\nraw_production_replay=PASS epochs="<<epochs<<" checkpoint_seed="<<seeded
        <<" target_quarantines="<<quarantines<<" primary_break_epochs="<<breaks
        <<" detection_delay_s="<<(detected<0 ? -1 : detected-7260)<<" graph_kf_replayed=0\n";
    return 0;
}
static std::map<double,std::vector<Row>> synthetic()
{
    std::map<double,std::vector<Row>> out;
    for(int epoch=0;epoch<20;epoch++)
        for(int sat=20;sat<=25;sat++)
        {
            const double iono=2+.01*epoch+.001*epoch*epoch;
            const double range=2e7;
            Row row;row.time=epoch*30;row.satellite=sat;row.range=range;row.el=45;
            row.l1=(range-iono)/wavelength[0]+(sat==25 && epoch>=18 ? 1 : 0);
            row.l2=(range-std::pow(wavelength[1]/wavelength[0],2)*iono)/wavelength[1];
            row.l5=(range-std::pow(wavelength[2]/wavelength[0],2)*iono)/wavelength[2];
            if(sat==25 && epoch==19)row.lli1=1;
            row.p1=range+iono;row.p2=range+std::pow(wavelength[1]/wavelength[0],2)*iono;
            row.phaseVar=4e-6;row.codeVar=.16;out[row.time].push_back(row);
        }
    return out;
}
int main(int argc,char** argv)
{
    try
    {
        setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
        setenv("ZHANG_R51_EVIDENCE_SLIP_QC","1",1);
        require(zhangRatioOnly(),"network ratio-only policy disabled");
        acsConfig.zhangFullRank.enable=true;acsConfig.epoch_interval=30;traceLevel=3;
        acsConfig.process_sys[E_Sys::GPS]=true;
        acsConfig.code_priorities[E_Sys::GPS]={E_ObsCode::L1C,E_ObsCode::L2W};
        auto& stationOptions=acsConfig.getRecOpts("MAL2");
        stationOptions.elevation_mask_deg=10;stationOptions.code_sigma=.4;stationOptions.phase_sigma=.002;
        stationOptions.error_model=E_NoiseModel::ELEVATION_DEPENDENT;
        for(const std::string code:{"L1C","L2W","L5Q"})
        {
            auto& options=acsConfig.getRecOpts("MAL2",{"GPS",code});
            options.code_sigma=.4;options.phase_sigma=.002;options.error_model=E_NoiseModel::ELEVATION_DEPENDENT;
        }
        if(argc==6 && std::string(argv[1])=="--raw")
            return rawReplay(argv[2],std::stoi(argv[3]),std::stoi(argv[4]),std::stoi(argv[5]));
        if(argc==3 && std::string(argv[1])=="--old-checkpoint")
        {
            ZhangCheckpointBundle bundle;
            const auto io=readZhangCheckpointBundle(argv[2],{},bundle);
            require(io.valid,io.failureReason);
            const auto& section=bundle.sections.at(ZHANG_RECEIVER_RUNTIME_CHECKPOINT_SECTION_NAME);
            require(section.schemaVersion==2,"fixture is not old receiver schema 2");
            ReceiverMap empty;ZhangReceiverRuntimeCheckpointRestorePlan plan;
            const auto result=preflightZhangReceiverRuntimeCheckpointSection(empty,
                bundle.manifest.runtimeId,section.payload,plan);
            // Decode and semantic validation occur BEFORE destination inventory
            // checking. This is explicitly NOT an import/restored filter claim.
            require(!result.valid && result.failureReason=="RECEIVER_RUNTIME_CHECKPOINT_RECEIVER_INVENTORY_MISMATCH",
                result.failureReason);
            require(result.receiverCount==180,"old checkpoint receiver count changed");
            std::cout<<"old_checkpoint_decode=PASS receivers="<<result.receiverCount
                <<" satellites="<<result.satelliteStatusCount<<" import_tested=0\n";
            return 0;
        }
        const bool selfTest=argc==1;
        auto epochs=selfTest ? synthetic() : std::map<double,std::vector<Row>>{};
        int n1=0,n2=0,n3=0;
        if(!selfTest)
        {
            require(argc==5,"usage: contract fixture n1 n2 n3");
            n1=std::stoi(argv[2]);n2=std::stoi(argv[3]);n3=std::stoi(argv[4]);
            std::ifstream input(argv[1]);require(bool(input),"fixture not found");std::string line;
            while(std::getline(input,line))
            {
                if(line.empty() || line[0]=='#')continue;
                Row row;std::istringstream stream(line);
                require(bool(stream>>row.time>>row.satellite>>row.range>>row.el>>row.l1>>row.l2>>row.l5
                    >>row.p1>>row.p2>>row.phaseVar>>row.codeVar>>row.lli1>>row.lli2>>row.lli5),"invalid raw fixture");
                if(row.satellite==25 && row.time>=7260){row.l1+=n1;row.l2+=n2;if(row.l5)row.l5+=n3;}
                epochs[row.time].push_back(row);
            }
        }
        ReceiverMap receivers;auto& receiver=receivers["MAL2"];receiver.id="MAL2";
        // Install real resolver-owned fixture records; the built-in dummy
        // pointers intentionally cannot satisfy a checkpoint import contract.
        auto& site=theSinex.mapsiteids["MAL2"];site.sitecode="MAL2";
        auto& equipment=theSinex.mapreceivers["MAL2"][GTime{}];equipment.sitecode="MAL2";
        auto& antenna=theSinex.mapantennas["MAL2"][GTime{}];antenna.sitecode="MAL2";
        auto& eccentricity=theSinex.mapeccentricities["MAL2"][GTime{}];eccentricity.sitecode="MAL2";
        eccentricity.ecc=Vector3d::Zero();
        receiver.snx.id_ptr=&site;receiver.snx.rec_ptr=&equipment;
        receiver.snx.ant_ptr=&antenna;receiver.snx.ecc_ptr=&eccentricity;
        receiver.aprioriPos=Vector3d(6378137,0,0);
        int quarantines=0,breaks=0;
        for(const auto& [time,rows]:epochs)
        {
            ObsList observations;
            for(const auto& row:rows)
            {
                auto obs=std::make_shared<GObs>();obs->time.bigTime=1405209600+time;obs->mount="MAL2";
                obs->Sat=SatSys(E_Sys::GPS,row.satellite);obs->rec_ptr=&receiver;
                obs->satStat_ptr=&receiver.satStatMap[obs->Sat];obs->satStat_ptr->el=row.el*D2R;
                obs->satNav_ptr=&nav.satNavMap[obs->Sat];
                obs->satNav_ptr->lamMap={{F1,wavelength[0]},{F2,wavelength[1]},{F5,wavelength[2]}};
                obs->rSatApc=Vector3d(6378137+row.range,0,0);obs->rSatCom=obs->rSatApc;obs->satClk=0;
                obs->ephPosValid=obs->ephClkValid=true;obs->clkSource=obs->posSource=E_Source::PRECISE;
                Sig one;one.code=E_ObsCode::L1C;one.L=row.l1;one.P=row.p1;one.LLI=row.lli1;
                one.phasVar=row.phaseVar;one.codeVar=row.codeVar;
                Sig two=one;two.code=E_ObsCode::L2W;two.L=row.l2;two.P=row.p2;two.LLI=row.lli2;
                obs->sigs={{F1,one},{F2,two}};obs->sigsLists[F1].push_back(one);obs->sigsLists[F2].push_back(two);
                if(row.l5){Sig five=one;five.code=E_ObsCode::L5Q;five.L=row.l5;five.LLI=row.lli5;obs->sigsLists[F5].push_back(five);}
                GObs bound=*obs;
                observations.push_back(static_cast<std::shared_ptr<GObs>>(bound));
            }
            clearSlips(observations);std::ostringstream trace;detectslips(trace,observations);
            if(!selfTest)std::cout<<trace.str();
            for(auto& item:observations)
            {
                auto& obs=dynamic_cast<GObs&>(*item);if(obs.Sat.prn!=25)continue;
                auto& first=obs.satStat_ptr->sigStatMap[ft2string(F1)];
                auto& second=obs.satStat_ptr->sigStatMap[ft2string(F2)];
                quarantines+=first.phaseQuarantined || second.phaseQuarantined;
                breaks+=first.slip.any!=0 || second.slip.any!=0;
                require(!first.slip.GF && !second.slip.GF && !first.slip.MW && !second.slip.MW,
                    "legacy GF/MW path overruled new evidence detector");
                if(selfTest && time==540)
                {
                    require(first.phaseQuarantined || first.slip.SCDIA,"real adapter silently used L1 injection");
                    require(!second.slip.any && !second.tracking.pendingBreak,"uncertain L1 injection retired healthy L2");
                }
                if(selfTest && time==570)
                {
                    require(first.slip.LLI && first.tracking.pendingBreak,"production adapter missed explicit L1 LLI");
                    require(!second.slip.any && !second.tracking.pendingBreak,"L1 LLI retired healthy L2");
                }
                if(!selfTest && n1==0 && n2==0 && n3==0)
                    require(!first.phaseQuarantined && !second.phaseQuarantined && !first.slip.any && !second.slip.any,
                        "MAL2 G25 replay generated false quarantine/retirement");
            }
        }
        std::cout<<"\nproduction_adapter=PASS quarantines="<<quarantines<<" primary_break_epochs="<<breaks<<'\n';
        // Verify schema-3 export can serialize the actual bounded detector
        // state. No configured-destination import is claimed here.
        if(selfTest)
        {
            std::string payload;
            const auto result=exportZhangReceiverRuntimeCheckpointSection(receivers,"SLIP-QC-TEST",payload);
            require(result.valid,result.failureReason);
            ReceiverMap empty;ZhangReceiverRuntimeCheckpointRestorePlan plan;
            const auto checked=preflightZhangReceiverRuntimeCheckpointSection(empty,"SLIP-QC-TEST",payload,plan);
            require(checked.failureReason=="RECEIVER_RUNTIME_CHECKPOINT_RECEIVER_INVENTORY_MISMATCH",checked.failureReason);
            ReceiverMap destination=receivers;
            destination.at("MAL2").satStatMap.clear();
            const auto prepared=preflightZhangReceiverRuntimeCheckpointSection(destination,"SLIP-QC-TEST",payload,plan);
            require(prepared.valid,prepared.failureReason);
            const auto imported=importZhangReceiverRuntimeCheckpointSection(destination,"SLIP-QC-TEST",plan);
            require(imported.valid,imported.failureReason);
            std::string reexported;
            const auto roundtrip=exportZhangReceiverRuntimeCheckpointSection(destination,"SLIP-QC-TEST",reexported);
            require(roundtrip.valid && reexported==payload,"schema3 full import did not preserve detector bytes");
            std::cout<<"schema3_export_import_roundtrip=PASS bytes="<<payload.size()<<'\n';
        }
        return 0;
    }
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
