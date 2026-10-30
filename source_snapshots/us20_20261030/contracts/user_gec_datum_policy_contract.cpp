#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangUserAppliedChartReceipt.hpp"
#include "common/zhangUserProductChart.hpp"
#include "pea/zhangPppAr.hpp"
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Production linked, not a copy of the adapter. All G/E/C and two bands are
// present in every trial; the event mask changes only the gauge event.
extern "C" int __wrap_main(int,char**)
{
 unsigned checks=0,trials=0;int active=-1;
 auto require=[&](bool value,const std::string& reason){++checks;if(!value)throw std::runtime_error(reason);};
 try {
  setenv("ZHANG_USER_QC","CAUSAL_TRACKING",1);setenv("ZHANG_USER_FLOAT_TRANSPORT","CONTRACT",1);
  setenv("ZHANG_USER_UNKNOWN_CHART","RECOVER",1);setenv("ZHANG_USER_AR_SEARCH","GONSIS_UDUC",1);
  setenv("ZHANG_USER_AR_FEEDBACK","EPOCH_ONLY",1);unsetenv("ZHANG_RESET_AUDIT_DIRECTORY");
  acsConfig.process_ppp=true;acsConfig.zhangPppAr.user_adapter=true;acsConfig.zhangPppAr.output_products=false;
  acsConfig.zhangPppAr.product_solution="PRODUCT_FIXED";acsConfig.zhangPppAr.user_use_full_product_covariance=false;
  acsConfig.zhangPppAr.fixed_lag_factor_capture_shadow=false;acsConfig.zhangPppAr.targeted_besd_capture_shadow=false;
  acsConfig.zhangFullRank.enable=false;acsConfig.phaseClockOsb.enable=false;
  acsConfig.zhangPppAr.integer_strategy="CANONICAL_USER_SD_WL_L1";
  acsConfig.pppOpts.ionoOpts.use_if_combo=false;
  const std::array<E_Sys,3> systems{E_Sys::GPS,E_Sys::GAL,E_Sys::BDS};
  const std::array<std::array<E_ObsCode,2>,3> codes{{{E_ObsCode::L1C,E_ObsCode::L2W},{E_ObsCode::L1C,E_ObsCode::L5Q},{E_ObsCode::L1P,E_ObsCode::L5P}}};
  const std::array<std::array<double,2>,3> frequencies{{{1575420000.,1227600000.},{1575420000.,1176450000.},{1575420000.,1176450000.}}};
  auto satAt=[&](int g,int sat){return SatSys(systems[g],sat+(g==2?18:0));};
  for(const auto policy:{"TRANSPORT_KNOWN","RESET_KNOWN_LOCAL"})for(int mask:{1,2,4,7})for(int scenario=0;scenario<14;++scenario)
  {
   ++trials;active=scenario;setenv("ZHANG_USER_DATUM_POLICY",policy,1);
   acsConfig.zhangPppAr.baseline_observables.clear();for(int g=0;g<3;++g)acsConfig.zhangPppAr.baseline_observables[systems[g]]={codes[g][0],codes[g][1]};
   auto integer=[&](int g,int sat,int band)->long long {
    if(!(mask&(1<<g)) || scenario==11)return 0;
    if(scenario==2)return band?2:5;
    if(scenario==1)return sat==1?(band?2:5):sat==2?(band?13:20):sat==4?(band?7:-3):0;
    return sat==2?(band?13:20):0;
   };
   auto endpoint=[&](int epoch,int g,int sat,int band){
    ZhangProductDatumEndpoint e;e.epoch=100000+30*epoch;e.counter=epoch;e.datum=1;
    e.satellite=satAt(g,sat).id();e.observable=enum_to_string(codes[g][band]);e.reference=satAt(g,1).id();
    e.segment=e.satellite+e.observable+"S"+std::to_string(epoch);e.physicalFunctional="FROZEN_PHYSICAL";
    e.component="GEC_"+enum_to_string(systems[g]);e.componentVersion=1;e.alignment=0;
    e.posterior="ALL_GEC_POSTERIOR_"+std::to_string(epoch);e.bundle="ALL_GEC_QP_"+std::to_string(epoch);return e;
   };
   const auto fixture=std::filesystem::current_path()/("gec_datum_"+std::to_string(trials)+".csv");
   std::ofstream out(fixture);out<<"synthetic CSV65\n"<<std::setprecision(17);
   for(int ep=0;ep<3;++ep)for(int g=0;g<3;++g)for(int sat=1;sat<=4;++sat)for(int band=0;band<2;++band)
   {
    if(scenario==5 && ep==1 && g==0 && sat==2)continue; // Missing intermediate file row, never bridged by numbers.
    ZhangProductDatumTransition t;t.target=endpoint(ep,g,sat,band);t.wavelength=CLIGHT/frequencies[g][band];
    t.hasSource=ep>0;if(ep)t.source=endpoint(ep-1,g,sat,band);
    t.correctionShiftCycles=ep==1?integer(g,sat,band):ep==2 && scenario==12?-integer(g,sat,band):0;t.correctionShiftMetres=t.wavelength*t.correctionShiftCycles;
    t.dualFrequencyBinding=t.target.bundle+"|"+t.target.satellite;t.integerProofIds={"SAME_TIME_PURE_PHASE_INTEGER","FULL_Q_MAP"};
    if(scenario==4 && ep==1 && g==0 && sat==2)t.source.posterior="WRONG_SOURCE";
    ZhangPublishedContinuityContract c;c.version=2;c.bundleId=t.target.bundle;c.jointCovarianceComplete=true;c.proofId=zhangDatumWireProof(t);
    if(ep) {
     c.previousEpochGpstSeconds=t.source.epoch;c.previousCounter=t.source.counter;c.previousDatumVersion=t.source.datum;c.previousPublicSegment=t.source.segment;
     if(!(scenario==6 && ep==1 && g==0 && sat==2))require(zhangAuthorizePublishedAffine(c,t.target.epoch,t.target.segment,t,true,true,true,true,true),"fixture certification");
     if(c.floatTransportValid && !t.correctionShiftCycles)c.kind="IDENTITY";
    }
    // IDENTITY means same target segment in the current wire contract.
    if(c.kind=="IDENTITY")c.kind="CERTIFIED_AFFINE";
    std::vector<std::string> row(52,"0");row[0]=std::to_string(100000+30*ep);row[1]="PRODUCT_FIXED";
    row[2]=t.target.satellite;row[3]=t.target.observable;row[5]=row[7]="0.01";row[10]="0.02";
    const double correction=.35*sat+.4*ep+(ep && !(ep==2 && scenario==12)?integer(g,sat,band)*t.wavelength:0);
    row[6]=std::to_string(-correction);row[9]=std::to_string(correction);row[11]=std::to_string(ep);row[14]="1";row[15]=row[0];row[16]="1";
    row[17]="SYNTHETIC_PROVEN_COORDINATE_CHANGE";row[18]="1";row[19]="CURRENT_ALIGNMENT_VALID";
    for(int f:{20,21,22,23,28,29,30,31,32,42,43,46})row[f]="1";
    row[24]=t.target.component;row[25]="DATUM";row[26]=row[27]=row[0];row[33]="NONE";row[34]="PHYSICAL0";row[35]="1";
    row[38]="4";row[39]=row[40]="3";row[44]="AR_VALID";row[47]=t.target.physicalFunctional;row[50]=c.floatTransportValid?"0":"1";row[51]=t.target.segment;
    for(int f=0;f<52;++f)out<<(f?",":"")<<row[f];out<<',';zhangWritePublishedContinuityCsv(out,c);out<<'\n';
   }
   out.close();acsConfig.zhangPppAr.product_filename=fixture.string();
   ReceiverMap receivers;auto& rec=receivers["GEC_DATUM_TEST"];rec.id="GEC_DATUM_TEST";rec.ready=true;
   for(int g=0;g<3;++g)for(int sat=1;sat<=4;++sat){
    auto satellite=satAt(g,sat);auto& status=rec.satStatMap[satellite];status.el=(80-10*sat)*D2R;
    GObs obs;obs.time.bigTime=100000;obs.Sat=satellite;obs.mount=rec.id;obs.satStat_ptr=&status;
    for(int b=0;b<2;++b){auto freq=code2Freq.at(systems[g]).at(codes[g][b]);auto& tracking=status.sigStatMap[ft2string(freq)].tracking;
     tracking.seen=true;tracking.lastCode=int(codes[g][b]);tracking.lastPhaseTime=100000;
     obs.sigs[freq].code=codes[g][b];obs.sigs[freq].P=2e7;obs.sigs[freq].L=1e8;}
    rec.obsList.push_back(static_cast<std::shared_ptr<GObs>>(obs));
   }
   KFState state;state.time.bigTime=100000;require(bindZhangCheckpointRuntimeId(state,"GEC_DATUM_"+std::to_string(trials)),"initial runtime");
   state.kfIndexMap.clear();state.stateTransitionMap.clear();state.kfIndexMap[state.oneKey]=0;state.stateTransitionMap[state.oneKey][state.oneKey][0]=1;
   state.x=VectorXd::Ones(1);state.dx=VectorXd::Zero(1);state.P=MatrixXd::Zero(1,1);
   std::ostringstream trace;updateZhangPppArUserReferences(trace,receivers,state);
   auto amb=[&](int g,int b,int sat){KFKey k;k.type=KF::AMBIGUITY;k.str=rec.id;k.Sat=satAt(g,sat);k.num=int(codes[g][b]);return k;};
   auto phase=[&](int g,int b){auto k=amb(g,b,1);k.type=KF::PHASE_BIAS;k.Sat=SatSys(systems[g],0);return k;};
   std::vector<KFKey> keys{state.oneKey};
   for(auto type:{KF::REC_POS,KF::TROP,KF::IONO_STEC}){KFKey k;k.type=type;k.str=rec.id;keys.push_back(k);}
   for(int g=0;g<3;++g)for(int b=0;b<2;++b){require(zhangPppArUserReference(state,rec.id,systems[g],codes[g][b])==satAt(g,1),"GEC initial references");if(scenario!=13){keys.push_back(phase(g,b));for(int sat=2;sat<=4;++sat)keys.push_back(amb(g,b,sat));}}
   const int n=keys.size();state.kfIndexMap.clear();state.stateTransitionMap.clear();
   state.x=VectorXd::Zero(n);state.dx=VectorXd::Zero(n);state.x(0)=1;
   MatrixXd factor=MatrixXd::Identity(n,n);for(int i=1;i<n;++i){state.x(i)=.7*i+.11;for(int j=1;j<n;++j)factor(i,j)+=.003*((i*3+j*7)%13);}
   state.P=factor*factor.transpose();state.P.row(0).setZero();state.P.col(0).setZero();
   for(int i=0;i<n;++i){state.kfIndexMap[keys[i]]=i;state.stateTransitionMap[keys[i]][keys[i]][0]=1;}
   persistZhangUserAppliedChartReceipts(state);
   if(scenario==7 || scenario==8){state.metaDataMap.erase(ZHANG_CHECKPOINT_RUNTIME_ID_METADATA);require(bindZhangCheckpointRuntimeId(state,"RESTORE_"+std::to_string(trials)),"runtime owner lost");
    if(scenario==8)std::erase_if(state.metaDataMap,[](const auto& x){return x.first.starts_with("USER_APPLIED_CHART_V1:");});}
   const auto before=state;
   state.time.bigTime=(scenario==3 || scenario==5 || scenario==6 || scenario==12)?100060:100030;
   for(auto& obs:only<GObs>(rec.obsList)){obs.time=state.time;for(auto& [name,sig]:obs.satStat_ptr->sigStatMap)sig.tracking.lastPhaseTime=state.time.bigTime;}
   if(scenario==10){for(auto it=rec.obsList.begin();it!=rec.obsList.end();)if(std::static_pointer_cast<GObs>(*it)->Sat==satAt(0,2))it=rec.obsList.erase(it);else ++it;}
   if(scenario==9)state.exactStateTransformCallback=[](auto&&...){return false;};
   const auto candidateBefore=state;bool threw=false;try{updateZhangPppArUserReferences(trace,receivers,state);}catch(const std::exception& e){threw=true;if(scenario!=9)std::cerr<<e.what()<<'\n'<<trace.str()<<'\n';}
   if(scenario==9){require(threw,"reject transaction callback");require(state.x==candidateBefore.x && state.P==candidateBefore.P && state.metaDataMap==candidateBefore.metaDataMap,"atomic rollback including all bands/systems");continue;}
   require(!threw,"production chart transaction succeeds");
   if(scenario==13){require(state.x==before.x && state.P==before.P && state.kfIndexMap==before.kfIndexMap,"a selected chart without historical phase states acknowledges a new chart without transporting or releasing information");continue;}
   if(scenario==8){for(int g=0;g<3;++g)for(int b=0;b<2;++b)require(!state.kfIndexMap.contains(phase(g,b)),"unowned state fails closed");continue;}
   // Independent matrix oracle includes every cross-system/TROP/IONO/POS block.
   std::map<KFKey,VectorXd> rows;
   for(const auto& k:keys){VectorXd row=VectorXd::Zero(n);row(before.kfIndexMap.at(k))=1;rows[k]=row;}
   for(int g=0;g<3;++g)for(int b=0;b<2;++b){
    const auto ref=zhangPppArUserReference(state,rec.id,systems[g],codes[g][b]);int r=0;for(int sat=1;sat<=4;++sat)if(ref==satAt(g,sat))r=sat;
    require(r>0,"valid final reference");const double l=CLIGHT/frequencies[g][b];
    VectorXd refrow=VectorXd::Zero(n);if(r!=1)refrow(before.kfIndexMap.at(amb(g,b,r)))=1;
    rows[phase(g,b)]+=l*refrow;rows[phase(g,b)](0)+=scenario==12?0:l*integer(g,r,b);
    rows.erase(amb(g,b,r));
    for(int sat=1;sat<=4;++sat)if(sat!=r){VectorXd row=VectorXd::Zero(n);if(sat!=1)row(before.kfIndexMap.at(amb(g,b,sat)))=1;
     row-=refrow;row(0)+=scenario==12?0:integer(g,sat,b)-integer(g,r,b);rows[amb(g,b,sat)]=row;}
    for(int sat=1;sat<=4;++sat)require(zhangUserArcGeneration(state.metaDataMap,zhangUserArcGenerationKey(true,rec.id,satAt(g,sat).id(),int(codes[g][b])))==0,"product policy never invents physical cycle slip");
    const bool unknown=(scenario==4 || scenario==5 || scenario==6) && g==0;
    if(unknown)require(!state.kfIndexMap.contains(amb(g,b,2)),"broken proof chain releases only affected direction");
    if(std::string(policy)=="RESET_KNOWN_LOCAL" && scenario!=2 && mask&(1<<g))for(int sat=1;sat<=4;++sat)
      if(integer(g,sat,b)!=0)require(!state.kfIndexMap.contains(amb(g,b,sat)),"known-history R release is local");
   }
   MatrixXd T(state.x.size(),n);for(const auto& [key,index]:state.kfIndexMap){require(rows.contains(key),"oracle knows actual key");T.row(index)=rows.at(key);}
   require((state.x-T*before.x).cwiseAbs().maxCoeff()<1e-10,"transport/survivor marginal means equal independent full oracle");
   require((state.P-T*before.P*T.transpose()).cwiseAbs().maxCoeff()<1e-11,"complete covariance incl cross-system and released marginal matches oracle");
   const auto once=state;updateZhangPppArUserReferences(trace,receivers,state);
   require(state.x==once.x && state.P==once.P && state.metaDataMap==once.metaDataMap,"repeated packet exactly idempotent");
  }
  std::cout<<"PASS GEC UDUC production datum policy trials="<<trials<<" checks="<<checks<<" full_cross_Q=1 two_frequency=1 endpoint_chain=1 local_release=1 restart_owner=1 rollback=1\n";return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL trial="<<trials<<" scenario="<<active<<" "<<e.what()<<'\n';return 1;}
}
