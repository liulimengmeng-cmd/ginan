#include "common/zhangCycleSlipEvidence.hpp"
#include "common/zhangProductPhysicalCatalogue.hpp"
#include <set>
#include <iostream>
#include <stdexcept>
static void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
static ZhangSlipSample sample(int epoch,int sat,int slip=0)
{
 const std::array<double,3> l{299792458./1575420000,299792458./1227600000,299792458./1176450000};
 const double range=1000+10*sat,clock=epoch>=60?300000:0,I=3+sat+.001*epoch+.00001*epoch*epoch;
 std::array<double,3> phase;for(int f=0;f<3;f++)phase[f]=range+clock-std::pow(l[f]/l[0],2)*I+(f==0?l[0]*slip:0);
 return zhangMakeSlipSample(epoch*30,123,l,phase,{4e-6,4e-6,4e-6},
  {range+clock+I,range+clock+std::pow(l[1]/l[0],2)*I},{.04,.04},false,true,range,true);
}
static void recovery(bool restore,int injection)
{
 std::array<ZhangSlipEvidenceState,6> states;int rebases=0,usedAfter=0;
 for(int epoch=0;epoch<230;epoch++)
 {
  const auto before=states;std::array<ZhangSlipSample,6> raw;
  for(int sat=0;sat<6;sat++)raw[sat]=sample(epoch,sat,sat==0&&epoch>=160?injection:0);
  for(int sat=0;sat<6;sat++)
  {
   std::vector<ZhangSlipClockPeer> peers;const auto* ref=zhangSlipClockReference(before[sat]);
   for(int p=0;p<6;p++)if(p!=sat && ref)
   {
    if(sat==0 && epoch>=40 && (epoch<160 || !restore) && peers.size()>=2)continue;
    const auto& history=before[p].screening;
    auto at=std::find_if(history.begin(),history.end(),[&](const auto& x){return x.time==ref->time && x.geometry;});
    if(at==history.end())continue;
    peers.push_back({at->time,raw[p].rawIf-at->rawIf,raw[p].covariance[0]+at->covariance[0],history,raw[p],true});
   }
   const auto out=zhangInspectSlipEvidence(states[sat],raw[sat],peers);
   const auto audit=zhangSlipRecoveryAudit;
   if(sat!=0)continue;
   if(epoch>=40 && epoch<160)require(!out.usePhase && !out.breakArc,"two-peer outage manufactured continuity/retirement");
   if(audit.recoverySamples)
   {
    require(restore && epoch>=160,"rebased before independent partners returned");
    require(audit.recoverySamples==8 && out.breakArc && out.retire[0] && out.retire[1] && !out.usePhase,
      "fresh baseline inherited old integer arc or used rebaseline-epoch phase");
    require(out.action==ZhangSlipAction::UNRESOLVED_STEP,"rebaseline masquerades as identified integer step");
    require(states[sat].accepted.front().time>=raw[sat].time-210,"rebaseline imported stale accepted history");
    rebases++;
   }
   if(out.action==ZhangSlipAction::INTEGER_STEP)require(false,"recovery invented a unique physical slip certificate");
   if(epoch>=175)usedAfter+=out.usePhase;
  }
 }
 if(restore)require(rebases==1 && usedAfter>=50,"stable fresh data did not recover once after peer-window loss");
 else require(rebases==0 && usedAfter==0,"insufficient independent partners released quarantine");
}
int main()
{
 struct Graph{bool connected=true;std::string rootReceiver="A";std::set<std::string> receivers{"A","B"};std::set<int> satellites{6,10,12};} graph;
 const auto all=zhangPhysicalProductSatelliteCatalogue(graph);
 require(all==std::set<int>({6,10,12}),"tree-only reference node lost from physical catalogue");
 graph.connected=false;require(zhangPhysicalProductSatelliteCatalogue(graph).empty(),"forest authorized one cross-component datum");
 graph.connected=true;graph.receivers.erase("A");require(zhangPhysicalProductSatelliteCatalogue(graph).empty(),"missing root accepted");
 recovery(true,0);recovery(false,0);recovery(true,5);recovery(true,-3);
 std::cout<<"GPS10_QC_RECOVERY_PASS dual_frequency=1 outage_epochs=120 exact_clock_jump=1 fresh_window=8 both_old_arcs_retired=1 insufficient_peers_stay_quarantined=1\n";
}
