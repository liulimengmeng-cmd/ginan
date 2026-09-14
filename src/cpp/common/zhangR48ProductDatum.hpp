#pragma once
#include "common/zhangFullRank.hpp"
#include "common/zhangR51Integrity.hpp"
struct ZhangR48DatumProposal {
 ZhangGraphBasis continued, replacement;
 bool oldHealthy=false, identityTransport=false, exactBasisTransport=false;
 bool transportProven() const {return identityTransport || exactBasisTransport;}
 bool samePhysicalDomain=false;
 std::string status="TRANSPORT_UNPROVEN";
};
// Pure proposal: no controller, frontend potential or Ledger is reachable.
inline ZhangR48DatumProposal proposeProductDatum(
 const ZhangGraphBasis& old,const ZhangGraphBasis& replacement,
 const std::set<ZhangGraphEdge>& hardInvalid,bool initialized) {
 ZhangR48DatumProposal p;p.replacement=replacement;
 p.samePhysicalDomain=initialized && old.connected && replacement.connected && old.edges==replacement.edges && hardInvalid.empty();
 p.oldHealthy=initialized && old.connected && old.rootReceiver==replacement.rootReceiver;
 for(const auto& e:old.treeEdges)
  p.oldHealthy &= replacement.edges.count(e)>0 && hardInvalid.count(e)==0;
 if(p.oldHealthy)
  p.continued=zhangBuildSpanningTree(replacement.edges,replacement.rootReceiver,old.treeEdges);
 return p;
}
// C2 control policy admits identity transport only. Hard rebuilds receive no
// inherited integer authorization; non-identity transport needs fresh proof.
inline void auditProductDatumTransport(ZhangR48DatumProposal& p,
 const ZhangGraphBasis& old) {
 p.identityTransport=p.oldHealthy && p.continued.connected;
 for(const auto& e:old.treeEdges)p.identityTransport &= p.continued.treeEdges.count(e)>0;
 if(zhangR51Enabled() && p.samePhysicalDomain && old.treeEdges!=p.replacement.treeEdges) {
  // Prove the two coordinate maps by exact physical expansion in both
  // directions. Integer RHS and frontend alignment still require their
  // physical-ledger proofs at consumption; graph connectivity grants no AR.
  auto roundtrip=[](const ZhangGraphBasis& a,const ZhangGraphBasis& b) {
   for(const auto& chord:a.edges) if(!a.treeEdges.contains(chord)) {
    const auto cycle=zhangFundamentalCycle(a,chord);
    std::map<ZhangGraphEdge,ZhangExactInteger> original,rebuilt;
    for(const auto& [edge,x]:cycle) original[edge]+=x;
    for(const auto& [edge,x]:original) if(!b.treeEdges.contains(edge) && x!=0)
     for(const auto& [term,y]:zhangFundamentalCycle(b,edge)) rebuilt[term]+=x*y;
    std::erase_if(original,[](const auto& item){return item.second==0;});
    std::erase_if(rebuilt,[](const auto& item){return item.second==0;});
    if(original!=rebuilt)return false;
   }
   return true;
  };
  p.exactBasisTransport=roundtrip(old,p.replacement) && roundtrip(p.replacement,old);
  if(p.exactBasisTransport){p.identityTransport=false;p.continued=p.replacement;}
 }
 p.status=p.exactBasisTransport?"EXACT_NONIDENTITY_BASIS_ROUNDTRIP":p.identityTransport?"IDENTITY_TRANSPORT_PROVEN":"TRANSPORT_UNPROVEN";
}
inline ZhangGraphBasis commitProductDatum(const ZhangR48DatumProposal& p) {
 return p.transportProven()?p.continued:p.replacement;
}
