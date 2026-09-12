#pragma once
#include "common/zhangFullRank.hpp"
struct ZhangR48DatumProposal {
 ZhangGraphBasis continued, replacement;
 bool oldHealthy=false, identityTransport=false;
 std::string status="TRANSPORT_UNPROVEN";
};
// Pure proposal: no controller, frontend potential or Ledger is reachable.
inline ZhangR48DatumProposal proposeProductDatum(
 const ZhangGraphBasis& old,const ZhangGraphBasis& replacement,
 const std::set<ZhangGraphEdge>& hardInvalid,bool initialized) {
 ZhangR48DatumProposal p;p.replacement=replacement;
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
 p.status=p.identityTransport?"IDENTITY_TRANSPORT_PROVEN":"TRANSPORT_UNPROVEN";
}
inline ZhangGraphBasis commitProductDatum(const ZhangR48DatumProposal& p) {
 return p.identityTransport?p.continued:p.replacement;
}
