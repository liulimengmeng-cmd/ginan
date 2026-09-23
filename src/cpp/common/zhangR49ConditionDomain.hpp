#pragma once
#include "common/zhangR47ProductDomain.hpp"
#include "common/zhangIntegerCandidateNis.hpp"
#include <sstream>
enum class ZhangR49ConditionSource {CURRENT_NETWORK,TRANSPORTED_HISTORY,GAUGE,CURRENT_PRODUCT_SEARCH,BRIDGE};
struct ZhangR49ConditionRecord {
 ZhangExactVector row;ZhangExactInteger value;
 ZhangR49ConditionSource source;
 std::vector<std::string> algebraicParentIds,conditioningParentIds;
};
struct ZhangR49ConditionDomain {
 std::string posteriorId,chartId,policyId;
 std::uint64_t version=0;
 std::vector<ZhangR49ConditionRecord> records;
 std::vector<std::string> searchTickets;
 // HNF and physical/product projection remain in the existing finalFrame.
 // This metadata moves atomically with the complete route, not with a graph.
 std::string scope() const {return posteriorId+"|"+chartId+"|"+std::to_string(version)+"|"+policyId;}
};
inline std::string zhangR49CanonicalTarget(const ZhangR47ProductSearchFrame& domain,
 const ZhangExactMatrix& targets,int dimension) {
 if(!domain.valid || !domain.affine || !domain.affine->valid || !zhangExactRectangularMatrix(targets,dimension))return {};
 std::vector<std::string> signatures;
 for(const auto& row:targets) {
  ZhangExactVector reduced(domain.affine->kernelBasis.size());
  for(int k=0;k<reduced.size();++k)for(int c=0;c<domain.columns.size();++c)
   reduced[k]+=row[domain.columns[c]]*domain.affine->kernelBasis[k][c];
  std::vector<bool> present(dimension);for(int c:domain.columns)present[c]=true;
  for(int c=0;c<dimension;++c)if(!present[c])reduced.push_back(row[c]);
  int sign=1;for(const auto& x:reduced)if(x!=0){sign=x<0?-1:1;break;}
  std::ostringstream text;for(const auto& x:reduced)text<<sign*x<<",";
  signatures.push_back(text.str());
 }
 // Only signed permutations and exact integer translations are identified.
 // Never divide coefficients by gcd: 2Z is not Z.
 std::sort(signatures.begin(),signatures.end());std::string key;
 for(const auto& row:signatures)key+="["+row+"]";return key;
}
