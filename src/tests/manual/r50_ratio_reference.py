from pathlib import Path
import subprocess,tempfile,sys
root=Path(__file__).resolve().parents[3]
s=(root/'src/cpp/ambres/GNSSambres.cpp').read_text()
a=s.index('static std::multimap<double, VectorXd> lambdaSearchReducedSuffixRatio(')
b=s.index('/** Lambda algorithm',a)
body=s[a:b]
legacy_start=s.index('static map<double, VectorXd> lambdaSearchReducedSuffix(')
body=s[legacy_start:a]+body
harness=r'''#include <Eigen/Dense>
#include <map>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include "common/zhangR50Validation.hpp"
using namespace Eigen;using namespace std;
#define ROUND(x) floor((x)+0.5)
struct GinAR_opt {double ratthr=3;int nset=2;};
struct GinAR_mtx {VectorXd Dtrs,zflt;MatrixXd Ltrs;};
'''+body+r'''
int main(){
 assert(zhangR50AssessRatio(1,3,3).accepted);
 assert(!zhangR50AssessRatio(1,2.99,3).accepted);
 assert(!zhangR50AssessRatio(1,1,3).accepted);
 assert(!zhangR50AssessRatio(0,0,3).accepted);
 assert(zhangR50AssessRatio(0,1,3).accepted);
 assert(!zhangR50AssessRatio(1,numeric_limits<double>::infinity(),3).valid);
 assert(!zhangR50AssessRatio(-1,2,3).valid);
 int incomplete=0;
 for(double x: {0.0,0.1,0.49,0.5,-0.5,1.4}) for(double y:{0.0,0.2,0.5}) {
  GinAR_mtx m;m.Dtrs=Vector2d(0.3,0.8);m.zflt=Vector2d(x,y);m.Ltrs=Matrix2d::Identity();m.Ltrs(1,0)=0.3;
  auto got=lambdaSearchReducedSuffixRatio(m,2,{});assert(got.size()>=2);
  vector<double> costs;
  for(int i=-6;i<=6;++i)for(int j=-6;j<=6;++j){double d1=y-j,d0=x-i-d1*0.3;costs.push_back(d0*d0/0.3+d1*d1/0.8);}
  sort(costs.begin(),costs.end());int k=0;bool exact=true;for(auto& [c,v]:got){if(k==2)break;if(abs(c-costs[k++])>=1e-10)exact=false;}if(!exact)++incomplete;
  // With unique costs, R50 preserves R49's enumeration and candidate-limit exit.
  auto old=lambdaSearchReducedSuffix(m,2,{});
  bool unique=true;double prev=-1;for(auto& [c,v]:got){if(c==prev)unique=false;prev=c;}
  if(unique){assert(old.size()==got.size());auto it=old.begin();for(auto& [c,v]:got){assert(c==it->first);assert(v==it->second);++it;}}
 }
 GinAR_mtx m;m.Dtrs=VectorXd::Ones(1);m.zflt=VectorXd::Constant(1,0.5);m.Ltrs=MatrixXd::Identity(1,1);
 auto ties=lambdaSearchReducedSuffixRatio(m,1,{});assert(ties.size()>=2);assert(ties.begin()->first==next(ties.begin())->first);
 cout<<"R49 unchanged candidate-limit search: non-global-top-two fixtures="<<incomplete<<"/18\n";
 cout<<"PASS R50 ratio boundary cases and 18 R49 candidate-enumeration equivalence fixtures with brute-force diagnostics; equal-cost candidates retained\n";
}
'''
cap=int(sys.argv[1]) if len(sys.argv)>1 else 2
assert cap>=2
harness=harness.replace('int nset=2;',f'int nset={cap};')
print(f'candidate_set_size={cap}',flush=True)
with tempfile.TemporaryDirectory(prefix='r50_ratio_') as d:
 p=Path(d);(p/'test.cpp').write_text(harness)
 subprocess.run(['g++','-std=c++20','-O1','-I'+str(root/'src/cpp'),'-I/home/rx/GINAN/vcpkg_installed/linux/x64-linux/include/eigen3',str(p/'test.cpp'),'-o',str(p/'test')],check=True)
 subprocess.run([str(p/'test')],check=True)
