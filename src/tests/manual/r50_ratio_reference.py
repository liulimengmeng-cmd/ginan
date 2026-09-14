from pathlib import Path
import subprocess,tempfile
root=Path(__file__).resolve().parents[3]
s=(root/'src/cpp/ambres/GNSSambres.cpp').read_text()
a=s.index('static std::multimap<double, VectorXd> lambdaSearchReducedSuffixRatio(')
b=s.index('/** Lambda algorithm',a)
body=s[a:b]
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
 for(double x: {0.0,0.1,0.49,0.5,-0.5,1.4}) for(double y:{0.0,0.2,0.5}) {
  GinAR_mtx m;m.Dtrs=Vector2d(0.3,0.8);m.zflt=Vector2d(x,y);m.Ltrs=Matrix2d::Identity();m.Ltrs(1,0)=0.3;
  auto got=lambdaSearchReducedSuffixRatio(m,2,{});assert(got.size()==2);
  vector<double> costs;
  for(int i=-6;i<=6;++i)for(int j=-6;j<=6;++j){double d1=y-j,d0=x-i-d1*0.3;costs.push_back(d0*d0/0.3+d1*d1/0.8);}
  sort(costs.begin(),costs.end());int k=0;for(auto& [c,v]:got){assert(abs(c-costs[k++])<1e-10);}
 }
 GinAR_mtx m;m.Dtrs=VectorXd::Ones(1);m.zflt=VectorXd::Constant(1,0.5);m.Ltrs=MatrixXd::Identity(1,1);
 auto ties=lambdaSearchReducedSuffixRatio(m,1,{});assert(ties.size()==2);assert(ties.begin()->first==next(ties.begin())->first);
 cout<<"PASS R50 ratio boundary cases and 18 top-two brute-force comparisons; equal-cost candidates retained\n";
}
'''
with tempfile.TemporaryDirectory(prefix='r50_ratio_') as d:
 p=Path(d);(p/'test.cpp').write_text(harness)
 subprocess.run(['g++','-std=c++20','-O1','-I'+str(root/'src/cpp'),'-I/home/rx/GINAN/vcpkg_installed/linux/x64-linux/include/eigen3',str(p/'test.cpp'),'-o',str(p/'test')],check=True)
 subprocess.run([str(p/'test')],check=True)
