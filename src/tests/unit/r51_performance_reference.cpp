#include "common/zhangR51PhysicalImage.hpp"
#include "common/zhangR49PosteriorWorkspace.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>

static void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
template<class F> double seconds(F f) {auto t=std::chrono::steady_clock::now();f();return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();}
static void sameImage(const ZhangR47ProductSearchFrame& a,const ZhangR47ProductSearchFrame& b) {
 require(a.valid==b.valid && a.reason==b.reason,"image validity/reason differs");
 if(!a.valid)return;
 require(a.columns==b.columns && a.searchRank==b.searchRank && a.conditionerRank==b.conditionerRank,"image dimensions/order differs");
 require(a.affine.particularSolution==b.affine.particularSolution && a.affine.kernelBasis==b.affine.kernelBasis && a.imageGenerators==b.imageGenerators,"exact affine/image differs");
}
static MatrixXd condition(ZhangR49PosteriorWorkspace& workspace, const MatrixXd& p,
 const MatrixXd& h,const std::vector<std::string>& order,bool retain) {
 require(workspace.factor(p,order,"epoch"),"prior factor failed");
 MatrixXd c=h*p*h.transpose();
 Eigen::SelfAdjointEigenSolver<MatrixXd> eigen((c+c.transpose())*.5);
 require(eigen.info()==Eigen::Success,"constraint eigensolver");
 MatrixXd b=h*workspace.squareRoot;
 MatrixXd basis=MatrixXd::Zero(p.rows(),h.rows());int cols=0;
 for(int i=0;i<h.rows();++i)if(eigen.eigenvalues()[i]>1e-10)
  basis.col(cols++)=b.transpose()*eigen.eigenvectors().col(i)/std::sqrt(eigen.eigenvalues()[i]);
 MatrixXd root=workspace.squareRoot;
 if(cols)root-=(workspace.squareRoot*basis.leftCols(cols))*basis.leftCols(cols).transpose();
 if(retain){require(workspace.formGram(std::move(root),order,"epoch"),"retained Gram failed");return workspace.covariance;}
 MatrixXd result=root*root.transpose();return (.5*(result+result.transpose())).eval();
}
int main(int argc,char** argv) {
 std::mt19937 rng(510915);
 for(int test=0;test<1200;++test) {
  int n=2+rng()%6,m=1+rng()%n,k=rng()%5;
  ZhangExactMatrix t(m,ZhangExactVector(n)),h(k,ZhangExactVector(n));ZhangExactVector v(k),x(n);
  for(auto& value:x)value=int(rng()%7)-3;
  for(auto& row:t)for(auto& value:row)value=int(rng()%5)-2;
  for(int r=0;r<k;++r)for(int c=0;c<n;++c){h[r][c]=int(rng()%7)-3;v[r]+=h[r][c]*x[c];}
  if(k && test%3==0)++v[0]; // also exercise infeasible integer systems
  auto full=zhangR47CompileProductSearchFrame(t,h,v,n);
  auto image=zhangR47CompileProductSearchFrame(t,h,v,n,ZhangProductFrameWork::IMAGE_GENERATORS_ONLY);
  sameImage(full,image);
 }
 for(int n=4;n<=30;n+=2) {
  MatrixXd l=MatrixXd::Random(n,n),p=l*l.transpose()+MatrixXd::Identity(n,n);
  std::vector<std::string> order;for(int i=0;i<n;++i)order.push_back(std::to_string(i));
  MatrixXd old=p,fast=p;ZhangR49PosteriorWorkspace a,b;
  for(int step=0;step<3;++step) {
   MatrixXd h=MatrixXd::Zero(1,n);h(0,step)=1;
   old=condition(a,old,h,order,false);fast=condition(b,fast,h,order,true);
   require((old-fast).norm()<1e-9*std::max(1.0,old.norm()),"sequential posterior mismatch");
   require((h*fast).norm()<1e-8,"constraint closure");
   Eigen::SelfAdjointEigenSolver<MatrixXd> e(fast);require(e.eigenvalues().minCoeff()>-1e-9,"Gram PSD");
  }
  require(a.decompositions==3 && b.decompositions==1 && b.hits==2,"sequential cache reuse");
  fast(n-1,n-1)+=1e-7;
  require(b.factor(fast,order,"epoch") && b.decompositions==2,"changed covariance must invalidate");
  std::swap(order[0],order[1]);require(b.factor(fast,order,"epoch") && b.decompositions==3,"order must invalidate");
  require(b.factor(fast,order,"next") && b.decompositions==4,"epoch must invalidate");
 }
 ZhangR49PosteriorWorkspace invalid;
 require(!invalid.formGram(MatrixXd::Zero(2,3),{"a","b"},"t"),"reject nonsquare root");
 MatrixXd nan=MatrixXd::Identity(2,2);nan(0,0)=std::numeric_limits<double>::quiet_NaN();
 require(!invalid.formGram(nan,{"a","b"},"t"),"reject NaN root");
 std::cout<<"PASS 1200 exact differential systems; 14 sequential PSD/closure/posterior/cache cases\n"<<std::flush;
 if(argc==3 && std::string(argv[1])=="--posterior") {
  std::ifstream input(argv[2],std::ios::binary);std::string magic,line;
  std::getline(input,magic);std::getline(input,line);
  require(magic=="R49_POSTERIOR_V1_COLUMN_MAJOR","posterior fixture magic");
  int n=std::stoi(line);require(n>6 && n<10000,"posterior fixture dimension");
  VectorXd mean(n);MatrixXd p(n,n);
  input.read(reinterpret_cast<char*>(mean.data()),n*sizeof(double));
  input.read(reinterpret_cast<char*>(p.data()),n*n*sizeof(double));require(bool(input),"posterior fixture truncated");
  std::vector<std::string> order;for(int i=0;i<n;++i)order.push_back(std::to_string(i));
  MatrixXd oldP=p,newP=p;ZhangR49PosteriorWorkspace a,b;
  auto run=[&](MatrixXd& q,ZhangR49PosteriorWorkspace& w,bool retain){for(int j=0;j<3;++j){MatrixXd h=MatrixXd::Zero(2,n);h(0,j*2)=1;h(1,j*2+1)=1;q=condition(w,q,h,order,retain);}};
  double old=seconds([&]{run(oldP,a,false);}),fast=seconds([&]{run(newP,b,true);});
  double relative=(oldP-newP).norm()/std::max(1.0,oldP.norm());require(relative<1e-9,"recorded posterior mismatch");
  std::cout<<"RECORDED_ROOT_BENCH n="<<n<<" old_s="<<old<<" optimized_s="<<fast<<" speedup="<<old/fast<<" relative_error="<<relative<<" synthetic_constraints=1 old_decompositions="<<a.decompositions<<" new_decompositions="<<b.decompositions<<"\n";
  return 0;
 }
 if(argc>1) {
  int m=std::stoi(argv[1]),n=2*m+2,k=std::min(m/4,52);
  ZhangExactMatrix t(m,ZhangExactVector(n)),h(k,ZhangExactVector(n));ZhangExactVector v(k);
  for(int i=0;i<m;++i){t[i][2*i]=1;t[i][2*i+1]=-1;t[i][n-2]=-1;t[i][n-1]=1;}
  for(int i=0;i<k;++i){h[i]=t[i];v[i]=i%13;}
  ZhangR47ProductSearchFrame a,b;
  double old=seconds([&]{a=zhangR47CompileProductSearchFrame(t,h,v,n);});
  double fast=seconds([&]{b=zhangR47CompileProductSearchFrame(t,h,v,n,ZhangProductFrameWork::IMAGE_GENERATORS_ONLY);});
  sameImage(a,b);
  std::cout<<"DOMAIN_BENCH targets="<<m<<" variables="<<n<<" history="<<k<<" old_s="<<old<<" optimized_s="<<fast<<" speedup="<<old/fast<<" exact_equal=1\n"<<std::flush;
  int size=argc>2?std::stoi(argv[2]):500;
  MatrixXd l=MatrixXd::Random(size,size),p=l*l.transpose()+MatrixXd::Identity(size,size);
  std::vector<std::string> order;for(int i=0;i<size;++i)order.push_back(std::to_string(i));
  MatrixXd oldP=p,newP=p;ZhangR49PosteriorWorkspace wa,wb;
  auto run=[&](MatrixXd& q,ZhangR49PosteriorWorkspace& w,bool retain){for(int j=0;j<3;++j){MatrixXd h=MatrixXd::Zero(2,size);h(0,j*2)=1;h(1,j*2+1)=1;q=condition(w,q,h,order,retain);}};
  old=seconds([&]{run(oldP,wa,false);});fast=seconds([&]{run(newP,wb,true);});
  const double rel=(oldP-newP).norm()/oldP.norm();require(rel<1e-9,"benchmark posterior mismatch");
  std::cout<<"ROOT_BENCH n="<<size<<" old_s="<<old<<" optimized_s="<<fast<<" speedup="<<old/fast<<" relative_error="<<rel<<" old_decompositions="<<wa.decompositions<<" new_decompositions="<<wb.decompositions<<"\n";
 }
}
