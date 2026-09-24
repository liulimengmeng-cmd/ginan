#include "common/zhangR51PhysicalImage.hpp"
#include "common/zhangR49PosteriorWorkspace.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>

static void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
template<class F> double seconds(F f) {auto t=std::chrono::steady_clock::now();f();return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();}
// Deliberately retain the old all-rows kernel construction as an independent
// oracle for the sparse touched-row projection in zhangIntegerAudit.hpp.
static ZhangExactSurvivingLattice referenceSurvivingLattice(
 const ZhangExactMatrix& rows,const ZhangExactVector& values,
 const std::vector<bool>& mask,bool track) {
 ZhangExactSurvivingLattice result;
 if(rows.size()!=values.size()) {result.consistent=false;return result;}
 std::vector<std::size_t> kept,removed;
 for(std::size_t c=0;c<mask.size();++c)(mask[c]?kept:removed).push_back(c);
 for(const auto& row:rows) {
  if(row.size()!=mask.size()) {result.consistent=false;return result;}
  for(auto c:removed)if(row[c]!=0){++result.touchedRows;break;}
 }
 ZhangExactMatrix transpose(removed.size(),ZhangExactVector(rows.size()));
 for(std::size_t i=0;i<removed.size();++i)
  for(std::size_t j=0;j<rows.size();++j)transpose[i][j]=rows[j][removed[i]];
 const auto combinations=zhangExactIntegerKernel(std::move(transpose),rows.size());
 result.combinationRank=combinations.size();
 ZhangExactMatrix projected;
 ZhangExactVector rhs;
 for(const auto& combination:combinations) {
  ZhangExactVector row(kept.size());ZhangExactInteger value=0;
  for(std::size_t j=0;j<rows.size();++j) {
   value+=combination[j]*values[j];
   for(std::size_t c=0;c<kept.size();++c)
    row[c]+=combination[j]*rows[j][kept[c]];
  }
  projected.push_back(std::move(row));rhs.push_back(value);
 }
 auto hnf=zhangExactRowHermiteNormalForm(std::move(projected),std::move(rhs),track);
 result.consistent=hnf.consistent;
 result.basis=std::move(hnf.basis);result.values=std::move(hnf.values);
 if(track) {
  result.rowTransform=zhangExactZeroMatrix(result.basis.size(),rows.size());
  for(std::size_t i=0;i<hnf.rowTransform.size();++i)
   for(std::size_t k=0;k<combinations.size();++k)
    for(std::size_t j=0;j<rows.size();++j)
     result.rowTransform[i][j]+=hnf.rowTransform[i][k]*combinations[k][j];
 }
 return result;
}
static void checkSurvivingTransform(const ZhangExactSurvivingLattice& result,
 const ZhangExactMatrix& rows,const ZhangExactVector& values,
 const std::vector<bool>& mask) {
 if(!result.consistent)return;
 require(result.rowTransform.size()==result.basis.size(),"projection transform rank");
 std::vector<std::size_t> kept;
 for(std::size_t c=0;c<mask.size();++c)if(mask[c])kept.push_back(c);
 for(std::size_t i=0;i<result.basis.size();++i) {
  require(result.rowTransform[i].size()==rows.size(),"projection transform width");
  ZhangExactVector reconstructed(mask.size());ZhangExactInteger value=0;
  for(std::size_t j=0;j<rows.size();++j) {
   value+=result.rowTransform[i][j]*values[j];
   for(std::size_t c=0;c<mask.size();++c)
    reconstructed[c]+=result.rowTransform[i][j]*rows[j][c];
  }
  require(value==result.values[i],"projection transform RHS");
  for(std::size_t c=0;c<mask.size();++c)
   if(!mask[c])require(reconstructed[c]==0,"projection removed coefficient");
  for(std::size_t c=0;c<kept.size();++c)
   require(reconstructed[kept[c]]==result.basis[i][c],"projection transform coefficient");
 }
}
static void testSparseSurvivingLattice() {
 std::mt19937 rng(510925);
 for(int trial=0;trial<500;++trial) {
  const int n=2+rng()%8,m=1+rng()%10;
  ZhangExactMatrix rows(m,ZhangExactVector(n));ZhangExactVector values(m);
  std::vector<bool> mask(n);
  for(int c=0;c<n;++c)mask[c]=(rng()%3!=0);
  for(int r=0;r<m;++r) {
   for(int c=0;c<n;++c)rows[r][c]=int(rng()%7)-3;
   values[r]=int(rng()%11)-5;
  }
  if(trial%5==0 && m>1) {rows[1]=rows[0];values[1]=values[0]+(trial%2);}
  if(trial%3==0)std::fill(mask.begin(),mask.end(),true);
  const bool track=trial%2==0;
  const auto expected=referenceSurvivingLattice(rows,values,mask,track);
  const auto actual=zhangExactSurvivingLattice(rows,values,mask,track);
  if(!(expected.consistent==actual.consistent &&
   (!expected.consistent || (expected.basis==actual.basis && expected.values==actual.values)) &&
   expected.touchedRows==actual.touchedRows &&
   expected.combinationRank==actual.combinationRank)) {
   std::cerr<<"projection mismatch trial="<<trial<<" n="<<n<<" m="<<m
    <<" track="<<track<<" expected_consistent="<<expected.consistent
    <<" actual_consistent="<<actual.consistent
    <<" expected_rank="<<expected.basis.size()
    <<" actual_rank="<<actual.basis.size()
    <<" expected_combinations="<<expected.combinationRank
    <<" actual_combinations="<<actual.combinationRank<<"\n";
   for(std::size_t r=0;r<rows.size();++r){std::cerr<<"input "<<r<<":";
    for(const auto& v:rows[r])std::cerr<<" "<<v;
    std::cerr<<" = "<<values[r]<<"\n";}
   for(std::size_t r=0;r<expected.basis.size();++r){std::cerr<<"expected "<<r<<":";
    for(const auto& v:expected.basis[r])std::cerr<<" "<<v;
    std::cerr<<" = "<<expected.values[r]<<"\n";}
   for(std::size_t r=0;r<actual.basis.size();++r){std::cerr<<"actual "<<r<<":";
    for(const auto& v:actual.basis[r])std::cerr<<" "<<v;
    std::cerr<<" = "<<actual.values[r]<<"\n";}
   throw std::runtime_error("sparse projection differs from all-rows integer kernel reference");
  }
  if(track)checkSurvivingTransform(actual,rows,values,mask);
 }
 const ZhangExactMatrix rows={{2,0,1,0},{0,1,1,0},{0,0,0,1}};
 const ZhangExactVector values={1,3,0};
 const std::vector<bool> mask={true,true,false,true};
 const auto parity=zhangExactSurvivingLattice(rows,values,mask,true);
 const auto parityReference=referenceSurvivingLattice(rows,values,mask,true);
 require(parity.basis==parityReference.basis && parity.values==parityReference.values,
  "nonprimitive parity relation changed");
 checkSurvivingTransform(parity,rows,values,mask);
 std::cout<<"PASS 500 independent all-rows integer projection comparisons and provenance audits\n";
}
static void sameImage(const ZhangR47ProductSearchFrame& a,const ZhangR47ProductSearchFrame& b) {
 require(a.valid==b.valid && a.reason==b.reason,"image validity/reason differs");
 if(!a.valid)return;
 require(a.columns==b.columns && a.searchRank==b.searchRank && a.conditionerRank==b.conditionerRank,"image dimensions/order differs");
 require(a.affine && b.affine &&
    a.affine->particularSolution==b.affine->particularSolution &&
    a.imageGenerators==b.imageGenerators,
    "exact particular/image differs");
 require(b.affine->kernelBasis.empty(),"generator-only path built a physical kernel");
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
 testSparseSurvivingLattice();
 if(argc>=2 && std::string(argv[1])=="--projection") {
  const int m=argc>=3?std::stoi(argv[2]):400;
  const int n=m/2+72;
  ZhangExactMatrix rows(m,ZhangExactVector(n));ZhangExactVector rhs(m);
  std::vector<bool> mask(n,true);
  for(int c=n-72;c<n;++c)mask[c]=false;
  for(int r=0;r<m;++r){rows[r][r%(n-72)]=1;rhs[r]=r%(n-72);}
  rows[0][n-1]=1;rows[1][n-1]=1;
  ZhangExactSurvivingLattice old,fast;
  const auto oldSeconds=seconds([&]{old=referenceSurvivingLattice(rows,rhs,mask,false);});
  const auto fastSeconds=seconds([&]{fast=zhangExactSurvivingLattice(rows,rhs,mask,false);});
  require(old.consistent==fast.consistent && old.basis==fast.basis &&
   old.values==fast.values,"sparse benchmark exact mismatch");
  std::cout<<"SPARSE_PROJECTION_BENCH rows="<<m<<" removed_columns=72 touched_rows=2"
   <<" old_s="<<oldSeconds<<" optimized_s="<<fastSeconds
   <<" speedup="<<oldSeconds/fastSeconds<<" exact_equal=1\n";
  return 0;
 }
 auto sharedDomain=zhangR47CompileProductSearchFrame(
    {{1,1,1}},{{1,0,0}},{2},3);
 auto sharedTarget=zhangR49CompileTargetOnDomain(
    {{0,1,0}},sharedDomain,3);
 require(sharedDomain.valid && sharedTarget.valid &&
    sharedTarget.affine.get()==sharedDomain.affine.get(),
    "target frame copied the affine domain");
 auto differentDomain=zhangR47CompileProductSearchFrame(
    {{1,1,1}},{{0,1,0}},{2},3);
 require(differentDomain.valid &&
    differentDomain.affine.get()!=sharedDomain.affine.get(),
    "different WL/L1 histories shared an affine domain");
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
  if(full.valid && image.valid && full.affine &&
     !full.affine->kernelBasis.empty()) {
   ZhangExactMatrix compact(t.size(),ZhangExactVector(full.columns.size()));
   for(int r=0;r<t.size();++r)for(int c=0;c<full.columns.size();++c)
    compact[r][c]=t[r][full.columns[c]];
   ZhangExactMatrix k(full.columns.size(),ZhangExactVector(full.affine->kernelBasis.size()));
   for(int r=0;r<full.affine->kernelBasis.size();++r)
    for(int c=0;c<full.columns.size();++c)k[c][r]=full.affine->kernelBasis[r][c];
   require(zhangExactProjectedKernelColumns(full.affine->deterministicBasis,
      compact,full.columns.size())==zhangExactMultiply(compact,k),
      "projected free columns differ from full integer kernel");
  }
  auto fullPhysical=zhangR51PhysicalImage(t,h,v,t,n,
      ZhangProductFrameWork::INTEGER_PROJECTOR);
  auto projectedPhysical=zhangR51PhysicalImage(t,h,v,t,n);
  require(fullPhysical.valid==projectedPhysical.valid &&
      fullPhysical.reason==projectedPhysical.reason,
      "physical image validity differs");
  if(fullPhysical.valid)
   require(fullPhysical.particularTarget==projectedPhysical.particularTarget &&
      fullPhysical.generators==projectedPhysical.generators &&
      fullPhysical.projector==projectedPhysical.projector &&
      fullPhysical.offsets==projectedPhysical.offsets,
      "physical image coordinates differ");
 }
 auto parity=zhangR51PhysicalImage({{2}},{},{},{{2}},1);
 require(parity.valid && parity.generators==ZhangExactMatrix{{2}},
    "nonprimitive 2Z image was saturated");
 auto hiddenParity=zhangR51PhysicalImage({{0,1}},{{2,1}},{1},{{0,1}},2);
 auto hiddenParityFull=zhangR51PhysicalImage({{0,1}},{{2,1}},{1},{{0,1}},2,
     ZhangProductFrameWork::INTEGER_PROJECTOR);
 require(hiddenParity.valid && hiddenParityFull.valid &&
    hiddenParity.generators==ZhangExactMatrix{{2}} &&
    hiddenParity.particularTarget==hiddenParityFull.particularTarget &&
    hiddenParity.projector==hiddenParityFull.projector &&
    hiddenParity.offsets==hiddenParityFull.offsets &&
    hiddenParity.particularTarget[0]%2!=0,
    "hidden integer parity/coset changed under projection");
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
