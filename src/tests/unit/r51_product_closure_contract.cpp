#include "common/zhangR47ProductDomain.hpp"
#include "common/zhangR48Marginal.hpp"
#include "common/zhangR48ValidationCache.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <sys/resource.h>

static void require(bool value,const char* reason){if(!value)throw std::runtime_error(reason);}
static void matrix(std::ostream& out,const ZhangExactMatrix& m){
    out<<m.size()<<" [";for(const auto& row:m){out<<"[";for(const auto& v:row)out<<v<<",";out<<"]";}out<<"]\n";
}
static void dumpVector(std::ostream& out,const ZhangExactVector& v){matrix(out,{v});}
static void snapshot(std::ostream& out,const ZhangR47ProductSearchFrame& f){
    out<<f.valid<<" "<<f.reason<<" "<<f.targetRank<<" "<<f.conditionerRank<<" "<<f.searchRank<<"\n";
    for(auto c:f.columns)out<<c<<",";out<<"\n";
    matrix(out,f.projector);matrix(out,f.imageGenerators);dumpVector(out,f.offsets);
    if(f.affine){
        const auto& a=*f.affine;out<<a.valid<<" "<<a.failureReason<<" "<<a.deterministicRank<<" "<<a.quotientRank<<"\n";
        matrix(out,a.deterministicBasis);dumpVector(out,a.deterministicValues);
        dumpVector(out,a.particularSolution);matrix(out,a.kernelBasis);dumpVector(out,a.smithInvariants);
    }
}
static ZhangExactMatrix transpose(const ZhangExactMatrix& m,int columns){
    ZhangExactMatrix t(columns,ZhangExactVector(m.size()));
    for(int r=0;r<m.size();++r)for(int c=0;c<columns;++c)t[c][r]=m[r][c];return t;
}
static void exactContracts(std::ostream& out){
    // Analytic oracle, independent of HNF/SNF/kernel functions: N0=z,
    // N1=1-2*z, N2=2. The product image is exactly 1+2Z, NOT Z.
    const ZhangExactMatrix h={{2,1,0},{0,0,3}},t={{0,1,0}};
    const auto f=zhangR47CompileProductSearchFrame(t,h,{1,6},3);
    require(f.valid && f.searchRank==1,"parity affine image invalid");
    const auto origin=zhangExactMatrixTimesColumn(t,f.affine->particularSolution);
    require(origin[0]%2!=0,"odd affine offset lost");
    require(f.imageGenerators.size()==1 && f.imageGenerators[0].size()==1 &&
        zhangExactAbs(f.imageGenerators[0][0])==2,"2Z incorrectly saturated");
    for(int z=-15;z<=15;++z){
        ZhangExactVector n={z,1-2*z,2};
        require(zhangExactMatrixTimesColumn(h,n)==ZhangExactVector({1,6}),"analytic oracle");
        const auto q=zhangExactMatrixTimesColumn(f.projector,n);
        require(f.imageGenerators[0][0]*(q[0]+f.offsets[0])==
            zhangExactMatrixTimesColumn(t,n)[0]-origin[0],"affine coordinate reconstruction");
    }
    require(!zhangR47CompileProductSearchFrame({{1,0}},{{2,0},{0,3}},{1,6},2).valid,"2N=1 accepted");
    require(!zhangR47CompileProductSearchFrame({{1,0}},{{0,0}},{1},2).valid,"zero contradiction accepted");
    require(!zhangR47CompileProductSearchFrame({{1,0}},{{1,0},{1,0}},{1,2},2).valid,"duplicate contradiction accepted");
    std::mt19937 gen(20261019);
    for(int trial=0;trial<220;++trial){
        const int n=2+gen()%7,r=gen()%7;
        ZhangExactMatrix a(r,ZhangExactVector(n)),targets(1+gen()%4,ZhangExactVector(n));
        ZhangExactVector n0(n);for(auto& v:n0)v=int(gen()%9)-4;
        for(auto& row:a)for(auto& v:row)v=int(gen()%7)-3;
        if(r && trial%4==0)a.back()=a.front();
        if(r && trial%11==0)a.back().assign(n,0);
        if(r && trial%13==0)a[0][0]=ZhangExactInteger(1)<<130;
        for(auto& row:targets)for(auto& v:row)v=int(gen()%5)-2;
        auto b=zhangExactMatrixTimesColumn(a,n0);
        if(r && trial%5==0)++b.back();
        for(auto work:{ZhangProductFrameWork::INTEGER_PROJECTOR,ZhangProductFrameWork::IMAGE_GENERATORS_ONLY}){
            const auto frame=zhangR47CompileProductSearchFrame(targets,a,b,n,work);snapshot(out,frame);
            if(!frame.valid)continue;
            const auto& domain=*frame.affine;
            ZhangExactVector fullOrigin(n);
            for(int c=0;c<frame.columns.size();++c)fullOrigin[frame.columns[c]]=domain.particularSolution[c];
            require(zhangExactMatrixTimesColumn(a,fullOrigin)==b,"invalid integer particular solution");
            for(const auto& k:domain.kernelBasis){
                ZhangExactVector full(n);for(int c=0;c<frame.columns.size();++c)full[frame.columns[c]]=k[c];
                require(zhangExactMatrixTimesColumn(a,full)==ZhangExactVector(a.size()),"invalid kernel direction");
            }
            if(!domain.kernelBasis.empty()){
                ZhangExactMatrix compact(targets.size(),ZhangExactVector(frame.columns.size()));
                for(int r=0;r<targets.size();++r)for(int c=0;c<frame.columns.size();++c)compact[r][c]=targets[r][frame.columns[c]];
                const auto fullImage=zhangExactMultiply(compact,transpose(domain.kernelBasis,frame.columns.size()));
                const auto imageOnly=zhangExactProjectedKernelColumns(domain.deterministicBasis,compact,frame.columns.size());
                require(fullImage==imageOnly,"projected kernel changed coordinates/order");
            }
        }
        // Membership and combination invariants are included in old/new byte
        // comparison, not just the final floating product values.
        const auto members=zhangIntegerRowLatticeContainsBatch(a,targets);
        for(int i=0;i<members.size();++i){
            const auto single=zhangIntegerRowLatticeContains(a,targets[i]);
            const auto tape=zhangIntegerRowLatticeContains(a,targets[i],true,true);
            require(single.contained==members[i].contained && single.rank==members[i].rank &&
                single.smithInvariants==members[i].smithInvariants && single.combination==members[i].combination &&
                tape.combination==single.combination,"single/batch/tape mismatch");
            out<<single.contained<<" "<<single.rank<<"\n";dumpVector(out,single.smithInvariants);dumpVector(out,single.combination);
        }
    }
}
static void numericContracts(){
    ZhangArThreading::ScopedThreads policy(true);
    constexpr int n=140,m=83;
    Eigen::MatrixXd a=Eigen::MatrixXd::NullaryExpr(n,n,[](auto r,auto c){return std::sin(.13*r+.23*c);});
    Eigen::MatrixXd p=a*a.transpose()+Eigen::MatrixXd::Identity(n,n);
    Eigen::VectorXd mean=Eigen::VectorXd::LinSpaced(n,-2,3);
    Eigen::MatrixXd h=Eigen::MatrixXd::Zero(m,n);for(int i=0;i<m;++i){h(i,i)=1;h(i,n-1)=.1;}
    Eigen::MatrixXd j=Eigen::MatrixXd::Random(12,n);
    Eigen::VectorXd v=h*mean+Eigen::VectorXd::Constant(m,.02);
    ZhangR48MarginalWorkspace w(mean,p);require(w.valid,"workspace factor failed");
    auto got=w.project(j,h,v);require(got.valid,"marginal rejected");
    // Independent dense LDLT conditional-Gaussian oracle.
    const Eigen::MatrixXd s=h*p*h.transpose(),cross=j*p*h.transpose();
    const Eigen::VectorXd expectedMean=j*mean+cross*s.ldlt().solve(v-h*mean);
    const Eigen::MatrixXd expectedCov=j*p*j.transpose()-cross*s.ldlt().solve(cross.transpose());
    require((got.mean-expectedMean).norm()/std::max(1.,expectedMean.norm())<1e-10,"full marginal mean mismatch");
    require((got.covariance-expectedCov).norm()/expectedCov.norm()<1e-10,"full marginal covariance mismatch");
    auto same=w.project(j,h,v);require(w.decompositions==1 && w.targetHits==1,"target cache missed");
    auto v2=v;v2(0)+=.1;auto changed=w.project(j,h,v2);
    require(changed.valid && w.decompositions==1 && w.hits==1 && (changed.mean-same.mean).norm()>1e-5,"RHS was cached");
    // Duplicate conditioner with consistent/inconsistent deterministic RHS.
    Eigen::MatrixXd hd(m+1,n);hd.topRows(m)=h;hd.bottomRows(1)=h.topRows(1);
    Eigen::VectorXd vd(m+1);vd.head(m)=v;vd(m)=v(0);
    require(w.project(j,hd,vd).valid,"consistent redundant conditioner rejected");
    vd(m)+=1;require(!w.project(j,hd,vd).valid,"deterministic contradiction missed");
#ifndef CLOSURE_BASELINE
    ZhangR48ValidationCache validation(mean,p);
    require(validation.assess(h,v,1e-6).valid && validation.assess(h,v,1e-6).valid && validation.hits==1,"validation reuse failed");
    validation.assess(h,v2,1e-6);require(validation.evaluations==2,"changed RHS reused validation");
    validation.assess(h,v2,1e-5);require(validation.evaluations==3,"changed alpha reused validation");
    auto hx=h;hx(0,n-2)=.1;validation.assess(hx,v2,1e-5);
    require(validation.evaluations==4,"same size changed rows reused validation");
    require(std::isnan(validation.assess(h,v,1e-6).nis),"NIS numeric calculation enabled");
#endif
}
static void benchmark(int n,std::ostream& out){
    const int rank=n-54;
    ZhangExactMatrix h(rank,ZhangExactVector(n)),t(44,ZhangExactVector(n));
    ZhangExactVector b(rank);
    for(int r=0;r<rank;++r){h[r][r]=1;b[r]=r%23-11;if(r%181==0)h[r][n-1]=2;}
    for(int r=0;r<44;++r){t[r][r]=1;t[r][rank+r]=2;}
    const auto start=std::chrono::steady_clock::now();
    const auto f=zhangR47CompileProductSearchFrame(t,h,b,n);require(f.valid,"benchmark frame failed");
    const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    snapshot(out,f);struct rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<"frame_seconds="<<seconds<<" dimension="<<n<<" rank="<<rank<<" peak_rss_kib="<<usage.ru_maxrss<<"\n";
    for(const auto& [phase,x]:zhangR48ExactMetrics)std::cout<<phase<<"="<<x.seconds<<" calls="<<x.calls<<"\n";
}
static void replay(const std::string& path,std::ostream& out){
    std::ifstream input(path);std::string magic;int n,work;input>>magic>>n>>work;
    require(magic=="R51_PRODUCT_FRAME_V1" && n>0 && (work==0||work==1),"replay format");
    auto readMatrix=[&](){
        int count;input>>count;ZhangExactMatrix m(count,ZhangExactVector(n));
        for(auto& row:m){int nnz;input>>nnz;for(int i=0;i<nnz;++i){int c;ZhangExactInteger v;input>>c>>v;require(c>=0 && c<n,"replay column");row[c]=v;}}return m;
    };
    const auto t=readMatrix(),h=readMatrix();int count;input>>count;ZhangExactVector b(count);for(auto& v:b)input>>v;
    require(bool(input),"replay truncated");
    const auto start=std::chrono::steady_clock::now();
    const auto f=zhangR47CompileProductSearchFrame(t,h,b,n,static_cast<ZhangProductFrameWork>(work));
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    snapshot(out,f);struct rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<"frame_seconds="<<seconds<<" dimension="<<n<<" input_rows="<<h.size()<<" peak_rss_kib="<<usage.ru_maxrss<<" valid="<<f.valid<<"\n";
    for(const auto& [phase,x]:zhangR48ExactMetrics)std::cout<<phase<<"="<<x.seconds<<" calls="<<x.calls<<"\n";
}
static void marginalBenchmark(int n,std::ostream& out){
    ZhangArThreading::ScopedThreads policy(true);
    const int m=n-117;
    Eigen::MatrixXd p=Eigen::MatrixXd::Identity(n,n);
    Eigen::VectorXd common=Eigen::VectorXd::LinSpaced(n,.01,.02);p.noalias()+=common*common.transpose();
    Eigen::VectorXd mean=Eigen::VectorXd::LinSpaced(n,-2,3);
    Eigen::MatrixXd h=Eigen::MatrixXd::Zero(m,n),j=Eigen::MatrixXd::Zero(44,n);
    for(int i=0;i<m;++i){h(i,i)=1;h(i,n-1)=.1;}
    for(int i=0;i<44;++i){j(i,i)=1;j(i,n-2-i)=2;}
    Eigen::VectorXd v=h*mean+Eigen::VectorXd::Constant(m,.02);
    ZhangR48MarginalWorkspace w(mean,p);require(w.valid,"marginal benchmark root");
    const auto start=std::chrono::steady_clock::now();
    const auto f=w.project(j,h,v);require(f.valid,"marginal benchmark rejected");
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    out<<std::setprecision(17)<<f.mean.transpose()<<"\n"<<f.covariance<<"\n";
    struct rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<"marginal_seconds="<<seconds<<" dimension="<<n<<" rows="<<m<<" peak_rss_kib="<<usage.ru_maxrss<<"\n";
}
int main(int argc,char** argv){try{
    std::ofstream out(argc>2?argv[2]:"/dev/null");
    if(argc>1 && std::string(argv[1])=="bench")benchmark(argc>3?std::stoi(argv[3]):800,out);
    else if(argc>3 && std::string(argv[1])=="replay")replay(argv[3],out);
    else if(argc>1 && std::string(argv[1])=="marginal")marginalBenchmark(argc>3?std::stoi(argv[3]):800,out);
    else {exactContracts(out);numericContracts();std::cout<<"closure contracts passed\n";}
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<"\n";return 1;}}
