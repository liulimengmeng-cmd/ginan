#include "common/zhangIntegerConditioner.hpp"
#ifndef CONDITION_BASELINE
#include "common/zhangIncrementalIntegerConditioner.hpp"
#endif
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

static void require(bool pass,const char* why)
{ if(!pass)throw std::runtime_error(why); }

struct Runner
{
#ifdef CONDITION_BASELINE
    VectorXd m; MatrixXd p;
    std::size_t calls=0;
    ZhangPosteriorEffectiveConditioningResult trial;
    Runner(VectorXd mean,MatrixXd covariance):m(std::move(mean)),p(std::move(covariance)){}
    const auto& mean()const{return m;} const auto& covariance()const{return p;}
    const auto& prepare(const MatrixXd& a,const VectorXd& z)
    {++calls;trial=zhangConditionPosteriorEffectiveIntegers(m,p,a,z);return trial;}
    void commit(){require(trial.valid,"invalid baseline commit");m=std::move(trial.mean);p=std::move(trial.covariance);}
    std::size_t builds()const{return calls;}
#else
    ZhangIncrementalIntegerConditioner workspace;
    ZhangIncrementalIntegerConditioner::Trial trial;
    Runner(VectorXd mean,MatrixXd covariance):workspace(std::move(mean),std::move(covariance)){}
    const auto& mean()const{return workspace.mean();} const auto& covariance()const{return workspace.covariance();}
    const auto& prepare(const MatrixXd& a,const VectorXd& z)
    {trial=workspace.prepare(a,z);return trial.result;}
    void commit(){require(workspace.commit(std::move(trial)),"stale conditioner commit");}
    std::size_t builds()const{return workspace.decompositions;}
#endif
};
template<class T> static void read(std::istream& s,T& value)
{s.read(reinterpret_cast<char*>(&value),sizeof(value));}
template<class T> static void write(std::ostream& s,const T& value)
{s.write(reinterpret_cast<const char*>(&value),sizeof(value));}
template<class T> static void dense(std::ostream& s,const T& value)
{s.write(reinterpret_cast<const char*>(value.data()),value.size()*sizeof(double));}

static void emit(std::ostream& out,const ZhangPosteriorEffectiveConditioningResult& result)
{
    int valid=result.valid,rank=result.conditionerRank,det=result.deterministicRows;
    int length=result.failureReason.size();
    write(out,valid);write(out,rank);write(out,det);write(out,length);
    out.write(result.failureReason.data(),length);
    if(result.valid){dense(out,result.mean);dense(out,result.covariance);}
}

static void replay(const std::string& input,const std::string& output)
{
    std::ifstream in(input,std::ios::binary);in.exceptions(std::ios::badbit|std::ios::failbit);
    int n;read(in,n);require(n>0&&n<20000,"invalid replay dimension");
    VectorXd mean(n);MatrixXd covariance(n,n);
    in.read(reinterpret_cast<char*>(mean.data()),n*sizeof(double));
    in.read(reinterpret_cast<char*>(covariance.data()),n*n*sizeof(double));
    Runner runner(std::move(mean),std::move(covariance));
    std::ofstream out(output,std::ios::binary);out.exceptions(std::ios::badbit|std::ios::failbit);write(out,n);
    int calls=0,commits=0;
    const auto started=std::chrono::steady_clock::now();
    while(in.peek()!=EOF)
    {
        char event;read(in,event);
        if(event=='P')
        {
            int count;read(in,count);require(count>0&&count<=256,"invalid replay rows");
            MatrixXd rows(count,n);VectorXd z(count);
            in.read(reinterpret_cast<char*>(rows.data()),rows.size()*sizeof(double));
            in.read(reinterpret_cast<char*>(z.data()),z.size()*sizeof(double));
            emit(out,runner.prepare(rows,z));++calls;
        }
        else if(event=='C'){runner.commit();++commits;}
        else throw std::runtime_error("invalid replay event");
    }
    std::cout<<"REPLAY n="<<n<<" calls="<<calls<<" commits="<<commits<<" builds="<<runner.builds()
        <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<'\n';
}

static void contracts(int n,int steps,int block)
{
    MatrixXd l(n,11);VectorXd mean(n);
    for(int i=0;i<n;++i){mean(i)=.02*std::sin(.13*i);for(int j=0;j<11;++j)l(i,j)=.03*std::sin(.11*i+.17*j);}
    const MatrixXd prior=l*l.transpose()+MatrixXd::Identity(n,n);
    Runner runner(mean,prior);
    const auto started=std::chrono::steady_clock::now();
    for(int step=0;step<steps;++step)
    {
        MatrixXd a=MatrixXd::Zero(block,n);a.middleCols(step*block,block).setIdentity();
        const auto& result=runner.prepare(a,VectorXd::Zero(block));
        require(result.valid&&result.conditionerRank==block,"sequential rank mismatch");runner.commit();
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    MatrixXd a=MatrixXd::Zero(steps*block,n);a.leftCols(steps*block).setIdentity();
    const MatrixXd cross=prior*a.transpose();const MatrixXd s=a*cross;
    // Independent one-shot Schur complement, not the sequential core.
    const VectorXd expectedMean=mean-cross*s.ldlt().solve(a*mean);
    const MatrixXd expectedP=prior-cross*s.ldlt().solve(cross.transpose());
    const double meanError=(runner.mean()-expectedMean).cwiseAbs().maxCoeff();
    const double covarianceError=(runner.covariance()-expectedP).cwiseAbs().maxCoeff();
    require(meanError<2e-10&&covarianceError<2e-10,"independent full-posterior oracle mismatch");
    std::cout<<std::setprecision(15)<<"CONDITION n="<<n<<" steps="<<steps<<" block="<<block<<" seconds="<<seconds
        <<" builds="<<runner.builds()<<" mean_error="<<meanError<<" covariance_error="<<covarianceError<<'\n';
    MatrixXd repeated=MatrixXd::Zero(2,n);repeated(0,0)=1;repeated(1,0)=2;
    auto same=runner.prepare(repeated,VectorXd::Zero(2));
    require(same.valid&&same.deterministicRows==2,"repeated equality rejected");
    VectorXd conflict=VectorXd::Zero(2);conflict(1)=1;
    require(!runner.prepare(repeated,conflict).valid,"deterministic affine conflict accepted");
    MatrixXd bad=MatrixXd::Identity(3,3);bad(0,1)=bad(1,0)=2;
    Runner nonPsd(VectorXd::Zero(3),bad);
    require(!nonPsd.prepare(MatrixXd::Identity(3,3),VectorXd::Zero(3)).valid,"indefinite covariance accepted");
#ifndef CONDITION_BASELINE
    ZhangIncrementalIntegerConditioner owned(VectorXd::Zero(3),MatrixXd::Identity(3,3));
    MatrixXd row=MatrixXd::Zero(1,3);row(0,0)=1;
    auto ticket=owned.prepare(row,VectorXd::Ones(1));auto stale=owned.prepare(row,VectorXd::Ones(1));
    require(owned.mean().isZero(),"prepare changed current mean");
    require(owned.commit(std::move(ticket))&&!owned.commit(std::move(stale)),"stale ticket accepted");
    require(!owned.commit(std::move(ticket)),"ticket reused");
    ZhangIncrementalIntegerConditioner other(VectorXd::Zero(3),MatrixXd::Identity(3,3));
    auto foreign=other.prepare(row,VectorXd::Ones(1));
    require(!owned.commit(std::move(foreign)),"foreign posterior ticket accepted");
    MatrixXd psd=MatrixXd::Identity(3,3);psd(2,2)=0;
    ZhangIncrementalIntegerConditioner singular(VectorXd::Zero(3),psd);
    row.setZero();row(0,2)=1;
    require(singular.prepare(row,VectorXd::Zero(1)).result.valid,"valid covariance-null row rejected");
    require(!singular.prepare(row,VectorXd::Ones(1)).result.valid,"null-space affine conflict accepted");
    Runner deterministic(VectorXd::Zero(3),MatrixXd::Zero(3,3));
    require(deterministic.prepare(MatrixXd::Identity(3,3),VectorXd::Zero(3)).valid,
        "zero-rank posterior rejected");
    require(!deterministic.prepare(MatrixXd::Identity(3,3),VectorXd::Ones(3)).valid,
        "zero-rank contradiction accepted");
    require(runner.builds()<=std::size_t(steps/16+2),"full decomposition repeated per block");
#endif
    std::cout<<"PASS repeated/contradictory constraints, PSD rejection, independent oracle, transactional root\n";
}

static void suite(const std::string& output)
{
    std::mt19937_64 random(20261020);
    std::normal_distribution<double> normal;
    std::ofstream out(output,std::ios::binary);
    int total=0;
    for(int sample=0;sample<160;++sample)
    {
        const int n=8+sample%33;
        MatrixXd root(n,n);VectorXd mean(n),truth(n);
        for(int i=0;i<n;++i){truth(i)=i%5-2;mean(i)=truth(i)+.03*normal(random);for(int j=0;j<n;++j)root(i,j)=.1*normal(random);}
        if(sample%3==0)root.rightCols(n/2).setZero();
        MatrixXd p=root*root.transpose();
        if(sample%3==1)p.diagonal().array()+=1;
        if(sample>=80)
        {
            VectorXd scale(n);
            for(int i=0;i<n;++i)scale(i)=std::pow(10.,-5.+10.*i/(n-1));
            p=scale.asDiagonal()*p*scale.asDiagonal();
            mean=truth;
        }
        if(sample>=120)
        {
            p=(root*root.transpose()+MatrixXd::Identity(n,n))*1e8;
            mean=truth;
        }
        if(sample%11==0){p.setIdentity();p(n-1,n-1)=1e-15;mean=truth;}
        Runner runner(mean,p);
        for(int step=0;step<14;++step)
        {
            const int count=1+step%4;
            MatrixXd a=MatrixXd::Zero(count,n);
            for(int row=0;row<count;++row){a(row,(step*3+row)%n)=1;a(row,(step*7+row+1)%n)+=row%2?-1:2;}
            if(step%4==2&&count>1)a.row(count-1)=2*a.row(0);
            VectorXd z=a*truth;
            if(step%7==6)z(count-1)+=1;
            const auto& result=runner.prepare(a,z);
            write(out,n);emit(out,result);++total;
            if(result.valid&&step%5!=3)runner.commit(); // abandoned proof too
        }
    }
    require(bool(out),"suite output failed");
    std::cout<<"SUITE cases="<<total<<" random systems=160 includes PSD, duplicate, conflict, abandoned proof, rank boundary, scaled full-rank\n";
}

int main(int argc,char** argv)
{
    ZhangArThreading::ScopedThreads threads(true);
    if(argc>1&&std::string(argv[1])=="suite"){suite(argv[2]);return 0;}
    if(argc>1&&std::string(argv[1])=="replay"){replay(argv[2],argv[3]);return 0;}
    contracts(argc>1?std::stoi(argv[1]):96,argc>2?std::stoi(argv[2]):20,argc>3?std::stoi(argv[3]):3);
}
