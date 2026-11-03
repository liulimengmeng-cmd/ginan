// Production header kernels linked with actual PEA objects. The controlled
// graph covers structural product compilation, not the entire closure stage.
#include "common/acsConfig.hpp"
#include "ambres/GNSSambres.hpp"
#include "common/algebra.hpp"
#include "common/blasThreading.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangProductRelationBasis.hpp"
#include "common/zhangR49PosteriorWorkspace.hpp"
#include "common/zhangR51QuotientCompatibility.hpp"
#include "pea/zhangReference.hpp"
#ifdef AR_OPTIMIZATION_PROFILE
#include "common/zhangArThreading.hpp"
#endif
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/resource.h>
#include <omp.h>

static void check(bool pass, const char* why)
{ if (!pass) throw std::runtime_error(why); }
static double elapsed(std::chrono::steady_clock::time_point start)
{ return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); }

static void rootBenchmark(int n, const std::string& output)
{
    MatrixXd l(n, 23);
    for (int j=0;j<23;++j) for(int i=0;i<n;++i) l(i,j)=0.02*std::sin(i*0.09+j*0.37);
    MatrixXd p=l*l.transpose();
    for(int i=0;i<n;++i)p(i,i)+=1+double(i)/n;
    std::vector<std::string> order;
    for(int i=0;i<n;++i)order.push_back("STATE_"+std::to_string(i));
    ZhangR49PosteriorWorkspace workspace;
#ifdef AR_OPTIMIZATION_PROFILE
    ZhangR49FloatRootScope pinned(p,order,"BENCH");
#endif
    const auto start=std::chrono::steady_clock::now();
    check(workspace.factor(p,order,"BENCH"),"FLOAT factor failed");
    MatrixXd root=workspace.squareRoot;
    root.col(0)*=0.5;
    check(workspace.formGram(std::move(root),order,"BENCH"),"conditioned Gram failed");
    check(workspace.factor(p,order,"BENCH"),"protected FLOAT restore failed");
    const double seconds=elapsed(start);
    const double reconstruction=(workspace.squareRoot*workspace.squareRoot.transpose()-p).norm()/p.norm();
    check(reconstruction<2e-12,"root reconstruction mismatch");
    std::cout<<"ROOT seconds="<<seconds<<" decompositions="<<workspace.decompositions
             <<" hits="<<workspace.hits<<" relative_reconstruction="<<reconstruction<<'\n';
    // Compare covariance, not a nonunique eigenvector/root basis.
    MatrixXd reconstructed=workspace.squareRoot*workspace.squareRoot.transpose();
    std::ofstream file(output+".root",std::ios::binary);
    file.write(reinterpret_cast<const char*>(reconstructed.data()),reconstructed.size()*sizeof(double));
}

static void compatibilityBenchmark(int n)
{
    ZhangR51QuotientCompatibility compatibility(n,true);
    const auto start=std::chrono::steady_clock::now();
    for(int first=0;first<n;first+=24)
    {
        const int count=std::min(24,n-first);
        ZhangExactMatrix rows(count,ZhangExactVector(count));
        ZhangExactVector values(count);
        std::vector<int> subset;
        for(int i=0;i<count;++i){rows[i][i]=1;values[i]=first+i;subset.push_back(first+i);}
        auto trial=compatibility.assess(rows,values,subset);
        check(trial.compatible() && compatibility.commit(trial),"independent block rejected");
    }
    check(!compatibility.assess({{1}},{-100},{0}).compatible(),"conflict accepted");
    std::cout<<"COMPATIBILITY seconds="<<elapsed(start)<<" queries="<<compatibility.queries
             <<" committed_rows="<<compatibility.acceptedRows()<<'\n';
}

static void matrixText(std::ostream& out,const char* name,const ZhangExactMatrix& matrix)
{
    out<<name<<' '<<matrix.size()<<'\n';
    for(const auto& row:matrix){for(const auto& x:row)out<<x<<',';out<<'\n';}
}
static void vectorText(std::ostream& out,const char* name,const ZhangExactVector& row)
{ out<<name<<' ';for(const auto& x:row)out<<x<<',';out<<'\n'; }

static void closureBenchmark(int stations,int satellites,const std::string& output)
{
    // Explicit controlled estimator/product trees exercise nonzero targets.
    // The normal startup controller may choose equal trees, producing an
    // empty product lattice; that is not a valid performance fixture.
    std::set<ZhangGraphEdge> edges,preferred;
    const SatSys reference(E_Sys::GPS,1);
    for(int r=0;r<stations;++r)for(int s=1;s<=satellites;++s)
        edges.insert({"R"+std::to_string(r),SatSys(E_Sys::GPS,s)});
    for(int r=0;r<stations;++r)preferred.insert({"R"+std::to_string(r),reference});
    for(int s=2;s<=satellites;++s)preferred.insert({"R1",SatSys(E_Sys::GPS,s)});
    const auto current=zhangBuildSpanningTree(edges,"R0");
    const auto product=zhangBuildSpanningTree(edges,"R0",preferred);
    const auto start=std::chrono::steady_clock::now();
    auto basis=ProductRelationBasisBuilder::build(current,product,reference,E_Sys::GPS,E_ObsCode::L1C);
    check(basis.valid,"structural product compiler rejected controlled graph");
    ZhangExactMatrix named;
    for(const auto& row:basis.namedRelations)named.push_back(row.currentCycleCoefficients);
    basis.wholeLattice=zhangCompileWholeProductLattice(named,ZhangExactVector(named.size()),
        std::vector<bool>(basis.currentChords.size(),true));
    check(basis.wholeLattice.valid,"whole product compiler rejected controlled graph");
    const double seconds=elapsed(start);
    std::cout<<"CLOSURE_STRUCTURAL seconds="<<seconds<<" chords="<<basis.currentChords.size()
             <<" product_rank="<<basis.exactRank<<" mappable="<<basis.mappableTargetRank<<'\n';
#ifdef AR_OPTIMIZATION_PROFILE
    std::cout<<"EXACT_WITNESS fast_batches="<<zhangPrivateUnitMembershipMetrics.fastBatches
        <<" fallback_batches="<<zhangPrivateUnitMembershipMetrics.fallbackBatches
        <<" actual_workers="<<zhangPrivateUnitMembershipMetrics.maximumWorkers<<'\n';
#endif
    std::ofstream file(output+".basis");
    file<<basis.valid<<' '<<basis.exactRank<<' '<<basis.mappableTargetRank<<' '
        <<basis.failureReason<<' '<<basis.exactHnf<<' '<<basis.saturationIndex<<'\n';
    matrixText(file,"exactRowBasis",basis.exactRowBasis);
    matrixText(file,"networkIntegerBasis",basis.networkIntegerBasis);
    matrixText(file,"networkContainmentTransform",basis.networkContainmentTransform);
    matrixText(file,"whole.searchNamedRows",basis.wholeLattice.searchNamedRows);
    matrixText(file,"whole.searchPosteriorRows",basis.wholeLattice.searchPosteriorRows);
    matrixText(file,"whole.namedRecovery",basis.wholeLattice.namedRecovery);
    vectorText(file,"offsets",basis.affineOffsets);
    vectorText(file,"whole.recoveryOffsets",basis.wholeLattice.namedRecoveryOffsets);
    file<<std::setprecision(17)<<"transform "<<basis.transform<<'\n';
}

int main(int argc,char** argv)
{
    omp_set_dynamic(0);omp_set_num_threads(4);BlasThreading::openblas_set_num_threads(4);
    std::cout<<std::setprecision(15);
#ifdef AR_OPTIMIZATION_PROFILE
    ZhangArThreading::ScopedThreads profile(true);
#endif
    const std::string mode=argc>1?argv[1]:"all",output=argc>2?argv[2]:"ar_hotspot";
    if(mode=="root" || mode=="all")rootBenchmark(argc>3?std::stoi(argv[3]):1800,output);
    if(mode=="compat" || mode=="all")compatibilityBenchmark(argc>3?std::stoi(argv[3]):768);
    if(mode=="closure" || mode=="all")closureBenchmark(argc>3?std::stoi(argv[3]):32,argc>4?std::stoi(argv[4]):20,output);
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<"RESOURCES peak_rss_kib="<<usage.ru_maxrss<<" omp="<<omp_get_max_threads()
             <<" blas="<<BlasThreading::openblas_get_num_threads()<<'\n';
}
