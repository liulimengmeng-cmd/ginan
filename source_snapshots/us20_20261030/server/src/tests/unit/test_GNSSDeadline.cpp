#include "ambres/GNSSambres.hpp"
#include "common/zhangR51PhysicalImage.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// The production GNSS_AR implementation is linked below; no mock search,
// reference re-entry, or relaxed ratio threshold is used in this test.
int traceLevel=0;
static void require(bool condition,const char* message)
{
    if(!condition)throw std::runtime_error(message);
}
static void enable(const char* name)
{
#ifdef _WIN32
    _putenv_s(name,"1");
#else
    setenv(name,"1",1);
#endif
}
static GinAR_opt options(int dimension)
{
    GinAR_opt out;out.mode=E_ARmode::LAMBDA;out.nset=2;
    out.min_lambda_fix_count=dimension;out.max_lambda_fix_count=dimension;
    out.ratthr=3;out.sucthr=0.9999;return out;
}
static GinAR_mtx one(double value=0.03)
{
    GinAR_mtx out;out.aflt=VectorXd::Constant(1,value);
    out.Paflt=MatrixXd::Constant(1,1,1e-4);return out;
}
static void close(double actual,double expected,const char* message)
{
    require(std::abs(actual-expected)<=1e-8*std::max(1.0,std::abs(expected)),message);
}
int main()
{
    try {
        enable("ZHANG_R51_ENABLE");enable("ZHANG_R51_RATIO_ONLY");
        std::ostringstream trace;
        {
            auto matrix=one();auto opt=options(1);
            matrix.Ztrs=MatrixXd::Identity(1,1);matrix.zfix=VectorXd::Ones(1);
            opt.search_deadline=std::chrono::steady_clock::now()-std::chrono::seconds(1);
            require(GNSS_AR(trace,matrix,opt)==0,"expired deadline returned a fix");
            require(matrix.searchDiagnostic.wallTimedOut,"deadline not recorded");
            require(matrix.searchDiagnostic.reason=="WALL_TIMEOUT","wrong deadline reason");
            require(!matrix.searchDiagnostic.topTwoComplete && !matrix.searchDiagnostic.ratioExecuted,
                "deadline produced an authorized ratio");
            require(matrix.Ztrs.rows()==0 && matrix.zfix.size()==0,"timeout leaked candidate rows");
        }
        {
            auto matrix=one();auto opt=options(1);opt.lambda_node_limit=1;
            require(GNSS_AR(trace,matrix,opt)==0,"node-limited search returned a fix");
            require(matrix.searchDiagnostic.nodeLimited,"node limit not recorded");
            require(matrix.searchDiagnostic.reason=="NODE_LIMIT","wrong node-limit reason");
            require(!matrix.searchDiagnostic.topTwoComplete && !matrix.searchDiagnostic.ratioExecuted,
                "incomplete top two passed ratio");
            require(matrix.Ztrs.rows()==0 && matrix.zfix.size()==0,"node limit leaked candidate rows");
        }
        {
            // Exact preparation and GNSS_AR use the SAME absolute deadline;
            // entering enumeration never grants a fresh one-second budget.
            ZhangExactWorkBudget preparation;
            preparation.deadline=std::chrono::steady_clock::now()-std::chrono::milliseconds(1);
            ZhangR51PhysicalImage image;
            require(!zhangRunWithinExactWorkBudget(preparation,[&] {
                return zhangR51PhysicalImage({{1}}, {}, {},{{1}},1);
            },image),"expired preparation budget returned an exact image");
            auto matrix=one();auto opt=options(1);opt.search_deadline=preparation.deadline;
            require(GNSS_AR(trace,matrix,opt)==0 && matrix.searchDiagnostic.wallTimedOut,
                "enumeration reset the exhausted preparation deadline");
            require(!matrix.searchDiagnostic.topTwoComplete && matrix.zfix.size()==0,
                "shared exhausted budget leaked an integer fix");
        }
        {
            auto matrix=one();auto opt=options(1);
            opt.search_deadline=std::chrono::steady_clock::now()+std::chrono::seconds(5);
            require(GNSS_AR(trace,matrix,opt)==1,"complete 1D search rejected");
            const auto& d=matrix.searchDiagnostic;
            require(d.ilsComplete && d.topTwoComplete && d.ratioExecuted,"1D top two not complete");
            require(!d.wallTimedOut && !d.nodeLimited,"stale timeout/node-limit flag");
            close(d.bestSquaredDistance,9,"1D best distance differs from enumeration");
            close(d.secondSquaredDistance,9409,"1D second distance differs from enumeration");
            close(d.ratio,9409.0/9,"1D ratio changed");
            require(d.totalSeconds>=d.factorSeconds+d.reductionSeconds+d.enumerationSeconds,
                "phase timings exceed total search time");
        }
        {
            GinAR_mtx matrix;matrix.aflt.resize(2);matrix.aflt<<0.01,-0.02;
            matrix.Paflt.resize(2,2);matrix.Paflt<<1e-4,2e-5,2e-5,2e-4;
            // Independent exhaustive reference, NOT a second GNSS_AR call.
            const MatrixXd q=matrix.Paflt;const VectorXd mean=matrix.aflt;
            Eigen::LDLT<MatrixXd> decomposition(q);
            std::vector<std::pair<double,VectorXd>> brute;
            for(int a=-4;a<=4;++a)for(int b=-4;b<=4;++b) {
                VectorXd n(2);n<<a,b;const VectorXd delta=n-mean;
                brute.emplace_back(delta.dot(decomposition.solve(delta)),n);
            }
            std::sort(brute.begin(),brute.end(),[](const auto& a,const auto& b){return a.first<b.first;});
            auto opt=options(2);
            require(GNSS_AR(trace,matrix,opt)==2,"complete 2D search rejected");
            require(matrix.searchDiagnostic.topTwoComplete && matrix.searchDiagnostic.ratioExecuted,
                "2D incomplete top two");
            close(matrix.searchDiagnostic.bestSquaredDistance,brute[0].first,"2D best distance mismatch");
            close(matrix.searchDiagnostic.secondSquaredDistance,brute[1].first,"2D second distance mismatch");
            const VectorXd fixed=matrix.Ztrs.fullPivLu().solve(matrix.zfix);
            require((fixed-brute[0].second).cwiseAbs().maxCoeff()<1e-10,"2D fixed integer mismatch");
        }
        {
            auto matrix=one(0.5);auto opt=options(1);
            require(GNSS_AR(trace,matrix,opt)==0,"equal top-two candidate distances passed ratio=3");
            require(matrix.searchDiagnostic.topTwoComplete && matrix.searchDiagnostic.ratioExecuted,
                "tied top-two enumeration was incomplete");
            close(matrix.searchDiagnostic.ratio,1,"tied candidates ratio is not one");
        }
        {
            auto matrix=one();matrix.Paflt(0,0)=-1;auto opt=options(1);
            require(GNSS_AR(trace,matrix,opt)==0,"negative covariance returned a fix");
            require(matrix.searchDiagnostic.reason=="INVALID_COVARIANCE","wrong covariance reason");
            require(!matrix.searchDiagnostic.ratioExecuted,"invalid covariance executed ratio");
        }
        {
            auto matrix=one();matrix.Paflt.resize(1,2);matrix.Paflt.setOnes();auto opt=options(1);
            require(GNSS_AR(trace,matrix,opt)==0,"wrong covariance shape returned a fix");
            require(matrix.searchDiagnostic.reason=="INVALID_COVARIANCE","wrong shape reason");
        }
        {
            const ZhangExactMatrix targets{{1,0},{0,1}};
            const auto zero=zhangR51PhysicalImage(targets,targets,{3,5},targets,2);
            require(zero.valid && zero.projector.empty(),"determined image is not rank zero");
            require(zero.particularTarget==ZhangExactVector({3,5}),"rank-zero target values changed");
            const auto oneDirection=zhangR51PhysicalImage(targets,{{1,0}},{3},targets,2);
            require(oneDirection.valid && oneDirection.projector.size()==1,"one free signal searched as two");
            const auto twoDirections=zhangR51PhysicalImage(targets,{}, {},targets,2);
            require(twoDirections.valid && twoDirections.projector.size()==2,"two free signals lost rank");
            const auto structural=zhangR51PhysicalImage({{0,0},{0,1}}, {}, {},{{0,0},{0,1}},2);
            require(structural.valid && structural.projector.size()==1,
                "structural-zero frequency became a new AR direction");
        }
        {
            // Independent parity reasoning: 2*Nold+Nlive=1 permits exactly
            // odd Nlive. Direct projection must preserve 2Z, not saturate it.
            const ZhangExactMatrix target{{0,1}},history{{2,1}};
            const auto image=zhangR51PhysicalImage(target,history,{1},target,2);
            require(image.valid && image.projector.size()==1 && image.generators.size()==1,
                "odd-image projection failed");
            require(boost::multiprecision::abs(image.generators[0][0])==2,
                "2Z image was silently saturated to Z");
            ZhangExactMatrix lifted;ZhangExactVector values;
            require(zhangR51LiftRationalRows({{1}},{1},image.projector,image.offsets,lifted,values),
                "rational lift failed");
            require(lifted.size()==1 && lifted[0][0]==0 && lifted[0][1]!=0,
                "lift escaped visible target");
            require(values[0]%lifted[0][1]==0,"lift RHS not divisible");
            const ZhangExactInteger live=values[0]/lifted[0][1];
            require(live%2!=0,"rational image lift admitted even Nlive");
            auto rows=history;rows.insert(rows.end(),lifted.begin(),lifted.end());
            ZhangExactVector rhs{1};rhs.insert(rhs.end(),values.begin(),values.end());
            require(zhangR47AffineIntegerFeasible(rows,rhs,2),"lift conflicts with parity constraint");
            require(!zhangR47AffineIntegerFeasible({{2,1},{0,1}},{1,2},2),
                "even Nlive passed exact feasibility");
        }
        std::cout<<"GNSS deadline, node-limit, complete top-two and independent 1D/2D tests passed\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
