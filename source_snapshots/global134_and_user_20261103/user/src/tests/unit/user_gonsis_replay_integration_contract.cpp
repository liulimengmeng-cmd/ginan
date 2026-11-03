#include "ambres/GNSSambres.hpp"
#include "common/zhangUserGonsisReplay.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{
int checks=0;
void require(bool condition,const char* reason)
{
    ++checks;
    if (!condition) throw std::runtime_error(reason);
}

struct OracleCandidate {double cost; VectorXd integer;};

// Independent exhaustive original-coordinate oracle. It does NOT call the
// production decorrelation, candidate recovery or LAMBDA enumeration helpers.
std::vector<OracleCandidate> oracle(const VectorXd& mean,const MatrixXd& Q)
{
    Eigen::LDLT<MatrixXd> solve(Q);
    require(solve.info()==Eigen::Success,"oracle covariance solve failed");
    std::vector<OracleCandidate> candidates;
    for (int a=-4;a<=4;++a) for (int b=-4;b<=4;++b)
        for (int c=-4;c<=4;++c) for (int d=-4;d<=4;++d)
        {
            VectorXd integer(4); integer<<a,b,c,d;
            const VectorXd delta=mean-integer;
            candidates.push_back({delta.dot(solve.solve(delta)),integer});
        }
    std::sort(candidates.begin(),candidates.end(),[](const auto& left,const auto& right)
        {return left.cost<right.cost;});
    Eigen::SelfAdjointEigenSolver<MatrixXd> spectrum(Q);
    const double outside=5-mean.cwiseAbs().maxCoeff();
    require(outside*outside/spectrum.eigenvalues().maxCoeff()>candidates[1].cost,
            "bounded oracle did not certify global top two");
    return candidates;
}

GinAR_mtx input()
{
    GinAR_mtx matrix;
    matrix.aflt.resize(4); matrix.aflt<<0.11,-0.27,0.33,0.17;
    matrix.Paflt.resize(4,4);
    matrix.Paflt<<.04,.015,.006,.002,
                  .015,.09,.01,.007,
                  .006,.01,.025,.008,
                  .002,.007,.008,.06;
    return matrix;
}
}

// Linked to ALL dependency-complete production objects via --wrap=main.
// Synthetic numeric contract only, not real-product correct-fix evidence.
extern "C" int __wrap_main()
{
    try
    {
        setenv("ZHANG_R51_RATIO_ONLY","1",1);
        setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS","0",1);
        GinAR_opt options;
        options.mode=E_ARmode::LAMBDA_ALT;
        options.nset=2; options.nitr=100;
        options.min_lambda_fix_count=4; options.max_lambda_fix_count=4;
        options.ratthr=2; options.force_bootstrap_success=false;
        options.maximum_lambda_conditional_sigma=0; options.lambda_candidate_nis_alpha=0;
        auto baseline=input();
        const auto source=baseline;
        const auto expected=oracle(source.aflt,source.Paflt);
        std::ostringstream baselineTrace;
        const int baselineRank=GNSS_AR(baselineTrace,baseline,options);
        require(baseline.searchDiagnostic.ilsComplete && baseline.searchDiagnostic.ratioExecuted,
                "production baseline top two was incomplete");
        ZhangUserGonsisLambdaTopTwo actual;
        int callbacks=0;
        options.lambdaReplayTopTwo=[&](const auto& topTwo){++callbacks; actual=topTwo;};
        auto diagnostic=input();
        std::ostringstream diagnosticTrace;
        const int diagnosticRank=GNSS_AR(diagnosticTrace,diagnostic,options);
        require(callbacks==1 && actual.inputCoordinatesValid,"complete full top two hook not emitted");
        require((actual.firstInput.array()==expected[0].integer.array()).all() &&
                (actual.secondInput.array()==expected[1].integer.array()).all(),
                "actual candidates differ from independent global original-coordinate oracle");
        const auto cost=[&](const VectorXd& integer)
        {const VectorXd delta=source.aflt-integer; return delta.dot(source.Paflt.ldlt().solve(delta));};
        require(std::abs(cost(actual.firstInput)-actual.firstSquaredCost)<1e-10 &&
                std::abs(cost(actual.secondInput)-actual.secondSquaredCost)<1e-10,
                "candidate costs do not match original mean/full covariance");
        require((actual.integerBasis*actual.firstInput-actual.firstReduced).norm()==0 &&
                (actual.integerBasis*actual.secondInput-actual.secondReduced).norm()==0,
                "basis direction or original coordinate labels incorrect");
        require(diagnosticRank==baselineRank && diagnostic.searchDiagnostic.reason==baseline.searchDiagnostic.reason &&
                diagnostic.searchDiagnostic.ratio==baseline.searchDiagnostic.ratio &&
                zhangUserGonsisMatrixBitsEqual(diagnostic.Ztrs,baseline.Ztrs) &&
                zhangUserGonsisMatrixBitsEqual(diagnostic.zfix,baseline.zfix),
                "diagnostic callback changed return/ratio/basis/candidate behavior");
        require(zhangUserGonsisMatrixBitsEqual(diagnostic.aflt,source.aflt) &&
                zhangUserGonsisMatrixBitsEqual(diagnostic.Paflt,source.Paflt),"diagnostic mutated FLOAT input");
        require(!diagnostic.searchDiagnostic.localNisExecuted && !diagnostic.lambda_candidate_nis_valid,
                "diagnostic unexpectedly calculated NIS");
        // The hook precedes the ratio veto: retain the actual rejected pair,
        // while preserving the solver's normal rejection and returned basis.
        callbacks=0;
        options.ratthr=actual.secondSquaredCost/actual.firstSquaredCost+1;
        auto rejected=input(); std::ostringstream rejectedTrace;
        require(GNSS_AR(rejectedTrace,rejected,options)==0 && callbacks==1 &&
                rejected.searchDiagnostic.reason=="RATIO_REJECTED" && actual.inputCoordinatesValid,
                "ratio-rejected complete top two not retained diagnostically");
        // A partial suffix cannot determine the full original integer vector.
        callbacks=0;
        options.min_lambda_fix_count=2; options.max_lambda_fix_count=2;
        auto partial=input(); std::ostringstream partialTrace;
        GNSS_AR(partialTrace,partial,options);
        require(callbacks==1 && actual.integerBasis.rows()==2 && actual.integerBasis.cols()==4 &&
                !actual.inputCoordinatesValid && actual.firstInput.size()==0,
                "partial reduced suffix was mislabelled as unique original coordinates");
        callbacks=0; options.mode=E_ARmode::OFF;
        auto off=input(); std::ostringstream offTrace;
        require(GNSS_AR(offTrace,off,options)==0 && callbacks==0,"no-enumeration path invoked top-two hook");
        std::cout<<"USER_GONSIS_REPLAY_INTEGRATION_CONTRACT PASS checks="<<checks
                 <<" production_GNSS_AR=1 independent_global_top2_oracle=1 original_coordinate_cost=1"
                 <<" ratio_veto_unchanged=1 partial_coordinate_claim_blocked=1 nis_evaluations=0"
                 <<" synthetic_contract_not_real_data_quality=1\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr<<"USER_GONSIS_REPLAY_INTEGRATION_CONTRACT FAIL "<<error.what()<<'\n';
        return EXIT_FAILURE;
    }
}
