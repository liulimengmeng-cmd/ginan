#include "common/zhangUserSharedPar.hpp"
#include "common/zhangUserResidualQuotient.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// Real GNSS_AR is linked from production objects. The independent candidate
// oracle below uses a literal 2x2 inverse and enumerates every integer in a
// bounded box; it does not call LAMBDA, LDLT, HNF or the production ratio gate.
extern "C" int __wrap_main()
{
    unsigned checks = 0;
    auto check = [&](bool ok, const char* label)
    {
        ++checks;
        if (!ok) throw std::runtime_error(std::string("FAIL ") + label);
    };
    try
    {
        setenv("ZHANG_R51_RATIO_ONLY", "1", 1);
        setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS", "0", 1);
        std::ostringstream trace;
        GinAR_opt options;
        options.ratthr = 2;
        options.nset = 2;
        options.sucthr = .9999;
        options.mode = E_ARmode::LAMBDA_ALT;
        options.min_lambda_fix_count = 1;
        options.lambda_candidate_nis_alpha = 0;
        for (int fixture = 0; fixture < 3; ++fixture)
        {
            GinAR_mtx input;
            input.aflt.resize(2); input.Paflt.resize(2,2);
            if (fixture == 0)
            {
                input.aflt << .01,.499;
                input.Paflt << .0001,0,0,1;
            }
            else if (fixture == 1)
            {
                input.aflt << .01,.499;
                input.Paflt << .0001,.003,.003,1;
            }
            else
            {
                // The precise coordinate is x+y, not x or y separately.
                input.aflt << -.489,.499;
                input.Paflt << 1.0001,-1,-1,1;
            }
            std::vector<double> distances;
            const double a=input.Paflt(0,0), b=input.Paflt(0,1), c=input.Paflt(1,1);
            const double determinant=a*c-b*b;
            for(int n0=-3;n0<=3;++n0) for(int n1=-3;n1<=3;++n1)
            {
                const double r0=n0-input.aflt(0),r1=n1-input.aflt(1);
                distances.push_back((c*r0*r0-2*b*r0*r1+a*r1*r1)/determinant);
            }
            std::sort(distances.begin(),distances.end());
            const double oracleJointRatio=distances[1]/distances[0];
            check(oracleJointRatio<2,"independent full integer ratio rejected");
            auto full=input;
            options.max_lambda_fix_count=2;
            check(GNSS_AR(trace,full,options)==0,"real full ILS rejected");
            check(full.searchDiagnostic.ratioExecuted && full.searchDiagnostic.ilsComplete,
                  "full ILS actually completed its ratio test");
            check(std::abs(full.searchDiagnostic.ratio-oracleJointRatio)<1e-7,
                  "joint ratio agrees with independent bounded exhaustive candidates");
            auto partial=input;
            ZhangUserSharedParAudit audit;
            check(zhangResolveUserSharedPar(trace,partial,options,&audit)==1,
                  "real PAR accepts reliable one-dimensional suffix");
            check(audit.testedDimensions==std::vector<int>({2,1}),"PAR actually tries full then partial");
            check(audit.acceptedRank==1 && audit.reason=="PARTIAL_ILS_ACCEPTED","PAR receipt truthful");
            check((partial.aflt-input.aflt).norm()==0 && (partial.Paflt-input.Paflt).norm()==0,
                  "every PAR search uses unchanged full float moments");
            check(partial.Ztrs.rows()==1 && partial.Ztrs.cols()==2 && partial.zfix.size()==1,
                  "accepted suffix maps to complete ambient coordinates");
            check(!partial.searchDiagnostic.localNisExecuted && !partial.lambda_candidate_nis_valid,
                  "no NIS evaluation performed");
            const VectorXd integerRow=partial.Ztrs.row(0).transpose();
            check((integerRow.array()-integerRow.array().round()).abs().maxCoeff()<1e-12,
                  "PAR row coefficients remain exact integers");
            if(fixture==2)
                check(std::abs(integerRow(0))==1 && std::abs(integerRow(1))==1,
                      "mixed reliable suffix is not replaced with principal covariance rows");
            const double scalarMean=integerRow.dot(input.aflt);
            const double scalarVariance=(integerRow.transpose()*input.Paflt*integerRow)(0,0);
            std::vector<std::pair<double,int>> scalarDistances;
            for(int integer=-5;integer<=5;++integer)
                scalarDistances.emplace_back((integer-scalarMean)*(integer-scalarMean)/scalarVariance,integer);
            std::sort(scalarDistances.begin(),scalarDistances.end());
            const double oraclePartialRatio=scalarDistances[1].first/scalarDistances[0].first;
            check(oraclePartialRatio>=2,"independent scalar ratio accepts");
            check(std::abs(partial.searchDiagnostic.ratio-oraclePartialRatio)<1e-5,
                  "partial ratio agrees with independent scalar integer enumeration");
            check(partial.zfix(0)==scalarDistances[0].second,"PAR best integer matches independent oracle");
            ZhangExactVector row{static_cast<long long>(std::llround(integerRow(0))),
                                 static_cast<long long>(std::llround(integerRow(1)))};
            const auto exact=zhangAuditUserIntegerUnion({row},
                {static_cast<long long>(std::llround(partial.zfix(0)))},2);
            check(exact.valid && exact.rank==1,"PAR mapped ambient integer union feasible");
        }
        GinAR_mtx singular;
        singular.aflt=VectorXd::Zero(2); singular.Paflt=MatrixXd::Zero(2,2);
        ZhangUserSharedParAudit bad;
        check(zhangResolveUserSharedPar(trace,singular,options,&bad)==0,"singular search rejects without jitter");
        check(bad.reason=="PAR_INPUT_NOT_POSITIVE_DEFINITE" && bad.testedDimensions.empty(),
              "singular input never reaches an integer search");
        GinAR_mtx half;
        half.aflt=VectorXd::Constant(2,.5); half.Paflt=MatrixXd::Identity(2,2);
        check(zhangResolveUserSharedPar(trace,half,options,&bad)==0,"all half-integer suffixes ratio-rejected");
        check(bad.testedDimensions==std::vector<int>({2,1}),"all rejected path genuinely attempts PAR");
        check(half.Ztrs.rows()==0 && half.zfix.size()==0,"failed candidates cannot escape into held ledger");
        GinAR_mtx precise;
        precise.aflt=VectorXd::Constant(2,.01); precise.Paflt=MatrixXd::Identity(2,2)*.0001;
        check(zhangResolveUserSharedPar(trace,precise,options,&bad)==2,"full ILS accepted when independently precise");
        check(bad.testedDimensions==std::vector<int>({2}) && bad.reason=="FULL_ILS_ACCEPTED",
              "full accepted path does not spend resources on unnecessary PAR");
        check(trace.str().find("nis_evaluations=0")!=std::string::npos,"auditable no-NIS receipts");
        // SAFE_RESEARCH must enforce independent uncertainty screening even
        // while global ratio-only remains on. The large scalar ratio is not a
        // substitute for absolute precision or a calibrated reliability claim.
        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE", "SAFE_RESEARCH", 1);
        MatrixXd support = MatrixXd::Zero(2,6);
        support(0,0)=-1; support(0,2)=1;
        support(1,0)=-1; support(1,4)=1;
        GinAR_mtx weak;
        weak.aflt=VectorXd::Constant(2,.01); weak.Paflt=MatrixXd::Identity(2,2);
        check(zhangResolveUserSharedPar(trace,weak,options,&bad,&support)==0,
            "ratio-only cannot bypass independent success gate");
        check(bad.independentSuccessGate && !bad.trialReasons.empty() &&
            bad.trialReasons[0]=="BEST_1D_BOOTSTRAP_TOO_LOW", "actual conditional-D success veto recorded");
        check(!weak.searchDiagnostic.localNisExecuted && !weak.lambda_candidate_nis_valid,
            "new success veto never restores NIS");
        GinAR_mtx precisionOnly;
        precisionOnly.aflt=VectorXd::Constant(2,.01);
        precisionOnly.Paflt=MatrixXd::Identity(2,2)*.04;
        options.sucthr=.9;
        check(zhangResolveUserSharedPar(trace,precisionOnly,options,&bad,&support)==0 &&
            bad.trialReasons[0]=="BEST_1D_CONDITIONAL_SIGMA_TOO_HIGH",
            "absolute conditional precision independent of ratio and permissive success threshold");
        options.sucthr=.9999;
        GinAR_mtx safePrecise;
        safePrecise.aflt=VectorXd::Constant(2,.01); safePrecise.Paflt=MatrixXd::Identity(2,2)*.0001;
        check(zhangResolveUserSharedPar(trace,safePrecise,options,&bad,&support)==2,
            "safe precise independent directions accepted");
        check(bad.trialBootstrapSuccess[0]>.9999 && bad.trialSupportSatellites[0]==3,
            "selected suffix success and physical support truthful");
        GinAR_mtx single;
        single.aflt=VectorXd::Constant(1,.01); single.Paflt=MatrixXd::Identity(1,1)*.0001;
        MatrixXd oneSupport=MatrixXd::Zero(1,4); oneSupport(0,0)=-1; oneSupport(0,2)=1;
        check(zhangResolveUserSharedPar(trace,single,options,&bad,&oneSupport)==0 &&
            bad.testedDimensions.empty(), "safe PAR cannot shrink to one dimension");
        GinAR_mtx budget;
        budget.aflt=VectorXd::Constant(6,.01); budget.Paflt=MatrixXd::Identity(6,6);
        MatrixXd budgetSupport=MatrixXd::Zero(6,14);
        for(int i=0;i<6;++i) {budgetSupport(i,0)=-1; budgetSupport(i,2*(i+1))=1;}
        check(zhangResolveUserSharedPar(trace,budget,options,&bad,&budgetSupport)==0 &&
            bad.reason=="PAR_TRIAL_BUDGET_EXHAUSTED" && bad.testedDimensions.size()==4,
            "trial budget exact no incomplete plan reported as accepted");
        check(trace.str().find("risk_claim=UNCALIBRATED_MODEL_CONDITIONAL_NOT_FFRT")!=std::string::npos,
            "screening does not claim fixed-failure-rate calibration");
        unsetenv("ZHANG_USER_AR_RELIABILITY_PROFILE");
        std::cout << "USER_SHARED_PAR_CONTRACT checks=" << checks
                  << " status=PASS real_GNSS_AR=1 independent_oracle=BOUNDED_INTEGER_ENUMERATION\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
