#include "common/zhangProductConditioningReceipt.hpp"
#include "common/zhangSmallProductResponse.hpp"
#include "common/zhangR48Marginal.hpp"
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

static void require(bool condition,const char* reason)
{ if(!condition)throw std::runtime_error(reason); }

static KFKey key(const char* receiver,int number)
{
    KFKey result;result.type=KF::AMBIGUITY;result.str=receiver;result.num=number;return result;
}

int main()
{
    try {
#ifdef _WIN32
        _putenv_s("ZHANG_R51_ENABLE","1");_putenv_s("ZHANG_R51_RATIO_ONLY","1");
#else
        setenv("ZHANG_R51_ENABLE","1",1);setenv("ZHANG_R51_RATIO_ONLY","1",1);
#endif
        const KFKey a=key("A",1),b=key("B",2),zero=key("unused",3);
        const std::vector<ZhangCurrentEpochNumericConstraint> recipe{
            {{{a,2},{b,1}},5},{{{b,1}},1}};
        ZhangProductConditioningIdentity identity;
        identity.runtimeId="run";identity.epochGpstSeconds=123456;
        identity.posteriorId="exact-ambiguity-moment";identity.chartId="signal-chart";
        identity.sourceRoot="whole-FLOAT-xP";identity.sourceStateOrder="full-order";
        identity.resultRoot="successful-whole-private-xP";identity.resultStateOrder="full-order";
        ZhangProductConditioningReceipt receipt;std::string reason;
        {
            bool violated=false,published=false,checkpointed=false;
            try {
                zhangRequireUnchangedFloatAuthority(false,violated);
                published=true;checkpointed=true;
            } catch(const std::runtime_error& error) {
                require(std::string(error.what())=="ZHANG_FLOAT_AUTHORITY_CHANGED_NO_PUBLICATION_OR_CHECKPOINT",
                    "FLOAT authority failure lost its terminal diagnostic");
            }
            require(violated && !published && !checkpointed,
                "damaged FLOAT proceeded to publication or checkpoint");
            // Simulate an outer catch trying a later FLOAT fallback. A PASS
            // argument must never clear the already latched violation.
            bool retried=false;
            try {zhangRequireUnchangedFloatAuthority(true,violated);retried=true;}
            catch(const std::runtime_error&) {}
            require(!retried,"caught FLOAT failure was later allowed to publish");
            bool clean=false;zhangRequireUnchangedFloatAuthority(true,clean);
            require(!clean,"unchanged FLOAT was falsely rejected");
        }
        require(receipt.bind(identity,recipe,reason),"valid successful-conditioning binding rejected");
        require(receipt.matches(identity,recipe,reason),"same complete numerical recipe rejected");
        {
            // Independent expected recipe: same two actual equations, merely
            // permuted/repeated and with an explicit zero-valued coefficient.
            auto reordered=recipe;std::swap(reordered[0],reordered[1]);
            reordered[0].rowCycles[zero]=-0.0;reordered.push_back(recipe[0]);
            reordered.push_back({{{zero,0}},0});
            require(receipt.matches(identity,reordered,reason),"canonical row equality failed");
            require(receipt.constraints.size()==2,"duplicates became extra fixed rows");
        }
        {
            auto changed=recipe;changed[0].value=6;
            require(!receipt.matches(identity,changed,reason),"same rank but changed RHS reused posterior");
            changed=recipe;changed.push_back({recipe[0].rowCycles,6});
            require(!receipt.matches(identity,changed,reason),"conflicting duplicate RHS reused posterior");
        }
        {
            auto changed=recipe;changed.push_back({{{a,1}},2});
            require(!receipt.matches(identity,changed,reason),"additional condition reused old posterior");
            changed=recipe;changed[0].rowCycles[a]=3;
            require(!receipt.matches(identity,changed,reason),"same rank but changed row reused posterior");
            changed=recipe;changed.erase(changed.begin());
            require(!receipt.matches(identity,changed,reason),"missing condition reused posterior");
            changed=recipe;for(auto& [coordinate,coefficient]:changed[0].rowCycles)coefficient*=2;
            changed[0].value*=2;
            require(!receipt.matches(identity,changed,reason),"rank/proportionality replaced complete row binding");
        }
        for(int field=0;field<8;++field) {
            auto changed=identity;
            switch(field) {
                case 0:changed.runtimeId="another-run";break;
                case 1:++changed.epochGpstSeconds;break;
                case 2:changed.posteriorId="changed-moment";break;
                case 3:changed.chartId="changed-chart";break;
                case 4:changed.sourceRoot="changed-whole-FLOAT-xP";break;
                case 5:changed.sourceStateOrder="permuted-source-order";break;
                case 6:changed.resultRoot="mutated-private-result-xP";break;
                case 7:changed.resultStateOrder="permuted-result-order";break;
            }
            require(!receipt.matches(changed,recipe,reason),"changed full root/order/scope reused posterior");
        }
        {
            auto changed=recipe;changed[0].rowCycles[a]=std::numeric_limits<double>::infinity();
            require(!receipt.matches(identity,changed,reason),"nonfinite row reused posterior");
            changed=recipe;changed.push_back({{{zero,0}},1});
            require(!receipt.matches(identity,changed,reason),"inconsistent structural-zero row reused posterior");
            auto badIdentity=identity;badIdentity.sourceStateOrder="different-order";
            ZhangProductConditioningReceipt invalid;
            require(!invalid.bind(badIdentity,recipe,reason),"conditioner changing order was certified");
            require(!invalid.matches(identity,recipe,reason),"invalid receipt remained usable");
        }
        {
            std::map<ZhangCoordinateRow,ZhangExactInteger> applied;
            for(const auto& row:receipt.constraints)applied.emplace(row.rowCycles,row.value);
            std::vector<ZhangCurrentEpochNumericConstraint> extras;
            auto current=recipe;current.push_back({{{a,1}},2});
            require(zhangProductResponseExtraConditions(applied,current,extras,reason) && extras.size()==1,
                "exact response recipe difference did not retain new condition");
            auto missing=recipe;missing.erase(missing.begin());
            require(!zhangProductResponseExtraConditions(applied,missing,extras,reason),
                "cached posterior reused after an applied condition disappeared");
            auto conflict=recipe;conflict[0].value=6;
            require(!zhangProductResponseExtraConditions(applied,conflict,extras,reason),
                "cached posterior reused with changed applied RHS");
            auto proportional=recipe;for(auto& [key,value]:proportional[0].rowCycles)value*=2;
            proportional[0].value*=2;
            require(!zhangProductResponseExtraConditions(applied,proportional,extras,reason),
                "rank/proportionality replaced exact applied-recipe presence");
            auto tooMany=recipe;
            for(int i=0;i<65;++i)tooMany.push_back({{{key("extra",100+i),1}},0});
            require(!zhangProductResponseExtraConditions(applied,tooMany,extras,reason),
                "small response exceeded its 64-extra-row bound");
        }
        {
            // Independent scalar Gaussian conditioning oracle. It does not
            // call a production conditioner, a matrix inverse, or a second
            // copy of the response implementation.
            const KFKey c=key("C",3);
            const long double m[3]={0.2L,-0.1L,0.3L};
            const long double q[3][3]={{4,1,0.5L},{1,3,0.75L},{0.5L,0.75L,2}};
            long double m0[3],q0[3][3],mFinal[3],qFinal[3][3];
            KFState base;base.kfIndexMap={{a,0},{b,1},{c,2}};
            base.x.resize(3);base.P.resize(3,3);
            for(int i=0;i<3;++i) {
                m0[i]=m[i]+q[i][1]/q[1][1]*(1-m[1]);
                base.x[i]=static_cast<double>(m0[i]);
                for(int j=0;j<3;++j) {
                    q0[i][j]=q[i][j]-q[i][1]*q[1][j]/q[1][1];
                    base.P(i,j)=static_cast<double>(q0[i][j]);
                }
            }
            for(int i=0;i<3;++i) {
                mFinal[i]=m0[i]+q0[i][2]/q0[2][2]*(2-m0[2]);
                for(int j=0;j<3;++j)qFinal[i][j]=q0[i][j]-q0[i][2]*q0[2][j]/q0[2][2];
            }
            const std::vector<ZhangCoordinateRow> target{{{a,1},{c,2}},{{b,-1},{c,1}}};
            const std::vector<ZhangCurrentEpochNumericConstraint> extra{{{{c,1}},2}};
            const auto response=zhangBuildSmallProductResponse(base,target,extra);
            require(response.valid && response.mean.size()==3 && response.extraSelector.rows()==1,
                "small response construction failed");
            const double targetExtraCross=static_cast<double>(q0[0][2]+2*q0[2][2]);
            require(std::abs(response.covariance(0,2)-targetExtraCross)<1e-12 &&
                std::abs(response.covariance(0,2))>1,
                "small response dropped target-extra cross covariance");
            ZhangR48MarginalWorkspace workspace(response.mean,response.covariance);
            const auto result=workspace.project(response.targetProjector,response.extraSelector,response.extraValues);
            require(result.valid,"small response marginal rejected valid scalar oracle fixture");
            const long double t[2][3]={{1,0,2},{0,-1,1}};
            for(int r=0;r<2;++r) {
                long double mean=0;for(int i=0;i<3;++i)mean+=t[r][i]*mFinal[i];
                require(std::abs(result.mean[r]-static_cast<double>(mean))<1e-10,"response mean differs from scalar oracle");
                for(int s=0;s<2;++s) {
                    long double covariance=0;
                    for(int i=0;i<3;++i)for(int j=0;j<3;++j)covariance+=t[r][i]*qFinal[i][j]*t[s][j];
                    require(std::abs(result.covariance(r,s)-static_cast<double>(covariance))<1e-10,
                        "response covariance differs from scalar oracle");
                }
            }
            auto invalidTarget=target;invalidTarget[0][zero]=1;
            require(!zhangBuildSmallProductResponse(base,invalidTarget,extra).valid,
                "small response assumed an unknown KF variable was zero");
            ZhangExactWorkBudget budget;budget.maximumCheckpoints=8;ZhangSmallProductResponse cancelled;
            require(!zhangRunWithinExactWorkBudget(budget,[&] {
                return zhangBuildSmallProductResponse(base,target,extra);
            },cancelled) && !cancelled.valid,"small response ignored cooperative deadline");
        }
        std::cout<<"Conditioning recipe identity, response set difference and independent scalar posterior oracle tests passed\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
