#include "common/zhangUserGonsisComponents.hpp"
#include "common/zhangUserGonsis.hpp"
#include "common/zhangUserGonsisOutput.hpp"

#include <iostream>
#include <stdexcept>

int main()
{
    unsigned checks=0;
    auto check=[&](bool valid,const char* label)
    {
        ++checks;
        if (!valid) throw std::runtime_error(std::string("FAIL ")+label);
    };
    try
    {
        ZhangUserGonsisComponentCertificate c0{true,true,"C0",28,0};
        auto c1=c0; c1.component="C1";
        check(zhangUserGonsisComponentPairCertified(c0,c0),"same current component pair authorised");
        check(!zhangUserGonsisComponentPairCertified(c0,c1),"cross component pair not integer authorised");
        for(int fault=0;fault<7;++fault)
        {
            auto changed=c0;
            if(fault==0) changed.pppArUsable=false;
            if(fault==1) changed.dualFrequencyArValid=false;
            if(fault==2) changed.component="UNRESOLVED";
            if(fault==3) changed.component="NONE";
            if(fault==4) changed.version++;
            if(fault==5) changed.alignmentGeneration++;
            if(fault==6) changed.version=-1;
            check(!zhangUserGonsisComponentPairCertified(c0,changed),"invalid or mismatched authority rejects pair");
        }
        auto secondBand=c0; secondBand.component="OTHER_BAND";
        const auto d0=zhangUserGonsisDualComponentDomain(c0,secondBand);
        auto changedSecond=secondBand; changedSecond.alignmentGeneration++;
        const auto d1=zhangUserGonsisDualComponentDomain(c1,secondBand);
        const auto d2=zhangUserGonsisDualComponentDomain(c0,changedSecond);
        check(!d0.empty() && d0!=d1 && d0!=d2,"both band component tuples govern domain");
        auto oddName=c0; oddName.component="a b:c|d";
        check(zhangUserGonsisComponentToken(oddName).find_first_of(" \t\n:|")==std::string::npos,
            "trace domain token has no provider whitespace delimiters");
        const auto groups=zhangUserGonsisPartitionDomains({d0,d1,d0,"",d2});
        check(groups.size()==3 && groups.at(d0)==std::vector<int>({0,2}),"partition preserves indices and rejects unknown domain");

        // Independent synthetic FLOAT chart. Each component has an unobserved
        // common fractional gauge relative to the chart reference. The SD
        // oracle is written explicitly, NOT by calling the partition/selector.
        const int count=12,dimension=2*count;
        VectorXd mean(dimension);
        MatrixXd factor=MatrixXd::Identity(dimension,dimension)*.01;
        std::vector<ZhangUserGonsisSatellite> satellites(count);
        for(int i=0;i<count;++i)
        {
            const int local=i%6;
            mean(2*i)=100+2*local+(i<6?.37:.62)+.001*local;
            mean(2*i+1)=85+local+(i<6?.11:-.43)+.0005*local;
            for(int j=0;j<dimension;++j)
                factor(2*i,j)+=.00003*(1+(i+j)%7);
            satellites[i]={"G"+std::to_string(i+1),"physical:"+std::to_string(i),
                1.2-.03*local,0,.0001,true};
        }
        const MatrixXd covariance=factor*factor.transpose();
        const auto frozenMean=mean;
        const auto frozenCovariance=covariance;
        check(std::abs(mean(0)-std::round(mean(0)))>.25,
            "individual nonchart coordinate deliberately not integer");
        MatrixXd allRows(20,dimension);
        VectorXd allValues(20),allVariances(20);
        MatrixXd independentSD=MatrixXd::Zero(10,dimension);
        for(int component=0;component<2;++component)
        {
            const int first=6*component;
            MatrixXd select=MatrixXd::Zero(12,dimension);
            std::vector<ZhangUserGonsisSatellite> localSatellites;
            for(int i=0;i<6;++i)
            {
                select(2*i,2*(first+i))=1;
                select(2*i+1,2*(first+i)+1)=1;
                localSatellites.push_back(satellites[first+i]);
            }
            const MatrixXd localP=select*covariance*select.transpose();
            auto solver=[](const VectorXd& values,const MatrixXd&)
            {
                ZhangUserGonsisIlsSolution result;
                result.complete=true;
                result.integers=values.unaryExpr([](double value){return std::floor(value+.5);});
                result.bestSquaredCost=1; result.secondSquaredCost=3;
                result.reason="SYNTHETIC_ACCEPTED_CONTRACT";
                return result;
            };
            const auto result=zhangResolveUserGonsis(select*mean,localP,localSatellites,solver);
            check(result.accepted && result.referenceIndex==0 && result.wideLaneRows.rows()==5,
                "component local reference preserves five authorised SD pairs");
            MatrixXd l1Oracle=MatrixXd::Zero(5,dimension);
            MatrixXd wlOracle=MatrixXd::Zero(5,dimension);
            for(int i=1;i<6;++i)
            {
                l1Oracle(i-1,2*(first+i))=1;
                l1Oracle(i-1,2*first)=-1;
                wlOracle.row(i-1)=l1Oracle.row(i-1);
                wlOracle(i-1,2*(first+i)+1)=-1;
                wlOracle(i-1,2*first+1)=1;
            }
            check((result.l1Rows*select-l1Oracle).norm()==0 &&
                (result.wideLaneRows*select-wlOracle).norm()==0,
                "production SD rows match independent hand specified oracle");
            const MatrixXd projected=result.l1Rows*localP*result.l1Rows.transpose();
            check((projected-l1Oracle*covariance*l1Oracle.transpose()).norm()<1e-14,
                "complete projected covariance equals independent full P oracle");
            independentSD.middleRows(5*component,5)=l1Oracle;
            allRows.middleRows(10*component,10)=result.constraintRows*select;
            allValues.segment(10*component,10)=result.constraintValues;
            allVariances.segment(10*component,10)=result.constraintVariances;
            localSatellites.resize(4);
            check(!zhangResolveUserGonsis((select*mean).head(8),localP.topLeftCorner(8,8),
                localSatellites,solver).accepted,"four satellites do not weaken four WL pair threshold");
        }
        const MatrixXd sdCovariance=independentSD*covariance*independentSD.transpose();
        check(sdCovariance.topRightCorner(5,5).norm()>1e-12,
            "independent components retain real cross covariance");
        const auto joint=zhangConditionUserGonsisOutput(mean,covariance,allRows,allValues,allVariances);
        check(joint.valid,"one full covariance joint output accepts both local decisions");
        // Analytic Gaussian conditioning oracle using a distinct full-pivot LU
        // solve and the subtraction covariance formula, not Joseph/helper reuse.
        const MatrixXd cross=covariance*allRows.transpose();
        MatrixXd innovation=allRows*cross;
        innovation.diagonal()+=allVariances;
        Eigen::FullPivLU<MatrixXd> independentSolve(innovation);
        const VectorXd expectedMean=mean+cross*independentSolve.solve(allValues-allRows*mean);
        const MatrixXd expectedCovariance=covariance-cross*independentSolve.solve(cross.transpose());
        check((joint.mean-expectedMean).norm()<1e-9 &&
            (joint.covariance-expectedCovariance).norm()<1e-10,
            "joint output mean complete P matches independent analytic oracle");
        check(joint.covariance.topRightCorner(12,12).norm()>1e-12,
            "joint output never zeros cross component covariance");
        check((mean.array()==frozenMean.array()).all() &&
            (covariance.array()==frozenCovariance.array()).all(),"all component operations preserve frozen FLOAT");
        std::cout<<"PASS user_gonsis_components_contract checks="<<checks
                 <<" independent_sd_oracle=1 complete_covariance_oracle=1\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return EXIT_FAILURE;
    }
}
