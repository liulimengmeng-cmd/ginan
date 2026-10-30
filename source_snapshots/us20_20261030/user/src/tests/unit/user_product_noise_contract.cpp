#include "common/acsConfig.hpp"
#include "pea/zhangPppAr.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>

// Production reader, independently specified clock/phase covariance. The
// expected residual covariance does not reuse the production square root.
extern "C" int __wrap_main()
{
    unsigned checks=0;
    auto check=[&](bool ok,const char* label){++checks;if(!ok)throw std::runtime_error(label);};
    try
    {
        const auto path=std::filesystem::current_path()/"user_product_noise_fixture.csv";
        const char* parameters[]{"CLOCK","PHASE","PHASE"};
        const char* codes[]{"NONE","L1C","L2W"};
        const double Q[3][3]{{.04,.01,-.002},{.01,.03,.004},{-.002,.004,.02}};
        std::ofstream csv(path);csv<<"synthetic full covariance, metres squared\n";
        for(int epoch : {100000,100030})
            for(int i=0;i<3;++i)for(int j=i;j<3;++j)
                csv<<epoch<<",PRODUCT_FIXED,G01,"<<parameters[i]<<','<<codes[i]
                   <<",G01,"<<parameters[j]<<','<<codes[j]<<','<<Q[i][j]<<'\n';
        for(int i=0;i<3;++i)
            csv<<"100060,PRODUCT_FIXED,G01,"<<parameters[i]<<','<<codes[i]
               <<",G01,"<<parameters[i]<<','<<codes[i]<<','<<Q[i][i]<<'\n';
        csv.close();
        acsConfig.zhangPppAr.user_use_full_product_covariance=true;
        acsConfig.zhangPppAr.product_solution="PRODUCT_FIXED";
        acsConfig.zhangPppAr.product_covariance_filename=path.string();
        GTime time;time.bigTime=100000;
        ZhangUserProductNoiseContract code,first,second;
        std::string failure;
        check(queryZhangInternalProductNoiseContract(time,SatSys("G01"),E_ObsCode::L1C,false,code,&failure),"production code noise contract");
        check(queryZhangInternalProductNoiseContract(time,SatSys("G01"),E_ObsCode::L1C,true,first,&failure),"production L1 noise contract");
        check(queryZhangInternalProductNoiseContract(time,SatSys("G01"),E_ObsCode::L2W,true,second,&failure),"production L2 noise contract");
        check(code.noiseId==first.noiseId && first.noiseId==second.noiseId,"same full joint source identity across code and both phase rows");
        check(first.parameterOrder.size()==3 && first.parameterWeights.size()==3 && first.units=="metre",
            "explicit covariance parameter order, measurement operator and units");
        check(first.statisticalTarget=="UNSPECIFIED_BY_LEGACY_CSV" && first.temporalParent=="UNAVAILABLE" &&
              first.stateProductCrossCovariance=="UNAVAILABLE","missing statistical and temporal contract remains unavailable");
        auto dot=[](const auto& a,const auto& b){double v=0;for(unsigned i=0;i<a.size();++i)v+=a[i]*b[i];return v;};
        check(std::abs(dot(first.factors,first.factors)-.05)<1e-12,"L1 uses full clock/phase cross term");
        check(std::abs(dot(second.factors,second.factors)-.064)<1e-12,"L2 uses full clock/phase cross term");
        check(std::abs(dot(first.factors,second.factors)-.036)<1e-12,"cross-frequency noise retained");
        check(std::abs(dot(code.factors,first.factors)-.03)<1e-12,"code/phase noise retained");
        time.bigTime=100030;ZhangUserProductNoiseContract next;
        check(queryZhangInternalProductNoiseContract(time,SatSys("G01"),E_ObsCode::L1C,true,next,&failure) &&
            next.noiseId!=first.noiseId,"distinct samples never collide as one white error");
        const auto otherPath=std::filesystem::current_path()/"user_product_noise_other_fixture.csv";
        std::ofstream other(otherPath);
        other<<"synthetic alternative source\n100030,PRODUCT_FIXED,G01,CLOCK,NONE,G01,CLOCK,NONE,.01\n";other.close();
        acsConfig.zhangPppAr.product_covariance_filename=otherPath.string();
        ZhangUserProductNoiseContract alternative;
        check(queryZhangInternalProductNoiseContract(time,SatSys("G01"),E_ObsCode::L1C,false,alternative,&failure) &&
            std::abs(dot(alternative.factors,alternative.factors)-.01)<1e-12 && alternative.noiseId!=next.noiseId,
            "reader source changes at the same epoch cannot reuse the old cache");
        acsConfig.zhangPppAr.product_covariance_filename=path.string();time.bigTime=100060;
        check(!queryZhangInternalProductNoiseContract(time,SatSys("G01"),E_ObsCode::L1C,true,next,&failure) &&
            failure=="COVARIANCE_CROSS_BLOCKS_UNSPECIFIED","absent cross blocks are not replaced with zero");
        std::filesystem::remove(path);
        std::filesystem::remove(otherPath);
        std::cout<<"USER_PRODUCT_NOISE_CONTRACT PASS checks="<<checks<<" production_reader=1 independent_Q_oracle=1 temporal_NA=1\n";
        return 0;
    }
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
