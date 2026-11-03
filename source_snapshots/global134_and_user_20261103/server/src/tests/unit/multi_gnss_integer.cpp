#include "common/zhangMultiGnssInteger.hpp"

#include <Eigen/Cholesky>
#include <iostream>
#include <random>
#include <stdexcept>

static void require(bool pass, const char* reason)
{
    if (!pass) throw std::runtime_error(reason);
}

int main()
{
    std::map<int, KFKey> keys;
    const std::array<std::pair<E_Sys, E_ObsCode>, 6> signals{{
        {E_Sys::GPS,E_ObsCode::L1C},{E_Sys::GPS,E_ObsCode::L2W},
        {E_Sys::GAL,E_ObsCode::L1C},{E_Sys::GAL,E_ObsCode::L5Q},
        {E_Sys::BDS,E_ObsCode::L1P},{E_Sys::BDS,E_ObsCode::L5P}}};
    for (int column = 0; column < signals.size(); ++column)
    {
        KFKey key;
        key.type=KF::AMBIGUITY;
        key.str="R";
        key.Sat=SatSys(signals[column].first,column<4?1:19);
        key.num=static_cast<int>(signals[column].second);
        keys[column]=key;
    }
    ZhangExactVector gps(6),gal(6),bds(6);
    gps[0]=1;gps[1]=-1;gal[2]=1;gal[3]=-1;bds[4]=1;bds[5]=-1;
    require(zhangMultiGnssIntegerRowsUseSystem({gps},keys,E_Sys::GPS),
        "foreign zero columns incorrectly removed or rejected");
    require(zhangMultiGnssIntegerRowsUseSystem({gal},keys,E_Sys::GAL),
        "Galileo L1C/L5Q integer row rejected");
    require(zhangMultiGnssIntegerRowsUseSystem({bds},keys,E_Sys::BDS),
        "BDS L1P/L5P integer row rejected");
    require(!zhangMultiGnssIntegerRowsUseSystem({gps,gal},keys,E_Sys::GPS),
        "foreign nonzero integer row accepted");
    auto mixed=gps;mixed[2]=1;
    require(!zhangMultiGnssIntegerRowsUseSystem({mixed},keys,E_Sys::GPS),
        "mixed GPS/Galileo integer relation accepted");
    auto missing=keys;missing.erase(0);
    require(!zhangMultiGnssIntegerRowsUseSystem({gps},missing,E_Sys::GPS),
        "missing coordinate was silently zero-filled");
    const std::array<E_ObsCode,2> gCodes{E_ObsCode::L1C,E_ObsCode::L2W};
    const std::array<E_ObsCode,2> eCodes{E_ObsCode::L1C,E_ObsCode::L5Q};
    const std::array<E_ObsCode,2> cCodes{E_ObsCode::L1P,E_ObsCode::L5P};
    using PhysicalRow=std::map<std::string,ZhangExactInteger>;
    require(zhangMultiGnssPhysicalRowsUseSystem({PhysicalRow{{"L1C|R|G01|V2",1}}},
        E_Sys::GPS,gCodes),"GPS physical identity rejected");
    require(zhangMultiGnssPhysicalRowsUseSystem({PhysicalRow{{"L5Q|R|E01|V3",1}}},
        E_Sys::GAL,eCodes),"Galileo physical identity rejected");
    require(zhangMultiGnssPhysicalRowsUseSystem({PhysicalRow{{"L1P|R|C19|V4",1}}},
        E_Sys::BDS,cCodes),"BDS physical identity rejected");
    require(!zhangMultiGnssPhysicalRowsUseSystem({PhysicalRow{{"L1C|R|E01|V2",1}}},
        E_Sys::GPS,gCodes),"same-PRN foreign physical identity accepted");
    require(!zhangMultiGnssPhysicalRowsUseSystem({PhysicalRow{{"L5X|R|E01|V3",1}}},
        E_Sys::GAL,eCodes),"tracking code silently substituted by X");
    require(!zhangMultiGnssPhysicalRowsUseSystem({PhysicalRow{{"L1P|R|C19|Vx",1}}},
        E_Sys::BDS,cCodes),"malformed physical arc accepted");

    // Independent full-rank Schur-complement reference: do not call the
    // production conditioner to construct the reference posterior.
    std::mt19937 random(20261021);
    std::normal_distribution<double> normal;
    for (int trial=0;trial<100;++trial)
    {
        constexpr int n=9;
        MatrixXd seed(n,n);
        VectorXd mean(n);
        for(int r=0;r<n;++r)
        {
            mean(r)=normal(random);
            for(int c=0;c<n;++c)seed(r,c)=normal(random);
        }
        const MatrixXd p=seed*seed.transpose()+MatrixXd::Identity(n,n);
        const VectorXd beforeMean=mean;
        const MatrixXd beforeP=p;
        MatrixXd a=MatrixXd::Zero(3,n);
        a(0,0)=1;a(0,1)=-1;
        a(1,2)=1;a(1,3)=-1;
        a(2,4)=1;a(2,5)=-1;
        VectorXd values(3);values<<2,-1,3;
        const MatrixXd ap=a*p;
        const MatrixXd gram=ap*a.transpose();
        const MatrixXd gain=p*a.transpose()*gram.ldlt().solve(MatrixXd::Identity(3,3));
        const VectorXd referenceMean=mean+gain*(values-a*mean);
        const MatrixXd referenceP=p-gain*ap;
        VectorXd currentMean=mean;
        MatrixXd currentP=p;
        for(int system=0;system<3;++system)
        {
            const MatrixXd row=a.row(system);
            const VectorXd value=values.segment(system,1);
            const auto conditioned=zhangConditionPosteriorEffectiveIntegers(
                currentMean,currentP,row,value);
            require(conditioned.valid,"sequential joint conditioner failed");
            currentMean=conditioned.mean;
            currentP=conditioned.covariance;
        }
        require((currentMean-referenceMean).cwiseAbs().maxCoeff()<1e-10,
            "sequential full-state mean differs from independent joint Schur reference");
        require((currentP-referenceP).cwiseAbs().maxCoeff()<1e-10,
            "sequential full-state covariance differs from independent joint Schur reference");
        require((mean-beforeMean).norm()==0 && (p-beforeP).norm()==0,
            "authoritative FLOAT input changed");
        require(std::abs(currentP(6,7))>1e-8 || std::abs(currentP(7,8))>1e-8,
            "joint nuisance covariance was block-diagonalized");
    }
    std::cout<<"PASS multi-GNSS integer scope, exact tracking identity, and 100 independent joint covariance references\n";
}
