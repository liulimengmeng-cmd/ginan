#include "common/zhangUserProductJointMath.hpp"
#include "common/zhangUserProductNoiseContract.hpp"
#include <random>
#include <iostream>

int main()
{
    unsigned checks=0;
    auto check=[&](bool ok,const char* label){++checks;if(!ok)throw std::runtime_error(label);};
    auto rejects=[&](auto call){try {call();return false;} catch(const std::runtime_error&){return true;}};
    try
    {
        using Eigen::MatrixXd; using Eigen::VectorXd;
        MatrixXd T=MatrixXd::Identity(1,1),E=MatrixXd::Ones(1,1),K=MatrixXd::Zero(1,2);
        K(0,0)=1; VectorXd shift=VectorXd::Constant(1,.37),k=VectorXd::Zero(1),mean(2);
        MatrixXd joint(2,2); joint<<.040009,-.04,-.04,.05; mean<<4,.2;
        auto moved=zhangUserTransportJointChart(mean,joint,T,E,shift,K,k,"SAME_INSTANT_CHART_DISPLACEMENT",true);
        check(std::abs(moved.covariance(0,0)-.010009)<1e-14,"negative state/displacement cross block is included");
        check(std::abs(moved.covariance(0,1)-.01)<1e-14 && moved.covariance(1,1)==.05,
            "displacement and transformed cross block are retained");
        check(std::abs((moved.integerMap*moved.mean+moved.integerOffset)(0)-4)<1e-14,
            "integer functional retains delta and fractional affine offset");
        check(rejects([&]{zhangUserTransportJointChart(mean,joint,T,E,shift,K,k,"CURRENT_PRODUCT_Q",true);}),
            "Qt cannot masquerade as same-instant Qdelta");
        check(rejects([&]{zhangUserTransportJointChart(mean,joint,T,E,shift,K,k,"SAME_INSTANT_CHART_DISPLACEMENT",false);}),
            "missing cross covariance is not silently assumed zero");
        MatrixXd bad=joint;bad(0,1)=bad(1,0)=.1;
        check(rejects([&]{zhangUserTransportJointChart(mean,bad,T,E,shift,K,k,"SAME_INSTANT_CHART_DISPLACEMENT",true);}),"indefinite joint posterior rejected");
        std::mt19937 generator(20261009); std::normal_distribution<double> normal;
        for (double cross : {-.04,0.0,.04})
        {
            joint<<.040009,cross,cross,.05;
            const auto result=zhangUserTransportJointChart(mean,joint,T,E,shift,K,k,"SAME_INSTANT_CHART_DISPLACEMENT",true);
            Eigen::SelfAdjointEigenSolver<MatrixXd> solver(joint);
            const MatrixXd root=solver.eigenvectors()*solver.eigenvalues().cwiseSqrt().asDiagonal();
            double sum=0,sum2=0; const int count=40000;
            for(int i=0;i<count;++i){VectorXd z(2);z<<normal(generator),normal(generator);
                const VectorXd sample=root*z; const double value=sample(0)+sample(1);sum+=value;sum2+=value*value;}
            const double observed=(sum2-sum*sum/count)/(count-1);
            const double expected=.040009+.05+2*cross;
            check(std::abs(result.covariance(0,0)-expected)<1e-14,"independent joint covariance formula");
            check(std::abs(observed-expected)<6*expected*std::sqrt(2.0/(count-1)),"Monte Carlo covariance within six sampling standard errors");
        }
        joint.setZero(); joint(0,0)=.03;
        for(int n=-3;n<=3;++n)
        {
            mean<<n,0;T(0,0)=2;
            const auto result=zhangUserTransportJointChart(mean,joint,T,E,shift,K,k,"SAME_INSTANT_CHART_DISPLACEMENT",true);
            check(std::abs((result.integerMap*result.mean+result.integerOffset)(0)-n)<1e-14 && result.integerMap(0,0)==.5,
                "non-unimodular real scaling preserves the restricted integer lattice, not arbitrary Z");
        }
        MatrixXd R(3,3);R<<2,.2,.1,.2,3,.4,.1,.4,1;
        MatrixXd B=MatrixXd::Ones(3,1);
        const MatrixXd projection=zhangUserUnknownChartProjection(R,B);
        check(projection.rows()==2 && (projection*B).norm()<1e-12,"unknown common chart direction exactly removed");
        check((projection*R*projection.transpose()-MatrixXd::Identity(2,2)).norm()<1e-12,"projected observation noise whitened");
        const MatrixXd weight=R.inverse()-R.inverse()*B*(B.transpose()*R.inverse()*B).inverse()*B.transpose()*R.inverse();
        check((projection.transpose()*projection-weight).norm()<1e-12,"independent nuisance elimination information oracle");
        check(zhangUserUnknownChartProjection(R,MatrixXd::Identity(3,3)).rows()==0,
            "one unconstrained chart offset per row contributes zero phase information");
        // WHITE augmented posterior and R+GQG^T must have the same x marginal.
        MatrixXd P(2,2);P<<2,.4,.4,1;MatrixXd H(3,2);H<<1,.2,-.3,1,.5,-.7;
        MatrixXd G(3,2);G<<1,0,.2,1,-.1,.3;MatrixXd Q(2,2);Q<<.4,.1,.1,.3;
        MatrixXd prior=MatrixXd::Zero(4,4);prior.topLeftCorner(2,2)=P;prior.bottomRightCorner(2,2)=Q;
        MatrixXd A(3,4);A<<H,G;VectorXd innovation(3);innovation<<.3,-.2,.1;
        const MatrixXd S=H*P*H.transpose()+R+G*Q*G.transpose();
        const MatrixXd gain=prior*A.transpose()*S.inverse();
        const MatrixXd posterior=prior-gain*A*prior;
        check((posterior.topLeftCorner(2,2)-(P-P*H.transpose()*S.inverse()*H*P)).norm()<1e-12,
            "WHITE augmentation equals correlated R marginal covariance");
        check(((gain*innovation).head(2)-P*H.transpose()*S.inverse()*innovation).norm()<1e-12,
            "WHITE augmentation equals correlated R marginal mean");
        ZhangUserProductNoiseContract packet;packet.noiseId="same-error";packet.sample="1";
        std::map<std::string,std::string> metadata;
        zhangUserClaimProductNoise(metadata,packet);zhangUserClaimProductNoise(metadata,packet);
        check(metadata.size()==2,"shared measurement rows do not duplicate the noise owner");
        packet.representation="AUGMENTED_STATE";
        check(rejects([&]{zhangUserClaimProductNoise(metadata,packet);}),"double R/state ownership rejected");
        packet.sample="2";zhangUserClaimProductNoise(metadata,packet);
        check(metadata.size()==2,"noise owner registry is bounded by one sample");
        std::cout<<"USER_PRODUCT_JOINT_MATH_CONTRACT PASS checks="<<checks
                 <<" independent_covariance_oracle=1 monte_carlo=1 unknown_nuisance_projection=1"
                 <<" white_augmentation_equivalence=1 integer_map=1 duplicate_noise_guard=1 NIS_computed=0\n";
        return 0;
    }
    catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
