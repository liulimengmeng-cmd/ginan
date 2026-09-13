#include "common/zhangR48Marginal.hpp"
#include "common/zhangR49PosteriorWorkspace.hpp"
#include <iostream>
#include <stdexcept>
void require(bool ok) {if(!ok)throw std::runtime_error("R49 posterior reference failed");}
int main() {
    for(int n=3;n<15;++n) {
        Eigen::MatrixXd l=Eigen::MatrixXd::Random(n,n);
        Eigen::MatrixXd p=l*l.transpose()+Eigen::MatrixXd::Identity(n,n);
        Eigen::VectorXd m=Eigen::VectorXd::Random(n);
        Eigen::MatrixXd h=Eigen::MatrixXd::Identity(2,n),j=Eigen::MatrixXd::Random(3,n);
        Eigen::VectorXd v=h*m+Eigen::VectorXd::Constant(2,0.02);
        ZhangR48MarginalWorkspace w(m,p);
        const Eigen::MatrixXd inv=(h*p*h.transpose()).inverse(),cross=j*p*h.transpose();
        auto a=w.project(j,h,v);
        require(a.valid && (a.mean-j*m-cross*inv*(v-h*m)).norm()<1e-10);
        require((a.covariance-(j*p*j.transpose()-cross*inv*cross.transpose())).norm()<1e-9);
        auto b=w.project(j,h,v);require(b.valid && w.targetHits==1);
        v(0)+=0.2;auto c=w.project(j,h,v);
        require(c.valid && w.decompositions==1 && w.hits==1);
        require((c.mean-j*m-cross*inv*(v-h*m)).norm()<1e-10);
    }
    Eigen::MatrixXd p=Eigen::MatrixXd::Identity(2,2);p(1,1)=0;
    ZhangR48MarginalWorkspace w(Eigen::VectorXd::Zero(2),p);
    Eigen::MatrixXd h(1,2);h<<0,1;
    require(!w.project(Eigen::MatrixXd::Identity(2,2),h,Eigen::VectorXd::Ones(1)).valid);
    require(w.project(Eigen::MatrixXd::Identity(2,2),h,Eigen::VectorXd::Zero(1)).valid);
    ZhangR49PosteriorWorkspace root;
    require(root.factor(p,{"a","b"},"t0"));require(root.factor(p,{"a","b"},"t0") && root.hits==1);
    require(root.factor(p,{"b","a"},"t0") && root.decompositions==2);
    require(root.factor(p,{"b","a"},"t1") && root.decompositions==3);
    p(0,0)=2;require(root.factor(p,{"b","a"},"t1") && root.decompositions==4);
    Eigen::VectorXd local(2);local<<std::sqrt(18.55717),0;
    auto gate=assessZhangIntegerCandidateNis(local,Eigen::MatrixXd::Identity(2,2),1e-6);
    require(gate.valid && gate.nis<gate.threshold);
    boost::math::chi_squared full(860);
    require(1065.08583+gate.nis>quantile(complement(full,1e-6)));
    std::cout<<"R49 posterior references PASS: means/covariance, RHS/order/epoch invalidation, null-space rejection, local-pass whole-reject\n";
}
