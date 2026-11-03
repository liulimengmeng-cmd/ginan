#include "common/zhangUserConditionalHold.hpp"
#include "common/zhangUserGonsisOutput.hpp"
#include <cassert>
#include <iostream>

int main()
{
    using Node=ZhangUserHoldNode;
    Node a{"A","arcA/product1","d"},b{"B","arcB/product1","d"},c{"C","arcC/product1","d"},d{"D","arcD/product1","d"};
    ZhangUserConditionalHold ledger; int added;
    assert(ledger.merge({{b,a,3,8},{c,a,5,12},{d,a,9,19}},1e-8,9e-6,added) && added==3);
    const auto pristine=ledger;
    assert(ledger.merge({{b,a,3,8},{c,b,2,4}},1e-8,9e-6,added) && added==0);
    assert(ledger.edges.size()==3 && ledger.noise==pristine.noise);
    assert(!ledger.merge({{c,b,2,5}},1e-8,9e-6,added));
    assert(ledger.edges.size()==3 && ledger.noise==pristine.noise);
    // Independent hand-written reference transform; deliberately no call to
    // ledger.path/project in the oracle. New rows C-B,D-B share the old B-A.
    const auto p=ledger.project({b,c,d});
    assert(p.edges.size()==2 && p.edges[0].wl==2 && p.edges[0].l1==4 && p.edges[1].wl==6 && p.edges[1].l1==11);
    Eigen::MatrixXd oracle=Eigen::MatrixXd::Zero(4,6);
    oracle(0,0)=-1;oracle(0,2)=1;oracle(1,1)=-1;oracle(1,3)=1;
    oracle(2,0)=-1;oracle(2,4)=1;oracle(3,1)=-1;oracle(3,5)=1;
    assert(p.transform==oracle && (p.noise-oracle*pristine.noise*oracle.transpose()).norm()<1e-20);
    assert(p.noise(0,2)==1e-8 && p.noise(1,3)==9e-6);
    auto retired=ledger;
    Node slippedA=a; slippedA.identity="arcA2/product1";
    assert(retired.retireReplaced({slippedA,b,c,d})==1);
    auto q=retired.project({b,c,d});assert(q.edges.size()==2 && (q.noise-p.noise).norm()<1e-20);
    assert(retired.project({slippedA,b}).edges.empty());
    Node slippedC=c;slippedC.identity="arcC2/product1";
    assert(ledger.retireReplaced({a,b,slippedC,d})==1);
    assert(ledger.project({a,b,d}).edges.size()==2 && ledger.project({slippedC,b}).edges.empty());
    ledger=pristine;
    Node newB=b;newB.identity="arcB/product2";newB.domain="d2";
    ledger.retireReplaced({a,newB,c,d});
    assert(ledger.project({a,c,d}).edges.size()==2 && ledger.project({a,newB}).edges.empty());
    // Missing endpoints suspend only; the exact identity can resume after a
    // product outage. No code treats a missing product as a zero measurement.
    ledger=pristine;auto missing=ledger.project({c,d});assert(missing.edges.size()==1);
    assert(ledger.project({a,b,c,d}).edges.size()==3);
    assert(ledger.project({}).edges.empty());
    // New domain cannot create cross-component continuity.
    assert(!ledger.merge({{newB,a,1,1}},1e-8,9e-6,added));
    // Independent direct Gaussian oracle and full covariance equivalence under
    // a reference change (not just matching positions or diagonal variances).
    Eigen::VectorXd x(6);x<<0.1,0.2,3.1,8.3,5.2,12.4;
    Eigen::MatrixXd t=Eigen::MatrixXd::Identity(6,6);t(0,4)=0.2;t(1,5)=0.3;
    Eigen::MatrixXd P=t*t.transpose();
    Eigen::MatrixXd H=Eigen::MatrixXd::Zero(4,6);
    H(0,0)=-1;H(0,1)=1;H(0,2)=1;H(0,3)=-1;
    H(1,0)=-.2;H(1,2)=.2;
    H(2,0)=-1;H(2,1)=1;H(2,4)=1;H(2,5)=-1;
    H(3,0)=-.2;H(3,4)=.2;
    Eigen::VectorXd z(4);z<<3,1.6,5,2.4;
    Eigen::MatrixXd R=Eigen::MatrixXd::Zero(4,4);R.diagonal()<<1e-8,9e-6,1e-8,9e-6;
    const auto xBefore=x;const auto pBefore=P;
    const auto output=zhangConditionUserGonsisOutput(x,P,H,z,R.diagonal(),&R);assert(output.valid);
    Eigen::MatrixXd S=H*P*H.transpose()+R;
    Eigen::MatrixXd K=P*H.transpose()*S.inverse();
    assert((output.mean-(x+K*(z-H*x))).norm()<1e-9);
    assert((output.covariance-(P-K*S*K.transpose())).norm()<1e-9);
    Eigen::MatrixXd B=Eigen::MatrixXd::Zero(4,4);
    B(0,0)=-1;B(1,1)=-1;B(2,0)=-1;B(2,2)=1;B(3,1)=-1;B(3,3)=1;
    Eigen::MatrixXd R2=B*R*B.transpose(),H2=B*H;Eigen::VectorXd z2=B*z;
    const auto changed=zhangConditionUserGonsisOutput(x,P,H2,z2,R2.diagonal(),&R2);assert(changed.valid);
    assert((output.mean-changed.mean).norm()<1e-9 && (output.covariance-changed.covariance).norm()<1e-9);
    // Simultaneous state-column reordering; stale numerical indices not stored.
    Eigen::MatrixXd J=Eigen::MatrixXd::Zero(6,6);for(int i=0;i<6;++i) J(i,5-i)=1;
    const auto reordered=zhangConditionUserGonsisOutput(J*x,J*P*J.transpose(),H*J.transpose(),z,R.diagonal(),&R);
    assert(reordered.valid && (J.transpose()*reordered.mean-output.mean).norm()<1e-9 &&
        (J.transpose()*reordered.covariance*J-output.covariance).norm()<1e-9);
    assert(x==xBefore && P==pBefore); // authoritative FLOAT is immutable
    assert(ZhangUserConditionalHold{}.edges.empty()); // cold restart has no integers
    bool overflow=false;try{zhangUserHoldChecked((__int128)std::numeric_limits<std::int64_t>::max()+1);}catch(...){overflow=true;}
    assert(overflow);
    bool invalid=false;try{zhangUserHoldInteger(2.01);}catch(...){invalid=true;}assert(invalid);
    // Repeated projection/merge does not accumulate evidence or shrink noise.
    ledger=pristine;
    for(int i=0;i<2000;++i) assert(ledger.merge({{c,b,2,4}},1e-8,9e-6,added)&&added==0);
    assert(ledger.edges.size()==3 && ledger.noise==pristine.noise);
    std::cout<<"PASS conditional hold: exact reference transport, full covariance oracle, selective slip/reset, gap suspension, dedup/conflict, immutable FLOAT, restart, overflow\n";
}
