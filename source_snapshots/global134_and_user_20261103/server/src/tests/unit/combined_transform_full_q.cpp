// Production applyStateTransform versus a literal dense same-epoch oracle.
// This tests a combined transform; it does NOT claim that two graph receipts
// have been composed or allow MULTIPLE_GRAPH_EVENTS_UNCOMPOSED to pass.
#include "common/algebra.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool ok,const char* why) { if (!ok) throw std::runtime_error(why); }
int main()
{
    KFState state;
    KFKey phase,clock,ambiguity;
    phase.type=KF::PHASE_BIAS; phase.str="phase";
    clock.type=KF::REC_CLOCK; clock.str="clock";
    ambiguity.type=KF::AMBIGUITY; ambiguity.str="ambiguity";
    state.kfIndexMap={{KFState::oneKey,0},{phase,1},{clock,2},{ambiguity,3}};
    state.x.resize(4); state.x << 1,2,3,5;
    state.dx=VectorXd::Zero(4);
    state.P.resize(4,4);
    state.P << 0,0,0,0, 0,4,1,2, 0,1,9,3, 0,2,3,16;
    const VectorXd beforeX=state.x;
    const MatrixXd beforeP=state.P;
    KFState sequential=state;
    // Event A: phase'=phase+2*ambiguity; event B: clock'=clock-phase'.
    // The literal combined rows are phase+2a, clock-phase-2a, and a.
    std::map<KFKey,std::map<KFKey,double>> combined{
        {KFState::oneKey,{{KFState::oneKey,1}}},
        {phase,{{phase,1},{ambiguity,2}}},
        {clock,{{clock,1},{phase,-1},{ambiguity,-2}}},
        {ambiguity,{{ambiguity,1}}}};
    std::ostringstream trace;
    require(sequential.applyStateTransform(trace,{
        {KFState::oneKey,{{KFState::oneKey,1}}},{phase,{{phase,1},{ambiguity,2}}},
        {clock,{{clock,1}}},{ambiguity,{{ambiguity,1}}}},"SAME_EPOCH_EVENT_A"),"first same-epoch transform failed");
    require(sequential.applyStateTransform(trace,{
        {KFState::oneKey,{{KFState::oneKey,1}}},{phase,{{phase,1}}},
        {clock,{{clock,1},{phase,-1}}},{ambiguity,{{ambiguity,1}}}},"SAME_EPOCH_EVENT_B"),"second same-epoch transform failed");
    require(state.applyStateTransform(trace,combined,"INDEPENDENT_COMBINED_EVENT_FULL_Q"),"production combined transform failed");
    const double literal[4][4]={{1,0,0,0},{0,1,0,2},{0,-1,1,-2},{0,0,0,1}};
    // Scalar loops form the oracle; they call neither production transform
    // nor the same Eigen T*P*T^T expression as the implementation.
    const std::vector<KFKey> keys={KFState::oneKey,phase,clock,ambiguity};
    for (int i=0;i<4;++i)
    {
        double mean=0;
        for(int k=0;k<4;++k) mean+=literal[i][k]*beforeX(k);
        const int actualI=state.kfIndexMap.at(keys[i]);
        require(std::abs(state.x(actualI)-mean)<1e-12,"combined event mean differs from literal oracle");
        require(std::abs(sequential.x(sequential.kfIndexMap.at(keys[i]))-mean)<1e-12,
            "same-epoch sequential means differ from combined literal oracle");
        for(int j=0;j<4;++j)
        {
            double covariance=0;
            for(int k=0;k<4;++k) for(int l=0;l<4;++l)
                covariance+=literal[i][k]*beforeP(k,l)*literal[j][l];
            require(std::abs(state.P(actualI,state.kfIndexMap.at(keys[j]))-covariance)<1e-12,
                "combined event full Q/phase-clock cross covariance differs from literal oracle");
            require(std::abs(sequential.P(sequential.kfIndexMap.at(keys[i]),sequential.kfIndexMap.at(keys[j]))-covariance)<1e-12,
                "same-epoch sequential full Q differs from independent combined literal oracle");
        }
    }
    std::cout<<"combined_transform_full_q: PASS\n";
}
