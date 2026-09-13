#include "common/zhangR47History.hpp"
#include <random>
#include <iostream>
#include <stdexcept>
void require(bool ok) {if(!ok)throw std::runtime_error("R49 selection equivalence failed");}
int main() {
    require(!zhangR47AffineIntegerFeasible({{2}},{1},1));
    require(zhangR47AffineIntegerFeasible({{2}},{6},1));
    require(!zhangR47AffineIntegerFeasible({{1,1},{1,-1}},{0,1},2));
    require(!zhangR47AffineIntegerFeasible({{2,1},{0,1}},{0,1},2));
    std::mt19937 random(4901);
    int cases=0;
    for(int dimension=1;dimension<=8;++dimension)for(int trial=0;trial<40;++trial) {
        const int count=2*dimension+2;
        ZhangExactVector truth(dimension);for(auto& v:truth)v=int(random()%7)-3;
        ZhangExactMatrix rows(count,ZhangExactVector(dimension));ZhangExactVector rhs(count);
        std::vector<ZhangDecisionProofs> parents;std::vector<double> gain(count);
        for(int i=0;i<count;++i) {
            for(int c=0;c<dimension;++c){rows[i][c]=int(random()%5)-2;rhs[i]+=rows[i][c]*truth[c];}
            if(trial%2 && random()%3==0)rhs[i]+=1;
            auto p=std::make_shared<ZhangIntegerDecisionProof>();
            p->id=std::to_string(i);p->originalStatement="fixture";p->observationProvenance="seed4901";
            p->conditionalFailureBound=(trial%3==0?2e-4:1e-6);
            parents.push_back({p});gain[i]=random()%4;
        }
        const auto old=zhangR49SelectHistorySubsetReference(rows,rhs,parents,{},gain,dimension,0.001,0.00025);
        const auto now=zhangR47SelectHistorySubset(rows,rhs,parents,{},gain,dimension,0.001,0.00025);
        require(old.valid==now.valid && old.rows==now.rows && old.values==now.values &&
            old.selected==now.selected && old.reasons==now.reasons && old.parents==now.parents && old.risk==now.risk);
        ++cases;
    }
    std::cout<<"R49 selection reference PASS cases="<<cases<<" plus 4 integer-feasibility regressions\n";
}
