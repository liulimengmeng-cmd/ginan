#include "common/zhangR47History.hpp"
#include <fstream>
#include <chrono>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
    if(argc!=2)return 2;
    std::ifstream file(argv[1]);int rows,dimension;file>>rows>>dimension;
    ZhangExactMatrix matrix(rows,ZhangExactVector(dimension));ZhangExactVector values(rows);
    for(int i=0;i<rows;++i) {int nz;file>>values[i]>>nz;for(int j=0;j<nz;++j){int c;ZhangExactInteger v;file>>c>>v;matrix[i][c]=v;}}
    if(!file)throw std::runtime_error("invalid fixture");
    auto proof=std::make_shared<ZhangIntegerDecisionProof>();proof->id="fixture";proof->originalStatement="recorded R48 00:30 domain";
    proof->observationProvenance="R48 final TRACE";proof->conditionalFailureBound=1e-6;
    for(int count:{64,128,rows}) {
        ZhangExactMatrix a(matrix.begin(),matrix.begin()+count);ZhangExactVector b(values.begin(),values.begin()+count);
        std::vector<ZhangDecisionProofs> parents(count,ZhangDecisionProofs{proof});std::vector<double> gains(count,1);
        auto start=std::chrono::steady_clock::now();
        auto now=zhangR47SelectHistorySubset(a,b,parents,{},gains,dimension,1e-3,2.5e-4);
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"new rows="<<count<<" selected="<<now.rows.size()<<" seconds="<<seconds<<std::endl;
        if(!now.valid)throw std::runtime_error("new selector rejected recorded domain");
        if(count<=128) {
            start=std::chrono::steady_clock::now();
            auto old=zhangR49SelectHistorySubsetReference(a,b,parents,{},gains,dimension,1e-3,2.5e-4);
            seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            if(now.rows!=old.rows || now.values!=old.values || now.parents!=old.parents || now.risk!=old.risk || now.selected!=old.selected || now.reasons!=old.reasons)
                throw std::runtime_error("recorded domain selection differs");
            std::cout<<"reference rows="<<count<<" seconds="<<seconds<<" identical=1"<<std::endl;
        }
    }
}
