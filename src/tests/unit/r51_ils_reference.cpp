#include "ambres/GNSSambres.cpp"
#include <random>
#include <iostream>
#include <stdexcept>
int main(){
    zhangR51SingleBlock=true;GinAR_opt options;
    std::mt19937 gen(51);std::uniform_real_distribution<double> z(-.45,.45),c(-.2,.2),d(.5,1.5);
    for(int test=0;test<120;++test){
        const int n=2+test%3;GinAR_mtx m;m.zflt.resize(n);m.Dtrs.resize(n);m.Ltrs=MatrixXd::Identity(n,n);
        for(int i=0;i<n;++i){m.zflt(i)=z(gen);m.Dtrs(i)=d(gen);for(int j=0;j<i;++j)m.Ltrs(i,j)=c(gen);}
        const auto found=lambdaSearchReducedSuffix(m,n,options);
        std::vector<double> reference;VectorXd integer(n);
        auto enumerate=[&](auto&& self,int k)->void{if(k<n){for(int x=-3;x<=3;++x){integer(k)=x;self(self,k+1);}return;}
            VectorXd residual=m.zflt-integer;
            for(int i=n-1;i>=0;--i)for(int j=i+1;j<n;++j)residual(i)-=m.Ltrs(j,i)*residual(j);
            reference.push_back((residual.array().square()/m.Dtrs.array()).sum());};
        enumerate(enumerate,0);std::sort(reference.begin(),reference.end());
        if(!m.searchDiagnostic.ilsComplete || found.size()!=2)throw std::runtime_error("incomplete search");
        int i=0;for(const auto& [cost,value]:found)if(std::abs(cost-reference[i++])>1e-10)throw std::runtime_error("ILS disagrees with exhaustive oracle");
    }
    GinAR_mtx limited;limited.zflt=VectorXd::Constant(3,.2);limited.Dtrs=VectorXd::Ones(3);limited.Ltrs=MatrixXd::Identity(3,3); if(!lambdaSearchReducedSuffix(limited,3,options,1).empty() || limited.searchDiagnostic.ilsComplete)throw std::runtime_error("node cap accepted incomplete search");
    std::cout<<"R51 ILS top two agree with independent exhaustive enumeration on 120 correlated 2-4D cases\n";
}
