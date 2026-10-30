#pragma once
// Include after algebra.hpp (Eigen's project plugin must be installed first).
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <cstdlib>
#include <boost/json.hpp>

inline bool zhangResetAuditEpoch(long double epoch)
{
    const char* root=std::getenv("ZHANG_RESET_AUDIT_DIRECTORY");
    const char* values=std::getenv("ZHANG_RESET_AUDIT_EPOCHS");
    if(!root || !*root || !values || !*values)return false;
    std::istringstream list(values);std::string value;
    while(std::getline(list,value,','))if(std::abs(std::stold(value)-epoch)<1e-6)return true;
    return false;
}
inline std::filesystem::path zhangResetAuditPath(long double epoch,const std::string& label)
{
    const char* root=std::getenv("ZHANG_RESET_AUDIT_DIRECTORY");
    if(!root || !*root)return {};
    auto path=std::filesystem::path(root)/std::to_string(static_cast<long long>(epoch))/label;
    std::filesystem::create_directories(path);return path;
}
template<class Derived> inline void zhangResetAuditMatrix(const std::filesystem::path& file,const Eigen::MatrixBase<Derived>& value)
{
    // Explicit dimensions and column-major float64. No rounding or unit rescale.
    const MatrixXd matrix=value;
    std::ofstream out(file,std::ios::binary);const std::uint64_t dims[2]={std::uint64_t(matrix.rows()),std::uint64_t(matrix.cols())};
    out.write(reinterpret_cast<const char*>(dims),sizeof(dims));
    out.write(reinterpret_cast<const char*>(matrix.data()),matrix.size()*sizeof(double));
    out.flush();if(!out)throw std::runtime_error("RESET_EXPERIMENT_AUDIT_WRITE_FAILED");
}
inline void zhangResetAuditState(const KFState& state,const std::string& label,bool event=false)
{
    if(!event && !zhangResetAuditEpoch(state.time.bigTime))return;
    const auto path=zhangResetAuditPath(state.time.bigTime,label);if(path.empty())return;
    zhangResetAuditMatrix(path/"x.f64",state.x);zhangResetAuditMatrix(path/"P.f64",state.P);zhangResetAuditMatrix(path/"dx.f64",state.dx);
    std::ofstream keys(path/"state_keys.txt");for(const auto& [key,index]:state.kfIndexMap)
        keys<<index<<' '<<int(key.type)<<' '<<std::quoted(key.str)<<' '<<std::quoted(key.Sat.id())<<' '<<key.num<<' '<<std::quoted(key.comment)<<'\n';
    boost::json::object metadata;for(const auto& [k,v]:state.metaDataMap)metadata[k]=v;
    std::ofstream file(path/"metadata.json");file<<boost::json::serialize(metadata)<<'\n';
}
struct ZhangResetFilterBlock
{
    int begX=0,numX=0,begH=0,numH=0;
    MatrixXd priorP,H,R,S,K,solveRightHandSide,solveFactor;VectorXd priorX,nu;
    int solverUsed=0;std::vector<int> pivots;
};
struct ZhangResetFilterCapture;
inline std::atomic<ZhangResetFilterCapture*> zhangResetActiveCapture{nullptr};
struct ZhangResetFilterCapture
{
    bool enabled=false;std::mutex mutex;std::vector<ZhangResetFilterBlock> blocks;
    explicit ZhangResetFilterCapture(const KFState& prior):enabled(zhangResetAuditEpoch(prior.time.bigTime))
    {if(enabled)zhangResetActiveCapture.store(this);}
    ~ZhangResetFilterCapture(){if(enabled)zhangResetActiveCapture.store(nullptr);}
    void capture(const KFState& prior,const KFMeas& meas,const MatrixXd& S,const MatrixXd& K,int begX,int numX,int begH,int numH,
        const MatrixXd& solveRightHandSide,const MatrixXd& solveFactor,int solverUsed,const std::vector<int>& pivots)
    {
        ZhangResetFilterBlock b;b.begX=begX;b.numX=numX;b.begH=begH;b.numH=numH;
        b.priorX=prior.x.segment(begX,numX);b.priorP=prior.P.block(begX,begX,numX,numX);
        b.H=meas.H.block(begH,begX,numH,numX);b.R=meas.R.block(begH,begH,numH,numH);b.nu=meas.V.segment(begH,numH);b.S=S;b.K=K;
        b.solveRightHandSide=solveRightHandSide;b.solveFactor=solveFactor;b.solverUsed=solverUsed;b.pivots=pivots;
        std::lock_guard<std::mutex> lock(mutex);blocks.push_back(std::move(b));
    }
    void committed(const KFState& prior,const KFState& posterior,const KFMeas& meas)
    {
        if(!enabled)return;
        zhangResetActiveCapture.store(nullptr); // AR/diagnostic solves cannot contaminate the PPP receipt.
        const auto path=zhangResetAuditPath(prior.time.bigTime,"committed_filter");
        zhangResetAuditState(prior,"committed_prior");zhangResetAuditState(posterior,"committed_posterior");
        zhangResetAuditMatrix(path/"H.f64",meas.H);zhangResetAuditMatrix(path/"R.f64",meas.R);zhangResetAuditMatrix(path/"nu.f64",meas.V);
        boost::json::array selected;std::set<std::pair<int,int>> seen;
        for(auto it=blocks.rbegin();it!=blocks.rend();++it)
        {
            const auto& b=*it;
            if(b.begX+b.numX>prior.x.size() || b.begH+b.numH>meas.V.size() || seen.contains({b.begX,b.begH}))continue;
            if(b.priorX!=prior.x.segment(b.begX,b.numX) || b.priorP!=prior.P.block(b.begX,b.begX,b.numX,b.numX) ||
               b.H!=meas.H.block(b.begH,b.begX,b.numH,b.numX) || b.R!=meas.R.block(b.begH,b.begH,b.numH,b.numH) || b.nu!=meas.V.segment(b.begH,b.numH))continue;
            seen.insert({b.begX,b.begH});const std::string id="block_"+std::to_string(b.begX)+"_"+std::to_string(b.begH);
            zhangResetAuditMatrix(path/(id+"_S.f64"),b.S);zhangResetAuditMatrix(path/(id+"_K.f64"),b.K);
            // These are the consumed HP and factorization, not a reconstruction
            // from rounded posterior matrices. Keep supplementary diagnostics
            // separate from the already frozen matrix comparison inventory.
            const auto numerical=path/"linear_solver"/id;std::filesystem::create_directories(numerical);
            zhangResetAuditMatrix(numerical/"actual_HP.f64",b.solveRightHandSide);
            zhangResetAuditMatrix(numerical/"actual_factor.f64",b.solveFactor);
            boost::json::array pivots;for(int p:b.pivots)pivots.push_back(p);
            boost::json::object solve{{"solver",b.solverUsed==1?"DPOTRF_DPOTRS":b.solverUsed==2?"DSYTRF_DSYTRS":"DGETRF_DGETRS"},
                {"solver_used",b.solverUsed},{"uplo","U"},{"factor_info",0},{"gain_solve_info",0},{"pivots",pivots},
                {"equilibration",false},{"equation","S * K_transpose = actual_HP"},{"actual_final_committed_block",true}};
            std::ofstream details(numerical/"solver.json");details<<boost::json::serialize(solve)<<'\n';
            selected.push_back(boost::json::object{{"id",id},{"begX",b.begX},{"numX",b.numX},{"begH",b.begH},{"numH",b.numH}});
        }
        boost::json::object result{{"scope","actual final committed PPP kFilter blocks"},{"candidates",blocks.size()},{"blocks",selected},
            {"postfit_substituted_for_innovation",false},{"gain_reconstructed_after_run",false}};
        std::ofstream out(path/"receipt.json");out<<boost::json::serialize(result)<<'\n';
        std::ofstream keys(path/"measurement_keys.txt");for(int i=0;i<meas.obsKeys.size();++i)
            keys<<i<<' '<<meas.obsKeys[i]<<'\n';
    }
};
