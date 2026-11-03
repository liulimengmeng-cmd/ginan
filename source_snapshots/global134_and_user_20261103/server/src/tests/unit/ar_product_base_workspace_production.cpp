// Production-linked P1/P2 contract. The golden conditioner core is the frozen
// original v6 body in its own namespace, never the candidate's core functions.
#include "common/acsConfig.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangProductBasePosteriorWorkspace.hpp"
#include "common/zhangR49PosteriorWorkspace.hpp"
#include "ambres/GNSSambres.hpp"
#include "fixtures/ar_p1_original_integer_conditioner.hpp"
#include <cstring>
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

using SharedDomain = ZhangProductBasePosteriorWorkspace<ZhangProductBaseColumnDescriptor>;
extern std::vector<ZhangProductBaseColumnDescriptor> zhangProductBaseOrderedColumns(const std::map<int,KFKey>&);
extern ZhangProductBasePosteriorIdentity zhangProductBaseSourceIdentityForAudit(
    const KFState&,const GinAR_mtx&);
extern bool zhangBuildProductBaseWorkspaceForAudit(Trace&,const KFState&,const KFState&,
    const GinAR_mtx&,E_Sys,const std::vector<GinAR_mtx>&,SharedDomain*,GinAR_mtx&,std::string&);
static int checks=0;
static void require(bool ok,const char* why)
{++checks;if(!ok)throw std::runtime_error(why);}
static bool rawEqual(const MatrixXd& a,const MatrixXd& b)
{return a.rows()==b.rows() && a.cols()==b.cols() &&
    (!a.size() || std::memcmp(a.data(),b.data(),std::size_t(a.size())*sizeof(double))==0);}
static bool rawEqual(const VectorXd& a,const VectorXd& b)
{return a.size()==b.size() && (!a.size() || std::memcmp(a.data(),b.data(),std::size_t(a.size())*sizeof(double))==0);}
static void sameResult(const ZhangPosteriorEffectiveConditioningResult& actual,
    const zhangP1Original::ZhangPosteriorEffectiveConditioningResult& golden)
{
    require(actual.valid==golden.valid,"validity differs from original conditioner");
    require(actual.failureReason==golden.failureReason,"failure reason differs from original conditioner");
    require(actual.conditioned==golden.conditioned,"conditioned flag differs");
    require(actual.conditionerRank==golden.conditionerRank,"retained rank differs");
    require(actual.deterministicRows==golden.deterministicRows,"deterministic row count differs");
    require(actual.maximumDeterministicResidual==golden.maximumDeterministicResidual,"deterministic residual differs");
    require(rawEqual(actual.mean,golden.mean),"complete mean raw bytes differ");
    require(rawEqual(actual.covariance,golden.covariance),"complete covariance raw bytes differ");
}
static ZhangDecisionProofPtr proof(const std::string& id)
{
    auto p=std::make_shared<ZhangIntegerDecisionProof>();
    p->id=id;p->originalStatement="unit CURRENT accepted row";
    p->observationProvenance="synthetic fixture only";p->conditionalFailureBound=1;
    return p;
}
struct Fixture
{
    KFState state;
    std::shared_ptr<const GinAR_mtx> source;
    std::vector<std::pair<int,KFKey>> order;
    ZhangProductBasePosteriorIdentity identity;
    Fixture()
    {
        state.time.bigTime=1405216800;state.factorCommitSequence=7;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="P1_FROZEN_GEC_UNIT";
        state.kfIndexMap.clear();state.x.resize(9);state.dx.resize(9);
        const char* sats[]={"G02","G03","G05","E02","E03","E05","C20","C21","C22"};
        for(int i=0;i<9;++i)
        {
            KFKey key;key.type=KF::AMBIGUITY;key.str="P1";key.Sat=SatSys(sats[i]);key.num=int(E_ObsCode::L1C);
            state.kfIndexMap[key]=i;order.push_back({i,key});
            state.x[i]=double(i+1)/8;state.dx[i]=double(9-i)/32;
        }
        VectorXd u(9);for(int i=0;i<9;++i)u[i]=double(i+1)/32;
        state.P=MatrixXd::Identity(9,9)+u*u.transpose();
        GinAR_mtx original;for(const auto& [i,key]:order)original.ambmap[i]=key;
        original.aflt=state.x;original.Paflt=state.P;
        source=std::make_shared<const GinAR_mtx>(std::move(original));
        identity=zhangProductBaseSourceIdentityForAudit(state,*source);
        require(!identity.owner.empty() && !identity.root.empty(),"fixture runtime/root is not bound");
    }
    std::shared_ptr<const MatrixXd> covariance() const
    {return {source,&source->Paflt};}
    GinAR_mtx receipt(int column) const
    {
        GinAR_mtx r;r.ambmap=source->ambmap;r.Ztrs=MatrixXd::Zero(1,9);r.Ztrs(0,column)=1;
        r.zfix=VectorXd::Constant(1,2);r.decisionProofs={proof("P1_ROW_"+std::to_string(column))};return r;
    }
};
static void originalConditionCases()
{
    std::vector<MatrixXd> covariances;
    MatrixXd p=MatrixXd::Identity(6,6);p(0,3)=p(3,0)=0.25;p(1,4)=p(4,1)=-0.125;
    covariances.push_back(p);
    p=MatrixXd::Identity(6,6);p(5,5)=0;covariances.push_back(p);
    p=MatrixXd::Identity(6,6);p(4,4)=1e-18;p(5,5)=0;covariances.push_back(p);
    p=MatrixXd::Identity(6,6);p(4,4)=-1e-4;covariances.push_back(p);
    p=MatrixXd::Identity(6,6);p(2,2)=std::numeric_limits<double>::quiet_NaN();covariances.push_back(p);
    for(std::size_t fixture=0;fixture<covariances.size();++fixture)
    {
        const auto frozen=std::make_shared<const MatrixXd>(covariances[fixture]);
        ZhangProductBasePosteriorIdentity identity{"OWNER","EPOCH","FROZEN_ROOT",7};
        ZhangProductBasePosteriorWorkspace<int> cache(identity,{0,1,2,3,4,5},frozen);
        const VectorXd mean=VectorXd::Zero(6);
        for(int kind=0;kind<4;++kind)
        {
            MatrixXd a=MatrixXd::Zero(2,6);a(0,0)=1;a(1,5)=1;
            if(kind==2)a.row(1)=a.row(0);
            if(kind==3)a.row(1)=2*a.row(0);
            VectorXd rhs=VectorXd::Zero(2);
            if(kind==1)rhs[1]=1;
            if(kind==3)rhs[1]=1; // dependent, inconsistent affine relation
            const auto golden=zhangP1Original::zhangConditionPosteriorEffectiveIntegers(mean,*frozen,a,rhs);
            const auto workspace=cache.acquire(identity,{0,1,2,3,4,5},*frozen);
            const auto actual=zhangConditionPosteriorEffectiveIntegers(mean,*frozen,a,rhs,1e-7,workspace,nullptr,true);
            sameResult(actual,golden);
        }
        require(cache.builds==(frozen->allFinite()?1U:0U),
            "fixed domain rebuilt or nonfinite source reached cache decomposition");
    }
    // A fully stochastic workspace must still execute the old conditional audit.
    std::ostringstream trace;
    {
        ZhangPhaseTimer timer(trace,"P1_STRICT_CHECK_TRACE");
        const MatrixXd p=MatrixXd::Identity(3,3),a=MatrixXd::Identity(3,3);
        const auto workspace=zhangBuildPosteriorEffectiveWorkspace(p);
        auto result=zhangConditionPosteriorEffectiveIntegers(VectorXd::Zero(3),p,a,
            VectorXd::Zero(3),1e-7,&workspace,nullptr,true);
        require(result.valid,"strict stochastic fixture failed");
    }
    require(trace.str().find("POSTERIOR_EFFECTIVE_CONDITIONAL_RESIDUAL")!=std::string::npos,
        "shared factor bypassed the original conditional residual audit");
}
static void domainInvalidation()
{
    auto frozen=std::make_shared<const MatrixXd>(MatrixXd::Identity(3,3));
    const ZhangProductBasePosteriorIdentity id{"OWNER","EPOCH","ROOT",7};
    ZhangProductBasePosteriorWorkspace<int> cache(id,{0,1,2},frozen);
    require(cache.acquire(id,{0,1,2},*frozen)!=nullptr,"valid domain not admitted");
    for(int change=0;change<4;++change)
    {
        auto bad=id;
        if(change==0)bad.owner+="_CLONE";
        if(change==1)bad.epoch+="_NEXT";
        if(change==2)bad.root+="_NEW";
        if(change==3)++bad.factorSequence;
        require(!cache.acquire(bad,{0,1,2},*frozen),"stale owner/epoch/root/sequence reused");
    }
    require(!cache.acquire(id,{1,0,2},*frozen),"same-dimension column reorder reused");
    MatrixXd changed=*frozen;changed(0,1)=changed(1,0)=0.125;
    require(!cache.acquire(id,{0,1,2},changed),"P mutation reused despite unchanged declared root");
    require(!cache.acquire(id,{0,1,2,3},MatrixXd::Identity(4,4)),"expanded nuisance domain reused");
    require(cache.builds==1 && cache.hits==0 && cache.fallbacks==7,"cache diagnostics wrong");
    // A reordered posterior must use its own original factor/checks, not the
    // old same-dimension basis. Compare the complete fallback Gaussian result.
    MatrixXd changedOrder=MatrixXd::Identity(3,3);
    changedOrder(0,2)=changedOrder(2,0)=0.25;
    MatrixXd a=MatrixXd::Zero(1,3);a(0,2)=1;
    const VectorXd mean=VectorXd::Zero(3),rhs=VectorXd::Ones(1);
    const auto noReuse=cache.acquire(id,{2,0,1},changedOrder);
    require(!noReuse,"reordered Gaussian used stale effective root");
    sameResult(zhangConditionPosteriorEffectiveIntegers(mean,changedOrder,a,rhs,1e-7,noReuse,nullptr,true),
        zhangP1Original::zhangConditionPosteriorEffectiveIntegers(mean,changedOrder,a,rhs));
    frozen.reset(); // cache keeps the exact immutable source alive
    require(cache.sourceCovarianceBytes()==9*sizeof(double),"source lifetime not retained");
}
static void descriptorCommentBinding()
{
    const auto p=std::make_shared<const MatrixXd>(MatrixXd::Identity(2,2));
    const ZhangProductBasePosteriorIdentity id{"OWNER","EPOCH","ORIGINAL_ROOT_TOKEN",7};
    const std::vector<ZhangProductBaseColumnDescriptor> columns{
        {0,1,2,3,4,"RX","old coordinate comment"},{1,1,2,5,6,"RX",""}};
    ZhangProductBasePosteriorWorkspace<ZhangProductBaseColumnDescriptor> cache(id,columns,p);
    require(cache.acquire(id,columns,*p,true)!=nullptr,"exact descriptor source rejected");
    auto changed=columns;changed[0].comment="new coordinate comment";
    require(!cache.acquire(id,changed,*p,true),"comment-only descriptor change reused the factor");
    const MatrixXd equalButDifferentObject=*p;
    require(!cache.acquire(id,columns,equalButDifferentObject,true),"unproven equal P object used frozen root token");
    require(cache.frozenSourceRoot()==id.root && cache.builds==1 && cache.hits==0 && cache.fallbacks==2,
        "frozen root token or descriptor diagnostics changed");
}
static void productionBaseAndCrossQ()
{
    Fixture f;const VectorXd x=f.state.x,dx=f.state.dx;const MatrixXd p=f.state.P;
    SharedDomain cache(f.identity,zhangProductBaseOrderedColumns(f.source->ambmap),f.covariance());
    const E_Sys systems[]={E_Sys::GPS,E_Sys::GAL,E_Sys::BDS};
    std::ostringstream trace;
    for(int s=0;s<3;++s)
    {
        const int column=3*s;GinAR_mtx actual;std::string reason;
        const auto r=f.receipt(column);
        require(zhangBuildProductBaseWorkspaceForAudit(trace,f.state,f.state,*f.source,
            systems[s],{r},&cache,actual,reason),"actual production base rejected");
        MatrixXd a=MatrixXd::Zero(1,9);a(0,column)=1;
        const auto golden=zhangP1Original::zhangConditionPosteriorEffectiveIntegers(
            f.source->aflt,f.source->Paflt,a,r.zfix);
        require(golden.valid,"golden base fixture rejected");
        require(rawEqual(actual.aflt,golden.mean),"production base mean differs from frozen old core");
        require(rawEqual(actual.Paflt,golden.covariance),"production base full covariance differs");
        require(actual.Paflt.rows()==9 && actual.Paflt(1,4)!=0,"cross-system covariance was dropped");
        require(rawEqual(f.state.x,x)&&rawEqual(f.state.P,p)&&rawEqual(f.state.dx,dx),"production base mutated FLOAT x/P/dx");
    }
    require(cache.builds==1 && cache.requests==3 && cache.hits==2,"three system bases did not share one decomposition");
    require(trace.str().find("preserve_uncached_checks=1")!=std::string::npos,"production strict gate not observed");
    auto changedColumns=f.source->ambmap;changedColumns.at(0).comment="comment-only coordinate change";
    require(!cache.acquire(f.identity,zhangProductBaseOrderedColumns(changedColumns),f.source->Paflt,true),
        "production column descriptor ignored a comment-only change");
    require(cache.fallbacks==1 && cache.builds==1 && cache.hits==2,"comment miss affected existing factor");
    auto inconsistent=f.receipt(0);
    inconsistent.Ztrs=MatrixXd::Zero(2,9);inconsistent.Ztrs(0,0)=inconsistent.Ztrs(1,0)=1;
    inconsistent.zfix.resize(2);inconsistent.zfix<<2,3;
    GinAR_mtx rejected;std::string rejectedReason;
    require(!zhangBuildProductBaseWorkspaceForAudit(trace,f.state,f.state,*f.source,
        E_Sys::GPS,{inconsistent},&cache,rejected,rejectedReason),"inconsistent base was accepted");
    require(rejectedReason=="NETWORK_WL_HELD_UNION_INCONSISTENT","inconsistent base failure changed");
    require(rawEqual(f.state.x,x)&&rawEqual(f.state.P,p)&&rawEqual(f.state.dx,dx),"failed base polluted original FLOAT");
    // Real base expansion: the receipt contains an active nuisance key absent
    // from the eight-column search posterior. The complete state still has it.
    GinAR_mtx narrow=*f.source;narrow.ambmap.erase(8);narrow.aflt=f.source->aflt.head(8);
    narrow.Paflt=f.source->Paflt.topLeftCorner(8,8);
    const auto narrowOwner=std::make_shared<const GinAR_mtx>(narrow);
    const auto identity=zhangProductBaseSourceIdentityForAudit(f.state,*narrowOwner);
    const auto order=zhangProductBaseOrderedColumns(narrow.ambmap);
    SharedDomain narrowCache(identity,order,{narrowOwner,&narrowOwner->Paflt});
    GinAR_mtx r=f.receipt(6);r.Ztrs(0,8)=2;r.zfix[0]=3;
    GinAR_mtx actual;std::string reason;
    require(zhangBuildProductBaseWorkspaceForAudit(trace,f.state,f.state,*narrowOwner,
        E_Sys::BDS,{r},&narrowCache,actual,reason),"production nuisance expansion rejected");
    MatrixXd a=MatrixXd::Zero(1,9);a(0,6)=1;a(0,8)=2;
    const auto golden=zhangP1Original::zhangConditionPosteriorEffectiveIntegers(
        f.source->aflt,f.source->Paflt,a,r.zfix);
    require(golden.valid,"golden nuisance fixture rejected");
    require(rawEqual(actual.aflt,VectorXd(golden.mean.head(8))),"nuisance expansion mean changed");
    require(rawEqual(actual.Paflt,MatrixXd(golden.covariance.topLeftCorner(8,8))),"nuisance full-Q projection changed");
    require(narrowCache.fallbacks==1 && narrowCache.builds==0,"expanded evidence did not use uncached fallback");
    require(rawEqual(f.state.x,x)&&rawEqual(f.state.P,p)&&rawEqual(f.state.dx,dx),"nuisance fallback mutated FLOAT");
}
static void symmetryAndAtomicWorkspace()
{
    std::vector<MatrixXd> matrices;
    MatrixXd a(3,3);a<<-0.0,0.25,-0.125,0.5,-0.0,0.75,-0.5,1,0;
    matrices.push_back(a);
    a=MatrixXd::Constant(3,3,std::numeric_limits<double>::denorm_min());a(0,1)=-0.0;matrices.push_back(a);
    a=MatrixXd::Constant(3,3,std::numeric_limits<double>::max()*0.75);matrices.push_back(a);
    a=MatrixXd::Identity(3,3);a(1,2)=std::numeric_limits<double>::quiet_NaN();matrices.push_back(a);
    for(const auto& input:matrices)
    {
        const MatrixXd golden=(0.5*(input+input.transpose())).eval();
        MatrixXd actual=input;zhangR49SymmetrizeCovarianceInPlace(actual);
        require(golden.allFinite()==actual.allFinite(),"symmetry reduced finite/overflow gate");
        if(golden.allFinite())require(rawEqual(actual,golden),"symmetry signed-zero/subnormal raw bytes changed");
    }
    MatrixXd root=MatrixXd::Identity(3,3);root(0,1)=0.125;
    MatrixXd golden=root*root.transpose();golden=(0.5*(golden+golden.transpose())).eval();
    ZhangR49PosteriorWorkspace workspace;
    const std::vector<std::string> order{"G","E","C"};
    require(workspace.formGram(root,order,"epoch"),"new Gram rejected valid root");
    require(rawEqual(workspace.covariance,golden),"Gram full covariance differs");
    require(workspace.factor(golden,order,"epoch") && workspace.hits==1,"conditioned-root reuse was lost");
    const auto previousDecompositions=workspace.decompositions;
    require(workspace.factor(golden,{"E","G","C"},"epoch") && workspace.hits==1 &&
        workspace.decompositions==previousDecompositions+1,"reordered cache got an exact-order hit");
    root=MatrixXd::Identity(3,3);root(0,0)=std::sqrt(std::numeric_limits<double>::max()*0.75);
    MatrixXd oldGram=root*root.transpose();oldGram=(0.5*(oldGram+oldGram.transpose())).eval();
    require(!oldGram.allFinite(),"overflow fixture did not exercise old diagonal overflow");
    require(!workspace.formGram(root,order,"epoch"),"new Gram accepted old overflowing covariance");
    root(0,0)=std::numeric_limits<double>::quiet_NaN();
    require(!workspace.formGram(root,order,"epoch"),"new Gram accepted nonfinite root");
}
int main()
{
    originalConditionCases();domainInvalidation();descriptorCommentBinding();productionBaseAndCrossQ();symmetryAndAtomicWorkspace();
    std::cout<<"AR_PRODUCT_BASE_WORKSPACE_CONTRACT_PASSED checks="<<checks
        <<" complete_x_P_dx=1 independent_frozen_core=1 current_integer_policy_changed=0\n";
}
