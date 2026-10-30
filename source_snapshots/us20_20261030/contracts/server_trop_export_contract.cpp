#include "common/algebra.hpp"
#include "common/zhangRegionalJointProduct.hpp"
#include "common/zhangRegionalTropBackground.hpp"
#include "common/zhangRegionalTropReceipt.hpp"
#include <iostream>
#include <stdexcept>
using std::vector;using std::string;
namespace contract_export {
#include "pea/zhangRegionalTropExport.inc"
}
void check(bool ok,const char* reason) {if(!ok)throw std::runtime_error(reason);}
int main()
{
 try {
    setenv("ZHANG_REGIONAL_JOINT_EXPORT","1",1);setenv("ZHANG_REGIONAL_TROPOSPHERE","GPT2_RESIDUAL_ENH",1);
    acsConfig.gpt2grid_files={"/mnt/d/GINAN_R20/inputData/products/tables/gpt_25.grd"};readgrid(acsConfig.gpt2grid_files.front());
    acsConfig.pppOpts.equate_tropospheres=false;
    GTime time;str2time("2026 10 01 02 00 00",0,19,time,E_TimeSys::GPST);
    KFState authority;authority.time=time;authority.x=Eigen::VectorXd::Zero(2);authority.P=Eigen::MatrixXd::Identity(2,2);
    authority.lastPppTransactionCommitted=authority.lastPppQcConverged=true;
    KFKey k0;k0.type=KF::TROP;k0.str="T0";KFKey k1=k0;k1.str="T1";
    authority.kfIndexMap.clear();authority.kfIndexMap[k0]=0;authority.kfIndexMap[k1]=1;
    KFMeas measurement;measurement.time=time;measurement.H=Eigen::MatrixXd::Zero(2,2);measurement.H(0,0)=2;
    measurement.R=Eigen::MatrixXd::Identity(2,2);KFKey obs;obs.type=KF::PHAS_MEAS;obs.str="T0";
    measurement.obsKeys={obs,obs};
    const auto px=authority.x;const auto pq=authority.P;
    zhangCaptureRegionalTropDesign(authority,measurement);
    check(authority.metaDataMap.contains(zhangRegionalTropReceiptPrefix+"T0") && !authority.metaDataMap.contains(zhangRegionalTropReceiptPrefix+"T1"),"only supported final physical rows authorize freshness");
    check(authority.outageLimitMap.empty() && (authority.x-px).norm()==0 && (authority.P-pq).norm()==0,"freshness bookkeeping introduces no estimator or outage policy mutation");
    authority.lastPppTransactionCommitted=false;zhangCaptureRegionalTropDesign(authority,measurement);
    check(!authority.metaDataMap.contains(zhangRegionalTropReceiptPrefix+"T0"),"failed PPP cannot keep current freshness");
    ZhangRegionalJointSnapshot s;s.valid=s.pppCommitted=s.pppQcValid=true;s.solution="PRODUCT_FIXED";s.root="SYNTHETIC_AUTHORITY";s.epoch=time.bigTime;
    s.mean=Eigen::VectorXd::Zero(23);s.mean.head(3)<<5,1,2;
    s.parameters={{"CLOCK","","G01","NONE",s.epoch},{"PHASE","","G01","L1C",s.epoch},{"IONO_L1_M","T0","G01","L1C",s.epoch}};
    std::vector<Eigen::Vector3d> xyz;Eigen::Vector3d origin=Eigen::Vector3d::Zero();
    for(int i=0;i<20;++i) {
        Eigen::Vector3d p(6378137+10*(i%4),13000*std::cos(i*1.17),13000*std::sin(i*1.17));xyz.push_back(p);origin+=p;
        auto& opt=acsConfig.getRecOpts("T"+std::to_string(i));opt.apriori_pos=p;opt.tropModel.enable=true;opt.tropModel.models={E_TropModel::GPT2};
        s.parameters.push_back({"TROP_ZTD_M","T"+std::to_string(i),"","NONE",s.epoch});
    }
    origin/=20;Eigen::MatrixXd design(20,4);Eigen::Vector4d coeff;coeff<<.03,.004,-.002,.01;
    for(int i=0;i<20;++i) {design.row(i)=zhangtrop::basis(origin,xyz[i]);s.mean(3+i)=zhangRegionalTropBackground(time,xyz[i])+design.row(i).dot(coeff);}
    Eigen::MatrixXd root=Eigen::MatrixXd::Identity(23,23)*.004;for(int i=3;i<23;++i) root(i,0)=.0005*(1+i%3);
    s.covariance=root*root.transpose();for(int i=0;i<23;++i)s.jacobian.push_back({{i,1},{-1,0}});
    const auto original=s;
    const auto metadata=contract_export::extendRegionalTroposphere(s);
    check(s.mean.size()==27 && metadata.find("trop_status=NONE")!=string::npos,"producer appended four derived coefficients");
    check((s.mean.tail(4)-coeff).norm()<1e-11,"source GPT2 residual fit, background subtracted once");
    check((s.mean.head(23)-original.mean).norm()==0 && (s.covariance.topLeftCorner(23,23)-original.covariance).norm()==0,"legacy means and Q exactly unchanged");
    Eigen::MatrixXd j=(design.transpose()*design).ldlt().solve(design.transpose());
    Eigen::MatrixXd projection=Eigen::MatrixXd::Zero(27,23);projection.topRows(23).setIdentity();projection.bottomRightCorner(4,20)=j;
    check((s.covariance-projection*original.covariance*projection.transpose()).norm()<1e-12,"full native Q projection including CP-I-TROP");
    for(int i=0;i<27;++i) {double value=0;for(auto [k,a]:s.jacobian[i])value+=a*(k<0?1:original.mean(k));
        check(std::abs(value-s.mean(i))<1e-12,"linear map plus affine GPT2 subtraction reconstructs each published mean");}
    auto insufficient=original;for(int i=3;i<18;++i)insufficient.parameters[i].sourceUpdated-=30;
    check(contract_export::extendRegionalTroposphere(insufficient).find("INSUFFICIENT_FRESH_REFERENCES")!=string::npos && insufficient.mean.size()==23,"stale references cannot create an extension");
    setenv("ZHANG_REGIONAL_TROPOSPHERE","OFF",1);auto off=original;
    check(contract_export::extendRegionalTroposphere(off).empty() && (off.mean-original.mean).norm()==0,"OFF readonly closure");
    // Native slant-path definition remains total ZTD with wet mapping derivative.
    VectorPos pos=ecef2pos(xyz[0]);AzEl azel;azel.az=.4;azel.el=.5;
    TropStates a,b;a.zenith=2.4;b.zenith=2.401;TropMapping da,db;double va=0,vb=0;std::ostringstream trace;
    const double pa=tropModel(trace,{E_TropModel::GPT2},time,pos,azel,a,da,va);
    const double pb=tropModel(trace,{E_TropModel::GPT2},time,pos,azel,b,db,vb);
    check(std::abs((pb-pa)-.001*da.wetMap)<1e-12,"native total ZTD derivative is wetMap, not hydrostaticMap");
    std::cout<<"SERVER_TROP_EXPORT_CONTRACT_PASS actual_exporter=1 original_CP_I_unchanged=1 full_Q=1 affine_mean=1 stale_rejected=1 total_ZTD_mapping=1\n";return 0;
 } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
