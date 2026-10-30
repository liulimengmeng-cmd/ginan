// Reuse the original GEC contract's synthetic products/helpers, not its
// implementation. New assertions exercise actual production adapter + KFMeas.
#define __wrap_main old_regional_v1_main
#include "../user/source/src/tests/unit/user_regional_gec_contract.cpp"
#undef __wrap_main
#include "common/zhangRegionalTropBackground.hpp"
#include <boost/json.hpp>

extern "C" int __wrap_main()
{
 try {
    Fixture fixture;
    const auto bundle=fixture.directory/"bundles/epoch_0";
    const long long epoch=fixture.times[0].bigTime;
    acsConfig.gpt2grid_files={"/mnt/d/GINAN_R20/inputData/products/tables/gpt_25.grd"};
    readgrid(acsConfig.gpt2grid_files.front());
    const Eigen::Vector3d position(6378137,0,0);
    const double background=zhangRegionalTropBackground(fixture.times[0],position);
    Eigen::MatrixXd rr=Eigen::MatrixXd::Zero(76,76);
    rr.topLeftCorner(72,72)=fixture.covariance.llt().matrixL();
    for(int i=72;i<76;++i) {rr(i,i)=.01;rr(i,0)=.002*(i-70);rr(i,1)=.001;}
    const Eigen::MatrixXd fullQ=rr*rr.transpose();
    std::ifstream original(bundle/"regional_parameters.csv");
    std::ostringstream params;params<<original.rdbuf();params<<std::setprecision(17);
    Eigen::Vector4d coefficients;coefficients<<.03,.002,-.001,0;
    for(int i=0;i<4;++i) params<<epoch<<",PRODUCT_FIXED,"<<72+i<<",TROP_COEFF,REGIONAL,,"<<i
        <<",m,"<<coefficients(i)<<','<<epoch<<",1,1,1,0,UNAVAILABLE,-1\n";
    write(bundle/"regional_parameters.csv",params.str());
    std::ostringstream cov;cov<<std::setprecision(17)<<"gpst_seconds,solution,row_index,column_index,covariance_m2\n";
    for(int i=0;i<76;++i) for(int j=i;j<76;++j) cov<<epoch<<",PRODUCT_FIXED,"<<i<<','<<j<<','<<fullQ(i,j)<<'\n';
    write(bundle/"regional_covariance.csv",cov.str());
    std::vector<Eigen::Vector3d> refs;
    for(const auto& xy:fixture.offsets) refs.emplace_back(6378137,1000*xy.x(),1000*xy.y());
    refs.emplace_back(6378137,14000,9000);refs.emplace_back(6378137,-13000,4000);refs.emplace_back(6378137,4000,-14000);
    std::ostringstream geometry;geometry<<std::setprecision(17)<<"receiver,role,x_m,y_m,z_m\n";
    for(int i=0;i<8;++i) geometry<<"R"<<i<<",REFERENCE,"<<refs[i].x()<<','<<refs[i].y()<<','<<refs[i].z()<<'\n';
    geometry<<user<<",USER,1,2,3\n";write(fixture.directory/"geometry.csv",geometry.str());
    std::ostringstream contract;contract<<std::setprecision(17)
        <<"schema=REGIONAL_JOINT_V2\nepoch_gpst_seconds="<<epoch<<"\nsolution=PRODUCT_FIXED\nsource_numeric_root=SYNTHETIC_ZTD_ROOT"
        <<"\nparameters=76\nall_parameters_unit=m\nfull_same_epoch_cross_covariance=1"
        <<"\nionosphere=NETWORK_ESTIMABLE_GF_FIRST_SIGNAL_METRES_NOT_ABSOLUTE_TEC"
        <<"\ntrop_status=NONE\ntrop_model="<<zhangtrop::modelName
        <<"\ntrop_background=GINAN_GPT2_DRY_PLUS_WET_V1\ntrop_height=ELLIPSOIDAL_WGS84_NOT_ENU_UP"
        <<"\ntrop_time=SOURCE_GTIME_GPST\ntrop_temporal_policy=WHITE_APPROX\ntrop_grid_sha256="<<zhangRegionalTropGridHash()
        <<"\ntrop_horizontal_scale_m=20000\ntrop_height_scale_m=1000\ntrop_coefficients=C0,E,N,H"
        <<"\ntrop_model_variance_m2=0.001\ntrop_origin=6378137,0,0\ntrop_references=R0,R1,R2,R3,R4,R5,R6,R7"
        <<"\ntrop_order=3\ntrop_leverage=0.125,0,0,0,0,0.3,0,0,0,0,0.3,0,0,0,0,0\n";
    for(int i=0;i<8;++i) contract<<"trop_ref_R"<<i<<'='<<refs[i].x()<<','<<refs[i].y()<<','<<refs[i].z()<<'\n';
    write(bundle/"regional_contract.txt",contract.str());
    boost::json::object files;
    for(const char* name:{"regional_parameters.csv","regional_covariance.csv","regional_contract.txt"}) {
        std::string fail;files[name]=boost::json::object{{"bytes",fs::file_size(bundle/name)},
            {"sha256",zhangCheckpointFileSha256((bundle/name).string(),&fail)}};
    }
    const auto manifest=boost::json::serialize(boost::json::object{{"epoch",fixture.times[0].to_string(0)},
        {"ar",true},{"root","SYNTHETIC_ZTD_ROOT"},{"regional_joint",boost::json::object{{"schema","REGIONAL_JOINT_V2"},{"files",files}}}});
    write(bundle/"COMMIT.json",manifest);write(bundle/"SHA256SUMS",zhangCheckpointSha256(manifest)+"  COMMIT.json\n");
    environment("ZHANG_USER_REGIONAL_MODE","SD_PLANE_JOINT");
    environment("ZHANG_USER_REGIONAL_BUNDLE_ROOT",(fixture.directory/"bundles").string());
    environment("ZHANG_USER_REGIONAL_GEOMETRY",(fixture.directory/"geometry.csv").string());
    environment("ZHANG_USER_PRODUCT_ERROR_POLICY","WHITE");
    environment("ZHANG_USER_TROP_CONSUMPTION_LEDGER",(fixture.directory/"ledger").string());
    acsConfig.zhangPppAr.user_adapter=true;acsConfig.zhangPppAr.user_use_full_product_covariance=true;
    acsConfig.zhangPppAr.output_products=false;acsConfig.zhangPppAr.product_solution="PRODUCT_FIXED";
    acsConfig.zhangPppAr.product_filename=(fixture.directory/"products.csv").string();
    acsConfig.zhangFullRank.enable=false;acsConfig.pppOpts.ionoOpts.common_ionosphere=true;
    acsConfig.pppOpts.equate_ionospheres=false;acsConfig.pppOpts.equate_tropospheres=false;
    acsConfig.pppOpts.ionoOpts.use_if_combo=acsConfig.pppOpts.ionoOpts.use_gf_combo=false;
    auto& opts=acsConfig.getRecOpts(user);
    opts.ion_stec.estimate={true};opts.ion_stec.sigma={100};opts.ion_stec.process_noise={0};
    opts.tropModel.models={E_TropModel::GPT2};opts.tropModel.enable=true;
    opts.trop.estimate={true};opts.trop.sigma={.3};opts.trop.process_noise={.003};
    opts.pos.estimate={true,true,true};opts.pos.sigma={100,100,100};opts.pos.process_noise={0,0,0};
    KFState state;state.time=fixture.times[0];KFMeasEntryList rows;std::string identity;
    KFKey ztd;ztd.type=KF::TROP;ztd.str=user;ztd.num=0;
    Eigen::MatrixXd d=Eigen::MatrixXd::Zero(28,76);
    int ordinal=0;
    for(int sat=0;sat<9;++sat) for(int band=0;band<3;++band) {
        std::vector<double> factors;std::string id,reason;
        require(zhangUserRegionalNoiseFactors(state.time,user,SatSys(satellites[sat].c_str()),codeAt(sat,band==2?2:1),band!=0,factors,id,nullptr,&reason),reason);
        identity=id;KFMeasEntry row(&state);row.obsKey.type=band?KF::PHAS_MEAS:KF::CODE_MEAS;
        row.obsKey.str=user;row.obsKey.Sat=SatSys(satellites[sat].c_str());row.obsKey.num=band;
        InitialState init=initialStateFromConfig(opts.trop,0);init.x=background;
        row.addDsgnEntry(ztd,1,init);row.usedValueMap[ztd]=background;
        for(int j=0;j<3;++j) {KFKey k;k.type=KF::REC_POS;k.str=user;k.num=j;
            auto pi=initialStateFromConfig(opts.pos,j);pi.x=position(j);row.addDsgnEntry(k,j==0?1:.2,pi);row.usedValueMap[k]=position(j);}
        sharedFactors(row,id,factors);rows.push_back(row);
        d(ordinal,3*sat)=1;if(band)d(ordinal,3*sat+band)=-1;++ordinal;
    }
    std::ostringstream trace;
    environment("ZHANG_USER_REGIONAL_TROP","OFF");zhangUserRegionalTropPseudoObs(trace,state,rows);
    require(rows.size()==27,"TROP OFF leaves physical rows unchanged");
    environment("ZHANG_USER_REGIONAL_IONO","OFF");zhangUserRegionalIonoPseudoObs(trace,state,rows);
    require(rows.size()==27,"IONO OFF changes no physical rows");
    environment("ZHANG_USER_REGIONAL_TROP","ON");
    for(auto& row:rows) for(auto& [key,value]:row.usedValueMap)
        if(key.type==KF::REC_POS && key.num==0)value+=1000;
    zhangUserRegionalTropPseudoObs(trace,state,rows);
    require(rows.size()==27 && trace.str().find("TROP_HEIGHT_OR_LEVERAGE_OUTSIDE_SUPPORT")!=std::string::npos,"height extrapolation is not promoted to a precise constraint");
    for(auto& row:rows) for(auto& [key,value]:row.usedValueMap)
        if(key.type==KF::REC_POS && key.num==0)value-=1000;
    zhangUserRegionalTropPseudoObs(trace,state,rows);
    require(rows.size()==28,"one common GEC TROP row added: "+trace.str());
    const auto& trop=rows.back();
    require(std::abs(trop.innov-.03)<1e-12,"total ZTD residual, background counted once");
    require(trop.metaDataMap.contains(ZHANG_USER_REGIONAL_TROP_MARKER),"AR marker attached");
    for(int j=0;j<3;++j) {Eigen::Vector3d plus=position,minus=position;plus(j)+=.1;minus(j)-=.1;
        auto predict=[&](const Eigen::Vector3d& xyz){return zhangRegionalTropBackground(state.time,xyz)+zhangtrop::basis(position,xyz).dot(coefficients);};
        KFKey k;k.type=KF::REC_POS;k.str=user;k.num=j;
        require(std::abs(trop.designEntryMap.at(k)+(predict(plus)-predict(minus))/.2)<1e-8,"full background+residual ECEF derivative sign and step check");}
    d(27,72)=1;Eigen::MatrixXd expected=d*fullQ*d.transpose();expected(27,27)+=.001;
    state.stateTransition(trace,state.time);KFMeas actual(state,rows,state.time);
    near(actual.R,expected,1e-10,"actual joint KFMeas CP/TROP cross covariance oracle");
    require(expected.topRightCorner(27,1).norm()>1e-5,"nonzero correlation retained");
    zhangUserRegionalTropPseudoObs(trace,state,rows);
    require(rows.size()==28 && trace.str().find("PRODUCT_ALREADY_ATTEMPTED")!=std::string::npos,"duplicate attempt cannot add information");
    require(fs::is_directory(fixture.directory/"ledger"/zhangCheckpointSha256(user+":"+std::to_string(epoch))),"persistent atomic claim exists");
    state.time.bigTime+=600;
    zhangUserRegionalTropPseudoObs(trace,state,rows);
    require(rows.size()==28 && trace.str().find("COMMITTED_EPOCH_UNAVAILABLE")!=std::string::npos,"missing exact epoch locally skips TROP; no stale reuse");
    std::cout<<"USER_TROP_RUNTIME_PASS production_adapter=1 actual_KFMeas=1 joint_cross_Q=1 derivative=1 total_ZTD=1 independent_switches=1 duplicate_attempt=1 truth_unused=1\n";
    write(fixture.directory/"runtime_trace.txt",trace.str());
    return 0;
 } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
