// Read-only postprocessing of a complete authoritative snapshot. No estimator
// changes, no user coordinates, no conditioning on an omitted product.
string extendRegionalTroposphere(ZhangRegionalJointSnapshot& s)
{
    if (!zhangRegionalTropEnabled()) return {};
    std::ostringstream metadata; metadata<<std::setprecision(17);
    metadata<<"trop_model="<<zhangtrop::modelName
        <<"\ntrop_background=GINAN_GPT2_DRY_PLUS_WET_V1"
        <<"\ntrop_height=ELLIPSOIDAL_WGS84_NOT_ENU_UP"
        <<"\ntrop_grid_sha256="<<zhangRegionalTropGridHash()
        <<"\ntrop_time=SOURCE_GTIME_GPST"
        <<"\ntrop_freshness=FINAL_COMMITTED_PPP_DESIGN_NONZERO_NOT_PER_LINK_QC"
        <<"\ntrop_horizontal_scale_m="<<zhangtrop::horizontalScale
        <<"\ntrop_height_scale_m="<<zhangtrop::heightScale
        <<"\ntrop_coefficients=C0,E,N,H\ntrop_temporal_policy=WHITE_APPROX"
        <<"\ntrop_model_variance_m2="<<zhangtrop::modelVarianceFloor
        <<"\ntrop_model_variance_calibrated=0\ntrop_loo=INTERNAL_NOT_INDEPENDENT_NETWORK_VALIDATION\n";
    vector<int> indices; vector<Eigen::Vector3d> positions;
    GTime time; time.bigTime=s.epoch;
    for (int i=0;i<s.mean.size();++i) {
        const auto& p=s.parameters[i];
        if (p.parameter!="TROP_ZTD_M" || p.sourceUpdated!=s.epoch ||
            !s.pppCommitted || !s.pppQcValid || s.mean(i)<0.5 || s.mean(i)>4 ||
            s.covariance(i,i)>0.05*0.05) continue;
        auto& options=acsConfig.getRecOpts(p.receiver);
        if (!options.tropModel.enable || options.tropModel.models!=vector<E_TropModel>{E_TropModel::GPT2})
            throw std::runtime_error("REGIONAL_TROP_REFERENCE_BACKGROUND_MISMATCH");
        const Eigen::Vector3d xyz=options.apriori_pos;
        if (!xyz.allFinite() || xyz.norm()<6e6 || xyz.norm()>7e6) continue;
        indices.push_back(i); positions.push_back(xyz);
    }
    if (indices.size()<8) { metadata<<"trop_status=INSUFFICIENT_FRESH_REFERENCES\n"; return metadata.str(); }
    Eigen::Vector3d origin=Eigen::Vector3d::Zero();
    for (const auto& xyz:positions) origin+=xyz;
    origin/=positions.size();
    const int n=indices.size(), old=s.mean.size();
    MatrixXd design(n,4), q(n,n); VectorXd residual(n);
    for (int i=0;i<n;++i) {
        design.row(i)=zhangtrop::basis(origin,positions[i]);
        residual(i)=s.mean(indices[i])-zhangRegionalTropBackground(time,positions[i]);
        for (int j=0;j<n;++j) q(i,j)=s.covariance(indices[i],indices[j]);
    }
    const auto fit=zhangtrop::fit(design,residual,q);
    metadata<<"trop_status="<<fit.reason<<"\ntrop_loo_rms_m="<<fit.looRms
        <<"\ntrop_geometry_condition="<<fit.condition<<"\ntrop_order="<<fit.order<<'\n';
    if (!fit.valid) return metadata.str();
    metadata<<"trop_origin="<<origin.x()<<','<<origin.y()<<','<<origin.z()<<"\ntrop_references=";
    for (int i=0;i<n;++i) metadata<<(i?",":"")<<s.parameters[indices[i]].receiver;
    metadata<<"\ntrop_leverage=";
    for (int i=0;i<4;++i) for(int j=0;j<4;++j) metadata<<(i||j?",":"")<<fit.leverage(i,j);
    metadata<<'\n';
    for (int i=0;i<n;++i) {
        metadata<<"trop_ref_"<<s.parameters[indices[i]].receiver<<'='<<positions[i].x()<<','<<positions[i].y()<<','<<positions[i].z()<<'\n';
        metadata<<"trop_loo_"<<s.parameters[indices[i]].receiver<<'='<<fit.leaveOneOutResiduals(i)<<','<<fit.leaveOneOutFormalVariances(i)<<'\n';
    }
    MatrixXd projection=MatrixXd::Zero(old+4,old);
    projection.topRows(old).setIdentity();
    for(int i=0;i<4;++i) for(int j=0;j<n;++j) projection(old+i,indices[j])=fit.jacobian(i,j);
    auto newRows=zhangRegionalTransformRows(projection,s.jacobian);
    VectorXd mean=projection*s.mean;
    for(int i=0;i<4;++i) {
        const double constant=fit.coefficients(i)-mean(old+i);
        newRows[old+i].push_back({-1,constant}); mean(old+i)=fit.coefficients(i);
        s.parameters.push_back({"TROP_COEFF","REGIONAL","",std::to_string(i),s.epoch});
    }
    const MatrixXd covariance=projection*s.covariance*projection.transpose();
    // Legacy block copied exactly. Matrix multiplication roundoff on that
    // unchanged block must not perturb an OFF experiment.
    MatrixXd finalQ=covariance; finalQ.topLeftCorner(old,old)=s.covariance;
    s.mean=std::move(mean);s.covariance=std::move(finalQ);s.jacobian=std::move(newRows);
    zhangRegionalValidate(s.mean,s.covariance);
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(s.covariance,Eigen::EigenvaluesOnly);
    s.minimumEigenvalue=eigen.eigenvalues().minCoeff();
    return metadata.str();
}
