// Included in the regional adapter TU: one common ZTD row, shared with the
// physical/iono rows' complete SAME-EPOCH covariance root.
void zhangUserRegionalTropPseudoObs(Trace& trace,KFState& state,KFMeasEntryList& entries)
{
    if (!regionalSwitch("ZHANG_USER_REGIONAL_TROP",false)) return;
    if (!zhangUserRegionalEnabled()) reject("TROP_REQUIRES_REGIONAL_JOINT_MODE");
    if (acsConfig.pppOpts.equate_tropospheres) reject("TROP_REQUIRES_PER_RECEIVER_STATE");
    const char* ledgerEnv=std::getenv("ZHANG_USER_TROP_CONSUMPTION_LEDGER");
    if (!ledgerEnv || !*ledgerEnv) reject("TROP_PERSISTENT_CONSUMPTION_LEDGER_REQUIRED");
    const char* intervalEnv=std::getenv("ZHANG_USER_TROP_INTERVAL_SECONDS");
    const std::string intervalText=intervalEnv?intervalEnv:"300";
    if (intervalText!="300" && intervalText!="900" && intervalText!="1800") reject("TROP_INTERVAL_MUST_BE_300_900_OR_1800");
    const int interval=std::stoi(intervalText);
    const long long epoch=epochSeconds(state.time);
    if (epoch%interval) return; // no held/previous-epoch TROP product
    const std::lock_guard<std::recursive_mutex> lock(regionalCacheMutex);
    const auto skip=[&](const std::string& receiver,const std::string& why) {
        trace<<"\nZHANG_USER_REGIONAL_TROP time="<<state.time.to_string(0)<<" receiver="<<receiver
            <<" status=SKIPPED reason="<<why<<" physical_rows_unchanged=1";
    };
    if (!loadEpoch(state.time)) { skip("ALL",cache.reason); return; }
    const auto& props=cache.properties;
    auto require=[&](const std::string& key,const std::string& value) {
        if (!props.contains(key) || props.at(key)!=value) throw std::runtime_error("TROP_CONTRACT:"+key);
    };
    auto numbers=[&](const std::string& key,int n) {
        std::string value=props.at(key); std::replace(value.begin(),value.end(),',',' ');
        std::istringstream in(value); Eigen::VectorXd out(n);
        for(int i=0;i<n;++i) if (!(in>>out(i)) || !std::isfinite(out(i))) throw std::runtime_error("TROP_VECTOR:"+key);
        std::string extra;if(in>>extra) throw std::runtime_error("TROP_VECTOR_EXTRA:"+key);
        return out;
    };
    Eigen::Vector3d origin; Eigen::Vector4d coefficients;
    Eigen::Matrix4d leverage; std::vector<Eigen::Vector3d> references;
    Eigen::MatrixXd selector=Eigen::MatrixXd::Zero(4,cache.mean.size());
    double spatialVariance=0,hmin=1e10,hmax=-1e10;
    try {
        require("schema","REGIONAL_JOINT_V2"); require("trop_status","NONE");
        require("trop_model",zhangtrop::modelName);
        require("trop_background","GINAN_GPT2_DRY_PLUS_WET_V1");
        require("trop_height","ELLIPSOIDAL_WGS84_NOT_ENU_UP");
        require("trop_time","SOURCE_GTIME_GPST"); require("trop_temporal_policy","WHITE_APPROX");
        require("trop_grid_sha256",zhangRegionalTropGridHash());
        require("trop_coefficients","C0,E,N,H");
        if (number(props.at("trop_horizontal_scale_m"))!=zhangtrop::horizontalScale ||
            number(props.at("trop_height_scale_m"))!=zhangtrop::heightScale) throw std::runtime_error("TROP_SCALE");
        spatialVariance=number(props.at("trop_model_variance_m2"));
        if (std::abs(spatialVariance-zhangtrop::modelVarianceFloor)>1e-15) throw std::runtime_error("TROP_FIXED_EXPERIMENT_VARIANCE");
        origin=numbers("trop_origin",3);
        const auto v=numbers("trop_leverage",16);
        for(int i=0;i<4;++i) for(int j=0;j<4;++j) leverage(i,j)=v(i*4+j);
        if ((leverage-leverage.transpose()).norm()>1e-9 ||
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d>(leverage).eigenvalues().minCoeff() < -1e-10)
            throw std::runtime_error("TROP_LEVERAGE_MATRIX");
        std::istringstream ids(props.at("trop_references")); std::string id; std::set<std::string> seen;
        while(std::getline(ids,id,',')) {
            const auto found=settings().geometry.find(id);
            if (!seen.insert(id).second || found==settings().geometry.end() || found->second.role!="REFERENCE") throw std::runtime_error("TROP_REFERENCE_ROLE");
            const Eigen::Vector3d xyz=numbers("trop_ref_"+id,3);
            if ((xyz-found->second.xyz).norm()>0.001) throw std::runtime_error("TROP_REFERENCE_GEOMETRY_MISMATCH");
            references.push_back(xyz);
            const double h=zhangtrop::geodetic(xyz).z(); hmin=std::min(hmin,h);hmax=std::max(hmax,h);
        }
        if (references.size()<8) throw std::runtime_error("TROP_REFERENCE_COUNT");
        for(int i=0;i<4;++i) {
            const int index=lookup(cache,"TROP_COEFF","REGIONAL","",std::to_string(i));
            if (index<0) throw std::runtime_error("TROP_CURRENT_COEFFICIENT_MISSING");
            selector(i,index)=1; coefficients(i)=cache.mean(index);
        }
    } catch(const std::exception& error) { skip("ALL",error.what()); return; }
    std::map<std::string,KFKey> supported;
    for(const auto& row:entries) if(row.valid && (row.obsKey.type==KF::CODE_MEAS || row.obsKey.type==KF::PHAS_MEAS))
        for(const auto& [key,coefficient]:row.designEntryMap)
            if (coefficient!=0 && key.type==KF::TROP && key.num==0 && !key.str.empty()) supported[key.str]=key;
    for(const auto& [receiver,key]:supported) {
        try {
            const auto user=settings().geometry.find(receiver);
            if (user==settings().geometry.end() || user->second.role!="USER") throw std::runtime_error("TROP_USER_ROLE");
            Eigen::Vector3d position; std::string positionOrigin;
            if (!physicalOperatingPosition(receiver,entries,state,position,positionOrigin)) throw std::runtime_error("TROP_CURRENT_POSITION_UNAVAILABLE");
            auto& options=acsConfig.getRecOpts(receiver);
            if (!options.tropModel.enable || options.tropModel.models!=std::vector<E_TropModel>{E_TropModel::GPT2}) throw std::runtime_error("TROP_USER_BACKGROUND_MISMATCH");
            InitialState init=initialStateFromConfig(options.trop,0);
            if (!init.estimate) throw std::runtime_error("TROP_STATE_NOT_ESTIMATED");
            bool exists=state.getKFValue(key,init.x)!=E_Source::NONE;
            if (!exists) for(const auto& row:entries) {
                auto it=row.usedValueMap.find(key);if(it!=row.usedValueMap.end()) {init.x=it->second;exists=true;break;}
            }
            if (!exists || !std::isfinite(init.x)) throw std::runtime_error("TROP_CURRENT_STATE_UNAVAILABLE");
            const auto b=zhangtrop::basis(origin,position);
            const double weightNorm=b.dot(leverage*b), height=zhangtrop::geodetic(position).z();
            std::vector<Eigen::Vector2d> horizontal;
            for(const auto& xyz:references) horizontal.push_back(horizontalKm(position,xyz));
            if (!zhangregional::detail::originInsideHull(horizontal)) throw std::runtime_error("TROP_OUTSIDE_REFERENCE_HULL");
            if (height<hmin-20 || height>hmax+20 || weightNorm>4 || weightNorm< -1e-10) throw std::runtime_error("TROP_HEIGHT_OR_LEVERAGE_OUTSIDE_SUPPORT");
            auto predict=[&](const Eigen::Vector3d& xyz) {
                return zhangRegionalTropBackground(state.time,xyz)+zhangtrop::basis(origin,xyz).dot(coefficients);
            };
            const double prediction=predict(position);
            Eigen::Vector3d derivative;
            for(int j=0;j<3;++j) {
                Eigen::Vector3d plus=position,minus=position;plus(j)+=1;minus(j)-=1;
                derivative(j)=(predict(plus)-predict(minus))/2;
            }
            if (!derivative.allFinite() || derivative.norm()>0.01 || prediction<0.5 || prediction>4) throw std::runtime_error("TROP_PREDICTION_OR_DERIVATIVE_INVALID");
            KFMeasEntry row(&state);
            row.obsKey.type=KF::TROP;row.obsKey.str=receiver;row.obsKey.num=0;
            row.obsKey.comment=std::string(ZHANG_USER_REGIONAL_TROP_LABEL)+cache.solution;
            row.addDsgnEntry(key,1,init);
            for(int j=0;j<3;++j) {
                auto posInit=initialStateFromConfig(options.pos,j);posInit.x=position(j);
                KFKey pk;pk.type=KF::REC_POS;pk.str=receiver;pk.num=j;
                row.addDsgnEntry(pk,-derivative(j),posInit);
            }
            row.setInnov(prediction-init.x);
            row.metaDataMap["pseudoObs"]=reinterpret_cast<void*>(true);
            row.metaDataMap[ZHANG_USER_REGIONAL_TROP_MARKER]=reinterpret_cast<void*>(true);
            const Eigen::RowVectorXd weights=b.transpose()*selector;
            const Eigen::RowVectorXd factors=weights*cache.root;
            addFactors(row,cache.noiseId,factors);
            Eigen::RowVectorXd space(1);space(0)=std::sqrt(spatialVariance);
            addFactors(row,cache.noiseId+":TROP_SPACE:"+receiver,space);
            // Persistent, atomic at-most-once ATTEMPT. A crash or QC rejection
            // may lose a constraint; it can never cause repeat information.
            // New independent experiment = new ledger. A restart reuses it.
            const fs::path ledger(ledgerEnv);fs::create_directories(ledger);
            const auto claim=ledger/zhangCheckpointSha256(receiver+":"+std::to_string(epoch));
            if (!fs::create_directory(claim)) {skip(receiver,"PRODUCT_ALREADY_ATTEMPTED");continue;}
            {std::ofstream out(claim/"source.txt",std::ios::binary);
                out<<cache.noiseId<<'\n'<<receiver<<'\n'<<epoch<<'\n';out.flush();
                if(!out) throw std::runtime_error("TROP_CONSUMPTION_LEDGER_WRITE");}
            entries.push_back(std::move(row));
            trace<<std::setprecision(17)<<"\nZHANG_USER_REGIONAL_TROP time="<<state.time.to_string(0)
                <<" receiver="<<receiver<<" status=APPENDED prediction_m="<<prediction<<" state_m="<<init.x
                <<" innovation_m="<<prediction-init.x<<" formal_variance_m2="<<factors.squaredNorm()
                <<" spatial_variance_m2="<<spatialVariance<<" weight_squared_norm="<<weightNorm
                <<" derivative_x="<<derivative.x()<<" derivative_y="<<derivative.y()<<" derivative_z="<<derivative.z()
                <<" model_order="<<props.at("trop_order")<<" geometry_origin="<<positionOrigin
                <<" shared_same_epoch_Q=1 temporal_policy=WHITE_APPROX interval_seconds="<<interval
                <<" user_truth_geometry_used=0 common_GEC_rows=1 update_pending_QC=1";
        } catch(const std::exception& error) {skip(receiver,error.what());}
    }
}
