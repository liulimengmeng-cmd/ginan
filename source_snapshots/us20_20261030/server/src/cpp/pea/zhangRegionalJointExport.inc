// Included in zhangPppAr.cpp after ProductCovarianceParameter. The snapshot
// lifetime is one writer transaction, never a global cache or historical KF copy.
void addRegionalIonosphereParameters(const KFState& state,
    vector<ProductCovarianceParameter>& parameters) {
    if (zhangRegionalTropEnabled()) {
        if (acsConfig.pppOpts.equate_tropospheres)
            throw std::runtime_error("REGIONAL_TROP_REQUIRES_DISTINCT_REFERENCE_ZTD");
        // Appended only after the legacy C/P rows; a readonly projection of
        // this SAME authoritative posterior, including all cross-covariance.
        for (const auto& [key,index]:state.kfIndexMap)
            if (key.type==KF::TROP && key.num==0 && !key.str.empty())
                parameters.push_back({key,SatSys{},"TROP_ZTD_M",E_ObsCode::NONE,{{index,1}}});
    }
    set<SatSys> satellites;
    for(const auto& p:parameters) satellites.insert(p.satellite);
    for(const auto& [key,index]:state.kfIndexMap) {
        if(key.type!=KF::IONO_STEC || key.str.empty() || !satellites.contains(key.Sat)) continue;
        const auto codes=acsConfig.zhangPppAr.baseline_observables.find(key.Sat.sys);
        if(codes==acsConfig.zhangPppAr.baseline_observables.end() || codes->second.empty()) continue;
        const double lambda=wavelength(key.Sat.sys,codes->second.front());
        if(!(lambda>0)) throw std::runtime_error("REGIONAL_IONO_REFERENCE_WAVELENGTH");
        const double alpha=40.3e16/std::pow(CLIGHT/lambda,2);
        // IONO_STEC is TECU; exported ionosphere and every covariance block are
        // metres at the first configured signal (GPS L1C in this experiment).
        parameters.push_back({key,key.Sat,"IONO_L1_M",codes->second.front(),{{index,alpha}}});
    }
}
ZhangRegionalJointSnapshot makeRegionalSnapshot(const KFState& state,const string& solution,
    const vector<ProductCovarianceParameter>& parameters,const VectorXd& mean,
    const MatrixXd& q,const MatrixXd& transform) {
    ZhangRegionalJointSnapshot out;
    out.epoch=static_cast<double>(state.time.bigTime);out.solution=solution;
    out.root=zhangR51NumericRoot(state.x,state.P);
    out.pppCommitted=state.lastPppTransactionCommitted;out.pppQcValid=state.lastPppQcConverged;
    out.mean=mean;out.covariance=q;
    vector<ZhangRegionalRow> raw;
    for(const auto& p:parameters) {
        // kfIndexMap keys do NOT refresh estimatedTime on observations.
        // outageLimitMap is the authoritative last-design-use receipt.
        const auto used=state.outageLimitMap.find(p.key);
        double updated=used==state.outageLimitMap.end()
            ? std::numeric_limits<double>::quiet_NaN()
            : static_cast<double>(used->first.estimatedTime.bigTime);
        if (p.parameter=="TROP_ZTD_M") {
            const auto receipt=state.metaDataMap.find(zhangRegionalTropReceiptPrefix+p.key.str);
            updated=receipt!=state.metaDataMap.end() && receipt->second==state.time.to_string(9)
                ?out.epoch:std::numeric_limits<double>::quiet_NaN();
        }
        out.parameters.push_back({p.parameter,p.key.str,p.parameter=="TROP_ZTD_M"?"":p.satellite.id(),enum_to_string(p.observable),
            updated});
        raw.push_back(p.stateTerms);
    }
    out.jacobian=zhangRegionalTransformRows(transform,raw);
    // Preserve affine integer/gauge constants as source index -1. They are
    // deterministic and add no covariance; do not silently omit them from Hx+c.
    for(int i=0;i<out.mean.size();++i) {
        double linearMean=0;
        for(auto [j,a]:out.jacobian[i]) linearMean+=a*state.x(j);
        out.jacobian[i].push_back({-1,out.mean(i)-linearMean});
    }
    zhangRegionalValidate(mean,q);
    Eigen::SelfAdjointEigenSolver<MatrixXd> eig(q,Eigen::EigenvaluesOnly);
    out.minimumEigenvalue=eig.eigenvalues().minCoeff();out.valid=true;
    return out;
}
#include "pea/zhangRegionalTropExport.inc"

void writeRegionalJointSnapshot(const ZhangRegionalJointSnapshot& input,
    const vector<ZhangInternalProduct>& products) {
    if(!zhangRegionalJointEnabled()) return;
    auto s=input;
    // A rejected private candidate must be allowed to fall back to the raw
    // FLOAT writer. If anyone tries to publish this incomplete candidate, the
    // mandatory joint-file check in OutputBundle still fails closed.
    if(!s.valid) return;
    if(s.parameters.size()!=std::size_t(s.mean.size()))
        throw std::runtime_error("REGIONAL_JOINT_SNAPSHOT_UNAVAILABLE");
    // Publication may apply a certified late integer gamma AFTER covariance
    // staging. Bind the means to the FINAL CSV rows, not the pre-gamma candidate.
    for(int i=0;i<s.mean.size();++i) {
        const auto& p=s.parameters[i];if(p.parameter=="IONO_L1_M" || p.parameter=="TROP_ZTD_M") continue;
        for(const auto& q:products) {
            if(q.solution!=s.solution || q.satellite.id()!=p.satellite ||
                (p.parameter=="PHASE" && q.observable!=string_to_enum<E_ObsCode>(p.observable))) continue;
            const double finalMean=p.parameter=="CLOCK"?q.clock_m:q.phase_m;
            const double delta=finalMean-s.mean(i);
            const SatSys sat(p.satellite.c_str());
            const double unit=p.parameter=="CLOCK"?1:wavelength(sat.sys,string_to_enum<E_ObsCode>(p.observable));
            if((p.parameter=="CLOCK" && std::abs(delta)>1e-8) ||
                std::abs(delta/unit-std::round(delta/unit))>1e-7)
                throw std::runtime_error("REGIONAL_LATE_MEAN_NOT_CERTIFIED_INTEGER_OFFSET");
            s.mean(i)=finalMean;
            for(auto& [index,a]:s.jacobian[i]) if(index==-1) a+=delta;
            break;
        }
    }
    const auto tropMetadata=extendRegionalTroposphere(s);
    const auto directory=std::filesystem::path(acsConfig.zhangPppAr.product_filename).parent_path();
    // This exporter must run inside the existing private atomic epoch bundle.
    if(!directory.filename().string().ends_with(".pending"))
        throw std::runtime_error("REGIONAL_JOINT_REQUIRES_PRIVATE_EPOCH_BUNDLE");
    std::ofstream values(directory/"regional_parameters.csv"),cov(directory/"regional_covariance.csv"),
        linear(directory/"regional_linear_map.csv"),quality(directory/"regional_quality.csv"),contract(directory/"regional_contract.txt");
    for(auto* stream:{&values,&cov,&linear,&quality,&contract}) *stream<<std::setprecision(17);
    values<<"gpst_seconds,solution,index,parameter,receiver,satellite,observable,unit,value,source_last_design_gpst,current_epoch_design_support,ppp_authority_valid,ppp_usable,pppar_usable,public_phase_segment,datum_version\n";
    cov<<"gpst_seconds,solution,row_index,column_index,covariance_m2\n";
    linear<<"gpst_seconds,solution,target_index,source_state_index,coefficient\n";
    quality<<"gpst_seconds,solution,receiver,satellite,observable,combined_m,combined_sigma_m,clock_phase_sigma_m,product_pppar_usable,spatial_error_included\n";
    std::map<string,int> clock;
    for(int i=0;i<s.mean.size();++i) if(s.parameters[i].parameter=="CLOCK") clock[s.parameters[i].satellite]=i;
    for(int i=0;i<s.mean.size();++i) {
        const auto& p=s.parameters[i];const ZhangInternalProduct* product=nullptr;
        for(const auto& q:products) if(q.solution==s.solution && q.satellite.id()==p.satellite &&
            (p.parameter=="CLOCK" || q.observable==string_to_enum<E_ObsCode>(p.observable))) {product=&q;break;}
        values<<s.epoch<<','<<s.solution<<','<<i<<','<<p.parameter<<','<<p.receiver<<','<<p.satellite<<','<<p.observable
            <<",m,"<<s.mean(i)<<','<<p.sourceUpdated<<','<<(p.sourceUpdated==s.epoch)<<','
            <<(s.pppCommitted&&s.pppQcValid)<<','<<((p.parameter=="TROP_COEFF" && s.pppCommitted && s.pppQcValid) || (product&&product->ppp_usable))<<','
            <<(product&&product->pppar_usable)<<','<<(product?product->public_phase_segment_id:"UNAVAILABLE")<<','
            <<(product?product->datum_version:-1)<<'\n';
        for(int j=i;j<s.mean.size();++j) cov<<s.epoch<<','<<s.solution<<','<<i<<','<<j<<','<<s.covariance(i,j)<<'\n';
        for(auto [j,a]:s.jacobian[i]) linear<<s.epoch<<','<<s.solution<<','<<i<<','<<j<<','<<a<<'\n';
        if(p.parameter!="IONO_L1_M" || !clock.contains(p.satellite)) continue;
        const int c=clock.at(p.satellite);
        const SatSys sat(p.satellite.c_str());
        const auto base=string_to_enum<E_ObsCode>(p.observable);
        for(int b=0;b<s.mean.size();++b) {
            const auto& bp=s.parameters[b];
            if(bp.parameter!="PHASE" || bp.satellite!=p.satellite) continue;
            const double mu=std::pow(wavelength(sat.sys,string_to_enum<E_ObsCode>(bp.observable))/wavelength(sat.sys,base),2);
            VectorXd row=VectorXd::Zero(s.mean.size());row(c)=1;row(b)=-1;row(i)=mu;
            const double v=(row.transpose()*s.covariance*row)(0,0);
            row(i)=0;const double cp=(row.transpose()*s.covariance*row)(0,0);row(i)=mu;
            if(v < -1e-8*std::max(1.0,s.covariance.diagonal().cwiseAbs().maxCoeff()))
                throw std::runtime_error("REGIONAL_COMBINATION_NEGATIVE_VARIANCE");
            bool pairUsable=false;
            for(const auto& q:products) if(q.solution==s.solution && q.satellite.id()==bp.satellite &&
                q.observable==string_to_enum<E_ObsCode>(bp.observable)) pairUsable=q.pppar_usable;
            quality<<s.epoch<<','<<s.solution<<','<<p.receiver<<','<<p.satellite<<','<<bp.observable<<','
                <<row.dot(s.mean)<<','<<std::sqrt(std::max(0.0,v))<<','<<std::sqrt(std::max(0.0,cp))<<','
                <<pairUsable<<",0\n";
        }
    }
    contract<<"schema="<<(zhangRegionalTropEnabled()?"REGIONAL_JOINT_V2":"REGIONAL_JOINT_V1")<<"\nepoch_gpst_seconds="<<s.epoch<<"\nsolution="<<s.solution
        <<"\nsource_numeric_root="<<s.root<<"\nparameters="<<s.mean.size()<<"\nminimum_eigenvalue_m2="<<s.minimumEigenvalue
        <<"\nall_parameters_unit=m\nionosphere=NETWORK_ESTIMABLE_GF_FIRST_SIGNAL_METRES_NOT_ABSOLUTE_TEC"
        <<"\nfull_same_epoch_cross_covariance=1\nspatial_prediction_error=NOT_PROVIDED"
        <<"\nlinear_map_source_minus_one=AFFINE_CONSTANT_METRES"
        <<"\ntemporal_cross_covariance=NOT_PROVIDED\nquality_location=REFERENCE_LINK_NOT_USER_LOCATION"
        <<"\ninteger_authority=ONLY_FINAL_PRODUCT_ROWS_NOT_STATE_PRESENCE\nnis_evaluations=0\n";
    contract<<tropMetadata;
    contract<<"freshness_receipt=LAST_DESIGN_USE_NOT_PER_LINK_QC_ACCEPTANCE\n";
    for(auto* stream:{&values,&cov,&linear,&quality,&contract}) {
        stream->flush();if(!*stream) throw std::runtime_error("REGIONAL_JOINT_WRITE_FAILED");
    }
}

void writeRegionalSourceAudit(const KFState& state) {
    const char* setting=std::getenv("ZHANG_REGIONAL_AUDIT_EPOCHS");
    if(!zhangRegionalJointEnabled() || !setting) return;
    const auto stamp=std::to_string(static_cast<long long>(std::llround(state.time.bigTime)));
    if((","+std::string(setting)+",").find(","+stamp+",")==std::string::npos) return;
    const auto directory=std::filesystem::path(acsConfig.zhangPppAr.product_filename).parent_path();
    if(!directory.filename().string().ends_with(".pending")) throw std::runtime_error("REGIONAL_AUDIT_NOT_PRIVATE");
    {std::ofstream out(directory/"regional_source_mean.bin",std::ios::binary);
     out.write(reinterpret_cast<const char*>(state.x.data()),state.x.size()*sizeof(double));
     out.flush();if(!out) throw std::runtime_error("REGIONAL_SOURCE_MEAN_WRITE");}
    {std::ofstream out(directory/"regional_source_covariance.bin",std::ios::binary);
     out.write(reinterpret_cast<const char*>(state.P.data()),state.P.size()*sizeof(double));
     out.flush();if(!out) throw std::runtime_error("REGIONAL_SOURCE_COVARIANCE_WRITE");}
    {std::ofstream out(directory/"regional_source_keys.csv");
     out<<"index,type,receiver,satellite,number\n";
     for(const auto& [key,index]:state.kfIndexMap) out<<index<<','<<enum_to_string(key.type)<<','<<key.str<<','<<key.Sat.id()<<','<<key.num<<'\n';
     out.flush();if(!out) throw std::runtime_error("REGIONAL_SOURCE_KEYS_WRITE");}
    std::ofstream contract(directory/"regional_contract.txt",std::ios::app);
    contract<<"source_audit_dimension="<<state.x.size()<<"\nsource_audit_matrix_order=column_major\n";
    for(const auto* name:{"regional_source_mean.bin","regional_source_covariance.bin","regional_source_keys.csv"})
        contract<<"source_audit_file="<<name<<" sha256="<<zhangR51BundleFileHash(directory/name)
            <<" bytes="<<std::filesystem::file_size(directory/name)<<'\n';
    contract.flush();if(!contract) throw std::runtime_error("REGIONAL_SOURCE_AUDIT_CONTRACT_WRITE");
}
