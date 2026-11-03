#include "common/zhangCheckpoint.hpp"
#include <algorithm>
#include <cmath>
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

static ZhangCheckpointBundle load(const std::filesystem::path& path,const std::string& binary)
{
    std::ifstream file(path.parent_path()/"checkpoint_manifest.json");
    std::ostringstream bytes;bytes<<file.rdbuf();
    const auto document=boost::json::parse(bytes.str());
    const auto& object=document.as_object();
    auto field=[&](const char* name){
        const auto& value=object.at(name).as_string();
        return std::string(value.data(),value.size());
    };
    ZhangCheckpointExpectations expected;
    expected.experimentMode=field("experiment_mode");
    expected.binarySha256=binary;
    expected.configSha256=zhangCheckpointSha256(field("config_text"));
    expected.inputManifestSha256=zhangCheckpointSha256(field("input_manifest_text"));
    expected.platformFingerprint=field("platform");
    expected.compilerFingerprint=field("compiler");
    expected.linearAlgebraFingerprint=field("linear_algebra");
    expected.endianness=field("endianness");
    if(expected.binarySha256!=field("binary_sha256") || expected.configSha256!=field("config_sha256") ||
       expected.inputManifestSha256!=field("input_manifest_sha256"))
        throw std::runtime_error("Manifest identity/hash mismatch");
    const auto envelope=verifyZhangCheckpointBundleEnvelope(path.string());
    if(!envelope.valid)throw std::runtime_error(envelope.failureReason);
    ZhangCheckpointBundle bundle;
    const auto decoded=readZhangCheckpointBundle(path.string(),expected,bundle);
    if(!decoded.valid)throw std::runtime_error(decoded.failureReason);
    if(bundle.sections.empty() || !bundle.kfCore.x.allFinite() || !bundle.kfCore.P.allFinite() ||
       bundle.kfCore.P.rows()!=bundle.kfCore.x.size() || bundle.kfCore.P.cols()!=bundle.kfCore.x.size())
        throw std::runtime_error("Decoded KF/sections invalid");
    std::cout<<"PASS checkpoint envelope, producer identity, config/input hashes, all sections and full KF decode"
             <<" dimension="<<bundle.kfCore.x.size()<<" sections="<<bundle.sections.size()
             <<" epoch="<<bundle.manifest.epoch<<" payload_sha256="<<decoded.payloadSha256<<'\n';
    return bundle;
}

int main(int argc,char** argv)
{
    if(argc!=3 && argc!=5)
        throw std::runtime_error("Usage: new.bundle new-producer-sha256 [reference.bundle reference-producer-sha256]");
    const auto current=load(argv[1],argv[2]);
    if(argc==5) {
        const auto reference=load(argv[3],argv[4]);
        auto order=[](const ZhangCheckpointKfCore& core) {
            std::vector<std::string> labels(core.x.size());
            for(const auto& [key,index]:core.kfIndexMap)
                labels.at(index)=std::to_string(static_cast<int>(key.type))+":"+
                    key.satellite.id()+":"+key.receiver+":"+std::to_string(key.number);
            return labels;
        };
        if(current.manifest.epoch!=reference.manifest.epoch || order(current.kfCore)!=order(reference.kfCore))
            throw std::runtime_error("Cannot compare different epochs/state identities/orders");
        const auto& a=current.kfCore;const auto& b=reference.kfCore;
        const double meanMaximum=(a.x-b.x).cwiseAbs().maxCoeff();
        const double covarianceMaximum=(a.P-b.P).cwiseAbs().maxCoeff();
        const double meanRelative=(a.x-b.x).norm()/std::max(1.,b.x.norm());
        const double covarianceRelative=(a.P-b.P).norm()/std::max(1.,b.P.norm());
        double meanComponentScaled=0,covarianceEntryScaled=0;
        for(int i=0;i<a.x.size();++i)
            meanComponentScaled=std::max(meanComponentScaled,
                std::abs(a.x(i)-b.x(i))/std::max(1.,std::abs(b.x(i))));
        for(int j=0;j<a.P.cols();++j)for(int i=0;i<a.P.rows();++i)
            covarianceEntryScaled=std::max(covarianceEntryScaled,
                std::abs(a.P(i,j)-b.P(i,j))/std::max(1.,
                    std::sqrt(std::abs(b.P(i,i)))*std::sqrt(std::abs(b.P(j,j)))));
        std::cout<<std::setprecision(17)<<"REAL_FLOAT_FULL_COMPARE entries="<<a.P.size()
            <<" mean_maximum_native="<<meanMaximum<<" covariance_maximum_native="<<covarianceMaximum
            <<" mean_relative_norm="<<meanRelative<<" covariance_relative_frobenius="<<covarianceRelative
            <<" mean_component_scaled="<<meanComponentScaled
            <<" covariance_entry_scaled="<<covarianceEntryScaled
            <<" state_order_identical=1 mixed_native_units=1\n"<<std::flush;
        if(meanRelative>2e-10 || covarianceRelative>2e-11 ||
           meanComponentScaled>2e-10 || covarianceEntryScaled>2e-11)
            throw std::runtime_error("Real full FLOAT mean/covariance relative equivalence failed");
    }
}
