#pragma once
#include "common/acsConfig.hpp"
#include "common/zhangR51Integrity.hpp"
#include "common/zhangR51CommittedBundle.hpp"
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include "common/zhangCheckpoint.hpp"
#include "common/zhangMultiGnssProductManifest.hpp"

/** A validated, renamed epoch directory is the only publication authority.
 * Aggregate CSVs are not updated because two separate files cannot be
 * atomically appended for an external user. Rejected private directories
 * remain available for audit. */
struct ZhangR51OutputBundle {
    std::string product,covariance,root,epoch,route;
    std::filesystem::path pending,committed;
    bool restored=false;
    ZhangR51OutputBundle(const std::string& time,const std::string& name,const std::string& rootId)
        : product(acsConfig.zhangPppAr.product_filename),covariance(acsConfig.zhangPppAr.product_covariance_filename),root(rootId),epoch(time),route(name) {
        auto base=std::filesystem::path(product+".epochs");
        std::filesystem::create_directories(base);
        std::string tag=time;for(auto& c:tag)if(c==':' || c==' ')c='_';
        static std::uint64_t sequence=0;
        tag+="_"+name+"_"+std::to_string(++sequence);
        pending=base/(tag+".pending");committed=base/tag;
        if(!std::filesystem::create_directory(pending))throw std::runtime_error("R51_BUNDLE_EXISTS");
        acsConfig.zhangPppAr.product_filename=(pending/"products.csv").string();
        acsConfig.zhangPppAr.product_covariance_filename=(pending/"covariance.csv").string();
    }
    void restore() {if(!restored){acsConfig.zhangPppAr.product_filename=product;acsConfig.zhangPppAr.product_covariance_filename=covariance;restored=true;}}
    ~ZhangR51OutputBundle(){restore();}
    static std::string fileHash(const std::filesystem::path& path) {
        return zhangR51BundleFileHash(path);
    }
    void publish(const std::string& certificate,bool ar) {
        restore();
        if(!std::filesystem::exists(pending/"products.csv"))throw std::runtime_error("R51_PRODUCT_FILE_MISSING");
        if(ar && !std::filesystem::exists(pending/"covariance.csv"))throw std::runtime_error("R51_COVARIANCE_FILE_MISSING");
        {std::ofstream out(pending/"certificate.txt",std::ios::binary);out<<certificate;out.flush();if(!out)throw std::runtime_error("R51_CERTIFICATE_WRITE_FAILED");}
        {boost::json::object manifest;
            manifest["epoch"]=epoch;manifest["route"]=route;
            manifest["root"]=root;manifest["ar"]=ar;
            ZhangMultiGnssManifestSignals systemSignals;
            for (const auto& [system, signals] : acsConfig.zhangPppAr.baseline_observables)
            {
                const auto enabled = acsConfig.process_sys.find(system);
                if (enabled == acsConfig.process_sys.end() || !enabled->second) continue;
                auto& codes = systemSignals[enum_to_string(system)];
                for (const auto code : signals) codes.push_back(enum_to_string(code));
            }
            if (systemSignals.size() > 1)
            {
                std::ifstream rows(pending/"products.csv");
                auto systems = zhangMultiGnssBuildSystemManifest(rows, certificate, systemSignals);
                manifest["schema_version"] = 2;
                manifest["systems"] = std::move(systems);
            }
            manifest["products_bytes"]=static_cast<std::int64_t>(std::filesystem::file_size(pending/"products.csv"));
            manifest["covariance_bytes"]=static_cast<std::int64_t>(std::filesystem::exists(pending/"covariance.csv")?std::filesystem::file_size(pending/"covariance.csv"):0);
            std::ofstream out(pending/"COMMIT.json");
            out<<boost::json::serialize(manifest)<<'\n';out.flush();
            if(!out)throw std::runtime_error("R51_MANIFEST_WRITE_FAILED");}
        {std::ofstream out(pending/"SHA256SUMS");
            for(const auto& name:{"products.csv","covariance.csv","certificate.txt","COMMIT.json"})
                if(std::filesystem::exists(pending/name))out<<fileHash(pending/name)<<"  "<<name<<"\n";
            out.flush();if(!out)throw std::runtime_error("R51_HASH_MANIFEST_WRITE_FAILED");}
        const auto verified=zhangR51VerifyBundleContents(pending,ar,true);
        if(!verified.valid || verified.epoch!=epoch || verified.route!=route ||
           verified.root!=root || verified.ar!=ar)
            throw std::runtime_error("R51_BUNDLE_PREFLIGHT_FAILED:"+verified.failureReason);
        std::filesystem::rename(pending,committed);
    }
};
