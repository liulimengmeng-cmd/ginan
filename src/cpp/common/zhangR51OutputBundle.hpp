#pragma once
#include "common/acsConfig.hpp"
#include "common/zhangR51Integrity.hpp"
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include "common/zhangCheckpoint.hpp"

/** Epoch directories are the publication authority. Aggregate CSVs are
 * derived convenience views and readers must use the committed manifest for
 * cross-file atomicity. Rejected private directories are retained for audit. */
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
    static void appendView(const std::filesystem::path& from,const std::string& to) {
        if(to.empty() || !std::filesystem::exists(from))return;
        const bool header=std::filesystem::exists(to) && std::filesystem::file_size(to)>0;
        std::ifstream in(from);std::ofstream out(to,std::ios::app);
        std::string line;if(header)std::getline(in,line);
        out<<in.rdbuf();out.flush();
        if(!out)throw std::runtime_error("R51_AGGREGATE_VIEW_WRITE_FAILED");
    }
    static std::string fileHash(const std::filesystem::path& path) {
        if(!std::filesystem::exists(path))return "";
        std::ifstream in(path,std::ios::binary);std::ostringstream bytes;bytes<<in.rdbuf();
        if(in.bad())throw std::runtime_error("R51_BUNDLE_READ_FAILED");return zhangCheckpointSha256(bytes.str());
    }
    void publish(const std::string& certificate,bool ar) {
        restore();
        if(!std::filesystem::exists(pending/"products.csv"))throw std::runtime_error("R51_PRODUCT_FILE_MISSING");
        if(ar && !std::filesystem::exists(pending/"covariance.csv"))throw std::runtime_error("R51_COVARIANCE_FILE_MISSING");
        {std::ofstream out(pending/"certificate.txt");out<<certificate;out.flush();if(!out)throw std::runtime_error("R51_CERTIFICATE_WRITE_FAILED");}
        {std::ofstream out(pending/"COMMIT.json");out<<"{\"epoch\":\""<<epoch<<"\",\"route\":\""<<route<<"\",\"root\":\""<<root<<"\",\"ar\":"<<(ar?"true":"false")<<",\"products_bytes\":"<<std::filesystem::file_size(pending/"products.csv")<<",\"covariance_bytes\":"<<(std::filesystem::exists(pending/"covariance.csv")?std::filesystem::file_size(pending/"covariance.csv"):0)<<"}\n";out.flush();if(!out)throw std::runtime_error("R51_MANIFEST_WRITE_FAILED");}
        {std::ofstream out(pending/"SHA256SUMS");
            for(const auto& name:{"products.csv","covariance.csv","certificate.txt","COMMIT.json"})
                if(std::filesystem::exists(pending/name))out<<fileHash(pending/name)<<"  "<<name<<"\n";
            out.flush();if(!out)throw std::runtime_error("R51_HASH_MANIFEST_WRITE_FAILED");}
        std::filesystem::rename(pending,committed);
        appendView(committed/"products.csv",product);appendView(committed/"covariance.csv",covariance);
    }
};
