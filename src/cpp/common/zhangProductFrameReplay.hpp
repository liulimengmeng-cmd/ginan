#pragma once
#include "common/zhangIntegerAudit.hpp"
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>

// Opt-in sparse input capture for strict old/new exact-backend replay. No
// observation data, covariance approximation or production decision cache.
inline void zhangCaptureProductFrame(const ZhangExactMatrix& targets,
    const ZhangExactMatrix& rows,const ZhangExactVector& rhs,int dimension,int work)
{
    const char* directory=std::getenv("ZHANG_R51_EXACT_REPLAY_DIRECTORY");
    if(!directory || !*directory)return;
    static std::atomic<unsigned long> serial{0};
    const auto path=std::filesystem::path(directory)/("frame_"+std::to_string(++serial)+".txt");
    std::filesystem::create_directories(directory);
    if(std::filesystem::exists(path))throw std::runtime_error("PRODUCT_FRAME_REPLAY_EXISTS");
    std::ofstream out(path);
    out<<"R51_PRODUCT_FRAME_V1 "<<dimension<<" "<<work<<"\n";
    auto matrix=[&](const ZhangExactMatrix& m){
        out<<m.size()<<"\n";
        for(const auto& row:m){
            std::size_t count=0;for(const auto& v:row)count+=v!=0;
            out<<count;for(int c=0;c<row.size();++c)if(row[c]!=0)out<<" "<<c<<" "<<row[c];out<<"\n";
        }
    };
    matrix(targets);matrix(rows);out<<rhs.size();for(const auto& v:rhs)out<<" "<<v;out<<"\n";
    out.flush();if(!out)throw std::runtime_error("PRODUCT_FRAME_REPLAY_WRITE_FAILED");
}
