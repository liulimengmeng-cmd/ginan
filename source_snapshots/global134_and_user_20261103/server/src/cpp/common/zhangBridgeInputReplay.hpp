#pragma once
#include "common/zhangIntegerAudit.hpp"
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>

/** Opt-in event input capture BEFORE expensive work. No observations or
 * previous integer decisions are imported by this diagnostic. Each new run
 * must use its own writable capture directory. */
inline void zhangCaptureBridgeSupportInput(const std::string& epoch,
    const std::string& satellite, std::uint64_t generation, const std::string& root,
    const std::vector<std::string>& ids, const ZhangExactMatrix& targets,
    const ZhangExactMatrix& rows, const ZhangExactVector& rhs,
    const std::set<std::size_t>& unavailable)
{
    const char* directory=std::getenv("ZHANG_R51_BRIDGE_REPLAY_DIRECTORY");
    if(!directory || !*directory)return;
    if(const char* filter=std::getenv("ZHANG_R51_BRIDGE_CAPTURE_EPOCHS");filter && *filter)
        if(epoch.size()<8 || std::string(filter).find(epoch.substr(epoch.size()-8))==std::string::npos)return;
    static std::atomic<unsigned long> serial{0};
    const auto number=++serial;
    if(number>2048)throw std::runtime_error("BRIDGE_REPLAY_CAPTURE_LIMIT");
    const auto path=std::filesystem::path(directory)/("bridge_"+std::to_string(number)+".txt");
    std::filesystem::create_directories(directory);
    if(std::filesystem::exists(path))throw std::runtime_error("BRIDGE_REPLAY_EXISTS");
    std::ofstream out(path);
    out<<"R51_BRIDGE_SUPPORT_V1 "<<ids.size()<<" "<<generation<<"\n"
       <<std::quoted(epoch)<<" "<<std::quoted(satellite)<<" "<<std::quoted(root)<<"\n";
    for(const auto& id:ids)out<<std::quoted(id)<<"\n";
    out<<unavailable.size();for(const auto c:unavailable)out<<" "<<c;out<<"\n";
    auto matrix=[&](const ZhangExactMatrix& input){
        out<<input.size()<<"\n";
        for(const auto& row:input){
            std::size_t count=0;for(const auto& value:row)count+=value!=0;
            out<<count;for(std::size_t c=0;c<row.size();++c)if(row[c]!=0)out<<" "<<c<<" "<<row[c];out<<"\n";
        }
    };
    matrix(targets);matrix(rows);out<<rhs.size();for(const auto& value:rhs)out<<" "<<value;out<<"\n";
    out.flush();if(!out)throw std::runtime_error("BRIDGE_REPLAY_CAPTURE_WRITE_FAILED");
}
