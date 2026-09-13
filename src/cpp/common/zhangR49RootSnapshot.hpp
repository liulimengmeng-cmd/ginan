#pragma once
#include <filesystem>
#include <vector>
#include <cstdint>
#include <atomic>
#include <fstream>
#include <cstdlib>
#include <stdexcept>
#include "common/eigenIncluder.hpp"
// Little-endian native IEEE754 doubles on the frozen x86_64 execution platform.
// One file per owned network posterior, not one copy for every bridge.
inline std::string zhangR49WriteRootSnapshot(const VectorXd& mean,const MatrixXd& covariance,
 const std::vector<std::string>& order,std::uint64_t serial) {
 const char* directory=std::getenv("ZHANG_R49_SNAPSHOT_DIRECTORY");if(!directory)return "DISABLED";
 std::filesystem::create_directories(directory);
 const auto base=std::filesystem::path(directory)/("root_"+std::to_string(serial));
 const auto file=base.string()+".bin",names=base.string()+".columns.txt";
 if(std::filesystem::exists(file)||std::filesystem::exists(names))throw std::runtime_error("R49_SNAPSHOT_EXISTS");
 if(covariance.rows()!=mean.size()||covariance.cols()!=mean.size()||order.size()!=mean.size())throw std::runtime_error("R49_SNAPSHOT_SHAPE");
 std::ofstream out(file,std::ios::binary);
 out<<"R49_POSTERIOR_V1_COLUMN_MAJOR\n"<<mean.size()<<"\n";
 out.write(reinterpret_cast<const char*>(mean.data()),mean.size()*sizeof(double));
 out.write(reinterpret_cast<const char*>(covariance.data()),covariance.size()*sizeof(double));
 out.flush();if(!out)throw std::runtime_error("R49_SNAPSHOT_WRITE_FAILED");
 std::ofstream columns(names);for(int i=0;i<order.size();++i)columns<<i<<" "<<order[i]<<"\n";
 columns.flush();if(!columns)throw std::runtime_error("R49_SNAPSHOT_COLUMNS_FAILED");return file;
}
