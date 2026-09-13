#pragma once
#include <filesystem>
#include <fstream>
#include <atomic>
#include <cstdlib>
#include "common/eigenIncluder.hpp"
inline void zhangDumpUserMeasurement(const VectorXd& x,const MatrixXd& p,
 const MatrixXd& h,const VectorXd& v,const MatrixXd& r,
 const VectorXd& xp,const MatrixXd& pp,const std::string& label) {
 const char* directory=std::getenv("ZHANG_USER_MEASUREMENT_DUMP");if(!directory)return;
 static std::atomic<unsigned long> serial{0};
 std::filesystem::create_directories(directory);
 const auto path=std::filesystem::path(directory)/("measurement_"+std::to_string(++serial)+".bin");
 if(std::filesystem::exists(path))throw std::runtime_error("USER_MEASUREMENT_DUMP_EXISTS");
 std::ofstream out(path,std::ios::binary);
 out<<"USER_MEASUREMENT_V1\n"<<label<<"\n";
 auto write=[&](const char* name,const MatrixXd& a){out<<name<<" "<<a.rows()<<" "<<a.cols()<<"\n";out.write(reinterpret_cast<const char*>(a.data()),a.size()*sizeof(double));out<<"\n";};
 write("x",x);write("P",p);write("H",h);write("v",v);write("R",r);write("xp",xp);write("Pp",pp);
 out.flush();if(!out)throw std::runtime_error("USER_MEASUREMENT_DUMP_FAILED");
}
