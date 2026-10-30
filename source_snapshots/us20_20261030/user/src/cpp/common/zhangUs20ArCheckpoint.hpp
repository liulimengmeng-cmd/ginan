#pragma once
#include <string>
std::string exportZhangUs20UserArCheckpoint(const std::string& runtime,double completed);
void importZhangUs20UserArCheckpoint(const std::string& payload,const std::string& runtime,double completed,bool dryRun);
