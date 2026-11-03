#pragma once
#include "common/gTime.hpp"
#include <filesystem>
#include <string>
struct Navigation;

namespace zhangRealtimeNavigation
{
bool enabled();
GTime selectionEpoch(GTime evaluationTime);
Navigation& orbitalContext(Navigation& live);
// Server: retain the incoming history, select one common orbit/clock family
// per satellite, then freeze BEFORE SPP, nominal orbit or PPP calculations.
void beginServer(const GTime& time);
void endServer();
void writeSnapshot(const std::filesystem::path& directory, const std::string& epoch);
// User: validates a complete private value before replacing only orbital
// inputs. It never imports network KF, receiver or ambiguity state.
void importSnapshot(const std::filesystem::path& path, const GTime& time);
void missingEpoch(const GTime& time);
const std::string& snapshotHash();
} // namespace zhangRealtimeNavigation
