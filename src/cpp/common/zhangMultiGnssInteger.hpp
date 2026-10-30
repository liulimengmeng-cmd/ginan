#pragma once

#include "common/algebra.hpp"
#include "common/zhangR47Candidate.hpp"
#include "common/zhangProductPhysicalCycleChart.hpp"

#include <array>
#include <map>
#include <string>

// Each constellation retains its own exact lattice, physical chart and
// certificate.  The delivery posterior below is nevertheless one full joint
// KF posterior: this map is never a request to block-diagonalise its P matrix.
using ZhangMultiGnssProductResults =
    std::map<E_Sys, ZhangProductRelationFixResult>;

struct ZhangMultiGnssProductSystemDelivery
{
    ZhangProductIntegerConstraintSet constraints;
    std::string sourcePosteriorId;
    std::string authoritativeCycleChartId;
};

struct ZhangMultiGnssProductDelivery
{
    bool valid = false;
    std::string runtimeId;
    long int epochGpstSeconds = 0;
    std::string floatRoot;
    std::string floatStateOrderId;
    std::string conditionedRoot;
    std::string conditionedStateOrderId;
    std::map<E_Sys, ZhangMultiGnssProductSystemDelivery> systems;
};

/** Positive current-epoch allow-list of the systems actually admitted to one
 * joint product posterior. A merely valid registry entry is not delivery
 * authority: a failed candidate may still leave diagnostic evidence there. */
struct ZhangMultiGnssCurrentEvidenceAuthority
{
    long int epochGpstSeconds=0;
    std::string floatRoot,floatStateOrderId;
    std::map<E_Sys,std::pair<std::string,std::string>> systemIdentities;
};

inline std::map<std::string,ZhangMultiGnssCurrentEvidenceAuthority>&
zhangMultiGnssCurrentEvidenceAuthorities()
{
    static std::map<std::string,ZhangMultiGnssCurrentEvidenceAuthority> registry;
    return registry;
}

inline bool zhangMultiGnssAuthorizeCurrentProductEvidence(
    const ZhangMultiGnssProductDelivery& delivery)
{
    auto& registry=zhangMultiGnssCurrentEvidenceAuthorities();
    registry.erase(delivery.runtimeId);
    if(!delivery.valid || delivery.runtimeId.empty() || delivery.epochGpstSeconds<=0 ||
        delivery.floatRoot.empty() || delivery.floatStateOrderId.empty() ||
        delivery.systems.empty())return false;
    ZhangMultiGnssCurrentEvidenceAuthority scope;
    scope.epochGpstSeconds=delivery.epochGpstSeconds;
    scope.floatRoot=delivery.floatRoot;scope.floatStateOrderId=delivery.floatStateOrderId;
    for(const auto& [system,entry]:delivery.systems)
    {
        if(system==E_Sys::NONE || entry.constraints.system!=system ||
            !entry.constraints.reliable || entry.sourcePosteriorId.empty() ||
            entry.authoritativeCycleChartId.empty())return false;
        scope.systemIdentities.emplace(system,std::make_pair(
            entry.sourcePosteriorId,entry.authoritativeCycleChartId));
    }
    registry.emplace(delivery.runtimeId,std::move(scope));
    return true;
}

/** A global state/ambiguity order may contain zero columns from the other
 * constellations.  A nonzero integer row may not cross this system boundary.
 * Reject an out-of-range or missing column instead of dropping its term. */
inline bool zhangMultiGnssIntegerRowsUseSystem(
    const ZhangExactMatrix& rows, const std::map<int, KFKey>& coordinates,
    E_Sys system)
{
    if (system == E_Sys::NONE ||
        !zhangExactRectangularMatrix(rows, coordinates.size())) return false;
    for (const auto& row : rows)
        for (std::size_t column = 0; column < row.size(); ++column)
        {
            if (row[column] == 0) continue;
            const auto found = coordinates.find(static_cast<int>(column));
            if (found == coordinates.end() || found->second.Sat.sys != system)
                return false;
        }
    return true;
}

/** Physical arc identity is signal|receiver|system-qualified-satellite|Vn.
 * Integer certificates do not authorise a relation to an arc of another
 * constellation, even when its PRN and/or observable number are the same. */
inline bool zhangMultiGnssPhysicalRowsUseSystem(
    const std::vector<std::map<std::string, ZhangExactInteger>>& rows,
    E_Sys system, const std::array<E_ObsCode, 2>& codes)
{
    if (system == E_Sys::NONE || codes[0] == E_ObsCode::NONE ||
        codes[1] == E_ObsCode::NONE || codes[0] == codes[1]) return false;
    for (const auto& row : rows)
        for (const auto& [identity, coefficient] : row)
        {
            if (coefficient == 0) continue;
            const auto a = identity.find('|');
            const auto b = a == std::string::npos ? a : identity.find('|', a + 1);
            const auto c = b == std::string::npos ? b : identity.find('|', b + 1);
            if (a == std::string::npos || b == std::string::npos ||
                c == std::string::npos || a == 0 || b == a + 1 ||
                c == b + 1 || c + 2 >= identity.size() ||
                identity[c + 1] != 'V' || identity.find('|', c + 1) != std::string::npos)
                return false;
            const auto code = identity.substr(0, a);
            if (code != enum_to_string(codes[0]) && code != enum_to_string(codes[1]))
                return false;
            const auto id = identity.substr(b + 1, c - b - 1);
            const SatSys satellite(id.c_str());
            if (satellite.sys != system || satellite.prn <= 0 || satellite.id() != id)
                return false;
            for (std::size_t i = c + 2; i < identity.size(); ++i)
                if (identity[i] < '0' || identity[i] > '9') return false;
        }
    return true;
}
