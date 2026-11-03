#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include "common/zhangIntegerAudit.hpp"

using ZhangR51PhysicalRow = std::map<std::string, ZhangExactInteger>;
using ZhangR51PhysicalRows = std::vector<ZhangR51PhysicalRow>;

struct ZhangR51PhysicalConsequenceMetrics
{
    std::size_t rows = 0, columns = 0, targets = 0;
    std::size_t builds = 0, directQueries = 0;
    std::size_t inputNonzeros = 0, unitPivots = 0, fillAdded = 0;
    std::size_t residualRows = 0, residualColumns = 0, residualNonzeros = 0;
    std::size_t basisRows = 0, basisNonzeros = 0;
    double buildSeconds = 0, sparseSeconds = 0, hnfSeconds = 0;
    double feasibilitySeconds = 0, querySeconds = 0;
};

// Production entry point, implemented in ppp_ambres.cpp and shared by
// PhysicalStore and the linked reference tests. A true result means every
// integer solution of this *feasible* domain satisfies all targets. It does
// not authorize an integer row-combination certificate or product publication.
// Unknown physical identities are free variables, never implicit zeros.
bool zhangR51PhysicalConsequence(
    const ZhangR51PhysicalRows& rows,
    const ZhangExactVector& values,
    bool domainValid,
    const ZhangR51PhysicalRows& targets,
    const ZhangExactVector& rhs,
    ZhangR51PhysicalConsequenceMetrics* metrics = nullptr);
