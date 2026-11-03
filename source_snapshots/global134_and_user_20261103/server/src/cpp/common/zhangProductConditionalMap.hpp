#pragma once

#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangUserIntegerFunctional.hpp"

/** A coordinate identity is not an integer decision.  This result only says
 * that the actual new product row differs from the transported old product
 * row by the stated versioned physical integer functional.  The writer must
 * independently certify that functional from current-epoch WL/L1 evidence. */
struct ZhangProductConditionalMapResult
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
    ZhangCoordinateRow bridgeStateRowCycles;
    std::vector<ZhangGraphEdge> physicalEdges;
    std::vector<int> physicalArcVersions;
    std::vector<int> physicalCoefficients;
};

// A complete versioned FLOAT functional may be recoverable from typed
// carriers even though none of its individual old arcs can be reconstructed.
// The caller supplies a row already resolved with an exact whole-span witness.
// This proves the numeric map only; it never certifies the integer RHS.
inline ZhangProductConditionalMapResult zhangProveWholeProductConditionalMap(
    const ZhangCoordinateRow& oldRow,const ZhangCoordinateRow& writerRow,
    const ZhangProductIntegerTransition& difference,double wavelength,
    bool oldJointProven,bool completeSpanProven,
    const ZhangCoordinateRow& completeDifferenceMetres,double offsetMetres)
{
    ZhangProductConditionalMapResult result;
    if(!oldJointProven || !completeSpanProven)
    {result.reason="WHOLE_PRODUCT_JOINT_MOMENT_OR_EXACT_SPAN_UNPROVEN";return result;}
    if(oldRow.empty() || writerRow.empty() || !difference.valid ||
        difference.physicalEdges.size()!=difference.physicalArcVersions.size() ||
        difference.physicalEdges.size()!=difference.coefficients.size() ||
        difference.affineOffsetCycles!=0 || offsetMetres!=0 ||
        !(wavelength>0) || !std::isfinite(wavelength))
    {result.reason="INVALID_WHOLE_PRODUCT_FUNCTIONAL_DIFFERENCE";return result;}
    for(size_t i=0;i<difference.coefficients.size();++i)
    {
        const auto& value=difference.coefficients[i];
        const int coefficient=value.convert_to<int>();
        if(difference.physicalArcVersions[i]<0 || ZhangExactInteger(coefficient)!=value)
        {result.reason="PRODUCT_BRIDGE_VERSION_OR_COEFFICIENT_OUT_OF_RANGE";return result;}
        result.physicalEdges.push_back(difference.physicalEdges[i]);
        result.physicalArcVersions.push_back(difference.physicalArcVersions[i]);
        result.physicalCoefficients.push_back(coefficient);
    }
    for(const auto& [key,value]:completeDifferenceMetres)
    {
        if(!std::isfinite(value))
        {result.reason="WHOLE_PRODUCT_BRIDGE_ROW_NONFINITE";return result;}
        if(value!=0)result.bridgeStateRowCycles.emplace(key,value/wavelength);
    }
    ZhangCoordinateRow residual=writerRow;
    zhangAddCoordinateRow(residual,oldRow,-1);
    zhangAddCoordinateRow(residual,completeDifferenceMetres,-1);
    for(const auto& [key,value]:residual)
        if(!std::isfinite(value) || std::abs(value)>1e-10)
        {result.reason="WHOLE_PRODUCT_FUNCTIONAL_BRIDGE_ROW_IDENTITY_FAILED";return result;}
    result.valid=true;result.reason="EXACT_WHOLE_VERSIONED_PRODUCT_FUNCTIONAL_MAP_NO_INTEGER_AUTHORITY";
    return result;
}

inline ZhangProductConditionalMapResult zhangProveProductConditionalMap(
    const ZhangCoordinateRow& transportedOldProductRow,
    const ZhangCoordinateRow& actualNewWriterRow,
    const ZhangProductIntegerTransition& versionedDifference,
    double wavelength,
    bool oldJointMomentTransportProven,
    const std::function<bool(const ZhangGraphEdge&, int,
        ZhangCoordinateRow&, std::string&)>& versionedPhysicalRow)
{
    ZhangProductConditionalMapResult result;
    if (!oldJointMomentTransportProven)
    {
        result.reason = "OLD_PRODUCT_JOINT_MOMENT_TRANSPORT_UNPROVEN";
        return result;
    }
    if (transportedOldProductRow.empty() || actualNewWriterRow.empty() ||
        !versionedDifference.valid ||
        versionedDifference.physicalEdges.size() !=
            versionedDifference.physicalArcVersions.size() ||
        versionedDifference.physicalEdges.size() !=
            versionedDifference.coefficients.size() ||
        versionedDifference.affineOffsetCycles != 0 ||
        !(wavelength > 0) || !std::isfinite(wavelength))
    {
        result.reason = "INVALID_PRODUCT_FUNCTIONAL_DIFFERENCE";
        return result;
    }
    for (size_t term = 0; term < versionedDifference.coefficients.size();
         term++)
    {
        const auto& edge = versionedDifference.physicalEdges[term];
        const int version = versionedDifference.physicalArcVersions[term];
        ZhangCoordinateRow physical;
        std::string rowReason;
        if (version < 0 ||
            !versionedPhysicalRow(edge, version, physical, rowReason) ||
            physical.empty())
        {
            result.reason = rowReason.empty()
                ? "VERSIONED_PRODUCT_PHYSICAL_ROW_MISSING" : rowReason;
            return result;
        }
        int exactCoefficient;
        double numericCoefficient;
        try
        {
            exactCoefficient =
                versionedDifference.coefficients[term].convert_to<int>();
            numericCoefficient =
                versionedDifference.coefficients[term].convert_to<double>();
        }
        catch (...)
        {
            result.reason = "PRODUCT_BRIDGE_COEFFICIENT_OUT_OF_RANGE";
            return result;
        }
        if (ZhangExactInteger(exactCoefficient) !=
                versionedDifference.coefficients[term] ||
            !std::isfinite(numericCoefficient))
        {
            result.reason = "PRODUCT_BRIDGE_COEFFICIENT_OUT_OF_RANGE";
            return result;
        }
        for (const auto& [key, coefficient] : physical)
            if (!std::isfinite(coefficient))
            {
                result.reason = "PRODUCT_BRIDGE_ROW_NONFINITE";
                return result;
            }
        zhangAddCoordinateRow(result.bridgeStateRowCycles, physical,
            numericCoefficient / wavelength);
        result.physicalEdges.push_back(edge);
        result.physicalArcVersions.push_back(version);
        result.physicalCoefficients.push_back(exactCoefficient);
    }
    ZhangCoordinateRow residual = actualNewWriterRow;
    zhangAddCoordinateRow(residual, transportedOldProductRow, -1);
    zhangAddCoordinateRow(residual, result.bridgeStateRowCycles, -wavelength);
    for (const auto& [key, coefficient] : residual)
        if (!std::isfinite(coefficient) || std::abs(coefficient) > 1e-10)
        {
            result.reason = "PRODUCT_FUNCTIONAL_BRIDGE_ROW_IDENTITY_FAILED";
            return result;
        }
    result.valid = true;
    result.reason = versionedDifference.coefficients.empty()
        ? "STRUCTURAL_ZERO_PRODUCT_FUNCTIONAL_MAP"
        : "EXACT_VERSIONED_PRODUCT_FUNCTIONAL_MAP";
    return result;
}

struct ZhangProductOverlapNumericRowCheck
{
    bool valid = false;
    std::string reason = "NOT_EVALUATED";
};

/** The exact versioned-edge witness and the posterior row used by LAMBDA
 * must describe the same bridge.  The former is checked separately as an
 * integer identity; this check compares every KF coordinate coefficient in
 * the numerical row, including retained old-arc auxiliaries.  A two-signal
 * event bridge also needs a currently admitted arc in each signal. */
inline ZhangProductOverlapNumericRowCheck zhangCheckProductOverlapNumericRows(
    const std::array<ZhangCoordinateRow, 2>& subjectBridgeCycles,
    const std::array<ZhangCoordinateRow, 2>& referenceBridgeCycles,
    const std::array<ZhangCoordinateRow, 2>& searchedRowsCycles,
    const std::array<bool, 2>& currentAdmittedArcSeen,
    bool correlatedOldAuxiliarySeen,
    bool retiredArcActuallyRequired = true)
{
    ZhangProductOverlapNumericRowCheck result;
    if (retiredArcActuallyRequired && !correlatedOldAuxiliarySeen)
    {
        result.reason = "RETIRED_ARC_AUXILIARY_REQUIRED";
        return result;
    }
    for (int signal = 0; signal < 2; signal++)
    {
        ZhangCoordinateRow expected = subjectBridgeCycles[signal];
        zhangAddCoordinateRow(expected, referenceBridgeCycles[signal], -1);
        bool nonzero = false;
        for (const auto& [key, value] : expected)
        {
            if (!std::isfinite(value))
            {
                result.reason = "PRODUCT_BRIDGE_NUMERIC_ROW_NONFINITE";
                return result;
            }
            nonzero = nonzero || std::abs(value) > 1e-10;
        }
        // A structural zero has no integer decision and needs no replacement
        // arc. Its actual search row must nevertheless be identically zero.
        if (nonzero && !currentAdmittedArcSeen[signal])
        {
            result.reason = signal == 0
                ? "L1_CURRENT_ADMITTED_ARC_REQUIRED"
                : "L2_CURRENT_ADMITTED_ARC_REQUIRED";
            return result;
        }
        if (nonzero && searchedRowsCycles[signal].empty())
        {
            result.reason = signal == 0
                ? "L1_PRODUCT_BRIDGE_NUMERIC_ROW_EMPTY"
                : "L2_PRODUCT_BRIDGE_NUMERIC_ROW_EMPTY";
            return result;
        }
        ZhangCoordinateRow residual = searchedRowsCycles[signal];
        zhangAddCoordinateRow(residual, expected, -1);
        for (const auto& [key, coefficient] : residual)
            if (!std::isfinite(coefficient) ||
                std::abs(coefficient) > 1e-10)
            {
                result.reason = signal == 0
                    ? "L1_PRODUCT_BRIDGE_NUMERIC_ROW_MISMATCH"
                    : "L2_PRODUCT_BRIDGE_NUMERIC_ROW_MISMATCH";
                return result;
            }
    }
    result.valid = true;
    result.reason = "DUAL_SIGNAL_VERSIONED_PRODUCT_BRIDGE_ROWS_MATCH";
    return result;
}
