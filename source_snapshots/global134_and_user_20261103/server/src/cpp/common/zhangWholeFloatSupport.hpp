#pragma once
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include <array>
#include <functional>

struct ZhangWholeFloatSupportResult
{
    bool valid=false,oldAux=false;
    std::string reason;
    std::vector<std::map<KFKey,double>> rows;
    ZhangBridgeSpanRow offsets;
    std::array<bool,2> current{},signals{};
};
using ZhangWholeFloatArcResolver = std::function<bool(
    const std::string&,std::map<KFKey,double>&,std::string&)>;
using ZhangWholeFloatCurrentVerifier = std::function<bool(std::string&)>;
struct ZhangPreparedWholeFloatSupport;

// Non-inline production assembly shared with the live event-overlap path.
// H/b belong to an already verified feasible CURRENT integer domain. The arc
// callback supplies actual Phi/lambda rows and observation-status diagnostics;
// these inputs do not themselves grant integer/publication authority. Typed C
// physicalScale/nuisance bindings are read from the actual KF metadata. This
// resolves numerical search rows without changing the integer image/frame.
ZhangWholeFloatSupportResult zhangResolveWholeFloatSupportProduction(
    const KFState& state,E_Sys system,const std::array<E_ObsCode,2>& codes,
    const std::vector<std::string>& physicalIds,
    const ZhangExactMatrix& currentH,const ZhangExactVector& currentValues,
    const ZhangBridgeSpanMatrix& queries,const ZhangBridgeSpanRow& offsets,
    const ZhangWholeFloatArcResolver& resolveArc,
    const ZhangWholeFloatCurrentVerifier& verifyCurrentBindings,
    std::ostream* trace=nullptr,const std::string& targetLabel={},
    ZhangPreparedWholeFloatSupport* prepared=nullptr);
