#pragma once

#include <array>
#include <memory>
#include <string>

#include "common/algebra.hpp"
#include "common/zhangGraphCoordinateTransport.hpp"
#include "common/zhangUserIntegerFunctional.hpp"

/** Issued only after the production factory has hashed the complete immutable
 * FLOAT x/P and checked its current evidence binding. The owner is borrowed:
 * use within a single-threaded preparation context, never mutate it or retain
 * this token across a filter commit. This is not an integer authority. */
class ZhangCurrentProductFloatOwnerToken
{
    const KFState* owner=nullptr;
    std::uint64_t commitSequence=0;
    long int epoch=0;
    Eigen::Index dimension=0;
    std::string runtime,posterior,chart,root,order;
    const KFState* conditionedOwner=nullptr;
    std::uint64_t conditionedSequence=0;
    std::string conditionedRoot,conditionedOrder,conditionedReason="NO_CONDITIONED_BASE_OWNER";
    std::map<ZhangCoordinateRow,ZhangExactInteger> conditionedRecipe;
    ZhangCurrentProductFloatOwnerToken()=default;
    friend std::shared_ptr<const ZhangCurrentProductFloatOwnerToken>
        zhangBindCurrentProductFloatOwner(const KFState&,E_Sys,long int,
            const std::string&,const std::string&,std::string&,const KFState*);
public:
    const std::string& floatRoot() const {return root;}
    const std::string& stateOrder() const {return order;}
    const std::string& conditionedBaseReason() const {return conditionedReason;}
    const std::map<ZhangCoordinateRow,ZhangExactInteger>& appliedRecipe() const {return conditionedRecipe;}
    const KFState* conditionedBase(const std::string& actualOrder) const
    {
        if(!conditionedOwner || conditionedRoot.empty() || conditionedRecipe.empty() ||
            conditionedOwner->factorCommitSequence!=conditionedSequence ||
            static_cast<long int>(conditionedOwner->time.bigTime)!=epoch ||
            conditionedOwner->x.size()!=dimension || conditionedOwner->P.rows()!=dimension ||
            conditionedOwner->P.cols()!=dimension || actualOrder!=conditionedOrder ||
            conditionedOrder!=order)return nullptr;
        return conditionedOwner;
    }
    const KFState* conditionedBaseOwner() const {return conditionedOwner;}
    bool matches(const KFState& state,const std::string& runtimeId,
        long int epochId,const std::string& posteriorId,
        const std::string& chartId,const std::string& actualOrder) const
    {
        return owner==&state && commitSequence==state.factorCommitSequence &&
            epoch==epochId && epoch==static_cast<long int>(state.time.bigTime) &&
            dimension==state.x.size() && dimension==state.P.rows() &&
            dimension==state.P.cols() && runtime==runtimeId &&
            posterior==posteriorId && chart==chartId && order==actualOrder &&
            !root.empty() && !order.empty();
    }
};

// Hashing is performed once, outside individual bridge budgets. A caller
// cannot obtain a token merely by supplying a purported root string.
std::shared_ptr<const ZhangCurrentProductFloatOwnerToken>
zhangBindCurrentProductFloatOwner(const KFState& state,E_Sys system,
    long int epoch,const std::string& posteriorId,const std::string& chartId,
    std::string& reason,const KFState* conditionedBaseOwner=nullptr);

/** Complete current physical targets, not merely the latest graph increment.
 * fullJointMapProven authorises the geometric map only; it supplies no integer
 * value. Both slots must identify their signal. An unavailable other segment
 * may use a proven native-to-same-native structural zero, but must not pretend
 * to authorise that signal's unavailable old segment. */
struct ZhangCurrentProductBridgePlan
{
    std::shared_ptr<const ZhangCurrentProductFloatOwnerToken> immutableFloat;
    SatSys satellite;
    std::array<ZhangProductIntegerTransition,2> completeSignalTargets;
    std::array<bool,2> fullJointMapProven{};
    std::array<std::string,2> fullJointMapWitness;
    int changeMask=0;
};
