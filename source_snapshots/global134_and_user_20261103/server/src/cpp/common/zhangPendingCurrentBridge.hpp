#pragma once
#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangCurrentBridgeWorkIdentity.hpp"
#include "common/zhangCurrentProductBridgePlan.hpp"
#include "common/zhangSharedProductEvent.hpp"
#include "common/zhangR51PhysicalImage.hpp"
#include "ambres/GNSSambres.hpp"

struct ZhangPreparedWholeFloatSupport;

// Work receipt only. exactAudited means the original integer-image/lift checks
// finished on the frozen parent, NOT that any relation has been authorised.
// No field is exported in a checkpoint or consumed by the product writer.
struct ZhangPendingCurrentBridge
{
    std::weak_ptr<const ZhangCurrentProductFloatOwnerToken> owner;
    const KFState* stateOwner=nullptr;
    std::uint64_t sequence=0;
    long int epoch=0;
    E_Sys system=E_Sys::NONE;
    std::array<E_ObsCode,2> codes{};
    std::array<ZhangSharedProductTargetKey,2> targets;
    std::array<bool,2> mapProven{};
    std::array<std::string,2> mapWitnesses;
    ZhangCurrentBridgeConditioningKey conditioning;
    std::string modelAndGraphIdentity;
    std::shared_ptr<ZhangCurrentEpochEventTransaction> transaction;
    std::shared_ptr<ZhangPreparedWholeFloatSupport> support;
    std::shared_ptr<const GinAR_mtx> candidate;
    Eigen::VectorXd candidateInputMean;
    Eigen::MatrixXd candidateInputCovariance;
    ZhangR51PhysicalImage candidateImage;
    bool exactAudited=false,numericAudited=false,proofsTraced=false,prepared=false;
    std::size_t stagedRows=0;
    std::vector<ZhangCurrentEpochPhysicalRow> physicalRows;
    std::vector<ZhangCoordinateRow> numericRows;
    ZhangExactVector values,targetValues;
    ZhangDecisionProofs proofs;
    std::string witness;
    SatSys satellite,referenceSatellite;
    int residualDimension=0,changeMask=0;
    bool oldAuxSeen=false;
    double ratio=0,threshold=0;
};

// The live event path and linked native regression use the same finalisation
// boundary. The receipt remains private until the base transaction commits.
bool zhangFinishPendingCurrentBridge(std::ostream&,const KFState&,
    ZhangCurrentEpochEvidence&,ZhangPendingCurrentBridge&,std::string&);

