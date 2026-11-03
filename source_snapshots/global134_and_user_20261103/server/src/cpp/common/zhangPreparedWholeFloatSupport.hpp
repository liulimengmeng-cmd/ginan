#pragma once
#include "common/zhangWholeFloatSupport.hpp"
#include "common/zhangCurrentProductBridgePlan.hpp"
#include "common/zhangCurrentEpochEvidence.hpp"
#include <type_traits>

// Ephemeral, non-serialised work. Nothing in this object grants integer or
// publication authority. The owner token and exact domain must BOTH survive;
// every caller additionally checks the actual owner sequence/root/order.
struct ZhangPreparedWholeFloatSupport
{
    using SymbolRow=std::map<std::string,ZhangBridgeSpanRational>;
    struct BoundRow {SymbolRow symbols;ZhangCoordinateRow numeric;bool oldAux=false;};
    std::weak_ptr<const ZhangCurrentProductFloatOwnerToken> owner;
    std::weak_ptr<const ZhangCurrentEpochCompiledDomain> domain;
    const KFState* stateOwner=nullptr;
    std::uint64_t sequence=0,generation=0;
    long int epoch=0;
    E_Sys system=E_Sys::NONE;
    std::array<E_ObsCode,2> codes{};
    std::string runtime,posterior,chart,root,order,modelAndGraphIdentity;
    std::vector<std::string> physicalIds;
    ZhangExactMatrix h;
    ZhangExactVector rhs;
    std::map<std::string,std::string> metadata;
    bool snapshotStarted=false,snapshotComplete=false,decodeStarted=false,decodeComplete=false,
        metadataValid=false,catalogueComplete=false,prepared=false,baseComplete=false;
    std::string snapshotCursor,decodeCursor,reason="PREPARED_SUPPORT_NOT_READY";
    std::vector<ZhangProductFloatTargetRequest> registered;
    std::map<std::string,ZhangProductFloatTargetRequest> carriers;
    std::set<std::string> ids;
    std::size_t idCollectCursor=0,structuredCursor=0,hCursor=0,sCursor=0;
    std::string arcCursor;
    bool arcStarted=false,arcsComplete=false;
    std::map<std::string,BoundRow> availablePhi;
    std::vector<BoundRow> structuredPhi;
    std::map<std::string,ZhangCoordinateRow> arcRows;
    std::map<std::string,std::string> arcReasons;
    std::map<std::string,bool> arcResolved;
    std::vector<SymbolRow> fullH,reducedH,fullS,reducedS;
    std::vector<std::size_t> hIndices,sIndices;
    std::uint64_t decodeEntries=0,decodeBytes=0,arcCalls=0,arcHits=0,queryCalls=0,
        scopeChecks=0,orderBuilds=0;
    ~ZhangPreparedWholeFloatSupport()
    {
        ZhangPhaseTimer timer("BRIDGE_SUPPORT_CONTEXT_DESTRUCT");
        // Release all owning containers/string allocations before the timer
        // ends; no shrink/copy of live data is introduced by this measurement.
        auto release=[](auto& value) {std::decay_t<decltype(value)> empty;value.swap(empty);};
        release(registered);release(carriers);release(metadata);release(availablePhi);
        release(structuredPhi);release(arcRows);release(arcReasons);release(arcResolved);
        release(fullH);release(reducedH);release(fullS);release(reducedS);release(h);release(rhs);
        release(physicalIds);release(ids);release(hIndices);release(sIndices);
        release(runtime);release(posterior);release(chart);release(root);release(order);
        release(modelAndGraphIdentity);release(snapshotCursor);release(decodeCursor);release(arcCursor);release(reason);
    }
};

// Production factory also used by the linked native regression. A caller
// still needs the verified owner token; this factory grants no authority.
std::shared_ptr<ZhangPreparedWholeFloatSupport> zhangAcquirePreparedWholeFloatSupport(
    const KFState&,E_Sys,const std::array<E_ObsCode,2>&,
    const std::shared_ptr<const ZhangCurrentProductFloatOwnerToken>&,
    const std::shared_ptr<const ZhangCurrentEpochCompiledDomain>&,
    const std::vector<std::string>&,const std::string&,long int,
    const std::string&,const std::string&);

