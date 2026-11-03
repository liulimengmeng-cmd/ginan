// Real updateZhangFullRankReferences/applyStateTransform, with an explicitly
// constructed CURRENT evidence fixture on a literal positive-semidefinite
// FLOAT posterior. This covers the production commit-boundary API contract;
// it is NOT a real atomic product writer or real-data AR acceptance test.
#include "common/acsConfig.hpp"
#include "common/common.hpp"
#include "common/constants.hpp"
#include "common/receiver.hpp"
#include "common/observations.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangCurrentEpochEvidence.hpp"
#include "common/zhangCertifiedTreeTransport.hpp"
#include "common/zhangR51Integrity.hpp"
#include "ambres/GNSSambres.hpp"
#include "pea/zhangReference.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

static void require(bool ok,const std::string& why) { if (!ok) throw std::runtime_error(why); }
static const std::vector<E_ObsCode> codes={E_ObsCode::L1C,E_ObsCode::L2W};
static KFKey clockKey(SatSys satellite) { KFKey key; key.type=KF::SAT_CLOCK; key.Sat=satellite; return key; }
static void addState(KFState& state,const KFKey& key,double mean)
{ InitialState initial; initial.estimate=true; initial.x=mean; initial.P=1; initial.Q=0; state.addKFState(key,initial); }

struct Fixture
{
    KFState state;
    ReceiverMap receivers;
    std::ostringstream trace;
    ZhangGraphIntegerContext before;
    explicit Fixture(bool reverse,const std::string& owner,bool sharedProductBranch=false)
    {
        // CURRENT provenance requires a real positive epoch, not noTime().
        state.time.bigTime=1405216800;
        state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]=owner;
        std::vector<ZhangGraphEdge> edges={
            {"A",SatSys("G01")},{"A",SatSys("G02")},{"B",SatSys("G01")},
            {"B",SatSys("G02")},{"C",SatSys("G01")},{"C",SatSys("G02")}};
        if(sharedProductBranch)
        {
            edges.push_back({"B",SatSys("G03")});
            edges.push_back({"B",SatSys("G04")});
        }
        if (reverse) std::reverse(edges.begin(),edges.end());
        for(const auto& edge:edges) addObservation(edge);
        updateZhangFullRankReferences(trace,receivers,state);
        require(zhangGraphIntegerContext(state,E_Sys::GPS,before),"production initial graph missing");
        auto expectedTree=std::set<ZhangGraphEdge>{edgesFor("A",1),edgesFor("A",2),edgesFor("B",1),edgesFor("C",1)};
        if(sharedProductBranch) { expectedTree.insert(edgesFor("B",3));expectedTree.insert(edgesFor("B",4)); }
        require(before.basis.treeEdges==expectedTree,
            "fixture initial tree differs from the literal oracle tree");
        for(const auto code:codes)
        {
            for(const auto& receiver:before.basis.receivers)
                if(receiver!="A") addState(state,zhangTransportReceiverKey(E_Sys::GPS,code,receiver),20);
            for(const auto& satellite:before.basis.satellites)
                addState(state,zhangTransportSatelliteKey(code,satellite),10+satellite.prn);
            addState(state,zhangTransportAmbiguityKey(code,edgesFor("B",2)),11.8);
            addState(state,zhangTransportAmbiguityKey(code,edgesFor("C",2)),18.4);
        }
        for(const auto& satellite:before.basis.satellites) addState(state,clockKey(satellite),100+satellite.prn);
        require(state.stateTransition(trace,state.time),"fixture state creation failed");
        VectorXd v=VectorXd::LinSpaced(state.x.size(),0.01,0.25);
        state.P=MatrixXd::Identity(state.x.size(),state.x.size())+v*v.transpose();
        const int one=state.kfIndexMap.at(KFState::oneKey);
        state.x(one)=1; state.P.row(one).setZero(); state.P.col(one).setZero();
        state.lastPppTransactionCommitted=true; state.lastPppQcConverged=true;
    }
    static ZhangGraphEdge edgesFor(const char* receiver,int satellite)
    { return {receiver,SatSys(E_Sys::GPS,satellite)}; }
    void addObservation(const ZhangGraphEdge& edge)
    {
        auto& receiver=receivers[edge.receiver]; receiver.id=edge.receiver; receiver.ready=true;
        auto observation=std::make_shared<GObs>(); observation->Sat=edge.satellite;
        auto& status=receiver.satStatMap[edge.satellite];
        // C-G02 is deliberately much higher quality than certified B-G02.
        status.el=edge.receiver=="B"?0.5:1.4; observation->satStat_ptr=&status;
        for(const auto code:codes)
        {
            auto& signal=observation->sigs[code2Freq[E_Sys::GPS][code]];
            signal.code=code; signal.L=1e8; signal.P=2e7;
            auto& options=acsConfig.getRecOpts(edge.receiver,{"GPS",enum_to_string(code)});
            options.ambiguity.sigma={1000}; options.phase_bias.sigma={1000};
        }
        receiver.obsList.push_back(observation);
    }
    void installCurrentEvidenceAndCommitContract()
    {
        const std::string owner=state.metaDataMap.at(ZHANG_CHECKPOINT_RUNTIME_ID_METADATA);
        ZhangProductPhysicalCycleChart chart; chart.columns=4;
        int column=0;
        for(const auto code:codes)
        {
            std::map<ZhangGraphEdge,int> versions;
            for(const auto& [edge,signals]:before.signalArcVersions) versions[edge]=signals.at(static_cast<int>(code));
            require(chart.add(column++,enum_to_string(code),edgesFor("B",2),before.basis,versions),"B cycle chart fixture failed");
            require(chart.add(column++,enum_to_string(code),edgesFor("C",2),before.basis,versions),"C cycle chart fixture failed");
        }
        const auto epoch=static_cast<long int>(state.time.bigTime);
        const auto root=zhangR51NumericRoot(state.x,state.P);
        const auto order=zhangCurrentEvidenceStateOrderId(state.kfIndexMap);
        auto& evidence=zhangCurrentEpochEvidenceRegistry()[{owner,E_Sys::GPS}];
        require(evidence.bind(owner,E_Sys::GPS,epoch,"FIXTURE_POSTERIOR","FIXTURE_CHART",chart,root),"current bind failed");
        // PSD is established independently by the literal I+vv^T constructor
        // on stochastic coordinates and deterministic zero-variance ONE.
        evidence.floatPosteriorPsdValidated=true; evidence.floatPosteriorPsdRootId=root;
        evidence.floatPosteriorPsdStateOrderId=order;
        column=0;
        for(const auto code:codes)
        {
            ZhangExactVector coordinate(4); coordinate[column]=1;
            ZhangCurrentEpochPhysicalRow physical;
            require(chart.expand(coordinate,physical),"physical expansion failed");
            auto proof=std::make_shared<ZhangIntegerDecisionProof>();
            proof->id=owner+"|FIXTURE_ACCEPTED|"+enum_to_string(code);
            proof->originalStatement="AUDITABLE_FIXTURE_CURRENT_B_G02_CYCLE_EQUALS_12_NOT_REAL_AR_PUBLICATION";
            proof->observationProvenance=owner+"|PPP_COMMITTED|E"+std::to_string(epoch);
            proof->conditionalFailureBound=0;
            require(evidence.add(physical,12,{proof}) && evidence.addNumericConstraint(
                {{zhangTransportAmbiguityKey(code,edgesFor("B",2)),1}},12,physical),"current accepted/bound row fixture failed");
            column+=2;
        }
        const auto originalMetadata=state.metaDataMap;
        const VectorXd x=state.x; const MatrixXd p=state.P;
        const std::map<E_Sys,std::pair<std::string,std::string>> scopes{
            {E_Sys::GPS,{"FIXTURE_POSTERIOR","FIXTURE_CHART"}}};
        ZhangPreparedCertifiedTreeTransportHistory prepared;
        std::string reason;
        auto privateCopy=state;
        privateCopy.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_BRANCH_ID_METADATA]=owner+"/branch/private";
        require(!zhangPrepareCertifiedTreeTransportHistory(trace,privateCopy,scopes,prepared,&reason),"private branch was promoted to base history");
        require(zhangPrepareCertifiedTreeTransportHistory(trace,state,scopes,prepared,&reason),"native current source prepare rejected: "+reason);
        require(state.metaDataMap==originalMetadata && (state.x-x).isZero() && (state.P-p).isZero(),"read-only staging changed FLOAT or metadata");
        auto wrongEpoch=prepared; ++wrongEpoch.sourceEpoch;
        require(!zhangCommitCertifiedTreeTransportHistory(trace,state,wrongEpoch,&reason),"staged source crossed epoch");
        // Boundary contract invocation only: this unit does not create a real
        // product bundle. End-to-end atomic publication provenance is a
        // separate native replay gate and must not be inferred from this test.
        require(zhangCommitCertifiedTreeTransportHistory(trace,state,prepared,&reason),"same-owner commit contract rejected: "+reason);
        require((state.x-x).isZero() && (state.P-p).isZero() && evidence.rows.size()==2,
            "transport commit conditioned FLOAT or altered CURRENT accepted rows");
    }
    ~Fixture() { eraseZhangGraphRuntime(state); }
};

static void runPostQcRerouteCase(bool history,bool cut,bool rejectCallback)
{
    Fixture fixture(false,"POST_QC_REROUTE_"+std::to_string(history)+std::to_string(cut)+std::to_string(rejectCallback),true);
    if(history) fixture.installCurrentEvidenceAndCommitContract();
    fixture.state.time=fixture.state.time+30;
    updateZhangFullRankReferences(fixture.trace,fixture.receivers,fixture.state);
    fixture.state.lastPppTransactionCommitted=true;fixture.state.lastPppQcConverged=true;
    const auto failed=Fixture::edgesFor("B",1);
    KFMeas measurement;measurement.time=fixture.state.time;
    unsigned int rejects=0;
    for(const auto& edge:fixture.before.basis.edges) for(const auto code:codes)
    {
        KFKey key;key.type=KF::PHAS_MEAS;key.str=edge.receiver;key.Sat=edge.satellite;key.num=int(code);
        measurement.obsKeys.push_back(key);
        std::map<std::string,void*> metadata;metadata["phaseRejectCount"]=&rejects;
        if(code==codes[1] && (edge==failed || (cut && edge==Fixture::edgesFor("B",2))))
            metadata["zhangPppDeweighted"]=&rejects;
        measurement.metaDataMaps.push_back(std::move(metadata));
    }
    const int rows=measurement.obsKeys.size();
    measurement.H=MatrixXd::Zero(rows,fixture.state.x.size());
    measurement.R=MatrixXd::Identity(rows,rows);
    std::string reason;
    require(zhangFinalizeGraphSignalObservationAdmission(fixture.state,E_Sys::GPS,measurement,reason),
        "post-QC fixture finalization failed: "+reason);
    ZhangGraphIntegerContext before;require(zhangGraphIntegerContext(fixture.state,E_Sys::GPS,before),"post-QC source missing");
    const auto oldIndex=fixture.state.kfIndexMap;const VectorXd oldX=fixture.state.x;const MatrixXd oldP=fixture.state.P;
    const auto oldMetadata=fixture.state.metaDataMap;const auto oldSequence=fixture.state.factorCommitSequence;
    if(rejectCallback) fixture.state.exactStateTransformCallback=[](const KFState&,GTime,
        const std::map<KFKey,int>&,const std::map<KFKey,int>&,const SparseMatrix<double>&,
        const std::string&,const VectorXd&,const MatrixXd&,std::uint64_t,std::uint64_t){return false;};
    const bool accepted=zhangTryPostQcCertifiedGraphReroute(fixture.trace,fixture.state,E_Sys::GPS,reason);
    ZhangGraphIntegerContext after;require(zhangGraphIntegerContext(fixture.state,E_Sys::GPS,after),"post-QC result missing");
    require(!after.signalObservationEdges.at(int(codes[1])).contains(failed),"deweighted signal was readmitted");
    require(after.signalObservationEdges.at(int(codes[0])).contains(failed),"passing signal was unnecessarily readmitted/removed");
    require(after.signalArcVersions==before.signalArcVersions,"temporary QC changed physical arc identities");
    require(zhangCurrentEpochEvidenceRegistry().empty(),"historical reroute granted CURRENT AR authority");
    if(!history || cut || rejectCallback)
    {
        require(!accepted && after.basis.treeEdges==before.basis.treeEdges &&
            fixture.state.kfIndexMap==oldIndex && fixture.state.metaDataMap==oldMetadata &&
            fixture.state.factorCommitSequence==oldSequence &&
            (fixture.state.x-oldX).isZero() && (fixture.state.P-oldP).isZero(),
            "failed post-QC candidate changed authoritative posterior, metadata or graph: "+reason);
        require(reason.find(!history?"UNCERTIFIED":cut?"CUT_HAS_NO_FULL_RANK":"TRANSACTION_REJECTED")!=std::string::npos,
            "post-QC rejection did not report the actual unavailable direction/cut/transaction: "+reason);
        return;
    }
    require(accepted && !after.basis.treeEdges.contains(failed) && after.basis.treeEdges.contains(Fixture::edgesFor("B",2)),
        "post-QC certified shared path not rerouted: "+reason);
    require(after.basis.edges==before.basis.edges && after.signalObservationEdges==before.signalObservationEdges,
        "post-QC exchange deleted FLOAT coordinates or fabricated current observations");
    // Literal independent oracle: changing B's datum subtracts lambda*N_B2
    // from EVERY B-only product satellite; preserve clock/phase cross-Q too.
    std::vector<std::pair<KFKey,VectorXd>> literal;
    for(const auto code:codes) for(const auto& satellite:before.basis.satellites)
    {
        const auto key=zhangTransportSatelliteKey(code,satellite);
        VectorXd row=VectorXd::Zero(oldX.size());row(oldIndex.at(key))=1;
        if(satellite.prn>=3)row(oldIndex.at(zhangTransportAmbiguityKey(code,Fixture::edgesFor("B",2))))=
            -CLIGHT/(code==codes[0]?FREQ1:FREQ2);
        literal.push_back({key,std::move(row)});
    }
    for(const auto& satellite:before.basis.satellites)
    {VectorXd row=VectorXd::Zero(oldX.size());const auto key=clockKey(satellite);row(oldIndex.at(key))=1;literal.push_back({key,std::move(row)});}
    for(const auto& [key,row]:literal)
    {
        const int i=fixture.state.kfIndexMap.at(key);
        require(std::abs(fixture.state.x(i)-row.dot(oldX))<1e-10,"post-QC mean violates literal product oracle");
        for(const auto& [other,otherRow]:literal)
            require(std::abs(fixture.state.P(i,fixture.state.kfIndexMap.at(other))-row.dot(oldP*otherRow))<1e-9,
                "post-QC complete product/clock cross-Q violates literal oracle");
    }
    require(fixture.trace.str().find("affected_products=4")!=std::string::npos,
        "shared-path multi-target impact was not diagnosed");
    require(!zhangTryPostQcCertifiedGraphReroute(fixture.trace,fixture.state,E_Sys::GPS,reason) &&
        reason=="POST_QC_COMMITTED_TREE_STILL_ADMITTED","repeating reroute did not remain idempotent");
}

static void runEquivalentHistoryCase()
{
    Fixture fixture(false,"CERTIFIED_TREE_EQUIVALENT_BASES");
    fixture.installCurrentEvidenceAndCommitContract();
    const std::string owner=fixture.state.metaDataMap.at(ZHANG_CHECKPOINT_RUNTIME_ID_METADATA);
    const auto source=zhangCurrentEpochEvidenceRegistry().at({owner,E_Sys::GPS});
    const auto oldMetadata=fixture.state.metaDataMap;
    for(const int multiplier:{2,3,5})
    {
        fixture.state.time=fixture.state.time+30;
        zhangBeginEpochOnlyIntegerTransaction(fixture.state);
        const auto epoch=static_cast<long int>(fixture.state.time.bigTime);
        const auto root=zhangR51NumericRoot(fixture.state.x,fixture.state.P);
        auto& evidence=zhangCurrentEpochEvidenceRegistry()[{owner,E_Sys::GPS}];
        require(evidence.bind(owner,E_Sys::GPS,epoch,"NEW_CURRENT_POSTERIOR","NEW_CURRENT_CHART",source.chart,root),
            "equivalent current basis bind failed");
        evidence.floatPosteriorPsdValidated=true;evidence.floatPosteriorPsdRootId=root;
        evidence.floatPosteriorPsdStateOrderId=zhangCurrentEvidenceStateOrderId(fixture.state.kfIndexMap);
        for(std::size_t i=0;i<source.rows.size();++i)
        {
            auto physical=source.rows[i];auto numeric=source.numericConstraints.at(i).rowCycles;
            for(auto& [id,value]:physical)value*=multiplier;
            for(auto& [key,value]:numeric)value*=multiplier;
            const auto value=source.values[i]*multiplier;
            auto proof=std::make_shared<ZhangIntegerDecisionProof>();
            proof->id=owner+"|EQUIVALENT_CURRENT|"+std::to_string(epoch)+"|"+std::to_string(i);
            proof->originalStatement="INDEPENDENT_FIXTURE_MULTIPLE_OF_OLD_CYCLE_NOT_REAL_AR";
            proof->observationProvenance=owner+"|PPP_COMMITTED|E"+std::to_string(epoch);
            proof->conditionalFailureBound=0;
            require(evidence.add(physical,value,{proof}) && evidence.addNumericConstraint(numeric,value,physical),
                "new equivalent basis fixture was not fully bound");
        }
        std::string reason;ZhangPreparedCertifiedTreeTransportHistory prepared;
        require(zhangPrepareCertifiedTreeTransportHistory(fixture.trace,fixture.state,
            {{E_Sys::GPS,{"NEW_CURRENT_POSTERIOR","NEW_CURRENT_CHART"}}},prepared,&reason),
            "equivalent basis staging failed: "+reason);
        require(prepared.sealedDirectories.empty() && fixture.state.metaDataMap==oldMetadata,
            "equivalent per-epoch bases grew facts/scopes/proof payloads or changed the old seal");
        require(!zhangCommitCertifiedTreeTransportHistory(fixture.trace,fixture.state,prepared,&reason) &&
            reason=="NO_PREPARED_TRANSPORT_FACTS" && fixture.state.metaDataMap==oldMetadata,
            "redundant old evidence was mislabeled as a fresh committed certificate");
    }
    // Keep the current physical/numeric binding internally consistent while
    // deliberately contradicting the old constant. The history import must
    // reject this synthetic bad accepted value, never replace the old seal.
    auto& evidence=zhangCurrentEpochEvidenceRegistry().at({owner,E_Sys::GPS});
    ++evidence.values.front();++evidence.numericConstraints.front().value;
    ++evidence.physicalNumericBindings.front().value;
    std::string reason;ZhangPreparedCertifiedTreeTransportHistory rejected;
    require(!zhangPrepareCertifiedTreeTransportHistory(fixture.trace,fixture.state,
        {{E_Sys::GPS,{"NEW_CURRENT_POSTERIOR","NEW_CURRENT_CHART"}}},rejected,&reason) &&
        reason=="TRANSPORT_COMMIT_CONFLICT_UNBOUND_PROOF_OR_CAPACITY" && fixture.state.metaDataMap==oldMetadata,
        "conflicting implied constant replaced accepted transport-only history");
}

static void runProductionCase(bool reverse)
{
    Fixture fixture(reverse,reverse?"CERTIFIED_TREE_PRODUCTION_REVERSE":"CERTIFIED_TREE_PRODUCTION_FORWARD");
    fixture.installCurrentEvidenceAndCommitContract();
    const auto oldKeys=fixture.state.kfIndexMap;
    const VectorXd oldX=fixture.state.x;
    const MatrixXd oldP=fixture.state.P;
    for(const auto code:codes)
        fixture.receivers.at("A").satStatMap[SatSys("G02")].sigStatMap[
            ft2string(code2Freq[E_Sys::GPS][code])].tracking.pendingBreak=true;
    fixture.state.time=fixture.state.time+30;
    updateZhangFullRankReferences(fixture.trace,fixture.receivers,fixture.state);
    ZhangGraphIntegerContext after;
    require(zhangGraphIntegerContext(fixture.state,E_Sys::GPS,after),"post-slip graph missing");
    require(after.basis.treeEdges.contains(Fixture::edgesFor("B",2)) &&
        !after.basis.treeEdges.contains(Fixture::edgesFor("C",2)),"real graph updater ignored the lower-quality certified replacement");
    require(zhangCurrentEpochEvidenceRegistry().empty(),"old RHS leaked through the epoch-only reset");
    std::vector<std::pair<KFKey,std::vector<double>>> productRows;
    for(const auto code:codes)
    {
        require(after.signalArcVersions.at(Fixture::edgesFor("A",2)).at(static_cast<int>(code)) >
            fixture.before.signalArcVersions.at(Fixture::edgesFor("A",2)).at(static_cast<int>(code)),"physical slip reused the old signal arc version");
        const double lambda=CLIGHT/(code==E_ObsCode::L1C?FREQ1:FREQ2);
        for(int satellite:{1,2})
        {
            const auto key=zhangTransportSatelliteKey(code,SatSys(E_Sys::GPS,satellite));
            std::vector<double> literal(oldX.size()); literal[oldKeys.at(key)]=1;
            if(satellite==2) literal[oldKeys.at(zhangTransportAmbiguityKey(code,Fixture::edgesFor("B",2)))]=lambda;
            productRows.push_back({key,std::move(literal)});
        }
    }
    for(int satellite:{1,2})
    {
        const auto key=clockKey(SatSys(E_Sys::GPS,satellite));
        std::vector<double> literal(oldX.size()); literal[oldKeys.at(key)]=1;
        productRows.push_back({key,std::move(literal)});
    }
    for(const auto& [key,row]:productRows)
    {
        double mean=0;
        for(int i=0;i<oldX.size();++i) mean+=row[i]*oldX(i);
        const int actual=fixture.state.kfIndexMap.at(key);
        require(std::abs(fixture.state.x(actual)-mean)<1e-10,"production product mean differs from literal pre-event oracle");
        for(const auto& [other,otherRow]:productRows)
        {
            double covariance=0;
            for(int i=0;i<oldX.size();++i) for(int j=0;j<oldX.size();++j)
                covariance+=row[i]*oldP(i,j)*otherRow[j];
            require(std::abs(fixture.state.P(actual,fixture.state.kfIndexMap.at(other))-covariance)<1e-10,
                "production complete phase/clock Q differs from literal pre-event oracle");
        }
    }
    // Crucial: the FLOAT transform used the stochastic 11.8, not the archived
    // accepted integer 12. History must not silently condition this direction.
    const int phase=fixture.state.kfIndexMap.at(zhangTransportSatelliteKey(E_ObsCode::L1C,SatSys("G02")));
    require(std::abs(fixture.state.x(phase)-(12+(CLIGHT/FREQ1)*12))>0.01,"history fixed an unfixed FLOAT direction");
    ZhangProductIntegerTransition oldCycle;
    oldCycle.valid=true; oldCycle.observable=E_ObsCode::L1C;
    oldCycle.physicalEdges={Fixture::edgesFor("B",2),Fixture::edgesFor("B",1),Fixture::edgesFor("A",2),Fixture::edgesFor("A",1)};
    oldCycle.coefficients={1,-1,-1,1};
    for(const auto& edge:oldCycle.physicalEdges)
        oldCycle.physicalArcVersions.push_back(fixture.before.signalArcVersions.at(edge).at(static_cast<int>(E_ObsCode::L1C)));
    const auto transport=zhangQueryCertifiedHistoricalTransport(fixture.state,E_Sys::GPS,oldCycle);
    require(transport.proven && transport.value==12 && !transport.currentArAuthority,
        "old physical certificate did not remain usable for transport only");
    oldCycle.physicalArcVersions[2]=after.signalArcVersions.at(Fixture::edgesFor("A",2)).at(static_cast<int>(E_ObsCode::L1C));
    require(!zhangQueryCertifiedHistoricalTransport(fixture.state,E_Sys::GPS,oldCycle).proven,
        "old relation was inherited by the new slipped physical arc");
    oldCycle.physicalEdges[2]=Fixture::edgesFor("D",2);
    require(!zhangQueryCertifiedHistoricalTransport(fixture.state,E_Sys::GPS,oldCycle).proven,
        "unconnected new receiver received an invented old endpoint");
}

static void runDisconnectedProductionCase()
{
    Fixture fixture(false,"CERTIFIED_TREE_PRODUCTION_SPLIT");
    fixture.installCurrentEvidenceAndCommitContract();
    for(const auto receiver:{"A","B","C"})
    {
        auto& observations=fixture.receivers.at(receiver).obsList;
        observations.erase(std::remove_if(observations.begin(),observations.end(),[](const auto& item)
        { const auto observation=std::dynamic_pointer_cast<GObs>(item); return observation && observation->Sat==SatSys("G02"); }),observations.end());
    }
    fixture.addObservation(Fixture::edgesFor("D",2));
    fixture.state.time=fixture.state.time+30;
    updateZhangFullRankReferences(fixture.trace,fixture.receivers,fixture.state);
    require(fixture.trace.str().find("all_old_endpoints_certified=0")!=std::string::npos,
        "disconnected real graph claimed certified continuity of all old endpoints");
    require(zhangCurrentEpochEvidenceRegistry().empty(),"component split restored old current RHS");
    ZhangProductIntegerTransition across;
    across.valid=true; across.observable=E_ObsCode::L1C;
    across.physicalEdges={Fixture::edgesFor("A",1),Fixture::edgesFor("D",2)};
    across.physicalArcVersions={0,0}; across.coefficients={1,-1};
    require(!zhangQueryCertifiedHistoricalTransport(fixture.state,E_Sys::GPS,across).proven,
        "real network split invented an integer relation to the new component");
}

static void runAllArcsRetiredFreshObservationCase(bool reverse)
{
    Fixture fixture(reverse,reverse?"CERTIFIED_TREE_ALL_RETIRED_REVERSE":"CERTIFIED_TREE_ALL_RETIRED_FORWARD");
    fixture.installCurrentEvidenceAndCommitContract();
    const auto traceStart=fixture.trace.str().size();
    // All observations remain present and QC-admissible. Every signal has a
    // genuine pending break; version numbers have not yet been consumed when
    // the production candidate selector runs. Old B-G02 is a historical chord
    // with a certified old cycle, but its fresh physical arc must not inherit it.
    for(const auto& edge:fixture.before.basis.edges)
        for(const auto code:codes)
            fixture.receivers.at(edge.receiver).satStatMap[edge.satellite].sigStatMap[
                ft2string(code2Freq[E_Sys::GPS][code])].tracking.pendingBreak=true;
    fixture.state.time=fixture.state.time+30;
    updateZhangFullRankReferences(fixture.trace,fixture.receivers,fixture.state);
    const auto eventTrace=fixture.trace.str().substr(traceStart);
    std::istringstream lines(eventTrace);std::string line;bool sawSelector=false;
    while(std::getline(lines,line))
        if(line.find("ZHANG_CERTIFIED_TREE_TRANSPORT ")!=std::string::npos)
        {
            sawSelector=true;
            require(line.find(" certified_replacements=0 ")!=std::string::npos,
                "pending-discontinuity chord was ranked as an old-certified replacement");
        }
    require(sawSelector,"all-retired production selector path was not exercised");
    require(eventTrace.find("reason=ALL_CURRENT_ARCS_RETIRED")!=std::string::npos,
        "FLOAT fresh-arc restart path was not exercised");
    ZhangGraphIntegerContext after;
    require(zhangGraphIntegerContext(fixture.state,E_Sys::GPS,after) && after.basis.connected &&
        after.basis.edges==fixture.before.basis.edges,
        "rejecting old certification removed legal fresh observations from FLOAT restart");
    require(fixture.state.x.allFinite() && fixture.state.P.allFinite(),"fresh FLOAT restart produced invalid posterior");
    require(zhangCurrentEpochEvidenceRegistry().empty(),"all-retired restart applied historical current RHS");
    for(const auto& edge:fixture.before.basis.edges)
        for(const auto code:codes)
            require(after.signalArcVersions.at(edge).at(static_cast<int>(code))>
                fixture.before.signalArcVersions.at(edge).at(static_cast<int>(code)),
                "all-retired restart failed to consume a physical signal break");
    ZhangProductIntegerTransition cycle;
    cycle.valid=true;cycle.observable=E_ObsCode::L1C;
    cycle.physicalEdges={Fixture::edgesFor("B",2),Fixture::edgesFor("B",1),Fixture::edgesFor("A",2),Fixture::edgesFor("A",1)};
    cycle.coefficients={1,-1,-1,1};
    for(const auto& edge:cycle.physicalEdges)
        cycle.physicalArcVersions.push_back(fixture.before.signalArcVersions.at(edge).at(static_cast<int>(E_ObsCode::L1C)));
    const auto old=zhangQueryCertifiedHistoricalTransport(fixture.state,E_Sys::GPS,cycle);
    require(old.proven && old.value==12 && !old.currentArAuthority,"restart destroyed valid old transport-only facts");
    for(std::size_t i=0;i<cycle.physicalEdges.size();++i)
        cycle.physicalArcVersions[i]=after.signalArcVersions.at(cycle.physicalEdges[i]).at(static_cast<int>(E_ObsCode::L1C));
    require(!zhangQueryCertifiedHistoricalTransport(fixture.state,E_Sys::GPS,cycle).proven,
        "all-retired fresh cycle inherited the old accepted integer RHS");
}

int main()
{
    setenv("ZHANG_R51_ENABLE","1",1); setenv("ZHANG_R51_RATIO_ONLY","1",1);
    setenv("ZHANG_R51_RATIO_THRESHOLD","2",1);
    setenv("ZHANG_R51_EPOCH_ONLY_AR","1",1); setenv("ZHANG_R51_SINGLE_BASIS_HANDOFF","1",1);
    setenv("ZHANG_R51_PRODUCT_REPAIR","1",1);
    acsConfig.zhangFullRank.product_basis_mode="SINGLE_S_BASIS_NATIVE";
    acsConfig.zhangPppAr.output_diagnostics=true;
    acsConfig.zhangFullRank.enable=true; acsConfig.process_sys[E_Sys::GPS]=true;
    auto& options=acsConfig.zhangFullRank.sysOpts[E_Sys::GPS];
    options.use_spanning_tree=true; options.reference_receiver="A";
    options.baseline_observables=codes; options.product_core_min_satellite_support=0;
    options.state_edge_grace_epochs=0;
    runProductionCase(false); runProductionCase(true); runDisconnectedProductionCase();
    runAllArcsRetiredFreshObservationCase(false);runAllArcsRetiredFreshObservationCase(true);
    runEquivalentHistoryCase();
    runPostQcRerouteCase(true,false,false);
    runPostQcRerouteCase(false,false,false);
    runPostQcRerouteCase(true,true,false);
    runPostQcRerouteCase(true,false,true);
    std::cout<<"CERTIFIED_TREE_TRANSPORT_PRODUCTION_CONTRACT_PASS FIXTURE_NOT_REAL_PRODUCT_PUBLICATION\n";
}
