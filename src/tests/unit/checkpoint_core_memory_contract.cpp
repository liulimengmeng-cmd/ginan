// Independently freezes the legacy full archive/equality as the reference.
// Link to the real candidate common/zhangCheckpoint.cpp, never a copied
// comparator implementation. This program does not invoke PEA or alter files.
#include "common/zhangCheckpoint.hpp"
#include "common/receiver.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
int assertions = 0;
void check(bool accepted, const std::string& reason)
{
    ++assertions;
    if (!accepted) throw std::runtime_error(reason);
}

// The exact pre-optimisation definition, intentionally independent of the
// new metadata comparator. Every field, order, key annotation and byte counts.
std::string legacyWithoutTime(const ZhangCheckpointKfCore& s)
{
    std::ostringstream out(std::ios::binary | std::ios::out);
    boost::archive::binary_oarchive ar(out, boost::archive::no_header);
    ar & s.x; ar & s.P; ar & s.dx; ar & s.prefitRatios; ar & s.postfitRatios;
    ar & s.kfIndexMap; ar & s.stateTransitionMap; ar & s.gaussMarkovTauMap;
    ar & s.gaussMarkovMuMap; ar & s.procNoiseMap; ar & s.initNoiseMap;
    ar & s.sigmaMaxMap; ar & s.outageLimitMap; ar & s.exponentialNoiseMap;
    ar & s.pseudoStateMap; ar & s.pseudoParentMap; ar & s.errorCountMap;
    ar & s.filterChunkMap; ar & s.metaDataMap; ar & s.lsqRequired;
    ar & s.sigmaPass; ar & s.chiQCPass; ar & s.chi2; ar & s.dof;
    ar & s.chi2PerDof; ar & s.qc; ar & s.id; ar & s.rtsBasename;
    ar & s.outputResiduals; ar & s.outputMongoMeasurements;
    ar & s.statisticsMap; ar & s.statisticsMapSum;
    return out.str();
}
bool legacyEqual(const ZhangCheckpointKfCore& a, const ZhangCheckpointKfCore& b)
{
    return a.time.bigTimeBytes == b.time.bigTimeBytes &&
        legacyWithoutTime(a) == legacyWithoutTime(b);
}
#ifdef ZHANG_TEST_BASELINE_NO_STREAMING
// An explicitly labelled baseline build links unchanged old production
// objects, which do not export the new streaming APIs. No source/object is
// patched or substituted; only this test main selects old public operations.
bool oldMatchesState(const ZhangCheckpointKfCore& core,const KFState& state)
{ return zhangCheckpointKfCoreBitwiseEqual(core,captureZhangCheckpointKfCore(state)); }
std::string oldCoreDigest(const ZhangCheckpointKfCore& core)
{ return zhangCheckpointSha256(serializeZhangCheckpointSectionPayload(core)); }
std::string oldStateDigest(const KFState& state)
{ return oldCoreDigest(captureZhangCheckpointKfCore(state)); }
#define zhangCheckpointKfCoreMatchesState oldMatchesState
#define zhangCheckpointKfCoreSha256Streaming oldCoreDigest
#define zhangCheckpointKfStateCoreSha256Streaming oldStateDigest
#endif
double fromBits(std::uint64_t bits)
{
    double value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

Receiver receiver;
KFState fixture(int n)
{
    KFState state;
    state.time.bigTime = 1405213230.L;
    state.x = VectorXd::LinSpaced(n, -2.75, 4.25);
    state.P.resize(n,n);
    for (int c=0; c<n; ++c) for (int r=0; r<n; ++r)
        state.P(r,c) = (r==c ? 2.0 : .125 / (1+std::abs(r-c)));
    state.dx = VectorXd::Constant(n, .03125);
    state.prefitRatios = VectorXd::Constant(3, .5);
    state.postfitRatios = VectorXd::Constant(5, .25);
    state.kfIndexMap.clear();
    receiver.id="RECEIVER_STABLE_ID";
    KFKey first;
    first.type=KF::AMBIGUITY; first.str="MAL2"; first.Sat=SatSys("G25");
    first.num=int(E_ObsCode::L1C); first.comment="FULL_COMMENT";
    first.estimatedTime.bigTime=state.time.bigTime-123.125L;
    first.rec_ptr=&receiver;
    auto second=first; second.Sat=SatSys("E11"); second.str="NICO";
    for (int i=0; i<n; ++i)
    {
        auto key=first; key.num=i; state.kfIndexMap[key]=i;
    }
    state.stateTransitionMap[first][second][0]=1.25;
    state.stateTransitionMap[first][first][2]=-.125;
    state.gaussMarkovTauMap[first]=7;
    state.gaussMarkovMuMap[first]=.125;
    state.procNoiseMap[first]=.25;
    state.initNoiseMap[first]=.5;
    state.sigmaMaxMap[first]=4;
    state.outageLimitMap[first]=90;
    state.exponentialNoiseMap[first]={.125,3600};
    state.pseudoStateMap[first][second]=.25;
    state.pseudoParentMap[first]=second;
    state.errorCountMap[first]=2;
    FilterChunk chunk;
    chunk.id="FULL_CHUNK"; chunk.begX=1; chunk.numX=2;
    chunk.begH=3; chunk.numH=4;
    state.filterChunkMap["main"]=chunk;
    state.metaDataMap[ZHANG_CHECKPOINT_RUNTIME_ID_METADATA]="LOCAL_MEMORY_NO_PEA";
    // Deliberately crosses the streaming buffer boundary.
    state.metaDataMap["large_opaque_metadata"]=std::string(140003,'q');
    state.lsqRequired=true; state.sigmaPass=true; state.chiQCPass=true;
    state.chi2=2.5; state.dof=7; state.chi2PerDof=.125; state.qc=.25;
    state.id="FULL_ID"; state.rts_basename="FULL_RTS";
    state.output_residuals=true; state.outputMongoMeasurements=true;
    state.statisticsMap["one"]=3; state.statisticsMapSum["total"]=5;
    return state;
}

void contract()
{
    auto state=fixture(7);
    const auto core=captureZhangCheckpointKfCore(state);
    const auto frozen=serializeZhangCheckpointSectionPayload(core);
    const auto legacyDigest=zhangCheckpointSha256(frozen);
    check(zhangCheckpointKfCoreMatchesState(core,state), "live core positive");
    check(zhangCheckpointKfCoreBitwiseEqual(core,core), "core positive");
    check(zhangCheckpointKfCoreSha256Streaming(core)==legacyDigest,"stream snapshot SHA differs");
    check(zhangCheckpointKfStateCoreSha256Streaming(state)==legacyDigest,"stream live SHA differs");

    using Mutation=std::pair<std::string,std::function<void(ZhangCheckpointKfCore&)>>;
    std::vector<Mutation> mutations={
        {"time",[](auto& s){s.time.bigTimeBytes[0]^=1;}},
        {"x",[](auto& s){s.x(1)+=.125;}},
        {"P diagonal",[](auto& s){s.P(1,1)+=.125;}},
        {"P cross system",[](auto& s){s.P(1,5)+=.125;}},
        {"P asymmetric half",[](auto& s){s.P(5,1)+=.125;}},
        {"dx",[](auto& s){s.dx(1)+=.125;}},
        {"prefit",[](auto& s){s.prefitRatios(1)+=.125;}},
        {"postfit",[](auto& s){s.postfitRatios(1)+=.125;}},
        {"vector dimension",[](auto& s){s.x.conservativeResize(6);}},
        {"matrix shape",[](auto& s){s.P.resize(1,49);}},
        {"index",[](auto& s){s.kfIndexMap.begin()->second++;}},
        {"key comment",[](auto& s){auto n=s.kfIndexMap.extract(s.kfIndexMap.begin());n.key().comment+="changed";s.kfIndexMap.insert(std::move(n));}},
        {"key estimate epoch",[](auto& s){auto n=s.kfIndexMap.extract(s.kfIndexMap.begin());n.key().estimatedTime.bigTimeBytes[1]^=1;s.kfIndexMap.insert(std::move(n));}},
        {"key receiver binding",[](auto& s){auto n=s.kfIndexMap.extract(s.kfIndexMap.begin());n.key().receiverPointerId+="changed";s.kfIndexMap.insert(std::move(n));}},
        {"transition",[](auto& s){s.stateTransitionMap.begin()->second.begin()->second.begin()->second+=.125;}},
        {"tau",[](auto& s){s.gaussMarkovTauMap.begin()->second+=1;}},
        {"mu",[](auto& s){s.gaussMarkovMuMap.begin()->second+=1;}},
        {"Q",[](auto& s){s.procNoiseMap.begin()->second+=1;}},
        {"init",[](auto& s){s.initNoiseMap.begin()->second+=1;}},
        {"sigma",[](auto& s){s.sigmaMaxMap.begin()->second+=1;}},
        {"outage",[](auto& s){s.outageLimitMap.begin()->second+=1;}},
        {"exponential value",[](auto& s){s.exponentialNoiseMap.begin()->second.value+=1;}},
        {"exponential tau",[](auto& s){s.exponentialNoiseMap.begin()->second.tau+=1;}},
        {"pseudo row",[](auto& s){s.pseudoStateMap.begin()->second.begin()->second+=1;}},
        {"pseudo parent",[](auto& s){s.pseudoParentMap.begin()->second.comment+="changed";}},
        {"error count",[](auto& s){s.errorCountMap.begin()->second++;}},
        {"chunk id",[](auto& s){s.filterChunkMap.begin()->second.id+="changed";}},
        {"chunk begX",[](auto& s){s.filterChunkMap.begin()->second.begX++;}},
        {"chunk numX",[](auto& s){s.filterChunkMap.begin()->second.numX++;}},
        {"chunk begH",[](auto& s){s.filterChunkMap.begin()->second.begH++;}},
        {"chunk numH",[](auto& s){s.filterChunkMap.begin()->second.numH++;}},
        {"metadata",[](auto& s){s.metaDataMap["added"]="changed";}},
        {"lsq",[](auto& s){s.lsqRequired=!s.lsqRequired;}},
        {"sigma pass",[](auto& s){s.sigmaPass=!s.sigmaPass;}},
        {"chi pass",[](auto& s){s.chiQCPass=!s.chiQCPass;}},
        {"chi2",[](auto& s){s.chi2+=.125;}},
        {"dof",[](auto& s){s.dof++;}},
        {"chi2 per dof",[](auto& s){s.chi2PerDof+=.125;}},
        {"qc",[](auto& s){s.qc+=.125;}},
        {"id",[](auto& s){s.id+="changed";}},
        {"rts",[](auto& s){s.rtsBasename+="changed";}},
        {"residual flag",[](auto& s){s.outputResiduals=!s.outputResiduals;}},
        {"mongo flag",[](auto& s){s.outputMongoMeasurements=!s.outputMongoMeasurements;}},
        {"statistics",[](auto& s){s.statisticsMap.begin()->second++;}},
        {"sum statistics",[](auto& s){s.statisticsMapSum.begin()->second++;}},
    };
    for (const auto& [name,mutate]:mutations)
    {
        auto changed=core;
        mutate(changed);
        check(!legacyEqual(core,changed),"bad independent negative "+name);
        check(!zhangCheckpointKfCoreBitwiseEqual(core,changed),"missed core field "+name);
        check(!zhangCheckpointKfCoreMatchesState(changed,state),"missed live field "+name);
        check(zhangCheckpointKfCoreSha256Streaming(changed)==
            zhangCheckpointSha256(serializeZhangCheckpointSectionPayload(changed)),
            "changed stream SHA "+name);
    }
    for (const auto bits:std::vector<std::uint64_t>{0,UINT64_C(0x8000000000000000),
        UINT64_C(0x7ff8000000000001),UINT64_C(0x7ff8000000000002),
        UINT64_C(0x7ff0000000000000),UINT64_C(0xfff0000000000000)})
    {
        state.x(1)=fromBits(bits);
        const auto exact=captureZhangCheckpointKfCore(state);
        check(zhangCheckpointKfCoreMatchesState(exact,state),"equal special bits rejected");
        check(zhangCheckpointKfCoreSha256Streaming(exact)==
            zhangCheckpointSha256(serializeZhangCheckpointSectionPayload(exact)),"special snapshot SHA");
        check(zhangCheckpointKfStateCoreSha256Streaming(state)==
            zhangCheckpointSha256(serializeZhangCheckpointSectionPayload(exact)),"special live SHA");
        auto changed=exact; changed.x(1)=fromBits(bits ^ 1);
        check(!legacyEqual(exact,changed),"special independent negative");
        check(!zhangCheckpointKfCoreBitwiseEqual(exact,changed),"special bits missed");
        check(!zhangCheckpointKfCoreMatchesState(changed,state),"special live bits missed");
    }
    check(frozen==serializeZhangCheckpointSectionPayload(core),"source archive mutated");
    auto empty=fixture(0);
    const auto emptyCore=captureZhangCheckpointKfCore(empty);
    check(zhangCheckpointKfCoreMatchesState(emptyCore,empty),"empty core compare");
    check(zhangCheckpointKfStateCoreSha256Streaming(empty)==
        zhangCheckpointSha256(serializeZhangCheckpointSectionPayload(emptyCore)),"empty SHA");
    std::cout << "{\"passed\":true,\"assertions\":" << assertions
        << ",\"fixture_dimension\":7,\"scientific_acceptance\":false}\n";
}

void runBenchmark(const std::string& mode,const KFState& state,
    const ZhangCheckpointKfCore& before,const std::string& artifact)
{
    bool equal=false;
    std::string beforeDigest,afterDigest;
    if (mode=="old")
    {
        const auto after=captureZhangCheckpointKfCore(state);
#ifdef ZHANG_TEST_BASELINE_NO_STREAMING
        equal=zhangCheckpointKfCoreBitwiseEqual(before,after);
#else
        // This optional counterfactual uses the frozen pre-change definition;
        // production old/new comparison should instead link separate binaries.
        equal=legacyEqual(before,after);
#endif
        beforeDigest=zhangCheckpointSha256(serializeZhangCheckpointSectionPayload(before));
        afterDigest=zhangCheckpointSha256(serializeZhangCheckpointSectionPayload(after));
    }
    else if (mode=="new")
    {
        equal=zhangCheckpointKfCoreMatchesState(before,state);
        beforeDigest=zhangCheckpointKfCoreSha256Streaming(before);
        afterDigest=zhangCheckpointKfStateCoreSha256Streaming(state);
    }
    else throw std::runtime_error("mode must be old or new");
    check(equal && beforeDigest.size()==64 && beforeDigest==afterDigest,"benchmark closure failed");
    if (!artifact.empty())
    {
        if (std::filesystem::exists(artifact)) throw std::runtime_error("refusing to overwrite artifact");
        std::ofstream output(artifact,std::ios::binary);
        output.exceptions(std::ios::badbit | std::ios::failbit);
        boost::archive::binary_oarchive archive(output,boost::archive::no_header);
        archive << before;
    }
    std::cout << "{\"passed\":true,\"mode\":\"" << mode
        << "\",\"dimension\":" << state.x.size() << ",\"core_sha256\":\"" << beforeDigest
        << "\",\"scientific_acceptance\":false}\n";
}

void benchmark(const std::string& mode,int n,const std::string& artifact)
{
    auto state=fixture(n);
    const auto before=captureZhangCheckpointKfCore(state);
    runBenchmark(mode,state,before,artifact);
}

void realCore(const std::string& mode,const std::string& inputPath,
    const std::string& artifact)
{
    ZhangCheckpointKfCore before;
    {
        std::ifstream input(inputPath,std::ios::binary);
        if (!input) throw std::runtime_error("real core open failed");
        boost::archive::binary_iarchive archive(input,boost::archive::no_header);
        archive >> before;
        if (input.peek()!=std::char_traits<char>::eof())
            throw std::runtime_error("real core trailing bytes");
    }
    if (before.x.size()<1 || before.x.size()>7196)
        throw std::runtime_error("real core dimension outside bounded test range");
    std::map<std::string,Receiver> receiverBindings;
    KFState state;
    std::string failure;
    const auto resolve=[&](const std::string& id)->Receiver*
    {
        auto& bound=receiverBindings[id];
        bound.id=id;
        return &bound;
    };
    check(restoreZhangCheckpointKfCoreWithReceiverResolver(before,state,resolve,&failure),
        "real core restore failed: "+failure);
    runBenchmark(mode,state,before,artifact);
}
}

int main(int argc,char** argv)
{
    try
    {
        if (argc==1) contract();
        else if ((argc==4 || argc==5) && std::string(argv[1])=="--benchmark")
        {
            const int n=std::stoi(argv[3]);
            if (n<1 || n>7196) throw std::runtime_error("bounded dimension must be 1..7196");
            benchmark(argv[2],n,argc==5 ? argv[4] : "");
        }
        else if ((argc==4 || argc==5) && std::string(argv[1])=="--real-core")
            realCore(argv[2],argv[3],argc==5 ? argv[4] : "");
        else throw std::runtime_error("usage: contract [--benchmark old|new dimension [full_core.bin]] or [--real-core old|new core.payload [full_core.bin]]");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "CONTRACT_FAILED:" << error.what() << '\n';
        return 1;
    }
}
