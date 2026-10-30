#include "common/zhangUserGonsis.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
int checks=0;
void require(bool condition,const std::string& reason)
{
    ++checks;
    if (!condition) throw std::runtime_error(reason);
}

struct Fixture
{
    VectorXd mean;
    MatrixXd Q;
    std::vector<ZhangUserGonsisSatellite> satellites;
};

Fixture fixture(int count)
{
    Fixture out;
    out.mean.resize(2*count);
    out.Q=MatrixXd::Identity(2*count,2*count)*0.04;
    for (int i=0;i<count;++i)
    {
        out.mean(2*i)=10+0.03*i;
        out.mean(2*i+1)=out.mean(2*i)-(2+i%3+0.02*i);
        out.satellites.push_back({"G"+std::to_string(i+1),"ARC_STATE_PRODUCT_"+std::to_string(i),
                                 1.2-0.1*i,1,double(i+1),true});
    }
    // Nonzero inter-satellite covariance; independent search-Q equality must
    // include it, not only the diagonal or a set of satellite names.
    out.Q(2,4)=out.Q(4,2)=0.003;
    return out;
}

ZhangUserGonsisIlsSolution rejectedComplete(const VectorXd& mean,const MatrixXd&)
{
    ZhangUserGonsisIlsSolution solution;
    solution.complete=true;
    solution.integers=mean.unaryExpr([](double value){return std::floor(value+0.5);});
    solution.bestSquaredCost=1;
    solution.secondSquaredCost=1.5;
    solution.reason="RATIO_REJECTED";
    return solution;
}

void testExactCacheKey()
{
    VectorXd mean(2),integers(2);
    mean<<0.1,-0.2; integers<<1,2;
    MatrixXd Q(2,2),wl(2,4),l1(2,4),ifr(2,4);
    Q<<0.04,0.01,0.01,0.08;
    wl<<1,-1,-1,1,0,0,1,-1;
    l1<<1,0,-1,0,0,0,1,0;
    ifr=0.1*l1+0.3*wl;
    ZhangUserGonsisExactIlsCache<ZhangUserGonsisIlsSolution> cache;
    const ZhangUserGonsisIlsKeyView key{mean,Q,wl,integers,l1,ifr};
    require(!cache.find(key),"empty cache hit");
    cache.remember(key,rejectedComplete(mean,Q));
    require(cache.size()==1 && cache.find(key),"exact complete key was not cached");
    for (int field=0;field<6;++field)
    {
        auto m=mean; auto q=Q; auto w=wl; auto z=integers; auto n=l1; auto f=ifr;
        if (field==0) m(0)=std::nextafter(m(0),1.0);
        if (field==1) q(0,1)=std::nextafter(q(0,1),1.0);
        if (field==2) w(0,0)=std::nextafter(w(0,0),2.0);
        if (field==3) z(0)+=1;
        if (field==4) n(0,0)+=1;
        if (field==5) f(0,0)=std::nextafter(f(0,0),2.0);
        require(!cache.find({m,q,w,z,n,f}),"modified full numeric/target key collided");
    }
    auto reversed=mean.reverse().eval();
    require(!cache.find({reversed,Q,wl,integers,l1,ifr}),"changed search order collided");
    auto signedMean=mean; signedMean(0)=0.0;
    ZhangUserGonsisExactIlsCache<ZhangUserGonsisIlsSolution> zeroCache;
    zeroCache.remember({signedMean,Q,wl,integers,l1,ifr},rejectedComplete(mean,Q));
    signedMean(0)=-0.0;
    require(!zeroCache.find({signedMean,Q,wl,integers,l1,ifr}),"signed-zero inputs were merged");
    ZhangUserGonsisExactIlsCache<ZhangUserGonsisIlsSolution> otherEpoch;
    require(!otherEpoch.find(key),"cache crossed frozen-call lifetime");
    auto incomplete=rejectedComplete(mean,Q); incomplete.complete=false;
    otherEpoch.remember(key,incomplete);
    require(otherEpoch.size()==0,"incomplete enumeration was cached");
}

void testRealDuplicateAndTargetedPar()
{
    auto f=fixture(7);
    f.mean(11)=f.mean(10)-4.4;
    f.mean(13)=f.mean(12)-2.45;
    f.satellites[6].metricVarianceMetres2=90;
    f.satellites[5].metricVarianceMetres2=80;
    f.satellites[4].metricVarianceMetres2=70;
    const auto originalMean=f.mean; const auto originalQ=f.Q;
    int calls=0;
    auto solver=[&](const VectorXd& mean,const MatrixXd& Q)
    {++calls; return rejectedComplete(mean,Q);};
    ZhangUserGonsisOptions oldPar;
    oldPar.targetFailedL1SubsetForPar=false;
    oldPar.reuseIdenticalIls=false;
    const auto uncached=zhangResolveUserGonsis(f.mean,f.Q,f.satellites,solver,oldPar);
    require(!uncached.accepted && calls==3 && uncached.trials.size()==3,"legacy duplicate fixture changed");
    require(uncached.trials[0].removalPool=="LEGACY_ALL_ACTIVE_L1_TARGET_IGNORED",
            "legacy pool falsely claims no actual L1 target");
    calls=0;
    oldPar.reuseIdenticalIls=true;
    int audits=0,reused=0;
    oldPar.replayTrial=[&](const ZhangUserGonsisTrial& trial,const ZhangUserGonsisReplayView& view)
    {
        ++audits; reused+=trial.ilsCacheReused;
        require(zhangUserGonsisMatrixBitsEqual(view.rawFrozenMean,originalMean),"replay source mean changed");
        require(zhangUserGonsisMatrixBitsEqual(view.rawFrozenCovariance,originalQ),"replay source Q changed");
        require(view.l1Mean.size()==4 && view.l1Covariance.rows()==4,"actual searched dimensions not captured");
        require(trial.wideLaneCandidateMeansCycles.size()==trial.wideLaneCandidateSatelliteIndices.size(),"actual WL dot means missing");
        require(view.wideLaneRows.rows()==4 && view.wlConditionedMean.size()==14,"WL/posterior view lost before move");
    };
    const auto cached=zhangResolveUserGonsis(f.mean,f.Q,f.satellites,solver,oldPar);
    require(!cached.accepted && calls==1 && cached.ilsCacheHits==2 && cached.ilsCalls==1,"duplicate ILS not reused exactly");
    require(audits==3 && reused==2,"each actual/reused trial was not replayed");
    require(cached.removedSatelliteIndices==uncached.removedSatelliteIndices,"cache changed PAR trajectory");
    for (int i=0;i<cached.trials.size();++i)
        require(cached.trials[i].ratio==uncached.trials[i].ratio && cached.trials[i].reason==uncached.trials[i].reason,
                "cached complete candidate changed decision");
    require(zhangUserGonsisMatrixBitsEqual(f.mean,originalMean) && zhangUserGonsisMatrixBitsEqual(f.Q,originalQ),"FLOAT was changed");
    ZhangUserGonsisOptions targeted;
    calls=0;
    const auto repaired=zhangResolveUserGonsis(f.mean,f.Q,f.satellites,solver,targeted);
    require(!repaired.accepted && !repaired.trials.empty(),"repaired trial fixture unavailable");
    require(repaired.trials[0].removedSatelliteIndex==4 && repaired.trials[0].removalPool=="ACTUAL_FAILED_WL_L1_ENDPOINTS",
            "PAR still removes already WL-rejected satellites");
    require(repaired.removedSatelliteIndices.front()!=5 && repaired.removedSatelliteIndices.front()!=6,
            "actual L1 target scope not respected");
    require(targeted.minimumWideLanePairs==4 && targeted.ratioThreshold==2 && targeted.wideLaneResidualCycles==0.25,
            "thresholds weakened");
    const auto secondCall=zhangResolveUserGonsis(f.mean,f.Q,f.satellites,solver,oldPar);
    require(secondCall.ilsCalls==1,"result/candidate history escaped one frozen call");
    // Insufficient WL has no failed L1 target; keep the explicit fallback.
    auto insufficient=fixture(6);
    for (int i=1;i<6;++i) insufficient.mean(2*i+1)=insufficient.mean(2*i)-2.4;
    const auto wlFailed=zhangResolveUserGonsis(insufficient.mean,insufficient.Q,insufficient.satellites,solver);
    require(wlFailed.trials[0].reason=="INSUFFICIENT_PASSED_WIDE_LANE_PAIRS" &&
            wlFailed.trials[0].removalPool=="ALL_ACTIVE_NO_FAILED_L1_TARGET","WL-only fallback mislabelled as L1 PAR");
}

void testExactCoordinateRecovery()
{
    MatrixXd basis(2,2); basis<<2,1,1,1;
    VectorXd first(2),second(2),answer;
    first<<-1,-4; second<<-11,-5;
    require(zhangUserGonsisRecoverInputInteger(basis,first,answer) && answer(0)==3 && answer(1)==-7,
            "reduced best not recovered in original input coordinates");
    require(zhangUserGonsisRecoverInputInteger(basis,second,answer) && answer(0)==-6 && answer(1)==1,
            "reduced second not recovered in original input coordinates");
    MatrixXd parity=MatrixXd::Identity(2,2); parity(0,0)=2; first<<1,1;
    require(!zhangUserGonsisRecoverInputInteger(parity,first,answer),"fractional lift rounded into integer candidate");
    require(!zhangUserGonsisRecoverInputInteger(basis.topRows(1),first.head(1),answer),"partial basis labelled unique input coordinates");
    MatrixXd singular=MatrixXd::Ones(2,2);
    require(!zhangUserGonsisRecoverInputInteger(singular,first,answer),"singular basis accepted");
    basis(0,0)=1.5;
    require(!zhangUserGonsisRecoverInputInteger(basis,first,answer),"noninteger basis accepted");
    basis=MatrixXd::Identity(2,2); basis(0,1)=9007199254740991.0; first<<-1,1;
    require(!zhangUserGonsisRecoverInputInteger(basis,first,answer),"unrepresentable integer input claimed exact double");
}

void testReplayDecimalRoundtripAndSelection()
{
    MatrixXd values(2,4);
    values<<-0.0,std::nextafter(1.0,2.0),std::numeric_limits<double>::denorm_min(),-1.2345678901234567,
            std::numeric_limits<double>::max(),1e-240,-1e180,0.0;
    std::ostringstream stream;
    stream<<std::scientific<<std::setprecision(5);
    const auto flags=stream.flags(); const auto precision=stream.precision();
    zhangUserGonsisReplayMatrix(stream,"time=2024-07-17T02:00:00 receiver=TEST trial=1","Q",values);
    require(stream.flags()==flags && stream.precision()==precision,"replay logging leaked stream formatting");
    const auto record=stream.str();
    require(record.find("rows=2 columns=4 layout=ROW_MAJOR decimal_digits=17")!=std::string::npos,"matrix shape/order/precision schema missing");
    std::istringstream csv(record.substr(record.find(" data=")+6));
    std::string token; MatrixXd decoded(2,4);
    for (int i=0;i<8;++i)
    {
        require(bool(std::getline(csv,token,',')),"roundtrip data truncated");
        // strtod handles subnormals without stod's spurious underflow exception.
        decoded(i/4,i%4)=std::strtod(token.c_str(),nullptr);
    }
    require(zhangUserGonsisMatrixBitsEqual(values,decoded),"17-digit matrix log is not bit-exact roundtrip");
    ZhangUserGonsisReplayControls controls;
    controls.receivers.insert("JPLM");
    require(!controls.select("R","JPLM","2024-07-17T02:00:00") && controls.counters.empty(),"default disabled replay created history");
    controls.enabled=true;
    require(!controls.select("R","OTHER","2024-07-17T02:00:00"),"unselected receiver produced matrices");
    require(controls.select("R","JPLM","2024-07-17T02:00:00") && controls.select("R","JPLM","2024-07-17T02:00:00"),"multiple domains double-count same epoch");
    require(controls.select("R","JPLM","2024-07-17T02:00:30") && !controls.select("R","JPLM","2024-07-17T02:01:00"),"default two-epoch log limit failed");
    controls.startGpst="2024-07-17T02:10:00";
    require(!controls.select("S","JPLM","2024-07-17T02:09:30"),"replay start gate failed");
    controls.maximumEpochs=0;
    for (int i=0;i<200;++i)
    {
        std::ostringstream epoch;
        const int seconds=2*3600+10*60+30*i;
        epoch<<"2024-07-17T"<<std::setfill('0')<<std::setw(2)<<seconds/3600
             <<':'<<std::setw(2)<<seconds/60%60<<':'<<std::setw(2)<<seconds%60;
        require(controls.select("S","JPLM",epoch.str()),"explicit whole-record replay stopped");
    }
    require(controls.counters.size()==2,"per-epoch replay history grew without bound");
}
}

int main()
{
    try
    {
        testExactCacheKey();
        testRealDuplicateAndTargetedPar();
        testExactCoordinateRecovery();
        testReplayDecimalRoundtripAndSelection();
        std::cout<<"USER_GONSIS_REPLAY_CACHE_CONTRACT PASS checks="<<checks
                 <<" full_numeric_cache_key=1 epoch_local_cache=1 targeted_PAR=1"
                 <<" exact_input_integer_recovery=1 replay_roundtrip_bits=1 default_bounded_selected_logging=1\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr<<"USER_GONSIS_REPLAY_CACHE_CONTRACT FAIL "<<error.what()<<'\n';
        return 1;
    }
}
