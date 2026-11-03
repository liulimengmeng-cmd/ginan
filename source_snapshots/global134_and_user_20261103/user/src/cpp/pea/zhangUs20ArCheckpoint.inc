// Only causal user AR bookkeeping. The previous epoch's disposable observation
// snapshot is never a future factor: time and commit-sequence guards require a
// new authoritative PPP capture before the next epoch's AR can consume it.
#include "common/zhangUs20ArCheckpoint.hpp"
#include <boost/json.hpp>
namespace
{
using ArObject=boost::json::object;
using ArArray=boost::json::array;
string arString(const ArObject& o,const char* k){return string(o.at(k).as_string());}
double arDouble(const ArObject& o,const char* k)
{
    const auto v=boost::json::value_to<double>(o.at(k));
    if(!std::isfinite(v))throw std::runtime_error("US20_AR_NONFINITE");return v;
}
std::uint64_t arCount(const ArObject& o,const char* k){return boost::json::value_to<std::uint64_t>(o.at(k));}
ArObject arPolicy(const ZhangUserArExperimentPolicy& p)
{return {{"seconds",p.minimumArcSeconds},{"samples",p.minimumArcSamples},{"gap",p.maximumGapSeconds},
    {"confirmation_epochs",p.confirmationEpochs},{"confirmation_interval",p.confirmationIntervalSeconds},{"time_tolerance",p.timeToleranceSeconds}};}
ZhangUserArExperimentPolicy arPolicy(const ArObject& o)
{
    if(o.size()!=6)throw std::runtime_error("US20_AR_POLICY_FIELDS");
    ZhangUserArExperimentPolicy p;p.minimumArcSeconds=arDouble(o,"seconds");p.minimumArcSamples=arCount(o,"samples");
    p.maximumGapSeconds=arDouble(o,"gap");p.confirmationEpochs=arCount(o,"confirmation_epochs");
    p.confirmationIntervalSeconds=arDouble(o,"confirmation_interval");p.timeToleranceSeconds=arDouble(o,"time_tolerance");return p;
}
ArObject arTracker(const ZhangUserArcTrackerCheckpoint& value)
{
    ArArray records;
    for(const auto& [k,a]:value.arcs)records.push_back(ArObject{{"runtime",k.runtime},{"receiver",k.receiver},{"object",k.object},
        {"valid",a.valid},{"eligible",a.eligible},{"first",a.first},{"last",a.last},{"age",a.age},{"count",a.count},
        {"identity",a.identity},{"reason",a.reason}});
    return {{"policy",arPolicy(value.policy)},{"records",records}};
}
ZhangUserArcTrackerCheckpoint arTracker(const ArObject& o)
{
    if(o.size()!=2 || o.at("records").as_array().size()>100000)throw std::runtime_error("US20_AR_TRACKER_FIELDS");
    ZhangUserArcTrackerCheckpoint value;value.policy=arPolicy(o.at("policy").as_object());
    for(const auto& record:o.at("records").as_array())
    {
        const auto& v=record.as_object();if(v.size()!=11)throw std::runtime_error("US20_AR_RECORD_FIELDS");
        ZhangUserArExperimentKey k{arString(v,"runtime"),arString(v,"receiver"),arString(v,"object")};
        ZhangUserArExperimentArc a;a.valid=v.at("valid").as_bool();a.eligible=v.at("eligible").as_bool();
        a.first=arDouble(v,"first");a.last=arDouble(v,"last");a.age=arDouble(v,"age");a.count=arCount(v,"count");
        a.identity=arString(v,"identity");a.reason=arString(v,"reason");
        if(!value.arcs.emplace(k,a).second)throw std::runtime_error("US20_AR_DUPLICATE_RECORD");
    }
    return value;
}
string arTime(const GTime& t)
{
    if(t==GTime::noTime())return "NONE";
    std::ostringstream out;out<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<t.bigTime;return out.str();
}
GTime arTime(const string& text,double completed)
{
    if(text=="NONE")return GTime::noTime();
    std::size_t used=0;GTime t;t.bigTime=std::stold(text,&used);
    if(used!=text.size() || !std::isfinite(t.bigTime) || t==GTime::noTime() || t.bigTime>completed || t.bigTime<0)
        throw std::runtime_error("US20_AR_TIME_INVALID");return t;
}
ArObject arKey(const ZhangUserGonsisReferenceKey& key)
{return {{"runtime",std::get<0>(key)},{"receiver",std::get<1>(key)},{"system",int(std::get<2>(key))},{"domain",std::get<3>(key)}};}
ZhangUserGonsisReferenceKey arKey(const ArObject& o,const string& runtime)
{
    const auto r=arString(o,"runtime"),receiver=arString(o,"receiver"),domain=arString(o,"domain");
    const auto sys=boost::json::value_to<int>(o.at("system"));
    if(!validZhangAmbresRuntimeId(r) || !zhangAmbresRuntimeBelongsTo(r,runtime) || receiver.empty() || domain.empty() ||
        (sys!=int(E_Sys::GPS) && sys!=int(E_Sys::GAL) && sys!=int(E_Sys::BDS)))throw std::runtime_error("US20_AR_GROUP_IDENTITY");
    return {r,receiver,static_cast<E_Sys>(sys),domain};
}
void arSatellite(const string& name,E_Sys sys)
{const SatSys s(name.c_str());if(s.id()!=name || s.prn<=0 || s.sys!=sys)throw std::runtime_error("US20_AR_SATELLITE_IDENTITY");}
}

static std::string encodeZhangUs20UserArCheckpoint(const std::string& runtime,double completed,
    const ZhangUserArExperimentArcTracker& legacy,
    const map<pair<double,std::size_t>,ZhangUserArSeparatedArcTracker>& separatedTrackers,
    const decltype(zhangUserGonsisReferences)& referenceRecords,
    const decltype(zhangUserGonsisFixedLocks)& lockRecords)
{
    if(!validZhangAmbresRuntimeId(runtime) || !std::isfinite(completed))throw std::runtime_error("US20_AR_EXPORT_IDENTITY");
    ArObject doc{{"schema","US20_USER_AR_CAUSAL_V1"},{"runtime",runtime},{"completed",completed},
        {"epoch_observations","REPLACED_BEFORE_NEXT_AR_AND_TIME_SEQUENCE_GUARDED"},
        {"legacy_arcs",arTracker(legacy.checkpoint(runtime))}};
    ArArray separated,refs,locks;
    for(const auto& [policy,tracker]:separatedTrackers)
    {
        const auto values=tracker.checkpoint(runtime);
        if(std::all_of(values.begin(),values.end(),[](const auto& x){return x.arcs.empty();}))continue;
        ArArray four;for(const auto& x:values)four.push_back(arTracker(x));
        separated.push_back(ArObject{{"minimum_seconds",policy.first},{"minimum_samples",policy.second},{"trackers",four}});
    }
    for(const auto& [key,value]:referenceRecords)if(zhangAmbresRuntimeBelongsTo(std::get<0>(key),runtime))
    {
        auto o=arKey(key);o["time"]=arTime(value.time);o["satellite"]=value.hint.satellite;
        o["identity"]=value.hint.identity;o["previous_published_fix"]=value.hint.previousPublishedFix;refs.push_back(o);
    }
    for(const auto& [key,sats]:lockRecords)if(zhangAmbresRuntimeBelongsTo(std::get<0>(key),runtime))
    {
        auto o=arKey(key);ArArray records;
        for(const auto& [sat,v]:sats)records.push_back(ArObject{{"satellite",sat},{"identity",v.identity},{"count",v.count},{"last",arTime(v.lastIncrement)}});
        o["records"]=records;locks.push_back(o);
    }
    doc["separated_arcs"]=separated;doc["references"]=refs;doc["published_fix_locks"]=locks;
    return boost::json::serialize(doc);
}

std::string exportZhangUs20UserArCheckpoint(const std::string& runtime,double completed)
{
    return encodeZhangUs20UserArCheckpoint(runtime,completed,zhangUserExperimentArcs(),
        zhangUserSeparatedArcTrackers(),zhangUserGonsisReferences,zhangUserGonsisFixedLocks);
}

void importZhangUs20UserArCheckpoint(const std::string& payload,const std::string& runtime,double completed,bool dryRun)
{
    if(payload.empty() || payload.size()>64*1024*1024 || !validZhangAmbresRuntimeId(runtime) || !std::isfinite(completed))
        throw std::runtime_error("US20_AR_IMPORT_BOUNDARY");
    const auto doc=boost::json::parse(payload).as_object();
    if(doc.size()!=8 || arString(doc,"schema")!="US20_USER_AR_CAUSAL_V1" || arString(doc,"runtime")!=runtime ||
        arDouble(doc,"completed")!=completed || arString(doc,"epoch_observations")!="REPLACED_BEFORE_NEXT_AR_AND_TIME_SEQUENCE_GUARDED")
        throw std::runtime_error("US20_AR_IMPORT_IDENTITY");
    auto legacy=zhangUserExperimentArcs();legacy.restoreCheckpoint(arTracker(doc.at("legacy_arcs").as_object()),runtime,completed);
    auto separated=zhangUserSeparatedArcTrackers();
    for(auto& [key,tracker]:separated)
    {
        auto empty=tracker.checkpoint(runtime);for(auto& e:empty)e.arcs.clear();tracker.restoreCheckpoint(empty,runtime,completed);
    }
    const auto reliability=zhangUserArReliabilityPolicy();
    set<pair<double,std::size_t>> seen;
    for(const auto& entry:doc.at("separated_arcs").as_array())
    {
        const auto& o=entry.as_object();const auto seconds=arDouble(o,"minimum_seconds");const auto samples=arCount(o,"minimum_samples");
        if(o.size()!=3 || seconds!=reliability.minimumArcSeconds || samples!=reliability.minimumArcSamples ||
            !seen.emplace(seconds,samples).second || o.at("trackers").as_array().size()!=4)throw std::runtime_error("US20_AR_SEPARATED_POLICY");
        ZhangUserArExperimentPolicy p;p.minimumArcSeconds=seconds;p.minimumArcSamples=samples;
        auto& tracker=separated.try_emplace(std::make_pair(seconds,samples),p).first->second;
        std::array<ZhangUserArcTrackerCheckpoint,4> values;
        for(int i=0;i<4;++i)values[i]=arTracker(o.at("trackers").as_array()[i].as_object());
        tracker.restoreCheckpoint(values,runtime,completed);
    }
    auto refs=zhangUserGonsisReferences;auto locks=zhangUserGonsisFixedLocks;
    auto clearFamily=[&](auto& records){for(auto it=records.begin();it!=records.end();)
        if(zhangAmbresRuntimeBelongsTo(std::get<0>(it->first),runtime))it=records.erase(it);else ++it;};
    clearFamily(refs);clearFamily(locks);
    for(const auto& entry:doc.at("references").as_array())
    {
        const auto& o=entry.as_object();if(o.size()!=8)throw std::runtime_error("US20_AR_REFERENCE_FIELDS");const auto key=arKey(o,runtime);
        ZhangUserGonsisPublishedReference v;v.time=arTime(arString(o,"time"),completed);v.hint.satellite=arString(o,"satellite");
        v.hint.identity=arString(o,"identity");v.hint.previousPublishedFix=o.at("previous_published_fix").as_bool();
        arSatellite(v.hint.satellite,std::get<2>(key));
        if(v.time==GTime::noTime() || v.hint.identity.empty() || !refs.emplace(key,v).second)throw std::runtime_error("US20_AR_REFERENCE_INVALID");
    }
    for(const auto& entry:doc.at("published_fix_locks").as_array())
    {
        const auto& o=entry.as_object();if(o.size()!=5)throw std::runtime_error("US20_AR_LOCK_FIELDS");const auto key=arKey(o,runtime);
        map<string,ZhangUserGonsisFixedLock> values;
        for(const auto& item:o.at("records").as_array())
        {
            const auto& v=item.as_object();if(v.size()!=4)throw std::runtime_error("US20_AR_LOCK_RECORD_FIELDS");
            const auto sat=arString(v,"satellite");arSatellite(sat,std::get<2>(key));
            ZhangUserGonsisFixedLock lock;lock.identity=arString(v,"identity");lock.count=boost::json::value_to<int>(v.at("count"));
            lock.lastIncrement=arTime(arString(v,"last"),completed);
            if(lock.identity.empty() || lock.count<0 || (lock.count==0)!=(lock.lastIncrement==GTime::noTime()) ||
                !values.emplace(sat,lock).second)throw std::runtime_error("US20_AR_LOCK_RECORD_INVALID");
        }
        if(!locks.emplace(key,std::move(values)).second)throw std::runtime_error("US20_AR_DUPLICATE_LOCK_GROUP");
    }
    // Check exact canonical representation on the private candidate, before
    // publishing any member of the transaction. A different numeric spelling
    // or field order must not throw only after updating live FIX counters.
    if(encodeZhangUs20UserArCheckpoint(runtime,completed,legacy,separated,refs,locks)!=payload)
        throw std::runtime_error("US20_AR_NONCANONICAL_PAYLOAD");
    if(dryRun)return;
    zhangUserExperimentArcs()=std::move(legacy);zhangUserSeparatedArcTrackers()=std::move(separated);
    zhangUserGonsisReferences=std::move(refs);zhangUserGonsisFixedLocks=std::move(locks);
    for(auto it=zhangUserGonsisObservationSnapshots.begin();it!=zhangUserGonsisObservationSnapshots.end();)
        if(zhangAmbresRuntimeBelongsTo(it->first,runtime))it=zhangUserGonsisObservationSnapshots.erase(it);else ++it;
}
