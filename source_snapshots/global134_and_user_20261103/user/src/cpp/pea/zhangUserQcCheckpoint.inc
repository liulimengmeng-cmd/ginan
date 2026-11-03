// US20-only supplement. The legacy receiver envelope and its bytes are unchanged.
namespace
{
struct Us20QcSample
{
    double time=0,value=0,sigma=0;
    template<class A> void serialize(A& a,const unsigned int) {a & time & value & sigma;}
};
struct Us20QcTrend
{
    std::vector<Us20QcSample> accepted;
    unsigned suspects=0;
    bool rawSeen=false;
    double lastRawTime=0,firstThreshold=0;
    template<class A> void serialize(A& a,const unsigned int)
    {a & accepted & suspects & rawSeen & lastRawTime & firstThreshold;}
};
struct Us20QcSignal
{
    bool physicalBreak=false;
    int state=0;
    template<class A> void serialize(A& a,const unsigned int) {a & physicalBreak & state;}
};
struct Us20QcSatellite
{
    Us20QcTrend gf,mw,codeGf;
    bool accepted=false,warmup=false,gap=false,available=false;
    double prediction=0,innovation=0,threshold=0;
    std::string reason;
    int frequencyA=0,frequencyB=0;
    std::map<std::string,Us20QcSignal> signals;
    template<class A> void serialize(A& a,const unsigned int)
    {a & gf & mw & codeGf & accepted & warmup & gap & available & prediction & innovation & threshold
       & reason & frequencyA & frequencyB & signals;}
};
struct Us20QcEnvelope
{
    unsigned version=1;
    std::string runtime;
    double completed=0;
    std::map<std::string,std::map<SatSys,Us20QcSatellite>> receivers;
    template<class A> void serialize(A& a,const unsigned int) {a & version & runtime & completed & receivers;}
};
const char* us20QcReason(const std::string& reason)
{
    // Stable literals: never retain a pointer into a temporary decoded string.
    for(const char* value:{"INVALID_SAMPLE","TREND_WARMUP","STALE_ACCEPTED_HISTORY_UNRESOLVED",
        "PERSISTENT_UNRESOLVED_INNOVATION","UNRESOLVED_INNOVATION","TREND_ACCEPTED"})
        if(reason==value)return value;
    throw std::runtime_error("US20_QC_UNKNOWN_DECISION_REASON:"+reason);
}
Us20QcTrend us20QcCaptureTrend(const ZhangUserQcTrend& input)
{
    Us20QcTrend result;
    for(const auto& s:input.accepted)result.accepted.push_back({s.time,s.value,s.sigma});
    result.suspects=input.consecutiveSuspects;result.rawSeen=input.rawSeen;
    result.lastRawTime=input.lastRawTime;result.firstThreshold=input.firstSuspectThreshold;
    return result;
}
ZhangUserQcTrend us20QcRestoreTrend(const Us20QcTrend& input,double completed)
{
    if(input.accepted.size()>10000 || !std::isfinite(input.lastRawTime) || !std::isfinite(input.firstThreshold) ||
       input.firstThreshold<0 || (input.rawSeen && input.lastRawTime>completed) ||
       (!input.rawSeen && (!input.accepted.empty() || input.suspects)))
        throw std::runtime_error("US20_QC_TREND_HEADER_INVALID");
    ZhangUserQcTrend result;
    double previous=-std::numeric_limits<double>::infinity();
    for(const auto& s:input.accepted)
    {
        if(!std::isfinite(s.time) || !std::isfinite(s.value) || !std::isfinite(s.sigma) || s.sigma<=0 ||
           s.time<=previous || s.time>completed || s.time>input.lastRawTime)
            throw std::runtime_error("US20_QC_TREND_SAMPLE_INVALID");
        result.accepted.push_back({s.time,s.value,s.sigma});previous=s.time;
    }
    result.consecutiveSuspects=input.suspects;result.rawSeen=input.rawSeen;
    result.lastRawTime=input.lastRawTime;result.firstSuspectThreshold=input.firstThreshold;
    return result;
}
Us20QcEnvelope us20QcDecode(const std::string& payload,const std::string& runtime,double completed)
{
    Us20QcEnvelope value;std::string reason;
    if(!deserializePayload(payload,value,reason,"US20_QC"))throw std::runtime_error(reason);
    if(value.version!=1 || value.runtime!=runtime || runtime.empty() || value.completed!=completed ||
       !std::isfinite(completed) || value.receivers.empty())throw std::runtime_error("US20_QC_IDENTITY_MISMATCH");
    for(const auto& [receiver,satellites]:value.receivers)
        for(const auto& [satellite,s]:satellites)
        {
            (void)us20QcRestoreTrend(s.gf,completed);(void)us20QcRestoreTrend(s.mw,completed);
            (void)us20QcRestoreTrend(s.codeGf,completed);(void)us20QcReason(s.reason);
            if(!std::isfinite(s.prediction) || !std::isfinite(s.innovation) || !std::isfinite(s.threshold) ||
               s.threshold<0 || s.frequencyA<0 || s.frequencyB<0)throw std::runtime_error("US20_QC_DECISION_INVALID");
            for(const auto& [signal,q]:s.signals)
                if(signal.empty() || q.state<0 || q.state>2)throw std::runtime_error("US20_QC_SIGNAL_INVALID");
        }
    return value;
}
}

std::string exportZhangUs20QcCheckpoint(const ReceiverMap& receivers,const std::string& runtime,double completed)
{
    Us20QcEnvelope value;value.runtime=runtime;value.completed=completed;
    for(const auto& [id,receiver]:receivers)
    {
        auto& satellites=value.receivers[id];
        for(const auto& [satellite,status]:receiver.satStatMap)
        {
            auto& s=satellites[satellite];
            s.gf=us20QcCaptureTrend(status.userAcceptedGf);s.mw=us20QcCaptureTrend(status.userAcceptedMw);
            s.codeGf=us20QcCaptureTrend(status.userAcceptedCodeGf);
            const auto& d=status.userGfDecision;s.accepted=d.accepted;s.warmup=d.warmup;s.gap=d.gap;
            s.prediction=d.prediction;s.innovation=d.innovation;s.threshold=d.threshold;s.reason=d.reason;
            s.available=status.userGfAvailable;s.frequencyA=status.userQcFrequencyA;s.frequencyB=status.userQcFrequencyB;
            for(const auto& [signal,q]:status.sigStatMap)s.signals[signal]={q.userPhysicalBreakThisEpoch,static_cast<int>(q.userQcState)};
        }
    }
    std::string payload,reason;
    if(!serializePayload(value,payload,reason,"US20_QC"))throw std::runtime_error(reason);
    (void)us20QcDecode(payload,runtime,completed);return payload;
}

void validateZhangUs20QcCheckpoint(const std::string& payload,const std::string& legacyPayload,
    const std::string& runtime,double completed)
{
    const auto value=us20QcDecode(payload,runtime,completed);
    ZhangReceiverRuntimeEnvelope legacy;ZhangReceiverRuntimeCheckpointResult result;
    if(!decodeAndValidateReceiverEnvelope(legacyPayload,runtime,legacy,result))throw std::runtime_error(result.failureReason);
    if(value.receivers.size()!=legacy.receivers.size())throw std::runtime_error("US20_QC_RECEIVER_COVERAGE");
    for(const auto& receiver:legacy.receivers)
    {
        const auto& satellites=value.receivers.at(receiver.mapKey);
        if(satellites.size()!=receiver.satelliteStatuses.size())throw std::runtime_error("US20_QC_SATELLITE_COVERAGE");
        for(const auto& [satellite,status]:receiver.satelliteStatuses)
        {
            const auto& signals=satellites.at(satellite).signals;
            if(signals.size()!=status.signalStatuses.size())throw std::runtime_error("US20_QC_SIGNAL_COVERAGE");
            for(const auto& [signal,ignored]:status.signalStatuses)(void)signals.at(signal);
        }
    }
}

void importZhangUs20QcCheckpoint(const std::string& payload,ReceiverMap& receivers,
    const std::string& runtime,double completed)
{
    const auto value=us20QcDecode(payload,runtime,completed);
    // Check every destination before writing any supplemental field.
    if(value.receivers.size()!=receivers.size())throw std::runtime_error("US20_QC_TARGET_RECEIVER_COVERAGE");
    for(const auto& [id,satellites]:value.receivers)
    {
        const auto& receiver=receivers.at(id);
        if(satellites.size()!=receiver.satStatMap.size())throw std::runtime_error("US20_QC_TARGET_SATELLITE_COVERAGE");
        for(const auto& [satellite,s]:satellites)
        {
            const auto& signals=receiver.satStatMap.at(satellite).sigStatMap;
            if(s.signals.size()!=signals.size())throw std::runtime_error("US20_QC_TARGET_SIGNAL_COVERAGE");
            for(const auto& [signal,q]:s.signals)(void)signals.at(signal);
        }
    }
    for(const auto& [id,satellites]:value.receivers)
        for(const auto& [satellite,s]:satellites)
        {
            auto& status=receivers.at(id).satStatMap.at(satellite);
            status.userAcceptedGf=us20QcRestoreTrend(s.gf,completed);
            status.userAcceptedMw=us20QcRestoreTrend(s.mw,completed);
            status.userAcceptedCodeGf=us20QcRestoreTrend(s.codeGf,completed);
            auto& d=status.userGfDecision;d.accepted=s.accepted;d.warmup=s.warmup;d.gap=s.gap;
            d.prediction=s.prediction;d.innovation=s.innovation;d.threshold=s.threshold;d.reason=us20QcReason(s.reason);
            status.userGfAvailable=s.available;status.userQcFrequencyA=s.frequencyA;status.userQcFrequencyB=s.frequencyB;
            for(const auto& [signal,q]:s.signals)
            {
                auto& target=status.sigStatMap.at(signal);target.userPhysicalBreakThisEpoch=q.physicalBreak;
                target.userQcState=static_cast<ZhangUserQcState>(q.state);
            }
        }
    if(exportZhangUs20QcCheckpoint(receivers,runtime,completed)!=payload)
        throw std::runtime_error("US20_QC_RESTORE_NOT_EXACT");
}
