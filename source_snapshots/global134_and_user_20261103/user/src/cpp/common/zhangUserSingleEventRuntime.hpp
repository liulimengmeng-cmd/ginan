#pragma once

#include "common/zhangUserSingleEvent.hpp"
#include "common/zhangUserResetCheckpoint.hpp"
#include "common/zhangUserDatumPolicy.hpp"
#include <optional>

// Predeclared, same-plan runtime integration. A later post-selection fork
// contract must still prove its distinct parent/execution identities; this
// path does not permit a different plan, binary, input set or scientific YAML.
namespace zhangSingleEventRuntime
{
inline const char* sectionName() { return zhangUs20Environment("ZHANG_US20_EVENT_PROTOCOL")=="2" ? "us20.event_fork_control.v2" : "us20.event_fork_control.v1"; }
using namespace zhangSingleEvent;
inline bool declared()
{
    const auto value=zhangUs20Environment("ZHANG_US20_EVENT_PROTOCOL");
    require(value.empty() || value=="0" || value=="1" || value=="2","SINGLE_EVENT_PROTOCOL_OPTION_INVALID");
    return value=="1" || value=="2";
}

inline std::string string(const boost::json::object& v, const char* name)
{ return std::string(v.at(name).as_string()); }
inline long long integer(const boost::json::object& v, const char* name)
{ return v.at(name).as_int64(); }

inline std::string actionName(Action action)
{ return action == Action::KEEP ? "KEEP" : "RESET_SELECTED_EVENT"; }
inline Action action(const std::string& value)
{
    require(value == "KEEP" || value == "RESET_SELECTED_EVENT", "SINGLE_EVENT_ACTION_INVALID");
    return value == "KEEP" ? Action::KEEP : Action::RESET_SELECTED_EVENT;
}
inline std::string statusName(Status value)
{
    switch(value)
    {
        case Status::PENDING: return "PENDING";
        case Status::ELIGIBLE: return "ELIGIBLE";
        case Status::INELIGIBLE: return "INELIGIBLE";
        case Status::MISSED: return "MISSED";
        default: throw std::runtime_error("SINGLE_EVENT_NONPERSISTENT_STATUS");
    }
}
inline Status status(const std::string& value)
{
    if(value == "PENDING") return Status::PENDING;
    if(value == "ELIGIBLE") return Status::ELIGIBLE;
    if(value == "INELIGIBLE") return Status::INELIGIBLE;
    if(value == "MISSED") return Status::MISSED;
    throw std::runtime_error("SINGLE_EVENT_STATUS_INVALID");
}

struct Context
{
    Plan plan;
    Action action = Action::KEEP;
    std::string group, parentPlanSha, originBranch, originEpoch, commitPath;
    std::string observedBranch;
    boost::json::object parentPlan;
    bool registered=false;
    std::string candidateId,externalManifestSha,expectedOriginSha,expectedOriginPath,registeredCommonPlanSha;
    long long stopEpoch=0;
};

inline bool registeredMode();
inline Context registeredContext(const boost::json::object&,const std::string&,const std::string&);

inline bool branchEnabled(const boost::json::object& plan, const std::string& name)
{
    if(!plan.contains("single_event") && !plan.contains("single_events")) return false;
    const auto& branch = plan.at("branches").as_object().at(name).as_object();
    return branch.contains("single_event_action") && string(branch,"single_event_action") != "DISABLED";
}

inline const Context* context()
{
    if(!declared() || !zhangUs20CheckpointEnabled() || !zhangUs20Environment("ZHANG_US20_DUMP_PROVENANCE").empty()) return nullptr;
    static const std::optional<Context> configured = []() -> std::optional<Context>
    {
        const auto [text, plan] = zhangUs20ForkPlan();
        const auto name = zhangUs20Environment("ZHANG_US20_FORK_BRANCH");
        if(!branchEnabled(plan,name)) return std::nullopt;
        require(zhangUserDatumPolicyRequested() && zhangUserDatumPolicy() == ZhangUserDatumPolicy::TRANSPORT_KNOWN,
            "SINGLE_EVENT_REQUIRES_GLOBAL_TRANSPORT_BASE");
        if(registeredMode())return registeredContext(plan,zhangCheckpointSha256(text),name);
        Context out;
        out.parentPlan = plan; out.parentPlanSha = zhangCheckpointSha256(text); out.observedBranch = name;
        const auto& branch = plan.at("branches").as_object().at(name).as_object();
        out.action = action(string(branch,"single_event_action"));
        out.group = string(branch,"single_event_group");
        require(!out.group.empty(), "SINGLE_EVENT_EXECUTION_GROUP_MISSING");
        const auto& spec = plan.at("single_event").as_object();
        require(string(spec,"schema") == "US20_PREDECLARED_SINGLE_EVENT_V1", "SINGLE_EVENT_PLAN_SCHEMA");
        out.originBranch = string(spec,"origin_branch"); out.originEpoch = string(spec,"origin_epoch");
        require(!branchEnabled(plan,out.originBranch), "SINGLE_EVENT_ORIGIN_ALREADY_HAS_CONTROL");
        out.commitPath = string(spec,"available_commit_path");
        out.plan.receiver = string(spec,"receiver"); out.plan.manifestSha = zhangCheckpointSha256(boost::json::serialize(spec));
        out.plan.experimentPlanSha = out.parentPlanSha; out.plan.streamSha = string(spec,"product_stream_sha256");
        out.plan.availableCommitSha = string(spec,"available_commit_sha256");
        out.plan.eventEpoch = integer(spec,"event_gpst_seconds");
        out.plan.availableEpoch = integer(spec,"available_gpst_seconds");
        out.plan.applyEpoch = integer(spec,"apply_gpst_seconds");
        for(const auto& value : spec.at("certificates").as_array()) out.plan.certificates.push_back(std::string(value.as_string()));
        // Parent SHA is established by the actual verified parent bundle, not
        // inserted into the common plan before that bundle exists.
        out.plan.parentCheckpointSha = std::string(64,'0'); validate(out.plan);
        require(zhangCheckpointFileSha256(zhangUs20Environment("ZHANG_RESET_PRODUCT_MANIFEST")) == out.plan.streamSha,
            "SINGLE_EVENT_COMPLETE_PRODUCT_STREAM_MISMATCH");
        return out;
    }();
    return configured ? &*configured : nullptr;
}

struct State
{
    std::string originPath, originSha;
    std::optional<Receipt> receipt;
};
inline std::optional<State> committed;

inline boost::json::object receiptJson(const Receipt& r)
{
    boost::json::array targets, owners;
    for(const auto& [sat,code] : r.targets) targets.push_back(boost::json::array{sat,code});
    for(const auto& owner : r.stateOwners) owners.push_back(boost::json::value(owner));
    return {{"manifest_sha256",r.manifestSha},{"experiment_plan_sha256",r.experimentPlanSha},
        {"parent_checkpoint_sha256",r.parentCheckpointSha},{"common_plan_sha256",r.commonPlanSha},
        {"common_state_sha256",r.commonStateSha},{"committed_state_sha256",r.committedStateSha},
        {"epoch",r.epoch},{"status",statusName(r.status)},{"action",actionName(r.action)},
        {"targets",targets},{"owners",owners},{"released_rank",static_cast<long long>(r.releasedRank)},{"reason",r.reason}};
}
inline Receipt receiptFromJson(const boost::json::object& v)
{
    Receipt r;
    r.manifestSha=string(v,"manifest_sha256");r.experimentPlanSha=string(v,"experiment_plan_sha256");
    r.parentCheckpointSha=string(v,"parent_checkpoint_sha256");r.commonPlanSha=string(v,"common_plan_sha256");
    r.commonStateSha=string(v,"common_state_sha256");r.committedStateSha=string(v,"committed_state_sha256");
    r.epoch=integer(v,"epoch");r.status=status(string(v,"status"));r.action=action(string(v,"action"));
    require(integer(v,"released_rank")>=0,"SINGLE_EVENT_NEGATIVE_RELEASE_RANK");r.releasedRank=integer(v,"released_rank");
    r.reason=string(v,"reason");
    for(const auto& item:v.at("targets").as_array())
    {
        const auto& a=item.as_array();require(a.size()==2,"SINGLE_EVENT_TARGET_TUPLE_SIZE");
        r.targets.emplace_back(std::string(a[0].as_string()),std::string(a[1].as_string()));
    }
    for(const auto& item:v.at("owners").as_array())r.stateOwners.push_back(std::string(item.as_string()));
    require(receiptJson(r)==v,"SINGLE_EVENT_RECEIPT_SCHEMA_NOT_EXACT");return r;
}
inline Plan boundPlan(const Context& c,const State& s)
{ auto p=c.plan;p.parentCheckpointSha=s.originSha;validate(p);return p; }

inline void validateOrigin(const Context& c,const State& s)
{
    require(sha256(s.originSha) && zhangCheckpointFileSha256((std::filesystem::path(s.originPath)/"checkpoint.bundle").string())==s.originSha,
        "SINGLE_EVENT_ORIGINAL_PARENT_CHANGED");
    if(c.registered && !c.expectedOriginSha.empty())
        require(s.originSha==c.expectedOriginSha && std::filesystem::weakly_canonical(s.originPath).string()==c.expectedOriginPath,
            "REGISTERED_FORK_ORIGINAL_PARENT_NOT_SELECTED");
    const auto manifest=boost::json::parse(zhangUs20ReadFile(std::filesystem::path(s.originPath)/"checkpoint_manifest.json")).as_object();
    require(string(manifest,"epoch")==c.originEpoch && string(manifest,"binary_sha256")==string(c.parentPlan,"binary_sha256"),
        "SINGLE_EVENT_ORIGIN_EPOCH_OR_BINARY_MISMATCH");
}

inline std::string encode(const Context& c,const State& s)
{
    boost::json::object value{{"schema","US20_EVENT_FORK_CONTROL_V1"},{"parent_plan_sha256",c.parentPlanSha},
        {"selector_sha256",c.plan.manifestSha},{"group",c.group},{"action",actionName(c.action)},
        {"origin_path",s.originPath},{"origin_sha256",s.originSha},{"receipt",s.receipt?boost::json::value(receiptJson(*s.receipt)):boost::json::value(nullptr)}};
    if(c.registered)
    {
        value["schema"]="US20_EVENT_FORK_CONTROL_V2";value["candidate_id"]=c.candidateId;
        value["event_fork_manifest_sha256"]=c.externalManifestSha;value["stop_gpst_seconds"]=c.stopEpoch;
    }
    return boost::json::serialize(value);
}
inline State decode(const Context& c,const std::string& payload,long long completed)
{
    const auto v=boost::json::parse(payload).as_object();
    require(string(v,"schema")== (c.registered?"US20_EVENT_FORK_CONTROL_V2":"US20_EVENT_FORK_CONTROL_V1") && string(v,"parent_plan_sha256")==c.parentPlanSha &&
        string(v,"selector_sha256")==c.plan.manifestSha && string(v,"group")==c.group && string(v,"action")==actionName(c.action),
        "SINGLE_EVENT_CONTROL_BINDING_MISMATCH");
    State s{string(v,"origin_path"),string(v,"origin_sha256"),{}};
    if(!v.at("receipt").is_null())s.receipt=receiptFromJson(v.at("receipt").as_object());
    require(encode(c,s)==payload,"SINGLE_EVENT_CONTROL_NONCANONICAL_OR_EXTRA_FIELDS");
    validateOrigin(c,s);
    if(s.receipt)
    {
        require(s.receipt->action==c.action,"SINGLE_EVENT_CONTROL_RECEIPT_ACTION_MISMATCH");
        CommonEvidence evidence;evidence.epoch=completed;
        require(evaluate(boundPlan(c,s),evidence,&*s.receipt).status==Status::ALREADY_TERMINAL,
            "SINGLE_EVENT_CONTROL_TERMINAL_RECEIPT_INVALID");
    }
    else require(completed<c.plan.applyEpoch,"SINGLE_EVENT_PENDING_AFTER_APPLICATION_BOUNDARY");
    return s;
}

inline bool sectionExpected(const std::string& ledgerPayload)
{
    if(!declared())return false;
    const auto [text,plan]=zhangUs20ForkPlan();
    const auto ledger=boost::json::parse(ledgerPayload).as_object();
    require(string(ledger,"US20_RESET_CHECKPOINT_PLAN_SHA256")==zhangCheckpointSha256(text),
        "SINGLE_EVENT_SOURCE_LEDGER_PLAN_MISMATCH");
    return branchEnabled(plan,string(ledger,"US20_RESET_CHECKPOINT_SOURCE_BRANCH"));
}

inline void restore(const std::string* payload,const std::string& ledgerPayload,long long completed,bool dryRun)
{
    const bool expected=sectionExpected(ledgerPayload);const auto* c=context();
    if(!c)
    { require(!expected && payload==nullptr,"SINGLE_EVENT_CONTROL_CANNOT_RESTORE_AS_ORDINARY_PARENT");return; }
    require(expected==(payload!=nullptr),"SINGLE_EVENT_MISSING_OR_UNEXPECTED_CONTROL_SECTION");
    State next;
    if(payload)next=decode(*c,*payload,completed);
    else
    {
        const auto ledger=boost::json::parse(ledgerPayload).as_object();
        require(string(ledger,"US20_RESET_CHECKPOINT_SOURCE_BRANCH")==c->originBranch && completed==c->plan.applyEpoch-30,
            "SINGLE_EVENT_INITIAL_PARENT_NOT_IMMEDIATELY_BEFORE_EVENT");
        next.originPath=std::filesystem::canonical(zhangUs20Environment("ZHANG_E29_RESTORE_BUNDLE")).string();
        next.originSha=zhangCheckpointFileSha256((std::filesystem::path(next.originPath)/"checkpoint.bundle").string());
        validateOrigin(*c,next);boundPlan(*c,next);
    }
    if(!dryRun)committed=std::move(next);
}

inline std::string capture(long long completed)
{
    const auto* c=context();require(c && committed.has_value(),"SINGLE_EVENT_CONTROL_NOT_INITIALIZED");
    const auto text=encode(*c,*committed);decode(*c,text,completed);return text;
}
inline bool collecting(long long epoch)
{
    const auto* c=context();if(!c)return false;
    require(committed.has_value(),"SINGLE_EVENT_MUST_RESTORE_EXACT_PARENT_FIRST");
    return !committed->receipt && epoch==c->plan.applyEpoch;
}
inline std::string stateSha(const KFState& state)
{ return zhangCheckpointSha256(serializeZhangCheckpointSectionPayload(captureZhangCheckpointKfCore(state))); }

// Verify the consumed flat files and the immutable publication containing the
// selected event, rather than copying an alleged COMMIT identity from the plan.
inline void verifyInputStream(const Context& c, const std::string& products, const std::string& covariance,
    const std::string& epoch)
{
    const auto path=zhangUs20Environment("ZHANG_RESET_PRODUCT_MANIFEST");
    const auto text=zhangUs20ReadFile(path);
    require(zhangCheckpointSha256(text)==c.plan.streamSha,"SINGLE_EVENT_STREAM_CHANGED");
    const auto manifest=boost::json::parse(text).as_object();
    require(manifest.at("passed").as_bool(),"SINGLE_EVENT_STREAM_NOT_ACCEPTED");
    const auto& files=manifest.at("files").as_object();
    const auto root=std::filesystem::path(path).parent_path();
    static std::set<std::string> verifiedStreams;
    const auto verificationKey=c.plan.streamSha+"\n"+products+"\n"+covariance;
    if(!registeredMode() || !verifiedStreams.contains(verificationKey))
    for(const auto& pair:std::vector<std::pair<std::string,std::string>>{{"products.csv",products},{"covariance.csv",covariance}})
    {
        const auto& file=files.at(pair.first).as_object();
        require(std::filesystem::equivalent(root/pair.first,pair.second) &&
            zhangCheckpointFileSha256(pair.second)==string(file,"sha256"),"SINGLE_EVENT_CONSUMED_INPUT_MISMATCH");
    }
    verifiedStreams.insert(verificationKey);
    std::size_t matches=0;
    for(const auto& value:manifest.at("source_bundles").as_array())
    {
        const auto& item=value.as_object();if(string(item,"epoch")!=epoch)continue;
        ++matches;
        const auto actual=std::filesystem::path(string(item,"directory"))/"COMMIT.json";
        require(std::filesystem::equivalent(actual,c.commitPath) && string(item,"commit_sha256")==c.plan.availableCommitSha &&
            zhangCheckpointFileSha256(actual.string())==c.plan.availableCommitSha,"SINGLE_EVENT_ACTUAL_PUBLICATION_MISMATCH");
        const auto commit=boost::json::parse(zhangUs20ReadFile(actual)).as_object();
        require(string(commit,"epoch")==epoch,"SINGLE_EVENT_PUBLICATION_EPOCH_MISMATCH");
    }
    require(matches==1,"SINGLE_EVENT_PUBLICATION_NOT_UNIQUE");
}
inline std::string planSha(const CommonEvidence& e)
{
    std::ostringstream out;
    for(const auto& [key,h]:e.held)
    {
        for(const auto* v:{&key.first,&key.second,&h.sourceEndpoint,&h.targetEndpoint,&h.referenceBefore,&h.referenceAfterCommon,
                &h.stateOwner,&h.physicalGenerationBefore,&h.physicalGenerationAfter,&h.statisticalGenerationBefore,&h.statisticalGenerationAfter})
            out<<std::quoted(*v)<<'\n';
        out<<h.existedBefore<<' '<<h.existsAfterCommon<<' '<<h.authorizedBefore<<' '<<h.authorizedAfterCommon<<'\n';
        for(const auto& v:h.consumedCertificates)out<<std::quoted(v)<<'\n';
    }
    for(const auto& v:e.noncommonCertificates)out<<std::quoted(v)<<'\n';
    return zhangCheckpointSha256(out.str());
}
// Included inside zhangSingleEventRuntime. Registration is experiment control;
// it neither changes product permissions nor defines a new scientific policy.
inline bool registeredMode()
{ return zhangUs20Environment("ZHANG_US20_EVENT_PROTOCOL")=="2"; }

inline std::string canonicalPath(const std::string& p)
{ require(!p.empty(),"REGISTERED_FORK_EMPTY_PATH");return std::filesystem::weakly_canonical(p).string(); }
inline bool containsPath(const std::string& a,const std::string& b)
{
    const auto x=std::filesystem::path(canonicalPath(a)), y=std::filesystem::path(canonicalPath(b));
    auto i=x.begin(),j=y.begin();for(;i!=x.end() && j!=y.end();++i,++j)if(*i!=*j)return false;
    return i==x.end();
}
inline void disjoint(const std::string& a,const std::string& b)
{ require(!containsPath(a,b) && !containsPath(b,a),"REGISTERED_FORK_OUTPUT_PATH_OVERLAP"); }

inline Context registeredSpec(const boost::json::object& p,const std::string& planSha,const std::string& id)
{
    Context out;out.parentPlan=p;out.parentPlanSha=planSha;out.candidateId=id;out.registered=true;
    const auto& s=p.at("single_events").as_object().at(id).as_object();
    require(string(s,"schema")=="US20_REGISTERED_SINGLE_EVENT_V2","REGISTERED_FORK_SPEC_SCHEMA");
    out.originBranch=string(s,"origin_branch");out.originEpoch=string(s,"origin_epoch");
    out.commitPath=string(s,"available_commit_path");out.stopEpoch=integer(s,"stop_gpst_seconds");
    auto& plan=out.plan;plan.receiver=string(s,"receiver");plan.manifestSha=zhangCheckpointSha256(boost::json::serialize(s));
    plan.experimentPlanSha=planSha;plan.streamSha=string(s,"product_stream_sha256");
    plan.availableCommitSha=string(s,"available_commit_sha256");
    plan.eventEpoch=integer(s,"event_gpst_seconds");plan.availableEpoch=integer(s,"available_gpst_seconds");
    plan.applyEpoch=integer(s,"apply_gpst_seconds");plan.parentCheckpointSha=std::string(64,'0');
    for(const auto& v:s.at("certificates").as_array())plan.certificates.emplace_back(v.as_string());
    validate(plan);
    require(out.stopEpoch>=plan.applyEpoch && (out.stopEpoch-plan.applyEpoch)%30==0,
        "REGISTERED_FORK_STOP_OUTSIDE_CALENDAR");
    const auto& origin=p.at("branches").as_object().at(out.originBranch).as_object();
    require(string(origin,"role")=="FULL_DAY_T" && string(origin,"single_event_action")=="DISABLED" &&
        string(origin,"policy")=="TRANSPORT_KNOWN","REGISTERED_FORK_ORIGIN_NOT_FULL_DAY_T");
    require(string(p,"station")==plan.receiver,"REGISTERED_FORK_WRONG_STATION");
    GTime t;t.bigTime=plan.applyEpoch-30;
    require(out.originEpoch==t.to_string(0),"REGISTERED_FORK_WRONG_PARENT_SLOT");
    require(zhangCheckpointFileSha256(zhangUs20Environment("ZHANG_RESET_PRODUCT_MANIFEST"))==plan.streamSha,
        "SINGLE_EVENT_COMPLETE_PRODUCT_STREAM_MISMATCH");
    return out;
}

inline boost::json::object evidenceJson(const CommonEvidence& e)
{
    boost::json::array held,noncommon;
    for(const auto& [k,h]:e.held)
    {
        boost::json::array consumed;for(const auto& v:h.consumedCertificates)consumed.emplace_back(v);
        held.push_back(boost::json::object{{"satellite",k.first},{"observable",k.second},
            {"existed_before",h.existedBefore},{"exists_after_common",h.existsAfterCommon},
            {"authorized_before",h.authorizedBefore},{"authorized_after_common",h.authorizedAfterCommon},
            {"source_endpoint",h.sourceEndpoint},{"target_endpoint",h.targetEndpoint},
            {"reference_before",h.referenceBefore},{"reference_after_common",h.referenceAfterCommon},
            {"state_owner",h.stateOwner},{"physical_generation_before",h.physicalGenerationBefore},
            {"physical_generation_after",h.physicalGenerationAfter},{"statistical_generation_before",h.statisticalGenerationBefore},
            {"statistical_generation_after",h.statisticalGenerationAfter},{"consumed_certificates",consumed}});
    }
    for(const auto& v:e.noncommonCertificates)noncommon.emplace_back(v);
    return {{"receiver",e.receiver},{"stream_sha256",e.streamSha},{"available_commit_sha256",e.availableCommitSha},
        {"parent_checkpoint_sha256",e.parentCheckpointSha},{"processing_stage",e.processingStage},
        {"common_plan_sha256",e.commonPlanSha},{"common_state_sha256",e.commonStateSha},
        {"epoch",e.epoch},{"held",held},{"noncommon_certificates",noncommon}};
}
inline CommonEvidence evidenceFromJson(const boost::json::object& v)
{
    CommonEvidence e;e.receiver=string(v,"receiver");e.streamSha=string(v,"stream_sha256");
    e.availableCommitSha=string(v,"available_commit_sha256");e.parentCheckpointSha=string(v,"parent_checkpoint_sha256");
    e.processingStage=string(v,"processing_stage");e.commonPlanSha=string(v,"common_plan_sha256");
    e.commonStateSha=string(v,"common_state_sha256");e.epoch=integer(v,"epoch");
    for(const auto& item:v.at("held").as_array())
    {
        const auto& h=item.as_object();Key k{string(h,"satellite"),string(h,"observable")};HeldEndpoint b;
        b.existedBefore=h.at("existed_before").as_bool();b.existsAfterCommon=h.at("exists_after_common").as_bool();
        b.authorizedBefore=h.at("authorized_before").as_bool();b.authorizedAfterCommon=h.at("authorized_after_common").as_bool();
        b.sourceEndpoint=string(h,"source_endpoint");b.targetEndpoint=string(h,"target_endpoint");
        b.referenceBefore=string(h,"reference_before");b.referenceAfterCommon=string(h,"reference_after_common");
        b.stateOwner=string(h,"state_owner");b.physicalGenerationBefore=string(h,"physical_generation_before");
        b.physicalGenerationAfter=string(h,"physical_generation_after");b.statisticalGenerationBefore=string(h,"statistical_generation_before");
        b.statisticalGenerationAfter=string(h,"statistical_generation_after");
        for(const auto& x:h.at("consumed_certificates").as_array())b.consumedCertificates.emplace_back(x.as_string());
        require(e.held.emplace(k,std::move(b)).second,"REGISTERED_FORK_DUPLICATE_EVIDENCE_KEY");
    }
    for(const auto& x:v.at("noncommon_certificates").as_array())e.noncommonCertificates.emplace_back(x.as_string());
    require(evidenceJson(e)==v && planSha(e)==e.commonPlanSha,"REGISTERED_FORK_COMMON_PLAN_EVIDENCE_CHANGED");
    return e;
}

inline State registeredOrigin(const Context& c)
{
    const auto& b=c.parentPlan.at("branches").as_object().at(c.originBranch).as_object();
    const auto root=std::filesystem::path(string(b,"checkpoint_root"));
    State s;std::size_t matches=0;
    require(std::filesystem::is_directory(root),"REGISTERED_FORK_PARENT_DIRECTORY_MISSING");
    for(const auto& entry:std::filesystem::directory_iterator(root))
    {
        const auto path=entry.path()/"checkpoint_manifest.json";if(!std::filesystem::is_regular_file(path))continue;
        const auto m=boost::json::parse(zhangUs20ReadFile(path)).as_object();if(string(m,"epoch")!=c.originEpoch)continue;
        ++matches;s.originPath=canonicalPath(entry.path().string());
        s.originSha=zhangCheckpointFileSha256((entry.path()/"checkpoint.bundle").string());
        require(string(m,"config_sha256")==zhangCheckpointSha256("US20_GEC_USER_CHECKPOINT_V1_APPROVED_FORKS\n"+
            zhangUs20ReadFile(zhangUs20Environment("ZHANG_US20_FORK_PLAN"))),"REGISTERED_FORK_PARENT_PLAN_CHANGED");
    }
    require(matches==1,"REGISTERED_FORK_PARENT_SLOT_NOT_UNIQUE");validateOrigin(c,s);return s;
}

inline void verifyExternalManifest(Context& c)
{
    const auto receiptPath=zhangUs20Environment("ZHANG_US20_RESTORE_RECEIPT");
    require(!receiptPath.empty(),"REGISTERED_FORK_RESTORE_RECEIPT_REQUIRED");
    require(zhangCheckpointFileSha256(receiptPath)==zhangUs20Environment("ZHANG_US20_RESTORE_RECEIPT_SHA256"),
        "REGISTERED_FORK_RESTORE_RECEIPT_CHANGED");
    const auto receipt=boost::json::parse(zhangUs20ReadFile(receiptPath)).as_object();
    const auto path=string(receipt,"event_fork_manifest_path");c.externalManifestSha=string(receipt,"event_fork_manifest_sha256");
    const auto text=zhangUs20ReadFile(path);require(zhangCheckpointSha256(text)==c.externalManifestSha,
        "REGISTERED_FORK_MANIFEST_CHANGED");
    const auto m=boost::json::parse(text).as_object();
    require(string(m,"schema")=="US20_SELECTED_REGISTERED_EVENT_V2" && string(m,"plan_sha256")==c.parentPlanSha &&
        string(m,"candidate_id")==c.candidateId && string(m,"selector_sha256")==c.plan.manifestSha &&
        string(m,"origin_branch")==c.originBranch && string(m,"origin_epoch")==c.originEpoch &&
        integer(m,"stop_gpst_seconds")==c.stopEpoch,"REGISTERED_FORK_MANIFEST_WIDENS_REGISTERED_EDGE");
    c.expectedOriginSha=string(m,"origin_sha256");c.expectedOriginPath=canonicalPath(string(m,"origin_path"));
    const auto& prefix=m.at("ordered_prefix").as_array();const auto& order=c.parentPlan.at("candidate_order").as_array();
    require(!prefix.empty() && prefix.size()<=order.size(),"REGISTERED_FORK_PREFIX_SIZE_INVALID");
    for(std::size_t i=0;i<prefix.size();++i)
    {
        const auto& ref=prefix[i].as_object();const auto id=string(ref,"candidate_id");
        require(id==std::string(order[i].as_string()),"REGISTERED_FORK_SKIPPED_CANDIDATE");
        const auto candidate=registeredSpec(c.parentPlan,c.parentPlanSha,id);
        const auto& origin=c.parentPlan.at("branches").as_object().at(candidate.originBranch).as_object();
        const auto evidencePath=std::filesystem::path(string(origin,"catalog_root"))/(id+".json");
        require(canonicalPath(string(ref,"path"))==canonicalPath(evidencePath.string()) &&
            zhangCheckpointFileSha256(evidencePath.string())==string(ref,"sha256"),"REGISTERED_FORK_SELECTION_EVIDENCE_CHANGED");
        const auto record=boost::json::parse(zhangUs20ReadFile(evidencePath)).as_object();
        require(string(record,"schema")=="US20_MAIN_T_CANDIDATE_EVIDENCE_V2" && string(record,"candidate_id")==id &&
            string(record,"plan_sha256")==c.parentPlanSha && string(record,"selector_sha256")==candidate.plan.manifestSha,
            "REGISTERED_FORK_SELECTION_IDENTITY_CHANGED");
        State s{string(record,"origin_path"),string(record,"origin_sha256"),{}};validateOrigin(candidate,s);
        const auto actualOrigin=registeredOrigin(candidate);
        require(actualOrigin.originSha==s.originSha && canonicalPath(actualOrigin.originPath)==canonicalPath(s.originPath),
            "REGISTERED_FORK_SELECTION_PARENT_CHANGED");
        const auto evidence=evidenceFromJson(record.at("evidence").as_object());
        require(evidence.epoch==candidate.plan.applyEpoch,"REGISTERED_FORK_MISSING_EXPECTED_EVIDENCE");
        const auto decision=evaluate(boundPlan(candidate,s),evidence);
        require(string(record,"status")==statusName(decision.status) && string(record,"reason")==decision.reason,
            "REGISTERED_FORK_SELECTION_DECISION_CHANGED");
        if(i+1<prefix.size())require(decision.status==Status::INELIGIBLE,"REGISTERED_FORK_SKIPPED_FIRST_ELIGIBLE");
        else
        {
            require(id==c.candidateId && decision.status==Status::ELIGIBLE && s.originSha==c.expectedOriginSha &&
                canonicalPath(s.originPath)==c.expectedOriginPath,"REGISTERED_FORK_SELECTED_EVENT_NOT_FIRST_ELIGIBLE");
            c.registeredCommonPlanSha=evidence.commonPlanSha;
        }
    }
    boost::json::object exact{{"schema","US20_SELECTED_REGISTERED_EVENT_V2"},{"plan_sha256",c.parentPlanSha},
        {"candidate_id",c.candidateId},{"selector_sha256",c.plan.manifestSha},{"origin_branch",c.originBranch},
        {"origin_epoch",c.originEpoch},{"origin_path",c.expectedOriginPath},{"origin_sha256",c.expectedOriginSha},
        {"stop_gpst_seconds",c.stopEpoch},{"ordered_prefix",prefix}};
    require(exact==m,"REGISTERED_FORK_MANIFEST_UNKNOWN_FIELDS");
}

inline Context registeredContext(const boost::json::object& p,const std::string& planSha,const std::string& name)
{
    const auto& branch=p.at("branches").as_object().at(name).as_object();
    auto c=registeredSpec(p,planSha,string(branch,"candidate_id"));c.observedBranch=name;
    c.action=action(string(branch,"single_event_action"));c.group=string(branch,"single_event_group");
    require(!c.group.empty(),"SINGLE_EVENT_EXECUTION_GROUP_MISSING");
    const auto& origin=p.at("branches").as_object().at(c.originBranch).as_object();
    require(string(branch,"input_manifest_sha256")==string(origin,"input_manifest_sha256"),"REGISTERED_FORK_INPUT_IDENTITY_CHANGED");
    const auto role=string(branch,"role");
    require(role==(c.action==Action::KEEP?"SINGLE_EVENT_KEEP":"SINGLE_EVENT_RESET"),"REGISTERED_FORK_ROLE_ACTION_MISMATCH");
    const auto source=string(branch,"restore_from");
    if(source==c.originBranch)require(string(branch,"restore_epoch")==c.originEpoch,"REGISTERED_FORK_EDGE_SLOT_MISMATCH");
    else
    {
        const auto& parent=p.at("branches").as_object().at(source).as_object();
        require(string(parent,"candidate_id")==c.candidateId && string(parent,"single_event_group")==c.group &&
            string(parent,"single_event_action")==actionName(c.action) && string(parent,"role")==role,
            "REGISTERED_FORK_UNREGISTERED_OR_OTHER_EVENT_PARENT");
    }
    verifyExternalManifest(c);return c;
}

inline std::vector<Context> observers(long long epoch)
{
    std::vector<Context> out;if(!registeredMode() || !zhangUs20Environment("ZHANG_US20_DUMP_PROVENANCE").empty())return out;
    static const auto configured=[](){const auto [text,p]=zhangUs20ForkPlan();
        std::vector<Context> all;const auto name=zhangUs20Environment("ZHANG_US20_FORK_BRANCH");
        const auto& branch=p.at("branches").as_object().at(name).as_object();
        if(string(branch,"role")!="FULL_DAY_T")return all;
        require(zhangUserDatumPolicy()==ZhangUserDatumPolicy::TRANSPORT_KNOWN,"REGISTERED_FORK_OBSERVER_REQUIRES_T");
        std::set<std::string> seen;long long previous=0;
        for(const auto& v:p.at("candidate_order").as_array())
        {
            const auto id=std::string(v.as_string());auto c=registeredSpec(p,zhangCheckpointSha256(text),id);
            require(seen.insert(id).second && c.plan.applyEpoch>=previous && c.originBranch==name,
                "REGISTERED_FORK_CATALOG_ORDER_OR_ORIGIN_INVALID");
            previous=c.plan.applyEpoch;c.observedBranch=name;all.push_back(std::move(c));
        }
        require(seen.size()==p.at("single_events").as_object().size(),"REGISTERED_FORK_CATALOG_NOT_COMPLETE");return all;
    }();
    for(const auto& c:configured)if(c.plan.applyEpoch==epoch)out.push_back(c);
    return out;
}

inline void preflightRegisteredPaths(const std::string& output,const std::string& checkpoints)
{
    if(!registeredMode() || !zhangUs20Environment("ZHANG_US20_DUMP_PROVENANCE").empty())return;
    const auto [text,p]=zhangUs20ForkPlan();const auto name=zhangUs20Environment("ZHANG_US20_FORK_BRANCH");
    const auto& current=p.at("branches").as_object().at(name).as_object();
    require(canonicalPath(output)==canonicalPath(string(current,"output_root")) &&
        canonicalPath(checkpoints)==canonicalPath(string(current,"checkpoint_root")),"REGISTERED_FORK_ACTUAL_OUTPUT_ROLE_CHANGED");
    std::vector<std::string> roots;
    for(const auto& item:p.at("branches").as_object())
    {
        const auto& b=item.value().as_object();const auto root=string(b,"output_root");
        for(const auto& previous:roots)disjoint(root,previous);roots.push_back(root);
        require(containsPath(root,string(b,"checkpoint_root")) && containsPath(root,string(b,"catalog_root")),
            "REGISTERED_FORK_ARTIFACT_OUTSIDE_ROLE_ROOT");
        disjoint(string(b,"checkpoint_root"),string(b,"catalog_root"));
        for(const auto& v:p.at("protected_paths").as_array())disjoint(root,std::string(v.as_string()));
    }
    // Force selected-manifest validation before creating branch output files.
    (void)context();
}
inline bool stopAfterCommitted(long long epoch)
{
    if(!registeredMode())return false;const auto* c=context();if(!c)return false;
    require(epoch<=c->stopEpoch,"REGISTERED_FORK_OVERRAN_REGISTERED_WINDOW");
    if(epoch!=c->stopEpoch)return false;
    require(committed && committed->receipt && committed->receipt->status==Status::ELIGIBLE,
        "REGISTERED_FORK_CANNOT_FINISH_WITHOUT_TARGET_RECEIPT");return true;
}

inline void publishRecord(const std::filesystem::path& path,const std::string& text)
{
    require(!std::filesystem::exists(path),"REGISTERED_FORK_DUPLICATE_CANDIDATE_RECORD");
    std::filesystem::create_directories(path.parent_path());auto tmp=path;tmp+=".pending";
    require(!std::filesystem::exists(tmp),"REGISTERED_FORK_UNRESOLVED_PENDING_RECORD");
    {std::ofstream f(tmp,std::ios::binary);f.exceptions(std::ios::badbit|std::ios::failbit);f<<text;f.close();}
    std::filesystem::rename(tmp,path);
}

} // namespace zhangSingleEventRuntime
