// A second native reader of archives produced by two actual PEA processes.
// No KF reconstruction, estimation, import, or "same helper twice" oracle.
#include "common/zhangCheckpoint.hpp"
#include "common/zhangUs20E240Continuation.hpp"
#include "common/zhangCertifiedTreeTransport.hpp"
#include "common/zhangProductFloatCarrierMetadata.hpp"
#include "pea/zhangReceiverCheckpoint.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/json.hpp>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
namespace fs=std::filesystem;
namespace json=boost::json;
static void require(bool value,const std::string& why)
{if(!value)throw std::runtime_error(why);}
static std::string bytes(const fs::path& path)
{std::ifstream in(path,std::ios::binary);require(bool(in),"OPEN:"+path.string());
 return {(std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>()};}
template<class T> static std::string wire(const T& value)
{std::ostringstream out(std::ios::binary);{boost::archive::binary_oarchive ar(out);ar<<value;}return out.str();}
static ZhangCheckpointBundle load(const fs::path& directory,const std::string& producer)
{
    const auto raw=bytes(directory/"checkpoint_manifest.json");
    const auto doc=json::parse(raw).as_object();
    auto field=[&](const char* name){return std::string(doc.at(name).as_string().c_str());};
    require(field("binary_sha256")==producer,"PRODUCER_HASH");
    ZhangCheckpointExpectations expected;
    expected.experimentMode=ZHANG_US20_GLOBAL134_MODE;expected.binarySha256=producer;
    expected.configSha256=zhangCheckpointSha256(field("config_text"));
    expected.inputManifestSha256=zhangCheckpointSha256(field("input_manifest_text"));
    expected.platformFingerprint=field("platform");expected.compilerFingerprint=field("compiler");
    expected.linearAlgebraFingerprint=field("linear_algebra");expected.endianness=field("endianness");
    require(field("config_sha256")==expected.configSha256 && field("input_manifest_sha256")==expected.inputManifestSha256,"TEXT_SHA");
    ZhangCheckpointBundle bundle;
    const auto read=readZhangCheckpointBundle((directory/"checkpoint.bundle").string(),expected,bundle);
    require(read.valid,"NATIVE_DECODE:"+read.failureReason);
    require(bundle.manifest.configText==field("config_text") && bundle.manifest.inputManifestText==field("input_manifest_text") &&
        bundle.manifest.checkpointId==field("checkpoint_id") && bundle.manifest.runtimeId==field("runtime_id") &&
        bundle.manifest.epoch==field("epoch"),"SIDECAR_NATIVE_IDENTITY");
    std::string failure;
    require(zhangValidateCertifiedTreeTransportMetadata(bundle.kfCore.metaDataMap,bundle.manifest.runtimeId,failure),"TRANSPORT_METADATA:"+failure);
    return bundle;
}
static bool savedFloatSequence(const ZhangCheckpointBundle& bundle,std::uint64_t& value)
{
    std::map<KFKey,int> indices;
    for(const auto& [source,index]:bundle.kfCore.kfIndexMap) {
        KFKey key;key.type=static_cast<KF>(source.type);key.Sat=source.satellite;
        key.str=source.receiver;key.num=source.number;key.comment=source.comment;
        indices.emplace(key,index);
    }
    std::string reason;
    return readZhangFloatTargetSavedFactorSequence(bundle.kfCore.metaDataMap,bundle.manifest.runtimeId,
        zhangFloatCarrierStateOrder(indices),value,reason);
}
int main(int argc,char**argv)
{
    try {
        require(argc==7,"Usage: before-directory after-directory before-producer-sha after-producer-sha report.json source-migration|same-producer|e242-uninterrupted-vs-restored");
        const fs::path beforePath(argv[1]),afterPath(argv[2]);
        const bool migration=std::string(argv[6])=="source-migration";
        const bool trajectory=std::string(argv[6])=="e242-uninterrupted-vs-restored";
        require(migration || trajectory || std::string(argv[6])=="same-producer","AUDIT_KIND");
        require(zhangUs20Global134Sha256WellFormed(argv[3]) && zhangUs20Global134Sha256WellFormed(argv[4]),"SHA_ARGUMENT");
        if(migration) {
            require(beforePath.string()==ZHANG_US20_E240_SOURCE_PATH && std::string(argv[3])==ZHANG_US20_E240_SOURCE_BINARY_SHA,"PINNED_SOURCE_PATH");
            require(zhangCheckpointFileSha256((beforePath/"checkpoint_manifest.json").string())==ZHANG_US20_E240_SOURCE_MANIFEST_SHA,"PINNED_SOURCE_MANIFEST");
            require(zhangCheckpointFileSha256((beforePath/"checkpoint.bundle").string())==ZHANG_US20_E240_SOURCE_BUNDLE_SHA,"PINNED_SOURCE_BUNDLE");
        } else require(std::string(argv[3])==argv[4],"SAME_PRODUCER_REQUIRED");
        if(trajectory) {
            const std::string id="E29-20240717020030-e242-r32f626de9463";
            const std::string root=ZHANG_US20_E240_VALIDATION_ROOT;
            require(beforePath.string()==root+"/first_pass_archive/"+id &&
                afterPath.string()==root+"/checkpoints/"+id,"E242_EXACT_COMPARISON_PATHS_REQUIRED");
        }
        const auto before=load(beforePath,argv[3]),after=load(afterPath,argv[4]);
        json::object report,fields;json::array failed,sections;
        auto check=[&](const std::string& name,bool value){fields[name]=value;if(!value)failed.emplace_back(name);};
        check("manifest.runtimeId",before.manifest.runtimeId==after.manifest.runtimeId && after.manifest.runtimeId==ZHANG_US20_GLOBAL134_RUNTIME_ID);
        check("manifest.mode",before.manifest.experimentMode==after.manifest.experimentMode);
        check("manifest.epoch",before.manifest.epoch==after.manifest.epoch);
        if(trajectory) {
            const std::string id="E29-20240717020030-e242-r32f626de9463";
            const std::string parent="E29-20240717020000-e241-r32f626de9463";
            // main.cpp updates e29LastCheckpointId after the successful e241
            // capture and sets it from the actual source ID on e241 restore.
            // Thus both complete e242 trajectories have exactly this parent.
            check("manifest.parent_and_id",before.manifest.checkpointId==id && after.manifest.checkpointId==id &&
                before.manifest.parentCheckpointId==parent && after.manifest.parentCheckpointId==parent &&
                before.manifest.epoch=="2024-07-17 02:00:30" && after.manifest.epoch=="2024-07-17 02:00:30");
        } else check("manifest.parent_and_id",after.manifest.parentCheckpointId==before.manifest.checkpointId && after.manifest.checkpointId==before.manifest.checkpointId+"-after-restore");
        check("manifest.input",before.manifest.inputManifestText==after.manifest.inputManifestText && before.manifest.inputManifestSha256==after.manifest.inputManifestSha256);
        check("manifest.ABI",before.manifest.platformFingerprint==after.manifest.platformFingerprint && before.manifest.compilerFingerprint==after.manifest.compilerFingerprint &&
            before.manifest.linearAlgebraFingerprint==after.manifest.linearAlgebraFingerprint && before.manifest.endianness==after.manifest.endianness);
        if(!migration)check("manifest.config",before.manifest.configText==after.manifest.configText);
        else {
            check("source.dimension",before.kfCore.x.size()==7274);
            check("source.config",before.manifest.configSha256==ZHANG_US20_E240_SOURCE_CONFIG_SHA);
            check("source.input",before.manifest.inputManifestSha256==ZHANG_US20_E240_SOURCE_INPUT_SHA);
            std::map<std::string,std::string> targetFields,environment;std::string reason;
            bool canonical=ZhangGec134MigrationDetail::canonicalFields(after.manifest.configText,targetFields,reason) &&
                ZhangUs20E240Detail::environmentDirectory(targetFields,environment,reason);
            if(canonical) {
                environment["ZHANG_US20_GLOBAL134_STAGE"]="AR";
                environment["ZHANG_E29_RESTORE_BUNDLE"]=ZHANG_US20_E240_SOURCE_PATH;
                canonical=validateZhangUs20E240CanonicalDelta(before.manifest.configText,after.manifest.configText,environment,reason);
            }
            check("source_to_target.explicit_canonical_allowlist",canonical);
        }
        const auto&a=before.kfCore;const auto&b=after.kfCore;
        if(hasZhangFloatTargetV2State(a.metaDataMap) || hasZhangFloatTargetV2State(b.metaDataMap)) {
            std::uint64_t sourceSequence=0,afterSequence=0;
            const bool valid=savedFloatSequence(before,sourceSequence) && savedFloatSequence(after,afterSequence);
            check("runtime.saved_float_factor_sequence_metadata",valid && sourceSequence==afterSequence);
            if(!trajectory) {
            const auto receipt=json::parse(bytes(afterPath/"restore_validation_receipt.json")).as_object();
            const auto* actual=receipt.if_contains("float_target_saved_factor_sequence");
            const auto* verified=receipt.if_contains("float_target_actual_runtime_counter_verified");
            const bool valueMatches=actual && ((actual->is_uint64() && actual->as_uint64()==afterSequence) ||
                (actual->is_int64() && actual->as_int64()>=0 && std::uint64_t(actual->as_int64())==afterSequence));
            check("runtime.actual_PEA_factor_counter_receipt",valid && valueMatches && verified && verified->is_bool() && verified->as_bool());
            }
        }
        auto matrix=[&](const char* name,const auto& x,const auto& y){
            const bool shape=x.rows()==y.rows() && x.cols()==y.cols();
            check(std::string("core.")+name,shape && (x.size()==0 || std::memcmp(x.data(),y.data(),std::size_t(x.size())*sizeof(double))==0));};
        matrix("x",a.x,b.x);matrix("P",a.P,b.P);matrix("dx",a.dx,b.dx);
        matrix("prefitRatios",a.prefitRatios,b.prefitRatios);matrix("postfitRatios",a.postfitRatios,b.postfitRatios);
#define FIELD(name) check("core." #name,wire(a.name)==wire(b.name))
        FIELD(time);FIELD(kfIndexMap);FIELD(stateTransitionMap);FIELD(gaussMarkovTauMap);FIELD(gaussMarkovMuMap);
        FIELD(procNoiseMap);FIELD(initNoiseMap);FIELD(sigmaMaxMap);FIELD(outageLimitMap);FIELD(exponentialNoiseMap);
        FIELD(pseudoStateMap);FIELD(pseudoParentMap);FIELD(errorCountMap);FIELD(filterChunkMap);FIELD(metaDataMap);
        FIELD(lsqRequired);FIELD(sigmaPass);FIELD(chiQCPass);FIELD(chi2);FIELD(dof);FIELD(chi2PerDof);FIELD(qc);
        FIELD(id);FIELD(rtsBasename);FIELD(outputResiduals);FIELD(outputMongoMeasurements);FIELD(statisticsMap);FIELD(statisticsMapSum);
#undef FIELD
        const std::map<std::string,unsigned> expected={
            {"GEC_AR_RUNTIME_V1",1},{"pea.controller.v1",1},{"pea.receiver_runtime.v1",3},
            {"pea.rinex_file_streams.v1",2},{"pea.rng_declaration.v1",1},{"pea.satellite_runtime.v1",1},
            {"zhang_ambres_runtime",3},{"zhang_graph_runtime",5},{"zhang_ppp_ar_runtime",6}};
        check("sections.count",before.sections.size()==9 && after.sections.size()==9);
        for(const auto&[name,schema]:expected) {
            const auto left=before.sections.find(name),right=after.sections.find(name);
            // The scientific section name remains stable. Its only admitted
            // wire schemas are legacy V1 and scoped V2; restoration must
            // recapture the exact source format, not normalize either side.
            const bool gec=name=="GEC_AR_RUNTIME_V1";
            const unsigned wanted=gec && left!=before.sections.end()?left->second.schemaVersion:schema;
            const bool admitted=!gec || wanted==1 || wanted==2;
            bool equal=admitted && left!=before.sections.end() && right!=after.sections.end() &&
                left->second.schemaVersion==wanted && right->second.schemaVersion==wanted;
            bool boundary=false;
            if(equal && !trajectory && name=="pea.receiver_runtime.v1") {
                const auto receipt=auditZhangReceiverCheckpointBoundaryResetEquivalence(before.manifest.runtimeId,left->second.payload,right->second.payload);
                equal=receipt.valid;boundary=true;
            } else if(equal) equal=left->second.payload==right->second.payload;
            check("section."+name,equal);
            sections.emplace_back(json::object{{"name",name},{"schema",wanted},{"entire_persistent_payload_equal",equal},
                {"declared_ready_workset_boundary_policy_only",boundary},
                {"before_payload_bytes",left==before.sections.end()?0:left->second.payload.size()},
                {"after_payload_bytes",right==after.sections.end()?0:right->second.payload.size()},
                {"before_payload_sha256",left==before.sections.end()?"MISSING":left->second.sha256},
                {"after_payload_sha256",right==after.sections.end()?"MISSING":right->second.sha256}});
        }
        report["checks"]=std::move(fields);report["failed_fields"]=std::move(failed);
        report["sections"]=std::move(sections);report["passed"]=report.at("failed_fields").as_array().empty();
        report["audit_kind"]=argv[6];report["state_dimension"]=a.x.size();report["P_entries_raw_IEEE"]=a.P.size();
        report["estimation_executed_by_auditor"]=false;report["whole_process_restore_must_be_proven_by_PEAs_native_receipt"]=true;
        report["actual_epoch_boundary_comparison_not_restore_only_recapture"]=trajectory;
        report["receiver_workset_normalization_used"]=!trajectory;
        report["scientific_accuracy_accepted"]=false;
        std::ofstream out(argv[5],std::ios::binary);out<<json::serialize(report)<<'\n';out.flush();require(bool(out),"REPORT_WRITE");
        std::cout<<"US20_E240_NINE_SECTION_AUDIT passed="<<report.at("passed")<<" core_fields=33 sections=9\n";
        return report.at("passed").as_bool()?0:1;
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
