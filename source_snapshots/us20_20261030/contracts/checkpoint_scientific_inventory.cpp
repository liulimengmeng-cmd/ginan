// Standalone read-only evidence extractor. Include the existing DTO definitions
// and link without the production receiver-checkpoint object. No solver entry
// point is executed, no capture/restore acceptance is inferred from extraction.
#include "pea/zhangReceiverCheckpoint.cpp"
#include <boost/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

extern "C" int __wrap_main(int argc,char** argv)
{
    try
    {
        namespace fs=std::filesystem;
        if(argc!=4)throw std::runtime_error("usage: extractor checkpoint_directory output_directory branch_output_root");
        const fs::path input=argv[1],output=argv[2];const std::string branchRoot=argv[3];
        if(!fs::create_directory(output))throw std::runtime_error("fresh output directory required");
        auto read=[](const fs::path& p){std::ifstream f(p,std::ios::binary);if(!f)throw std::runtime_error("read "+p.string());
            return std::string(std::istreambuf_iterator<char>(f),{});};
        const auto manifest=boost::json::parse(read(input/"checkpoint_manifest.json")).as_object();
        auto field=[&](const char* key){return std::string(manifest.at(key).as_string());};
        ZhangCheckpointExpectations expected{field("experiment_mode"),field("binary_sha256"),field("config_sha256"),
            field("input_manifest_sha256"),field("platform"),field("compiler"),field("linear_algebra"),field("endianness")};
        ZhangCheckpointBundle bundle;
        const auto result=readZhangCheckpointBundle((input/"checkpoint.bundle").string(),expected,bundle);
        if(!result.valid)throw std::runtime_error(result.failureReason);
        boost::json::object inventory{{"input_bundle_sha256",zhangCheckpointFileSha256((input/"checkpoint.bundle").string())},
            {"binary_sha256",field("binary_sha256")},{"config_sha256",field("config_sha256")},
            {"input_manifest_sha256",field("input_manifest_sha256")},{"epoch",field("epoch")}};
        boost::json::array mapped;
        auto artifact=[&](std::string& value,const std::string& role)
        {
            if(value.empty())return;
            if(!value.starts_with(branchRoot+"/"))throw std::runtime_error("unmapped artifact path "+role+":"+value);
            mapped.push_back(boost::json::object{{"role",role},{"original",value}});
            value="<BRANCH_OUTPUT>"+value.substr(branchRoot.size());
        };
        boost::json::object payloads;
        auto write=[&](const std::string& name,const std::string& bytes)
        {
            if(name.find('/')!=std::string::npos)throw std::runtime_error("bad section name");
            std::ofstream f(output/(name+".bin"),std::ios::binary);f.write(bytes.data(),bytes.size());
            if(!f)throw std::runtime_error("write failed");
            payloads[name]=boost::json::object{{"sha256",zhangCheckpointSha256(bytes)},{"bytes",bytes.size()}};
        };
        artifact(bundle.kfCore.rtsBasename,"kfCore.rtsBasename");
        for(const char* role:{"traceFilename","posFilenameMTA1","posFilenamePKRD","posFilenameP800"})
        {
            const auto found=bundle.kfCore.metaDataMap.find(role);
            if(found!=bundle.kfCore.metaDataMap.end())artifact(found->second,string("kfCore.metadata.")+role);
        }
        write("kfCore",serializeZhangCheckpointSectionPayload(bundle.kfCore));
        for(const auto& [name,section]:bundle.sections)
        {
            std::string bytes=section.payload,failure;
            if(name==ZHANG_RECEIVER_RUNTIME_CHECKPOINT_SECTION_NAME)
            {
                ZhangReceiverRuntimeEnvelope envelope;
                if(!deserializePayload(bytes,envelope,failure,"INVENTORY"))throw std::runtime_error(failure);
                for(auto& receiver:envelope.receivers)
                {
                    boost::json::object bindings{
                        {"receiver",receiver.id},
                        {"site",receiver.sinex.siteBindingFingerprint},
                        {"receiver_binding",receiver.sinex.receiverBindingFingerprint},
                        {"antenna",receiver.sinex.antennaBindingFingerprint},
                        {"eccentricity",receiver.sinex.eccentricityBindingFingerprint},
                        {"default_site",siteBindingFingerprint(&dummySiteid)},
                        {"default_receiver",receiverBindingFingerprint(&dummyReceiver)},
                        {"default_antenna",antennaBindingFingerprint(&dummyAntenna)},
                        {"default_eccentricity",eccentricityBindingFingerprint(&dummySiteEcc)}};
                    std::ofstream(output/(receiver.id+"_sinex_binding_evidence.json"))<<boost::json::serialize(bindings)<<'\n';
                    const auto prefix="receiver."+receiver.mapKey+".";
                    artifact(receiver.traceFilename,prefix+"traceFilename");
                    artifact(receiver.jsonTraceFilename,prefix+"jsonTraceFilename");
                    artifact(receiver.sppOutputFile,prefix+"sppOutputFile");
                    artifact(receiver.solution.sppState.rtsBasename,prefix+"solution.sppState.rtsBasename");
                }
                if(!serializePayload(envelope,bytes,failure,"INVENTORY"))throw std::runtime_error(failure);
            }
            else if(name==ZHANG_SATELLITE_RUNTIME_CHECKPOINT_SECTION_NAME)
            {
                ZhangSatelliteRuntimeEnvelope envelope;
                if(!deserializePayload(bytes,envelope,failure,"INVENTORY"))throw std::runtime_error(failure);
                for(auto& satellite:envelope.satellites)
                {
                    artifact(satellite.traceFilename,"satellite."+satellite.satellite.id()+".traceFilename");
                    artifact(satellite.jsonTraceFilename,"satellite."+satellite.satellite.id()+".jsonTraceFilename");
                }
                if(!serializePayload(envelope,bytes,failure,"INVENTORY"))throw std::runtime_error(failure);
            }
            else if(name=="us20.regional_runtime.v1")
            {
                std::istringstream stream(bytes,std::ios::binary);boost::archive::binary_iarchive archive(stream);
                std::string tag,settingsIdentity;archive & tag & settingsIdentity;
                if(tag!="US20_REGIONAL_RUNTIME_V1" || settingsIdentity.size()!=64)throw std::runtime_error("regional header");
                const auto offset=static_cast<std::size_t>(stream.tellg());
                inventory["regional_readonly_settings_identity"]=settingsIdentity;
                // Exact remaining payload includes completion time, sequence,
                // actual cache mean/Q/root, noise identity and all properties.
                // The excluded header is retained and must be justified by the
                // separate immutable future-product transformation manifest.
                bytes=bytes.substr(offset);
            }
            else if(name=="us20.trop_ledger.v1")
            {
                const auto ledger=boost::json::parse(bytes).as_object();
                inventory["trop_ledger_original"]=ledger;
                // Retain the raw payload here. An explicit JSON-field audit,
                // not this extractor, decides allowed plan/branch identities.
            }
            write(name,bytes);
        }
        inventory["mapped_artifact_fields"]=mapped;inventory["payloads"]=payloads;
        inventory["scope"]="Read-only structured inventory; no scientific equivalence decision";
        std::ofstream(output/"inventory.json")<<boost::json::serialize(inventory)<<'\n';
        std::cout<<"CHECKPOINT_INVENTORY_EXTRACTED\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<"INVENTORY_FAILED "<<e.what()<<'\n';return 1;}
}
