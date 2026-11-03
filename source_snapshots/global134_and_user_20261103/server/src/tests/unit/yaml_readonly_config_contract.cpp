#ifdef NDEBUG
#error "yaml_readonly_config_contract requires assertions; compile with -UNDEBUG"
#endif
#include "common/zhangYamlReadOnlyLookup.hpp"
#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

using Path = std::vector<std::string>;
struct Default
{
    std::string value,comment,type;
    bool found=false;
    bool operator==(const Default& other) const
    {return std::tie(value,comment,type,found)==std::tie(other.value,other.comment,other.type,other.found);}
};
struct Registry
{
    std::map<std::string,Default> defaults;
    std::map<std::string,bool> available;
    bool operator==(const Registry& other) const
    {return defaults==other.defaults && available==other.available;}
};
void addAvailable(Registry& registry,const std::string& stack)
{
    auto& available=registry.available[stack];
    if(available)return;
    available=true;
    auto pos=stack.find_last_of(':',stack.size()-2);
    if(pos!=std::string::npos)addAvailable(registry,stack.substr(0,pos+1));
}
YAML::Node query(YAML::Node base,const Path& path,Registry& registry,
    bool candidate,const std::string& comment="comment",const std::string& initial="42")
{
    YAML::Node current;current.reset(base);std::string stack;
    for(std::size_t i=0;i<path.size();++i)
    {
        const auto& key=path[i];if(key.empty())continue;
        bool test=false;
        if(candidate)
        {
            const auto child=zhangYamlReadOnlyChildPreservingFallback(current,key);
            test=static_cast<bool>(child);
            current.reset(child);
        }
        else
        {
            if(current[key])test=true;
            current.reset(current[key]);
        }
        stack+=key+":";
        if(registry.defaults.find(stack)==registry.defaults.end())
        {
            registry.defaults[stack]=i==path.size()-1 ? Default{initial,comment,"int",test} : Default{"","","",test};
        }
        if(i==path.size()-1 && registry.defaults[stack].comment.empty() && !comment.empty())
            registry.defaults[stack].comment=comment;
    }
    addAvailable(registry,stack);
    return current;
}
struct Outcome
{
    int value=42,warnings=0;
    bool found=false;
    std::string error;
    bool operator==(const Outcome& other) const
    {return std::tie(value,warnings,found,error)==std::tie(other.value,other.warnings,other.found,other.error);}
};
Outcome tryInt(YAML::Node root,const Path& path,Registry& registry,bool candidate)
{
    Outcome out;
    try
    {
        auto node=query(root,path,registry,candidate);
        try {out.value=node.as<int>();out.found=true;}
        catch(...)
        {
            if(node.IsSequence() || (node.IsScalar() && !node.Scalar().empty()))++out.warnings;
        }
    }
    catch(const YAML::BadSubscript& error) {out.error="BadSubscript:"+std::string(error.what());}
    catch(const YAML::Exception& error) {out.error="YAMLException:"+std::string(error.what());}
    return out;
}
int main()
{
    struct Case {const char* name;std::string yaml;Path path;};
    const std::vector<Case> cases{
        {"present","a: {b: {value: 2}}",{"a","b","value"}},
        {"first_missing","a: 1",{"missing","b","value"}},
        {"middle_missing","a: {}",{"a","missing","value"}},
        {"last_missing","a: {b: {}}",{"a","b","missing"}},
        {"null_parent","a: null",{"a","b","value"}},
        {"root_null","null",{"a","b"}},
        {"scalar_parent","a: 1",{"a","b","value"}},
        {"root_scalar","hello",{"a"}},
        {"leaf_scalar_bad_conversion","a: nope",{"a"}},
        {"leaf_sequence_bad_conversion","a: [1,2]",{"a"}},
        {"sequence_zero_string","a: [3,4]",{"a","0"}},
        {"sequence_one_string","a: [3,4]",{"a","1"}},
        {"sequence_leading_zero","a: [3,4]",{"a","01"}},
        {"sequence_unknown_key","a: [3,4]",{"a","missing"}},
        {"sequence_out_of_range","a: [3,4]",{"a","3"}},
        {"sequence_scalar_parent_error","a: [3,4]",{"a","1","bad"}},
        {"nested_sequence","a: [[3,4],[5,6]]",{"a","1","0"}},
        {"empty_descriptor","a: {b: 7}",{"a","","b"}}
    };
    std::size_t assertions=0;
    for(const auto& test:cases)
    {
        auto old=YAML::Load(test.yaml),candidate=YAML::Load(test.yaml);
        const auto dump=YAML::Dump(candidate);
        Registry oldRegistry,candidateRegistry;
        const auto expected=tryInt(old,test.path,oldRegistry,false);
        const auto actual=tryInt(candidate,test.path,candidateRegistry,true);
        if(!(expected==actual))
        {
            std::cerr<<"Mismatch "<<test.name<<" old="<<expected.value<<"/"<<expected.error
                <<" new="<<actual.value<<"/"<<actual.error<<"\n";
            return 1;
        }
        assert(expected==actual);++assertions;
        assert(oldRegistry==candidateRegistry);++assertions;
        assert(YAML::Dump(candidate)==YAML::Dump(old));++assertions;
        if(std::string(test.name).find("sequence")==std::string::npos &&
           std::string(test.name).find("null")==std::string::npos)
        {assert(YAML::Dump(candidate)==dump);++assertions;}
    }
    auto missing=zhangYamlReadOnlyChildPreservingFallback(YAML::Load("a: 1"),"missing");
    assert(!missing.IsDefined());++assertions;
    assert(missing.Type()==YAML::NodeType::Undefined);++assertions;
    assert(!missing.IsSequence() && !missing.IsScalar());++assertions;
    bool conversionThrown=false;
    try {(void)missing.as<int>();}catch(const YAML::BadConversion&){conversionThrown=true;}
    assert(conversionThrown);++assertions;
    Registry repeated;
    auto stable=YAML::Load("a: {value: 2}");const auto original=YAML::Dump(stable);
    for(int i=0;i<1000;++i)
    {
        const auto result=tryInt(stable,{"a","unknown_"+std::to_string(i),"value"},repeated,true);
        assert(result.value==42 && !result.found && result.error.empty());++assertions;
    }
    assert(YAML::Dump(stable)==original);++assertions;
    // Metadata stays registered for every missing descriptor, including the
    // middle and final path. Existing comment/default registration is retained.
    assert(repeated.defaults.at("a:unknown_0:value:").comment=="comment");++assertions;
    assert(repeated.defaults.at("a:unknown_0:value:").value=="42");++assertions;
    assert(repeated.available.at("a:unknown_0:value:"));++assertions;
    assert(repeated.available.at("a:unknown_0:"));++assertions;
    assert(repeated.available.at("a:"));++assertions;
    // Existing non-map mutation is intentionally retained, not simulated.
    // It must affect a later whole-vector read exactly as the original API.
    auto oldSequence=YAML::Load("a: [3,4]"),candidateSequence=YAML::Load("a: [3,4]");
    Registry oldSequenceRegistry,candidateSequenceRegistry;
    assert(tryInt(oldSequence,{"a","0"},oldSequenceRegistry,false)==
           tryInt(candidateSequence,{"a","0"},candidateSequenceRegistry,true));++assertions;
    bool oldVector=false,candidateVector=false;
    std::vector<int> oldValues{99},candidateValues{99};
    try {oldValues=oldSequence["a"].as<std::vector<int>>();oldVector=true;}catch(const YAML::Exception&){}
    try {candidateValues=candidateSequence["a"].as<std::vector<int>>();candidateVector=true;}catch(const YAML::Exception&){}
    assert(oldVector==candidateVector && oldValues==candidateValues);++assertions;
    assert(YAML::Dump(oldSequence)==YAML::Dump(candidateSequence));++assertions;
    assert(oldSequenceRegistry==candidateSequenceRegistry);++assertions;
    Registry oldComments,candidateComments;
    auto oldCommentRoot=YAML::Load("a: {}"),candidateCommentRoot=YAML::Load("a: {}");
    (void)query(oldCommentRoot,{"a","missing"},oldComments,false,"","original");
    (void)query(candidateCommentRoot,{"a","missing"},candidateComments,true,"","original");
    (void)query(oldCommentRoot,{"a","missing"},oldComments,false,"updated","replacement");
    (void)query(candidateCommentRoot,{"a","missing"},candidateComments,true,"updated","replacement");
    assert(oldComments==candidateComments);++assertions;
    assert(candidateComments.defaults.at("a:missing:").comment=="updated");++assertions;
    assert(candidateComments.defaults.at("a:missing:").value=="original");++assertions;
    std::cout<<"{\"passed\":true,\"assertions\":"<<assertions
        <<",\"fixture_cases\":"<<cases.size()<<",\"missing_query_repetitions\":1000"
        <<",\"defaults_comments_available_registration_equal\":true"
        <<",\"scalar_parent_bad_subscript_equal\":true,\"invalid_node_fallback_safe\":true"
        <<",\"map_missing_source_dump_unchanged\":true,\"nonmap_legacy_mutations_equal\":true"
        <<",\"subsequent_vector_decode_equal\":true,\"uses_production_lookup_header\":true"
        <<",\"calls_production_acsconfig_parser\":false"
        <<",\"pea_run\":false}\n";
}
