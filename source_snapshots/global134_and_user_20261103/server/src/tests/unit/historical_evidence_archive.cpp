#include "common/zhangHistoricalEvidenceArchive.hpp"
#include <openssl/evp.h>
#include <iostream>
#include <stdexcept>

// Same SHA-256 primitive as checkpoint production, with no KF/PEA linkage.
#ifndef ZHANG_ARCHIVE_PRODUCTION_LINK
std::string zhangCheckpointSha256(const std::string& bytes)
{
    unsigned char digest[EVP_MAX_MD_SIZE];unsigned int size=0;
    if(!EVP_Digest(bytes.data(),bytes.size(),digest,&size,EVP_sha256(),nullptr))return {};
    std::ostringstream out;out<<std::hex<<std::setfill('0');
    for(unsigned int i=0;i<size;++i)out<<std::setw(2)<<unsigned(digest[i]);
    return out.str();
}
#endif

static unsigned int checks=0;
static void require(bool value,const std::string& why)
{++checks;if(!value)throw std::runtime_error(why);}

int main(int argc,char** argv)
try
{
    using namespace zhang_history_archive;
    require(argc==2,"requires new evidence-output directory");
    const auto root=std::filesystem::absolute(argv[1]);
    require(!std::filesystem::exists(root) && std::filesystem::create_directories(root),"output already exists");
    const std::string parentBytes="COMPLETE_SNAPSHOT_FIXTURE_WITH_NONZERO_CROSS_Q_NOT_PZZ_ONLY";
    const auto parent=root/"full-checkpoint.fixture";
    require(writeImmutable(parent,parentBytes),"durable parent write failed");
    const auto collision=root/"somebody-else.tmp";
    require(writeImmutable(collision,"OTHER_WRITER_BYTES"),"temp collision fixture setup failed");
    require(!writeImmutableWithTemporary(root/"collision-final","MUST_NOT_TRUNCATE",collision),
        "exclusive temp collision overwrote somebody else's file");
    std::string collisionBytes;
    require(readBytes(collision,1024,collisionBytes) && collisionBytes=="OTHER_WRITER_BYTES" &&
        !std::filesystem::exists(root/"collision-final"),"temp collision mutated existing temp/final");
    require(!writeImmutableWithTemporary(parent,"MUST_NOT_REPLACE",root/"late-install.tmp"),
        "no-replace installation overwrote existing immutable target");
    require(readBytes(parent,1024,collisionBytes) && collisionBytes==parentBytes,
        "immutable target bytes changed after rejected installation");
    require(!readBytes(parent,parentBytes.size()-1,collisionBytes),"bounded reader ignored declared byte limit");
    const Scope scope{"AUTHORITATIVE_OWNER","2024-07-17 02:27:30","CHECKPOINT_291","CHECKPOINT_290",zhangCheckpointSha256(parentBytes)};
    auto verify=[&](const std::string& expected,std::string& reason)
    {
        std::string actual;
        const bool valid=readBytes(parent,maximumObjectBytes,actual) && zhangCheckpointSha256(actual)==expected;
        reason=valid?"NONE":"MISSING_OR_CORRUPT_FULL_SNAPSHOT";return valid;
    };
    const std::vector<Evidence> evidence{
        {"TRANSPORT_ONLY_CERTIFICATE","CERT/BDS","complete physical versions, RHS, scopes and proof closure",{scope.snapshotPayloadSha256}},
        {"FLOAT_TARGET_DEFINITION","FLOAT/C27-L5","FALK-C28-L5P-V0 whole target definition; no Gaussian substitution",{scope.snapshotPayloadSha256}}};
    std::string reason;
    const auto positive=root/"positive";
    require(write(positive,scope,evidence,verify,reason),"archive write failed: "+reason);
    std::vector<Evidence> restored;
    require(read(positive,scope,verify,restored,reason) && restored==evidence,"archive lost exact payload/parents/scope");
    require(reason=="HISTORY_VERIFIED_TRANSPORT_ONLY_NO_CURRENT_AR_AUTHORITY","read exposed an authorization role");
    require(!write(positive,scope,evidence,verify,reason),"immutable archive was overwritten");
    auto unchanged=restored;
    auto wrong=scope;wrong.epoch="2024-07-17 02:28:00";
    require(!read(positive,wrong,verify,restored,reason) && restored==unchanged,"wrong epoch read committed output");
    wrong=scope;wrong.owner="OTHER_OWNER";
    require(!read(positive,wrong,verify,restored,reason) && restored==unchanged,"wrong owner read committed output");
    wrong=scope;wrong.parentCheckpointId="DIFFERENT_PARENT_CHECKPOINT";
    require(!read(positive,wrong,verify,restored,reason),"wrong parent identity accepted");
    auto missingParent=verify;
    std::filesystem::rename(parent,root/"parent.unavailable");
    require(!read(positive,scope,missingParent,restored,reason) && restored==unchanged,"missing snapshot parent became evidence");
    std::filesystem::rename(root/"parent.unavailable",parent);
    {
        std::ofstream output(parent,std::ios::binary|std::ios::app);output<<"CORRUPTION";
    }
    require(!read(positive,scope,verify,restored,reason),"corrupt complete snapshot parent accepted");
    {std::ofstream output(parent,std::ios::binary|std::ios::trunc);output<<parentBytes;}
    const auto partial=root/"partial";
    std::filesystem::create_directories(partial/"objects");
    const auto object=encode(evidence.front());const auto id=zhangCheckpointSha256(object);
    require(writeImmutable(partial/"objects"/(id+".payload"),object),"partial fixture setup failed");
    require(!read(partial,scope,verify,restored,reason),"objects without committed manifest accepted");
    require(writeImmutable(partial/"manifest",std::string("partially written")),"partial manifest setup failed");
    require(!read(partial,scope,verify,restored,reason),"partial manifest accepted");
    const auto corrupted=root/"corrupted";
    require(write(corrupted,scope,evidence,verify,reason),"corruption fixture write failed");
    {std::ofstream output(corrupted/"objects"/(id+".payload"),std::ios::binary|std::ios::trunc);output<<object.substr(0,object.size()/2);}
    require(!read(corrupted,scope,verify,restored,reason) && restored==unchanged,"partial object committed evidence");
    auto forged=evidence;forged.front().parents.front()=std::string(64,'0');
    require(!write(root/"wrong-object-parent",scope,forged,verify,reason),"object with wrong full-parent hash committed");
    auto invalid=scope;invalid.parentCheckpointId=invalid.checkpointId;
    require(!write(root/"self-parent",invalid,evidence,verify,reason),"self-parent checkpoint identity accepted");
    unsigned calls=0;
    auto changesBeforeCommit=[&](const std::string&,std::string& why){why="CANCELLED_OR_PARENT_CHANGED";return ++calls==1;};
    require(!write(root/"interrupted-before-manifest",scope,evidence,changesBeforeCommit,reason),
        "interrupted parent verification committed manifest");
    require(!std::filesystem::exists(root/"interrupted-before-manifest"/"manifest"),"unverified candidate obtained committed manifest");
    std::cout<<"HISTORICAL_EVIDENCE_ARCHIVE_PASS checks="<<checks
        <<" complete_parent_required=1 partial_corrupt_scope_parent_negatives=1 transport_only=1 no_solver=1\n";
    return 0;
}
catch(const std::exception& error)
{std::cerr<<"HISTORICAL_EVIDENCE_ARCHIVE_FAIL checks="<<checks<<" reason="<<error.what()<<'\n';return 1;}
