// Required only in US20 checkpoint mode. Store the actual eigenfactor, rather
// than recomputing a numerically different square root after restoration.
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

namespace
{
std::string us20RegionalSettingsIdentity()
{
    const auto& s=settings();std::ostringstream text;text<<std::setprecision(17);
    text<<s.enabled<<' '<<std::quoted(s.root.string())<<'\n';
    for(const auto& [id,g]:s.geometry)text<<std::quoted(id)<<' '<<std::quoted(g.role)<<' '
        <<g.xyz(0)<<' '<<g.xyz(1)<<' '<<g.xyz(2)<<'\n';
    for(const auto& [epoch,e]:s.epochs)text<<epoch<<' '<<std::quoted(e.directory.string())<<' '<<e.commitHash<<'\n';
    return zhangCheckpointSha256(text.str());
}
template<class A> void us20RegionalMatrix(A& archive,MatrixXd& value)
{
    std::uint64_t rows=value.rows(),cols=value.cols();archive & rows & cols;
    if(rows>4096 || cols>4096)reject("US20_CACHE_MATRIX_DIMENSIONS");
    if constexpr(A::is_loading::value)value.resize(rows,cols);
    for(std::uint64_t i=0;i<rows*cols;++i)archive & value.data()[i];
    if(!value.allFinite())reject("US20_CACHE_MATRIX_NONFINITE");
}
template<class A> void us20RegionalCache(A& archive,EpochCache& value)
{
    archive & value.epoch & value.solution & value.noiseId & value.reason & value.numericRoot
        & value.sourceDimension & value.loadSequence & value.minimumScaledEigenvalue & value.valid & value.properties;
    std::uint64_t count=value.parameters.size();archive & count;
    if(count>4096)reject("US20_CACHE_PARAMETER_COUNT");
    if constexpr(A::is_loading::value)value.parameters.resize(count);
    for(auto& p:value.parameters)archive & p.parameter & p.receiver & p.satellite & p.observable & p.fileIndex & p.mean & p.eligible;
    MatrixXd mean=value.mean;
    us20RegionalMatrix(archive,mean);us20RegionalMatrix(archive,value.covariance);us20RegionalMatrix(archive,value.root);
    if(mean.cols()!=1)reject("US20_CACHE_MEAN_DIMENSIONS");
    if constexpr(A::is_loading::value)
    {
        value.mean=mean;value.index.clear();
        for(int i=0;i<static_cast<int>(value.parameters.size());++i)
        {
            const auto& p=value.parameters[i];
            if(!value.index.emplace(std::make_tuple(p.parameter,p.receiver,p.satellite,p.observable),i).second)
                reject("US20_CACHE_DUPLICATE_PARAMETER");
        }
    }
    else
    {
        if(value.index.size()!=count)reject("US20_CACHE_INDEX_COUNT");
        for(int i=0;i<static_cast<int>(count);++i)
        {
            const auto& p=value.parameters[i];
            if(value.index.at({p.parameter,p.receiver,p.satellite,p.observable})!=i)reject("US20_CACHE_INDEX_ORDER");
        }
    }
    if(!std::isfinite(value.minimumScaledEigenvalue) || value.sourceDimension<0 || value.sourceDimension>4096)
        reject("US20_CACHE_SOURCE_DIMENSION");
    if(value.valid && (value.mean.size()!=count || value.covariance.rows()!=count || value.covariance.cols()!=count ||
        value.root.rows()!=count || value.noiseId.empty() || value.reason!="NONE"))reject("US20_CACHE_VALID_SHAPE");
}
}

std::string exportZhangUs20RegionalCheckpoint(double completed)
{
    const std::lock_guard<std::recursive_mutex> lock(regionalCacheMutex);
    if(!settings().enabled)reject("US20_CACHE_REGIONAL_REQUIRED");
    std::ostringstream stream(std::ios::binary);boost::archive::binary_oarchive archive(stream);
    const std::string tag="US20_REGIONAL_RUNTIME_V1",identity=us20RegionalSettingsIdentity();
    archive & tag & identity & completed & regionalLoadSequence;
    // Capture never writes to the live cache, counters, filter or ledger.
    EpochCache copy=cache;us20RegionalCache(archive,copy);return stream.str();
}

void importZhangUs20RegionalCheckpoint(const std::string& payload,double completed,bool dryRun)
{
    const std::lock_guard<std::recursive_mutex> lock(regionalCacheMutex);
    if(payload.empty() || payload.size()>512*1024*1024)reject("US20_CACHE_PAYLOAD_SIZE");
    std::istringstream stream(payload,std::ios::binary);boost::archive::binary_iarchive archive(stream);
    std::string tag,identity;double savedCompleted=0;std::uint64_t sequence=0;
    archive & tag & identity & savedCompleted & sequence;
    if(tag!="US20_REGIONAL_RUNTIME_V1" || identity!=us20RegionalSettingsIdentity() || savedCompleted!=completed)
        reject("US20_CACHE_IDENTITY_MISMATCH");
    EpochCache candidate;us20RegionalCache(archive,candidate);
    if(stream.peek()!=std::char_traits<char>::eof() || candidate.epoch>completed ||
       candidate.loadSequence!=sequence)reject("US20_CACHE_BOUNDARY_OR_TRAILING_BYTES");
    if(candidate.valid)
    {
        const auto location=settings().epochs.find(candidate.epoch);
        if(location==settings().epochs.end() || candidate.solution!=acsConfig.zhangPppAr.product_solution)
            reject("US20_CACHE_SOURCE_EPOCH_OR_SOLUTION");
        const auto& source=location->second;const auto bytes=readText(source.directory/"COMMIT.json");
        if(zhangCheckpointSha256(bytes)!=source.commitHash)reject("US20_CACHE_COMMIT_CHANGED");
        const auto doc=object(bytes);verifyHashes(source,doc);
        if(candidate.noiseId!="REGIONAL_JOINT:"+std::to_string(candidate.epoch)+":"+candidate.solution+":"+
            jsonString(doc,"root")+":"+source.commitHash)reject("US20_CACHE_NOISE_ID");
    }
    if(dryRun)return;
    cache=std::move(candidate);regionalLoadSequence=sequence;
    if(exportZhangUs20RegionalCheckpoint(completed)!=payload)reject("US20_CACHE_RESTORE_NOT_EXACT");
}
