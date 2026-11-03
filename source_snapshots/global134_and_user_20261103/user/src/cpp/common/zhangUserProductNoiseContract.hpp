#pragma once
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

struct ZhangUserProductNoiseContract
{
    std::string noiseId;
    std::string source;
    std::string units = "metre";
    std::string covarianceTarget = "PUBLISHED_CLOCK_PHASE_COVARIANCE";
    std::string statisticalTarget = "UNSPECIFIED_BY_LEGACY_CSV";
    std::string temporalParent = "UNAVAILABLE";
    std::string stateProductCrossCovariance = "UNAVAILABLE";
    std::string representation = "OBSERVATION_COVARIANCE";
    std::string sample;
    // Parameter descriptors are in the actual covariance row order, not a
    // guessed CLOCK/L1/L2 order. Factors share one identity across all rows.
    std::vector<std::string> parameterOrder;
    std::vector<double> parameterWeights;
    std::vector<double> factors;
    int numericalRank = 0;
};

// Multiple measurement rows may share an error, but may not charge the same
// error once in R and again as an augmented state. Registry lifetime is one
// measurement sample, so storage does not grow with elapsed history.
inline void zhangUserClaimProductNoise(std::map<std::string,std::string>& metadata,
    const ZhangUserProductNoiseContract& packet)
{
    if (packet.noiseId.empty() || packet.sample.empty() ||
        (packet.representation!="OBSERVATION_COVARIANCE" && packet.representation!="AUGMENTED_STATE"))
        throw std::runtime_error("INVALID_PRODUCT_NOISE_CONTRACT");
    const std::string prefix="USER_PRODUCT_NOISE_OWNER:";
    const auto sample=metadata.find(prefix+"SAMPLE");
    if (sample==metadata.end() || sample->second!=packet.sample)
    {
        for (auto it=metadata.begin();it!=metadata.end();)
            if (it->first.compare(0,prefix.size(),prefix)==0) it=metadata.erase(it); else ++it;
        metadata[prefix+"SAMPLE"]=packet.sample;
    }
    const auto key=prefix+packet.noiseId;
    const auto owner=metadata.find(key);
    if (owner!=metadata.end() && owner->second!=packet.representation)
        throw std::runtime_error("PRODUCT_NOISE_DOUBLE_OWNERSHIP_R_AND_STATE");
    metadata[key]=packet.representation;
}
