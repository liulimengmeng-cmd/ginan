#include "common/zhangRatioOnly.hpp"
#include "common/zhangRatioGate.hpp"
#include "common/zhangR51FloatReuse.hpp"

#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
static void setSetting(const char* name, const char* value)
{
#ifdef _WIN32
    _putenv_s(name, value ? value : "");
#else
    if (value) setenv(name, value, 1);
    else unsetenv(name);
#endif
}
struct RestoreEnvironment
{
    std::map<std::string, std::pair<bool, std::string>> settings;
    RestoreEnvironment()
    {
        for (const char* name : {"ZHANG_R51_ENABLE", "ZHANG_R51_RATIO_ONLY", "ZHANG_R51_RATIO_THRESHOLD"})
        {
            const char* value = std::getenv(name);
            settings[name] = {value != nullptr, value ? value : ""};
        }
    }
    ~RestoreEnvironment()
    {
        for (const auto& [name, setting] : settings)
            setSetting(name.c_str(), setting.first ? setting.second.c_str() : nullptr);
    }
};
int main()
{
    RestoreEnvironment restore;
    try
    {
        setSetting("ZHANG_R51_ENABLE", "1");
        setSetting("ZHANG_R51_RATIO_ONLY", "1");
        setSetting("ZHANG_R51_RATIO_THRESHOLD", nullptr);
        require(zhangR51RatioPolicyValid(3), "unset override rejected configured ratio=3");
        require(zhangR51EffectiveRatioThreshold(3) == 3, "unset override changed configured ratio=3");
        require(zhangR51EffectiveRatioThreshold(4) == 4, "unset override changed arbitrary configured ratio");
        setSetting("ZHANG_R51_RATIO_THRESHOLD", "2");
        require(zhangR51RatioPolicyValid(3), "explicit ratio=2 rejected");
        const double threshold = zhangR51EffectiveRatioThreshold(3);
        require(threshold == 2, "effective ratio did not use explicit override");
        require(zhangEvaluateRatioGate(true, true, 1, 2.1640035246414073, threshold).accepted,
            "observed G18 ratio=2.164 rejected at threshold=2");
        require(zhangEvaluateRatioGate(true, true, 1, 2.0917795579036786, threshold).accepted,
            "observed G18 ratio=2.092 rejected at threshold=2");
        require(!zhangEvaluateRatioGate(true, true, 1, 1.2679336606487015, threshold).accepted,
            "observed G08 ratio=1.268 accepted at threshold=2");
        require(!zhangEvaluateRatioGate(true, true, 1, 2.1640035246414073, 3).accepted,
            "configured ratio=3 reference did not reject ratio=2.164");
        require(!zhangEvaluateRatioGate(false, true, 1, 20, threshold).accepted,
            "incomplete top two accepted under ratio=2");
        require(!zhangEvaluateRatioGate(true, false, 1, 20, threshold).accepted,
            "non-distinct candidates accepted under ratio=2");
        require(!zhangEvaluateRatioGate(true, true, 1, 1, threshold).accepted,
            "tied candidates accepted under ratio=2");
        require(zhangEvaluateRatioGate(true, true, 0, 1, threshold).accepted,
            "valid zero best distance lost infinite-ratio semantics");

        for (const char* invalid : {"", "nan", "inf", "-2", "1", "0", "2junk", " 2", "2 ", "1e999"})
        {
            setSetting("ZHANG_R51_RATIO_THRESHOLD", invalid);
            require(!zhangR51RatioPolicyValid(3), "invalid override silently accepted");
            require(std::isnan(zhangR51EffectiveRatioThreshold(3)), "invalid override silently fell back");
            require(!zhangEvaluateRatioGate(true, true, 1, 20,
                zhangR51EffectiveRatioThreshold(3)).accepted, "invalid policy produced a ratio pass");
        }
        setSetting("ZHANG_R51_ENABLE", "0");
        require(zhangR51RatioPolicyValid(3) && zhangR51EffectiveRatioThreshold(3) == 3,
            "non-R51 path consumed experiment override");
        setSetting("ZHANG_R51_ENABLE", "1");
        setSetting("ZHANG_R51_RATIO_ONLY", "0");
        require(zhangR51RatioPolicyValid(3) && zhangR51EffectiveRatioThreshold(3) == 3,
            "non-ratio-only path consumed experiment override");

        std::map<std::string, std::string> checkpointPolicy;
        require(zhangR51FixtureRatioPolicyAllowed(checkpointPolicy, false),
            "pinned ratio=3 fixture default rejected");
        checkpointPolicy["ZHANG_R51_RATIO_THRESHOLD"] = "2";
        checkpointPolicy["ZHANG_R51_EFFECTIVE_RATIO_THRESHOLD"] = "2";
        require(zhangR51FixtureRatioPolicyAllowed(checkpointPolicy, true),
            "unfixed FLOAT source rejected explicit new ratio=2 experiment");
        require(!zhangR51FixtureRatioPolicyAllowed(checkpointPolicy, false),
            "fixed ratio=3 source migrated as same-policy ratio=2");
        checkpointPolicy["ZHANG_R51_EFFECTIVE_RATIO_THRESHOLD"] = "3";
        require(!zhangR51FixtureRatioPolicyAllowed(checkpointPolicy, true),
            "checkpoint recorded threshold disagreed with override");
        checkpointPolicy["ZHANG_R51_RATIO_THRESHOLD"] = "nan";
        require(!zhangR51FixtureRatioPolicyAllowed(checkpointPolicy, true),
            "invalid checkpoint policy allowed FLOAT migration");
        std::cout << "PASS R51 ratio=2 policy, unchanged defaults, invalid settings and checkpoint policy isolation\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
