// Calls the real production parser; intentionally never starts estimation.
#include "common/acsConfig.hpp"
#include "common/zhangSystemClockDatum.hpp"
#include <iostream>
#include <exception>

int main(int argc, char** argv)
{
    try
    {
        if (!configure(argc, argv)) return 2;
        std::cout << "MULTI_GNSS_SERVER_CONFIG_PROBE full_rank="
                  << acsConfig.zhangFullRank.enable
                  << " products=" << acsConfig.zhangPppAr.output_products
                  << " system_clocks=" << zhangUsesSystemReceiverClocks(acsConfig)
                  << " configured_systems=" << zhangSystemReceiverClockSystems(acsConfig).size()
                  << " no_estimation_started=1\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "CONFIG_REJECTED " << error.what() << '\n';
        return 3;
    }
}
