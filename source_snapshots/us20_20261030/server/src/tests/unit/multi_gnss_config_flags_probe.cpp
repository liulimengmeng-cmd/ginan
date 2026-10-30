// Read the production parser's effective flags; never starts PEA estimation.
#include "common/acsConfig.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    if (!configure(argc, argv)) return 2;
    std::cout << "MULTI_GNSS_CONFIG_FLAGS_PROBE"
              << " multi_epoch_relink_shadow="
              << acsConfig.zhangPppAr.multi_epoch_relink_shadow
              << " whitened_wl_fixed_lag_shadow="
              << acsConfig.zhangPppAr.whitened_wl_fixed_lag_shadow
              << " fixed_lag_factor_capture_shadow="
              << acsConfig.zhangPppAr.fixed_lag_factor_capture_shadow
              << " targeted_besd_capture_shadow="
              << acsConfig.zhangPppAr.targeted_besd_capture_shadow
              << " user_adapter=" << acsConfig.zhangPppAr.user_adapter
              << " no_estimation_started=1\n";
    return 0;
}
