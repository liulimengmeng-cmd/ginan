#include "common/zhangDiagnosticScalar.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>

int main()
{
    int checks = 0;
    const std::array<double, 8> values{
        17.675078751701641, -.11299442828851625, 0., -0.,
        86399.5, 1e-200, 1e200, std::nextafter(1., 2.)};
    try
    {
        for (double value : values)
        {
            const auto text = zhangDiagnosticDouble(value);
            if (std::bit_cast<std::uint64_t>(std::stod(text)) != std::bit_cast<std::uint64_t>(value))
                throw std::runtime_error("binary64 diagnostic did not round trip");
            ++checks;
            for (int precision : {6, 17})
            for (auto format : {std::ios_base::fmtflags(0), std::ios_base::fixed, std::ios_base::scientific})
            {
                std::ostringstream trace;
                trace.precision(precision); trace.setf(format, std::ios_base::floatfield);
                const auto flags = trace.flags();
                trace << zhangDiagnosticDouble(value);
                if (trace.str() != text || trace.precision() != precision || trace.flags() != flags)
                    throw std::runtime_error("diagnostic depended on or changed external stream flags");
                ++checks;
            }
        }
        std::cout << "PASS checks=" << checks << " scope=diagnostic_roundtrip_and_stream_isolation\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "FAIL checks=" << checks << " " << e.what() << '\n';
        return 1;
    }
}
