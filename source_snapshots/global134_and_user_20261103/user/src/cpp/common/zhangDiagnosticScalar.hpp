#pragma once

#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

// Machine-readable diagnostics must not inherit precision/floatfield flags
// left on the TRACE stream by an earlier epoch or a checkpoint capture.
// A private stream also leaves the caller's presentation state untouched.
inline std::string zhangDiagnosticDouble(double value)
{
    std::ostringstream text;
    text << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return text.str();
}
