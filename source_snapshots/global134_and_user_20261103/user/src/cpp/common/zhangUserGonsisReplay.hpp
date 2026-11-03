#pragma once
#include "common/eigenIncluder.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Diagnostic-only payload. Basis rows map ORIGINAL input integer coordinates
// to reduced coordinates. Never label a reduced candidate as an input vector.
struct ZhangUserGonsisLambdaTopTwo
{
    MatrixXd integerBasis;
    VectorXd firstReduced;
    VectorXd secondReduced;
    VectorXd firstInput;
    VectorXd secondInput;
    double firstSquaredCost = std::numeric_limits<double>::quiet_NaN();
    double secondSquaredCost = std::numeric_limits<double>::quiet_NaN();
    bool inputCoordinatesValid = false;
};

// Exact, diagnostic-only recovery for the full square integer basis used by
// GONSIS. Partial suffixes cannot uniquely recover all original coordinates.
// No floating inverse/rounding is permitted to certify an input candidate.
inline bool zhangUserGonsisRecoverInputInteger(
    const MatrixXd& basis, const VectorXd& reduced, VectorXd& input)
{
    using Integer = boost::multiprecision::cpp_int;
    using Rational = boost::rational<Integer>;
    const int n = basis.rows();
    input.resize(0);
    constexpr double exactLimit = 9007199254740991.0;
    if (n <= 0 || basis.cols() != n || reduced.size() != n ||
        !basis.allFinite() || !reduced.allFinite()) return false;
    std::vector<std::vector<Rational>> augmented(n, std::vector<Rational>(n + 1));
    std::vector<std::vector<Integer>> original(n, std::vector<Integer>(n));
    std::vector<Integer> rhs(n);
    auto integer = [&](double value, Integer& out)
    {
        if (std::abs(value) > exactLimit || value != std::floor(value)) return false;
        out = static_cast<long long>(value);
        return true;
    };
    for (int row = 0; row < n; ++row)
    {
        if (!integer(reduced(row), rhs[row])) return false;
        augmented[row][n] = Rational(rhs[row]);
        for (int column = 0; column < n; ++column)
        {
            if (!integer(basis(row, column), original[row][column])) return false;
            augmented[row][column] = Rational(original[row][column]);
        }
    }
    for (int column = 0; column < n; ++column)
    {
        int pivot = column;
        while (pivot < n && augmented[pivot][column].numerator() == 0) ++pivot;
        if (pivot == n) return false;
        if (pivot != column) std::swap(augmented[pivot], augmented[column]);
        const Rational value = augmented[column][column];
        for (int j = column; j <= n; ++j) augmented[column][j] /= value;
        for (int row = 0; row < n; ++row)
        {
            if (row == column || augmented[row][column].numerator() == 0) continue;
            const Rational multiple = augmented[row][column];
            for (int j = column; j <= n; ++j)
                augmented[row][j] -= multiple * augmented[column][j];
        }
    }
    const Integer limit = static_cast<long long>(exactLimit);
    std::vector<Integer> answer(n);
    VectorXd recovered(n);
    for (int row = 0; row < n; ++row)
    {
        const auto& value = augmented[row][n];
        if (value.denominator() != 1 || value.numerator() > limit || value.numerator() < -limit)
            return false;
        answer[row] = value.numerator();
        recovered(row) = answer[row].convert_to<double>();
    }
    for (int row = 0; row < n; ++row)
    {
        Integer check = 0;
        for (int column = 0; column < n; ++column) check += original[row][column] * answer[column];
        if (check != rhs[row]) return false;
    }
    input = std::move(recovered);
    return true;
}

template<typename Left, typename Right>
inline bool zhangUserGonsisMatrixBitsEqual(const Eigen::MatrixBase<Left>& left,
                                         const Eigen::MatrixBase<Right>& right)
{
    if (left.rows() != right.rows() || left.cols() != right.cols()) return false;
    if constexpr (Left::IsRowMajor != Right::IsRowMajor) return false;
    // All cache-owned/key matrices below are contiguous MatrixXd/VectorXd,
    // with identical storage order. Signed zero is deliberately not merged.
    return left.size() == 0 || std::memcmp(left.derived().data(), right.derived().data(),
                                          static_cast<std::size_t>(left.size()) * sizeof(double)) == 0;
}

struct ZhangUserGonsisIlsKeyView
{
    const VectorXd& mean;
    const MatrixXd& covariance;
    const MatrixXd& wideLaneRows;
    const VectorXd& wideLaneIntegers;
    const MatrixXd& l1Rows;
    const MatrixXd& ifRows;
};

// Local to ONE frozen domain call. No state, product identities or integers
// are cached across epochs/components; incomplete enumeration is not cached.
template<typename Solution>
class ZhangUserGonsisExactIlsCache
{
    struct Entry
    {
        VectorXd mean;
        MatrixXd covariance;
        MatrixXd wideLaneRows;
        VectorXd wideLaneIntegers;
        MatrixXd l1Rows;
        MatrixXd ifRows;
        Solution solution;
    };
    std::vector<Entry> entries;
public:
    const Solution* find(const ZhangUserGonsisIlsKeyView& key) const
    {
        for (const auto& entry : entries)
            if (zhangUserGonsisMatrixBitsEqual(entry.mean, key.mean) &&
                zhangUserGonsisMatrixBitsEqual(entry.covariance, key.covariance) &&
                zhangUserGonsisMatrixBitsEqual(entry.wideLaneRows, key.wideLaneRows) &&
                zhangUserGonsisMatrixBitsEqual(entry.wideLaneIntegers, key.wideLaneIntegers) &&
                zhangUserGonsisMatrixBitsEqual(entry.l1Rows, key.l1Rows) &&
                zhangUserGonsisMatrixBitsEqual(entry.ifRows, key.ifRows)) return &entry.solution;
        return nullptr;
    }
    void remember(const ZhangUserGonsisIlsKeyView& key, const Solution& solution)
    {
        if (!solution.complete || !std::isfinite(solution.bestSquaredCost) ||
            !std::isfinite(solution.secondSquaredCost)) return;
        entries.push_back({key.mean, key.covariance, key.wideLaneRows, key.wideLaneIntegers,
                           key.l1Rows, key.ifRows, solution});
    }
    std::size_t size() const { return entries.size(); }
};

inline bool zhangUserGonsisReplaySwitch(const char* name, bool fallback)
{
    const char* value = std::getenv(name);
    if (!value || !*value) return fallback;
    if (std::string(value) == "1") return true;
    if (std::string(value) == "0") return false;
    throw std::runtime_error(std::string("INVALID_GONSIS_BOOLEAN_CONTROL:") + name);
}

struct ZhangUserGonsisReplayControls
{
    bool enabled = false;
    std::set<std::string> receivers;
    std::size_t maximumEpochs = 2; // zero explicitly requests the whole selected record
    std::string startGpst;
    struct Counter { std::string lastEpoch; std::size_t epochs = 0; };
    std::map<std::string, Counter> counters;
    static ZhangUserGonsisReplayControls fromEnvironment()
    {
        ZhangUserGonsisReplayControls out;
        out.enabled = zhangUserGonsisReplaySwitch("ZHANG_USER_GONSIS_REPLAY", false);
        if (!out.enabled) return out;
        if (const char* value = std::getenv("ZHANG_USER_GONSIS_REPLAY_RECEIVERS"))
        {
            std::istringstream list(value);
            std::string receiver;
            while (std::getline(list, receiver, ','))
            {
                if (receiver.empty() || receiver.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos)
                    throw std::runtime_error("INVALID_GONSIS_REPLAY_RECEIVER_LIST");
                out.receivers.insert(receiver);
            }
        }
        if (const char* value = std::getenv("ZHANG_USER_GONSIS_REPLAY_MAX_EPOCHS"))
        {
            const std::string text(value);
            if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("INVALID_GONSIS_REPLAY_EPOCH_LIMIT");
            out.maximumEpochs = std::stoull(text);
            if (out.maximumEpochs > 1000000) throw std::runtime_error("GONSIS_REPLAY_EPOCH_LIMIT_TOO_LARGE");
        }
        if (const char* value = std::getenv("ZHANG_USER_GONSIS_REPLAY_START_GPST"))
        {
            out.startGpst = value;
            const auto& text = out.startGpst;
            bool valid = text.size() == 19;
            if (valid)
                for (std::size_t index = 0; index < text.size(); ++index)
                {
                    const char separator = index == 4 || index == 7 ? '-' :
                        index == 10 ? 'T' : index == 13 || index == 16 ? ':' : '\0';
                    valid &= separator ? text[index] == separator :
                        text[index] >= '0' && text[index] <= '9';
                }
            if (valid)
                valid = text.substr(5,2) >= "01" && text.substr(5,2) <= "12" &&
                    text.substr(8,2) >= "01" && text.substr(8,2) <= "31" &&
                    text.substr(11,2) <= "23" && text.substr(14,2) <= "59" && text.substr(17,2) <= "59";
            if (!valid)
                throw std::runtime_error("INVALID_GONSIS_REPLAY_START_GPST");
        }
        return out;
    }
    bool select(const std::string& runtime, const std::string& receiver, const std::string& epoch)
    {
        if (!enabled || !receivers.count(receiver) || (!startGpst.empty() && epoch < startGpst)) return false;
        auto& counter = counters[runtime + ":" + receiver];
        if (counter.lastEpoch != epoch)
        {
            counter.lastEpoch = epoch;
            ++counter.epochs;
        }
        return maximumEpochs == 0 || counter.epochs <= maximumEpochs;
    }
};

inline std::string zhangUserGonsisReplayHex(const std::string& value)
{
    constexpr char digits[] = "0123456789ABCDEF";
    std::string out;
    for (unsigned char character : value)
    {
        out.push_back(digits[character >> 4]);
        out.push_back(digits[character & 15]);
    }
    return out.empty() ? "EMPTY" : out;
}

template<typename Matrix>
inline void zhangUserGonsisReplayMatrix(std::ostream& stream, const std::string& context,
                                      const char* name, const Eigen::MatrixBase<Matrix>& value)
{
    const auto flags = stream.flags();
    const auto precision = stream.precision();
    const auto locale = stream.getloc();
    stream.imbue(std::locale::classic());
    stream << "\nZHANG_USER_GONSIS_REPLAY_MATRIX " << context << " name=" << name
           << " rows=" << value.rows() << " columns=" << value.cols()
           << " layout=ROW_MAJOR decimal_digits=" << std::numeric_limits<double>::max_digits10
           << " data=" << std::defaultfloat << std::setprecision(std::numeric_limits<double>::max_digits10);
    bool separator = false;
    for (Eigen::Index row = 0; row < value.rows(); ++row)
        for (Eigen::Index column = 0; column < value.cols(); ++column)
        {
            if (separator) stream << ',';
            stream << value(row, column);
            separator = true;
        }
    if (!separator) stream << "EMPTY";
    stream.flags(flags);
    stream.precision(precision);
    stream.imbue(locale);
}
