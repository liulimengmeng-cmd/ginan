#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/receiver.hpp"
#include "common/zhangSystemClockDatum.hpp"
#include "pea/zhangPppAr.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void updateRecClocks(Trace&, ReceiverMap&, KFState&);
void updateSppPppClkOffsets(ReceiverMap&, KFState&);
// Bound by the test runner to a copied object containing the actual private
// production writer. No production symbol or binary is changed for testing.
extern "C" void zhang_multi_gnss_test_append_product(const ZhangInternalProduct&);

namespace
{
void require(bool condition, const std::string& reason)
{
    if (!condition) throw std::runtime_error(reason);
}
struct Parameter { const char* satellite; const char* kind; const char* code; };
const Parameter parameters[] = {
    {"G01", "CLOCK", "NONE"}, {"G01", "PHASE", "L1C"}, {"G01", "PHASE", "L2W"},
    {"E01", "CLOCK", "NONE"}, {"E01", "PHASE", "L1C"}, {"E01", "PHASE", "L5Q"},
    {"C19", "CLOCK", "NONE"}, {"C19", "PHASE", "L1P"}, {"C19", "PHASE", "L5P"}};

std::vector<std::string> fixtureFields(const std::string& row)
{
    std::vector<std::string> values;
    std::istringstream input(row);
    std::string value;
    while (std::getline(input, value, ',')) values.push_back(value);
    return values;
}

void fixtureRow(std::ostream& stream, const std::vector<std::string>& fields)
{
    for (std::size_t i = 0; i < fields.size(); ++i)
        stream << (i ? "," : "") << fields[i];
    stream << '\n';
}

void actualProductRoundTrip(const std::filesystem::path& output)
{
    acsConfig.zhangPppAr.user_adapter = true;
    acsConfig.zhangPppAr.product_solution = "PRODUCT_FIXED";
    const auto written = output / "actual_writer_products.csv";
    std::filesystem::remove(written); // exact dedicated test fixture, never production input
    acsConfig.zhangPppAr.product_filename = written.string();
    ZhangInternalProduct original;
    original.time.bigTime = 19000;
    original.valid_from.bigTime = 18970;
    original.solution = "PRODUCT_FIXED";
    original.satellite = SatSys("E01");
    original.observable = E_ObsCode::L5Q;
    original.clock_m = 4.5;
    original.phase_m = 1.25;
    original.correction_m = 3.25;
    original.clock_sigma_m = 0.2;
    original.phase_sigma_m = 0.1;
    original.correction_sigma_m = 0.15;
    original.numeric_valid = original.branch_valid = original.continuity_valid = true;
    original.ppp_usable = original.pppar_usable = original.ar_valid = true;
    original.integer_valid = original.dual_frequency_ar_valid = true;
    original.integer_structure_valid = original.integer_datum_continuous = true;
    original.integer_precision_valid = original.cycle_closure_valid = true;
    original.discontinuity_counter = 2;
    original.datum_version = 3;
    original.product_state = "AR_VALID";
    original.same_segment_continuous = true;
    original.public_phase_segment_id = "PUBLIC_E01_L5Q_1";
    auto& contract = original.published_continuity;
    contract.previousEpochGpstSeconds = 18970;
    contract.previousCounter = 2;
    contract.previousDatumVersion = 3;
    contract.previousPublicSegment = original.public_phase_segment_id;
    contract.floatTransportValid = true;
    contract.floatCorrectionShiftMetres = 0;
    contract.integerTransportValid = true;
    contract.integerTransportShiftCycles = 0;
    contract.proofId = "INDEPENDENT_READER_CONTRACT_PROOF";
    contract.kind = "IDENTITY";
    contract.jointCovarianceComplete = true;
    contract.bundleId = "READER_CONTRACT_BUNDLE";
    zhang_multi_gnss_test_append_product(original);
    ZhangInternalProduct actual;
    require(queryZhangInternalProduct(original.time, original.satellite, original.observable, actual),
            "actual 65-column writer product was not consumed");
    require(actual.clock_m == original.clock_m && actual.phase_m == original.phase_m &&
            actual.public_phase_segment_id == original.public_phase_segment_id &&
            actual.published_continuity.previousEpochGpstSeconds == 18970 &&
            actual.published_continuity.previousCounter == 2 &&
            actual.published_continuity.previousDatumVersion == 3 &&
            actual.published_continuity.previousPublicSegment == original.public_phase_segment_id &&
            actual.published_continuity.floatTransportValid &&
            actual.published_continuity.floatCorrectionShiftMetres == 0 &&
            actual.published_continuity.integerTransportValid &&
            actual.published_continuity.integerTransportShiftCycles == 0 &&
            actual.published_continuity.proofId == contract.proofId &&
            actual.published_continuity.kind == "IDENTITY" &&
            actual.published_continuity.jointCovarianceComplete &&
            actual.published_continuity.bundleId == contract.bundleId,
            "writer/reader round trip dropped a continuity field");
    std::ifstream input(written);
    std::string headerLine, productLine;
    std::getline(input, headerLine);
    std::getline(input, productLine);
    const auto header = fixtureFields(headerLine);
    const auto row = fixtureFields(productLine);
    require(header.size() == 65 && row.size() == 65,
            "fixture did not come from the actual 65-column writer");
    auto queryFixture = [&](const std::string& name, const std::vector<std::string>& names,
                            const std::vector<std::string>& values, bool duplicate = false)
    {
        const auto filename = output / name;
        std::ofstream stream(filename);
        fixtureRow(stream, names);
        fixtureRow(stream, values);
        if (duplicate) fixtureRow(stream, values);
        stream.close();
        acsConfig.zhangPppAr.product_filename = filename.string();
        return queryZhangInternalProduct(original.time, original.satellite, original.observable, actual);
    };
    auto reorderedHeader = header, reorderedRow = row;
    std::reverse(reorderedHeader.begin(), reorderedHeader.end());
    std::reverse(reorderedRow.begin(), reorderedRow.end());
    require(queryFixture("reordered_products.csv", reorderedHeader, reorderedRow),
            "reader is still position-dependent for the current schema");
    auto affineRow = row;
    affineRow[58] = "0.375";
    affineRow[60] = "-3";
    affineRow[62] = "CERTIFIED_AFFINE";
    require(queryFixture("affine_products.csv", header, affineRow) &&
            actual.published_continuity.floatCorrectionShiftMetres == 0.375 &&
            actual.published_continuity.integerTransportShiftCycles == -3,
            "affine transport was flattened to a fabricated identity");
    auto unknownRow = row;
    unknownRow[57] = "0";
    unknownRow[58] = "nan";
    unknownRow[59] = "0";
    unknownRow[60] = "nan";
    unknownRow[61] = "UNAVAILABLE";
    unknownRow[62] = "UNKNOWN";
    require(queryFixture("unknown_products.csv", header, unknownRow) &&
            !actual.published_continuity.floatTransportValid &&
            !actual.published_continuity.integerTransportValid &&
            std::isnan(actual.published_continuity.floatCorrectionShiftMetres),
            "unknown transport was rejected or fabricated instead of preserved");
    auto bad = row;
    bad[52] = "2";
    require(!queryFixture("unsupported_contract.csv", header, bad), "unknown contract version accepted");
    bad = row; bad[57] = "2";
    require(!queryFixture("invalid_contract_flag.csv", header, bad), "non-Boolean transport flag accepted");
    bad = row; bad[63] = "0";
    require(!queryFixture("usable_without_q.csv", header, bad), "usable current product without complete Q accepted");
    auto badHeader = header;
    badHeader.back() = "unknown_continuity_extension";
    require(!queryFixture("unknown_contract_column.csv", badHeader, row), "unknown schema column accepted");
    badHeader = header; badHeader.pop_back();
    bad = row; bad.pop_back();
    require(!queryFixture("incomplete_contract.csv", badHeader, bad), "truncated current schema accepted");
    require(!queryFixture("duplicate_products.csv", header, row, true), "duplicate product identity accepted");
    const std::vector<std::string> legacyHeader(header.begin(), header.begin() + 52);
    const std::vector<std::string> legacyRow(row.begin(), row.begin() + 52);
    require(queryFixture("legacy_52_products.csv", legacyHeader, legacyRow) &&
            actual.public_phase_segment_id == original.public_phase_segment_id &&
            !actual.published_continuity.floatTransportValid,
            "legacy 52-column products were broken or granted fabricated transport");
}

double independentQ(int row, int column)
{
    return (row == column ? 1.0 + row * 0.1 : 0.0) +
        0.03 * (row + 1) * (column + 1);
}

void writeSidecar(const std::filesystem::path& path, long epoch,
                  bool missingCross, bool conflicting, bool negativeDiagonal)
{
    std::ofstream stream(path);
    stream << "epoch,solution,row_satellite,row_parameter,row_observable,"
              "column_satellite,column_parameter,column_observable,covariance\n";
    auto write = [&](int row, int column, double value)
    {
        const auto& a = parameters[row];
        const auto& b = parameters[column];
        stream << epoch << ",FLOAT," << a.satellite << ',' << a.kind << ',' << a.code
               << ',' << b.satellite << ',' << b.kind << ',' << b.code << ','
               << std::setprecision(17) << value << '\n';
    };
    for (int row = 0; row < 9; ++row)
        for (int column = row; column < 9; ++column)
        {
            if (missingCross && row == 0 && column == 3) continue;
            write(row, column, negativeDiagonal && row == 0 && column == 0
                ? -1.0 : independentQ(row, column));
        }
    if (conflicting) write(0, 3, independentQ(0, 3) + 0.5);
    require(stream.good(), "could not write independent covariance fixture");
}

void bindFixture(const std::filesystem::path& filename)
{
    acsConfig.zhangPppAr.user_adapter = true;
    acsConfig.zhangPppAr.user_use_full_product_covariance = true;
    acsConfig.zhangPppAr.product_solution = "FLOAT";
    acsConfig.zhangPppAr.product_filename = (filename.parent_path() / "no_product_input").string();
    acsConfig.zhangPppAr.product_covariance_filename = filename.string();
}

void bindJointClockModel()
{
    acsConfig.zhangFullRank.enable = true;
    acsConfig.zhangPppAr.user_adapter = false;
    acsConfig.process_sys.clear();
    for (const auto system : {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS})
    {
        acsConfig.process_sys[system] = true;
        acsConfig.zhangFullRank.sysOpts[system].baseline_observables =
            {E_ObsCode::L1C, E_ObsCode::L2W};
    }
}

void actualClockAdjustment()
{
    bindJointClockModel();
    acsConfig.adjust_rec_clocks_by_spp = true;
    acsConfig.adjust_clocks_for_jumps_only = false;
    ReceiverMap receivers;
    auto& receiver = receivers["CLOCK_CONTRACT"];
    receiver.id = "CLOCK_CONTRACT";
    receiver.sol.status = E_Solution::SINGLE;
    receiver.sol.clkRefSys = E_Sys::GPS;
    receiver.sol.sppClk = 20;
    receiver.sol.sppPppClkOffset = 2;
    receiver.sol.sppPppClkOffsetSys = E_Sys::GPS;
    receiver.sol.clkAdjustReady = true;
    auto& options = acsConfig.recOptsMap[receiver.id];
    options._initialised = true;
    options.clk.estimate = {true};
    options.clk.sigma = {1};
    options.clk_rate.estimate = {true};
    options.clk_rate.sigma = {1};
    KFState state;
    state.time.bigTime = 1000;
    tsync.bigTime = 1030;
    state.kfIndexMap.clear();
    state.kfIndexMap[KFState::oneKey] = 0;
    state.x = VectorXd::Zero(7);
    state.x[0] = 1;
    state.P = MatrixXd::Identity(7, 7);
    state.P(0, 0) = 0;
    const E_Sys systems[] = {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS};
    for (int i = 0; i < 3; ++i)
    {
        const auto clock = zhangReceiverClockKey(acsConfig, receiver.id, systems[i]);
        const auto rate = zhangReceiverClockKey(acsConfig, receiver.id, systems[i],
                                               0, KF::REC_CLOCK_RATE);
        state.kfIndexMap[clock] = i + 1;
        state.kfIndexMap[rate] = i + 4;
        state.x[i + 1] = 10 * (i + 1);
        state.x[i + 4] = 0.1 * (i + 1);
        // Production keeps a persistent STM for existing states. Merely
        // inserting kfIndexMap/x is not a valid initialized KF lifecycle:
        // setKFTrans would create a "new" state and add init.x-currentX.
        state.addKFState(clock, {.x = state.x[i + 1], .P = 1});
        state.addKFState(rate, {.x = state.x[i + 4], .P = 1});
        state.setKFTransRate(clock, rate, 1);
    }
    std::ostringstream trace;
    const auto beforeUnverified = state.stateTransitionMap;
    const auto beforeMean = state.x;
    const auto beforeCovariance = state.P;
    updateRecClocks(trace, receivers, state);
    // This intentionally contains only the old scalar fields, not a real SPP
    // solution or an external-chart receipt. It must no longer authorize the
    // former +9 m adjustment. Actual successful per-system SPP/postfilter and
    // physical-jump/rate paths belong to multi_gnss_clock_lifecycle.cpp, which
    // calls production spp instead of fabricating chart metadata here.
    require(state.stateTransitionMap == beforeUnverified &&
            (state.x-beforeMean).norm() == 0 && (state.P-beforeCovariance).norm() == 0,
            "unverified legacy scalar adjusted an effective receiver clock");
    const auto beforeChangedReference = state.stateTransitionMap;
    receiver.sol.clkRefSys = E_Sys::GAL;
    updateRecClocks(trace, receivers, state);
    require(state.stateTransitionMap == beforeChangedReference,
            "changed SPP reference reused stale GPS offset");
    updateSppPppClkOffsets(receivers, state);
    require(receiver.sol.sppPppClkOffsetSys == E_Sys::GPS &&
            std::abs(receiver.sol.sppPppClkOffset-2) < 1e-12,
            "postfilter invented new source-chart authority without actual SPP evidence");
    updateRecClocks(trace, receivers, state);
    require(state.stateTransitionMap == beforeChangedReference,
            "unverified changed-system scalar bypassed the source-chart boundary");
}

void actualClockOutputGauge()
{
    bindJointClockModel();
    acsConfig.reference_clock = "A";
    acsConfig.reference_bias = "NO_REFERENCE";
    KFState state;
    state.kfIndexMap.clear();
    state.kfIndexMap[KFState::oneKey] = 0;
    state.x = VectorXd::Zero(10);
    state.x[0] = 1;
    state.dx = VectorXd::Zero(10);
    MatrixXd factor = MatrixXd::Zero(10, 10);
    for (int row = 1; row < 10; ++row)
        for (int column = 1; column < 10; ++column)
            factor(row, column) = (row == column ? 2.0 : 0.0) +
                0.02 * (row + 1) * (column + 2);
    state.P = factor * factor.transpose();
    const E_Sys systems[] = {E_Sys::GPS, E_Sys::GAL, E_Sys::BDS};
    for (int i = 0; i < 3; ++i)
    {
        state.kfIndexMap[zhangReceiverClockKey(acsConfig, "A", systems[i])] = 1 + i;
        state.kfIndexMap[zhangReceiverClockKey(acsConfig, "B", systems[i])] = 4 + i;
        KFKey satellite;
        satellite.type = KF::SAT_CLOCK;
        satellite.Sat = SatSys(systems[i], 10);
        state.kfIndexMap[satellite] = 7 + i;
    }
    for (int i = 1; i < 10; ++i) state.x[i] = i * i;
    std::ostringstream trace;
    const auto output = propagateUncertainty(trace, state);
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
        {
            const int row = output.kfIndexMap.at(zhangReceiverClockKey(acsConfig, "B", systems[a]));
            const int column = output.kfIndexMap.at(zhangReceiverClockKey(acsConfig, "B", systems[b]));
            const double expectedMean = state.x[4 + a] - state.x[1 + a];
            const double expectedCovariance = state.P(4 + a, 4 + b) -
                state.P(1 + a, 4 + b) - state.P(4 + a, 1 + b) + state.P(1 + a, 1 + b);
            require(std::abs(output.x[row] - expectedMean) < 1e-12,
                    "actual output gauge subtracted another constellation reference");
            require(std::abs(output.P(row, column) - expectedCovariance) < 1e-11,
                    "actual output gauge dropped cross-constellation covariance");
        }
}
}

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2, "supply a dedicated fixture output directory");
        std::filesystem::path output(argv[1]);
        std::filesystem::create_directories(output);
        // Test the actual production reader, not a second matrix loader.
        unsetenv("ZHANG_R51_ENABLE");
        actualClockAdjustment();
        actualClockOutputGauge();
        acsConfig.zhangFullRank.enable = false;
        actualProductRoundTrip(output);
        const auto complete = output / "complete.csv";
        writeSidecar(complete, 12345, false, false, false);
        bindFixture(complete);
        std::vector<std::vector<double>> rows;
        const SatSys satellites[] = {SatSys("G01"), SatSys("E01"), SatSys("C19")};
        const E_ObsCode codes[] = {E_ObsCode::L1C, E_ObsCode::L1C, E_ObsCode::L1P};
        GTime time;
        time.bigTime = 12345;
        for (int system = 0; system < 3; ++system)
        {
            std::vector<double> row;
            std::string reason;
            int rank = 0;
            require(queryZhangInternalProductNoiseFactors(time, satellites[system],
                codes[system], true, row, &rank, &reason), "valid joint Q rejected: " + reason);
            require(rank == 9 && row.size() == 9, "global Q rank was split by system");
            rows.push_back(std::move(row));
        }
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
            {
                double reconstructed = 0;
                for (int factor = 0; factor < 9; ++factor)
                    reconstructed += rows[a][factor] * rows[b][factor];
                const int ca = 3 * a, pa = ca + 1, cb = 3 * b, pb = cb + 1;
                const double expected = independentQ(ca, cb) - independentQ(ca, pb) -
                                        independentQ(pa, cb) + independentQ(pa, pb);
                require(std::abs(reconstructed - expected) < 1e-11,
                        "production factor rows lost clock/phase cross-system Q");
            }

        struct BadCase { const char* file; const char* reason; bool missing, conflict, negative; };
        const BadCase cases[] = {
            {"missing_cross.csv", "COVARIANCE_INCOMPLETE_JOINT_TRIANGLE", true, false, false},
            {"conflicting.csv", "COVARIANCE_CONFLICTING_CANONICAL_PAIR", false, true, false},
            {"non_psd.csv", "COVARIANCE_NOT_POSITIVE_SEMIDEFINITE", false, false, true}};
        // First change the input file without changing epoch/solution. The
        // production reader must not reuse the previous file's cached Q.
        long epoch = 12345;
        for (const auto& bad : cases)
        {
            const auto filename = output / bad.file;
            writeSidecar(filename, epoch, bad.missing, bad.conflict, bad.negative);
            bindFixture(filename);
            time.bigTime = epoch++;
            std::vector<double> row;
            std::string reason;
            require(!queryZhangInternalProductNoiseFactors(time, satellites[0],
                codes[0], true, row, nullptr, &reason), "invalid Q accepted");
            require(reason == bad.reason, "unexpected covariance rejection: " + reason);
            require(row.empty(), "rejected covariance returned usable noise factors");
        }
        std::cout << "multi_gnss_clock_product_covariance: PASS actual_production_reader=1 "
                     "joint_primitive_dimension=9 cross_system_clock_phase=1\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "multi_gnss_clock_product_covariance: FAIL " << error.what() << '\n';
        return 1;
    }
}
