#include "common/zhangUserArExperiment.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>

static unsigned assertions = 0;

#define CHECK(condition) do { ++assertions; if (!(condition)) { \
    std::cerr << "FAIL line " << __LINE__ << ": " #condition "\n"; return EXIT_FAILURE; \
} } while (false)

int main()
{
    // Independent boundary fixtures; the helper only receives copies of scalar
    // geometry, so candidate rejection cannot edit the FLOAT observations.
    constexpr double rad = 3.14159265358979323846 / 180.0;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(zhangUserArElevationDecision(true, 15, 15 * rad, true).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 14.999 * rad, true).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, nan, true).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 30 * rad, false).eligible);
    CHECK(!zhangUserArElevationDecision(true, nan, 30 * rad, true).eligible);
    CHECK(!zhangUserArElevationDecision(true, -1, 30 * rad, true).eligible);
    CHECK(!zhangUserArElevationDecision(true, 91, 30 * rad, true).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 91 * rad, true).eligible);
    CHECK(zhangUserArElevationDecision(true, 0, 0, true).eligible);
    CHECK(zhangUserArElevationDecision(false, 15, 14.999 * rad, true).eligible);
    CHECK(zhangUserArElevationDecision(false, 15, nan, false).eligible);
    CHECK(!zhangUserArElevationDecision(false, 15, nan, false).enabled);
    CHECK(zhangUserArElevationDecision(true, 15, 30 * rad, true, 0.5, 0.5).eligible);
    CHECK(zhangUserArElevationDecision(true, 15, 30 * rad, true, -0.5, 0.5).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 30 * rad, true, 0.500001, 0.5).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 30 * rad, true, -30, 0.5).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 30 * rad, true, nan, 0.5).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 30 * rad, true, 0, nan).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 30 * rad, true, 0, -1).eligible);
    CHECK(zhangUserArElevationDecision(false, 15, 30 * rad, true, -30, 0.5).eligible);
    CHECK(std::string(zhangUserArElevationDecision(true, 15, 30 * rad, true, -30, 0.5).reason)
        == "STALE_OBSERVATION_GEOMETRY");
    // The native SD candidate gate consumes the same scalar contract before
    // choosing its reference. Missing/currently excluded observations cannot
    // become candidates merely because a retained KF ambiguity still exists.
    CHECK(std::string(zhangUserArElevationDecision(true, 15, nan, false, nan, 0.5).reason)
        == "INVALID_OBSERVATION_GEOMETRY");
    CHECK(std::string(zhangUserArElevationDecision(true, 15, 14.999 * rad, true, 0, 0.5).reason)
        == "BELOW_AR_ELEVATION_MASK");
    CHECK(std::string(zhangUserArElevationDecision(true, 15, 15 * rad, true, 0, 0.5).reason)
        == "AR_ELEVATION_ACCEPTED");
    CHECK(!zhangUserArElevationDecision(true, 15, 30 * rad, true, -0.500001, 0.5).eligible);
    CHECK(!zhangUserArElevationDecision(true, 15, 30 * rad, true,
        std::numeric_limits<double>::infinity(), 0.5).eligible);
    std::vector<double> floatObservationElevations{14.999 * rad, 15 * rad, 35 * rad};
    const auto untouchedFloatObservations = floatObservationElevations;
    std::vector<std::size_t> arCandidates;
    for (std::size_t i = 0; i < floatObservationElevations.size(); ++i)
        if (zhangUserArElevationDecision(true, 15, floatObservationElevations[i], true).eligible)
            arCandidates.push_back(i);
    CHECK(floatObservationElevations == untouchedFloatObservations);
    CHECK((arCandidates == std::vector<std::size_t>{1, 2}));
    // A low IF pivot cannot be implicitly reintroduced into an accepted set.
    const auto lowPivot = zhangUserArElevationDecision(true, 15, floatObservationElevations[0], true);
    CHECK(!lowPivot.eligible);
    CHECK(std::find(arCandidates.begin(), arCandidates.end(), 0) == arCandidates.end());

    const ZhangUserArExperimentKey arcKey{"run-a", "NICO", "G20"};
    const ZhangUserArExperimentKey gateKey{"run-a", "NICO", "GPS"};
    const auto identity = zhangUserArArcIdentity("l1-arc-1", "p1-v1", "l2-arc-1", "p2-v1");
    CHECK(!identity.empty());
    CHECK(zhangUserArArcIdentity("", "p1", "l2", "p2").empty());
    CHECK(zhangUserArArcIdentity("a:b", "c", "d", "e") !=
        zhangUserArArcIdentity("a", "b:c", "d", "e"));

    ZhangUserArExperimentArcTracker arcs;
    auto arc = arcs.update(arcKey, 0, identity);
    CHECK(arc.valid && !arc.eligible && arc.count == 1 && arc.age == 0);
    CHECK(arcs.update(arcKey, 0, identity).count == 1);
    for (int i = 1; i <= 4; ++i) arc = arcs.update(arcKey, i * 30, identity);
    CHECK(arc.eligible && arc.count == 5 && arc.age == 120 && arc.first == 0 && arc.last == 120);
    CHECK(arcs.update(arcKey, 120, identity).count == 5);
    arc = arcs.update(arcKey, 150, identity + "changed");
    CHECK(arc.count == 1 && !arc.eligible && arc.first == 150 && arc.reason == "IDENTITY_CHANGED");
    arc = arcs.update(arcKey, 196, identity + "changed");
    CHECK(arc.count == 1 && arc.reason == "OBSERVATION_GAP");
    arc = arcs.update(arcKey, 190, identity + "changed");
    CHECK(!arc.valid && arc.reason == "TIME_REVERSED" && arcs.size() == 0);
    arcs.update(arcKey, 210, identity);
    CHECK(!arcs.update(arcKey, 240, identity, false).valid && arcs.size() == 0);
    arcs.update(arcKey, 270, identity);
    CHECK(!arcs.update(arcKey, 300, "").valid && arcs.size() == 0);
    arcs.update(arcKey, 300, identity);
    CHECK(!arcs.update(arcKey, std::numeric_limits<double>::quiet_NaN(), identity).valid);
    CHECK(arcs.size() == 0);

    std::vector<ZhangUserArNamedInteger> rows = {
        {"G20", identity, "+0002", "-0003"},
        {"G25", identity, "0", "123456789012345678901234567890"}};
    const auto signature = zhangUserArCandidateSignature("G11", identity, rows);
    CHECK(!signature.empty());
    std::reverse(rows.begin(), rows.end());
    CHECK(signature == zhangUserArCandidateSignature("G11", identity, rows));
    rows[1].wideLane = "2";
    rows[1].first = "-3";
    CHECK(signature == zhangUserArCandidateSignature("G11", identity, rows));
    auto changedRows = rows;
    changedRows[0].first = "4";
    const auto changed = zhangUserArCandidateSignature("G11", identity, changedRows);
    CHECK(signature != changed);
    CHECK(signature != zhangUserArCandidateSignature("G12", identity, rows));
    CHECK(signature != zhangUserArCandidateSignature("G11", identity + "x", rows));
    changedRows = rows;
    changedRows[0].identity += "x";
    CHECK(signature != zhangUserArCandidateSignature("G11", identity, changedRows));
    changedRows = rows;
    changedRows.pop_back();
    CHECK(signature != zhangUserArCandidateSignature("G11", identity, changedRows));
    CHECK(zhangUserArCandidateSignature("G11", identity, {}).empty());
    changedRows = rows;
    changedRows.push_back(rows.front());
    CHECK(zhangUserArCandidateSignature("G11", identity, changedRows).empty());
    changedRows = rows;
    changedRows[0].first = "1.0";
    CHECK(zhangUserArCandidateSignature("G11", identity, changedRows).empty());
    std::string integer;
    CHECK(zhangUserArCanonicalInteger("-000", integer) && integer == "0");
    CHECK(!zhangUserArCanonicalInteger("+", integer));
    CHECK(!zhangUserArCanonicalInteger(" 1", integer));

    ZhangUserArFirstFeedbackGate gate;
    auto decision = gate.observe(gateKey, 0, signature);
    CHECK(decision.count == 1 && !decision.eligible && !decision.released);
    CHECK(gate.observe(gateKey, 0, signature).count == 1);
    CHECK(gate.observe(gateKey, 30, signature).count == 2);
    decision = gate.observe(gateKey, 60, signature);
    CHECK(decision.count == 3 && decision.eligible && !gate.released(gateKey));
    CHECK(gate.observe(gateKey, 60, signature).count == 3);
    CHECK(!gate.commit(gateKey, 60, signature, false) && !gate.released(gateKey));
    CHECK(gate.observe(gateKey, 90, signature).count == 1);
    CHECK(gate.observe(gateKey, 120, changed).count == 1);
    CHECK(gate.observe(gateKey, 150, changed).count == 2);
    CHECK(gate.observe(gateKey, 150, "").count == 0);
    CHECK(gate.observe(gateKey, 180, signature).count == 1);
    CHECK(gate.observe(gateKey, 210, signature).count == 2);
    CHECK(gate.observe(gateKey, 260, signature).count == 1);
    CHECK(gate.observe(gateKey, 270, signature).count == 1);
    CHECK(gate.observe(gateKey, 300, signature).count == 2);
    CHECK(gate.observe(gateKey, 290, signature).count == 0);
    CHECK(gate.observe(gateKey, 330, signature).count == 1);
    CHECK(gate.observe(gateKey, 360, signature).count == 2);
    CHECK(gate.observe(gateKey, 390, signature).eligible);
    CHECK(!gate.commit(gateKey, 390, changed, true) && !gate.released(gateKey));
    CHECK(gate.observe(gateKey, 420, signature).count == 1);
    gate.observe(gateKey, 450, signature);
    CHECK(gate.observe(gateKey, 480, signature).eligible);
    CHECK(!gate.commit(gateKey, 450, signature, true));
    gate.observe(gateKey, 510, signature);
    gate.observe(gateKey, 540, signature);
    CHECK(gate.observe(gateKey, 570, signature).eligible);
    CHECK(gate.commit(gateKey, 570, signature, true) && gate.released(gateKey));
    decision = gate.observe(gateKey, 900, changed);
    CHECK(decision.released && decision.eligible);
    decision = gate.observe(gateKey, 930, "", false);
    CHECK(decision.released && !decision.eligible);

    const ZhangUserArExperimentKey otherReceiver{"run-a", "DYNG", "GPS"};
    const ZhangUserArExperimentKey otherRuntime{"run-b", "NICO", "GPS"};
    CHECK(!gate.released(otherReceiver) && !gate.released(otherRuntime));
    CHECK(gate.observe(otherReceiver, 0, signature).count == 1);
    CHECK(gate.observe(otherRuntime, 0, signature).count == 1);
    CHECK(gate.size() == 3);
    ZhangUserArFirstFeedbackGate independent;
    CHECK(independent.observe(gateKey, 0, signature).count == 1 && !independent.released(gateKey));
    gate.clearRuntime("run-a");
    CHECK(gate.size() == 1 && !gate.released(gateKey));
    CHECK(gate.observe(otherRuntime, 30, signature).count == 2);
    arcs.update(arcKey, 0, identity);
    arcs.update({"run-b", "NICO", "G20"}, 0, identity);
    arcs.clearRuntime("run-a");
    CHECK(arcs.size() == 1);
    CHECK(arcs.update(arcKey, 30, identity).count == 1);
    CHECK(arcs.update({"run-b", "NICO", "G20"}, 30, identity).count == 2);
    gate.clear();
    arcs.clear();
    CHECK(gate.size() == 0 && arcs.size() == 0);

    std::cout << "PASS user_ar_experiment_contract assertions=" << assertions << '\n';
    return EXIT_SUCCESS;
}
