#include "common/zhangRealtimeNavigation.hpp"
#include "common/acsConfig.hpp"
#include "common/navigation.hpp"
#include "common/zhangCheckpoint.hpp"
#include <boost/json.hpp>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace zhangRealtimeNavigation
{
namespace j = boost::json;
using EphMap = decltype(nav.ephMap);
using SsrMap = std::map<SatSys, SSRMaps>;
namespace
{
GTime currentEpoch;
bool active = false, serverBackup = false;
bool selectingServerInputs = false;
EphMap previousEph;
SsrMap previousSsr;
ERP previousErp;
Navigation frozenNavigation;
std::string frozenBytes, frozenHash;

[[noreturn]] void reject(const std::string& reason)
{ throw std::runtime_error("REALTIME_NAV:" + reason); }

std::string setting(const char* key)
{
    const char* value = std::getenv(key);
    if (!value || !*value) reject(std::string("SETTING_REQUIRED:") + key);
    return value;
}

std::string timeText(GTime time)
{
    std::ostringstream text;
    text << std::setprecision(std::numeric_limits<long double>::max_digits10) << time.bigTime;
    return text.str();
}

GTime readTime(const j::value& value)
{
    if (!value.is_string()) reject("TIME_TYPE");
    const std::string text(value.as_string());
    std::size_t used = 0;
    GTime result;
    result.bigTime = std::stold(text, &used);
    if (used != text.size() || !std::isfinite(result.bigTime)) reject("TIME_VALUE");
    return result;
}

double real(const j::value& value)
{
    const double result = value.to_number<double>();
    if (!std::isfinite(result)) reject("NONFINITE_VALUE");
    return result;
}

j::array erpValues(const ERPValues& p)
{
    return {timeText(p.time), p.xp, p.yp, p.ut1Utc, p.lod, p.xpr, p.ypr,
        p.xpSigma, p.ypSigma, p.xprSigma, p.yprSigma, p.ut1UtcSigma, p.lodSigma,
        p.isPredicted, p.isFiltered};
}

ERPValues readErpValues(const j::array& p)
{
    if (p.size() != 15) reject("ERP_WIDTH");
    ERPValues result;
    result.time = readTime(p[0]);
    double* fields[] = {&result.xp, &result.yp, &result.ut1Utc, &result.lod, &result.xpr, &result.ypr,
        &result.xpSigma, &result.ypSigma, &result.xprSigma, &result.yprSigma,
        &result.ut1UtcSigma, &result.lodSigma};
    for (std::size_t i = 0; i < 12; ++i) *fields[i] = real(p[i + 1]);
    result.isPredicted = p[13].as_bool(); result.isFiltered = p[14].as_bool();
    return result;
}

j::array metadata(const SSRMeta& m)
{
    // receivedTime in the upstream decoder is the SSR MESSAGE epoch. It is
    // deliberately NOT labelled a wall-clock arrival timestamp here.
    return {m.epochTime1s, timeText(m.receivedTime), m.updateIntIndex,
        m.multipleMessage, m.referenceDatum, m.provider, m.solution, m.numSats};
}

j::array attitude(const AttStatus& p)
{
    j::array result{timeText(p.startTime), p.startSign, p.startYaw, p.startYawRate,
        timeText(p.excludeTime), p.nominalYaw, p.modelYaw, timeText(p.modelYawTime), p.modelYawValid};
    for (const auto* axis : {&p.eXBody, &p.eYBody, &p.eZBody, &p.eXAnt, &p.eYAnt, &p.eZAnt})
        for (int i = 0; i < 3; ++i) result.push_back(p.modelYawValid ? (*axis)(i) : 0.0);
    return result;
}

AttStatus readAttitude(const j::array& p)
{
    if (p.size() != 27) reject("ATTITUDE_WIDTH");
    AttStatus r;
    r.startTime = readTime(p[0]); r.startSign = real(p[1]); r.startYaw = real(p[2]);
    r.startYawRate = real(p[3]); r.excludeTime = readTime(p[4]);
    r.nominalYaw = real(p[5]); r.modelYaw = real(p[6]); r.modelYawTime = readTime(p[7]);
    r.modelYawValid = p[8].as_bool();
    int index = 9;
    for (auto* axis : {&r.eXBody, &r.eYBody, &r.eZBody, &r.eXAnt, &r.eYAnt, &r.eZAnt})
        for (int i = 0; i < 3; ++i) (*axis)(i) = real(p[index++]);
    return r;
}

SSRMeta readMetadata(const j::array& m)
{
    if (m.size() != 8) reject("SSR_METADATA_WIDTH");
    SSRMeta r;
    r.epochTime1s = m[0].to_number<int>(); r.receivedTime = readTime(m[1]);
    r.updateIntIndex = m[2].to_number<int>(); r.multipleMessage = m[3].to_number<int>();
    r.referenceDatum = m[4].to_number<unsigned>(); r.provider = m[5].to_number<unsigned>();
    r.solution = m[6].to_number<unsigned>(); r.numSats = m[7].to_number<unsigned>();
    return r;
}

j::array ephemeris(const Eph& e)
{
    // Explicit mathematical input DTO for LEGACY_GEC; no opaque object bytes,
    // pointers, unused/uninitialised debug fields or receiver filter state.
    return {e.iode, e.iodc, static_cast<int>(e.svh), e.code, e.flag,
        timeText(e.toe), timeText(e.toc), timeText(e.ttm),
        e.A, e.e, e.i0, e.OMG0, e.omg, e.M0, e.deln, e.OMGd, e.idot,
        e.crc, e.crs, e.cuc, e.cus, e.cic, e.cis, e.toes,
        e.f0, e.f1, e.f2, e.tgd[0], e.tgd[1], e.tgd[2], e.tgd[3],
        e.ura[0], e.ura[1], e.ura[2], e.ura[3]};
}

Eph readEphemeris(const j::array& values, SatSys satellite, E_NavMsgType type)
{
    if (values.size() != 35) reject("EPHEMERIS_WIDTH");
    Eph e{};
    e.Sat = satellite; e.type = type;
    e.iode = values[0].to_number<int>(); e.iodc = values[1].to_number<int>();
    e.svh = static_cast<E_Svh>(values[2].to_number<int>());
    e.code = values[3].to_number<int>(); e.flag = values[4].to_number<int>();
    e.toe = readTime(values[5]); e.toc = readTime(values[6]); e.ttm = readTime(values[7]);
    double* fields[] = {&e.A, &e.e, &e.i0, &e.OMG0, &e.omg, &e.M0, &e.deln, &e.OMGd, &e.idot,
        &e.crc, &e.crs, &e.cuc, &e.cus, &e.cic, &e.cis, &e.toes,
        &e.f0, &e.f1, &e.f2, &e.tgd[0], &e.tgd[1], &e.tgd[2], &e.tgd[3],
        &e.ura[0], &e.ura[1], &e.ura[2], &e.ura[3]};
    for (std::size_t i = 0; i < 27; ++i) *fields[i] = real(values[i + 8]);
    if (!(e.A > 1e7) || e.e < 0 || e.e >= 1 || e.iode < 0) reject("EPHEMERIS_VALUE");
    return e;
}

bool sameProvider(const SSRMeta& a, const SSRMeta& b)
{ return a.provider == b.provider && a.solution == b.solution; } // Clock headers carry no CRS bit.

void apply(EphMap eph, SsrMap ssr, ERP erp)
{
    nav.ephMap.swap(eph);
    // Preserve SatNav object addresses, antenna/attitude state and receiver
    // pointers. A missing satellite gets NO old SSR fallback.
    for (auto& [sat, data] : nav.satNavMap)
    {
        data.receivedSSR = {};
        data.satPos0.posTime = GTime::noTime();
    }
    for (auto& [sat, corrections] : ssr) nav.satNavMap[sat].receivedSSR = std::move(corrections);
    nav.erp = std::move(erp);
}

void decodeAndApply(const std::string& bytes, const GTime& time)
{
    if (bytes.size() > 16 * 1024 * 1024) reject("SNAPSHOT_SIZE_LIMIT");
    const auto document = j::parse(bytes).as_object();
    if (document.at("schema").as_string() != "REALTIME_NAV_LEGACY_GEC_V1" ||
        document.at("epoch").as_string() != time.to_string(0) ||
        document.at("generation").as_string() != setting("ZHANG_REALTIME_GENERATION") ||
        document.at("model_sha256").as_string() != setting("ZHANG_REALTIME_MODEL_SHA256") ||
        document.at("mode").as_string() != setting("ZHANG_REALTIME_NAV_MODE")) reject("SNAPSHOT_IDENTITY");
    if (document.at("ssr_offset").to_number<int>() != static_cast<int>(acsConfig.ssr_input_antenna_offset))
        reject("SSR_REFERENCE_POINT_MISMATCH");
    if (document.at("igs_bds_epoch_policy").as_string() != setting("ZHANG_REALTIME_IGS_BDS_TIME"))
        reject("BDS_WIRE_TIME_POLICY_MISMATCH");
    EphMap eph;
    SsrMap ssr;
    std::map<SatSys, AttStatus> attitudes;
    const auto& records = document.at("satellites").as_array();
    if (records.size() > 192) reject("SATELLITE_LIMIT");
    for (const auto& record : records)
    {
        const auto& r = record.as_object();
        const SatSys satellite(std::string(r.at("satellite").as_string()).c_str());
        if (satellite.sys != E_Sys::GPS && satellite.sys != E_Sys::GAL && satellite.sys != E_Sys::BDS)
            reject("UNSUPPORTED_SYSTEM");
        const auto type = static_cast<E_NavMsgType>(r.at("nav_type").to_number<int>());
        if (acsConfig.used_nav_types.at(satellite.sys) != type) reject("NAV_TYPE_MISMATCH");
        const auto key = readTime(r.at("ephemeris_key"));
        if (eph.contains(satellite)) reject("DUPLICATE_SATELLITE");
        auto value = readEphemeris(r.at("ephemeris").as_array(), satellite, type);
        if (!eph[satellite][type].emplace(key, std::move(value)).second) reject("DUPLICATE_SATELLITE");
        if (r.contains("orbit"))
        {
            const auto& p = r.at("orbit").as_array();
            const auto& c = r.at("clock").as_array();
            if (p.size() != 13 || c.size() != 8) reject("SSR_RECORD_WIDTH");
            SSREph pos; SSRClk clk;
            const GTime posKey = readTime(p[0]), clkKey = readTime(c[0]);
            pos.t0 = readTime(p[1]); pos.udi = real(p[2]); pos.iod = p[3].to_number<int>();
            pos.iode = p[4].to_number<int>(); pos.iodcrc = p[5].to_number<int>();
            for (int i = 0; i < 3; ++i) { pos.deph(i) = real(p[6 + i]); pos.ddeph(i) = real(p[9 + i]); }
            pos.ssrMeta = readMetadata(p[12].as_array());
            clk.t0 = readTime(c[1]); clk.udi = real(c[2]); clk.iod = c[3].to_number<int>();
            for (int i = 0; i < 3; ++i) clk.dclk[i] = real(c[4 + i]);
            clk.ssrMeta = readMetadata(c[7].as_array());
            if (posKey > time || clkKey > time || pos.iod != clk.iod ||
                !sameProvider(pos.ssrMeta, clk.ssrMeta) || pos.udi <= 0 || clk.udi <= 0 ||
                pos.ssrMeta.referenceDatum != 0 ||
                (satellite.sys == E_Sys::GAL ?
                    (pos.iode & 255) != (eph.at(satellite).at(type).begin()->second.iode & 255) :
                    pos.iode != eph.at(satellite).at(type).begin()->second.iode)) reject("SSR_PAIR_IDENTITY");
            ssr[satellite].ssrEph_map.emplace(posKey, pos);
            ssr[satellite].ssrClk_map.emplace(clkKey, clk);
        }
        attitudes.emplace(satellite, readAttitude(r.at("attitude").as_array()));
    }
    ERP erp;
    erp.filterValues = readErpValues(document.at("erp_filter").as_array());
    const auto& maps = document.at("erp_maps").as_array();
    if (maps.size() > 8) reject("ERP_MAP_LIMIT");
    for (const auto& history : maps)
    {
        std::map<GTime, ERPValues> map;
        if (history.as_array().size() > 512) reject("ERP_HISTORY_LIMIT");
        for (const auto& entry : history.as_array())
        {
            const auto& pair = entry.as_array();
            if (pair.size() != 2 || !map.emplace(readTime(pair[0]), readErpValues(pair[1].as_array())).second)
                reject("ERP_DUPLICATE");
        }
        erp.erpMaps.push_back(std::move(map));
    }
    if (document.at("leap_seconds").to_number<int>() != static_cast<int>(leapSeconds(time)))
        reject("LEAP_SECONDS_MISMATCH");
    frozenNavigation.ephMap = eph;
    frozenNavigation.satNavMap.clear();
    for (const auto& [satellite, corrections] : ssr)
        frozenNavigation.satNavMap[satellite].receivedSSR = corrections;
    frozenNavigation.erp = erp;
    apply(std::move(eph), std::move(ssr), std::move(erp));
    for (auto& [satellite, state] : attitudes) nav.satNavMap[satellite].attStatus = state;
    currentEpoch = time; active = true;
}
} // namespace

bool enabled()
{
    const char* value = std::getenv("ZHANG_REALTIME_NAV_MODE");
    if (!value || std::string(value) == "OFF") return false;
    if (std::string(value) != "FROZEN_SSR_GEC" && std::string(value) != "FROZEN_BROADCAST_GEC")
        reject("NAV_MODE_INVALID");
    return true;
}

GTime selectionEpoch(GTime evaluationTime)
{
    if (!enabled()) return evaluationTime;
    if (!active) reject("CONTEXT_NOT_PREPARED");
    if (std::abs((evaluationTime - currentEpoch).to_double()) > 2)
        reject("EVALUATION_OUTSIDE_FROZEN_WINDOW");
    return currentEpoch;
}

Navigation& orbitalContext(Navigation& live)
{
    if (!enabled()) return live;
    if (selectingServerInputs) return live;
    if (!active) reject("ORBITAL_CONTEXT_NOT_PREPARED");
    return frozenNavigation;
}

void beginServer(const GTime& time)
{
    if (!enabled()) return;
    if (serverBackup) reject("NESTED_SERVER_EPOCH");
    selectingServerInputs = true;
    struct EndSelection { ~EndSelection() { selectingServerInputs = false; } } endSelection;
    j::object document;
    document["schema"] = "REALTIME_NAV_LEGACY_GEC_V1";
    document["epoch"] = time.to_string(0);
    document["generation"] = setting("ZHANG_REALTIME_GENERATION");
    document["model_sha256"] = setting("ZHANG_REALTIME_MODEL_SHA256");
    document["mode"] = setting("ZHANG_REALTIME_NAV_MODE");
    document["ssr_offset"] = static_cast<int>(acsConfig.ssr_input_antenna_offset);
    document["leap_seconds"] = static_cast<int>(leapSeconds(time));
    document["evaluation_window_seconds"] = 2;
    document["satellite_reference_point"] = "COM";
    document["terrestrial_frame"] = setting("ZHANG_REALTIME_NAV_MODE") == "FROZEN_SSR_GEC"
        ? "CAS01 provider frame; no frame fit or conversion" : "broadcast native frame; no frame fit or conversion";
    document["external_bias_policy"] = "SSR code/phase bias and VTEC messages excluded by relay";
    document["high_rate_clock_policy"] = "not consumed by this frozen context";
    document["ura_policy"] = "native default, no received SSR URA used";
    const auto bdsPolicy = setting("ZHANG_REALTIME_IGS_BDS_TIME");
    if (bdsPolicy != "IGS_GPST" && bdsPolicy != "CAS01_BDT_COMPAT") reject("BDS_WIRE_TIME_POLICY_INVALID");
    document["igs_bds_epoch_policy"] = bdsPolicy;
    document["frozen_at_unix_seconds"] = std::chrono::duration<double>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    j::array records;
    for (const auto& [sat, types] : nav.ephMap)
    {
        if (sat.sys != E_Sys::GPS && sat.sys != E_Sys::GAL && sat.sys != E_Sys::BDS) continue;
        const auto enabledSystem = acsConfig.process_sys.find(sat.sys);
        if (enabledSystem == acsConfig.process_sys.end() || !enabledSystem->second) continue;
        const auto type = acsConfig.used_nav_types.at(sat.sys);
        const auto foundType = types.find(type);
        if (foundType == types.end() || foundType->second.empty()) continue;
        const SSREph* orbit = nullptr; const SSRClk* clock = nullptr;
        GTime orbitKey, clockKey;
        int iode = ANY_IODE;
        if (setting("ZHANG_REALTIME_NAV_MODE") == "FROZEN_SSR_GEC")
        {
            const auto found = nav.satNavMap.find(sat);
            if (found == nav.satNavMap.end()) continue;
            const auto& inputs = found->second.receivedSSR;
            for (auto p = inputs.ssrEph_map.lower_bound(time); p != inputs.ssrEph_map.end() && !orbit; ++p)
                for (auto c = inputs.ssrClk_map.lower_bound(time); c != inputs.ssrClk_map.end(); ++c)
                {
                    if (p->second.iod != c->second.iod || !sameProvider(p->second.ssrMeta, c->second.ssrMeta)) continue;
                    // The chosen polynomials must cover the full evaluation
                    // window, including finite-difference velocity queries.
                    const auto covers = [&](const auto& value)
                    { return value.udi > 0 && std::abs((time - value.t0).to_double()) + 2 <=
                        value.udi * acsConfig.validity_interval_factor; };
                    if (!covers(p->second) || !covers(c->second)) continue;
                    orbit = &p->second; clock = &c->second;
                    orbitKey = p->first; clockKey = c->first; iode = orbit->iode;
                    break;
                }
            if (!orbit) continue;
        }
        const auto* selected = seleph<Eph>(nullStream, time, sat, type, iode, nav);
        if (!selected) continue;
        const auto key = std::find_if(foundType->second.begin(), foundType->second.end(),
            [&](const auto& entry) { return &entry.second == selected; });
        if (key == foundType->second.end()) reject("SELECTED_EPHEMERIS_NOT_IN_MAP");
        j::object r;
        r["satellite"] = sat.id(); r["nav_type"] = static_cast<int>(type);
        r["ephemeris_key"] = timeText(key->first); r["ephemeris"] = ephemeris(*selected);
        r["attitude"] = attitude(nav.satNavMap[sat].attStatus);
        if (orbit)
        {
            r["orbit"] = j::array{timeText(orbitKey), timeText(orbit->t0), orbit->udi, orbit->iod,
                orbit->iode, orbit->iodcrc, orbit->deph(0), orbit->deph(1), orbit->deph(2),
                orbit->ddeph(0), orbit->ddeph(1), orbit->ddeph(2), metadata(orbit->ssrMeta)};
            r["clock"] = j::array{timeText(clockKey), timeText(clock->t0), clock->udi, clock->iod,
                clock->dclk[0], clock->dclk[1], clock->dclk[2], metadata(clock->ssrMeta)};
        }
        records.push_back(std::move(r));
    }
    j::object counts;
    for (const auto& record : records)
    {
        const std::string satellite(record.as_object().at("satellite").as_string());
        const std::string system(1, satellite.front());
        if (!counts.contains(system)) counts[system] = 0;
        counts[system] = counts[system].as_int64() + 1;
    }
    std::size_t orbitRecords=0, clockRecords=0;
    for (const auto& [satellite, data] : nav.satNavMap)
    {
        orbitRecords += data.receivedSSR.ssrEph_map.size();
        clockRecords += data.receivedSSR.ssrClk_map.size();
    }
    std::cerr << "ZHANG_REALTIME_NAV_FREEZE epoch=" << time.to_string(0)
        << " broadcast_satellites=" << nav.ephMap.size()
        << " ssr_orbits=" << orbitRecords << " ssr_clocks=" << clockRecords
        << " selected=" << j::serialize(counts) << std::endl;
    document["satellites"] = std::move(records);
    document["erp_filter"] = erpValues(nav.erp.filterValues);
    j::array maps;
    for (const auto& history : nav.erp.erpMaps)
    {
        j::array entries;
        for (const auto& [key, p] : history) entries.push_back(j::array{timeText(key), erpValues(p)});
        maps.push_back(std::move(entries));
    }
    document["erp_maps"] = std::move(maps);
    frozenBytes = j::serialize(document) + '\n';
    frozenHash = zhangCheckpointSha256(frozenBytes);
    previousEph = nav.ephMap; previousErp = nav.erp;
    previousSsr.clear();
    for (const auto& [sat, data] : nav.satNavMap) previousSsr.emplace(sat, data.receivedSSR);
    serverBackup = true;
    decodeAndApply(frozenBytes, time);
}

void endServer()
{
    if (!serverBackup) return;
    apply(std::move(previousEph), std::move(previousSsr), std::move(previousErp));
    serverBackup = false; active = false;
}

void writeSnapshot(const std::filesystem::path& directory, const std::string& epoch)
{
    if (!enabled()) return;
    if (!active || !serverBackup || epoch != currentEpoch.to_string(0) ||
        !directory.filename().string().ends_with(".pending")) reject("WRITER_CONTEXT");
    std::ofstream output(directory / "navigation.json", std::ios::binary);
    output << frozenBytes; output.flush();
    if (!output) reject("SNAPSHOT_WRITE_FAILED");
}

void importSnapshot(const std::filesystem::path& path, const GTime& time)
{
    if (!enabled()) reject("USER_NAV_MODE_REQUIRED");
    if (!std::filesystem::is_regular_file(path) || std::filesystem::file_size(path) > 16 * 1024 * 1024)
        reject("USER_SNAPSHOT_UNAVAILABLE_OR_OVERSIZE");
    std::ifstream input(path, std::ios::binary);
    std::ostringstream text; text << input.rdbuf();
    if (!input || input.bad()) reject("USER_SNAPSHOT_READ_FAILED");
    const auto bytes = text.str();
    decodeAndApply(bytes, time);
    frozenHash = zhangCheckpointSha256(bytes);
}

void missingEpoch(const GTime& time)
{
    if (!enabled()) return;
    // No orbital measurement can be formed from a stale snapshot. Keeping ERP
    // alone permits the ordinary prediction/output path to advance its clock.
    apply({}, {}, nav.erp);
    frozenNavigation.ephMap.clear(); frozenNavigation.satNavMap.clear();
    currentEpoch = time; active = true; frozenHash.clear();
}

const std::string& snapshotHash() { return frozenHash; }
} // namespace zhangRealtimeNavigation
