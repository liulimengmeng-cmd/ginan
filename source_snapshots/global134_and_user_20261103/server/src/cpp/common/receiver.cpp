#include "common/receiver.hpp"
#include <set>
#include "common/acsConfig.hpp"
#include "common/streamParser.hpp"
#include "common/streamRinex.hpp"
#include "common/streamRtcm.hpp"

ReceiverMap receiverMap;

void extractReceiverMetadata(Receiver& rec, Parser& parser, ObsList* obsList)
{
    auto& recOpts = acsConfig.getRecOpts(rec.id);
    rec.metadata.setPriority(recOpts.meta_priority);

    string parserType = parser.parserType();

    if (parserType == "RinexParser")
    {
        auto& rinexParser = static_cast<RinexParser&>(parser);
        rec.metadata.ingestRinex(rinexParser.rnxRec);
    }
    else if (parserType == "RtcmParser")
    {
        auto& rtcmParser = static_cast<RtcmParser&>(parser);

        if (auto* info_ptr = selectRtcmStationInfoForMetadata(
                rtcmParser.stationInfoMap,
                rtcmParser.lastReferenceStationId
            ))
        {
            rec.metadata.ingestRtcm(*info_ptr);
        }
    }

    (void)obsList;

    syncReceiverMetadata(rec);
}

void extractTrackedSignals(Receiver& rec, Parser& parser, ObsList* obsList)
{
    // Extract tracked signals from RINEX header or accumulate from RTCM observations
    // This enables receiver-specific signal tracking mode detection
    string parserType = parser.parserType();

    if (parserType == "RinexParser" && rec.trackedSignals.empty())
    {
        auto& rinexParser = static_cast<RinexParser&>(parser);

        for (auto& [sys, codeTypeMap] : rinexParser.sysCodeTypes)
        {
            vector<E_ObsCode> signals;
            for (auto& [idx, codeType] : codeTypeMap)
            {
                if (codeType.code != E_ObsCode::NONE)
                {
                    signals.push_back(codeType.code);
                }
            }
            // RINEX2 has an unqualified header type catalogue, stored under
            // NONE, whose code3 fields are empty.  That is absence of tracking
            // metadata, not a tracked constellation.  Do not persist a fake
            // NONE system or infer tracking suffixes from legacy L1/L2 names.
            // Nonempty malformed system records remain visible to validation.
            if (sys == E_Sys::NONE && signals.empty())
                continue;
            rec.trackedSignals[sys] = signals;
        }
    }
    else if (parserType == "RtcmParser" && obsList != nullptr)
    {
        // For RTCM streams, accumulate tracked signals from observations
        // This builds up the signal list dynamically as we receive data

        // Extract signals from the provided observation list
        for (auto& obs : only<GObs>(*obsList))
        {
            E_Sys          sys = obs.Sat.sys;
            set<E_ObsCode> signalSet;

            // If we already have signals for this system, start with those
            if (rec.trackedSignals.count(sys))
            {
                signalSet.insert(rec.trackedSignals[sys].begin(), rec.trackedSignals[sys].end());
            }

            // Add signals from this observation
            for (auto& [ftype, sigsList] : obs.sigsLists)
            {
                for (auto& sig : sigsList)
                {
                    if (sig.code != E_ObsCode::NONE)
                    {
                        signalSet.insert(sig.code);
                    }
                }
            }

            // Update the tracked signals with the expanded set
            rec.trackedSignals[sys] = vector<E_ObsCode>(signalSet.begin(), signalSet.end());
        }
    }
    else
    {
        // Signal extraction not supported for other parser types (UBX, Custom, etc.)
    }
}
