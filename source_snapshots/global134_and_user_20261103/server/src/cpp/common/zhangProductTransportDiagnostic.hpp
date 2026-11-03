#pragma once
#include <string>

// Failure classification is descriptive only. It cannot authorize a product,
// change a segment, synthesize an old endpoint, or release a safety gate.
inline std::string zhangProductTransportFailureClass(const std::string& reason)
{
    auto has = [&](const char* token) { return reason.find(token) != std::string::npos; };
    if (has("NO_OLD_ENDPOINT")) return "NO_OLD_ENDPOINT";
    if (has("BUDGET") || has("UNFINISHED") || has("DEADLINE")) return "BUDGET_INCOMPLETE";
    if (has("RATIO")) return "RATIO_REJECTED";
    if (has("FLOAT_TARGET") || has("FLOAT_ROW") || has("MISSING_POSTERIOR_COLUMN") || has("RETAINED_FLOAT"))
        return "OLD_FLOAT_TARGET_UNAVAILABLE";
    if (has("RECEIPT_MISSING") || has("UNCOMPOSED") || has("TRANSPORT_CHAIN_BROKEN"))
        return "IDENTITY_RECEIPT_MISSING";
    if (has("NOT_ATTEMPTED") || has("NO_CURRENT_EPOCH_EVIDENCE") || has("NOT_EVALUATED"))
        return "SEARCH_NOT_ATTEMPTED";
    if (has("UNCERTIFIED") || has("NOT_ENTAILED") || has("NO_CURRENT") || has("NOT_PROVEN"))
        return "CURRENT_INTEGER_UNCERTIFIED";
    return "CERTIFICATE_NOT_APPLICABLE";
}
