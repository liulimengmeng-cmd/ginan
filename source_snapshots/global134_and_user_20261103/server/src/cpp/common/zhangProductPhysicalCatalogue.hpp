#pragma once

/** Product nodes belong to the connected physical graph, not its independent
 * chord coordinate list. This helper grants no integer or arc authority. */
template<class Basis>
inline auto zhangPhysicalProductSatelliteCatalogue(const Basis& basis)
    -> decltype(basis.satellites)
{
    if (!basis.connected || !basis.receivers.contains(basis.rootReceiver)) return {};
    return basis.satellites;
}
