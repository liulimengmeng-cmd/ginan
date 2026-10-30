// #pragma GCC optimize ("O0")

#include "common/linearCombo.hpp"
#include <cmath>
#include <initializer_list>
#include "common/acsQC.hpp"
#include "common/debug.hpp"
#include "common/navigation.hpp"
#include "common/observations.hpp"
#include "common/satStat.hpp"
#include "common/trace.hpp"

namespace
{
double phaseMetresOrMissing(double cycles, double wavelength)
{
    if (!std::isfinite(cycles) || !std::isfinite(wavelength) || wavelength <= 0) return 0;
    const double metres = cycles * wavelength;
    return std::isfinite(metres) ? metres : 0;
}

double codeMetresOrMissing(double metres)
{
    return std::isfinite(metres) && metres > 0 ? metres : 0;
}

bool validLcNumbers(const S_LC& lc)
{
    if (lc.lam_A <= 0 || lc.lam_B <= 0 || lc.lam_A == lc.lam_B ||
        lc.lam_WL == 0 || lc.lam_NL <= 0) return false;
    for (const double value : {lc.GF_Phas_m, lc.GF_Code_m,
            lc.WL_Phas_m, lc.WL_Phas_c, lc.WL_Code_m, lc.WL_Code_c,
            lc.NL_Phas_m, lc.NL_Phas_c, lc.NL_Code_m, lc.NL_Code_c,
            lc.IF_Phas_m, lc.IF_Code_m, lc.MW_m, lc.MW_c,
            lc.lam_A, lc.lam_B, lc.lam_WL, lc.lam_NL})
        if (!std::isfinite(value)) return false;
    return true;
}
}

/** Create combinations between specific observation values
 */
S_LC getLC(
    double  L_A_m,   ///< Phase measurement A (in meters)
    double  L_B_m,   ///< Phase measurement B (in meters)
    double  P_A_m,   ///< Code measurement A (in meters)
    double  P_B_m,   ///< Code measurement B (in meters)
    double  lamA,    ///< Wavelength A
    double  lamB,    ///< Wavelength B
    double* c1_out,  ///< Ionosphere free coefficient 1
    double* c2_out   ///< Ionosphere free coefficient 2
)
{
    S_LC lc = {};

    if (c1_out) *c1_out = 0;
    if (c2_out) *c2_out = 0;

    // A duplicate configured band is not a second physical frequency. In
    // particular, satFreqs() can leave its default third band equal to F5
    // when only F1/F5 are configured. Never form an IF/WL denominator of 0.
    for (const double value : {L_A_m, L_B_m, P_A_m, P_B_m, lamA, lamB})
        if (!std::isfinite(value)) return lc;
    if (lamA <= 0 || lamB <= 0 || lamA == lamB) return lc;

    if (P_A_m <= 0 || P_B_m <= 0 || L_A_m == 0 || L_B_m == 0)
    {
        // 		tracepde(lv, fppde, "PDE, code observation insufficient\n");
        return lc;
    }

    lc.lam_A = lamA;
    lc.lam_B = lamB;

    double L_A_c = L_A_m / lamA;
    double L_B_c = L_B_m / lamB;
    double P_A_c = P_A_m / lamA;
    double P_B_c = P_B_m / lamB;

    /* phase gf, wl, mw */
    double c1 = lamB * lamB / (lamB * lamB - lamA * lamA); /* IF */
    double c2 = lamA * lamA / (lamB * lamB - lamA * lamA);
    // 	double c3 = lamB / (lamB - lamA); /* WL */
    // 	double c4 = lamA / (lamB - lamA);

    // 	lc.GF_Phas_m =      L_A_m -      L_B_m;
    // 	lc.IF_Phas_m = c1 * L_A_m - c2 * L_B_m;
    // 	lc.WL_Phas_m = c3 * L_A_m - c4 * L_B_m;

    // 	c3 = lamB / (lamB + lamA); /* MW */
    // 	c4 = lamA / (lamB + lamA);
    //
    // 	lc.GF_Code_m =      P_A_m -      P_B_m;	/* geometry-free codes are independent from phase */
    // 	lc.IF_Code_m = c1 * P_A_m - c2 * P_B_m;
    // 	lc.NL_Code_m  = c3 * P_A_m + c4 * P_B_m;
    //
    // 	double lamw = lamA * lamB / (lamB - lamA);
    // 	lc.MW_c  = (lc.WL_Phas_m - lc.NL_Code_m) / lamw; /* cycle */

    lc.lam_WL = lamA * lamB / (lamB - lamA);
    lc.lam_NL = lamA * lamB / (lamB + lamA);

    lc.WL_Phas_c = L_A_c - L_B_c;
    lc.WL_Code_c = P_A_c - P_B_c;

    lc.NL_Phas_c = L_A_c + L_B_c;
    lc.NL_Code_c = P_A_c + P_B_c;

    lc.GF_Phas_m = L_A_m - L_B_m;
    lc.GF_Code_m = P_A_m - P_B_m;

    lc.IF_Phas_m = c1 * L_A_m - c2 * L_B_m;
    lc.IF_Code_m = c1 * P_A_m - c2 * P_B_m;

    lc.WL_Phas_m = lc.WL_Phas_c * lc.lam_WL;
    lc.WL_Code_m = lc.WL_Code_c * lc.lam_NL;

    lc.NL_Phas_m = lc.NL_Phas_c * lc.lam_WL;
    lc.NL_Code_m = lc.NL_Code_c * lc.lam_NL;

    lc.MW_m = lc.WL_Phas_m - lc.NL_Code_m;
    lc.MW_c = lc.MW_m / lc.lam_WL;

    // Also reject overflow in otherwise finite but nonphysical inputs.
    if (!std::isfinite(c1) || !std::isfinite(c2) || !validLcNumbers(lc)) return {};
    if (c1_out) *c1_out = c1;
    if (c2_out) *c2_out = c2;
    lc.valid = true;
    return lc;
}

/** Get combinations from pre-computed values, or return an empty value
 */
S_LC& getLC(
    lc_t&   lcBase,  ///< Linear combination base object
    E_FType fA,      ///< Frequency type A
    E_FType fB       ///< Frequency type B
)
{
    if (fA > fB)
        std::swap(fA, fB);

    // try to get existing LC from the observation's satStat object
    auto& lc = lcBase.lcMap[{fA, fB}];
    if (fA == fB || !validLcNumbers(lc)) lc = {};
    return lc;
}

/** Get/calculate linear combination values for an observation
 */
S_LC& getLC(
    GObs&   obs,     ///< Observation to compute values form
    lc_t&   lcBase,  ///< Linear combination base object
    E_FType fA,      ///< Frequency type A
    E_FType fB       ///< Frequency type B
)
{
    // try to get existing LC from the observation's satStat object
    S_LC& lc = getLC(lcBase, fA, fB);

    if (fA == fB)
    {
        lc = {};
        return lc;
    }

    if (lc.valid)
    {
        return lc;
    }

    // make a new linear combination from the observation
    lcBase.time = obs.time;
    for (E_FType f : {fA, fB})
    {
        // These are derived LC caches, not the original observations. A
        // nonfinite derived value denotes missing data, never a valid zero
        // observation; scalar getLC still requires both nonzero phases/codes.
        if (!std::isfinite(lcBase.L_m[f])) lcBase.L_m[f] = 0;
        lcBase.P[f] = codeMetresOrMissing(lcBase.P[f]);
        if (lcBase.L_m[f] == 0)
        {
            // no L measurement, try to get from observation
            lcBase.L_m[f] = phaseMetresOrMissing(obs.sigs[f].L, obs.satNav_ptr->lamMap[f]);
            lcBase.P[f]   = codeMetresOrMissing(obs.sigs[f].P);
        }
        if (lcBase.L_m[f] == 0)
        {
            // still no measurement, give up
            return lc;
        }
    }

    double L_A  = lcBase.L_m[fA];
    double L_B  = lcBase.L_m[fB];
    double P_A  = lcBase.P[fA];
    double P_B  = lcBase.P[fB];
    double lamA = obs.satNav_ptr->lamMap[fA];
    double lamB = obs.satNav_ptr->lamMap[fB];

    lc = getLC(L_A, L_B, P_A, P_B, lamA, lamB, nullptr, nullptr);
    if (!lc.valid) return lc;

    // special cases
    if (fB == F5 && (obs.Sat.sys == E_Sys::GAL || obs.Sat.sys == E_Sys::BDS))
        lc.MW_c *= -1; /* cycle */

    if (fA == F1 && fB == F2)
    {
        lcBase.mp[F1] = P_A - L_A - 2.0 * lamA * lamA / (lamB * lamB - lamA * lamA) * (L_A - L_B);
        lcBase.mp[F2] = P_B - L_B - 2.0 * lamB * lamB / (lamB * lamB - lamA * lamA) * (L_A - L_B);
    }
    else if (fB == F5 && lcBase.mp[F5] == 0)
    {
        lcBase.mp[F5] = P_B - L_B - 2.0 * lamB * lamB / (lamB * lamB - lamA * lamA) * (L_A - L_B);
    }

    return lc;
}

/** Prepare a base object for linear combinations using observation data
 */
void lcPrepareBase(
    GObs& obs,    ///< Observation data to use
    lc_t& lcBase  ///< Linear combination base object to prepare
)
{
    lcBase.time = obs.time;
    lcBase.Sat  = obs.Sat;

    for (auto& [ft, sig] : obs.sigs)
    {
        // populate variables for later use.
        lcBase.L_m[ft] = phaseMetresOrMissing(sig.L, obs.satNav_ptr->lamMap[ft]);
        lcBase.P[ft]   = codeMetresOrMissing(sig.P);
    }
}

/** Function to prepare some predefined linear combinations from an observation
 */
void obs2lc(
    Trace& trace,  ///< Trace to output to
    GObs&  obs,    ///< Observation to prepare combinations for
    lc_t&  lcBase  ///< Linear combination base object to use
)
{
    E_Sys sys = obs.Sat.sys;

    E_FType frq1;
    E_FType frq2;
    E_FType frq3;

    bool pass = satFreqs(obs.Sat.sys, frq1, frq2, frq3);

    if (pass == false)
        return;

    char strprefix[64];
    snprintf(
        strprefix,
        sizeof(strprefix),
        "%3s sat=%4s",
        obs.time.to_string().c_str(),
        obs.Sat.id().c_str()
    );

    lcPrepareBase(obs, lcBase);

    // iterate pairwise over the frequencies.
    S_LC& lc12 = getLC(obs, lcBase, frq1, frq2);
    S_LC& lc15 = getLC(obs, lcBase, frq1, frq3);
    S_LC& lc25 = getLC(obs, lcBase, frq2, frq3);

    string frq1Str  = enum_to_string(frq1);
    string frq2Str  = enum_to_string(frq2);
    string frq3Str  = enum_to_string(frq3);
    string frq12Str = frq1Str + frq2Str;
    string frq13Str = frq1Str + frq3Str;
    string frq23Str = frq2Str + frq3Str;

    tracepdeex(
        3,
        trace,
        "%s zd L -- %-3s =%14.4f %-3s =%14.4f %-3s =%14.4f\n",
        strprefix,
        frq1Str.c_str(),
        lcBase.L_m[frq1],
        frq2Str.c_str(),
        lcBase.L_m[frq2],
        frq3Str.c_str(),
        lcBase.L_m[frq3]
    );
    tracepdeex(
        3,
        trace,
        "%s zd P -- %-3s =%14.4f %-3s =%14.4f %-3s =%14.4f\n",
        strprefix,
        frq1Str.c_str(),
        lcBase.P[frq1],
        frq2Str.c_str(),
        lcBase.P[frq2],
        frq3Str.c_str(),
        lcBase.P[frq3]
    );
    tracepdeex(
        3,
        trace,
        "%s mp P -- %-3s =%14.4f %-3s =%14.4f %-3s =%14.4f\n",
        strprefix,
        frq1Str.c_str(),
        lcBase.mp[frq1],
        frq2Str.c_str(),
        lcBase.mp[frq2],
        frq3Str.c_str(),
        lcBase.mp[frq3]
    );
    tracepdeex(
        3,
        trace,
        "%s gf L -- %-6s=%14.4f %-6s=%14.4f %-6s=%14.4f\n",
        strprefix,
        frq12Str.c_str(),
        lc12.GF_Phas_m,
        frq13Str.c_str(),
        lc15.GF_Phas_m,
        frq23Str.c_str(),
        lc25.GF_Phas_m
    );
    tracepdeex(
        3,
        trace,
        "%s gf P -- %-6s=%14.4f %-6s=%14.4f %-6s=%14.4f\n",
        strprefix,
        frq12Str.c_str(),
        lc12.GF_Code_m,
        frq13Str.c_str(),
        lc15.GF_Code_m,
        frq23Str.c_str(),
        lc25.GF_Code_m
    );
    tracepdeex(
        3,
        trace,
        "%s mw L -- %-6s=%14.4f %-6s=%14.4f %-6s=%14.4f\n",
        strprefix,
        frq12Str.c_str(),
        lc12.MW_c,
        frq13Str.c_str(),
        lc15.MW_c,
        frq23Str.c_str(),
        lc25.MW_c
    );
    tracepdeex(
        3,
        trace,
        "%s wl L -- %-6s=%14.4f %-6s=%14.4f %-6s=%14.4f\n",
        strprefix,
        frq12Str.c_str(),
        lc12.WL_Phas_m,
        frq13Str.c_str(),
        lc15.WL_Phas_m,
        frq23Str.c_str(),
        lc25.WL_Phas_m
    );
    tracepdeex(
        3,
        trace,
        "%s if L -- %-6s=%14.4f %-6s=%14.4f %-6s=%14.4f\n",
        strprefix,
        frq12Str.c_str(),
        lc12.IF_Phas_m,
        frq13Str.c_str(),
        lc15.IF_Phas_m,
        frq23Str.c_str(),
        lc25.IF_Phas_m
    );
    tracepdeex(
        3,
        trace,
        "%s if P -- %-6s=%14.4f %-6s=%14.4f %-6s=%14.4f\n",
        strprefix,
        frq12Str.c_str(),
        lc12.IF_Code_m,
        frq13Str.c_str(),
        lc15.IF_Code_m,
        frq23Str.c_str(),
        lc25.IF_Code_m
    );

    traceJson(
        5,
        trace,
        obs.time,
        {{"data", "linearCombos"}, {"Sat", obs.Sat.id()}},
        {{frq1Str, lcBase.L_m[frq1]}, {frq2Str, lcBase.L_m[frq2]}}
    );
}

/** Function to prepare some predefined linear combinations from a list of observations
 */
void obs2lcs(
    Trace&   trace,   ///< Trace to output to
    ObsList& obsList  ///< List of bservation to prepare combinations for
)
{
    int lv = 3;

    if (obsList.empty())
    {
        return;
    }

    tracepdeex(
        lv,
        trace,
        "\n   *-------- PDE form LC %s             --------*\n",
        obsList.front()->time.to_string().c_str()
    );

    for (auto& obs : only<GObs>(obsList))
    {
        if (obs.exclude)
        {
            continue;
        }

        lc_t& lc = obs.satStat_ptr->lc_new;
        obs2lc(trace, obs, lc);
    }
}
