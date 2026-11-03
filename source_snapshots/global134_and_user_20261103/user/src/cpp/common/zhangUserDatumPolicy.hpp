#pragma once
#include <cstdlib>
#include <stdexcept>
#include <string>

// The validity of a product map is computed before this experimental choice.
// Both policies use the same UNKNOWN recovery and physical tracking policy.
enum class ZhangUserDatumPolicy { TRANSPORT_KNOWN, RESET_KNOWN_LOCAL };
inline bool zhangUserDatumPolicyRequested()
{const char* p=std::getenv("ZHANG_USER_DATUM_POLICY");return p && *p;}
inline ZhangUserDatumPolicy zhangUserDatumPolicy()
{
    const char* p=std::getenv("ZHANG_USER_DATUM_POLICY");
    const std::string value=p?p:"TRANSPORT_KNOWN";
    if(value.empty() || value=="TRANSPORT_KNOWN")return ZhangUserDatumPolicy::TRANSPORT_KNOWN;
    if(value=="RESET_KNOWN_LOCAL")return ZhangUserDatumPolicy::RESET_KNOWN_LOCAL;
    throw std::runtime_error("UNKNOWN_ZHANG_USER_DATUM_POLICY");
}
inline const char* zhangUserDatumPolicyName()
{return zhangUserDatumPolicy()==ZhangUserDatumPolicy::TRANSPORT_KNOWN?"TRANSPORT_KNOWN":"RESET_KNOWN_LOCAL";}
