#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// Wire v2: an endpoint-bound deterministic PUBLIC-to-native coordinate map.
// Shift signs are in C=clock-phase, not in the phase-bias column itself.
// This record contains publisher attestation, not independent integer truth.
struct ZhangProductDatumEndpoint
{
    double epoch=std::numeric_limits<double>::quiet_NaN();
    int counter=-1,datum=-1,componentVersion=-1,alignment=-1;
    std::string satellite,observable,reference,segment,physicalFunctional,component,posterior,bundle;
};

struct ZhangProductDatumTransition
{
    bool hasSource=false;
    ZhangProductDatumEndpoint source,target;
    long long correctionShiftCycles=0;
    double wavelength=0,correctionShiftMetres=0;
    std::vector<std::string> integerProofIds;
    std::string dualFrequencyBinding;
};

inline bool zhangDatumToken(const std::string& x)
{ return !x.empty() && x!="UNAVAILABLE" && x!="UNRESOLVED" && x!="NONE" && x.find_first_of(",\r\n")==std::string::npos; }

inline bool zhangDatumEndpointValid(const ZhangProductDatumEndpoint& e)
{
    return std::isfinite(e.epoch) && e.counter>=0 && e.datum>=0 && e.componentVersion>=0 && e.alignment>=0 &&
        zhangDatumToken(e.satellite) && zhangDatumToken(e.observable) && zhangDatumToken(e.reference) &&
        zhangDatumToken(e.segment) && zhangDatumToken(e.physicalFunctional) && zhangDatumToken(e.component) &&
        zhangDatumToken(e.posterior) && zhangDatumToken(e.bundle);
}

inline void zhangWriteDatumEndpoint(std::ostream& s,const ZhangProductDatumEndpoint& e)
{
    s<<std::setprecision(17)<<e.epoch<<' '<<e.counter<<' '<<e.datum<<' '<<e.componentVersion<<' '<<e.alignment;
    for(const auto* x:{&e.satellite,&e.observable,&e.reference,&e.segment,&e.physicalFunctional,&e.component,&e.posterior,&e.bundle})s<<' '<<std::quoted(*x);
}

inline bool zhangReadDatumEndpoint(std::istream& s,ZhangProductDatumEndpoint& e)
{
    // Non-finite default endpoints are represented by a finite zero when no
    // source is present; the hasSource bit owns their interpretation.
    if(!(s>>e.epoch>>e.counter>>e.datum>>e.componentVersion>>e.alignment))return false;
    for(auto* x:{&e.satellite,&e.observable,&e.reference,&e.segment,&e.physicalFunctional,&e.component,&e.posterior,&e.bundle})
        if(!(s>>std::quoted(*x)))return false;
    return true;
}

inline std::string zhangDatumHex(const std::string& x)
{
    static const char digits[]="0123456789abcdef";std::string out;out.reserve(2*x.size());
    for(unsigned char c:x){out.push_back(digits[c>>4]);out.push_back(digits[c&15]);}return out;
}

inline bool zhangDatumUnhex(const std::string& x,std::string& out)
{
    if(x.empty() || x.size()%2 || x.size()>1024*1024)return false;
    auto digit=[](char c){return c>='0'&&c<='9'?c-'0':c>='a'&&c<='f'?c-'a'+10:-1;};
    out.clear();out.reserve(x.size()/2);
    for(std::size_t i=0;i<x.size();i+=2){int a=digit(x[i]),b=digit(x[i+1]);if(a<0||b<0)return false;out.push_back(char(16*a+b));}
    return true;
}

inline std::string zhangEncodeDatumTransition(const ZhangProductDatumTransition& t)
{
    std::ostringstream s;s<<1<<' '<<t.hasSource<<' ';auto source=t.source;
    if(!t.hasSource)source.epoch=0;
    zhangWriteDatumEndpoint(s,source);s<<' ';zhangWriteDatumEndpoint(s,t.target);
    s<<' '<<t.correctionShiftCycles<<' '<<std::setprecision(17)<<t.wavelength<<' '<<t.correctionShiftMetres<<' '
        <<std::quoted(t.dualFrequencyBinding)<<' '<<t.integerProofIds.size();
    for(const auto& p:t.integerProofIds)s<<' '<<std::quoted(p);
    return zhangDatumHex(s.str());
}

inline bool zhangDecodeDatumTransition(const std::string& text,ZhangProductDatumTransition& out)
{
    std::string decoded;if(!zhangDatumUnhex(text,decoded))return false;
    std::istringstream s(decoded);int version=0,hasSource=0;std::size_t n=0;ZhangProductDatumTransition t;
    if(!(s>>version>>hasSource) || version!=1 || (hasSource!=0&&hasSource!=1) ||
        !zhangReadDatumEndpoint(s,t.source) || !zhangReadDatumEndpoint(s,t.target) ||
        !(s>>t.correctionShiftCycles>>t.wavelength>>t.correctionShiftMetres>>std::quoted(t.dualFrequencyBinding)>>n) || n>10000)return false;
    t.hasSource=hasSource;
    for(std::size_t i=0;i<n;++i){std::string p;if(!(s>>std::quoted(p)) || !zhangDatumToken(p))return false;t.integerProofIds.push_back(std::move(p));}
    s>>std::ws;if(!s.eof())return false;out=std::move(t);return true;
}

inline std::string zhangDatumProofId(const ZhangProductDatumTransition& t)
{
    // Content fingerprint (not a cryptographic signature/authentication).
    std::uint64_t h=1469598103934665603ULL;
    for(unsigned char c:zhangEncodeDatumTransition(t)){h^=c;h*=1099511628211ULL;}
    std::ostringstream s;s<<"R51AFFINE-V2-"<<std::hex<<h;return s.str();
}

inline bool zhangDatumTransitionPureIntegerValid(const ZhangProductDatumTransition& t)
{
    if(!t.hasSource || !zhangDatumEndpointValid(t.source) || !zhangDatumEndpointValid(t.target) ||
        !(t.source.epoch<t.target.epoch) || t.source.satellite!=t.target.satellite || t.source.observable!=t.target.observable ||
        !(t.wavelength>0) || !std::isfinite(t.wavelength) || !std::isfinite(t.correctionShiftMetres) ||
        !zhangDatumToken(t.dualFrequencyBinding) || t.integerProofIds.empty() ||
        std::abs(static_cast<long double>(t.correctionShiftCycles))>4503599627370496.L)return false;
    for(const auto& proof:t.integerProofIds)if(!zhangDatumToken(proof))return false;
    const double expected=t.wavelength*static_cast<double>(t.correctionShiftCycles);
    return std::isfinite(expected) && std::abs(expected-t.correctionShiftMetres)<=
        16*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(expected));
}

inline std::string zhangDatumEndpointIdentity(const ZhangProductDatumEndpoint& e)
{
    std::ostringstream s;zhangWriteDatumEndpoint(s,e);return s.str();
}

inline std::string zhangDatumWireProof(const ZhangProductDatumTransition& t)
{ return "R51DATUM2:"+zhangEncodeDatumTransition(t); }

inline bool zhangReadDatumWireProof(const std::string& proof,ZhangProductDatumTransition& t)
{
    const std::string prefix="R51DATUM2:";
    return proof.compare(0,prefix.size(),prefix)==0 &&
        zhangDecodeDatumTransition(proof.substr(prefix.size()),t);
}

// The sign is checked with two representations at ONE physical instant.
// Neither an adjacent-epoch product difference nor a rounded float is proof.
inline bool zhangSameTimeNativeCorrectionShift(double clock,double phase,
    double lambda,long long gamma,double& shift)
{
    if(!std::isfinite(clock) || !std::isfinite(phase) || !(lambda>0) ||
        !std::isfinite(lambda) || std::abs(static_cast<long double>(gamma))>4503599627370496.L)return false;
    const double delta=lambda*static_cast<double>(gamma);
    const double native=clock-phase,aligned=clock-(phase+delta);
    const double error=(native-aligned)-delta;
    const double bound=16*std::numeric_limits<double>::epsilon()*
        std::max({1.,std::abs(clock),std::abs(phase),std::abs(delta)});
    if(!std::isfinite(delta) || !std::isfinite(aligned) || !std::isfinite(error) || std::abs(error)>bound)return false;
    shift=delta;return true;
}

inline bool zhangDatumDualPairValid(const ZhangProductDatumTransition& a,
    const ZhangProductDatumTransition& b)
{
    auto pair=[](const auto& x,const auto& y) {
        return x.epoch==y.epoch && x.satellite==y.satellite && x.observable!=y.observable &&
            x.reference==y.reference && x.posterior==y.posterior && x.bundle==y.bundle;
    };
    return zhangDatumTransitionPureIntegerValid(a) && zhangDatumTransitionPureIntegerValid(b) &&
        pair(a.source,b.source) && pair(a.target,b.target) &&
        a.dualFrequencyBinding==b.dualFrequencyBinding;
}

struct ZhangUserFloatMap
{
    bool valid=false,randomMapAvailable=false;
    ZhangProductDatumEndpoint source,target;
    double deterministicShiftMetres=0;
    std::string reason="FLOAT_MAP_UNKNOWN";
};
struct ZhangUserIntegerMap
{
    bool valid=false,identityCanonical=true;
    long long correctionShiftCycles=0;
    std::string sourcePhysicalIdentity,targetPhysicalIdentity,proofId,reason="INTEGER_MAP_UNKNOWN";
};
