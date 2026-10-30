// Test the real private receiver DTO and validation code. Link without its
// separate production object; all other objects are the frozen current build.
#include "pea/zhangReceiverCheckpoint.cpp"
#include <bit>
#include <iostream>

extern "C" int __wrap_main(int,char**)
{
    int checks=0;
    auto require=[&](bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);};
    try
    {
        for(const char* sat:{"E02","C20"})
        {
            lc_t native;native.Sat=SatSys(sat);native.time.bigTime=1474855230;
            native.L_m[F5]=23557868.628;native.P[F5]=23557842.306;
            auto value=getLC(native.L_m[F5],native.L_m[F5],native.P[F5],native.P[F5],
                0.25482804879085386,0.25482804879085386,nullptr,nullptr);
            value.MW_c*=-1;native.lcMap[{F5,F5}]=value;
            SatNav satNav;satNav.lamMap[F5]=value.lam_A;
            GObs observation;observation.Sat=native.Sat;observation.satNav_ptr=&satNav;
            lc_t legacy;legacy.Sat=native.Sat;legacy.L_m=native.L_m;legacy.P=native.P;
            getLC(observation,legacy,F5,F5);
            native.mp=legacy.mp;
            require(value.valid && std::isinf(value.lam_WL) && std::isnan(value.IF_Phas_m),"actual legacy singular combination");
            auto dto=captureLcState(native);std::string reason,bytes,reencoded;
            unsetenv("ZHANG_US20_RESET_CHECKPOINT");
            require(!validLcState(dto,reason),"legacy E29 must still reject nonfinite cache");
            setenv("ZHANG_US20_RESET_CHECKPOINT","1",1);
            std::feclearexcept(FE_ALL_EXCEPT);std::feraiseexcept(FE_INEXACT);
            const int exceptions=std::fetestexcept(FE_ALL_EXCEPT);
            require(validLcState(dto,reason),"exact US20 derived representation");
            require(std::fetestexcept(FE_ALL_EXCEPT)==exceptions,"capture validation preserves FP flags");
            require(serializePayload(dto,bytes,reason,"TEST"),"serialize all bits");
            ZhangRuntimeLcState decoded;
            require(deserializePayload(bytes,decoded,reason,"TEST") && validLcState(decoded,reason),"decode complete derived cache");
            require(serializePayload(captureLcState(restoreLcState(decoded)),reencoded,reason,"TEST") && bytes==reencoded,"live DTO roundtrip all bits");
            auto bad=dto;bad.satellite=SatSys("G01");require(!validLcState(bad,reason),"GPS outside supported degenerate scope");
            bad=dto;auto item=bad.combinations.extract({F5,F5});item.key()={F1,F5};bad.combinations.insert(std::move(item));
            require(!validLcState(bad,reason),"nonfinite actual dual frequency rejected");
            bad=dto;bad.phaseMetres[F5]=std::numeric_limits<double>::quiet_NaN();require(!validLcState(bad,reason),"nonfinite physical sample rejected");
            bad=dto;bad.combinations.at({F5,F5}).wavelengthB+=0.001;require(!validLcState(bad,reason),"nonidentical wavelength rejected");
            bad=dto;bad.combinations.at({F5,F5}).gfPhaseMetres=1;require(!validLcState(bad,reason),"corrupt finite derived field rejected");
            bad=dto;bad.combinations.at({F5,F5}).valid=false;require(!validLcState(bad,reason),"unexpected valid flag rejected");
            bad=dto;auto& nan=bad.combinations.at({F5,F5}).mwCycles;
            nan=std::bit_cast<double>(std::bit_cast<std::uint64_t>(nan)^std::uint64_t{1});
            require(!validLcState(bad,reason),"different NaN payload rejected");
            bad=dto;bad.codeMetres.erase(F5);require(!validLcState(bad,reason),"missing source code rejected");
            bad=dto;bad.phaseMetres[F1]=1;require(!validLcState(bad,reason),"multipath exception requires absent first phase");
            bad=dto;auto& multipath=bad.multipathMetres.at(F5);
            multipath=std::bit_cast<double>(std::bit_cast<std::uint64_t>(multipath)^std::uint64_t{1});
            require(!validLcState(bad,reason),"different derived multipath NaN rejected");
            setenv("ZHANG_US20_RESET_CHECKPOINT","0",1);require(!validLcState(dto,reason),"inactive mode remains strict");
        }
        setenv("ZHANG_US20_RESET_CHECKPOINT","1",1);
        const auto site=siteBindingFingerprint(&dummySiteid);
        const auto receiver=receiverBindingFingerprint(&dummyReceiver);
        const auto antenna=antennaBindingFingerprint(&dummyAntenna);
        const auto ecc=eccentricityBindingFingerprint(&dummySiteEcc);
        require(resolveSiteBinding(site)==&dummySiteid,"exact default site binding");
        require(resolveReceiverBinding(receiver)==&dummyReceiver,"exact default receiver binding");
        require(resolveAntennaBinding(antenna)==&dummyAntenna,"exact default antenna binding");
        require(resolveEccentricityBinding(ecc)==&dummySiteEcc,"exact default eccentricity binding");
        require(!resolveSiteBinding(site+"0") && !resolveAntennaBinding(ecc),"changed or wrong-type identity rejected");
        const auto saved=dummySiteid.desc;dummySiteid.desc="changed default contents";
        require(!resolveSiteBinding(site),"placeholder contents must match exactly");dummySiteid.desc=saved;
        SinexSiteId copy=dummySiteid;
        require(siteBindingFingerprint(&copy)!=site,"equal record values cannot impersonate default object identity");
        copy.sitecode="TEST";theSinex.mapsiteids.emplace("TEST",copy);
        const auto real=siteBindingFingerprint(&theSinex.mapsiteids.at("TEST"));
        require(resolveSiteBinding(real)==&theSinex.mapsiteids.at("TEST"),"normal configured identity unchanged");
        unsetenv("ZHANG_US20_RESET_CHECKPOINT");
        require(!resolveSiteBinding(site) && !resolveReceiverBinding(receiver)
            && !resolveAntennaBinding(antenna) && !resolveEccentricityBinding(ecc),"US20 placeholders unavailable to legacy E29");
        require(resolveSiteBinding(real)==&theSinex.mapsiteids.at("TEST"),"legacy normal configured identity unchanged");
        std::cout<<"PASS exact singular LC and default SINEX checkpoint representation checks="<<checks<<'\n';return 0;
    }
    catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
