#include "user_chart_fixture.hpp"
#include "common/acsConfig.hpp"
#include "common/algebra.hpp"
#include "common/constants.hpp"
#include "common/observations.hpp"
#include "common/receiver.hpp"
#include "common/satStat.hpp"
#include "common/zhangCheckpoint.hpp"
#include "common/zhangUserArcGeneration.hpp"
#include "common/zhangUserGonsisEpoch.hpp"
#include "common/zhangUserProductChart.hpp"
#include "pea/zhangPppAr.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
bool sameMoments(const KFState& a,const KFState& b)
{
    return a.x.size()==b.x.size() && a.P.rows()==b.P.rows() && a.P.cols()==b.P.cols() &&
        a.dx.size()==b.dx.size() &&
        std::memcmp(a.x.data(),b.x.data(),a.x.size()*sizeof(double))==0 &&
        std::memcmp(a.P.data(),b.P.data(),a.P.size()*sizeof(double))==0 &&
        std::memcmp(a.dx.data(),b.dx.data(),a.dx.size()*sizeof(double))==0;
}

struct Fixture
{
    ReceiverMap receivers;
    KFState state;
    KFMeas measurement;
    unsigned callbackCalls=0;

    explicit Fixture(unsigned scenario)
    {
        auto& receiver=receivers["COMPONENT_TEST"];
        receiver.id="COMPONENT_TEST";
        receiver.ready=true;
        state.time.bigTime=100000+30*scenario;
        if(!bindZhangCheckpointRuntimeId(state,"gonsis-component-hook-"+std::to_string(scenario)))
            throw std::runtime_error("component fixture owner bind failed");
        state.kfIndexMap.clear();
        state.kfIndexMap[state.oneKey]=0;
        auto add=[&](KFKey key)
        {
            const int column=state.kfIndexMap.size();
            state.kfIndexMap[key]=column;
            state.stateTransitionMap[key][key][0]=1;
        };
        for(int axis=0;axis<3;++axis)
        {
            KFKey key; key.type=KF::REC_POS; key.str=receiver.id;
            key.num=axis; key.rec_ptr=&receiver; add(key);
        }
        for(auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
        {
            KFKey key; key.type=KF::PHASE_BIAS; key.str=receiver.id;
            key.Sat=SatSys(E_Sys::GPS,0); key.num=static_cast<int>(code);
            key.rec_ptr=&receiver; add(key);
        }
        for(int number=1;number<=12;++number)
        {
            const SatSys sat(E_Sys::GPS,number);
            auto& status=receiver.satStatMap[sat];
            status.el=(85-3*(number-1))*D2R;
            for(auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
            {
                const auto frequency=code2Freq.at(E_Sys::GPS).at(code);
                auto& tracking=status.sigStatMap[ft2string(frequency)];
                tracking.tracking.seen=true;
                tracking.tracking.lastCode=static_cast<int>(code);
                tracking.tracking.lastPhaseTime=double(state.time.bigTime);
                status.sigStatMap[enum_to_string(code)]={};
                if(number==1) continue; // Actual omitted KF chart reference.
                KFKey key; key.type=KF::AMBIGUITY; key.str=receiver.id;
                key.Sat=sat; key.num=static_cast<int>(code); key.rec_ptr=&receiver;
                add(key);
            }
            GObs observation; observation.Sat=sat; observation.time=state.time;
            observation.mount=receiver.id; observation.satStat_ptr=&status;
            for(auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
            {
                const auto frequency=code2Freq.at(E_Sys::GPS).at(code);
                observation.sigs[frequency].code=code;
                observation.sigs[frequency].P=2e7;
                observation.sigs[frequency].L=frequency==F1?1e8:8e7;
            }
            receiver.obsList.push_back(static_cast<std::shared_ptr<GObs>>(observation));
        }
        const int dimension=state.kfIndexMap.size();
        state.x=VectorXd::Zero(dimension); state.x(0)=1;
        state.dx=VectorXd::Zero(dimension);
        MatrixXd factor=MatrixXd::Identity(dimension,dimension)*.01;
        for(int i=1;i<dimension;++i)
            for(int j=1;j<dimension;++j) factor(i,j)+=.00003*(1+(i+2*j)%7);
        factor.row(0).setZero(); factor.col(0).setZero();
        state.P=factor*factor.transpose();
        for(const auto& [key,column] : state.kfIndexMap)
        {
            if(key.type==KF::REC_POS) state.x(column)=key.num==0?6378137:0;
            if(key.type!=KF::AMBIGUITY) continue;
            const int local=(key.Sat.prn-1)%6;
            const bool second=key.num==static_cast<int>(E_ObsCode::L2W);
            state.x(column)=(second?local:2*local)+(second?.0005:.001)*local;
            if(key.Sat.prn>=7) state.x(column)+=second?85.11:100.37;
        }
        std::ostringstream referenceTrace;
        initializeFixtureChart(referenceTrace,receivers,state);
        if (scenario==3)
            for (const auto& [key,column] : state.kfIndexMap)
                if (key.type==KF::AMBIGUITY && key.Sat.prn>=7)
                {
                    const double offset=(key.num==static_cast<int>(E_ObsCode::L1C)?.37:-.23)*(key.Sat.prn-6);
                    state.x(column)+=offset;
                    zhangUserSetChartOffset(state,receiver.id,key.Sat,key.num,offset);
                }
        if (scenario==4)
            for (const auto& [position,row] : state.kfIndexMap)
                if (position.type==KF::REC_POS)
                    for (const auto& [other,column] : state.kfIndexMap)
                        if (other.type!=KF::REC_POS)
                        { state.P(row,column)*=.5; state.P(column,row)*=.5; }
        if(zhangPppArUserReference(state,receiver.id,E_Sys::GPS,E_ObsCode::L1C)!=SatSys("G01") ||
            zhangPppArUserReference(state,receiver.id,E_Sys::GPS,E_ObsCode::L2W)!=SatSys("G01"))
            throw std::runtime_error("fixture did not establish actual G01 chart reference");
        // Scope changes occur only AFTER establishing the real chart. C1
        // references are ordinary state-backed rows, never invented zeroes.
        if(scenario>0)
            for(auto& obs : only<GObs>(receiver.obsList))
                if(obs.Sat.prn>=2 && obs.Sat.prn<=6) obs.excludeConfig=true;
        if(scenario==2)
            for(int number=11;number<=12;++number)
                receiver.satStatMap[SatSys(E_Sys::GPS,number)].el=14*D2R;

        state.lastPppTransactionCommitted=true;
        state.lastPppQcConverged=true;
        state.factorCommitSequence=90+scenario;
        KFMeasEntryList entries;
        for(const auto& obs : only<GObs>(receiver.obsList))
        {
            if(obs.exclude) continue;
            for(auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                for(auto type : {KF::CODE_MEAS,KF::PHAS_MEAS})
                {
                    KFMeasEntry entry(&state);
                    entry.obsKey.type=type; entry.obsKey.str=receiver.id;
                    entry.obsKey.Sat=obs.Sat; entry.obsKey.num=static_cast<int>(code);
                    entry.obsKey.rec_ptr=&receiver;
                    for(const auto& [key,column] : state.kfIndexMap)
                    {
                        if(key.type==KF::REC_POS) entry.designEntryMap[key]=.2+.1*key.num;
                        if(type==KF::PHAS_MEAS && key.type==KF::PHASE_BIAS &&
                            key.num==static_cast<int>(code)) entry.designEntryMap[key]=1;
                        if(type==KF::PHAS_MEAS && key.type==KF::AMBIGUITY &&
                            key.Sat==obs.Sat && key.num==static_cast<int>(code))
                            entry.designEntryMap[key]=genericWavelength.at(code2Freq.at(E_Sys::GPS).at(code));
                    }
                    entry.addNoiseEntry(entry.obsKey,1,type==KF::PHAS_MEAS?.003*.003:.3*.3);
                    entry.setValue(0); entries.push_back(entry);
                }
        }
        measurement=KFMeas(state,entries,state.time);
        measurement.V.setZero(); measurement.VV.setZero();
        state.acceptedMeasurementFactorCallback=[&](auto&&...){++callbackCalls;return true;};
        state.stateTransitionFactorCallback=[&](auto&&...){++callbackCalls;return true;};
        state.exactStateTransformCallback=[&](auto&&...){++callbackCalls;return true;};
        state.authoritativeFactorTransactionCallback=[&](auto&&...){++callbackCalls;return true;};
    }
};

unsigned occurrences(const std::string& text,const std::string& token)
{
    unsigned count=0;
    for(std::size_t position=0;(position=text.find(token,position))!=std::string::npos;position+=token.size()) ++count;
    return count;
}
}

extern "C" int __wrap_main()
{
    unsigned checks=0;
    auto check=[&](bool valid,const char* label)
    {
        ++checks;
        if(!valid) throw std::runtime_error(std::string("FAIL ")+label);
    };
    try
    {
        for(const char* name : {"ZHANG_USER_AR_MIN_ARC_SECONDS","ZHANG_USER_AR_MIN_ARC_SAMPLES",
            "ZHANG_USER_AR_CONFIRM_EPOCHS","ZHANG_USER_AR_MAX_CONDITIONAL_SIGMA_CYCLES",
            "ZHANG_USER_AR_PAR_MIN_RANK","ZHANG_USER_AR_PAR_MAX_TRIALS"}) unsetenv(name);
        setenv("ZHANG_USER_AR_RELIABILITY_PROFILE","GONSIS_UDUC",1);
        setenv("ZHANG_USER_AR_SEARCH","GONSIS_UDUC",1);
        setenv("ZHANG_USER_QC","CAUSAL_TRACKING",1);
        setenv("ZHANG_USER_AR_FEEDBACK","EPOCH_ONLY",1);
        setenv("ZHANG_USER_ENFORCE_AR_ELEVATION","1",1);
        setenv("ZHANG_R51_RATIO_ONLY","1",1);
        setenv("ZHANG_R51_EXPENSIVE_DIAGNOSTICS","0",1);
        acsConfig.process_ppp=true;
        acsConfig.zhangPppAr.user_adapter=true;
        acsConfig.zhangPppAr.output_products=false;
        acsConfig.zhangPppAr.product_solution="PRODUCT_FIXED";
        acsConfig.zhangPppAr.integer_strategy="CANONICAL_USER_SD_WL_L1";
        acsConfig.zhangPppAr.user_use_full_product_covariance=false;
        acsConfig.zhangFullRank.enable=false; acsConfig.phaseClockOsb.enable=false;
        acsConfig.pppOpts.ionoOpts.use_if_combo=false;
        acsConfig.zhangPppAr.baseline_observables[E_Sys::GPS]={E_ObsCode::L1C,E_ObsCode::L2W};
        acsConfig.solve_amb_for[E_Sys::GPS]=true;
        acsConfig.ambrOpts.mode=E_ARmode::LAMBDA_ALT;
        acsConfig.ambrOpts.ratioThres=2; acsConfig.ambrOpts.succsThres=.999;
        acsConfig.ambrOpts.elevation_mask_deg=15;
        acsConfig.ambrOpts.lambda_set=2; acsConfig.ambrOpts.AR_max_itr=100;
        acsConfig.assign_closest_epoch=false; acsConfig.epoch_tolerance=.5;
        const auto path=std::filesystem::current_path()/"user_gonsis_components_fixture.csv";
        std::ofstream products(path);
        products<<"gpst_seconds,solution,satellite,observable,clock_m,clock_sigma_m,phase_m,phase_sigma_m,clock_phase_covariance_m2,correction_m,correction_sigma_m,discontinuity_counter,integer_shift_cycles,fractional_shift_cycles,datum_version,valid_from_gpst_seconds,product_iod,reset_reason,persistent_relation_known,current_alignment_state,integer_structure_valid,integer_datum_continuous,integer_precision_valid,integer_valid,integer_component_id,integer_datum_id,solution_interval_start_gpst_seconds,solution_interval_end_gpst_seconds,numeric_valid,branch_valid,continuity_valid,ppp_usable,pppar_usable,invalid_reason,phase_product_segment_id,integer_component_version,integer_alignment_generation,real_gauge_generation,integer_component_size,integer_component_rank,certified_relation_count,redundant_relation_count,cycle_closure_valid,dual_frequency_ar_valid,product_state,discontinuity,ar_valid,support_segment_fingerprint,backend_s_basis_generation\n";
        for(int scenario=0;scenario<5;++scenario)
            for(int number=1;number<=12;++number)
                for(auto code : {E_ObsCode::L1C,E_ObsCode::L2W})
                {
                    const int time=100000+30*scenario;
                    products<<time<<",PRODUCT_FIXED,"<<SatSys(E_Sys::GPS,number).id()<<','<<enum_to_string(code)
                        <<",0,.01,0,.01,0,0,.02,0,0,0,0,"<<time
                        <<",1,NONE,1,CURRENT_ALIGNMENT_VALID,1,1,1,1,PRODUCT-C"<<(number<=6?0:1)
                        <<",DATUM,"<<time<<','<<time
                        <<",1,1,1,1,1,NONE,SEG0,28,0,0,6,5,5,0,1,1,PRODUCT_FIXED,0,1,FIXTURE,0\n";
                }
        products.close();
        acsConfig.zhangPppAr.product_filename=path.string();
        VectorXd firstPositionCorrection;
        for(unsigned scenario=0;scenario<5;++scenario)
        {
            Fixture fixture(scenario);
            const auto& receiver=fixture.receivers.at("COMPONENT_TEST");
            check(!zhangPppArUserAmbiguityIntegerValid(fixture.state,receiver.id,SatSys("G07"),E_ObsCode::L1C),
                "strict old getter still refuses C1 coordinate relative to C0 chart");
            const KFState frozen=fixture.state;
            captureZhangUserGonsisFloatObservations(fixture.state,fixture.measurement);
            KFState output=fixture.state;
            std::ostringstream trace;
            fixEpochOnlyUserAmbiguities(trace,output,fixture.state);
            const std::string evidence=trace.str();
            check(sameMoments(fixture.state,frozen) && fixture.state.kfIndexMap==frozen.kfIndexMap &&
                fixture.state.metaDataMap==frozen.metaDataMap &&
                fixture.state.factorCommitSequence==frozen.factorCommitSequence && fixture.callbackCalls==0,
                "public production entry preserves authoritative FLOAT bitwise and callbacks");
            check(occurrences(evidence,"\nZHANG_USER_GONSIS_AR_SUMMARY ")==1 &&
                occurrences(evidence,"\nZHANG_USER_GONSIS_OUTPUT_TRANSACTION ")==1,
                "multiple components publish one summary and one joint transaction");
            check(evidence.find("component_ils=INDEPENDENT")!=std::string::npos &&
                evidence.find("component_domain=")!=std::string::npos &&
                evidence.find("nis_evaluations=0")!=std::string::npos,
                "production evidence states independent component searches and zero NIS");
            if(scenario==2)
            {
                check(sameMoments(output,frozen) && evidence.find("committed_rank=0 reason=PASS")!=std::string::npos,
                    "four C1 satellites never bypass minimum four WL pairs");
                check(evidence.find("ZHANG_USER_GONSIS_L1_SEARCH")==std::string::npos,
                    "insufficient legal pairs do not invoke ILS");
                continue;
            }
            const int components=scenario==0?2:1;
            const int rows=10*components;
            check(evidence.find("committed_rank="+std::to_string(rows)+" reason=PASS")!=std::string::npos,
                "actual LAMBDA and joint gates accept all component WL L1 rows");
            check(occurrences(evidence," reason=ACCEPTED next_removed_satellite=")==components,
                "exactly one accepted PAR trial per supported component");
            MatrixXd independent=MatrixXd::Zero(rows,frozen.x.size());
            VectorXd values(rows),variances(rows);
            const double l1=genericWavelength.at(F1),l2=genericWavelength.at(F2);
            const double nl=l1*l2/(l1+l2),wl=l1*l1*l2/(l2*l2-l1*l1);
            int offset=0;
            for(int component=scenario==0?0:1;component<2;++component)
            {
                const int reference=1+6*component;
                auto addCoefficient=[&](int row,int number,E_ObsCode code,double value)
                {
                    if(number==1) return;
                    KFKey key; key.type=KF::AMBIGUITY; key.str=receiver.id;
                    key.Sat=SatSys(E_Sys::GPS,number); key.num=static_cast<int>(code);
                    independent(row,frozen.kfIndexMap.at(key))+=value;
                    if (scenario==3)
                        independent(row,frozen.kfIndexMap.at(frozen.oneKey))-=value*
                            (code==E_ObsCode::L1C?.37:-.23)*(number-6);
                };
                for(int local=1;local<=5;++local)
                {
                    const int rover=reference+local,row=offset+local-1;
                    addCoefficient(row,rover,E_ObsCode::L1C,1);
                    addCoefficient(row,rover,E_ObsCode::L2W,-1);
                    addCoefficient(row,reference,E_ObsCode::L1C,-1);
                    addCoefficient(row,reference,E_ObsCode::L2W,1);
                    independent.row(row+5)=wl*independent.row(row);
                    addCoefficient(row+5,rover,E_ObsCode::L1C,nl);
                    addCoefficient(row+5,reference,E_ObsCode::L1C,-nl);
                    values(row)=local; values(row+5)=nl*(2*local)+wl*local;
                    variances(row)=1e-8; variances(row+5)=.003*.003;
                }
                offset+=10;
            }
            const MatrixXd cross=frozen.P*independent.transpose();
            MatrixXd innovation=independent*cross;
            innovation.diagonal()+=variances;
            Eigen::FullPivLU<MatrixXd> solve(innovation);
            const VectorXd expectedMean=frozen.x+cross*solve.solve(values-independent*frozen.x);
            const MatrixXd expectedP=frozen.P-cross*solve.solve(cross.transpose());
            check((output.x-expectedMean).norm()<1e-8 && (output.P-expectedP).norm()<1e-10,
                "production complete mean P matches independent direct Gaussian oracle");
            if (scenario==1 || scenario==4)
            {
                VectorXd correction(3);
                for (const auto& [key,index] : frozen.kfIndexMap)
                    if (key.type==KF::REC_POS) correction(key.num)=output.x(index)-frozen.x(index);
                if (scenario==1) firstPositionCorrection=correction;
                else check(firstPositionCorrection.norm()>1e-9 &&
                    (correction-.5*firstPositionCorrection).norm()<2e-9,
                    "same ambiguity search marginal but different current Cxa gives different fixed position");
            }
            if(scenario==0)
            {
                KFKey left,right; left.type=right.type=KF::AMBIGUITY;
                left.str=right.str=receiver.id; left.num=right.num=static_cast<int>(E_ObsCode::L1C);
                left.Sat=SatSys("G02"); right.Sat=SatSys("G08");
                check(std::abs(output.P(frozen.kfIndexMap.at(left),frozen.kfIndexMap.at(right)))>1e-12,
                    "real joint output retains cross component posterior covariance");
            }
            if(scenario==1)
            {
                // Same product identity and same physical tracking serial,
                // but an internal coordinate has been retired/re-created.
                // Exercise the REAL adapter via the production AR entry, not
                // a copied identity-building helper.
                const auto physicalKey=zhangUserArcGenerationKey(true,receiver.id,"G07",
                    static_cast<int>(E_ObsCode::L1C));
                const auto stateKey=zhangUserArcGenerationKey(false,receiver.id,"G07",
                    static_cast<int>(E_ObsCode::L1C));
                const auto physicalGeneration=zhangUserArcGeneration(fixture.state.metaDataMap,physicalKey);
                zhangIncrementUserArcGeneration(fixture.state.metaDataMap,stateKey);
                const KFState reinitialised=fixture.state;
                captureZhangUserGonsisFloatObservations(fixture.state,fixture.measurement);
                KFState repeatOutput=fixture.state;
                std::ostringstream repeatTrace;
                fixEpochOnlyUserAmbiguities(repeatTrace,repeatOutput,fixture.state);
                const auto repeatEvidence=repeatTrace.str();
                const auto lineStart=repeatEvidence.find(" receiver="+receiver.id+" satellite=G07 provider=");
                check(lineStart!=std::string::npos,"state-only reset has an actual G07 endpoint receipt");
                const auto endpointLine=repeatEvidence.substr(lineStart,repeatEvidence.find('\n',lineStart)-lineStart);
                check(endpointLine.find("state_reason=IDENTITY_CHANGED")!=std::string::npos &&
                    endpointLine.find("physical_reason=DUPLICATE_EPOCH")!=std::string::npos &&
                    endpointLine.find("product_reason=DUPLICATE_EPOCH")!=std::string::npos,
                    "STATE-only generation changes state identity without product or physical change");
                check(zhangUserArcGeneration(fixture.state.metaDataMap,physicalKey)==physicalGeneration &&
                    sameMoments(fixture.state,reinitialised) && fixture.callbackCalls==0,
                    "state-only adapter evidence does not invent a physical break or mutate FLOAT");
            }
        }
        std::filesystem::remove(path);
        std::cout<<"PASS user_gonsis_components_integration_contract checks="<<checks
                 <<" affine_integer_chart_search=1 current_Cxa_output=1"
                 <<" real_product_parser_reference_ar_hooks=1 independent_complete_P_oracle=1"
                 <<" synthetic_contract_not_real_data_quality=1\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return EXIT_FAILURE;
    }
}
