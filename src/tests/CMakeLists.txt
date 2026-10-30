set(UNIT_TEST_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/../cpp
    ${EIGEN3_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

function(configure_unit_test target)
    target_include_directories(${target} PRIVATE
        ${UNIT_TEST_INCLUDE_DIRS}
    )

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
    endif()

    add_test(NAME ${target} COMMAND ${target})
endfunction()

add_executable(zhang_float_rank_diagnostic_policy_tests unit/float_rank_diagnostic_policy.cpp)
configure_unit_test(zhang_float_rank_diagnostic_policy_tests)

add_executable(zhang_float_matrix_parallel_helper_tests unit/float_matrix_parallel_helpers.cpp)
configure_unit_test(zhang_float_matrix_parallel_helper_tests)
if(BLA_VENDOR STREQUAL "OpenBLAS")
    target_compile_definitions(zhang_float_matrix_parallel_helper_tests PRIVATE GINAN_USE_OPENBLAS=1)
    target_link_libraries(zhang_float_matrix_parallel_helper_tests PRIVATE ${BLAS_LIBRARIES})
endif()
set_tests_properties(zhang_float_matrix_parallel_helper_tests PROPERTIES
    ENVIRONMENT "OPENBLAS_NUM_THREADS=2;OMP_NUM_THREADS=4;MKL_NUM_THREADS=1" TIMEOUT 60)

add_executable(zhang_multi_gnss_clock_tests unit/multi_gnss_clock_contract.cpp)
configure_unit_test(zhang_multi_gnss_clock_tests)
target_link_libraries(zhang_multi_gnss_clock_tests PRIVATE Boost::program_options yaml-cpp)

add_executable(zhang_multi_gnss_integer_tests
    unit/multi_gnss_integer.cpp unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_multi_gnss_integer_tests)

add_executable(zhang_multi_gnss_product_manifest_tests
    unit/multi_gnss_product_manifest.cpp unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp ../cpp/common/zhangCheckpoint.cpp)
configure_unit_test(zhang_multi_gnss_product_manifest_tests)
target_link_libraries(zhang_multi_gnss_product_manifest_tests PRIVATE
    Boost::program_options yaml-cpp Boost::json Boost::serialization OpenSSL::Crypto)

set(RECEIVER_METADATA_TEST_SOURCES
    ../cpp/common/receiverMetadata.cpp
    unit/test_ReceiverMetaData.cpp
)

set(SANITY_CHECKER_TEST_SOURCES
    ../cpp/common/sanityCheckers/ConfigSanityManager.cpp
    ../cpp/common/sanityCheckers/EphemerisTimeDelayChecker.cpp
    ../cpp/common/sanityCheckers/EpochToleranceChecker.cpp
    ../cpp/common/sanityCheckers/IonosphericFreeComboChecker.cpp
    ../cpp/common/sanityCheckers/IonosphericOutageChecker.cpp
    ../cpp/common/sanityCheckers/RequiredSiteEccentricityChecker.cpp
    ../cpp/common/sanityCheckers/SbasSanityChecker.cpp
)

function(add_sanity_checker_test target test_source)
    add_executable(${target}
        ${SANITY_CHECKER_TEST_SOURCES}
        ${test_source}
    )

    target_link_libraries(${target} PRIVATE
        Boost::log
        Boost::log_setup
        Boost::program_options
        Boost::thread
        Boost::unit_test_framework
    )

    configure_unit_test(${target})
endfunction()

add_executable(receiver_metadata_tests
    ${RECEIVER_METADATA_TEST_SOURCES}
)

add_executable(trace_file_cache_tests
    unit/test_TraceFileCache.cpp
)

add_executable(phase_clock_osb_tests
    unit/test_PhaseClockOsb.cpp
)

add_executable(zhang_full_rank_tests
    unit/test_ZhangRatioOnly.cpp
    unit/test_ZhangR47.cpp
    unit/test_ZhangR48.cpp
    unit/test_ZhangSignalArc.cpp
	unit/test_ZhangSharedProductEvent.cpp
    unit/test_ZhangR46.cpp
    ../cpp/common/satSys.cpp
    unit/test_ZhangFullRank.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/zhangCheckpoint.cpp
)

# Keep the deterministic restart boundary tests in a small target.  Pulling
# these into zhang_full_rank_tests would make failures in cursor/controller
# checkpointing difficult to distinguish from the mathematical regression
# suite, and linking the complete PEA executable would hide missing test
# seams behind unrelated production dependencies.
add_executable(zhang_checkpoint_infra0_tests
    unit/test_ZhangCheckpointInfra0.cpp
    unit/test_ZhangCheckpointSupport.cpp
    unit/test_ZhangInputCheckpointSupport.cpp
    ../cpp/common/zhangCheckpoint.cpp
    ../cpp/pea/zhangInputCheckpoint.cpp
)

add_sanity_checker_test(config_sanity_manager_tests
    unit/sanityCheckers/test_ConfigSanityManager.cpp
)
add_sanity_checker_test(epoch_tolerance_checker_tests
    unit/sanityCheckers/test_EpochToleranceChecker.cpp
)
add_sanity_checker_test(ephemeris_time_delay_checker_tests
    unit/sanityCheckers/test_EphemerisTimeDelayChecker.cpp
)
add_sanity_checker_test(ionospheric_free_combo_checker_tests
    unit/sanityCheckers/test_IonosphericFreeComboChecker.cpp
)
add_sanity_checker_test(ionospheric_outage_checker_tests
    unit/sanityCheckers/test_IonosphericOutageChecker.cpp
)
add_sanity_checker_test(required_site_eccentricity_checker_tests
    unit/sanityCheckers/test_RequiredSiteEccentricityChecker.cpp
)
add_sanity_checker_test(sbas_sanity_checker_tests
    unit/sanityCheckers/test_SbasSanityChecker.cpp
)

add_custom_target(unit_tests
    DEPENDS zhang_float_rank_diagnostic_policy_tests
            zhang_float_matrix_parallel_helper_tests
            receiver_metadata_tests
            trace_file_cache_tests
            phase_clock_osb_tests
            zhang_full_rank_tests
            zhang_checkpoint_infra0_tests
            zhang_r51_float_reuse_contract_tests
            config_sanity_manager_tests
            epoch_tolerance_checker_tests
            ephemeris_time_delay_checker_tests
            ionospheric_free_combo_checker_tests
            ionospheric_outage_checker_tests
            required_site_eccentricity_checker_tests
            sbas_sanity_checker_tests
)

target_link_libraries(receiver_metadata_tests PRIVATE
    Boost::unit_test_framework
)

target_link_libraries(trace_file_cache_tests PRIVATE
    Boost::log
    Boost::unit_test_framework
)

target_link_libraries(phase_clock_osb_tests PRIVATE
    Boost::unit_test_framework
)

target_link_libraries(zhang_full_rank_tests PRIVATE
	Boost::json
    Boost::serialization
    Boost::unit_test_framework
    OpenSSL::Crypto
)

# Exercise the same Eigen backend as PEA.  Invalid BLAS arguments must fail a
# test, rather than merely printing XERBLA warnings and reporting success.
target_compile_definitions(zhang_full_rank_tests PRIVATE EIGEN_USE_BLAS=1)
target_link_libraries(zhang_full_rank_tests PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(UNIX AND NOT APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_sources(zhang_full_rank_tests PRIVATE unit/zhangBlasArgumentAudit.cpp)
    target_link_options(zhang_full_rank_tests PRIVATE "-Wl,--wrap=dgemv_")
endif()

target_link_libraries(zhang_checkpoint_infra0_tests PRIVATE
    Boost::log
    Boost::log_setup
    Boost::serialization
    Boost::unit_test_framework
    OpenSSL::Crypto
)

if(WIN32)
    target_link_libraries(receiver_metadata_tests PRIVATE ws2_32 wsock32)
endif()

configure_unit_test(receiver_metadata_tests)
configure_unit_test(trace_file_cache_tests)
configure_unit_test(phase_clock_osb_tests)
configure_unit_test(zhang_full_rank_tests)
configure_unit_test(zhang_checkpoint_infra0_tests)

add_executable(zhang_r51_performance_tests unit/r51_performance_reference.cpp)
add_executable(zhang_epoch_stop_contract_tests unit/r51_epoch_stop_contract.cpp)
configure_unit_test(zhang_epoch_stop_contract_tests)
add_dependencies(unit_tests zhang_epoch_stop_contract_tests)
configure_unit_test(zhang_r51_performance_tests)
target_compile_definitions(zhang_r51_performance_tests PRIVATE EIGEN_USE_BLAS=1)
target_link_libraries(zhang_r51_performance_tests PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(UNIX AND NOT APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_sources(zhang_r51_performance_tests PRIVATE unit/zhangBlasArgumentAudit.cpp)
    target_link_options(zhang_r51_performance_tests PRIVATE "-Wl,--wrap=dgemv_")
endif()

add_executable(zhang_r51_entailment_tests unit/r51_entailment_reference.cpp)
configure_unit_test(zhang_r51_entailment_tests)

add_executable(zhang_r51_compatibility_tests unit/r51_compatibility_reference.cpp)
configure_unit_test(zhang_r51_compatibility_tests)

add_executable(zhang_r51_float_reuse_contract_tests unit/r51_float_reuse_contract.cpp)
configure_unit_test(zhang_r51_float_reuse_contract_tests)

add_executable(zhang_r51_current_epoch_evidence_cache_tests
    unit/r51_current_epoch_evidence_cache.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_r51_current_epoch_evidence_cache_tests)

add_executable(zhang_r51_proof_dictionary_tests unit/r51_proof_dictionary.cpp)
configure_unit_test(zhang_r51_proof_dictionary_tests)

add_executable(zhang_r51_final_product_posterior_tests
    unit/r51_final_product_posterior.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_r51_final_product_posterior_tests)
target_compile_definitions(zhang_r51_final_product_posterior_tests PRIVATE EIGEN_USE_BLAS=1)
target_link_libraries(zhang_r51_final_product_posterior_tests PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(UNIX AND NOT APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_sources(zhang_r51_final_product_posterior_tests PRIVATE unit/zhangBlasArgumentAudit.cpp)
    target_link_options(zhang_r51_final_product_posterior_tests PRIVATE "-Wl,--wrap=dgemv_")
endif()

add_executable(zhang_r51_full_joint_handoff_tests
    unit/r51_full_joint_handoff.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_r51_full_joint_handoff_tests)

add_executable(zhang_r51_handoff_bundle_reader_tests
    unit/r51_handoff_bundle_reader.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp
    ../cpp/common/zhangCheckpoint.cpp)
configure_unit_test(zhang_r51_handoff_bundle_reader_tests)
target_link_libraries(zhang_r51_handoff_bundle_reader_tests PRIVATE
    Boost::json Boost::serialization OpenSSL::Crypto)

# Compile the production LAMBDA entry point in isolation. The exhaustive
# small-dimensional reference is independent of the production enumerator.
add_executable(zhang_gnss_deadline_tests
    unit/test_GNSSDeadline.cpp
    ../cpp/ambres/GNSSambres.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_gnss_deadline_tests)
target_compile_definitions(zhang_gnss_deadline_tests PRIVATE EIGEN_USE_BLAS=1)
target_link_libraries(zhang_gnss_deadline_tests PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(UNIX AND NOT APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Unused counterfactual beam/oracle entry points depend on RTKLIB helpers;
    # discard them, rather than replacing the tested search with a mock.
    target_compile_options(zhang_gnss_deadline_tests PRIVATE -ffunction-sections -fdata-sections)
    target_link_options(zhang_gnss_deadline_tests PRIVATE "-Wl,--gc-sections")
endif()

add_executable(zhang_product_conditioning_receipt_tests
    unit/test_ZhangProductConditioningReceipt.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_product_conditioning_receipt_tests)

add_executable(zhang_r51_private_graph_state_tests
    unit/r51_private_graph_state.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_r51_private_graph_state_tests)

add_executable(zhang_exact_work_budget_tests
    unit/test_ZhangExactWorkBudget.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_exact_work_budget_tests)

add_executable(zhang_r51_handoff_search_budget_tests
    unit/r51_handoff_search_budget.cpp)
configure_unit_test(zhang_r51_handoff_search_budget_tests)

add_executable(zhang_r51_ratio_policy_tests unit/r51_ratio_policy.cpp)
configure_unit_test(zhang_r51_ratio_policy_tests)

add_executable(zhang_r51_ratio_reuse_contract_tests unit/r51_ratio_reuse_contract.cpp)
configure_unit_test(zhang_r51_ratio_reuse_contract_tests)

add_executable(zhang_r51_handoff_work_identity_tests
    unit/r51_handoff_work_identity.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_r51_handoff_work_identity_tests)

add_executable(zhang_r51_verified_domain_image_tests
    unit/r51_verified_domain_image.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_r51_verified_domain_image_tests)

add_executable(zhang_r51_current_bridge_representation_tests
    unit/r51_current_bridge_representation.cpp)
configure_unit_test(zhang_r51_current_bridge_representation_tests)

add_executable(zhang_r51_bridge_support_analysis_tests unit/r51_bridge_support_analysis.cpp)
configure_unit_test(zhang_r51_bridge_support_analysis_tests)
add_executable(zhang_r51_exact_column_plan_tests unit/r51_exact_column_plan.cpp)
configure_unit_test(zhang_r51_exact_column_plan_tests)
add_dependencies(unit_tests zhang_r51_bridge_support_analysis_tests zhang_r51_exact_column_plan_tests)
add_executable(zhang_r51_structural_publication_gap_tests
    unit/r51_structural_publication_gap.cpp
    unit/test_ZhangCheckpointSupport.cpp
    ../cpp/common/satSys.cpp)
configure_unit_test(zhang_r51_structural_publication_gap_tests)
add_dependencies(unit_tests zhang_r51_structural_publication_gap_tests)
add_executable(zhang_r51_handoff_event_budget_tests unit/r51_handoff_event_budget.cpp)
configure_unit_test(zhang_r51_handoff_event_budget_tests)
add_dependencies(unit_tests zhang_r51_handoff_event_budget_tests)
add_executable(zhang_r51_bridge_event_replay_tests unit/r51_bridge_event_replay.cpp)
configure_unit_test(zhang_r51_bridge_event_replay_tests)
add_dependencies(unit_tests zhang_r51_bridge_event_replay_tests)

add_executable(zhang_cycle_slip_evidence_tests unit/test_ZhangCycleSlipEvidence.cpp)
configure_unit_test(zhang_cycle_slip_evidence_tests)
target_link_libraries(zhang_cycle_slip_evidence_tests PRIVATE Boost::serialization)
add_dependencies(unit_tests zhang_cycle_slip_evidence_tests)
add_test(NAME zhang_cycle_slip_no_slip_holdout
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo 1000 500)
add_test(NAME zhang_cycle_slip_injected_noise_readiness
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo-slip 1000 10000)
add_test(NAME zhang_cycle_slip_doppler_no_slip
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo 500 70000)
add_test(NAME zhang_cycle_slip_doppler_readiness
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo-slip 500 60000)
add_test(NAME zhang_cycle_slip_storm_no_slip
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo 500 80000)
add_test(NAME zhang_cycle_slip_storm_injected
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo-slip 500 90000)
add_test(NAME zhang_cycle_slip_doppler_blunder
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo 500 100000)
add_test(NAME zhang_cycle_slip_smooth_phase_trend
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo 500 150000)
add_test(NAME zhang_cycle_slip_no_doppler_iono_trend
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo 500 240000)
add_test(NAME zhang_cycle_slip_no_doppler_phase_trend
    COMMAND zhang_cycle_slip_evidence_tests --monte-carlo 500 190000)
set_tests_properties(zhang_cycle_slip_no_doppler_iono_trend
    PROPERTIES ENVIRONMENT "ZHANG_SLIP_TEST_IONO_STRESS=1")
set_tests_properties(zhang_cycle_slip_no_doppler_phase_trend
    PROPERTIES ENVIRONMENT "ZHANG_SLIP_TEST_PHASE_TREND=1")
set_tests_properties(zhang_cycle_slip_smooth_phase_trend
    PROPERTIES ENVIRONMENT "ZHANG_SLIP_TEST_DOPPLER=1;ZHANG_SLIP_TEST_PHASE_TREND=1")
set_tests_properties(zhang_cycle_slip_doppler_no_slip zhang_cycle_slip_doppler_readiness
    PROPERTIES ENVIRONMENT "ZHANG_SLIP_TEST_DOPPLER=1")
set_tests_properties(zhang_cycle_slip_storm_no_slip zhang_cycle_slip_storm_injected
    PROPERTIES ENVIRONMENT "ZHANG_SLIP_TEST_DOPPLER=1;ZHANG_SLIP_TEST_IONO_STRESS=1")
set_tests_properties(zhang_cycle_slip_doppler_blunder
    PROPERTIES ENVIRONMENT "ZHANG_SLIP_TEST_DOPPLER=1;ZHANG_SLIP_TEST_DOPPLER_BLUNDER=1")
# Do not invert or waive the readiness test: a reproducible weak-slip miss must
# remain a failing gate, even if deterministic fixtures and other regressions pass.

# Keep the public aggregate complete: a successful unit_tests build must
# produce every R51 regression executable registered with CTest.
add_dependencies(unit_tests
    zhang_multi_gnss_clock_tests
    zhang_multi_gnss_integer_tests
    zhang_multi_gnss_product_manifest_tests
    zhang_r51_performance_tests
    zhang_r51_entailment_tests
    zhang_r51_compatibility_tests
    zhang_r51_current_epoch_evidence_cache_tests
    zhang_r51_proof_dictionary_tests
    zhang_r51_final_product_posterior_tests
    zhang_r51_full_joint_handoff_tests
    zhang_r51_handoff_bundle_reader_tests
    zhang_gnss_deadline_tests
    zhang_r51_private_graph_state_tests
    zhang_product_conditioning_receipt_tests
    zhang_exact_work_budget_tests
    zhang_r51_handoff_search_budget_tests
    zhang_r51_ratio_policy_tests
    zhang_r51_ratio_reuse_contract_tests
    zhang_r51_handoff_work_identity_tests
    zhang_r51_verified_domain_image_tests
    zhang_r51_current_bridge_representation_tests)

# These lifecycle contracts exercise the real graph/KF registries and the real
# transaction code. Reuse the already-built PEA objects; do not replace them
# with header-only mocks. The harness currently supports GNU/Linux Makefiles.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    find_program(ZHANG_CONTRACT_PYTHON NAMES python3 REQUIRED)
    add_custom_target(zhang_clock_lifecycle_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_clock_lifecycle_contract.py ${CMAKE_BINARY_DIR} --build-only
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_clock_lifecycle_contract_build)
    add_test(NAME zhang_clock_lifecycle_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/clock_lifecycle_contract/clock_lifecycle_contract
            ${CMAKE_BINARY_DIR}/clock_lifecycle_contract/ctest.TRACE)
    set_tests_properties(zhang_clock_lifecycle_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;ZHANG_R51_RATIO_ONLY=1;ZHANG_R51_RATIO_THRESHOLD=2;ZHANG_R51_EPOCH_ONLY_AR=1;ZHANG_CLOCK_LIFECYCLE_DIAGNOSTICS=1"
        TIMEOUT 120)
    add_custom_target(zhang_multi_gnss_covariance_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/multi_gnss_clock_compile.py ${CMAKE_BINARY_DIR} covariance
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_multi_gnss_covariance_contract_build)
    add_test(NAME zhang_multi_gnss_covariance_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/multi_gnss_clock_contract/product_covariance_contract
            ${CMAKE_BINARY_DIR}/multi_gnss_clock_contract/fixtures)
    set_tests_properties(zhang_multi_gnss_covariance_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1" TIMEOUT 120)
    add_custom_target(zhang_multi_gnss_runtime_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_multi_gnss_runtime.py ${CMAKE_BINARY_DIR}
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_multi_gnss_runtime_contract_build)
    add_test(NAME zhang_multi_gnss_runtime_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/multi_gnss_runtime/multi_gnss_runtime
            ${CMAKE_BINARY_DIR}/multi_gnss_runtime/ctest.TRACE)
    set_tests_properties(zhang_multi_gnss_runtime_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;ZHANG_R51_ENABLE=1;ZHANG_R51_RATIO_ONLY=1;ZHANG_R51_RATIO_THRESHOLD=2;ZHANG_R51_EPOCH_ONLY_AR=1"
        TIMEOUT 120)
    add_custom_target(zhang_r51_slip_qc_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/run_r51_slip_qc.py ${CMAKE_BINARY_DIR}
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_r51_slip_qc_contract_build)
    add_test(NAME zhang_r51_slip_qc_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/slip_qc_contract/slip_qc_contract)
    set_tests_properties(zhang_r51_slip_qc_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
    add_custom_target(zhang_r51_graph_runtime_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_r51_graph_runtime.py ${CMAKE_BINARY_DIR}
        DEPENDS pea
        VERBATIM)
    add_dependencies(unit_tests zhang_r51_graph_runtime_contract_build)
    add_custom_target(zhang_r51_handoff_float_lifecycle_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_r51_handoff_float_lifecycle.py ${CMAKE_BINARY_DIR}
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_r51_handoff_float_lifecycle_build)
    add_test(NAME zhang_r51_handoff_float_lifecycle_tests
        COMMAND ${CMAKE_BINARY_DIR}/handoff_float_lifecycle/handoff_float_lifecycle)
    set_tests_properties(zhang_r51_handoff_float_lifecycle_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;ZHANG_R51_ENABLE=1")
    add_test(NAME zhang_r51_graph_runtime_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/graph_runtime_contract/graph_runtime_contract)
    set_tests_properties(zhang_r51_graph_runtime_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;ZHANG_R51_ENABLE=1")
    add_custom_target(zhang_r51_float_qc_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_r51_float_qc.py ${CMAKE_BINARY_DIR}
        DEPENDS pea
        VERBATIM)
    add_dependencies(unit_tests zhang_r51_float_qc_contract_build)
    add_test(NAME zhang_r51_float_qc_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/float_qc_contract/float_qc_contract)
    set_tests_properties(zhang_r51_float_qc_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
    add_custom_target(zhang_r51_float_parallel_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_r51_float_parallel.py ${CMAKE_BINARY_DIR}
            ${CMAKE_BINARY_DIR}/float_parallel_contract
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_r51_float_parallel_contract_build)
    add_test(NAME zhang_r51_float_parallel_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/float_parallel_contract/float_parallel_benchmark all 17 31)
    set_tests_properties(zhang_r51_float_parallel_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=4;OMP_NUM_THREADS=4;MKL_NUM_THREADS=1"
        TIMEOUT 60)
    add_custom_target(zhang_r51_ar_optimization_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_r51_ar_hotspots.py ${CMAKE_BINARY_DIR}
            ${CMAKE_BINARY_DIR}/ar_optimization_contract contracts
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_r51_ar_optimization_contract_build)
    add_custom_target(zhang_r51_product_closure_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_r51_product_closure_tests.py ${CMAKE_BINARY_DIR}
            ${CMAKE_BINARY_DIR}/product_closure_contract
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_r51_product_closure_contract_build)
    add_custom_target(zhang_r51_incremental_condition_contract_build
        COMMAND ${ZHANG_CONTRACT_PYTHON}
            ${CMAKE_CURRENT_SOURCE_DIR}/run_r51_incremental_condition.py ${CMAKE_BINARY_DIR}
            ${CMAKE_BINARY_DIR}/incremental_condition_contract
        DEPENDS pea VERBATIM)
    add_dependencies(unit_tests zhang_r51_incremental_condition_contract_build)
    add_test(NAME zhang_r51_incremental_condition_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/incremental_condition_contract/incremental_condition_contract 96 20 3)
    set_tests_properties(zhang_r51_incremental_condition_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1;ZHANG_R51_AR_THREADS=all"
        TIMEOUT 120)
    add_test(NAME zhang_r51_product_closure_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/product_closure_contract/product_closure_contract tests)
    set_tests_properties(zhang_r51_product_closure_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=4;OMP_NUM_THREADS=4;MKL_NUM_THREADS=1;ZHANG_R51_AR_THREADS=all;ZHANG_R51_RATIO_ONLY=1"
        TIMEOUT 180)
    add_test(NAME zhang_r51_ar_optimization_contract_tests
        COMMAND ${CMAKE_BINARY_DIR}/ar_optimization_contract/ar_optimization_contract)
    set_tests_properties(zhang_r51_ar_optimization_contract_tests PROPERTIES
        ENVIRONMENT "OPENBLAS_NUM_THREADS=4;OMP_NUM_THREADS=4;MKL_NUM_THREADS=1;ZHANG_R51_AR_THREADS=all"
        TIMEOUT 120)
endif()

add_executable(zhang_product_publication_boundary_tests unit/product_publication_boundary.cpp)
configure_unit_test(zhang_product_publication_boundary_tests)
add_dependencies(unit_tests zhang_product_publication_boundary_tests)
set_tests_properties(zhang_product_publication_boundary_tests PROPERTIES TIMEOUT 30)
