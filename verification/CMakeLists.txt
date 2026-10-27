option(APRL_ENABLE_PHYSICS_ANALYSIS
  "Register optional non-gating particle/field energy-budget examples" OFF)

option(APRL_BUILD_WINDOW_FEASIBILITY
  "Build bounded moving-window research probes; no production window changes" OFF)
if(APRL_BUILD_WINDOW_FEASIBILITY)
  foreach(probe window_streaming pic_dependency window_analytic)
    add_executable(aprl_${probe}_probe ${probe}_probe.cpp)
    target_link_libraries(aprl_${probe}_probe PRIVATE direct_eb_core)
    target_compile_options(aprl_${probe}_probe PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>)
  endforeach()
endif()

add_executable(aprl_core_checks core_checks.cpp)
target_link_libraries(aprl_core_checks PRIVATE direct_eb_core)
target_compile_options(aprl_core_checks PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>
)

add_executable(aprl_output_checks output_checks.cpp)
target_link_libraries(aprl_output_checks PRIVATE direct_eb_core)
target_compile_options(aprl_output_checks PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>
)

add_executable(aprl_shutdown_watchdog_probe shutdown_watchdog_probe.cpp)
target_link_libraries(aprl_shutdown_watchdog_probe PRIVATE direct_eb_core)
target_compile_options(aprl_shutdown_watchdog_probe PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>
)

# Delay only a test output's H5Fclose to exercise the non-root abort/root
# close race without adding a fault-injection hook to production code.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  add_library(aprl_duration_close_delay SHARED duration_close_delay.cpp)
  target_include_directories(aprl_duration_close_delay PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(aprl_duration_close_delay PRIVATE
    MPI::MPI_CXX ${HDF5_LIBRARIES} ${CMAKE_DL_LIBS})
  target_compile_definitions(aprl_duration_close_delay PRIVATE
    OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX)
  target_compile_options(aprl_duration_close_delay PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/check_duration_close.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/check_duration_close.cmake" @ONLY)

  add_library(aprl_duration_barrier_fault SHARED duration_barrier_fault.cpp)
  target_include_directories(aprl_duration_barrier_fault PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(aprl_duration_barrier_fault PRIVATE
    MPI::MPI_CXX ${HDF5_LIBRARIES} ${CMAKE_DL_LIBS})
  target_compile_definitions(aprl_duration_barrier_fault PRIVATE
    OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX)
  target_compile_options(aprl_duration_barrier_fault PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/check_duration_barrier_fault.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/check_duration_barrier_fault.cmake" @ONLY)
endif()

add_executable(aprl_frozen_plane_scan frozen_plane_scan.cpp)
target_link_libraries(aprl_frozen_plane_scan PRIVATE direct_eb_core)
target_compile_options(aprl_frozen_plane_scan PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>
)

set(APRL_VERIFY_WORK "${CMAKE_CURRENT_BINARY_DIR}/work")

add_executable(aprl_threaded_fields threaded_fields.cpp)
target_link_libraries(aprl_threaded_fields PRIVATE direct_eb_core)
add_executable(aprl_mpi_overlap_checks mpi_overlap_checks.cpp)
target_link_libraries(aprl_mpi_overlap_checks PRIVATE direct_eb_core)
add_executable(aprl_detector_sampling_checks detector_sampling_checks.cpp)
target_link_libraries(aprl_detector_sampling_checks PRIVATE direct_eb_core)
target_compile_options(aprl_detector_sampling_checks PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>)
foreach(ranks 1 2 4)
  add_test(NAME required.detector.sampling-${ranks}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl_detector_sampling_checks>
      "${APRL_VERIFY_WORK}/detector-sampling-${ranks}" ${MPIEXEC_POSTFLAGS})
  set_tests_properties(required.detector.sampling-${ranks} PROPERTIES
    LABELS "required;functional;mpi;detector" TIMEOUT 120)
endforeach()
add_test(NAME required.parallel.fields COMMAND aprl_threaded_fields)
set_tests_properties(required.parallel.fields PROPERTIES
  LABELS "required;functional;parallel" TIMEOUT 120)
foreach(ranks 2 3)
  add_test(NAME required.parallel.halo-${ranks}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl_mpi_overlap_checks>
      ${MPIEXEC_POSTFLAGS})
  set_tests_properties(required.parallel.halo-${ranks} PROPERTIES
    LABELS "required;functional;parallel;mpi" TIMEOUT 120)
endforeach()
set(APRL_VERIFY_CARDS "${CMAKE_CURRENT_SOURCE_DIR}/cards")
foreach(directory
    smoke teardown detector mpi-one mpi-two ledger invalid-key invalid-stop
    invalid-watchdog
    physics-k0 physics-signal hdf5 sase auto-z-one auto-z-two
    auto-z-reference ballistic-injection ballistic-injection-two
    ballistic-injection-mixed ballistic-cpml-initial
    ballistic-cpml-initial-two)
  file(MAKE_DIRECTORY "${APRL_VERIFY_WORK}/${directory}")
endforeach()
file(MAKE_DIRECTORY
  "${APRL_VERIFY_WORK}/output-safety"
  "${APRL_VERIFY_WORK}/frozen-mpi-release"
  "${APRL_VERIFY_WORK}/postprocess/signal"
  "${APRL_VERIFY_WORK}/postprocess/baseline"
  "${APRL_VERIFY_WORK}/postprocess/ballistic")

foreach(card reconstruction field-analysis trajectory power closure)
  configure_file(
    "${APRL_VERIFY_CARDS}/postprocess-${card}.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-${card}.yaml" @ONLY)
endforeach()
configure_file(
  "${APRL_VERIFY_CARDS}/postprocess-reconstruction-protected.yaml.in"
  "${CMAKE_CURRENT_BINARY_DIR}/postprocess-reconstruction-protected.yaml"
  @ONLY)
configure_file(
  "${APRL_VERIFY_CARDS}/postprocess-ballistic-reconstruction.yaml.in"
  "${CMAKE_CURRENT_BINARY_DIR}/postprocess-ballistic-reconstruction.yaml"
  @ONLY)
configure_file(
  "${APRL_VERIFY_CARDS}/postprocess-ballistic-analysis.yaml.in"
  "${CMAKE_CURRENT_BINARY_DIR}/postprocess-ballistic-analysis.yaml"
  @ONLY)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/check_teardown_trace.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/check_teardown_trace.cmake"
  @ONLY)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/check_shutdown_watchdog.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/check_shutdown_watchdog.cmake"
  @ONLY)
configure_file(
  "${APRL_VERIFY_CARDS}/frozen-mpi-release.yaml.in"
  "${CMAKE_CURRENT_BINARY_DIR}/frozen-mpi-release.yaml"
  @ONLY)

function(aprl_required_test name)
  add_test(NAME ${name} COMMAND ${ARGN})
  set_tests_properties(${name} PROPERTIES
    LABELS "required;functional"
    TIMEOUT 60
  )
endfunction()

function(aprl_required_mpi_test name)
  add_test(NAME ${name} COMMAND ${ARGN})
  set_tests_properties(${name} PROPERTIES
    LABELS "required;functional;mpi"
    TIMEOUT 60
  )
endfunction()

add_test(NAME required.runtime.smoke
  COMMAND $<TARGET_FILE:aprl>
          "${APRL_VERIFY_CARDS}/runtime-smoke.yaml")
set_tests_properties(required.runtime.smoke PROPERTIES
  LABELS "required;smoke"
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/smoke"
  TIMEOUT 60
)

add_test(NAME required.runtime.teardown-trace
  COMMAND ${CMAKE_COMMAND}
          "-DAPRL_EXECUTABLE=$<TARGET_FILE:aprl>"
          "-DTRACE_CARD=${APRL_VERIFY_CARDS}/runtime-smoke.yaml"
          "-DTRACE_WORK=${APRL_VERIFY_WORK}/teardown"
          -P "${CMAKE_CURRENT_BINARY_DIR}/check_teardown_trace.cmake")
set_tests_properties(required.runtime.teardown-trace PROPERTIES
  LABELS "required;smoke;mpi"
  TIMEOUT 60
)

add_test(NAME required.runtime.shutdown-watchdog
  COMMAND ${CMAKE_COMMAND}
          "-DWATCHDOG_PROBE=$<TARGET_FILE:aprl_shutdown_watchdog_probe>"
          -P "${CMAKE_CURRENT_BINARY_DIR}/check_shutdown_watchdog.cmake")
set_tests_properties(required.runtime.shutdown-watchdog PROPERTIES
  LABELS "required;smoke;mpi"
  TIMEOUT 60
)

if(TARGET aprl_duration_close_delay)
  add_test(NAME required.runtime.duration-close
    COMMAND ${CMAKE_COMMAND}
      "-DAPRL_EXECUTABLE=$<TARGET_FILE:aprl>"
      "-DDURATION_DELAY_LIBRARY=$<TARGET_FILE:aprl_duration_close_delay>"
      "-DOUTPUT_CHECKER=$<TARGET_FILE:aprl_output_checks>"
      "-DDURATION_CARD=${APRL_VERIFY_CARDS}/duration-close.yaml"
      "-DDURATION_WORK=${APRL_VERIFY_WORK}/duration-close"
      -P "${CMAKE_CURRENT_BINARY_DIR}/check_duration_close.cmake")
  set_tests_properties(required.runtime.duration-close PROPERTIES
    LABELS "required;functional;mpi" TIMEOUT 65 PROCESSORS 4 RUN_SERIAL TRUE)

  add_test(NAME required.runtime.duration-barrier-watchdog
    COMMAND ${CMAKE_COMMAND}
      "-DAPRL_EXECUTABLE=$<TARGET_FILE:aprl>"
      "-DMPI_FAULT_LIBRARY=$<TARGET_FILE:aprl_duration_barrier_fault>"
      "-DOUTPUT_CHECKER=$<TARGET_FILE:aprl_output_checks>"
      "-DDURATION_CARD=${APRL_VERIFY_CARDS}/duration-close.yaml"
      "-DMPI_FAULT_WORK=${APRL_VERIFY_WORK}/duration-barrier-watchdog"
      -P "${CMAKE_CURRENT_BINARY_DIR}/check_duration_barrier_fault.cmake")
  set_tests_properties(required.runtime.duration-barrier-watchdog PROPERTIES
    LABELS "required;functional;mpi" TIMEOUT 50 PROCESSORS 4 RUN_SERIAL TRUE)
endif()

aprl_required_mpi_test(required.lorentz-elegant
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl_core_checks> lorentz-elegant
  ${MPIEXEC_POSTFLAGS})

aprl_required_test(required.boris-phase
  $<TARGET_FILE:aprl_core_checks> boris-phase)

add_test(NAME required.frozen-plane.free-drift
  COMMAND $<TARGET_FILE:aprl_frozen_plane_scan>
    "${APRL_VERIFY_WORK}/frozen-plane")
set_tests_properties(required.frozen-plane.free-drift PROPERTIES
  LABELS "required;functional" TIMEOUT 180)

aprl_required_mpi_test(required.frozen-plane.free-drift-mpi
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl_frozen_plane_scan>
  "${APRL_VERIFY_WORK}/frozen-plane-mpi2" ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.frozen-plane.free-drift-mpi PROPERTIES
  TIMEOUT 180)

aprl_required_mpi_test(required.frozen-plane.activation-slab-transfer
  ${CMAKE_COMMAND}
  "-DGENERATOR=$<TARGET_FILE:uniform_lab_plane_bunch_generator>"
  "-DAPRL=$<TARGET_FILE:aprl>"
  "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
  "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
  "-DCARD=${CMAKE_CURRENT_BINARY_DIR}/frozen-mpi-release.yaml"
  "-DINPUT=${APRL_VERIFY_WORK}/frozen-mpi-release/particles.h5"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/check_frozen_mpi_release.cmake")
set_tests_properties(required.frozen-plane.activation-slab-transfer PROPERTIES
  TIMEOUT 90)

aprl_required_test(required.longitudinal-layout
  $<TARGET_FILE:aprl_core_checks> longitudinal-layout)

aprl_required_test(required.runtime.defaults
  $<TARGET_FILE:aprl_core_checks> runtime-defaults)

aprl_required_test(required.initial-field.transverse-limit
  $<TARGET_FILE:aprl_core_checks> initial-field-transverse
  "${APRL_VERIFY_CARDS}/initial-field-transverse.yaml"
  "${APRL_VERIFY_CARDS}/initial-field-transverse-invalid.yaml")

aprl_required_test(required.initial-field.cpml-particles
  $<TARGET_FILE:aprl_core_checks> initial-field-cpml-particles)

aprl_required_test(required.charge-continuity
  $<TARGET_FILE:aprl_core_checks> charge-continuity)

aprl_required_test(required.cowan-dispersion
  $<TARGET_FILE:aprl_core_checks> cowan-dispersion)

aprl_required_mpi_test(required.current-filter
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl_core_checks> current-filter
  ${MPIEXEC_POSTFLAGS})

aprl_required_test(required.cpml-reflection
  $<TARGET_FILE:aprl_core_checks> cpml-reflection)

aprl_required_mpi_test(required.hdf5-versions
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl_core_checks> hdf5-versions
  "${APRL_VERIFY_WORK}/hdf5" ${MPIEXEC_POSTFLAGS})

if(APRL_BUILD_SASE_BUNCH_GENERATOR)
  aprl_required_test(required.sase-generator.self-test
    $<TARGET_FILE:sase_bunch_generator> --self-test)

  aprl_required_test(required.sase-generator.cleanup
    ${CMAKE_COMMAND} -E rm -f
    "${APRL_VERIFY_WORK}/sase/generated.h5")

  aprl_required_test(required.sase-generator.output
    $<TARGET_FILE:sase_bunch_generator>
    --output "${APRL_VERIFY_WORK}/sase/generated.h5"
    --electrons 1.0e8
    --length-m 8.0e-6
    --resonant-wavelength-m 1.0e-6
    --gamma 50.0
    --sigma-gamma 0.02
    --sigma-x-m 2.0e-5
    --sigma-y-m 3.0e-5
    --emittance-n-x-m 1.0e-6
    --emittance-n-y-m 1.5e-6
    --harmonics 3
    --beamlets-per-wavelength 8
    --seed 17)
  set_tests_properties(required.sase-generator.output PROPERTIES
    DEPENDS required.sase-generator.cleanup)

  aprl_required_mpi_test(required.sase-generator.aprl-read
    ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl_core_checks>
    sase-generator-file "${APRL_VERIFY_WORK}/sase/generated.h5"
    ${MPIEXEC_POSTFLAGS})
  set_tests_properties(required.sase-generator.aprl-read PROPERTIES
    DEPENDS required.sase-generator.output)

  aprl_required_test(required.sase-generator.refuses-overwrite
    $<TARGET_FILE:sase_bunch_generator>
    --output "${APRL_VERIFY_WORK}/sase/generated.h5"
    --electrons 1.0e8
    --length-m 8.0e-6
    --resonant-wavelength-m 1.0e-6
    --gamma 50.0
    --sigma-x-m 2.0e-5
    --sigma-y-m 3.0e-5)
  set_tests_properties(required.sase-generator.refuses-overwrite PROPERTIES
    DEPENDS required.sase-generator.output
    WILL_FAIL TRUE)
endif()

aprl_required_test(required.detector-stop.run
  $<TARGET_FILE:aprl> "${APRL_VERIFY_CARDS}/detector-stop.yaml")
set_tests_properties(required.detector-stop.run PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/detector")

aprl_required_test(required.detector-stop.check
  $<TARGET_FILE:aprl_output_checks> detector
  "${APRL_VERIFY_WORK}/detector/detectors/field.h5"
  "${APRL_VERIFY_WORK}/detector/detectors/particles.h5"
  4 2.0e-16 "${APRL_VERIFY_WORK}/detector/run-manifest.yaml")
set_tests_properties(required.detector-stop.check PROPERTIES
  DEPENDS required.detector-stop.run)

aprl_required_test(required.stop.invalid-empty-beamline
  $<TARGET_FILE:aprl> "${APRL_VERIFY_CARDS}/invalid-empty-stop.yaml")
set_tests_properties(required.stop.invalid-empty-beamline PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/invalid-stop"
  WILL_FAIL TRUE)

aprl_required_test(required.yaml.unknown-key
  $<TARGET_FILE:aprl> "${APRL_VERIFY_CARDS}/invalid-unknown-key.yaml")
set_tests_properties(required.yaml.unknown-key PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/invalid-key"
  WILL_FAIL TRUE)

aprl_required_test(required.runtime.shutdown-watchdog-invalid-timeout
  $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/invalid-shutdown-watchdog.yaml")
set_tests_properties(
  required.runtime.shutdown-watchdog-invalid-timeout PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/invalid-watchdog"
  WILL_FAIL TRUE)

aprl_required_test(required.runtime.shutdown-watchdog-invalid-finalize-timeout
  $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/invalid-shutdown-watchdog-finalize.yaml")
set_tests_properties(
  required.runtime.shutdown-watchdog-invalid-finalize-timeout PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/invalid-watchdog"
  WILL_FAIL TRUE)

aprl_required_test(required.output-safety.cleanup
  ${CMAKE_COMMAND} -E rm -f
  "${APRL_VERIFY_WORK}/output-safety/protected-manifest.yaml")

aprl_required_test(required.output-safety.first-run
  $<TARGET_FILE:aprl> "${APRL_VERIFY_CARDS}/output-safety.yaml")
set_tests_properties(required.output-safety.first-run PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/output-safety"
  DEPENDS required.output-safety.cleanup)

aprl_required_test(required.output-safety.refuses-overwrite
  $<TARGET_FILE:aprl> "${APRL_VERIFY_CARDS}/output-safety.yaml")
set_tests_properties(required.output-safety.refuses-overwrite PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/output-safety"
  DEPENDS required.output-safety.first-run
  WILL_FAIL TRUE)

aprl_required_mpi_test(required.mpi.one-rank
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/mpi-consistency.yaml" ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.mpi.one-rank PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/mpi-one")

aprl_required_mpi_test(required.ballistic-boundary-injection
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/ballistic-boundary-injection.yaml"
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.ballistic-boundary-injection PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/ballistic-injection"
  PASS_REGULAR_EXPRESSION "physical-aperture entries=[1-9]")

aprl_required_mpi_test(required.ballistic-boundary-injection-two-rank
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/ballistic-boundary-injection.yaml"
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.ballistic-boundary-injection-two-rank PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/ballistic-injection-two"
  PASS_REGULAR_EXPRESSION "physical-aperture entries=[1-9]")

# This check itself is serial, but its inputs are produced by the one- and
# two-rank MPI tests.  Exclude it from verify_required_no_mpi as a unit.
aprl_required_mpi_test(required.ballistic-boundary-ledger
  $<TARGET_FILE:aprl_output_checks> boundary-ledger
  "${APRL_VERIFY_WORK}/ballistic-injection/ledger/energy-ledger.h5"
  "${APRL_VERIFY_WORK}/ballistic-injection-two/ledger/energy-ledger.h5")
set_tests_properties(required.ballistic-boundary-ledger PROPERTIES
  DEPENDS "required.ballistic-boundary-injection;required.ballistic-boundary-injection-two-rank")

aprl_required_mpi_test(required.ballistic-boundary-mixed-poisson
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/ballistic-boundary-mixed.yaml"
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.ballistic-boundary-mixed-poisson PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/ballistic-injection-mixed"
  PASS_REGULAR_EXPRESSION "Boundary-injection initial Poisson source: [1-9]")

aprl_required_mpi_test(required.ballistic-boundary-cpml-initial
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/ballistic-boundary-cpml-initial.yaml"
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.ballistic-boundary-cpml-initial PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/ballistic-cpml-initial"
  PASS_REGULAR_EXPRESSION "including [1-9][0-9]* in the transverse CPML")

aprl_required_mpi_test(required.ballistic-boundary-cpml-initial-two-rank
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/ballistic-boundary-cpml-initial.yaml"
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.ballistic-boundary-cpml-initial-two-rank PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/ballistic-cpml-initial-two"
  PASS_REGULAR_EXPRESSION "including [1-9][0-9]* in the transverse CPML")

aprl_required_mpi_test(required.mpi.two-rank
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/mpi-consistency.yaml" ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.mpi.two-rank PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/mpi-two")

aprl_required_mpi_test(required.mpi.compare
  $<TARGET_FILE:aprl_output_checks> compare-particles
  "${APRL_VERIFY_WORK}/mpi-one/detectors/particles.h5"
  "${APRL_VERIFY_WORK}/mpi-two/detectors/particles.h5" 8)
set_tests_properties(required.mpi.compare PROPERTIES
  DEPENDS "required.mpi.one-rank;required.mpi.two-rank")

aprl_required_mpi_test(required.auto-z.one-rank
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/auto-longitudinal.yaml" ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.auto-z.one-rank PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/auto-z-one")

aprl_required_mpi_test(required.auto-z.two-rank
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/auto-longitudinal.yaml" ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.auto-z.two-rank PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/auto-z-two")

# This comparison consumes manifests produced by the one- and two-rank MPI
# runs above.  Mark it as MPI-coupled as well so --skip-mpi-run does not run
# the comparator after deliberately excluding both producers in a clean build.
aprl_required_mpi_test(required.auto-z.rank-independent
  ${CMAKE_COMMAND}
  -DONE_MANIFEST=${APRL_VERIFY_WORK}/auto-z-one/run-manifest.yaml
  -DTWO_MANIFEST=${APRL_VERIFY_WORK}/auto-z-two/run-manifest.yaml
  -P "${CMAKE_CURRENT_SOURCE_DIR}/check_auto_longitudinal.cmake")
set_tests_properties(required.auto-z.rank-independent PROPERTIES
  DEPENDS "required.auto-z.one-rank;required.auto-z.two-rank")

aprl_required_mpi_test(required.auto-z.reference-center
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:aprl>
  "${APRL_VERIFY_CARDS}/auto-longitudinal-reference.yaml"
  ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.auto-z.reference-center PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/auto-z-reference")

aprl_required_test(required.energy-ledger-k0.run
  $<TARGET_FILE:aprl> "${APRL_VERIFY_CARDS}/energy-ledger-k0.yaml")
set_tests_properties(required.energy-ledger-k0.run PROPERTIES
  WORKING_DIRECTORY "${APRL_VERIFY_WORK}/ledger")

aprl_required_test(required.energy-ledger-k0.check
  $<TARGET_FILE:aprl_output_checks> ledger
  "${APRL_VERIFY_WORK}/ledger/ledger/energy-ledger.h5")
set_tests_properties(required.energy-ledger-k0.check PROPERTIES
  DEPENDS required.energy-ledger-k0.run)

if(APRL_BUILD_POSTPROCESSORS)
  foreach(tool trajectory_radiation field_reconstruction
      field_plane_analysis field_power_compare energy_closure)
    aprl_required_test(required.yaml.unknown-key.${tool}
      $<TARGET_FILE:${tool}> "${APRL_VERIFY_CARDS}/postprocess-invalid.yaml")
    set_tests_properties(required.yaml.unknown-key.${tool} PROPERTIES
      WILL_FAIL TRUE)
  endforeach()

  aprl_required_test(required.postprocess.signal-run
    $<TARGET_FILE:aprl> "${APRL_VERIFY_CARDS}/postprocess-run.yaml")
  set_tests_properties(required.postprocess.signal-run PROPERTIES
    WORKING_DIRECTORY "${APRL_VERIFY_WORK}/postprocess/signal")

  aprl_required_test(required.postprocess.baseline-run
    $<TARGET_FILE:aprl> "${APRL_VERIFY_CARDS}/postprocess-run.yaml")
  set_tests_properties(required.postprocess.baseline-run PROPERTIES
    WORKING_DIRECTORY "${APRL_VERIFY_WORK}/postprocess/baseline")

  aprl_required_test(required.postprocess.reconstruction
    $<TARGET_FILE:field_reconstruction>
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-reconstruction.yaml")
  set_tests_properties(required.postprocess.reconstruction PROPERTIES
    DEPENDS required.postprocess.signal-run)

  aprl_required_test(required.postprocess.refuses-overwrite
    $<TARGET_FILE:field_reconstruction>
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-reconstruction-protected.yaml")
  set_tests_properties(required.postprocess.refuses-overwrite PROPERTIES
    DEPENDS required.postprocess.reconstruction
    WILL_FAIL TRUE)

  aprl_required_test(required.postprocess.field-analysis
    $<TARGET_FILE:field_plane_analysis>
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-field-analysis.yaml")
  set_tests_properties(required.postprocess.field-analysis PROPERTIES
    DEPENDS "required.postprocess.signal-run;required.postprocess.baseline-run")

  aprl_required_test(required.postprocess.ballistic-run
    $<TARGET_FILE:aprl>
    "${APRL_VERIFY_CARDS}/postprocess-ballistic-run.yaml")
  set_tests_properties(required.postprocess.ballistic-run PROPERTIES
    WORKING_DIRECTORY "${APRL_VERIFY_WORK}/postprocess/ballistic")

  aprl_required_test(required.postprocess.ballistic-reconstruction
    $<TARGET_FILE:field_reconstruction>
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-ballistic-reconstruction.yaml")
  set_tests_properties(required.postprocess.ballistic-reconstruction PROPERTIES
    DEPENDS required.postprocess.ballistic-run)

  aprl_required_test(required.postprocess.ballistic-field-analysis
    $<TARGET_FILE:field_plane_analysis>
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-ballistic-analysis.yaml")
  set_tests_properties(required.postprocess.ballistic-field-analysis PROPERTIES
    DEPENDS required.postprocess.ballistic-reconstruction)

  aprl_required_test(required.postprocess.field-routes-check
    $<TARGET_FILE:aprl_output_checks> field-routes
    "${APRL_VERIFY_WORK}/postprocess/field-analysis.h5"
    "${APRL_VERIFY_WORK}/postprocess/ballistic-reconstructed.h5"
    "${APRL_VERIFY_WORK}/postprocess/ballistic-field-analysis.h5")
  set_tests_properties(required.postprocess.field-routes-check PROPERTIES
    DEPENDS "required.postprocess.field-analysis;required.postprocess.ballistic-field-analysis")

  aprl_required_test(required.postprocess.trajectory-radiation
    $<TARGET_FILE:trajectory_radiation>
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-trajectory.yaml")
  set_tests_properties(required.postprocess.trajectory-radiation PROPERTIES
    DEPENDS required.postprocess.signal-run)

  aprl_required_test(required.postprocess.power-comparison
    $<TARGET_FILE:field_power_compare>
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-power.yaml")
  set_tests_properties(required.postprocess.power-comparison PROPERTIES
    DEPENDS "required.postprocess.signal-run;required.postprocess.baseline-run;required.postprocess.trajectory-radiation")

  aprl_required_test(required.postprocess.energy-closure
    $<TARGET_FILE:energy_closure>
    "${CMAKE_CURRENT_BINARY_DIR}/postprocess-closure.yaml")
  set_tests_properties(required.postprocess.energy-closure PROPERTIES
    DEPENDS "required.postprocess.signal-run;required.postprocess.field-analysis")

  aprl_required_test(required.postprocess.check
    $<TARGET_FILE:aprl_output_checks> postprocess
    "${APRL_VERIFY_WORK}/postprocess/reconstructed.h5"
    "${APRL_VERIFY_WORK}/postprocess/field-analysis.h5"
    "${APRL_VERIFY_WORK}/postprocess/trajectory-radiation.h5"
    "${APRL_VERIFY_WORK}/postprocess/power.h5"
    "${APRL_VERIFY_WORK}/postprocess/energy-closure.yaml")
  set_tests_properties(required.postprocess.check PROPERTIES
    DEPENDS "required.postprocess.reconstruction;required.postprocess.field-analysis;required.postprocess.trajectory-radiation;required.postprocess.power-comparison;required.postprocess.energy-closure")
endif()

if(APRL_ENABLE_PHYSICS_ANALYSIS)
  add_test(NAME analysis.energy-budget.k0
    COMMAND $<TARGET_FILE:aprl>
            "${APRL_VERIFY_CARDS}/physics-energy-k0.yaml")
  set_tests_properties(analysis.energy-budget.k0 PROPERTIES
    LABELS "physics;analysis"
    WORKING_DIRECTORY "${APRL_VERIFY_WORK}/physics-k0"
    TIMEOUT 120)

  add_test(NAME analysis.energy-budget.signal
    COMMAND $<TARGET_FILE:aprl>
            "${APRL_VERIFY_CARDS}/physics-energy-signal.yaml")
  set_tests_properties(analysis.energy-budget.signal PROPERTIES
    LABELS "physics;analysis"
    WORKING_DIRECTORY "${APRL_VERIFY_WORK}/physics-signal"
    TIMEOUT 120)

  add_test(NAME analysis.energy-budget.k0-report
    COMMAND $<TARGET_FILE:energy_ledger_report>
      "${APRL_VERIFY_WORK}/physics-k0/ledger/energy-ledger.h5")
  set_tests_properties(analysis.energy-budget.k0-report PROPERTIES
    LABELS "physics;analysis"
    DEPENDS analysis.energy-budget.k0)

  add_test(NAME analysis.energy-budget.signal-report
    COMMAND $<TARGET_FILE:energy_ledger_report>
      "${APRL_VERIFY_WORK}/physics-signal/ledger/energy-ledger.h5")
  set_tests_properties(analysis.energy-budget.signal-report PROPERTIES
    LABELS "physics;analysis"
    DEPENDS analysis.energy-budget.signal)
endif()

# Keep this full-pipeline check together: its two four-rank runs must consume
# exactly the same input, then compare all scientific datasets without h5py.
if(APRL_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
  add_executable(aprl_hybrid_output_checks hybrid_output_checks.cpp)
  target_include_directories(aprl_hybrid_output_checks PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(aprl_hybrid_output_checks PRIVATE MPI::MPI_CXX ${HDF5_LIBRARIES})
  target_compile_definitions(aprl_hybrid_output_checks PRIVATE
    OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX)
  configure_file("${APRL_VERIFY_CARDS}/hybrid-held.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/hybrid-held.yaml" @ONLY)
  add_test(NAME required.parallel.held-end-to-end
    COMMAND ${CMAKE_COMMAND}
      "-DGENERATOR=$<TARGET_FILE:uniform_lab_plane_bunch_generator>"
      "-DAPRL=$<TARGET_FILE:aprl>"
      "-DCOMPARATOR=$<TARGET_FILE:aprl_hybrid_output_checks>"
      "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
      "-DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
      "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
      "-DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}"
      "-DCARD=${CMAKE_CURRENT_BINARY_DIR}/hybrid-held.yaml"
      "-DWORK=${APRL_VERIFY_WORK}/hybrid-held"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_hybrid_held.cmake")
  set_tests_properties(required.parallel.held-end-to-end PROPERTIES
    LABELS "required;functional;parallel;mpi" TIMEOUT 420 PROCESSORS 8
    RUN_SERIAL TRUE)
endif()

add_custom_target(verify_required
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L required
  DEPENDS aprl aprl_core_checks aprl_output_checks
    uniform_lab_plane_bunch_generator
    aprl_shutdown_watchdog_probe energy_ledger_report
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  USES_TERMINAL
)

add_custom_target(verify_required_no_mpi
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L required -LE mpi
  DEPENDS aprl aprl_core_checks aprl_output_checks
    uniform_lab_plane_bunch_generator
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  USES_TERMINAL
  COMMENT "Running required non-MPI APRL verification tests"
)

if(APRL_BUILD_SASE_BUNCH_GENERATOR)
  add_dependencies(verify_required sase_bunch_generator)
  add_dependencies(verify_required_no_mpi sase_bunch_generator)
endif()

add_dependencies(verify_required aprl_threaded_fields aprl_mpi_overlap_checks)
add_dependencies(verify_required aprl_detector_sampling_checks)
if(TARGET aprl_duration_close_delay)
  add_dependencies(verify_required aprl_duration_close_delay)
  add_dependencies(verify_required aprl_duration_barrier_fault)
endif()
add_dependencies(verify_required_no_mpi aprl_threaded_fields)
if(TARGET aprl_hybrid_output_checks)
  add_dependencies(verify_required aprl_hybrid_output_checks)
endif()

if(APRL_BUILD_POSTPROCESSORS)
  add_dependencies(verify_required trajectory_radiation field_reconstruction
    field_plane_analysis field_power_compare energy_closure)
  add_dependencies(verify_required_no_mpi trajectory_radiation
    field_reconstruction field_plane_analysis field_power_compare
    energy_closure)
endif()

if(APRL_ENABLE_PHYSICS_ANALYSIS)
  add_custom_target(verify_physics
    COMMAND ${CMAKE_CTEST_COMMAND} --verbose --output-on-failure -L physics
    DEPENDS aprl energy_ledger_report
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    USES_TERMINAL
  )
endif()
