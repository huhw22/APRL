option(FEL_ENABLE_PHYSICS_ANALYSIS
  "Register optional non-gating particle/field energy-budget examples" OFF)

add_executable(fel_core_checks core_checks.cpp)
target_link_libraries(fel_core_checks PRIVATE direct_eb_core)
target_compile_options(fel_core_checks PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>
)

add_executable(fel_output_checks output_checks.cpp)
target_link_libraries(fel_output_checks PRIVATE direct_eb_core)
target_compile_options(fel_output_checks PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>
)

set(FEL_VERIFY_WORK "${CMAKE_CURRENT_BINARY_DIR}/work")
set(FEL_VERIFY_CARDS "${CMAKE_CURRENT_SOURCE_DIR}/cards")
foreach(directory
    smoke detector mpi-one mpi-two ledger invalid-key invalid-stop
    physics-k0 physics-signal hdf5)
  file(MAKE_DIRECTORY "${FEL_VERIFY_WORK}/${directory}")
endforeach()

function(fel_required_test name)
  add_test(NAME ${name} COMMAND ${ARGN})
  set_tests_properties(${name} PROPERTIES
    LABELS "required;functional"
    TIMEOUT 60
  )
endfunction()

add_test(NAME required.runtime.smoke
  COMMAND $<TARGET_FILE:simulator>
          "${FEL_VERIFY_CARDS}/runtime-smoke.yaml")
set_tests_properties(required.runtime.smoke PROPERTIES
  LABELS "required;smoke"
  WORKING_DIRECTORY "${FEL_VERIFY_WORK}/smoke"
  TIMEOUT 60
)

fel_required_test(required.lorentz-elegant
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:fel_core_checks> lorentz-elegant
  ${MPIEXEC_POSTFLAGS})

fel_required_test(required.boris-phase
  $<TARGET_FILE:fel_core_checks> boris-phase)

fel_required_test(required.charge-continuity
  $<TARGET_FILE:fel_core_checks> charge-continuity)

fel_required_test(required.cowan-dispersion
  $<TARGET_FILE:fel_core_checks> cowan-dispersion)

fel_required_test(required.cpml-reflection
  $<TARGET_FILE:fel_core_checks> cpml-reflection)

fel_required_test(required.hdf5-versions
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:fel_core_checks> hdf5-versions
  "${FEL_VERIFY_WORK}/hdf5" ${MPIEXEC_POSTFLAGS})

fel_required_test(required.detector-stop.run
  $<TARGET_FILE:simulator> "${FEL_VERIFY_CARDS}/detector-stop.yaml")
set_tests_properties(required.detector-stop.run PROPERTIES
  WORKING_DIRECTORY "${FEL_VERIFY_WORK}/detector")

fel_required_test(required.detector-stop.check
  $<TARGET_FILE:fel_output_checks> detector
  "${FEL_VERIFY_WORK}/detector/detectors/field.h5"
  "${FEL_VERIFY_WORK}/detector/detectors/particles.h5"
  4 2.0e-16)
set_tests_properties(required.detector-stop.check PROPERTIES
  DEPENDS required.detector-stop.run)

fel_required_test(required.stop.invalid-empty-beamline
  $<TARGET_FILE:simulator> "${FEL_VERIFY_CARDS}/invalid-empty-stop.yaml")
set_tests_properties(required.stop.invalid-empty-beamline PROPERTIES
  WORKING_DIRECTORY "${FEL_VERIFY_WORK}/invalid-stop"
  WILL_FAIL TRUE)

fel_required_test(required.yaml.unknown-key
  $<TARGET_FILE:simulator> "${FEL_VERIFY_CARDS}/invalid-unknown-key.yaml")
set_tests_properties(required.yaml.unknown-key PROPERTIES
  WORKING_DIRECTORY "${FEL_VERIFY_WORK}/invalid-key"
  WILL_FAIL TRUE)

fel_required_test(required.mpi.one-rank
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:simulator>
  "${FEL_VERIFY_CARDS}/mpi-consistency.yaml" ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.mpi.one-rank PROPERTIES
  WORKING_DIRECTORY "${FEL_VERIFY_WORK}/mpi-one")

fel_required_test(required.mpi.two-rank
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:simulator>
  "${FEL_VERIFY_CARDS}/mpi-consistency.yaml" ${MPIEXEC_POSTFLAGS})
set_tests_properties(required.mpi.two-rank PROPERTIES
  WORKING_DIRECTORY "${FEL_VERIFY_WORK}/mpi-two")

fel_required_test(required.mpi.compare
  $<TARGET_FILE:fel_output_checks> compare-particles
  "${FEL_VERIFY_WORK}/mpi-one/detectors/particles.h5"
  "${FEL_VERIFY_WORK}/mpi-two/detectors/particles.h5" 8)
set_tests_properties(required.mpi.compare PROPERTIES
  DEPENDS "required.mpi.one-rank;required.mpi.two-rank")

fel_required_test(required.energy-ledger-k0.run
  $<TARGET_FILE:simulator> "${FEL_VERIFY_CARDS}/energy-ledger-k0.yaml")
set_tests_properties(required.energy-ledger-k0.run PROPERTIES
  WORKING_DIRECTORY "${FEL_VERIFY_WORK}/ledger")

fel_required_test(required.energy-ledger-k0.check
  $<TARGET_FILE:fel_output_checks> ledger
  "${FEL_VERIFY_WORK}/ledger/ledger/energy-ledger.h5")
set_tests_properties(required.energy-ledger-k0.check PROPERTIES
  DEPENDS required.energy-ledger-k0.run)

if(FEL_ENABLE_PHYSICS_ANALYSIS)
  add_test(NAME analysis.energy-budget.k0
    COMMAND $<TARGET_FILE:simulator>
            "${FEL_VERIFY_CARDS}/physics-energy-k0.yaml")
  set_tests_properties(analysis.energy-budget.k0 PROPERTIES
    LABELS "physics;analysis"
    WORKING_DIRECTORY "${FEL_VERIFY_WORK}/physics-k0"
    TIMEOUT 120)

  add_test(NAME analysis.energy-budget.signal
    COMMAND $<TARGET_FILE:simulator>
            "${FEL_VERIFY_CARDS}/physics-energy-signal.yaml")
  set_tests_properties(analysis.energy-budget.signal PROPERTIES
    LABELS "physics;analysis"
    WORKING_DIRECTORY "${FEL_VERIFY_WORK}/physics-signal"
    TIMEOUT 120)

  add_test(NAME analysis.energy-budget.k0-report
    COMMAND $<TARGET_FILE:energy_ledger_report>
      "${FEL_VERIFY_WORK}/physics-k0/ledger/energy-ledger.h5")
  set_tests_properties(analysis.energy-budget.k0-report PROPERTIES
    LABELS "physics;analysis"
    DEPENDS analysis.energy-budget.k0)

  add_test(NAME analysis.energy-budget.signal-report
    COMMAND $<TARGET_FILE:energy_ledger_report>
      "${FEL_VERIFY_WORK}/physics-signal/ledger/energy-ledger.h5")
  set_tests_properties(analysis.energy-budget.signal-report PROPERTIES
    LABELS "physics;analysis"
    DEPENDS analysis.energy-budget.signal)
endif()

add_custom_target(verify_required
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -L required
  DEPENDS simulator fel_core_checks fel_output_checks energy_ledger_report
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  USES_TERMINAL
)

if(FEL_ENABLE_PHYSICS_ANALYSIS)
  add_custom_target(verify_physics
    COMMAND ${CMAKE_CTEST_COMMAND} --verbose --output-on-failure -L physics
    DEPENDS simulator energy_ledger_report
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    USES_TERMINAL
  )
endif()
