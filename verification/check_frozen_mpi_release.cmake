foreach(required GENERATOR APRL MPIEXEC MPIEXEC_NUMPROC_FLAG CARD INPUT)
  if(NOT DEFINED ${required})
    message(FATAL_ERROR "Missing ${required} for frozen MPI release test")
  endif()
endforeach()

execute_process(
  COMMAND "${GENERATOR}" --output "${INPUT}" --current-A 0.05 --overwrite
  RESULT_VARIABLE generate_result
  OUTPUT_VARIABLE generate_output
  ERROR_VARIABLE generate_error)
if(NOT generate_result EQUAL 0)
  message(FATAL_ERROR
    "Synthetic MPI release fixture failed: ${generate_output}${generate_error}")
endif()

execute_process(
  COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" 4 "${APRL}" "${CARD}"
  RESULT_VARIABLE run_result
  OUTPUT_VARIABLE run_output
  ERROR_VARIABLE run_error)
if(NOT run_result EQUAL 0)
  message(FATAL_ERROR
    "Frozen MPI release simulation failed: ${run_output}${run_error}")
endif()

string(REGEX MATCH
  "cross-rank activation pushes=([0-9]+)"
  cross_rank_count "${run_output}")
if(NOT cross_rank_count OR CMAKE_MATCH_1 LESS 1)
  message(FATAL_ERROR
    "Fixture did not exercise release across a z slab: ${run_output}")
endif()
message(STATUS "Frozen MPI release crossed slabs for ${CMAKE_MATCH_1} particles")
