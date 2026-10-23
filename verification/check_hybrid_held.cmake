foreach(required GENERATOR APRL COMPARATOR MPIEXEC MPIEXEC_NUMPROC_FLAG CARD WORK)
  if(NOT DEFINED ${required})
    message(FATAL_ERROR "Missing ${required} for hybrid held test")
  endif()
endforeach()
file(MAKE_DIRECTORY "${WORK}/threads-1" "${WORK}/threads-2")
execute_process(
  COMMAND "${GENERATOR}" --output "${WORK}/particles.h5"
    --current-A 0.05 --overwrite
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Cannot generate hybrid fixture: ${output}${error}")
endif()
foreach(threads 1 2)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
      "OMP_NUM_THREADS=${threads}" "OMP_DYNAMIC=FALSE" "APRL_PROFILE_INTERVAL=1"
      "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" 4 ${MPIEXEC_PREFLAGS}
      "${APRL}" ${MPIEXEC_POSTFLAGS} "${CARD}"
    WORKING_DIRECTORY "${WORK}/threads-${threads}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error
    TIMEOUT 180)
  file(WRITE "${WORK}/threads-${threads}/run.log" "${output}${error}")
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Hybrid ${threads}-thread run failed: ${output}${error}")
  endif()
  if(NOT output MATCHES "Compute threads per MPI rank: ${threads};")
    message(FATAL_ERROR "Requested compute threads were not used: ${output}")
  endif()
  string(REGEX MATCH "cross-rank activation pushes=([0-9]+)" crossed "${output}")
  if(NOT crossed OR CMAKE_MATCH_1 LESS 1)
    message(FATAL_ERROR "Hybrid run did not exercise cross-rank release: ${output}")
  endif()
  set(crossings "${CMAKE_MATCH_1}")
  string(REGEX MATCHALL "particles_max=[0-9]+" workloads "${output}")
  set(maximum 0)
  foreach(workload IN LISTS workloads)
    string(REPLACE "particles_max=" "" count "${workload}")
    if(count GREATER maximum)
      set(maximum "${count}")
    endif()
  endforeach()
  if(maximum LESS_EQUAL 2048)
    message(FATAL_ERROR "Hybrid fixture never exceeded the preparation batch: ${output}")
  endif()
  message(STATUS "Hybrid coverage: threads=${threads} ranks=4 max_local_particles=${maximum} cross_rank_release=${crossings}")
endforeach()
execute_process(COMMAND "${COMPARATOR}" "${WORK}/threads-1" "${WORK}/threads-2"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "Hybrid scientific outputs differ: ${output}${error}")
endif()
message(STATUS "${output}")
