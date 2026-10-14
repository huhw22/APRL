if(NOT DEFINED ONE_MANIFEST OR NOT DEFINED TWO_MANIFEST)
  message(FATAL_ERROR "ONE_MANIFEST and TWO_MANIFEST are required")
endif()

foreach(manifest IN ITEMS "${ONE_MANIFEST}" "${TWO_MANIFEST}")
  if(NOT EXISTS "${manifest}")
    message(FATAL_ERROR "Automatic longitudinal manifest is missing: ${manifest}")
  endif()
endforeach()

file(READ "${ONE_MANIFEST}" one_text)
file(READ "${TWO_MANIFEST}" two_text)

set(cells_pattern
  "resolved_mesh:[\r\n]+  cells: \\[([0-9]+), ([0-9]+), ([0-9]+)\\]")
string(REGEX MATCH "${cells_pattern}" one_match "${one_text}")
if(NOT one_match)
  message(FATAL_ERROR "One-rank manifest has no resolved_mesh cells")
endif()
set(one_nz "${CMAKE_MATCH_3}")

string(REGEX MATCH "${cells_pattern}" two_match "${two_text}")
if(NOT two_match)
  message(FATAL_ERROR "Two-rank manifest has no resolved_mesh cells")
endif()
set(two_nz "${CMAKE_MATCH_3}")

if(one_nz LESS 3)
  message(FATAL_ERROR "Automatic longitudinal sizing resolved invalid Nz=${one_nz}")
endif()
if(NOT one_nz STREQUAL two_nz)
  message(FATAL_ERROR
    "Automatic longitudinal Nz depends on MPI ranks: one=${one_nz}, two=${two_nz}")
endif()

# The self-field-enabled card exercises the represented-charge mean/RMS
# reductions as well as the extrema.  Compare the complete derived placement
# fields, not only the rounded cell count, so changing MPI decomposition
# cannot silently move the physical interval.
foreach(field IN ITEMS
    origin_m outer_bounds_m physical_bounds_m translation_m
    initial_bounds_m ballistic_bounds_m radiation_upper_m
    self_field_bounds_m required_bounds_m minimum_physical_cells
    alignment_extra_cells right_rounding_padding_m)
  set(field_pattern "    ${field}: ([^\r\n]+)")
  string(REGEX MATCH "${field_pattern}" one_field_match "${one_text}")
  if(NOT one_field_match)
    message(FATAL_ERROR "One-rank manifest lacks ${field}")
  endif()
  set(one_field_value "${CMAKE_MATCH_1}")
  string(REGEX MATCH "${field_pattern}" two_field_match "${two_text}")
  if(NOT two_field_match)
    message(FATAL_ERROR "Two-rank manifest lacks ${field}")
  endif()
  set(two_field_value "${CMAKE_MATCH_1}")
  if(NOT one_field_value STREQUAL two_field_value)
    message(FATAL_ERROR
      "Automatic ${field} depends on MPI ranks: one=${one_field_value}; two=${two_field_value}")
  endif()
endforeach()

foreach(text IN ITEMS "${one_text}" "${two_text}")
  if(NOT text MATCHES "mode: observation-bounded")
    message(FATAL_ERROR
      "Automatic longitudinal manifest lacks observation-bounded provenance")
  endif()
  if(NOT text MATCHES "initial_self_field_enabled: true")
    message(FATAL_ERROR
      "Automatic longitudinal manifest lost initial-self-field provenance")
  endif()
  if(NOT text MATCHES "self_field_padding_rms: 0")
    message(FATAL_ERROR
      "Automatic longitudinal manifest lost box-clipped self-field padding")
  endif()
endforeach()

message(STATUS "Automatic longitudinal Nz is rank-independent: ${one_nz}")
