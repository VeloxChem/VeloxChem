#.rst:
#
# Sets the alignment
#
# Variables modified::
#
#   VLX_ALIGN

option_with_default(VLX_ALIGN "Alignment" 64)
list(APPEND _valid_align 32 64)
if(DEFINED VLX_ALIGN AND NOT VLX_ALIGN IN_LIST _valid_align)
  message(STATUS "${VLX_ALIGN} not a valid alignment, resetting to default value 64")
  set(VLX_ALIGN 64 CACHE STRING "Alignment" FORCE)
endif()
