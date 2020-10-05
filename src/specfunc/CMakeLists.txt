target_sources(vlxobjs
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/BoysFunction.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
