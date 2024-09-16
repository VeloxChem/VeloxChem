target_sources(vlxobjs
  PRIVATE
    DftFunc.cpp
    GtoValues.cpp
    GtoValuesRecD.cpp
    GtoValuesRecP.cpp
    GtoValuesRecS.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
