target_sources(vlxobjs
  PRIVATE
    GtoValuesRecS.cpp
    GtoValuesRecP.cpp
    GtoValuesRecD.cpp
    GtoValuesRecF.cpp
    DftFunc.cpp 
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
