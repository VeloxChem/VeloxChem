target_sources(vlxobjs
  PRIVATE
    DftFunc.cpp 
    GtoValues.cpp
    GtoValuesRecS.cpp
    GtoValuesRecP.cpp
    GtoValuesRecD.cpp
    GtoValuesRecF.cpp
    Prescreener.cpp 
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
