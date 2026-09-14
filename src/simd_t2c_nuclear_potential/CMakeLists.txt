target_sources(vlxobjs
  PRIVATE
    SimdNuclearPotentialDriver.cpp
    SimdNuclearPotentialFunc.cpp
    SimdNuclearPotentialRecSS.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
