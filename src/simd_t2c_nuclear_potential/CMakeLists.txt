target_sources(vlxobjs
  PRIVATE
    SimdNuclearPotentialCtrVrrSS.cpp
    SimdNuclearPotentialDriver.cpp
    SimdNuclearPotentialFunc.cpp
    SimdNuclearPotentialRecPP.cpp
    SimdNuclearPotentialRecPS.cpp
    SimdNuclearPotentialRecSP.cpp
    SimdNuclearPotentialRecSS.cpp
    SimdNuclearPotentialVrrRecPS.cpp
    SimdNuclearPotentialVrrRecSD.cpp
    SimdNuclearPotentialVrrRecSP.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
