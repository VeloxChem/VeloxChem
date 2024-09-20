target_sources(vlxobjs
  PRIVATE
    AngularMomentumIntegrals.cpp
    GtoInfo.cpp
    LinearMomentumIntegrals.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
