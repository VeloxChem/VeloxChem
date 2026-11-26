target_sources(vlxobjs
  PRIVATE
    SimdTransformD.cpp
    SimdTransformF.cpp
    SimdTransformG.cpp
    SimdTransformH.cpp
    SimdTransformI.cpp
    SimdTransformK.cpp
    SimdTransformL.cpp
    SimdTransformP.cpp
    SimdTransformS.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
