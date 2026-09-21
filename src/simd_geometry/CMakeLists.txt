target_sources(vlxobjs
  PRIVATE
    SimdGeometryD1.cpp
    SimdGeometryF1.cpp
    SimdGeometryG1.cpp
    SimdGeometryH1.cpp
    SimdGeometryI1.cpp
    SimdGeometryK1.cpp
    SimdGeometryL1.cpp
    SimdGeometryP1.cpp
    SimdGeometryS1.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
