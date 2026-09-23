target_sources(vlxobjs
  PRIVATE
    EmbeddingRegion.cpp
    PermanentMultipoles.cpp
    PolarizableEmbedding.cpp
    PolarizableSite.cpp
    PolarizableForceField.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
