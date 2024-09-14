target_sources(vlxobjs
  PRIVATE
    Codata.cpp
    ErrorHandler.cpp
    MpiFunc.cpp
    OpenMPFunc.cpp
    StringFormat.cpp
    TensorLabels.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
