target_sources(vlxobjs
  PRIVATE
    Codata.cpp
    ErrorHandler.cpp
    MpiFunc.cpp
    MultiTimer.cpp
    OpenMPFunc.cpp
    StringFormat.cpp
    TensorLabels.cpp
    Timer.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
