target_sources(vlxobjs
  PRIVATE
    ErrorHandler.cpp
    MultiTimer.cpp
    StringFormat.cpp
    OpenMPFunc.cpp
    MpiFunc.cpp
    Timer.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )
