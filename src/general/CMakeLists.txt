target_sources(vlxobjs
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/Codata.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ErrorHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FourIndexes.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MemAlloc.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MpiFunc.cpp
    ${CMAKE_CURRENT_LIST_DIR}/OMPTasks.cpp
    ${CMAKE_CURRENT_LIST_DIR}/StringFormat.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ThreeIndexes.cpp
    ${CMAKE_CURRENT_LIST_DIR}/TwoIndexes.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
