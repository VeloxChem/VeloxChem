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

if(VLX_LA_VENDOR MATCHES MKL)
  if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(_mkl_interface_layer_ "MKL_INTERFACE_LP64+MKL_INTERFACE_GNU")
    set(_mkl_threading_layer_ "MKL_THREADING_GNU")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(_mkl_interface_layer_ "MKL_INTERFACE_LP64")
    set(_mkl_threading_layer_ "MKL_THREADING_INTEL")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(_mkl_interface_layer_ "MKL_INTERFACE_LP64")
    set(_mkl_threading_layer_ "MKL_THREADING_INTEL")
  endif()

  # generate MKL configuration header file
  configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/ConfigMKL.hpp.in
    ${CMAKE_CURRENT_LIST_DIR}/ConfigMKL.hpp
    @ONLY
    )
endif()
