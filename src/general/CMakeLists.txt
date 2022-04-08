target_sources(vlxobjs
  PRIVATE
    Codata.cpp
    ErrorHandler.cpp
    FourIndexes.cpp
    MpiFunc.cpp
    OMPTasks.cpp
    StringFormat.cpp
    ThreeIndexes.cpp
    TwoIndexes.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
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
    ConfigMKL.hpp.in
    ConfigMKL.hpp
    @ONLY
    )
endif()
