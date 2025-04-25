if(DEFINED ENV{EIGEN_INCLUDE_DIR})

  message(STATUS "Searching for Eigen in $ENV{EIGEN_INCLUDE_DIR}")

  find_path(_eigen_include_dirs
    NAMES
      Eigen/Dense
    PATHS
      $ENV{EIGEN_INCLUDE_DIR}
    )

  if(_eigen_include_dirs)
    message(STATUS "Found Eigen in $ENV{EIGEN_INCLUDE_DIR}")
    add_library(Math::LA INTERFACE IMPORTED)
    target_include_directories(Math::LA
      SYSTEM
      INTERFACE
        "${_eigen_include_dirs}"
      )
  else()
    message(STATUS "Could not find Eigen in $ENV{EIGEN_INCLUDE_DIR}")
    message(FATAL_ERROR "Please double check environment variable EIGEN_INCLUDE_DIR")
  endif()

else()

  message(STATUS "Searching for Eigen")

  find_path(_eigen_include_dirs
    NAMES
      eigen3/Eigen/Dense
    HINTS
      ${CMAKE_PREFIX_PATH}/include
    )

  if(_eigen_include_dirs)
    cmake_path(APPEND _eigen_include_dirs "eigen3")
    message(STATUS "Found Eigen in ${_eigen_include_dirs}")
    add_library(Math::LA INTERFACE IMPORTED)
    target_include_directories(Math::LA
      SYSTEM
      INTERFACE
        "${_eigen_include_dirs}"
      )
  else()
    message(STATUS "Could not find Eigen.")
    message(FATAL_ERROR "Please define environment variable EIGEN_INCLUDE_DIR")
  endif()

endif()
