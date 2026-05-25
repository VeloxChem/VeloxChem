if(DEFINED ENV{BOOST_INCLUDE_DIR})

  message(STATUS "Searching for Boost Math in $ENV{BOOST_INCLUDE_DIR}")

  find_path(_boost_include_dirs
    NAMES
      boost/math
    PATHS
      $ENV{BOOST_INCLUDE_DIR}
    )

  if(_boost_include_dirs)
    message(STATUS "Found Boost in $ENV{BOOST_INCLUDE_DIR}")
    add_library(Boost::headers INTERFACE IMPORTED)
    target_include_directories(Boost::headers
      SYSTEM
      INTERFACE
        "${_boost_include_dirs}"
      )
  else()
    message(STATUS "Could not find Boost in $ENV{BOOST_INCLUDE_DIR}")
    message(FATAL_ERROR "Please double check environment variable BOOST_INCLUDE_DIR")
  endif()

else()

  message(STATUS "Searching for Boost")

  find_path(_boost_include_dirs
    NAMES
      boost/math
    HINTS
      ${CMAKE_PREFIX_PATH}/include
    )

  if(_boost_include_dirs)
    message(STATUS "Found Boost in ${_boost_include_dirs}")
    add_library(Boost::headers INTERFACE IMPORTED)
    target_include_directories(Boost::headers
      SYSTEM
      INTERFACE
        "${_boost_include_dirs}"
      )
  else()
    message(STATUS "Could not find Boost.")
    message(FATAL_ERROR "Please define environment variable BOOST_INCLUDE_DIR")
  endif()

endif()
