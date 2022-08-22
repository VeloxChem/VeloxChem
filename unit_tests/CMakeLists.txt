find_package(GTest QUIET)
if(TARGET GTest::GTest)
  message(STATUS "Found GTest: ${GTEST_LIBRARIES}")
else()
  message(STATUS "Suitable Google Test could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(_gtest
    QUIET
    URL https://github.com/google/googletest/archive/release-1.12.1.tar.gz
  )
  # for windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  # prevent Google Test from using PThreads
  set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)
  # do not build Google Mock
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(_gtest)

  # Provide an alias, so linking to GTest looks the same regardless if it was
  # found on the system or if it was fetched at configuration
  add_library(GTest::GTest ALIAS gtest)
endif()

add_executable(utests UnitTestVeloxChem.cpp)

target_link_libraries(utests
  PRIVATE
    GTest::GTest
    vlxobjs
  )

add_subdirectory(dft)
add_subdirectory(eridrv)
add_subdirectory(general)
add_subdirectory(math)
add_subdirectory(moldata)
add_subdirectory(oneints)
add_subdirectory(oneints_amom)
add_subdirectory(oneints_dip)
add_subdirectory(oneints_ef)
add_subdirectory(oneints_kin)
add_subdirectory(oneints_lmom)
add_subdirectory(oneints_npot)
add_subdirectory(oneints_ovl)
add_subdirectory(orbdata)
add_subdirectory(response)
add_subdirectory(specfunc)
add_subdirectory(tools)
add_subdirectory(twoints)
add_subdirectory(visualization)

set_target_properties(utests
  PROPERTIES
    OUTPUT_NAME UnitTestVeloxChem.x
  )

include(GoogleTest)
gtest_discover_tests(utests
  PROPERTIES
    LABELS "unit"
  DISCOVERY_TIMEOUT  # how long to wait (in seconds) before crashing
    240
  DISCOVERY_MODE  # when to discover the tests
    PRE_TEST
  )