find_package(GTest QUIET)
if(GTEST_FOUND)
  message(STATUS "Found GTest: ${GTEST_LIBRARIES}")
else()
  message(STATUS "Suitable Google Test could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(gtest_sources
    QUIET
    URL https://github.com/google/googletest/archive/release-1.10.0.tar.gz
  )

  FetchContent_GetProperties(gtest_sources)

  if(NOT gtest_sources_POPULATED)
    FetchContent_Populate(gtest_sources)

    # Prevent Google Test from using PThreads
    set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)

    # adds the targers: gtest, gtest_main, gmock, gmock_main
    add_subdirectory(
      ${gtest_sources_SOURCE_DIR}
      ${gtest_sources_BINARY_DIR}
      )
  endif()
endif()

add_executable(utests UnitTestVeloxChem.cpp)

target_link_libraries(utests
  PRIVATE
    GTest::GTest
    vlxobjs
    ${PYTHON_LIBRARIES}
  )

add_subdirectory(dft)
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
  )
