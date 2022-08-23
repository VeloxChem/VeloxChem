find_package(GTest QUIET)
if(NOT TARGET GTest::GTest)
  message(STATUS "Suitable Google Test could not be located. Unit tests target won't be available!")
else()
  message(STATUS "Found GTest: ${GTEST_LIBRARIES}")
  file(RELATIVE_PATH _build_dir ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
  message(STATUS "Compile unit tests with: cmake --build ${_build_dir} --target utests")

  add_executable(utests 
    EXCLUDE_FROM_ALL
    UnitTestVeloxChem.cpp
    )
  
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
endif()

