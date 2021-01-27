find_package(GTest CONFIG REQUIRED)

# Add targets mapping the same library names as defined in
# older versions of CMake's FindGTest
if(NOT TARGET GTest::GTest)
  add_library(GTest::GTest INTERFACE IMPORTED)
  target_link_libraries(GTest::GTest INTERFACE GTest::gtest)
endif()

add_executable(utests UnitTestVeloxChem.cpp)

target_link_libraries(utests
  PRIVATE
    GTest::GTest
    vlxobjs
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
