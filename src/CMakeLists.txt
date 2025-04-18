# we compile everything as object files, then we link to the Python module and
# the unit test executable
add_library(vlxobjs OBJECT)

if(WIN32)
  target_compile_definitions(vlxobjs
    PUBLIC
      _USE_MATH_DEFINES
      _CRT_SECURE_NO_WARNINGS
      _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
    )
endif()

# flags to prepend
target_compile_options(vlxobjs
  BEFORE
  PRIVATE
    "${VLX_CXX_FLAGS}"
  PUBLIC
    "${ARCH_FLAG}"
  )

# flags to append
target_compile_options(vlxobjs
  PRIVATE
    "${EXTRA_CXXFLAGS}"
  )

# FIXME the includes should be relative to the source folder:
# that is: #include "dft/AOKohnShamMatrix.hpp", rather than: #include "AOKohnShamMatrix.hpp"
target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

target_link_libraries(vlxobjs
  PUBLIC
    OpenMP::OpenMP_CXX
    Math::LA
    Libxc::xc
  )

add_subdirectory(dft_func)
add_subdirectory(general)
add_subdirectory(math)
add_subdirectory(moldata)
add_subdirectory(onee_ints)
add_subdirectory(orbdata)
add_subdirectory(pdft_func)
add_subdirectory(pymodule)
add_subdirectory(python)
add_subdirectory(t2c_electric_dipole_momentum)
add_subdirectory(t2c_electric_dipole_momentum_geom_X00)
add_subdirectory(t2c_electron_repulsion)
add_subdirectory(t2c_electron_repulsion_geom_X00)
add_subdirectory(t2c_func)
add_subdirectory(t2c_geometrical_derivatives)
add_subdirectory(t2c_kinetic_energy)
add_subdirectory(t2c_kinetic_energy_geom_X00)
add_subdirectory(t2c_kinetic_energy_geom_X0Y)
add_subdirectory(t2c_nuclear_potential)
add_subdirectory(t2c_nuclear_potential_geom_0X0)
add_subdirectory(t2c_nuclear_potential_geom_X00)
add_subdirectory(t2c_nuclear_potential_geom_X0Y)
add_subdirectory(t2c_nuclear_potential_geom_XY0)
add_subdirectory(t2c_overlap)
add_subdirectory(t2c_overlap_geom_X00)
add_subdirectory(t2c_overlap_geom_X0Y)
add_subdirectory(t3c_electron_repulsion)
add_subdirectory(t3c_electron_repulsion_geom_0X0)
add_subdirectory(t3c_electron_repulsion_geom_X00)
add_subdirectory(t3c_func)
add_subdirectory(t3c_overlap)
add_subdirectory(t4c_electron_repulsion)
add_subdirectory(t4c_electron_repulsion_geom_X000)
add_subdirectory(t4c_electron_repulsion_geom_XY00)
add_subdirectory(t4c_electron_repulsion_geom_X0Y0)
add_subdirectory(t4c_func)
add_subdirectory(t4c_geometrical_derivatives)
add_subdirectory(visualization)
