add_library(vlxobjs OBJECT "")

target_compile_definitions(vlxobjs
  PUBLIC
    # TODO should probably be namespaced, like VLX_HAS_MPI or similar
    $<BUILD_INTERFACE:ENABLE_MPI=$<BOOL:$<TARGET_EXISTS:MPI::MPI_CXX>>>
    $<BUILD_INTERFACE:VLX_ALIGN=${VLX_ALIGN}>
    # TODO should probably be namespaced, like VLX_HAS_MKL or similar
    $<BUILD_INTERFACE:ENABLE_MKL=$<BOOL:$<TARGET_EXISTS:MKL>>>
  )

# NOTE we only add the OpenMP flag for the selected compiler, but do not link
# against the runtime explicitly.  This avoids conflicts between what the
# compiled and MKL think
target_compile_options(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${OpenMP_CXX_FLAGS}>
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
    MPI::MPI_CXX
    MKL
  )

set_target_properties(vlxobjs
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )

add_subdirectory(dft)
add_subdirectory(dftfuncs)
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
add_subdirectory(python)
add_subdirectory(response)
add_subdirectory(specfunc)
add_subdirectory(twoints)
add_subdirectory(visualization)

# TODO figure out how to handle the Python modules...
# add_subdirectory(pymodule)

# TODO Install rules
