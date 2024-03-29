# we compile everything as object files, then we link to the Python module and
# the unit test executable
add_library(vlxobjs OBJECT)

target_compile_definitions(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:ENABLE_MPI=$<BOOL:$<TARGET_EXISTS:MPI::MPI_CXX>>>
    $<BUILD_INTERFACE:VLX_ALIGN=${VLX_ALIGN}>
  )

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
    MPI::MPI_CXX
    Math::LA
    Libxc::xc
  )

add_subdirectory(device)
add_subdirectory(dft)
add_subdirectory(dftfuncs)
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
add_subdirectory(pymodule)
add_subdirectory(python)
add_subdirectory(response)
add_subdirectory(specfunc)
add_subdirectory(twoints)
add_subdirectory(visualization)
add_subdirectory(xtb)
