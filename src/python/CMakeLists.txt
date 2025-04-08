# find mpi4py
include(${PROJECT_SOURCE_DIR}/cmake/custom/FindPythonModule.cmake)

find_python_module(mpi4py REQUIRED)

# create Python module with pybind11
# * use pybind11 header files with -isystem
# * explicitly disable LTO, even for release builds
pybind11_add_module(veloxchemlib
  NO_EXTRAS
    ExportDft.cpp
    ExportGeneral.cpp
    ExportMath.cpp
    ExportMoldata.cpp
    ExportOneElecInts.cpp
    ExportOrbdata.cpp
    ExportVisualization.cpp
    ExportT2CIntegrals.cpp
    ExportT3CIntegrals.cpp
    ExportT4CIntegrals.cpp
    ExportVeloxChemLib.cpp
  )

target_link_libraries(veloxchemlib
  PRIVATE
    vlxobjs
  )

file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${PYMOD_INSTALL_FULLDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(VLX_RPATH "@loader_path/${_rel}")
else()
  set(VLX_RPATH "\$ORIGIN/${_rel}")
endif()

set_target_properties(veloxchemlib
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${VLX_RPATH}${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH ON
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${PYMOD_INSTALL_FULLDIR}
  )

install(
  TARGETS
    veloxchemlib
  DESTINATION
    ${PYMOD_INSTALL_FULLDIR}
  )
