# TODO: move dftd4 search function to cmake/custom/
message(STATUS "Searching for DFTD4...")

find_path(_dftd4_include_dirs
  NAMES
    dftd4.h
  PATHS
    $ENV{DFTD4_INCLUDE_DIR}
    $ENV{DFTD4_HOME}/include
  HINTS
    ${CMAKE_PREFIX_PATH}/include
  )

if(_eigen_include_dirs)
  find_library(_dftd4_libraries
    NAMES
      "dftd4"
    PATHS
      $ENV{DFTD4_LIB_DIR}
      $ENV{DFTD4_HOME}/lib
    HINTS
      ${CMAKE_PREFIX_PATH}/lib
    )

  if(_dftd4_libraries)
    message(STATUS "Found DFTD4: ${_dftd4_libraries}")
    target_compile_definitions(vlxobjs
      PUBLIC
        ENABLE_DFTD4
      )
    target_include_directories(vlxobjs
      SYSTEM
      PRIVATE
        ${_dftd4_include_dirs}
      )
    target_link_libraries(vlxobjs
      PUBLIC
        ${_dftd4_libraries}
      )
  endif()
endif()

if (NOT _dftd4_libraries)
  message(STATUS "Could not find DFTD4.")
endif()

target_sources(vlxobjs
  PRIVATE
    AtomicPartialChargesModel.cpp
    AtomicRadii.cpp
    ChemicalElement.cpp
    CoordinationNumber.cpp
    DispersionData.cpp
    DispersionModel.cpp
    DispersionParameters.cpp
    Molecule.cpp
    NewDispersionModel.cpp
    PartialCharges.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
