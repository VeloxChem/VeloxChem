target_sources(vlxobjs
  PRIVATE
    ChemicalElement.cpp
    CoordinationNumber.cpp
    DispersionData.cpp
    DispersionModel.cpp
    DispersionParameters.cpp
    Molecule.cpp
    PartialCharges.cpp
    VdwRadii.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
