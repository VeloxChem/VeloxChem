target_sources(vlxobjs
  PRIVATE
    AtomicRadii.cpp
    ChemicalElement.cpp
    CommonNeighbors.cpp
    CoordinationNumber.cpp
    DispersionData.cpp
    DispersionModel.cpp
    DispersionParameters.cpp
    Molecule.cpp
    PartialCharges.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
