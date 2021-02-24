target_sources(vlxobjs
  PRIVATE
    ChemicalElement.cpp
    Molecule.cpp
    PartialCharges.cpp
    VdwRadii.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
