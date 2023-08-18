target_sources(vlxobjs
  PRIVATE
    ChemicalElement.cpp
    CoordinationNumber.cpp
    Molecule.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
