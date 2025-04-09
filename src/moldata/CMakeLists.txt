target_sources(vlxobjs
  PRIVATE
    AtomicPartialChargesModel.cpp
    AtomicRadii.cpp
    ChemicalElement.cpp
    Molecule.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
