target_sources(utests
  PRIVATE
    ChemicalElementTest.cpp
    CommonNeighborsTest.cpp
    MoleculeSetter.cpp
    MoleculeTest.cpp
  )

target_include_directories(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
