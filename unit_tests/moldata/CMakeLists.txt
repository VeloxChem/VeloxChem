target_sources(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/ChemicalElementTest.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MoleculeSetter.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MoleculeTest.cpp
  )

target_include_directories(utests
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )
