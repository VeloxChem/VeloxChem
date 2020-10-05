target_sources(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/ChemicalElement.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Molecule.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PartialCharges.cpp
    ${CMAKE_CURRENT_LIST_DIR}/VdwRadii.cpp
  )

target_include_directories(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )