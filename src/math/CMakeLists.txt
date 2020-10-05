target_sources(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/DenseDiagonalizer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/DenseLinearAlgebra.cpp
    ${CMAKE_CURRENT_LIST_DIR}/DenseMatrix.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MathFunc.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SparseMatrix.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SparseVector.cpp
  )

target_include_directories(veloxchemlib
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
  )