target_sources(vlxobjs
  PRIVATE
    DenseDiagonalizer.cpp
    DenseLinearAlgebra.cpp
    DenseMatrix.cpp
    Matrices.cpp
    Matrix.cpp
    MatrixFunc.cpp
    SubMatrix.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
