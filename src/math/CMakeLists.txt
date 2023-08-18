target_sources(vlxobjs
  PRIVATE
    DenseDiagonalizer.cpp
    DenseLinearAlgebra.cpp
    DenseMatrix.cpp
    Matrices.cpp
    MatrixFunc.cpp
    MatricesFunc.cpp	
    Matrix.cpp
    AngularMomentum.cpp
    SubMatrix.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
