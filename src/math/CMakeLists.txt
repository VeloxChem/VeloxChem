target_sources(vlxobjs
  PRIVATE
    Dense4DTensor.cpp
    DenseMatrix.cpp
    Matrices.cpp
    Matrix.cpp
    MatrixFunc.cpp
    SerialDenseLinearAlgebra.cpp
    SubMatrix.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
