target_sources(vlxobjs
  PRIVATE
    Dense4DTensor.cpp
    DenseDiagonalizer.cpp
    DenseLinearAlgebra.cpp
    DenseMatrix.cpp
    MathFunc.cpp
    SparseMatrix.cpp
    SparseVector.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
