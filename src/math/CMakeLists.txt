target_sources(vlxobjs
  PRIVATE
    GeomThreeMatrices.cpp
    GeomFourMatrices.cpp
    Matrices.cpp
    MatrixFunc.cpp
    Matrix.cpp
    AngularMomentum.cpp
    GeomTwoMatrices.cpp
    SubMatrix.cpp
  )

target_include_directories(vlxobjs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )
