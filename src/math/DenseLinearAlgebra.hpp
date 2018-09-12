//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef DenseLinearAlgebra_hpp
#define DenseLinearAlgebra_hpp

#include "DenseMatrix.hpp"

namespace denblas { // denblas namespace
    
    /**
     Computes matrix multiplication: A * B.

     @param matrixA the matrix A.
     @param matrixB the matrix B
     @return the matrix A * B.
     */
    CDenseMatrix multAB(const CDenseMatrix& matrixA,
                        const CDenseMatrix& matrixB);
    
    /**
     Computes matrix multiplication: A * B^T.
     
     @param matrixA the matrix A.
     @param matrixB the matrix B
     @return the matrix A * B^T.
     */
    CDenseMatrix multABt(const CDenseMatrix& matrixA,
                         const CDenseMatrix& matrixB);
    
    /**
     Computes matrix multiplication: A^T * B.
     
     @param matrixA the matrix A.
     @param matrixB the matrix B
     @return the matrix A^T * B.
     */
    CDenseMatrix multAtB(const CDenseMatrix& matrixA,
                         const CDenseMatrix& matrixB);

    /**
     Computes matrix substraction: A - B.
     
     @param matrixA the matrix A.
     @param matrixB the matrix B
     @return the matrix A - B.
     */
    CDenseMatrix subAB(const CDenseMatrix& matrixA,
                       const CDenseMatrix& matrixB);
    
    /**
     Computes matrix addition: A + factor * B.
     
     @param matrixA the matrix A.
     @param matrixB the matrix B
     @param factor the scaling factor of matrix B.
     @return the matrix A +  factor * B.
     */
    CDenseMatrix addAB(const CDenseMatrix& matrixA,
                       const CDenseMatrix& matrixB,
                       const double        factor);
    
} // denblas namespace

#endif /* DenseLinearAlgebra_hpp */
