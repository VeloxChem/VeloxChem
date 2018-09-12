//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "DenseLinearAlgebra.hpp"

#include "mkl.h"

namespace denblas { // denblas namespace
    
    CDenseMatrix
    multAB(const CDenseMatrix& matrixA,
           const CDenseMatrix& matrixB)
    {
        // set up dimensions of matrix A
        
        auto narow = matrixA.getNumberOfRows();
        
        auto nacol = matrixA.getNumberOfColumns();
        
        // set up dimensions of matrix B
        
        auto nbcol = matrixB.getNumberOfColumns();
        
        // allocate dense matrix
        
        CDenseMatrix mat(narow, nbcol);
        
        // compute matrix-matrix multiplication
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    narow, nbcol, nacol,
                    1.0,
                    matrixA.values(), nacol,
                    matrixB.values(), nbcol,
                    0.0,
                    mat.values(), nbcol);
        
        return mat;
    }
    
    CDenseMatrix
    multABt(const CDenseMatrix& matrixA,
            const CDenseMatrix& matrixB)
    {
        // set up dimensions of matrix A
        
        auto narow = matrixA.getNumberOfRows();
        
        auto nacol = matrixA.getNumberOfColumns();
        
        // set up dimensions of matrix B
        
        auto nbrow = matrixB.getNumberOfRows();
        
        auto nbcol = matrixB.getNumberOfColumns();
        
        // allocate dense matrix
        
        CDenseMatrix mat(narow, nbrow);
        
        // compute matrix-matrix multiplcation
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    narow, nbrow, nacol,
                    1.0,
                    matrixA.values(), nacol,
                    matrixB.values(), nbcol,
                    0.0,
                    mat.values(), nbrow);
        
        return mat;
    }
    
    CDenseMatrix
    multAtB(const CDenseMatrix& matrixA,
            const CDenseMatrix& matrixB)
    {
        // set up dimensions of matrix A
        
        auto narow = matrixA.getNumberOfRows();
        
        auto nacol = matrixA.getNumberOfColumns();
        
        // set up dimensions of matrix B
        
        auto nbcol = matrixB.getNumberOfColumns();
        
        // allocate dense matrix
        
        CDenseMatrix mat(nacol, nbcol);

        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                    nacol, nbcol, narow,
                    1.0,
                    matrixA.values(), nacol,
                    matrixB.values(), nbcol,
                    0.0,
                    mat.values(), nbcol);

        return mat;
    }
    
    CDenseMatrix
    subAB(const CDenseMatrix& matrixA,
          const CDenseMatrix& matrixB)
    {
        // copy matrix
        
        CDenseMatrix mat = matrixA;
        
        // substract matrix
        
        cblas_daxpy(mat.getNumberOfElements(), -1.0, matrixB.values(), 1,
                    mat.values(), 1);
        
        return mat;
    }

    CDenseMatrix
    addAB(const CDenseMatrix& matrixA,
          const CDenseMatrix& matrixB,
          const double        factor)
    {
        // copy matrix
        
        CDenseMatrix mat = matrixA;
        
        // add scaled matrix
        
        cblas_daxpy(mat.getNumberOfElements(), factor, matrixB.values(), 1,
                    mat.values(), 1);
        
        return mat;
    }

    
} // denblas namespace
