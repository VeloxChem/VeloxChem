//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DenseLinearAlgebra.hpp"
#include "ErrorHandler.hpp"

#ifdef ENABLE_MKL
#include "mkl.h"
#else
#include "cblas.h"
#endif

namespace denblas {  // denblas namespace

CDenseMatrix
multAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB)
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();

    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();

    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(nacol == nbrow, "denblas::multAB - Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(narow, nbcol);

    // compute matrix-matrix multiplication

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                narow,
                nbcol,
                nacol,
                1.0,
                matrixA.values(),
                nacol,
                matrixB.values(),
                nbcol,
                0.0,
                mat.values(),
                nbcol);

    return mat;
}

CDenseMatrix
multABt(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB)
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();

    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();

    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(nacol == nbcol, "denblas::multABt - Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(narow, nbrow);

    // compute matrix-matrix multiplcation

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasTrans,
                narow,
                nbrow,
                nacol,
                1.0,
                matrixA.values(),
                nacol,
                matrixB.values(),
                nbcol,
                0.0,
                mat.values(),
                nbrow);

    return mat;
}

CDenseMatrix
multAtB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB)
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();

    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();

    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(narow == nbrow, "denblas::multAtB - Inconsistent sizes in matrix multiplication");

    // allocate dense matrix

    CDenseMatrix mat(nacol, nbcol);

    cblas_dgemm(CblasRowMajor,
                CblasTrans,
                CblasNoTrans,
                nacol,
                nbcol,
                narow,
                1.0,
                matrixA.values(),
                nacol,
                matrixB.values(),
                nbcol,
                0.0,
                mat.values(),
                nbcol);

    return mat;
}

CDenseMatrix
multDiagByA(const CMemBlock<double>& diagonal, const CDenseMatrix& matrix)
{
    // set up dimensions of matrix

    auto nrow = matrix.getNumberOfRows();

    auto ncol = matrix.getNumberOfColumns();

    // allocate results matrix

    CDenseMatrix mat(nrow, ncol);

    // set up pointers to matrices

    auto mval = mat.values();

    auto sval = matrix.values();

    auto dval = diagonal.data();

    // compute matrix multiplication

    for (int32_t i = 0; i < nrow; i++)
    {
        // set up local pointers to rows

        auto cmval = &mval[i * ncol];

        auto csval = &sval[i * ncol];

        // fetch value of diagonal

        auto f = dval[i];

        #pragma omp simd
        for (int32_t j = 0; j < ncol; j++)
        {
            cmval[j] = f * csval[j];
        }
    }

    return mat;
}

CDenseMatrix
multDiagByAt(const CMemBlock<double>& diagonal, const CDenseMatrix& matrix)
{
    // set up dimensions of matrix

    auto nrow = matrix.getNumberOfRows();

    auto ncol = matrix.getNumberOfColumns();

    // allocate results matrix

    CDenseMatrix mat(ncol, nrow);

    // set up pointers to matrices

    auto mval = mat.values();

    auto sval = matrix.values();

    auto dval = diagonal.data();

    // compute matrix multiplication

    for (int32_t i = 0; i < ncol; i++)
    {
        auto ioff = i * nrow;

        for (int32_t j = 0; j < nrow; j++)
        {
            mval[ioff + j] = dval[i] * sval[j * ncol + i];
        }
    }

    return mat;
}

CDenseMatrix
subAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB)
{
    errors::assertMsgCritical(matrixA.getNumberOfElements() == matrixB.getNumberOfElements(),
                              "denblas::subAB - Inconsistent sizes in matrix subtraction");

    // copy matrix

    CDenseMatrix mat = matrixA;

    // substract matrix

    cblas_daxpy(mat.getNumberOfElements(), -1.0, matrixB.values(), 1, mat.values(), 1);

    return mat;
}

CDenseMatrix
addAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB, const double factor)
{
    errors::assertMsgCritical(matrixA.getNumberOfElements() == matrixB.getNumberOfElements(),
                              "denblas::addAB - Inconsistent sizes in matrix addition");

    // copy matrix

    CDenseMatrix mat = matrixA;

    // add scaled matrix

    cblas_daxpy(mat.getNumberOfElements(), factor, matrixB.values(), 1, mat.values(), 1);

    return mat;
}

void
multABt(CDenseMatrix& matrixC, const double alpha, const CDenseMatrix& matrixA, const CDenseMatrix& matrixB)
{
    // set up dimensions of matrix A

    auto narow = matrixA.getNumberOfRows();

    auto nacol = matrixA.getNumberOfColumns();

    // set up dimensions of matrix B

    auto nbrow = matrixB.getNumberOfRows();

    auto nbcol = matrixB.getNumberOfColumns();

    errors::assertMsgCritical(nacol == nbcol, "denblas::multABt - Inconsistent sizes in matrix multiplication");

    errors::assertMsgCritical(narow == matrixC.getNumberOfRows(), "denblas::multABt - Inconsistent sizes in matrix multiplication");

    errors::assertMsgCritical(nbrow == matrixC.getNumberOfColumns(), "denblas::multABt - Inconsistent sizes in matrix multiplication");

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasTrans,
                narow,
                nbrow,
                nacol,
                alpha,
                matrixA.values(),
                nacol,
                matrixB.values(),
                nbcol,
                1.0,
                matrixC.values(),
                nbrow);
}

void
multAtB(double* matrixC, const double alpha, const CDenseMatrix& matrixA, const CDenseMatrix& matrixB)
{
    // set up dimensions of matrix A
    
    auto narow = matrixA.getNumberOfRows();
    
    auto nacol = matrixA.getNumberOfColumns();
    
    // set up dimensions of matrix B
    
    auto nbrow = matrixB.getNumberOfRows();
    
    auto nbcol = matrixB.getNumberOfColumns();
    
    errors::assertMsgCritical(narow == nbrow, "denblas::multAtB - Inconsistent sizes in matrix multiplication");
    
    // allocate dense matrix
    
    CDenseMatrix mat(nacol, nbcol);
    
    cblas_dgemm(CblasRowMajor,
                CblasTrans,
                CblasNoTrans,
                nacol,
                nbcol,
                narow,
                alpha,
                matrixA.values(),
                nacol,
                matrixB.values(),
                nbcol,
                1.0,
                matrixC,
                nbcol);
}
    
double
dot(const CMemBlock<double>& vectorA, const CMemBlock<double>& vectorB)
{
    errors::assertMsgCritical(vectorA.size() == vectorB.size(), "denblas::dot - Inconsistent sizes in dot product of vectors");

    return cblas_ddot(vectorA.size(), vectorA.data(), 1, vectorB.data(), 1);
}

double
dot(const CMemBlock<double>& vectorA, const CDenseMatrix& matrixB)
{
    errors::assertMsgCritical(vectorA.size() == matrixB.getNumberOfElements(),
                              "denblas::dot - Inconsistent sizes in dot product of vector and column matrix");

    return cblas_ddot(vectorA.size(), vectorA.data(), 1, matrixB.values(), 1);
}

double
trace(const CDenseMatrix& matrix)
{
    errors::assertMsgCritical(matrix.getNumberOfColumns() == matrix.getNumberOfRows(), "denblas::trace - Non-square matrix");

    auto pvals = matrix.values();

    auto mdim = matrix.getNumberOfRows();

    double fsum = 0.0;

    for (int32_t i = 0; i < mdim; i++)
        fsum += pvals[i * mdim + i];

    return fsum;
}

double
trace(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB)
{
    return denblas::trace(denblas::multAB(matrixA, matrixB));
}

}  // namespace denblas
