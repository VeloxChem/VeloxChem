//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "DenseDiagonalizer.hpp"
#include "DenseLinearAlgebra.hpp"

#include <cmath>

#ifdef ENABLE_MKL
#include "mkl.h"
#else
#include "lapacke.h"
#endif

CDenseDiagonalizer::CDenseDiagonalizer()

    : _state(true)

    , _isSolved(false)
{
}

CDenseDiagonalizer::~CDenseDiagonalizer()
{
}

void
CDenseDiagonalizer::diagonalize(const CDenseMatrix& matrix)
{
    // copy matrix into temporary storage

    _matrix = matrix;

    // determine dimensions of matrix

    int32_t ndim = _matrix.getNumberOfRows();

    // initialize eigenvalues and eigenvectors

    _eigenVectors = CDenseMatrix(ndim);

    _eigenValues = CMemBlock<double>(ndim);

    // set up pointers to matrices and vectors

    auto mat = _matrix.values();

    auto evecs = _eigenVectors.values();

    auto evals = _eigenValues.data();

    // temporary array for pivot data

    CMemBlock<int32_t> idx(2 * ndim);

    // initialize number of eigenvalues

    int32_t nval = 0;

    // diagonalize matrix

    auto st = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', ndim, mat, ndim, 0.0, 0.0, 0, 0, 1.0e-13, &nval, evals, evecs, ndim, idx.data());

    // update state of diagonalizer

    _state = (st == 0);

    // set eigenvales and eigenvectors availabilty flag

    if (_state) _isSolved = true;
}

bool
CDenseDiagonalizer::getState() const
{
    return _state;
}

bool
CDenseDiagonalizer::isLinearlyDependentBasis(const double threshold) const
{
    if (getNumberOfEigenValues(threshold) != _eigenValues.size())
    {
        return true;
    }

    return false;
}

CDenseMatrix
CDenseDiagonalizer::getEigenVectors() const
{
    return _eigenVectors;
}

CDenseMatrix
CDenseDiagonalizer::getEigenVectors(const double threshold) const
{
    auto rdim = getNumberOfEigenValues(threshold);

    auto nrow = _eigenVectors.getNumberOfRows();

    if (rdim != nrow)
    {
        auto spos = nrow - rdim;

        return _eigenVectors.slice(0, spos, nrow, rdim);
    }

    return _eigenVectors;
}

CMemBlock<double>
CDenseDiagonalizer::getEigenValues() const
{
    return _eigenValues;
}

CMemBlock<double>
CDenseDiagonalizer::getEigenValues(const double threshold) const
{
    auto rdim = getNumberOfEigenValues(threshold);

    auto ndim = _eigenValues.size();

    if (rdim != ndim)
    {
        auto spos = ndim - rdim;

        return _eigenValues.slice(spos, rdim);
    }

    return _eigenValues;
}

CDenseMatrix
CDenseDiagonalizer::getInvertedSqrtMatrix() const
{
    if (_isSolved)
    {
        // set up temporary e^-1/2 vector

        auto evec = _eigenValues;

        // set up pointers and dimensions of e^-1/2 vector

        auto fvals = evec.data();

        auto ndim = evec.size();

        // compute e^-1/2 vector

        #pragma omp simd aligned(fvals: VLX_ALIGN)
        for (int32_t i = 0; i < ndim; i++)
        {
            fvals[i] = 1.0 / std::sqrt(fvals[i]);
        }

        // construct A^-1/2 matrix

        auto mat = denblas::multDiagByAt(evec, _eigenVectors);

        return denblas::multAB(_eigenVectors, mat);
    }

    return CDenseMatrix();
}

CDenseMatrix
CDenseDiagonalizer::getInvertedMatrix() const
{
    if (_isSolved)
    {
        // set up temporary e^-1 vector

        auto evec = _eigenValues;

        // set up pointers and dimensions of e^-1 vector

        auto fvals = evec.data();

        auto ndim = evec.size();

        // compute e^-1 vector

        #pragma omp simd aligned(fvals: VLX_ALIGN)
        for (int32_t i = 0; i < ndim; i++)
        {
            fvals[i] = 1.0 / fvals[i];
        }

        // construct A^-1 matrix

        auto mat = denblas::multDiagByAt(evec, _eigenVectors);

        return denblas::multAB(_eigenVectors, mat);
    }

    return CDenseMatrix();
}

int32_t
CDenseDiagonalizer::getNumberOfEigenValues(const double threshold) const
{
    // NOTE: dserv stores eigenvalues in ascending order

    if (_isSolved)
    {
        for (int32_t i = 0; i < _eigenValues.size(); i++)
        {
            if (_eigenValues.at(i) > threshold)
            {
                return _eigenValues.size() - i;
            }
        }
    }

    return 0;
}
