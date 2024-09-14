//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "DenseDiagonalizer.hpp"

#include <cmath>

#include "DenseLinearAlgebra.hpp"

#ifdef ENABLE_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

CDenseDiagonalizer::CDenseDiagonalizer()

    : _state(true)

    , _isSolved(false)
{
}

auto
CDenseDiagonalizer::diagonalize(const CDenseMatrix& matrix) -> void
{
    // copy matrix into temporary storage

    _matrix = CDenseMatrix(matrix);

    // determine dimensions of matrix

    auto ndim = _matrix.getNumberOfRows();

    // initialize eigenvalues and eigenvectors

    _eigenVectors = CDenseMatrix(ndim, ndim);

    _eigenValues = std::vector<double>(ndim);

    // set up pointers to matrices and vectors

    auto mat = _matrix.values();

    auto evecs = _eigenVectors.values();

    auto evals = _eigenValues.data();

    // temporary array for pivot data

    auto ndim_int32 = static_cast<int32_t>(ndim);

    std::vector<int32_t> idx_int32(2 * ndim_int32);

    // initialize number of eigenvalues

    int32_t nval_int32 = 0;

    // diagonalize matrix

    auto st = LAPACKE_dsyevr(LAPACK_ROW_MAJOR,
                             'V',
                             'A',
                             'U',
                             ndim_int32,
                             mat,
                             ndim_int32,
                             0.0,
                             0.0,
                             0,
                             0,
                             1.0e-13,
                             &nval_int32,
                             evals,
                             evecs,
                             ndim_int32,
                             idx_int32.data());

    // update state of diagonalizer

    _state = (st == 0);

    // set eigenvales and eigenvectors availabilty flag

    if (_state) _isSolved = true;
}

auto
CDenseDiagonalizer::getState() const -> bool
{
    return _state;
}

auto
CDenseDiagonalizer::getEigenVectors() const -> CDenseMatrix
{
    return _eigenVectors;
}

auto
CDenseDiagonalizer::getEigenValues() const -> std::vector<double>
{
    return _eigenValues;
}

auto
CDenseDiagonalizer::getInvertedSqrtMatrix() const -> CDenseMatrix
{
    if (_isSolved)
    {
        // set up temporary e^-1/2 vector

        auto evec = _eigenValues;

        // set up pointers and dimensions of e^-1/2 vector

        auto fvals = evec.data();

        auto ndim = evec.size();

        // compute e^-1/2 vector

        for (size_t i = 0; i < ndim; i++)
        {
            fvals[i] = 1.0 / std::sqrt(fvals[i]);
        }

        // construct A^-1/2 matrix

        auto mat = denblas::multDiagByAt(evec, _eigenVectors);

        return denblas::multAB(_eigenVectors, mat);
    }

    return CDenseMatrix();
}

auto
CDenseDiagonalizer::getInvertedMatrix() const -> CDenseMatrix
{
    if (_isSolved)
    {
        // set up temporary e^-1 vector

        auto evec = _eigenValues;

        // set up pointers and dimensions of e^-1 vector

        auto fvals = evec.data();

        auto ndim = evec.size();

        // compute e^-1 vector

        for (size_t i = 0; i < ndim; i++)
        {
            fvals[i] = 1.0 / fvals[i];
        }

        // construct A^-1 matrix

        auto mat = denblas::multDiagByAt(evec, _eigenVectors);

        return denblas::multAB(_eigenVectors, mat);
    }

    return CDenseMatrix();
}
