//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#include "OverlapMatrix.hpp"

#include <cmath>

#include "DenseDiagonalizer.hpp"
#include "StringFormat.hpp"

COverlapMatrix::COverlapMatrix()
{
}

COverlapMatrix::COverlapMatrix(const CDenseMatrix& matrix)

    : _matrix(matrix)
{
}

COverlapMatrix::COverlapMatrix(const COverlapMatrix& source)

    : _matrix(source._matrix)
{
}

COverlapMatrix::COverlapMatrix(COverlapMatrix&& source) noexcept

    : _matrix(std::move(source._matrix))
{
}

COverlapMatrix::~COverlapMatrix()
{
}

COverlapMatrix&
COverlapMatrix::operator=(const COverlapMatrix& source)
{
    if (this == &source) return *this;

    _matrix = source._matrix;

    return *this;
}

COverlapMatrix&
COverlapMatrix::operator=(COverlapMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _matrix = std::move(source._matrix);

    return *this;
}

bool
COverlapMatrix::operator==(const COverlapMatrix& other) const
{
    if (_matrix != other._matrix) return false;

    return true;
}

bool
COverlapMatrix::operator!=(const COverlapMatrix& other) const
{
    return !(*this == other);
}

std::string
COverlapMatrix::getString() const
{
    return _matrix.getString();
}

int32_t
COverlapMatrix::getNumberOfRows() const
{
    return _matrix.getNumberOfRows();
}

int32_t
COverlapMatrix::getNumberOfColumns() const
{
    return _matrix.getNumberOfColumns();
}

int32_t
COverlapMatrix::getNumberOfElements() const
{
    return _matrix.getNumberOfElements();
}

const double*
COverlapMatrix::values() const
{
    return _matrix.values();
}

CDenseMatrix
COverlapMatrix::getOrthogonalizationMatrix(const double threshold) const
{
    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(_matrix);

    auto omat = diagdrv.getEigenVectors(threshold);

    auto eigs = diagdrv.getEigenValues(threshold);

    // set up dimensions

    auto ncol = eigs.size();

    auto nrow = omat.getNumberOfRows();

    // loop over matrix columns

    auto mdat = omat.values();

    for (int32_t i = 0; i < ncol; i++)
    {
        auto fact = 1.0 / std::sqrt(eigs.at(i));

        for (int32_t j = 0; j < nrow; j++)
        {
            mdat[j * ncol + i] *= fact;
        }
    }

    return omat;
}

std::ostream&
operator<<(std::ostream& output, const COverlapMatrix& source)
{
    output << std::endl;

    output << "[COverlapMatrix (Object):" << &source << "]" << std::endl;

    output << "_matrix: " << source._matrix << std::endl;

    return output;
}
