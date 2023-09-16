//
//                           VELOXCHEM 1.0-RC2
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

#include "DftSubMatrix.hpp"

#include "StringFormat.hpp"

namespace dftsubmat {  // dftsubmat namespace

auto
getSubDensityMatrix(const CDenseMatrix& densityMatrix, const std::vector<int64_t>& aoIndices) -> CDenseMatrix
{
    const auto naos = densityMatrix.getNumberOfRows();

    const auto aocount = static_cast<int64_t>(aoIndices.size());

    if (aocount <= naos)
    {
        CDenseMatrix sub_dens(aocount, aocount);

        for (int64_t i = 0; i < aocount; i++)
        {
            auto sub_dens_row = sub_dens.row(i);

            auto dens_row = densityMatrix.row(aoIndices[i]);

            for (int64_t j = 0; j < aocount; j++)
            {
                sub_dens_row[j] = dens_row[aoIndices[j]];
            }
        }

        return sub_dens;
    }

    return CDenseMatrix();
}

auto
distributeSubMatrixToDenseMatrix(CDenseMatrix& matrix, const CDenseMatrix& subMatrix, const std::vector<int64_t>& aoIndices) -> void
{
    const auto naos = matrix.getNumberOfRows();

    const auto aocount = static_cast<int64_t>(aoIndices.size());

    if (aocount <= naos)
    {
        for (int64_t row = 0; row < subMatrix.getNumberOfRows(); row++)
        {
            auto row_orig = aoIndices[row];

            auto mat_row_orig = matrix.row(row_orig);

            auto submat_row = subMatrix.row(row);

            for (int64_t col = 0; col < subMatrix.getNumberOfColumns(); col++)
            {
                auto col_orig = aoIndices[col];

                mat_row_orig[col_orig] += submat_row[col];
            }
        }
    }
}

}  // namespace dftsubmat
