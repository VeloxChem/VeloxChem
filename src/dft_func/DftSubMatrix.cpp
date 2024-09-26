//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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
getSubDensityMatrix(const double* densityPointer, const std::vector<int>& aoIndices, const int naos) -> CDenseMatrix
{
    const auto aocount = static_cast<int>(aoIndices.size());

    if (aocount <= naos)
    {
        CDenseMatrix sub_dens(aocount, aocount);

        for (int i = 0; i < aocount; i++)
        {
            auto sub_dens_row = sub_dens.row(i);

            auto dens_row = densityPointer + aoIndices[i] * naos;

            for (int j = 0; j < aocount; j++)
            {
                sub_dens_row[j] = dens_row[aoIndices[j]];
            }
        }

        return sub_dens;
    }

    return CDenseMatrix();
}

auto
distributeSubMatrixToKohnSham(CAOKohnShamMatrix& aoKohnShamMatrix, const CDenseMatrix& subMatrix, const std::vector<int>& aoIndices) -> void
{
    const auto naos = aoKohnShamMatrix.getNumberOfRows();

    const auto aocount = static_cast<int>(aoIndices.size());

    if (aocount <= naos)
    {
        for (int row = 0; row < subMatrix.getNumberOfRows(); row++)
        {
            auto row_orig = aoIndices[row];

            auto ksmat_row_orig = aoKohnShamMatrix.alphaValues() + row_orig * naos;

            auto submat_row = subMatrix.row(row);

            for (int col = 0; col < subMatrix.getNumberOfColumns(); col++)
            {
                auto col_orig = aoIndices[col];

                ksmat_row_orig[col_orig] += submat_row[col];
            }
        }
    }
}

auto
distributeSubMatrixToKohnSham(CAOKohnShamMatrix&               aoKohnShamMatrix,
                              const std::vector<CDenseMatrix>& subMatrices,
                              const std::vector<int>&          aoIndices) -> void
{
    const auto naos = aoKohnShamMatrix.getNumberOfRows();

    const auto aocount = static_cast<int>(aoIndices.size());

    if (aocount <= naos)
    {
        auto nrows = subMatrices[0].getNumberOfRows();

        auto ncols = subMatrices[0].getNumberOfColumns();

        for (int row = 0; row < nrows; row++)
        {
            auto row_orig = aoIndices[row];

            auto ksmat_a_row_orig = aoKohnShamMatrix.alphaValues() + row_orig * naos;

            auto ksmat_b_row_orig = aoKohnShamMatrix.betaValues() + row_orig * naos;

            auto submat_a_row = subMatrices[0].row(row);

            auto submat_b_row = subMatrices[1].row(row);

            for (int col = 0; col < ncols; col++)
            {
                auto col_orig = aoIndices[col];

                ksmat_a_row_orig[col_orig] += submat_a_row[col];

                ksmat_b_row_orig[col_orig] += submat_b_row[col];
            }
        }
    }
}

auto
distributeSubMatrixToFock(const std::vector<double*>& aoFockPointers,
                          const int                   fockIndex,
                          const CDenseMatrix&         subMatrix,
                          const std::vector<int>&     aoIndices,
                          const int                   naos) -> void
{
    const auto aocount = static_cast<int>(aoIndices.size());

    if (aocount <= naos)
    {
        for (int row = 0; row < subMatrix.getNumberOfRows(); row++)
        {
            auto row_orig = aoIndices[row];

            auto fock_row_orig = aoFockPointers[fockIndex] + row_orig * naos;

            auto submat_row = subMatrix.row(row);

            for (int col = 0; col < subMatrix.getNumberOfColumns(); col++)
            {
                auto col_orig = aoIndices[col];

                fock_row_orig[col_orig] += submat_row[col];
            }
        }
    }
}

}  // namespace dftsubmat
