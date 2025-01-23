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
getSubAODensityMatrix(const std::vector<const double*>& densityPointers, const std::vector<int>& aoIndices, const int naos) -> CAODensityMatrix
{
    const auto aocount = static_cast<int>(aoIndices.size());

    if (aocount > naos) return CAODensityMatrix();

    std::vector<CDenseMatrix> submatrices;

    auto numdens = static_cast<int>(densityPointers.size());

    for (int idens = 0; idens < numdens; idens++)
    {
        CDenseMatrix sub_dens(aocount, aocount);

        auto dens = densityPointers[idens];

        for (int i = 0; i < aocount; i++)
        {
            auto sub_dens_row = sub_dens.row(i);

            auto dens_row = dens + aoIndices[i] * naos;

            for (int j = 0; j < aocount; j++)
            {
                sub_dens_row[j] = dens_row[aoIndices[j]];
            }
        }

        submatrices.push_back(sub_dens);
    }

    return CAODensityMatrix(submatrices, denmat::rest);
}

auto
getSubMatrixByColumnSlicing(const CDenseMatrix& denseMatrix, const std::vector<int>& aoIndices, const int naos) -> CDenseMatrix
{
    const auto aocount = static_cast<int>(aoIndices.size());

    const auto nrows = denseMatrix.getNumberOfRows();

    if ((aocount <= naos) && (nrows > 0))
    {
        CDenseMatrix sub_matrix(nrows, aocount);

        auto dense = denseMatrix.values();

        for (int i = 0; i < nrows; i++)
        {
            auto sub_matrix_row = sub_matrix.row(i);

            auto dense_row = dense + i * naos;

            for (int j = 0; j < aocount; j++)
            {
                sub_matrix_row[j] = dense_row[aoIndices[j]];
            }
        }

        return sub_matrix;
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
                              const CDenseMatrix&              subMatrix_a,
                              const CDenseMatrix&              subMatrix_b,
                              const std::vector<int>&          aoIndices) -> void
{
    const auto naos = aoKohnShamMatrix.getNumberOfRows();

    const auto aocount = static_cast<int>(aoIndices.size());

    if (aocount <= naos)
    {
        auto nrows = subMatrix_a.getNumberOfRows();

        auto ncols = subMatrix_a.getNumberOfColumns();

        for (int row = 0; row < nrows; row++)
        {
            auto row_orig = aoIndices[row];

            auto ksmat_a_row_orig = aoKohnShamMatrix.alphaValues() + row_orig * naos;

            auto ksmat_b_row_orig = aoKohnShamMatrix.betaValues() + row_orig * naos;

            auto submat_a_row = subMatrix_a.row(row);

            auto submat_b_row = subMatrix_b.row(row);

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
distributeSubMatrixToKohnSham(CAOKohnShamMatrix&               aoKohnShamMatrix,
                              const std::vector<CDenseMatrix>& subMatrices,
                              const std::vector<int>&          aoIndices) -> void
{
    distributeSubMatrixToKohnSham(aoKohnShamMatrix, subMatrices[0], subMatrices[1], aoIndices);
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

void
distributeSubMatrixToDenseMatrix(CDenseMatrix&           matrix,
                                 const CDenseMatrix&     subMatrix,
                                 const std::vector<int>& aoIndices,
                                 const int               naos)
{
    const auto aocount = static_cast<int>(aoIndices.size());

    if (aocount <= naos)
    {
        for (int row = 0; row < subMatrix.getNumberOfRows(); row++)
        {
            auto row_orig = aoIndices[row];

            auto mat_row_orig = matrix.row(row_orig);

            auto submat_row = subMatrix.row(row);

            for (int col = 0; col < subMatrix.getNumberOfColumns(); col++)
            {
                auto col_orig = aoIndices[col];

                mat_row_orig[col_orig] += submat_row[col];
            }
        }
    }
}

auto
distributeSubmatrixTo4DTensor(CDense4DTensor& fullTensor, const CDenseMatrix& subMatrix, const std::vector<int>& aoIndices) -> void
{
    const auto aocount = static_cast<int>(aoIndices.size());

    auto w_full = fullTensor.values();

    auto w_small = subMatrix.values();

    auto nAct3 = subMatrix.getNumberOfColumns();

    for (int i = 0; i < aocount; i++)
    {
        auto irow = nAct3 * i;

        auto irow_full = nAct3 * aoIndices[i];

        for (int jkl = 0; jkl < nAct3; jkl++)
        {
            w_full[irow_full + jkl] += w_small[irow + jkl];
        }
    }
}

}  // namespace dftsubmat
