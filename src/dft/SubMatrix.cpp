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

#include "SubMatrix.hpp"
#include "StringFormat.hpp"

namespace submat {  // submat namespace

CDenseMatrix
getSubDensityMatrix(const CAODensityMatrix&     densityMatrix,
                    const int32_t               densityIndex,
                    const std::string&          densitySpin,
                    const std::vector<int32_t>& aoIndices,
                    const int32_t               aoCount,
                    const int32_t               nAOs)
{
    if (aoCount <= nAOs)
    {
        CDenseMatrix sub_dens(aoCount, aoCount);

        bool alphaspin = (fstr::upcase(densitySpin) == std::string("ALPHA"));

        auto dens = alphaspin ? densityMatrix.alphaDensity(densityIndex) : densityMatrix.betaDensity(densityIndex);

        for (int32_t i = 0; i < aoCount; i++)
        {
            auto sub_dens_row = sub_dens.row(i);

            auto dens_row = dens + aoIndices[i] * nAOs;

            for (int32_t j = 0; j < aoCount; j++)
            {
                sub_dens_row[j] = dens_row[aoIndices[j]];
            }
        }

        return sub_dens;
    }

    return CDenseMatrix();
}

CDenseMatrix
getSubMatrixByColumnSlicing(const CDenseMatrix&         denseMatrix,
                            const std::vector<int32_t>& aoIndices,
                            const int32_t               aoCount,
                            const int32_t               nAOs)
{
    auto nrows = denseMatrix.getNumberOfRows();

    if ((aoCount <= nAOs) && (nrows > 0))
    {
        CDenseMatrix sub_matrix(nrows, aoCount);

        auto dense = denseMatrix.values();

        for (int32_t i = 0; i < nrows; i++)
        {
            auto sub_matrix_row = sub_matrix.row(i);

            auto dense_row = dense + i * nAOs;

            for (int32_t j = 0; j < aoCount; j++)
            {
                sub_matrix_row[j] = dense_row[aoIndices[j]];
            }
        }

        return sub_matrix;
    }

    return CDenseMatrix();
}

CAODensityMatrix
getSubDensityMatrix(const CAODensityMatrix&     densityMatrix,
                    const std::vector<int32_t>& aoIndices,
                    const int32_t               aoCount)
{
    auto naos = densityMatrix.getNumberOfRows(0);

    if (aoCount > naos) return CAODensityMatrix();

    std::vector<CDenseMatrix> submatrices;

    auto numdens = densityMatrix.getNumberOfDensityMatrices();

    for (int32_t idens = 0; idens < numdens; idens++)
    {
        CDenseMatrix sub_dens(aoCount, aoCount);

        auto dens = densityMatrix.getReferenceToDensity(idens);

        for (int32_t i = 0; i < aoCount; i++)
        {
            auto sub_dens_row = sub_dens.row(i);

            auto dens_row = dens.row(aoIndices[i]);

            for (int32_t j = 0; j < aoCount; j++)
            {
                sub_dens_row[j] = dens_row[aoIndices[j]];
            }
        }

        submatrices.push_back(sub_dens);
    }

    return CAODensityMatrix(submatrices, denmat::rest);
}

void
distribute4DSubTo4DFull(CDense4DTensor&          fullMatrix,
                                          const CDense4DTensor&       subMatrix,
                                          const std::vector<int32_t>& aoIndices,
                                          const int32_t               aoCount)
{
    auto W_full = fullMatrix.values();

    auto W_small = subMatrix.values();

    auto nAct =  subMatrix.getjIndex();

    int32_t nAct3 = nAct * nAct * nAct;

    for (int32_t i = 0; i < aoCount; i++) // looping over AO
    {
        auto irow = nAct3 * i;

        auto irow_full = nAct3 * aoIndices[i];

        for (int32_t jkl = 0; jkl < nAct3; jkl++)
        {

            W_full[irow_full + jkl] += W_small[irow + jkl];

        }
    }

}



void
distributeSubMatrixToKohnSham(CAOKohnShamMatrix&          aoKohnShamMatrix,
                              const CDenseMatrix&         subMatrix,
                              const std::vector<int32_t>& aoIndices,
                              const int32_t               aoCount,
                              const int32_t               nAOs)
{
    if (aoCount <= nAOs)
    {
        for (int32_t row = 0; row < subMatrix.getNumberOfRows(); row++)
        {
            auto row_orig = aoIndices[row];

            auto ksmat_row_orig = aoKohnShamMatrix.getMatrix(0) + row_orig * nAOs;

            auto submat_row = subMatrix.row(row);

            for (int32_t col = 0; col < subMatrix.getNumberOfColumns(); col++)
            {
                auto col_orig = aoIndices[col];

                ksmat_row_orig[col_orig] += submat_row[col];
            }
        }
    }
}

void
distributeSubMatrixToKohnSham(CAOKohnShamMatrix&               aoKohnShamMatrix,
                              const std::vector<CDenseMatrix>& subMatrices,
                              const std::vector<int32_t>&      aoIndices,
                              const int32_t                    aoCount,
                              const int32_t                    nAOs)
{
    if (aoCount <= nAOs)
    {
        auto nrows = subMatrices[0].getNumberOfRows();

        auto ncols = subMatrices[0].getNumberOfColumns();

        for (int32_t row = 0; row < nrows; row++)
        {
            auto row_orig = aoIndices[row];

            auto ksmat_a_row_orig = aoKohnShamMatrix.getKohnSham(false) + row_orig * nAOs;

            auto ksmat_b_row_orig = aoKohnShamMatrix.getKohnSham(true) + row_orig * nAOs;

            auto submat_a_row = subMatrices[0].row(row);

            auto submat_b_row = subMatrices[1].row(row);

            for (int32_t col = 0; col < ncols; col++)
            {
                auto col_orig = aoIndices[col];

                ksmat_a_row_orig[col_orig] += submat_a_row[col];

                ksmat_b_row_orig[col_orig] += submat_b_row[col];
            }
        }
    }
}

void
distributeSubMatrixToFock(CAOFockMatrix&              aoFockMatrix,
                          const int32_t               fockIndex,
                          const CDenseMatrix&         subMatrix,
                          const std::vector<int32_t>& aoIndices,
                          const int32_t               aoCount,
                          const int32_t               nAOs)
{
    if (aoCount <= nAOs)
    {
        for (int32_t row = 0; row < subMatrix.getNumberOfRows(); row++)
        {
            auto row_orig = aoIndices[row];

            auto fock_row_orig = aoFockMatrix.getFock(fockIndex, "ALPHA") + row_orig * nAOs;

            auto submat_row = subMatrix.row(row);

            for (int32_t col = 0; col < subMatrix.getNumberOfColumns(); col++)
            {
                auto col_orig = aoIndices[col];

                fock_row_orig[col_orig] += submat_row[col];
            }
        }
    }
}

void
distributeSubMatrixToDenseMatrix(CDenseMatrix&               matrix,
                                 const CDenseMatrix&         subMatrix,
                                 const std::vector<int32_t>& aoIndices,
                                 const int32_t               aoCount,
                                 const int32_t               nAOs)
{
    if (aoCount <= nAOs)
    {
        for (int32_t row = 0; row < subMatrix.getNumberOfRows(); row++)
        {
            auto row_orig = aoIndices[row];

            auto mat_row_orig = matrix.row(row_orig);

            auto submat_row = subMatrix.row(row);

            for (int32_t col = 0; col < subMatrix.getNumberOfColumns(); col++)
            {
                auto col_orig = aoIndices[col];

                mat_row_orig[col_orig] += submat_row[col];
            }
        }
    }
}

}  // namespace submat
