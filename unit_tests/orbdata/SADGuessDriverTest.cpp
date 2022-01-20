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

#include "SADGuessDriverTest.hpp"

#include <mpi.h>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "SADGuessDriver.hpp"
#include "CheckFunctions.hpp"

TEST_F(CSADGuessDriverTest, AtomIdxForAO)
{
    CSADGuessDriver saddrv(MPI_COMM_WORLD);

    auto h2o = vlxmol::getMoleculeH2O();

    auto min_basis = vlxbas::getMinimalBasisForH2O();

    auto ao_inds = saddrv.getAOIndicesOfAtoms(h2o, min_basis);

    std::vector<std::vector<int32_t>> ref_inds;

    ref_inds.push_back(std::vector<int32_t>({0, 1, 4, 5, 6}));

    ref_inds.push_back(std::vector<int32_t>({2}));

    ref_inds.push_back(std::vector<int32_t>({3}));

    ASSERT_EQ(ao_inds, ref_inds);
}

TEST_F(CSADGuessDriverTest, OccupationNumbers)
{
    std::vector<std::vector<double>> refocc;

    refocc.push_back(std::vector<double>());

    refocc.push_back(std::vector<double>({0.5}));

    refocc.push_back(std::vector<double>({1}));

    refocc.push_back(std::vector<double>({1, 0.5}));

    refocc.push_back(std::vector<double>({1, 1}));

    refocc.push_back(std::vector<double>({1, 0.375, 0.375, 0.375, 0.375}));

    refocc.push_back(std::vector<double>({1, 0.5, 0.5, 0.5, 0.5}));

    refocc.push_back(std::vector<double>({1, 0.625, 0.625, 0.625, 0.625}));

    refocc.push_back(std::vector<double>({1, 0.75, 0.75, 0.75, 0.75}));

    refocc.push_back(std::vector<double>({1, 0.875, 0.875, 0.875, 0.875}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 0.5, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 0.375, 1, 0.375, 1, 0.375, 1, 0.375}));

    refocc.push_back(std::vector<double>({1, 1, 0.5, 1, 0.5, 1, 0.5, 1, 0.5}));

    refocc.push_back(std::vector<double>({1, 1, 0.625, 1, 0.625, 1, 0.625, 1, 0.625}));

    refocc.push_back(std::vector<double>({1, 1, 0.75, 1, 0.75, 1, 0.75, 1, 0.75}));

    refocc.push_back(std::vector<double>({1, 1, 0.875, 1, 0.875, 1, 0.875, 1, 0.875}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 0.5, 1, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.2, 0.2, 0.2, 0.2, 0.2}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.3, 0.3, 0.3, 0.3, 0.3}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.4, 0.4, 0.4, 0.4, 0.4}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.6, 0.6, 0.6, 0.6, 0.6}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.8, 0.8, 0.8, 0.8, 0.8}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.9, 0.9, 0.9, 0.9, 0.9}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 0.375, 1, 1, 0.375, 1, 1, 0.375, 1, 1, 0.375, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 0.5, 1, 1, 0.5, 1, 1, 0.5, 1, 1, 0.5, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 0.625, 1, 1, 0.625, 1, 1, 0.625, 1, 1, 0.625, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 0.75, 1, 1, 0.75, 1, 1, 0.75, 1, 1, 0.75, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 0.875, 1, 1, 0.875, 1, 1, 0.875, 1, 1, 0.875, 1, 1, 1, 1, 1}));

    refocc.push_back(std::vector<double>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}));

    CSADGuessDriver saddrv(MPI_COMM_WORLD);

    for (int32_t i = 0; i < static_cast<int32_t>(refocc.size()); i++)
    {
        auto qocc = saddrv.getOccupationNumbersForElement(i, 0.0);

        ASSERT_EQ(refocc[i].size(), qocc.size());

        vlxtest::compare(refocc[i], qocc.data());
    }
}

TEST_F(CSADGuessDriverTest, InitialGuess)
{
    CSADGuessDriver saddrv(MPI_COMM_WORLD);

    auto h2o = vlxmol::getMoleculeH2O();

    auto min_basis = vlxbas::getMinimalBasisForH2O();

    auto ao_basis = vlxbas::getMolecularBasisForH2O();

    auto dsad = saddrv.compute(h2o, min_basis, ao_basis, "restricted");

    ASSERT_EQ(1, dsad.getNumberOfDensityMatrices());

    int32_t nrows = dsad.getNumberOfRows(0);

    int32_t ncols = dsad.getNumberOfColumns(0);

    ASSERT_EQ(nrows * ncols, dsad.getNumberOfElements(0));

    ASSERT_EQ(dsad.alphaDensity(0), dsad.betaDensity(0));

    // out of bounds test

    ASSERT_EQ(0, dsad.getNumberOfRows(1));

    ASSERT_EQ(0, dsad.getNumberOfColumns(1));

    ASSERT_EQ(0, dsad.getNumberOfElements(1));

    ASSERT_EQ(nullptr, dsad.alphaDensity(1));

    ASSERT_EQ(nullptr, dsad.betaDensity(1));
}
