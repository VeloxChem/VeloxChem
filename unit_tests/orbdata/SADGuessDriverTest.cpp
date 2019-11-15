//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "SADGuessDriverTest.hpp"

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "OverlapIntegralsDriver.hpp"
#include "OverlapMatrix.hpp"
#include "SADGuessDriver.hpp"

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

TEST_F(CSADGuessDriverTest, InitialGuess)
{
    COverlapIntegralsDriver ovldrv(MPI_COMM_WORLD);

    CSADGuessDriver saddrv(MPI_COMM_WORLD);

    auto h2o = vlxmol::getMoleculeH2O();

    auto min_basis = vlxbas::getMinimalBasisForH2O();

    auto ao_basis = vlxbas::getMolecularBasisForH2O();

    auto S12 = ovldrv.compute(h2o, min_basis, ao_basis);

    auto S22 = ovldrv.compute(h2o, ao_basis);

    auto dsad = saddrv.compute(h2o, min_basis, ao_basis, S12, S22, true);

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
