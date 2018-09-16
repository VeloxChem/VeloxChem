//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SADGuessTest.hpp"

#include "SADGuessDriver.hpp"
#include "DenseMatrix.hpp"
#include "OverlapMatrix.hpp"
#include "OverlapIntegralsDriver.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"

TEST_F(CSADGuessTest, AtomIdxForAO)
{
    CSADGuessDriver saddrv(mpi::rank(MPI_COMM_WORLD),
                           mpi::nodes(MPI_COMM_WORLD),
                           MPI_COMM_WORLD);

    auto h2o = vlxmol::getMoleculeH2O();
    
    auto min_basis = vlxbas::getMinimalBasisForH2O();

    std::vector< std::vector<int32_t> > ao_inds = 
        
        saddrv.getAOIndicesOfAtoms(h2o, min_basis);

    std::vector< std::vector<int32_t> > ref_inds;

    ref_inds.push_back(std::vector<int32_t>({0,1,4,5,6}));

    ref_inds.push_back(std::vector<int32_t>({2}));

    ref_inds.push_back(std::vector<int32_t>({3}));

    ASSERT_EQ(ao_inds, ref_inds);
}

TEST_F(CSADGuessTest, InitialGuess)
{
    COverlapIntegralsDriver ovldrv(mpi::rank(MPI_COMM_WORLD),
                                   mpi::nodes(MPI_COMM_WORLD),
                                   MPI_COMM_WORLD);
    
    CSADGuessDriver saddrv(mpi::rank(MPI_COMM_WORLD),
                           mpi::nodes(MPI_COMM_WORLD),
                           MPI_COMM_WORLD);

    auto h2o = vlxmol::getMoleculeH2O();
    
    auto min_basis = vlxbas::getMinimalBasisForH2O();

    auto ao_basis  = vlxbas::getMolecularBasisForH2O();

    COutputStream ost(std::string("dummy.out"));
    
    COverlapMatrix S12 = ovldrv.compute(h2o, min_basis, ao_basis, ost, MPI_COMM_WORLD);

    COverlapMatrix S22 = ovldrv.compute(h2o, ao_basis, ost, MPI_COMM_WORLD);

    CDenseMatrix dsad = saddrv.compute(h2o, min_basis, ao_basis, S12, S22, ost, MPI_COMM_WORLD);

    std::vector<double> intvals{  1.057352923440807,  0.129815238298234,  0.111536835372263,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.129815238298234,  0.338924489241330,  0.215629473540853,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.111536835372263,  0.215629473540853,  0.138018642725352,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.172758446767233,  0.147269526277504,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.147269526277504,  0.125541261662431,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.172758446767233,
                                  0.147269526277504,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.147269526277504,
                                  0.125541261662431,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.351956162526701,  0.237186325053813,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.237186325053813,  0.159841931417425,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.351956162526700,
                                  0.237186325053813,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.237186325053813,
                                  0.159841931417425,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.351956162526700,  0.237186325053813,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.237186325053813,  0.159841931417425,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000,  0.000000000000000,
                                 -0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000, -0.000000000000000, -0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                 -0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000,
                                  0.000000000000000,  0.000000000000000,  0.000000000000000};

    ASSERT_EQ(dsad.getNumberOfElements(), intvals.size());

    CDenseMatrix m (intvals, dsad.getNumberOfRows(), dsad.getNumberOfColumns());

    ASSERT_EQ(dsad, m);
}
