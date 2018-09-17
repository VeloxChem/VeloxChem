//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SADGuessDriver.hpp"

#include "SystemClock.hpp"
#include "StringFormat.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DenseDiagonalizer.hpp"

CSADGuessDriver::CSADGuessDriver(const int32_t  globRank,
                                 const int32_t  globNodes,
                                       MPI_Comm comm)

    : _globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CSADGuessDriver::~CSADGuessDriver()
{
    
}

std::vector< std::vector<double> >
CSADGuessDriver::_buildQocc() const
{
    std::vector< std::vector<double> > qocc;

    // dummy atom

    qocc.push_back(std::vector<double> ());

    // H,He

    //                                    1s

    qocc.push_back(std::vector<double> ({ 0.5 }));

    qocc.push_back(std::vector<double> ({ 1.0 }));

    // Li,Be

    //                                    1s   2s

    qocc.push_back(std::vector<double> ({ 1.0, 0.5 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0 }));

    // B,C,N,O,F,Ne

    //                                    1s   2s     2p-1   2p0    2p+1

    qocc.push_back(std::vector<double> ({ 1.0, 0.375, 0.375, 0.375, 0.375 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.500, 0.500, 0.500, 0.500 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.625, 0.625, 0.625, 0.625 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.750, 0.750, 0.750, 0.750 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.875, 0.875, 0.875, 0.875 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.000, 1.000, 1.000, 1.000 }));

    // Na,Mg

    //                                    1s   2s   3s     2p-1 2p0  2p+1

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 0.500, 1.0, 1.0, 1.0 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.000, 1.0, 1.0, 1.0 }));

    // Al,Si,P,S,Cl,Ar

    //                                    1s   2s   3s     2p-1 3p-1   2p0  3p0    2p+1 3p+1

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 0.375, 1.0, 0.375, 1.0, 0.375, 1.0, 0.375 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 0.500, 1.0, 0.500, 1.0, 0.500, 1.0, 0.500 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 0.625, 1.0, 0.625, 1.0, 0.625, 1.0, 0.625 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 0.750, 1.0, 0.750, 1.0, 0.750, 1.0, 0.750 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 0.875, 1.0, 0.875, 1.0, 0.875, 1.0, 0.875 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.000, 1.0, 1.000, 1.0, 1.000, 1.0, 1.000 }));

    return qocc;
}

std::vector< std::vector<int32_t> >
CSADGuessDriver::getAOIndicesOfAtoms(const CMolecule&       molecule,
                                     const CMolecularBasis& basis) const
{
    std::vector< std::vector<int32_t> > aoinds_atoms;

    int32_t natoms = molecule.getNumberOfAtoms();

    for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
    {
        aoinds_atoms.push_back(std::vector<int32_t>());
    }

    int32_t max_angl = basis.getMolecularMaxAngularMomentum(molecule);

    for (int32_t aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        for (int32_t s = -angl; s <= angl; s++)
        {
            for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
            {
                int32_t idelem = molecule.getIdsElemental()[atomidx];

                int32_t nao = basis.getNumberOfBasisFunctions(idelem, angl);

                for (int32_t i = 0; i < nao; i++, aoidx++)
                {
                    aoinds_atoms[atomidx].push_back(aoidx);
                }
            }
        }
    }

    return aoinds_atoms;
}

CDenseMatrix
CSADGuessDriver::compute(const CMolecule&       molecule,
                         const CMolecularBasis& basis_1,
                         const CMolecularBasis& basis_2,
                         const COverlapMatrix&  S12,
                         const COverlapMatrix&  S22,
                               COutputStream&   oStream,
                               MPI_Comm         comm) const 
{
    CSystemClock timer;
    
    CDenseMatrix dsad;
    
    if (_locRank == mpi::master())
    {
        // generate SAD guess
        
        dsad = _compSADGuess(molecule, basis_1, basis_2, S12, S22);
    }
    
    _printComputationTime(timer, oStream);
    
    return dsad;
}

CDenseMatrix
CSADGuessDriver::_compSADGuess(const CMolecule&       molecule,
                               const CMolecularBasis& basis_1,
                               const CMolecularBasis& basis_2,
                               const COverlapMatrix&  S12,
                               const COverlapMatrix&  S22) const
{
    int32_t natoms = molecule.getNumberOfAtoms();

    int32_t nao_1 = S12.getNumberOfRows();

    int32_t nao_2 = S12.getNumberOfColumns();

    // AO indices for atoms

    std::vector< std::vector<int32_t> >

        aoinds_atoms_1 = getAOIndicesOfAtoms(molecule, basis_1);

    std::vector< std::vector<int32_t> >
        
        aoinds_atoms_2 = getAOIndicesOfAtoms(molecule, basis_2);

    // occupation numbers

    std::vector< std::vector<double> > qocc = _buildQocc();

    // C_SAD matrix

    CDenseMatrix csad (nao_2, nao_1);

    csad.zero();

    #pragma omp parallel for schedule(dynamic)
    for (int atomidx = 0; atomidx < natoms; atomidx++) {

        // AO indices for this atom

        const std::vector<int32_t>& aoinds_1 = aoinds_atoms_1[atomidx];

        const std::vector<int32_t>& aoinds_2 = aoinds_atoms_2[atomidx];
        
        // set up AO indices dimensions
        
        auto naodim_1 = static_cast<int32_t>(aoinds_1.size());
        
        auto naodim_2 = static_cast<int32_t>(aoinds_2.size());
        
        // atomic block of AOs

        CDenseMatrix block_12 (naodim_1, naodim_2);

        CDenseMatrix block_22 (naodim_2, naodim_2);

        block_12.zero();

        block_22.zero();

        for (int32_t i = 0; i < naodim_1; i++)
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                block_12.values()[i * naodim_2 + j] = 
                    
                    S12.values()[aoinds_1[i] * nao_2 + aoinds_2[j]];
            }
        }

        for (int32_t i = 0; i < naodim_2; i++)
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                block_22.values()[i * naodim_2 + j] = 
                    
                    S22.values()[aoinds_2[i] * nao_2 + aoinds_2[j]];
            }
        }

        // A = S12' C1(identity)

        CDenseMatrix mat_c1 (naodim_1, naodim_1);

        mat_c1.zero();

        for (int32_t i = 0; i < naodim_1; i++)
        {
            mat_c1.values()[i * naodim_1 + i] = 1.0;
        }

        CDenseMatrix mat_a = denblas::multAtB(block_12, mat_c1);

        // S22^-1

        CDenseDiagonalizer diagdrv;

        diagdrv.diagonalize(block_22);

        if (! diagdrv.getState())
        {
            throw "DenseMatrix diagonalization failed in getSADInitialGuess!";
        }

        CDenseMatrix block_22_inv = diagdrv.getInvertedMatrix();

        // M = A' S22^-1 A

        CDenseMatrix prod = denblas::multAB(block_22_inv, mat_a);

        CDenseMatrix mat_m = denblas::multAtB(mat_a, prod);

        // M^-1/2

        diagdrv.diagonalize(mat_m);

        if (! diagdrv.getState())
        {
            throw "DenseMatrix diagonalization failed in getSADInitialGuess!";
        }

        CDenseMatrix mat_m_invsqrt = diagdrv.getInvertedSqrtMatrix();

        // C2 = S22^-1 A M^-1/2

        prod = denblas::multAB(mat_a, mat_m_invsqrt);

        CDenseMatrix mat_c2 = denblas::multAB(block_22_inv, prod);

        // update csad

        const int32_t idelem = molecule.getIdsElemental()[atomidx];

        for (int32_t j = 0; j < naodim_2; j++)
        {
            for (int32_t i = 0; i < naodim_1; i++)
            {
                csad.values()[aoinds_2[j] * nao_1 + aoinds_1[i]] = 

                    mat_c2.values()[j * naodim_1 + i] * sqrt(qocc[idelem][i]);
            }
        }
    }

    // D_SAD density matrix

    CDenseMatrix dsad = denblas::multABt(csad, csad);

    return dsad;
}

void
CSADGuessDriver::_printComputationTime(const CSystemClock&  timer,
                                             COutputStream& oStream) const
{
    auto tsec = timer.getElapsedTimeInSeconds();
    
    if (_isLocalMode)
    {
        // FIX ME: we need tags for each driver to be implemented to manage
        //         MPI send/receive cycle.
    }
    
    if (_globRank == mpi::master())
    {
        oStream << fmt::info << "SAD initial guess computed in ";
        
        oStream << fstr::to_string(tsec, 2) << " sec.";
        
        oStream << fmt::end << fmt::blank;
    }
}
