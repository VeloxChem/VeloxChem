//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SADGuess.hpp"

#include "DenseLinearAlgebra.hpp"
#include "DenseDiagonalizer.hpp"

namespace sad_guess { // sad_guess namespace

std::vector< std::vector<double> > buildQocc()
{
    std::vector< std::vector<double> > qocc;

    // dummy atom

    qocc.push_back(std::vector<double> ());

    // H,He

    qocc.push_back(std::vector<double> ({ 0.5 }));

    qocc.push_back(std::vector<double> ({ 1.0 }));

    // Li,Be,B,C,N,O,F,Ne

    qocc.push_back(std::vector<double> ({ 1.0, 0.5 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.375, 0.375, 0.375, 0.375 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.500, 0.500, 0.500, 0.500 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.625, 0.625, 0.625, 0.625 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.750, 0.750, 0.750, 0.750 }));

    qocc.push_back(std::vector<double> ({ 1.0, 0.875, 0.875, 0.875, 0.875 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.000, 1.000, 1.000, 1.000 }));

    // Na,Mg,Al,Si,P,S,Cl,Ar

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.0, 1.0, 1.0, 0.500 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.0, 1.0, 1.0, 1.000 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.0, 1.0, 1.0, 0.375, 0.375, 0.375, 0.375 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.0, 1.0, 1.0, 0.500, 0.500, 0.500, 0.500 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.0, 1.0, 1.0, 0.625, 0.625, 0.625, 0.625 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.0, 1.0, 1.0, 0.750, 0.750, 0.750, 0.750 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.0, 1.0, 1.0, 0.875, 0.875, 0.875, 0.875 }));

    qocc.push_back(std::vector<double> ({ 1.0, 1.0, 1.0, 1.0, 1.0, 1.000, 1.000, 1.000, 1.000 }));

    return qocc;
}

std::vector< std::vector<int32_t> >
getAOIndicesOfAtoms(const CMolecule&       molecule,
                    const CMolecularBasis& basis)
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
getSADInitialGuess(const CMolecule&       molecule,
                   const CMolecularBasis& basis_1,
                   const CMolecularBasis& basis_2,
                   const COverlapMatrix&  S12,
                   const COverlapMatrix&  S22)
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

    std::vector< std::vector<double> > qocc = buildQocc();

    // C_SAD matrix

    CDenseMatrix csad (nao_2, nao_1);

    csad.zero();

    for (int atomidx = 0; atomidx < natoms; atomidx++) {

        // AO indices for this atom

        const std::vector<int32_t>& aoinds_1 = aoinds_atoms_1[atomidx];

        const std::vector<int32_t>& aoinds_2 = aoinds_atoms_2[atomidx];

        // atomic block of AOs

        CDenseMatrix block_12 (aoinds_1.size(), aoinds_2.size());

        CDenseMatrix block_22 (aoinds_2.size(), aoinds_2.size());

        block_12.zero();

        block_22.zero();

        for (int32_t i = 0; i < aoinds_1.size(); i++)
        {
            for (int32_t j = 0; j < aoinds_2.size(); j++)
            {
                block_12.values()[i * aoinds_2.size() + j] = 
                    
                    S12.values()[aoinds_1[i] * nao_2 + aoinds_2[j]];
            }
        }

        for (int32_t i = 0; i < aoinds_2.size(); i++)
        {
            for (int32_t j = 0; j < aoinds_2.size(); j++)
            {
                block_22.values()[i * aoinds_2.size() + j] = 
                    
                    S22.values()[aoinds_2[i] * nao_2 + aoinds_2[j]];
            }
        }

        // A = S12' C1(identity)

        CDenseMatrix mat_c1 (aoinds_1.size(), aoinds_1.size());

        mat_c1.zero();

        for (int32_t i = 0; i < aoinds_1.size(); i++)
        {
            mat_c1.values()[i * aoinds_1.size() + i] = 1.0;
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

        for (int32_t j = 0; j < aoinds_2.size(); j++)
        {
            for (int32_t i = 0; i < aoinds_1.size(); i++)
            {
                csad.values()[aoinds_2[j] * nao_1 + aoinds_1[i]] = 

                    mat_c2.values()[j * aoinds_1.size() + i] * sqrt(qocc[idelem][i]);
            }
        }
    }

    // D_SAD density matrix

    CDenseMatrix dsad = denblas::multABt(csad, csad);

    return dsad;
}

} // sad_guess namespace
