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
    std::vector< std::vector<int32_t> > aoIndicesOfAtoms;

    int32_t numAtoms = molecule.getNumberOfAtoms();

    for (int32_t atomIdx = 0; atomIdx < numAtoms; atomIdx++)
    {
        aoIndicesOfAtoms.push_back(std::vector<int32_t>());
    }

    int32_t maxAngMom = basis.getMolecularMaxAngularMomentum(molecule);

    for (int32_t aoIdx = 0, angMom = 0; angMom <= maxAngMom; angMom++)
    {
        for (int32_t s = -angMom; s <= angMom; s++)
        {
            for (int32_t atomIdx = 0; atomIdx < numAtoms; atomIdx++)
            {
                int32_t idElem = molecule.getIdsElemental()[atomIdx];

                int32_t numAOs = basis.getNumberOfBasisFunctions(idElem, angMom);

                for (int32_t i = 0; i < numAOs; i++, aoIdx++)
                {
                    aoIndicesOfAtoms[atomIdx].push_back(aoIdx);
                }
            }
        }
    }

    return aoIndicesOfAtoms;
}

CDenseMatrix
getSADInitialGuess(const CMolecule&       molecule,
                   const CMolecularBasis& basis_1,
                   const CMolecularBasis& basis_2,
                   const COverlapMatrix&  S12,
                   const COverlapMatrix&  S22)
{
    int32_t numAtoms = molecule.getNumberOfAtoms();

    int32_t numAO_1 = S12.getNumberOfRows();

    int32_t numAO_2 = S12.getNumberOfColumns();

    // AO indices for atoms

    std::vector< std::vector<int32_t> >

        aoIndicesOfAtoms_1 = getAOIndicesOfAtoms(molecule, basis_1);

    std::vector< std::vector<int32_t> >
        
        aoIndicesOfAtoms_2 = getAOIndicesOfAtoms(molecule, basis_2);

    // occupation numbers

    std::vector< std::vector<double> > qocc = buildQocc();

    // C_SAD matrix

    CDenseMatrix C_SAD (numAO_2, numAO_1);

    C_SAD.zero();

    for (int atomIdx = 0; atomIdx < numAtoms; atomIdx++) {

        // AO indices for this atom

        const std::vector<int32_t>& aoIdx_1 = aoIndicesOfAtoms_1[atomIdx];

        const std::vector<int32_t>& aoIdx_2 = aoIndicesOfAtoms_2[atomIdx];

        // atomic block of AOs

        CDenseMatrix block_12 (aoIdx_1.size(), aoIdx_2.size());

        CDenseMatrix block_22 (aoIdx_2.size(), aoIdx_2.size());

        block_12.zero();

        block_22.zero();

        for (int32_t i = 0; i < aoIdx_1.size(); i++)
        {
            for (int32_t j = 0; j < aoIdx_2.size(); j++)
            {
                block_12.values()[i * aoIdx_2.size() + j] = 
                    
                    S12.values()[aoIdx_1[i] * numAO_2 + aoIdx_2[j]];
            }
        }

        for (int32_t i = 0; i < aoIdx_2.size(); i++)
        {
            for (int32_t j = 0; j < aoIdx_2.size(); j++)
            {
                block_22.values()[i * aoIdx_2.size() + j] = 
                    
                    S22.values()[aoIdx_2[i] * numAO_2 + aoIdx_2[j]];
            }
        }

        // A = S12' C1(identity)

        CDenseMatrix C1 (aoIdx_1.size(), aoIdx_1.size());

        C1.zero();

        for (int32_t i = 0; i < aoIdx_1.size(); i++)
        {
            C1.values()[i * aoIdx_1.size() + i] = 1.0;
        }

        CDenseMatrix A = denblas::multAtB(block_12, C1);

        // S22^-1

        CDenseDiagonalizer diagdrv;

        diagdrv.diagonalize(block_22);

        if (! diagdrv.getState())
        {
            throw "DenseMatrix diagonalization failed in getSADInitialGuess!";
        }

        CDenseMatrix block_22_inv = diagdrv.getInvertedMatrix();

        // M = A' S22^-1 A

        CDenseMatrix prod = denblas::multAB(block_22_inv, A);

        CDenseMatrix M = denblas::multAtB(A, prod);

        // M^-1/2

        diagdrv.diagonalize(M);

        if (! diagdrv.getState())
        {
            throw "DenseMatrix diagonalization failed in getSADInitialGuess!";
        }

        CDenseMatrix M_invsqrt = diagdrv.getInvertedSqrtMatrix();

        // C2 = S22^-1 A M^-1/2

        prod = denblas::multAB(A, M_invsqrt);

        CDenseMatrix C2 = denblas::multAB(block_22_inv, prod);

        // update C_SAD

        const int32_t idElem = molecule.getIdsElemental()[atomIdx];

        for (int32_t j = 0; j < aoIdx_2.size(); j++)
        {
            for (int32_t i = 0; i < aoIdx_1.size(); i++)
            {
                C_SAD.values()[aoIdx_2[j] * numAO_1 + aoIdx_1[i]] = 

                    C2.values()[j * aoIdx_1.size() + i] * sqrt(qocc[idElem][i]);
            }
        }
    }

    // density matrix

    CDenseMatrix D_SAD = denblas::multABt(C_SAD, C_SAD);

    return D_SAD;
}

} // sad_guess namespace
