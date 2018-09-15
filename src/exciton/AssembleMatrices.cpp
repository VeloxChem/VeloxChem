//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "AssembleMatrices.hpp"

#include <algorithm>
#include <vector>
#include <string>

namespace dimerfunc { // dimerfunc namespace

CDenseMatrix
assembleDenseMatrices(const CMolecule&       mol_1,
                      const CMolecule&       mol_2,
                      const CMolecularBasis& basis_1,
                      const CMolecularBasis& basis_2,
                      const CDenseMatrix&    S11,
                      const CDenseMatrix&    S22,
                      const CDenseMatrix&    S12,
                      const CDenseMatrix&    S21)
{
    // get indices of atomic orbitals located on each molecule
    
    // angular momentum based AO ordering
    // S, P-1, P0, P+1, D-2, D-1, D0, D+1, D+2, ...
    
    // AO_type:  S S S S S S S S P-1 P-1 P0 P0 P+1 P+1 ...
    // Molecule: A A A A B B B B A   B   A  B  A   B   ...

    int32_t maxAngMom_1 = basis_1.getMolecularMaxAngularMomentum(mol_1);

    int32_t maxAngMom_2 = basis_2.getMolecularMaxAngularMomentum(mol_2);

    int32_t maxAngMom = std::max(maxAngMom_1, maxAngMom_2);

    std::vector< std::vector<int32_t> > aoIndicesOfMolecules;

    aoIndicesOfMolecules.push_back(std::vector<int32_t>());

    aoIndicesOfMolecules.push_back(std::vector<int32_t>());

    for (int32_t aoIdx = 0, angMom = 0; angMom <= maxAngMom; angMom++)
    {
        int32_t numAO_1 = basis_1.getNumberOfBasisFunctions(mol_1, angMom);

        int32_t numAO_2 = basis_2.getNumberOfBasisFunctions(mol_2, angMom);

        for (int32_t s = -angMom; s <= angMom; s++)
        {
            for (int32_t i = 0; i < numAO_1; i++, aoIdx++)
            {
                aoIndicesOfMolecules[0].push_back(aoIdx);
            }

            for (int32_t i = 0; i < numAO_2; i++, aoIdx++)
            {
                aoIndicesOfMolecules[1].push_back(aoIdx);
            }
        }
    }

    // find out the AO index mapping from monomer to dimer

    const std::vector<int32_t>& aoIdx_1 = aoIndicesOfMolecules[0];

    const std::vector<int32_t>& aoIdx_2 = aoIndicesOfMolecules[1];

    // form the four blocks of dimer matrix

    const int32_t numAO_1 = S11.getNumberOfRows();

    const int32_t numAO_2 = S22.getNumberOfRows();

    const int32_t numAO = numAO_1 + numAO_2;

    CDenseMatrix smat (numAO, numAO);

    smat.zero();

    // [1,1] block

    for (int32_t i = 0; i < numAO_1; i++)
    {
        for (int32_t j = 0; j < numAO_1; j++)
        {
            smat.values()[aoIdx_1[i] * numAO + aoIdx_1[j]] = 
                
                S11.values()[i * numAO_1 + j];
        }
    }

    // [2,2] block

    for (int32_t i = 0; i < numAO_2; i++)
    {
        for (int32_t j = 0; j < numAO_2; j++)
        {
            smat.values()[aoIdx_2[i] * numAO + aoIdx_2[j]] = 
                
                S22.values()[i * numAO_2 + j];
        }
    }

    // [1,2] block

    for (int32_t i = 0; i < numAO_1; i++)
    {
        for (int32_t j = 0; j < numAO_2; j++)
        {
            smat.values()[aoIdx_1[i] * numAO + aoIdx_2[j]] = 
                
                S12.values()[i * numAO_2 + j];
        }
    }

    // [2,1] block

    for (int32_t i = 0; i < numAO_2; i++)
    {
        for (int32_t j = 0; j < numAO_1; j++)
        {
            smat.values()[aoIdx_2[i] * numAO + aoIdx_1[j]] = 
                
                S21.values()[i * numAO_1 + j];
        }
    }

    return smat;
}

COverlapMatrix
assembleOverlapMatrices(const CMolecule&       mol_1,
                        const CMolecule&       mol_2,
                        const CMolecularBasis& basis_1,
                        const CMolecularBasis& basis_2,
                        const COverlapMatrix&  S11,
                        const COverlapMatrix&  S22,
                        const COverlapMatrix&  S12,
                        const COverlapMatrix&  S21)
{
    std::vector<double> v11 (S11.values(), S11.values() + S11.getNumberOfElements());

    std::vector<double> v22 (S22.values(), S22.values() + S22.getNumberOfElements());

    std::vector<double> v12 (S12.values(), S12.values() + S12.getNumberOfElements());

    std::vector<double> v21 (S21.values(), S21.values() + S21.getNumberOfElements());

    CDenseMatrix m11 (v11, S11.getNumberOfRows(), S11.getNumberOfColumns());

    CDenseMatrix m22 (v22, S22.getNumberOfRows(), S22.getNumberOfColumns());

    CDenseMatrix m12 (v12, S12.getNumberOfRows(), S12.getNumberOfColumns());

    CDenseMatrix m21 (v21, S21.getNumberOfRows(), S21.getNumberOfColumns());

    CDenseMatrix m = assembleDenseMatrices(mol_1, mol_2,
                                           basis_1, basis_2,
                                           m11, m22, m12, m21);

    return COverlapMatrix(m);
}

CKineticEnergyMatrix
assembleKineticEnergyMatrices(const CMolecule&             mol_1,
                              const CMolecule&             mol_2,
                              const CMolecularBasis&       basis_1,
                              const CMolecularBasis&       basis_2,
                              const CKineticEnergyMatrix&  S11,
                              const CKineticEnergyMatrix&  S22,
                              const CKineticEnergyMatrix&  S12,
                              const CKineticEnergyMatrix&  S21)
{
    std::vector<double> v11 (S11.values(), S11.values() + S11.getNumberOfElements());

    std::vector<double> v22 (S22.values(), S22.values() + S22.getNumberOfElements());

    std::vector<double> v12 (S12.values(), S12.values() + S12.getNumberOfElements());

    std::vector<double> v21 (S21.values(), S21.values() + S21.getNumberOfElements());

    CDenseMatrix m11 (v11, S11.getNumberOfRows(), S11.getNumberOfColumns());

    CDenseMatrix m22 (v22, S22.getNumberOfRows(), S22.getNumberOfColumns());

    CDenseMatrix m12 (v12, S12.getNumberOfRows(), S12.getNumberOfColumns());

    CDenseMatrix m21 (v21, S21.getNumberOfRows(), S21.getNumberOfColumns());

    CDenseMatrix m = assembleDenseMatrices(mol_1, mol_2,
                                           basis_1, basis_2,
                                           m11, m22, m12, m21);

    return CKineticEnergyMatrix(m);
}

CNuclearPotentialMatrix
assembleNuclearPotentialMatrices(const CMolecule&               mol_1,
                                 const CMolecule&               mol_2,
                                 const CMolecularBasis&         basis_1,
                                 const CMolecularBasis&         basis_2,
                                 const CNuclearPotentialMatrix& S11,
                                 const CNuclearPotentialMatrix& S22,
                                 const CNuclearPotentialMatrix& S12,
                                 const CNuclearPotentialMatrix& S21)
{
    std::vector<double> v11 (S11.values(), S11.values() + S11.getNumberOfElements());

    std::vector<double> v22 (S22.values(), S22.values() + S22.getNumberOfElements());

    std::vector<double> v12 (S12.values(), S12.values() + S12.getNumberOfElements());

    std::vector<double> v21 (S21.values(), S21.values() + S21.getNumberOfElements());

    CDenseMatrix m11 (v11, S11.getNumberOfRows(), S11.getNumberOfColumns());

    CDenseMatrix m22 (v22, S22.getNumberOfRows(), S22.getNumberOfColumns());

    CDenseMatrix m12 (v12, S12.getNumberOfRows(), S12.getNumberOfColumns());

    CDenseMatrix m21 (v21, S21.getNumberOfRows(), S21.getNumberOfColumns());

    CDenseMatrix m = assembleDenseMatrices(mol_1, mol_2,
                                           basis_1, basis_2,
                                           m11, m22, m12, m21);

    return CNuclearPotentialMatrix(m);
}

} // dimerfunc namespace
