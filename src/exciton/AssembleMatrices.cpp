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

int32_t
getMolecularMaxAngularMomentum(const CMolecule&       mol,
                               const CMolecularBasis& basis)
{
    int32_t maxAM = 0;

    for (int32_t i = 0; i < mol.getNumberOfAtoms(); i++)
    {
        int32_t idElem = mol.getIdsElemental()[i];

        int32_t elemMaxAM = basis.getMaxAngularMomentum(idElem);

        if (elemMaxAM > maxAM)
        {
            maxAM = elemMaxAM;
        }
    }

    return maxAM;
}

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
    // assign the AOs to molecules
    
    // angular momentum based AO ordering
    // S, P-1, P0, P+1, D-2, D-1, D0, D+1, D+2, ...
    
    // AO_type:  S S S S S S S S P-1 P-1 P0 P0 P+1 P+1 ...
    // Molecule: A A A A B B B B A   B   A  B  A   B   ...

    int32_t maxAM_1 = getMolecularMaxAngularMomentum(mol_1, basis_1);

    int32_t maxAM_2 = getMolecularMaxAngularMomentum(mol_2, basis_2);

    std::vector<std::string> molIdx;

    for (int32_t i = 0; i <= std::max(maxAM_1, maxAM_2); i++)
    {
        for (int32_t s = -i; s <= i; s++)
        {
            if (i <= maxAM_1)
            {
                for (int32_t k = 0; k < basis_1.getNumberOfBasisFunctions(mol_1, i); k++)
                {
                    molIdx.push_back("A");
                }
            }

            if (i <= maxAM_2)
            {
                for (int32_t k = 0; k < basis_2.getNumberOfBasisFunctions(mol_2, i); k++)
                {
                    molIdx.push_back("B");
                }
            }
        }
    }

    // find out the AO index mapping from monomer to dimer

    std::vector<int32_t> aoIdx_1;

    std::vector<int32_t> aoIdx_2;

    for (int32_t i = 0; i < molIdx.size(); i++)
    {
        if (molIdx[i] == std::string("A"))
        {
            aoIdx_1.push_back(i);
        }

        if (molIdx[i] == std::string("B"))
        {
            aoIdx_2.push_back(i);
        }
    }

    // form the four blocks of dimer matrix

    const int32_t nAO_1 = S11.getNumberOfRows();

    const int32_t nAO_2 = S22.getNumberOfRows();

    const int32_t nAO = nAO_1 + nAO_2;

    CDenseMatrix smat (nAO, nAO);

    smat.zero();

    // [1,1] block

    for (int32_t i = 0; i < nAO_1; i++)
    {
        for (int32_t j = 0; j < nAO_1; j++)
        {
            smat.values()[aoIdx_1[i] * nAO + aoIdx_1[j]] = S11.values()[i * nAO_1 + j];
        }
    }

    // [2,2] block

    for (int32_t i = 0; i < nAO_2; i++)
    {
        for (int32_t j = 0; j < nAO_2; j++)
        {
            smat.values()[aoIdx_2[i] * nAO + aoIdx_2[j]] = S22.values()[i * nAO_2 + j];
        }
    }

    // [1,2] block

    for (int32_t i = 0; i < nAO_1; i++)
    {
        for (int32_t j = 0; j < nAO_2; j++)
        {
            smat.values()[aoIdx_1[i] * nAO + aoIdx_2[j]] = S12.values()[i * nAO_2 + j];
        }
    }

    // [2,1] block

    for (int32_t i = 0; i < nAO_2; i++)
    {
        for (int32_t j = 0; j < nAO_1; j++)
        {
            smat.values()[aoIdx_2[i] * nAO + aoIdx_1[j]] = S21.values()[i * nAO_1 + j];
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

