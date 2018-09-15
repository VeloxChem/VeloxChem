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

    int32_t max_angl_1 = basis_1.getMolecularMaxAngularMomentum(mol_1);

    int32_t max_angl_2 = basis_2.getMolecularMaxAngularMomentum(mol_2);

    int32_t max_angl = std::max(max_angl_1, max_angl_2);

    std::vector< std::vector<int32_t> > aoinds_mol;

    aoinds_mol.push_back(std::vector<int32_t>());

    aoinds_mol.push_back(std::vector<int32_t>());

    for (int32_t aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        int32_t nao_1 = basis_1.getNumberOfBasisFunctions(mol_1, angl);

        int32_t nao_2 = basis_2.getNumberOfBasisFunctions(mol_2, angl);

        for (int32_t s = -angl; s <= angl; s++)
        {
            for (int32_t i = 0; i < nao_1; i++, aoidx++)
            {
                aoinds_mol[0].push_back(aoidx);
            }

            for (int32_t i = 0; i < nao_2; i++, aoidx++)
            {
                aoinds_mol[1].push_back(aoidx);
            }
        }
    }

    // find out the AO index mapping from monomer to dimer

    const std::vector<int32_t>& aoinds_1 = aoinds_mol[0];

    const std::vector<int32_t>& aoinds_2 = aoinds_mol[1];

    // form the four blocks of dimer matrix

    const int32_t nao_1 = S11.getNumberOfRows();

    const int32_t nao_2 = S22.getNumberOfRows();

    const int32_t nao = nao_1 + nao_2;

    CDenseMatrix smat (nao, nao);

    smat.zero();

    // [1,1] block

    for (int32_t i = 0; i < nao_1; i++)
    {
        for (int32_t j = 0; j < nao_1; j++)
        {
            smat.values()[aoinds_1[i] * nao + aoinds_1[j]] = 
                
                S11.values()[i * nao_1 + j];
        }
    }

    // [2,2] block

    for (int32_t i = 0; i < nao_2; i++)
    {
        for (int32_t j = 0; j < nao_2; j++)
        {
            smat.values()[aoinds_2[i] * nao + aoinds_2[j]] = 
                
                S22.values()[i * nao_2 + j];
        }
    }

    // [1,2] block

    for (int32_t i = 0; i < nao_1; i++)
    {
        for (int32_t j = 0; j < nao_2; j++)
        {
            smat.values()[aoinds_1[i] * nao + aoinds_2[j]] = 
                
                S12.values()[i * nao_2 + j];
        }
    }

    // [2,1] block

    for (int32_t i = 0; i < nao_2; i++)
    {
        for (int32_t j = 0; j < nao_1; j++)
        {
            smat.values()[aoinds_2[i] * nao + aoinds_1[j]] = 
                
                S21.values()[i * nao_1 + j];
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
