//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AOIndices.hpp"

namespace aoindices {  // aoindices namespace

std::vector<std::vector<int32_t>>
getDimerAOIndices(const CMolecule& mol_1, const CMolecule& mol_2, const CMolecularBasis& basis_1, const CMolecularBasis& basis_2)
{
    // get indices of atomic orbitals located on each molecule

    // angular momentum based AO ordering
    // S, P-1, P0, P+1, D-2, D-1, D0, D+1, D+2, ...

    // AO_type:  S S S S S S S S P-1 P-1 P0 P0 P+1 P+1 ...
    // Molecule: A A A A B B B B A   B   A  B  A   B   ...

    int32_t max_angl_1 = basis_1.getMolecularMaxAngularMomentum(mol_1);

    int32_t max_angl_2 = basis_2.getMolecularMaxAngularMomentum(mol_2);

    int32_t max_angl = std::max(max_angl_1, max_angl_2);

    std::vector<std::vector<int32_t>> aoinds_mol;

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

    return aoinds_mol;
}

}  // namespace aoindices
