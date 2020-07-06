//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef AOIndices_hpp
#define AOIndices_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace aoindices {  // aoindices namespace

/**
 Gets AO indices of the two molecules in a molecular dimer.

 @param mol_1 the first molecule.
 @param mol_2 the second molecule.
 @param basis_1 the basis set for the first molecule.
 @param basis_2 the basis set for the first molecule.
 @return a vector of vector containing the AO indices of the two molecules.

 @author X. Li
 */
std::vector<std::vector<int32_t>> getDimerAOIndices(const CMolecule&       mol_1,
                                                    const CMolecule&       mol_2,
                                                    const CMolecularBasis& basis_1,
                                                    const CMolecularBasis& basis_2);

}  // namespace aoindices

#endif /* AOIndices_hpp */
