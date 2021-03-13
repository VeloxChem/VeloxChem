//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef GtoTransform_hpp
#define GtoTransform_hpp

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace gtotra {  // gtotra namespace

/**
 Transforms AO matrix from Dalton to VeloxChem format.

 @param matrix the AO matrix in Dalton format.
 @param basis the molecular basis set.
 @param molecule the molecule.
 @return the AO matrix in VeloxChem format.
 */
CDenseMatrix to_veloxchem(const CDenseMatrix& matrix, const CMolecularBasis& basis, const CMolecule& molecule);

/**
 Transforms AO matrix from VeloxChem to Dalton format.

 @param matrix the AO matrix in VeloxChem format.
 @param basis the molecular basis set.
 @param molecule the molecule.
 @return the AO matrix in Dalton format.
 */
CDenseMatrix to_dalton(const CDenseMatrix& matrix, const CMolecularBasis& basis, const CMolecule& molecule);

/**
 Gets basis function indices for an atom.

 @param molecule the molecule.
 @param basis the molecular basis set.
 @param atomIdx the (0-based) index of the atom in the molecule.
 @return the list of basis function indices.
 */
std::vector<std::vector<int32_t>> getBasisFunctionIndicesForAtom(const CMolecule& molecule, const CMolecularBasis& basis, const int32_t atomIdx);

}  // namespace gtotra

#endif /* GtoTransform_hpp */
