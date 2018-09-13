//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef SADGuess_hpp
#define SADGuess_hpp

#include <vector>

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "DenseDiagonalizer.hpp"
#include "DenseLinearAlgebra.hpp"
#include "OverlapMatrix.hpp"

namespace sad_guess { // sad_guess namespace

/**
 Gets the occupation numbers for SAD initial guess.

 @return a 2D vector of occupation numbers for SAD initial guess.
 */
std::vector< std::vector<double> > buildQocc();

/**
 Gets the atomic indices associated with each atomic orbital.
 
 @param molecule the molecule.
 @param basis the basis set.
 @return a vector of atomic indices for atomic oribtals.
 */
std::vector<int32_t>
getAtomIdxForAO(const CMolecule&       molecule,
                const CMolecularBasis& basis);

/**
 Form the density matrix of SAD initial guess.
 
 @param molecule the molecule.
 @param basis_1 the minimal (smaller) basis set.
 @param basis_2 the normal (larger) basis set.
 @param S12 the crossing overlap matrix between basis_1 and basis_2.
 @param S22 the overlap matrix from basis_2.
 @return the density matrix of SAD initial guess.
 */
CDenseMatrix
getSADInitialGuess(const CMolecule&       molecule,
                   const CMolecularBasis& basis_1,
                   const CMolecularBasis& basis_2,
                   const COverlapMatrix&  S12,
                   const COverlapMatrix&  S22);

} // sad_guess namespace

#endif /* SADGuess_hpp */
