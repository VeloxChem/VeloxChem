//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef AssembleMatrices_hpp
#define AssembleMatrices_hpp

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "OverlapMatrix.hpp"
#include "KineticEnergyMatrix.hpp"
#include "NuclearPotentialMatrix.hpp"

namespace dimerfunc { // dimerfunc namespace

/**
 Gets the maximum angular momentum of a molecule & basis set.
 
 @param mol the molecule.
 @param basis the basis set for the molecule.
 */
int32_t
getMolecularMaxAngularMomentum(const CMolecule&       mol,
                               const CMolecularBasis& basis);

/**
 Assembles dense matrices of two molecules into one dense matrix of the
 molecular dimer.
 
 @param mol_1 the first molecule.
 @param mol_2 the second molecule.
 @param basis_1 the basis set for the first molecule.
 @param basis_2 the basis set for the first molecule.
 @param S11 the "top left" block of dense matrix (from mol_1).
 @param S22 the "bottom right" block of dense matrix (from mol_2).
 @param S12 the "top right" block of dense matrix (crossing terms).
 @param S21 the "bottom left" block of dense matrix (crossing terms).
 */
CDenseMatrix
assembleDenseMatrices(const CMolecule&       mol_1,
                      const CMolecule&       mol_2,
                      const CMolecularBasis& basis_1,
                      const CMolecularBasis& basis_2,
                      const CDenseMatrix&    S11,
                      const CDenseMatrix&    S22,
                      const CDenseMatrix&    S12,
                      const CDenseMatrix&    S21);

/**
 Assembles overlap matrices of two molecules into one overlap matrix of the
 molecular dimer.
 
 @param mol_1 the first molecule.
 @param mol_2 the second molecule.
 @param basis_1 the basis set for the first molecule.
 @param basis_2 the basis set for the first molecule.
 @param S11 the "top left" block of overlap matrix (from mol_1).
 @param S22 the "bottom right" block of overlap matrix (from mol_2).
 @param S12 the "top right" block of overlap matrix (crossing terms).
 @param S21 the "bottom left" block of overlap matrix (crossing terms).
 */
COverlapMatrix
assembleOverlapMatrices(const CMolecule&       mol_1,
                        const CMolecule&       mol_2,
                        const CMolecularBasis& basis_1,
                        const CMolecularBasis& basis_2,
                        const COverlapMatrix&  S11,
                        const COverlapMatrix&  S22,
                        const COverlapMatrix&  S12,
                        const COverlapMatrix&  S21);

/**
 Assembles kinetic energy matrices of two molecules into one kinetic energy
 matrix of the molecular dimer.
 
 @param mol_1 the first molecule.
 @param mol_2 the second molecule.
 @param basis_1 the basis set for the first molecule.
 @param basis_2 the basis set for the first molecule.
 @param S11 the "top left" block of kinetic energy matrix (from mol_1).
 @param S22 the "bottom right" block of kinetic energy matrix (from mol_2).
 @param S12 the "top right" block of kinetic energy matrix (crossing terms).
 @param S21 the "bottom left" block of kinetic energy matrix (crossing terms).
 */
CKineticEnergyMatrix
assembleKineticEnergyMatrices(const CMolecule&             mol_1,
                              const CMolecule&             mol_2,
                              const CMolecularBasis&       basis_1,
                              const CMolecularBasis&       basis_2,
                              const CKineticEnergyMatrix&  S11,
                              const CKineticEnergyMatrix&  S22,
                              const CKineticEnergyMatrix&  S12,
                              const CKineticEnergyMatrix&  S21);

/**
 Assembles nuclear potential matrices of two molecules into one nuclear
 potential matrix of the molecular dimer.
 
 @param mol_1 the first molecule.
 @param mol_2 the second molecule.
 @param basis_1 the basis set for the first molecule.
 @param basis_2 the basis set for the first molecule.
 @param S11 the "top left" block of nuclear potential matrix (from mol_1).
 @param S22 the "bottom right" block of nuclear potential matrix (from mol_2).
 @param S12 the "top right" block of nuclear potential matrix (crossing terms).
 @param S21 the "bottom left" block of nuclear potential matrix (crossing terms).
 */
CNuclearPotentialMatrix
assembleNuclearPotentialMatrices(const CMolecule&                mol_1,
                                 const CMolecule&                mol_2,
                                 const CMolecularBasis&          basis_1,
                                 const CMolecularBasis&          basis_2,
                                 const CNuclearPotentialMatrix&  S11,
                                 const CNuclearPotentialMatrix&  S22,
                                 const CNuclearPotentialMatrix&  S12,
                                 const CNuclearPotentialMatrix&  S21);

} // dimerfunc namespace

#endif /* AssembleMatrices_hpp */
