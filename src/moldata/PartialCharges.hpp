//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef PartialCharges_hpp
#define PartialCharges_hpp

#include "DenseMatrix.hpp"
#include "Molecule.hpp"

namespace parchg {  // parchg namespace

/**
 Creates atomic partial charges (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @param netcharge net charge of the molecule.
 @return a vector of atomic partial charges for a molecule.
 */
std::vector<double> getPartialCharges(const CMolecule& molecule, const double netcharge);

/**
 Creates atomic partial charges (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @param netcharge net charge of the molecule.
 @param dqdr the derivative matrix of dimension (3N,N+1).
 @return a vector of atomic partial charges for a molecule.
 */
std::vector<double> getPartialCharges(const CMolecule& molecule, const double netcharge, CDenseMatrix& dqdr);

}  // namespace parchg

#endif /* PartialCharges_hpp */
