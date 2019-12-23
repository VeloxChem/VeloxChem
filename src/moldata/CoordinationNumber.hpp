//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef CoordinationNumber_hpp
#define CoordinationNumber_hpp

#include "DenseMatrix.hpp"
#include "Molecule.hpp"

#include <vector>

namespace coordnum {  // coordnum namespace

/**
 Creates covalent radii (reference: dftd4 v2.4.0).

 @return a vector of covalent radii with nuclear charge as index.
 */
std::vector<double> getCovalentRadius();

/**
 Creates atomic coordination numbers (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @return a vector of atomic coordination numbers for a molecule.
 */
std::vector<double> getCoordinationNumber(const CMolecule& molecule);

/**
 Creates atomic coordination numbers (reference: dftd4 v2.4.0).

 @param molecule the molecule.
 @param dcndr the derivative matrix of dimension (3N,N).
 @return a vector of atomic coordination numbers for a molecule.
 */
std::vector<double> getCoordinationNumber(const CMolecule& molecule, CDenseMatrix& dcndr);

}  // namespace coordnum

#endif /* CoordinationNumber_hpp */
