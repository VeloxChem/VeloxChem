//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef PartialCharges_hpp
#define PartialCharges_hpp

#include "Molecule.hpp"

#include <vector>

namespace parchg {  // parchg namespace

/**
 Creates atomic partial charges.

 @param molecule the molecule.
 @param netcharge net charge of the molecule.
 @return a vector of atomic partial charges for a molecule.
 */
std::vector<double> getPartialCharges(const CMolecule& molecule, const double netcharge);

/**
 Creates atomic electronegativity.

 @return a vector of atomic electronegativity with nuclear charge as index.
 */
std::vector<double> getElectronegativity();

/**
 Creates atomic hardness.

 @return a vector of atomic hardness with nuclear charge as index.
 */
std::vector<double> getHardness();

}  // namespace parchg

#endif /* PartialCharges_hpp */
