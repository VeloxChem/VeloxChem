//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef CubeGenerator_hpp
#define CubeGenerator_hpp

#include <vector>

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularOrbitals.hpp"

namespace cubes { // cubes namespace

/**
 Builds the components of Cartesian angular momentum for a shell.
 
 @param angl the angular momentum of the shell.
 @return a vector of vector, e.g. {{1,0,0}, {0,1,0}, {0,0,1}} for p-shell.
 */
const std::vector<std::vector<int32_t>>
buildCartesianAngularMomentum(int angl);

const std::vector<double>
getPhiAtomicOrbitals(const CMolecule&       molecule,
                     const CMolecularBasis& basis,
                     const double           xp,
                     const double           yp,
                     const double           zp);

/**
 Calculates the psi value of a molecular orbital at a give grid point.
 
 @param molecule the molecule.
 @param basis the basis set for the molecule.
 @param molorb the molecular orbitals of the molecule.
 @param moidx the index of the molecular orbital (0-based).
 @param xp the X coordinate of the grid point.
 @param yp the Y coordinate of the grid point.
 @param zp the Z coordinate of the grid point.
 @return psi value of the molecular orbital.
 */
double
getPsiMolecularOrbital(const CMolecule&          molecule,
                       const CMolecularBasis&    basis,
                       const CMolecularOrbitals& molorb,
                       const int32_t             moidx,
                       const std::string&        mospin,
                       const double              xp,
                       const double              yp,
                       const double              zp);

double
getPsiDensity(const CMolecule&        molecule,
              const CMolecularBasis&  basis,
              const CAODensityMatrix& density,
              const int32_t           densityIndex,
              const std::string&      densitySpin,
              const double            xp,
              const double            yp,
              const double            zp);

} // cubes namespace

#endif /* CubeGenerator_hpp */
