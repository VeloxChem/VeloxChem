//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef PartitionFunc_hpp
#define PartitionFunc_hpp

#include <cstdint>

namespace partfunc { // partfunc namespace
    
    /**
     Applies SSF partitioning scheme to grid weights.
     Reference: R. E. Stratmann, G. E. Scuseria, and M. J. Frisch, Chem. Phys.
     Lett., 213 (257), 1996.

     @param gridCoordsX the vector of Cartesian X coordinates of grid points.
     @param gridCoordsY the vector of Cartesian Y coordinates of grid points.
     @param gridCoordsZ the vector of Cartesian Z coordinates of grid points.
     @param gridWeights the vector of weights of grid points.
     @param nGridPoints the number of grid points.
     @param molCoordsX the vector of Cartesian X coordinates of atoms in molecule.
     @param molCoordsY the vector of Cartesian Y coordinates of atoms in molecule.
     @param molCoordsZ the vector of Cartesian Z coordinates of atoms in molecule.
     @param nAtoms the number of atoms in molecule.
     @param partWeights the temporary SSF weights.
     @param idAtom the identifier of atom associated with grid points.
     */
    void ssf(double* gridCoordsX, double* gridCoordsY, double* gridCoordsZ,
             double* gridWeights, const int32_t nGridPoints,
             const double* molCoordsX, const double* molCoordsY,
             const double* molCoordsZ, const int32_t nAtoms,
             double* partWeights, const int32_t idAtom);
    
    /**
     Computes polynomial weight function in eliptical coordinates. See Eq. 14
     in R. E. Stratmann, G. E. Scuseria, and M. J. Frisch, Chem. Phys. Lett.,
     213 (257), 1996.

     @param eRadius the eliptical coordinate.
     @return the polynomial weight.
     */
    double zeta(const double eRadius);
    
} // partfunc namespace

#endif /* PartitionFunc_hpp */
