//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ThreeCenterHrrFunc_hpp
#define ThreeCenterHrrFunc_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "GtoPairsBlock.hpp"
#include "VecIndexes.hpp"

namespace t3hrrfunc { // t3hrrfunc namespace
 
    /**
     Computes batch of half transformed (X|g(r,r')|PP) electron repulsion
     integrals and stores results in integrals buffer.
     
     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braAngularMomentum the angular momentum of bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForXPP(      CMemBlock2D<double>&  contrBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  cdDistances,
                                     const int32_t               braAngularMomentum,
                                     const CGtoPairsBlock&       ketGtoPairsBlock);
    
} // t3hrrfunc namespace

#endif /* ThreeCenterHrrFunc_hpp */
