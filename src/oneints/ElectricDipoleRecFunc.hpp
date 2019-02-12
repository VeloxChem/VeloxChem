//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ElectricDipoleRecFunc_hpp
#define ElectricDipoleRecFunc_hpp

#include "MemBlock2D.hpp"
#include "GtoBlock.hpp"
#include "VecIndexes.hpp"

namespace ediprecfunc { // ediprecfunc namespace
    
    /**
     Computes batch of primitive (S|M|S) electric dipole integrals and stores
     results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectricDipoleForSS(      CMemBlock2D<double>&  primBuffer,
                                 const CVecThreeIndexes&     recPattern,
                                 const std::vector<int32_t>& recIndexes,
                                 const CMemBlock2D<double>&  osFactors,
                                 const CMemBlock2D<double>&  abDistances,
                                 const CMemBlock2D<double>&  pcDistances,
                                 const CGtoBlock&            braGtoBlock,
                                 const CGtoBlock&            ketGtoBlock,
                                 const int32_t               iContrGto);
    
} // ediprecfunc namespace

#endif /* ElectricDipoleRecFunc_hpp */
