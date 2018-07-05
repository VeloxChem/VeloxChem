//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef GenFunc_hpp
#define GenFunc_hpp

#include <cstdint>

#include "MemBlock2D.hpp"
#include "SphericalMomentum.hpp"

namespace genfunc { // genfunc namespace
    
    /**
     Contracts primitive data vectors to contracted data vectors.

     @param contrData the contracted data vectors.
     @param contrIndex the index of first contracted data vector.
     @param primData the primitive data vectors.
     @param primIndex the index of first contracted data vector.
     @param startPositions the vector of start positions in contraction pattern.
     @param endPositions the vector of end positions in contractrion pattern.
     @param nElements the number of elements in individual contracted data
            vector.
     @param nBlocks the number of contracted vectors.
     */
    void contract(      CMemBlock2D<double>& contrData,
                  const int32_t              contrIndex,
                  const CMemBlock2D<double>& primData,
                  const int32_t              primIndex,
                  const int32_t*             startPositions,
                  const int32_t*             endPositions,
                  const int32_t              nElements,
                  const int32_t              nBlocks);
    
    /**
     Transforms Cartesian data vectors to spherical data vectors.

     @param spherData the spherical data vectors.
     @param cartData the Cartesian data vectors.
     @param spherMomentum the spherical momentum objecct.
     @param nElements the number elements in individual data vector.
     @param nBlocks the number of data vectors per spherical momentum component.
     */
    void transform(     CMemBlock2D<double>& spherData,
                  const CMemBlock2D<double>& cartData,
                  const CSphericalMomentum&  spherMomentum,
                  const int32_t              nElements,
                  const int32_t              nBlocks);
    

} // genfunc namespace

#endif /* GenFunc_hpp */
