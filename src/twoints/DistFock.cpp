//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DistFock.hpp"

namespace distfock { // distfock namespace

    void
    distRestJK(      double*              fockMatrix,
               const int32_t              nFockColumns,
               const double*              densityMatrix,
               const int32_t              nDensityColumns,
               const CMemBlock2D<double>& spherInts,
               const CGtoPairsBlock&      braGtoPairsBlock,
               const CGtoPairsBlock&      ketGtoPairsBlock,
               const bool                 isBraEqualKet,
               const int32_t              nKetContrPairs,
               const int32_t              iContrPair)
    {
        
    }
    
    
} // distfock namespace
