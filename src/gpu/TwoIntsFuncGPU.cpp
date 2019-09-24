//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TwoIntsFuncGPU.hpp"

namespace twointsgpu { // twointsgpu namespace
    
    void
    compDistancesPQ(        double*         pqDistancesData,
                    const   size_t          pitchOfDistancesData,
                    const   double*         braGtoPairsData,
                    const   size_t          pitchOfBraGtoPairsData,
                    const   double*         ketGtoPairsData,
                    const   size_t          pitchOfKetGtoPairsData,
                    const   CGtoPairsBlock& braGtoPairsBlock,
                    const   int32_t         nKetPrimPairs,
                    const   int32_t         iContrPair,
                    const   CCudaDevices*   cudaDevices)
    {
        
    }
    
} // intsfunc namespace
