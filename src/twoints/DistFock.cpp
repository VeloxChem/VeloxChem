//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DistFock.hpp"

#include "AngularMomentum.hpp"

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
        // set up angular momentum on bra
        
        auto aang = braGtoPairsBlock.getBraAngularMomentum();
        
        auto bang = braGtoPairsBlock.getKetAngularMomentum();
        
        // set up angular momentum on ket
        
        auto cang = ketGtoPairsBlock.getBraAngularMomentum();
        
        auto dang = ketGtoPairsBlock.getKetAngularMomentum();
        
        // set up number of angular components on bra
        
        auto acomp = angmom::to_SphericalComponents(aang);
        
        auto bcomp = angmom::to_SphericalComponents(bang);
        
        // set up number of angular components on ket
        
        auto ccomp = angmom::to_SphericalComponents(cang);
        
        auto dcomp = angmom::to_SphericalComponents(dang);
        
        // set up symmetry on bra side
        
        auto idzi = braGtoPairsBlock.getBraIdentifiers(0);
        
        auto idzj = braGtoPairsBlock.getKetIdentifiers(0);
        
        bool symbra = idzi[iContrPair] == idzj[iContrPair];
        
        // set up pointers to zero identifiers on ket side
        
        auto idzk = braGtoPairsBlock.getBraIdentifiers(0);
        
        auto idzl = braGtoPairsBlock.getKetIdentifiers(0);
        
        // loop over angular components on bra side
        
        int32_t bracomp = 0;
        
        for (int32_t i = 0; i < acomp; i++)
        {
            // apply symmetry on bra side
            
            auto joff = symbra ? i : 0;
            
            // set up index
            
            for (int32_t j = joff; j < bcomp; j++)
            {
                // loop over angular components on ket side
                
                int32_t ketcomp = 0;
                
                for (int32_t k = 0; k < ccomp; k++)
                {
                    for (int32_t l = 0; l < dcomp; l++)
                    {
                        // loop over pairs on ket side
                        
                        for (int32_t m = 0; m < nKetContrPairs; m++)
                        {
                            // apply symmetry of ket side
                            
                            if ((idzk[m] == idzl[m]) && (l < k)) continue;
                            
                            // set up indexes on ket side
                            
                            // scale integral values
                        }
                        
                        // update ket componets counter
                        
                        ketcomp++;
                    }
                }
                
                // update bra components counter
                
                bracomp++;
            }
        }
        
    }
    
    
} // distfock namespace
