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
        
        // determine symmetry of angular components on bra side
        
        auto refp = (braGtoPairsBlock.getBraIdentifiers(0))[iContrPair];
        
        auto refq = (braGtoPairsBlock.getKetIdentifiers(0))[iContrPair];
        
        bool symbra = (refp == refq);
        
        // set up pointers to reference indexes on ket side
        
        auto prefk = ketGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefl = ketGtoPairsBlock.getKetIdentifiers(0);
        
        // loop over angular components on bra side
        
        int32_t bracomp = 0;
        
        for (int32_t i = 0; i < acomp; i++)
        {
            // set up index P for bra side
            
            int32_t idp = (braGtoPairsBlock.getBraIdentifiers(i))[iContrPair];
         
            // starting angular index of Q
            
            auto jstart = (symbra) ? i : 0;
            
            for (int32_t j = jstart; j < bcomp; j++)
            {
                // set up index Q for bra side
                
                int32_t idq = (braGtoPairsBlock.getKetIdentifiers(j))[iContrPair];
                
                // angular offset from bra side
                
                auto braoff = (i * bcomp + j) * ccomp * dcomp;
                
                // loop over angular components on ket side
                
                int32_t ketcomp = 0;
                
                for (int32_t k = 0; k < ccomp; k++)
                {
                    // set up pointer to R indexes on ket side
                    
                    auto idxk = ketGtoPairsBlock.getBraIdentifiers(k);
                    
                    for (int32_t l = 0; l < dcomp; l++)
                    {
                        // set up pointer to S indexes on ket side
                        
                        auto idxl = ketGtoPairsBlock.getKetIdentifiers(l);
                        
                        // set up pointer to integrals
                        
                        auto pints = spherInts.data(braoff + k * dcomp + l);
                        
                        // loop over pairs on ket side
                        
                        for (int32_t m = 0; m < nKetContrPairs; m++)
                        {
                            // symmetry restriction for ket angular components
                            
                            auto refr = prefk[m];
                            
                            auto refs = prefl[m];
                            
                            if ((refr == refs) && (l < k)) continue;
                            
                            // semmetry restriction for bra/ket angular componets
                            
                            bool braeqket = (refp == refr) && (refq == refs);
                            
                            bool redbraket = (braeqket) ? (bracomp <= ketcomp) : true;
                            
                            if (redbraket)
                            {
                                // set up S and R indexes
                            
                                auto idr = idxk[m];
                            
                                auto ids = idxl[m];
                            
                                // scale integral value
                            
                                auto fval = pints[m];
                            
                                if (idp == idq) fval *= 0.5;
                            
                                if (idr == ids) fval *= 0.5;
                            
                                if ((idp == idr) && (idq == ids)) fval *= 0.5;
                            
                                // Coulomb contributions
                            
                                fockMatrix[idp * nFockColumns + idq] += 4.0 * fval * densityMatrix[idr * nDensityColumns + ids];
                            
                                fockMatrix[idr * nFockColumns + ids] += 4.0 * fval * densityMatrix[idp * nDensityColumns + idq];
                            
                                // exchange contributions
                            
                                fockMatrix[idp * nFockColumns + idr] -= fval * densityMatrix[idq * nDensityColumns + ids];
                            
                                fockMatrix[idp * nFockColumns + ids] -= fval * densityMatrix[idq * nDensityColumns + idr];
                            
                                fockMatrix[idq * nFockColumns + idr] -= fval * densityMatrix[idp * nDensityColumns + ids];
                            
                                fockMatrix[idq * nFockColumns + ids] -= fval * densityMatrix[idp * nDensityColumns + idr];
                            }
                        }
                        
                        // update angular components counter for ket
                        
                        ketcomp++;
                    }
                }
                
                // update angular components counter for bra
                
                bracomp++;
            }
        }
    }
    
    
} // distfock namespace
