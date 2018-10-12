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
        
        // loop over angular components on bra side
        
        for (int32_t i = 0; i < acomp; i++)
        {
            // apply symmetry on bra side
            
            //auto joff = symbra ? i : 0;
            
            // set up index P for bra side
            
            int32_t idp = (braGtoPairsBlock.getBraIdentifiers(i))[iContrPair];
            
            for (int32_t j = 0; j < bcomp; j++)
            {
                // set up index Q for bra side
                
                int32_t idq = (braGtoPairsBlock.getKetIdentifiers(j))[iContrPair];
                
                // angular offset from bra side
                
                auto braoff = (i * bcomp + j) * ccomp * dcomp;
                
                // loop over angular components on ket side
                
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
                            // set up S and R indexes
                            
                            auto idr = idxk[m];
                            
                            auto ids = idxl[m];
                            
                            // scale integral value
                            
                            auto fval = pints[m];
                            
                            if (idp == idq) fval *= 0.5;
                            
                            if (idr == ids) fval *= 0.5;
                            
                            if ((idp == idr) && (idq == ids)) fval *= 0.5;
                            
                            // Coulomb contribution
                            
                            auto fpq = 4.0 * fval * densityMatrix[idr * nDensityColumns + ids];
                            
                            fockMatrix[idp * nFockColumns + idq] += fpq;
                            
                            //fockMatrix[idq * nFockColumns + idp] += fpq;
                            
                            auto frs = 4.0 * fval * densityMatrix[idp * nDensityColumns + idq];
                            
                            fockMatrix[idr * nFockColumns + ids] += frs;
                            
                            //fockMatrix[ids * nFockColumns + idr] += frs;
                            
                            // exchange contribution
                            
                            auto fpr = fval * densityMatrix[idq * nDensityColumns + ids];
                            
                            fockMatrix[idp * nFockColumns + idr] -= fpr;
                            
                            //fockMatrix[idr * nFockColumns + idp] -= fpr;
                            
                            auto fps = fval * densityMatrix[idq * nDensityColumns + idr];
                            
                            fockMatrix[idp * nFockColumns + ids] -= fps;
                            
                            //fockMatrix[ids * nFockColumns + idp] -= fps;
                            
                            auto fqr = fval * densityMatrix[idp * nDensityColumns + ids];
                            
                            fockMatrix[idq * nFockColumns + idr] -= fqr;
                            
                            //fockMatrix[idr * nFockColumns + idq] -= fqr;
                            
                            auto fqs = fval * densityMatrix[idp * nDensityColumns + idr];
                            
                            fockMatrix[idq * nFockColumns + ids] -= fqs;
                            
                            //fockMatrix[ids * nFockColumns + idq] -= fqs;
                        }
                    }
                }
            }
        }
    }
    
    
} // distfock namespace
