//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DistMaxDensity.hpp"

#include <cmath>

#include "AngularMomentum.hpp"

namespace distmaxden { // distmaxden namespace

    void
    getMaxRestDenJK(      CMemBlock<double>&   maxDensityElements,
                    const double*              densityMatrix,
                    const int32_t              nDensityColumns,
                    const CGtoPairsBlock&      braGtoPairsBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
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
        
        auto prefp = braGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefq = braGtoPairsBlock.getKetIdentifiers(0);
        
        //bool symbra = (refp == refq);
        
        // set up pointers to reference indexes on ket side
        
        auto prefk = ketGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefl = ketGtoPairsBlock.getKetIdentifiers(0);
        
        // set up pointer to vector max. density elements
        
        auto mdenvec = maxDensityElements.data();
        
        // set up pointers to contraction pattern on bra side
        
        auto cbspos = braGtoPairsBlock.getContrStartPositions();
        
        auto cbepos = braGtoPairsBlock.getContrEndPositions();
        
        // set up pointers to contraction pattern on ket side
        
        auto ckspos = ketGtoPairsBlock.getContrStartPositions();
        
        auto ckepos = ketGtoPairsBlock.getContrEndPositions();
        
        for (int32_t i = cbspos[iContrPair]; i < cbepos[iContrPair]; i++)
        {
            auto refp = prefp[i];
            
            auto refq = prefq[i];
            
            bool symbra = (refp == refq);
            
            // loop over angular components on bra side
            
            for (int32_t j = 0; j < acomp; j++)
            {
                // set up index P for bra side
                
                int32_t idp  = (braGtoPairsBlock.getBraIdentifiers(j))[i];
                
                // starting angular index of Q
                
                auto kstart = (symbra) ? j : 0;
                
                for (int32_t k = kstart; k < bcomp; k++)
                {
                    // set up index Q for bra side
                    
                    int32_t idq = (braGtoPairsBlock.getKetIdentifiers(k))[i];
                    
                    // loop over angular components on ket side
                    
                    for (int32_t l = 0; l < ccomp; l++)
                    {
                        // set up pointer to R indexes on ket side
                        
                        auto idxk = ketGtoPairsBlock.getBraIdentifiers(l);
                        
                        for (int32_t m = 0; m < dcomp; m++)
                        {
                            // set up pointer to S indexes on ket side
                            
                            auto idxl = ketGtoPairsBlock.getKetIdentifiers(m);
                            
                            // loop over pairs on ket side
                            
                            for (int32_t n = 0; n < nKetContrPairs; n++)
                            {
                                double mden = 0.0;
                                
                                for (int32_t o = ckspos[n]; o < ckepos[n]; o++)
                                {
                                    // symmetry restriction for ket angular components
                                    
                                    auto refr = prefk[o];
                                    
                                    auto refs = prefl[o];
                                    
                                    if ((refr == refs) && (m < l)) continue;
                                    
                                    // symmetry restriction for bra/ket angular componets
                                    
                                    bool braeqket = (refp == refr) && (refq == refs);
                                    
                                    if  (((l * dcomp + m) < (j * bcomp + k)) && braeqket) continue;
                                    
                                    // set up S and R indexes
                                    
                                    auto idr = idxk[o];
                                    
                                    auto ids = idxl[o];
                                    
                                    // Coulomb contributions
                                    
                                    auto cden = 4.0 * std::fabs(densityMatrix[idr * nDensityColumns + ids]);
                                    
                                    if (cden > mden) mden = cden;
                                    
                                    cden = 4.0 * std::fabs(densityMatrix[idp * nDensityColumns + idq]);
                                    
                                    if (cden > mden) mden = cden;
                                    
                                    // exchange contributions
                                    
                                    cden = std::fabs(densityMatrix[idq * nDensityColumns + ids]);
                                    
                                    if (cden > mden) mden = cden;
                                    
                                    cden = std::fabs(densityMatrix[idq * nDensityColumns + idr]);
                                    
                                    if (cden > mden) mden = cden;
                                    
                                    cden = std::fabs(densityMatrix[idp * nDensityColumns + ids]);
                                    
                                    if (cden > mden) mden = cden;
                                    
                                    cden = std::fabs(densityMatrix[idp * nDensityColumns + idr]);
                                    
                                    if (cden > mden) mden = cden;
                                }
                                
                                // copy max. density element
                                
                                if (mden > mdenvec[n]) mdenvec[n] = mden;
                            }
                        }
                    }
                }
            }
        }
    }
    
} // distmaxden namespace
