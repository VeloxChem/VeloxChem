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
    distRestJK(      CFockContainer&      fockContainer,
               const int32_t              iFockMatrix,
               const double*              densityMatrix,
               const int32_t              nDensityColumns,
               const CMemBlock2D<double>& spherInts,
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
        
        // set up pointers to reference indexes on bra side
        
        auto prefp = braGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefq = braGtoPairsBlock.getKetIdentifiers(0);
        
        // set up pointers to reference indexes on ket side
        
        auto prefk = ketGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefl = ketGtoPairsBlock.getKetIdentifiers(0);
        
        // set up pointers to start positions for A, B, C, D centers
        
        auto apos = fockContainer.getStartPositionsA(iFockMatrix);
        
        auto bpos = fockContainer.getStartPositionsB(iFockMatrix);
        
        auto cpos = fockContainer.getStartPositionsC(iFockMatrix);
        
        auto dpos = fockContainer.getStartPositionsD(iFockMatrix);
        
        // set up dimensions of submatrices for A, B, C, D centers
        
        auto bdim = fockContainer.getDimensionsB(iFockMatrix);
        
        auto cdim = fockContainer.getDimensionsC(iFockMatrix);
        
        auto ddim = fockContainer.getDimensionsD(iFockMatrix);
        
        // set up pointers to contraction pattern on bra side
        
        auto cbspos = braGtoPairsBlock.getContrStartPositions();
        
        auto cbepos = braGtoPairsBlock.getContrEndPositions();
        
        // set up pointers to contraction pattern on ket side
        
        auto ckspos = ketGtoPairsBlock.getContrStartPositions();
        
        auto ckepos = ketGtoPairsBlock.getContrEndPositions();
        
        // dimensions of ket side
        
        auto ketdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs(nKetContrPairs - 1);
        
        // loop over integrals
        
        int32_t iblock = 0;
        
        for (int32_t i = cbspos[iContrPair]; i < cbepos[iContrPair]; i++)
        {
            auto refp = prefp[i];
            
            auto refq = prefq[i];
            
            bool symbra = (refp == refq);
            
            auto intoff = iblock * ketdim;
            
            // loop over angular components on bra side
            
            for (int32_t j = 0; j < acomp; j++)
            {
                // set up index P for bra side
                
                int32_t idp  = (braGtoPairsBlock.getBraIdentifiers(j))[i];
                
                int32_t redp = idp - apos[j];
                
                auto kstart = (symbra) ? j : 0;
                
                for (int32_t k = kstart; k < bcomp; k++)
                {
                    // set up index Q for bra side
                    
                    int32_t idq = (braGtoPairsBlock.getKetIdentifiers(k))[i];
                    
                    int32_t redq = idq - bpos[k];
                    
                    // angular offset from bra side
                    
                    auto braoff = (j * bcomp + k) * ccomp * dcomp;
                    
                    // loop over angular components on ket side
                    
                    for (int32_t l = 0; l < ccomp; l++)
                    {
                        // set up pointer to R indexes on ket side
                        
                        auto idxk = ketGtoPairsBlock.getBraIdentifiers(l);
                        
                        for (int32_t m = 0; m < dcomp; m++)
                        {
                            // set up pointer to S indexes on ket side
                            
                            auto idxl = ketGtoPairsBlock.getKetIdentifiers(m);
                            
                            // set up pointer to integrals
                            
                            auto pints = spherInts.data(braoff + l * dcomp + m);
                            
                            // set up pointers to submatrices
                            
                            auto submat_pq = fockContainer.getSubMatrixData(iFockMatrix, 0, j * bcomp + k);
                            
                            auto submat_rs = fockContainer.getSubMatrixData(iFockMatrix, 1, l * dcomp + m);
                            
                            auto submat_pr = fockContainer.getSubMatrixData(iFockMatrix, 2, j * ccomp + l);
                            
                            auto submat_ps = fockContainer.getSubMatrixData(iFockMatrix, 3, j * dcomp + m);
                            
                            auto submat_qr = fockContainer.getSubMatrixData(iFockMatrix, 4, k * ccomp + l);
                            
                            auto submat_qs = fockContainer.getSubMatrixData(iFockMatrix, 5, k * dcomp + m);
                            
                            // loop over pairs on ket side
                            
                            for (int32_t n = 0; n < nKetContrPairs; n++)
                            {
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
                                    
                                    auto redr = idr - cpos[l];
                                    
                                    auto reds = ids - dpos[m];
                                    
                                    // scale integral value
                                    
                                    auto fval = pints[o + intoff];
                                    
                                    if (idp == idq) fval *= 0.5;
                                    
                                    if (idr == ids) fval *= 0.5;
                                    
                                    if ((idp == idr) && (idq == ids)) fval *= 0.5;
                                    
                                    // Coulomb contributions
                                    
                                    submat_pq[redp * bdim + redq] += 4.0 * fval * densityMatrix[idr * nDensityColumns + ids];
                                    
                                    submat_rs[redr * ddim + reds] += 4.0 * fval * densityMatrix[idp * nDensityColumns + idq];
                                    
                                    // exchange contributions
                                    
                                    submat_pr[redp * cdim + redr] -= fval * densityMatrix[idq * nDensityColumns + ids];
                                    
                                    submat_ps[redp * ddim + reds] -= fval * densityMatrix[idq * nDensityColumns + idr];
                                    
                                    submat_qr[redq * cdim + redr] -= fval * densityMatrix[idp * nDensityColumns + ids];
                                    
                                    submat_qs[redq * ddim + reds] -= fval * densityMatrix[idp * nDensityColumns + idr];
                                }
                            }
                        }
                    }
                }
            }
          
            iblock++;
        }
        
//        // loop over angular components on bra side
//
//        for (int32_t i = 0; i < acomp; i++)
//        {
//            // set up index P for bra side
//
//            int32_t idp  = (braGtoPairsBlock.getBraIdentifiers(i))[iContrPair];
//
//            int32_t redp = idp - apos[i];
//
//            // starting angular index of Q
//
//            auto jstart = (symbra) ? i : 0;
//
//            for (int32_t j = jstart; j < bcomp; j++)
//            {
//                // set up index Q for bra side
//
//                int32_t idq = (braGtoPairsBlock.getKetIdentifiers(j))[iContrPair];
//
//                int32_t redq = idq - bpos[j];
//
//                // angular offset from bra side
//
//                auto braoff = (i * bcomp + j) * ccomp * dcomp;
//
//                // loop over angular components on ket side
//
//                for (int32_t k = 0; k < ccomp; k++)
//                {
//                    // set up pointer to R indexes on ket side
//
//                    auto idxk = ketGtoPairsBlock.getBraIdentifiers(k);
//
//                    for (int32_t l = 0; l < dcomp; l++)
//                    {
//                        // set up pointer to S indexes on ket side
//
//                        auto idxl = ketGtoPairsBlock.getKetIdentifiers(l);
//
//                        // set up pointer to integrals
//
//                        auto pints = spherInts.data(braoff + k * dcomp + l);
//
//                        // set up pointers to submatrices
//
//                        auto submat_pq = fockContainer.getSubMatrixData(iFockMatrix, 0, i * bcomp + j);
//
//                        auto submat_rs = fockContainer.getSubMatrixData(iFockMatrix, 1, k * dcomp + l);
//
//                        auto submat_pr = fockContainer.getSubMatrixData(iFockMatrix, 2, i * ccomp + k);
//
//                        auto submat_ps = fockContainer.getSubMatrixData(iFockMatrix, 3, i * dcomp + l);
//
//                        auto submat_qr = fockContainer.getSubMatrixData(iFockMatrix, 4, j * ccomp + k);
//
//                        auto submat_qs = fockContainer.getSubMatrixData(iFockMatrix, 5, j * dcomp + l);
//
//                        // loop over pairs on ket side
//
//                        for (int32_t m = 0; m < nKetContrPairs; m++)
//                        {
//                            // symmetry restriction for ket angular components
//
//                            auto refr = prefk[m];
//
//                            auto refs = prefl[m];
//
//                            if ((refr == refs) && (l < k)) continue;
//
//                            // symmetry restriction for bra/ket angular componets
//
//                            bool braeqket = (refp == refr) && (refq == refs);
//
//                            if  (((k * dcomp + l) < (i * bcomp + j)) && braeqket) continue;
//
//                            // set up S and R indexes
//
//                            auto idr = idxk[m];
//
//                            auto ids = idxl[m];
//
//                            auto redr = idr - cpos[k];
//
//                            auto reds = ids - dpos[l];
//
//                            // scale integral value
//
//                            auto fval = pints[m];
//
//                            if (idp == idq) fval *= 0.5;
//
//                            if (idr == ids) fval *= 0.5;
//
//                            if ((idp == idr) && (idq == ids)) fval *= 0.5;
//
//                            // Coulomb contributions
//
//                            submat_pq[redp * bdim + redq] += 4.0 * fval * densityMatrix[idr * nDensityColumns + ids];
//
//                            submat_rs[redr * ddim + reds] += 4.0 * fval * densityMatrix[idp * nDensityColumns + idq];
//
//                            // exchange contributions
//
//                            submat_pr[redp * cdim + redr] -= fval * densityMatrix[idq * nDensityColumns + ids];
//
//                            submat_ps[redp * ddim + reds] -= fval * densityMatrix[idq * nDensityColumns + idr];
//
//                            submat_qr[redq * cdim + redr] -= fval * densityMatrix[idp * nDensityColumns + ids];
//
//                            submat_qs[redq * ddim + reds] -= fval * densityMatrix[idp * nDensityColumns + idr];
//                        }
//                    }
//                }
//            }
//        }
    }
    
} // distfock namespace
