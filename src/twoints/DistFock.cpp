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
        
        // determine symmetry of angular components on bra side
        
        auto refp = (braGtoPairsBlock.getBraIdentifiers(0))[iContrPair];
        
        auto refq = (braGtoPairsBlock.getKetIdentifiers(0))[iContrPair];
        
        bool symbra = (refp == refq);
        
        // set up pointers to reference indexes on ket side
        
        auto prefk = ketGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefl = ketGtoPairsBlock.getKetIdentifiers(0);
        
        // set up pointers to start positions for A, B, C, D centers
        
        auto apos = fockContainer.getStartPositionsA(iFockMatrix);
        
        auto bpos = fockContainer.getStartPositionsB(iFockMatrix);
        
        auto cpos = fockContainer.getStartPositionsC(iFockMatrix);
        
        auto dpos = fockContainer.getStartPositionsD(iFockMatrix);
        
        // set up dimensions of submatrices for B, C, D centers
        
        auto bdim = fockContainer.getDimensionsB(iFockMatrix);
        
        auto cdim = fockContainer.getDimensionsC(iFockMatrix);
        
        auto ddim = fockContainer.getDimensionsD(iFockMatrix);
        
        // loop over angular components on bra side
        
        for (int32_t i = 0; i < acomp; i++)
        {
            // set up index P for bra side
            
            int32_t idp  = (braGtoPairsBlock.getBraIdentifiers(i))[iContrPair];
            
            int32_t redp = idp - apos[i];
            
            // starting angular index of Q
            
            auto jstart = (symbra) ? i : 0;
            
            for (int32_t j = jstart; j < bcomp; j++)
            {
                // set up bra angular index
                
                auto braidx = i * bcomp + j;
                
                auto braoff = braidx * ccomp * dcomp;
                
                // set up index Q for bra side
                
                int32_t idq = (braGtoPairsBlock.getKetIdentifiers(j))[iContrPair];
                
                int32_t redq = idq - bpos[j];
            
                // loop over angular components on ket side
                
                for (int32_t k = 0; k < ccomp; k++)
                {
                    // set up pointer to R indexes on ket side
                    
                    auto idxk = ketGtoPairsBlock.getBraIdentifiers(k);
                    
                    for (int32_t l = 0; l < dcomp; l++)
                    {
                        // set up ket angular index
                        
                        auto ketidx = k * dcomp + l;
                        
                        // set up pointer to S indexes on ket side
                        
                        auto idxl = ketGtoPairsBlock.getKetIdentifiers(l);
                        
                        // set up pointer to integrals
                        
                        auto pints = spherInts.data(braoff + ketidx);
                        
                        // set up pointers to submatrices
                        
                        auto submat_pq = fockContainer.getSubMatrixData(iFockMatrix, 0, i * bcomp + j);
                        
                        auto submat_rs = fockContainer.getSubMatrixData(iFockMatrix, 1, k * dcomp + l);
                        
                        auto submat_pr = fockContainer.getSubMatrixData(iFockMatrix, 2, i * ccomp + k);
                        
                        auto submat_ps = fockContainer.getSubMatrixData(iFockMatrix, 3, i * dcomp + l);
                        
                        auto submat_qr = fockContainer.getSubMatrixData(iFockMatrix, 4, j * ccomp + k);
                        
                        auto submat_qs = fockContainer.getSubMatrixData(iFockMatrix, 5, j * dcomp + l);
                        
                        // loop over pairs on ket side
                        
                        for (int32_t m = 0; m < nKetContrPairs; m++)
                        {
                            // symmetry restriction for ket angular components
                            
                            auto refr = prefk[m];
                            
                            auto refs = prefl[m];
                            
                            if ((refr == refs) && (l < k)) continue;
                            
                            // symmetry restriction for bra/ket angular componets
                            
                            bool braeqket = (refp == refr) && (refq == refs);
                            
                            if  ((ketidx < braidx) && braeqket) continue;
                
                            // set up S and R indexes
                            
                            auto idr = idxk[m];
                            
                            auto ids = idxl[m];
                            
                            auto redr = idr - cpos[k];
                            
                            auto reds = ids - dpos[l];
                            
                            // scale integral value
                            
                            auto fval = pints[m];
                                
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
    
    void
    distRestJ(      CFockContainer&      fockContainer,
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
        
        // determine symmetry of angular components on bra side
        
        auto refp = (braGtoPairsBlock.getBraIdentifiers(0))[iContrPair];
        
        auto refq = (braGtoPairsBlock.getKetIdentifiers(0))[iContrPair];
        
        bool symbra = (refp == refq);
        
        // set up pointers to reference indexes on ket side
        
        auto prefk = ketGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefl = ketGtoPairsBlock.getKetIdentifiers(0);
        
        // set up pointers to start positions for A, B, C, D centers
        
        auto apos = fockContainer.getStartPositionsA(iFockMatrix);
        
        auto bpos = fockContainer.getStartPositionsB(iFockMatrix);
        
        auto cpos = fockContainer.getStartPositionsC(iFockMatrix);
        
        auto dpos = fockContainer.getStartPositionsD(iFockMatrix);
        
        // set up dimensions of submatrices for B, D centers
        
        auto bdim = fockContainer.getDimensionsB(iFockMatrix);
        
        auto ddim = fockContainer.getDimensionsD(iFockMatrix);
        
        // loop over angular components on bra side
        
        for (int32_t i = 0; i < acomp; i++)
        {
            // set up index P for bra side
            
            int32_t idp  = (braGtoPairsBlock.getBraIdentifiers(i))[iContrPair];
            
            int32_t redp = idp - apos[i];
            
            // starting angular index of Q
            
            auto jstart = (symbra) ? i : 0;
            
            for (int32_t j = jstart; j < bcomp; j++)
            {
                // set up bra angular index
                
                auto braidx = i * bcomp + j;
                
                auto braoff = braidx * ccomp * dcomp;
                
                // set up index Q for bra side
                
                int32_t idq = (braGtoPairsBlock.getKetIdentifiers(j))[iContrPair];
                
                int32_t redq = idq - bpos[j];
                
                // loop over angular components on ket side
                
                for (int32_t k = 0; k < ccomp; k++)
                {
                    // set up pointer to R indexes on ket side
                    
                    auto idxk = ketGtoPairsBlock.getBraIdentifiers(k);
                    
                    for (int32_t l = 0; l < dcomp; l++)
                    {
                        // set up ket angular index
                        
                        auto ketidx = k * dcomp + l;
                        
                        // set up pointer to S indexes on ket side
                        
                        auto idxl = ketGtoPairsBlock.getKetIdentifiers(l);
                        
                        // set up pointer to integrals
                        
                        auto pints = spherInts.data(braoff + ketidx);
                        
                        // set up pointers to submatrices
                        
                        auto submat_pq = fockContainer.getSubMatrixData(iFockMatrix, 0, i * bcomp + j);
                        
                        auto submat_rs = fockContainer.getSubMatrixData(iFockMatrix, 1, k * dcomp + l);
                        
                        // loop over pairs on ket side
                        
                        for (int32_t m = 0; m < nKetContrPairs; m++)
                        {
                            // symmetry restriction for ket angular components
                            
                            auto refr = prefk[m];
                            
                            auto refs = prefl[m];
                            
                            if ((refr == refs) && (l < k)) continue;
                            
                            // symmetry restriction for bra/ket angular componets
                            
                            bool braeqket = (refp == refr) && (refq == refs);
                            
                            if  ((ketidx < braidx) && braeqket) continue;
                            
                            // set up S and R indexes
                            
                            auto idr = idxk[m];
                            
                            auto ids = idxl[m];
                            
                            auto redr = idr - cpos[k];
                            
                            auto reds = ids - dpos[l];
                            
                            // scale integral value
                            
                            auto fval = pints[m];
                            
                            if (idp == idq) fval *= 0.5;
                            
                            if (idr == ids) fval *= 0.5;
                            
                            if ((idp == idr) && (idq == ids)) fval *= 0.5;
                            
                            // Coulomb contributions
                            
                            submat_pq[redp * bdim + redq] += 2.0 * fval * densityMatrix[idr * nDensityColumns + ids];
                            
                            submat_rs[redr * ddim + reds] += 2.0 * fval * densityMatrix[idp * nDensityColumns + idq];
                        }
                    }
                }
            }
        }
    }
    
    void
    distRestK(      CFockContainer&      fockContainer,
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
        
        // determine symmetry of angular components on bra side
        
        auto refp = (braGtoPairsBlock.getBraIdentifiers(0))[iContrPair];
        
        auto refq = (braGtoPairsBlock.getKetIdentifiers(0))[iContrPair];
        
        bool symbra = (refp == refq);
        
        // set up pointers to reference indexes on ket side
        
        auto prefk = ketGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefl = ketGtoPairsBlock.getKetIdentifiers(0);
        
        // set up pointers to start positions for A, B, C, D centers
        
        auto apos = fockContainer.getStartPositionsA(iFockMatrix);
        
        auto bpos = fockContainer.getStartPositionsB(iFockMatrix);
        
        auto cpos = fockContainer.getStartPositionsC(iFockMatrix);
        
        auto dpos = fockContainer.getStartPositionsD(iFockMatrix);
        
        // set up dimensions of submatrices for C, D centers
        
        auto cdim = fockContainer.getDimensionsC(iFockMatrix);
        
        auto ddim = fockContainer.getDimensionsD(iFockMatrix);
        
        // loop over angular components on bra side
        
        for (int32_t i = 0; i < acomp; i++)
        {
            // set up index P for bra side
            
            int32_t idp  = (braGtoPairsBlock.getBraIdentifiers(i))[iContrPair];
            
            int32_t redp = idp - apos[i];
            
            // starting angular index of Q
            
            auto jstart = (symbra) ? i : 0;
            
            for (int32_t j = jstart; j < bcomp; j++)
            {
                // set up bra angular index
                
                auto braidx = i * bcomp + j;
                
                auto braoff = braidx * ccomp * dcomp;
                
                // set up index Q for bra side
                
                int32_t idq = (braGtoPairsBlock.getKetIdentifiers(j))[iContrPair];
                
                int32_t redq = idq - bpos[j];
                
                // loop over angular components on ket side
                
                for (int32_t k = 0; k < ccomp; k++)
                {
                    // set up pointer to R indexes on ket side
                    
                    auto idxk = ketGtoPairsBlock.getBraIdentifiers(k);
                    
                    for (int32_t l = 0; l < dcomp; l++)
                    {
                        // set up ket angular index
                        
                        auto ketidx = k * dcomp + l;
                        
                        // set up pointer to S indexes on ket side
                        
                        auto idxl = ketGtoPairsBlock.getKetIdentifiers(l);
                        
                        // set up pointer to integrals
                        
                        auto pints = spherInts.data(braoff + ketidx);
                        
                        // set up pointers to submatrices
                        
                        auto submat_pr = fockContainer.getSubMatrixData(iFockMatrix, 0, i * ccomp + k);
                        
                        auto submat_ps = fockContainer.getSubMatrixData(iFockMatrix, 1, i * dcomp + l);
                        
                        auto submat_qr = fockContainer.getSubMatrixData(iFockMatrix, 2, j * ccomp + k);
                        
                        auto submat_qs = fockContainer.getSubMatrixData(iFockMatrix, 3, j * dcomp + l);
                        
                        // loop over pairs on ket side
                        
                        for (int32_t m = 0; m < nKetContrPairs; m++)
                        {
                            // symmetry restriction for ket angular components
                            
                            auto refr = prefk[m];
                            
                            auto refs = prefl[m];
                            
                            if ((refr == refs) && (l < k)) continue;
                            
                            // symmetry restriction for bra/ket angular componets
                            
                            bool braeqket = (refp == refr) && (refq == refs);
                            
                            if  ((ketidx < braidx) && braeqket) continue;
                            
                            // set up S and R indexes
                            
                            auto idr = idxk[m];
                            
                            auto ids = idxl[m];
                            
                            auto redr = idr - cpos[k];
                            
                            auto reds = ids - dpos[l];
                            
                            // scale integral value
                            
                            auto fval = pints[m];
                            
                            if (idp == idq) fval *= 0.5;
                            
                            if (idr == ids) fval *= 0.5;
                            
                            if ((idp == idr) && (idq == ids)) fval *= 0.5;
                            
                            // exchange contributions
                            
                            submat_pr[redp * cdim + redr] += fval * densityMatrix[idq * nDensityColumns + ids];
                            
                            submat_ps[redp * ddim + reds] += fval * densityMatrix[idq * nDensityColumns + idr];
                            
                            submat_qr[redq * cdim + redr] += fval * densityMatrix[idp * nDensityColumns + ids];
                            
                            submat_qs[redq * ddim + reds] += fval * densityMatrix[idp * nDensityColumns + idr];
                        }
                    }
                }
            }
        }
    }
    
    void
    distRestGenJ(      CFockContainer&      fockContainer,
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
        
        // determine symmetry of angular components on bra side
        
        auto refp = (braGtoPairsBlock.getBraIdentifiers(0))[iContrPair];
        
        auto refq = (braGtoPairsBlock.getKetIdentifiers(0))[iContrPair];
        
        bool symbra = (refp == refq);
        
        // set up pointers to reference indexes on ket side
        
        auto prefk = ketGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefl = ketGtoPairsBlock.getKetIdentifiers(0);
        
        // set up pointers to start positions for A, B, C, D centers
        
        auto apos = fockContainer.getStartPositionsA(iFockMatrix);
        
        auto bpos = fockContainer.getStartPositionsB(iFockMatrix);
        
        auto cpos = fockContainer.getStartPositionsC(iFockMatrix);
        
        auto dpos = fockContainer.getStartPositionsD(iFockMatrix);
        
        // set up dimensions of submatrices for A, B, D, C centers
        
        auto adim = fockContainer.getDimensionsA(iFockMatrix);
        
        auto bdim = fockContainer.getDimensionsB(iFockMatrix);
        
        auto cdim = fockContainer.getDimensionsC(iFockMatrix);
        
        auto ddim = fockContainer.getDimensionsD(iFockMatrix);
        
        // loop over angular components on bra side
        
        for (int32_t i = 0; i < acomp; i++)
        {
            // set up index P for bra side
            
            int32_t idp  = (braGtoPairsBlock.getBraIdentifiers(i))[iContrPair];
            
            int32_t redp = idp - apos[i];
            
            // starting angular index of Q
            
            auto jstart = (symbra) ? i : 0;
            
            for (int32_t j = jstart; j < bcomp; j++)
            {
                // set up bra angular index
                
                auto braidx = i * bcomp + j;
                
                auto braoff = braidx * ccomp * dcomp;
                
                // set up index Q for bra side
                
                int32_t idq = (braGtoPairsBlock.getKetIdentifiers(j))[iContrPair];
                
                int32_t redq = idq - bpos[j];
                
                // loop over angular components on ket side
                
                for (int32_t k = 0; k < ccomp; k++)
                {
                    // set up pointer to R indexes on ket side
                    
                    auto idxk = ketGtoPairsBlock.getBraIdentifiers(k);
                    
                    for (int32_t l = 0; l < dcomp; l++)
                    {
                        // set up ket angular index
                        
                        auto ketidx = k * dcomp + l;
                        
                        // set up pointer to S indexes on ket side
                        
                        auto idxl = ketGtoPairsBlock.getKetIdentifiers(l);
                        
                        // set up pointer to integrals
                        
                        auto pints = spherInts.data(braoff + ketidx);
                        
                        // set up pointers to submatrices
                        
                        auto submat_pq = fockContainer.getSubMatrixData(iFockMatrix, 0, i * bcomp + j);
                        
                        auto submat_qp = fockContainer.getSubMatrixData(iFockMatrix, 1, j * acomp + i);
                        
                        auto submat_rs = fockContainer.getSubMatrixData(iFockMatrix, 2, k * dcomp + l);
                        
                        auto submat_sr = fockContainer.getSubMatrixData(iFockMatrix, 3, l * ccomp + k);
                        
                        // loop over pairs on ket side
                        
                        for (int32_t m = 0; m < nKetContrPairs; m++)
                        {
                            // symmetry restriction for ket angular components
                            
                            auto refr = prefk[m];
                            
                            auto refs = prefl[m];
                            
                            if ((refr == refs) && (l < k)) continue;
                            
                            // symmetry restriction for bra/ket angular componets
                            
                            bool braeqket = (refp == refr) && (refq == refs);
                            
                            if  ((ketidx < braidx) && braeqket) continue;
                            
                            // set up S and R indexes
                            
                            auto idr = idxk[m];
                            
                            auto ids = idxl[m];
                            
                            auto redr = idr - cpos[k];
                            
                            auto reds = ids - dpos[l];
                            
                            // scale integral value
                            
                            auto fval = pints[m];
                            
                            if (idp == idq) fval *= 0.5;
                            
                            if (idr == ids) fval *= 0.5;
                            
                            if ((idp == idr) && (idq == ids)) fval *= 0.5;
                            
                            // Coulomb contributions
                            
                            auto f2rs = fval * (densityMatrix[idr * nDensityColumns + ids] +
                                                
                                                densityMatrix[ids * nDensityColumns + idr]);
                            
                            submat_pq[redp * bdim + redq] += f2rs;
                            
                            submat_qp[redq * adim + redp] += f2rs;
                            
                            auto f2pq = fval * (densityMatrix[idp * nDensityColumns + idq] +
                                                
                                                densityMatrix[idq * nDensityColumns + idp]);
                            
                            submat_rs[redr * ddim + reds] += f2pq;
                            
                            submat_sr[reds * cdim + redr] += f2pq;
                        }
                    }
                }
            }
        }
    }
    
    void
    distRestGenK(      CFockContainer&      fockContainer,
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
        
        // determine symmetry of angular components on bra side
        
        auto refp = (braGtoPairsBlock.getBraIdentifiers(0))[iContrPair];
        
        auto refq = (braGtoPairsBlock.getKetIdentifiers(0))[iContrPair];
        
        bool symbra = (refp == refq);
        
        // set up pointers to reference indexes on ket side
        
        auto prefk = ketGtoPairsBlock.getBraIdentifiers(0);
        
        auto prefl = ketGtoPairsBlock.getKetIdentifiers(0);
        
        // set up pointers to start positions for A, B, C, D centers
        
        auto apos = fockContainer.getStartPositionsA(iFockMatrix);
        
        auto bpos = fockContainer.getStartPositionsB(iFockMatrix);
        
        auto cpos = fockContainer.getStartPositionsC(iFockMatrix);
        
        auto dpos = fockContainer.getStartPositionsD(iFockMatrix);
        
        // set up dimensions of submatrices for A, B, C, D centers
        
        auto adim = fockContainer.getDimensionsA(iFockMatrix);
        
        auto bdim = fockContainer.getDimensionsB(iFockMatrix);
        
        auto cdim = fockContainer.getDimensionsC(iFockMatrix);
        
        auto ddim = fockContainer.getDimensionsD(iFockMatrix);
        
        // loop over angular components on bra side
        
        for (int32_t i = 0; i < acomp; i++)
        {
            // set up index P for bra side
            
            int32_t idp  = (braGtoPairsBlock.getBraIdentifiers(i))[iContrPair];
            
            int32_t redp = idp - apos[i];
            
            // starting angular index of Q
            
            auto jstart = (symbra) ? i : 0;
            
            for (int32_t j = jstart; j < bcomp; j++)
            {
                // set up bra angular index
                
                auto braidx = i * bcomp + j;
                
                auto braoff = braidx * ccomp * dcomp;
                
                // set up index Q for bra side
                
                int32_t idq = (braGtoPairsBlock.getKetIdentifiers(j))[iContrPair];
                
                int32_t redq = idq - bpos[j];
                
                // loop over angular components on ket side
                
                for (int32_t k = 0; k < ccomp; k++)
                {
                    // set up pointer to R indexes on ket side
                    
                    auto idxk = ketGtoPairsBlock.getBraIdentifiers(k);
                    
                    for (int32_t l = 0; l < dcomp; l++)
                    {
                        // set up ket angular index
                        
                        auto ketidx = k * dcomp + l;
                        
                        // set up pointer to S indexes on ket side
                        
                        auto idxl = ketGtoPairsBlock.getKetIdentifiers(l);
                        
                        // set up pointer to integrals
                        
                        auto pints = spherInts.data(braoff + ketidx);
                        
                        // set up pointers to submatrices
                        
                        auto submat_pr = fockContainer.getSubMatrixData(iFockMatrix, 0, i * ccomp + k);
                        
                        auto submat_rp = fockContainer.getSubMatrixData(iFockMatrix, 1, k * acomp + i);
                        
                        auto submat_ps = fockContainer.getSubMatrixData(iFockMatrix, 2, i * dcomp + l);
                        
                        auto submat_sp = fockContainer.getSubMatrixData(iFockMatrix, 3, l * acomp + i);
                        
                        auto submat_qr = fockContainer.getSubMatrixData(iFockMatrix, 4, j * ccomp + k);
                        
                        auto submat_rq = fockContainer.getSubMatrixData(iFockMatrix, 5, k * bcomp + j);
                        
                        auto submat_qs = fockContainer.getSubMatrixData(iFockMatrix, 6, j * dcomp + l);
                        
                        auto submat_sq = fockContainer.getSubMatrixData(iFockMatrix, 7, l * bcomp + j);
                        
                        // loop over pairs on ket side
                        
                        for (int32_t m = 0; m < nKetContrPairs; m++)
                        {
                            // symmetry restriction for ket angular components
                            
                            auto refr = prefk[m];
                            
                            auto refs = prefl[m];
                            
                            if ((refr == refs) && (l < k)) continue;
                            
                            // symmetry restriction for bra/ket angular componets
                            
                            bool braeqket = (refp == refr) && (refq == refs);
                            
                            if  ((ketidx < braidx) && braeqket) continue;
                            
                            // set up S and R indexes
                            
                            auto idr = idxk[m];
                            
                            auto ids = idxl[m];
                            
                            auto redr = idr - cpos[k];
                            
                            auto reds = ids - dpos[l];
                            
                            // scale integral value
                            
                            auto fval = pints[m];
                            
                            if (idp == idq) fval *= 0.5;
                            
                            if (idr == ids) fval *= 0.5;
                            
                            if ((idp == idr) && (idq == ids)) fval *= 0.5;
                            
                            // exchange contributions
                            
                            submat_pr[redp * cdim + redr] += fval * densityMatrix[idq * nDensityColumns + ids];
                            
                            submat_rp[redr * adim + redp] += fval * densityMatrix[ids * nDensityColumns + idq];
                            
                            submat_ps[redp * ddim + reds] += fval * densityMatrix[idq * nDensityColumns + idr];
                            
                            submat_sp[reds * adim + redp] += fval * densityMatrix[idr * nDensityColumns + idq];
                            
                            submat_qr[redq * cdim + redr] += fval * densityMatrix[idp * nDensityColumns + ids];
                            
                            submat_rq[redr * bdim + redq] += fval * densityMatrix[ids * nDensityColumns + idp];
                            
                            submat_qs[redq * ddim + reds] += fval * densityMatrix[idp * nDensityColumns + idr];
                            
                            submat_sq[reds * bdim + redq] += fval * densityMatrix[idr * nDensityColumns + idp];
                        }
                    }
                }
            }
        }
    }
    
} // distfock namespace
