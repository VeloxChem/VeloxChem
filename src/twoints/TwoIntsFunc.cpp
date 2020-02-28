//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "TwoIntsFunc.hpp"

namespace twointsfunc { // twointsfunc namespace
    
    void
    compDistancesAQ(      CMemBlock2D<double>& aqDistances,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto brx = braGtoBlock.getCoordinatesX();
        
        auto bry = braGtoBlock.getCoordinatesY();
        
        auto brz = braGtoBlock.getCoordinatesZ();
        
        auto spos = braGtoBlock.getStartPositions();
        
        // compute distances
        
        mathfunc::distances(aqDistances.data(0),
                            aqDistances.data(1),
                            aqDistances.data(2),
                            brx[spos[iContrGto]],
                            bry[spos[iContrGto]],
                            brz[spos[iContrGto]],
                            ketGtoPairsBlock.getCoordinatesPX(),
                            ketGtoPairsBlock.getCoordinatesPY(),
                            ketGtoPairsBlock.getCoordinatesPZ(),
                            ketGtoPairsBlock.getNumberOfScreenedPrimPairs());
    }

    void
    compFactorsForThreeCenterElectronRepulsion(      CMemBlock2D<double>& osFactors,
                                               const CGtoBlock&           braGtoBlock,
                                               const CGtoPairsBlock&      ketGtoPairsBlock,
                                               const int32_t              iContrGto)
    {
        // set up angular momentum data
        
        auto bang = braGtoBlock.getAngularMomentum();
        
        auto kang = ketGtoPairsBlock.getBraAngularMomentum()
        
                  + ketGtoPairsBlock.getKetAngularMomentum();
        
        // set up pointers to primitives data on bra side
        
        auto bexp = braGtoBlock.getExponents();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitive pairs data on ket side
        
        auto kfxi = ketGtoPairsBlock.getFactorsXi();
        
        auto koxi = ketGtoPairsBlock.getFactorsOneOverXi();
        
        auto nprim = ketGtoPairsBlock.getNumberOfScreenedPrimPairs();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(5 * idx);
            
            auto fz = osFactors.data(5 * idx + 1);
            
            auto fb = bexp[i];
            
            #pragma omp simd aligned(fx, fz, kfxi: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fx[j] = 1.0 / (fb + kfxi[j]);
                
                fz[j] = fb * kfxi[j] * fx[j];
            }
            
            if (bang > 1)
            {
                auto ga = osFactors.data(5 * idx + 2);
                
                auto ta = osFactors.data(5 * idx + 3);
                
                auto fga = 1.0 / fb;
                
                #pragma omp simd aligned(ga, ta, fz: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    ga[j] = fga;
                    
                    ta[j] = fga * fz[j];
                }
            }
            
            if (kang > 1)
            {
                auto td = osFactors.data(5 * idx + 4);
                
                #pragma omp simd aligned(td, koxi, fz: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    td[j] = koxi[j] * fz[j];
                }
            }
            
            idx++;
        }
    }
    
    void
    compCoordinatesForW(      CMemBlock2D<double>& wCoordinates,
                        const CMemBlock2D<double>& osFactors,
                        const int32_t              nFactors,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoPairsBlock&      ketGtoPairsBlock,
                        const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto brx = braGtoBlock.getCoordinatesX();
        
        auto bry = braGtoBlock.getCoordinatesY();
        
        auto brz = braGtoBlock.getCoordinatesZ();
        
        auto bexp = braGtoBlock.getExponents();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto kqx = ketGtoPairsBlock.getCoordinatesPX();
        
        auto kqy = ketGtoPairsBlock.getCoordinatesPY();
        
        auto kqz = ketGtoPairsBlock.getCoordinatesPZ();
        
        auto kfxi = ketGtoPairsBlock.getFactorsXi();
        
        auto nprim = ketGtoPairsBlock.getNumberOfScreenedPrimPairs();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to prefactors
            
            auto fx = osFactors.data(nFactors * idx);
            
            // set up primitive GTO data on bra side
            
            auto fax = bexp[i] * brx[i];
            
            auto fay = bexp[i] * bry[i];
            
            auto faz = bexp[i] * brz[i];
            
            // set up pointers to coordinates of P
            
            auto wx = wCoordinates.data(3 * idx);
            
            auto wy = wCoordinates.data(3 * idx + 1);
            
            auto wz = wCoordinates.data(3 * idx + 2);
            
            #pragma omp simd aligned(fx, wx, wy, wz, kfxi, kqx, kqy,\
                                     kqz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = fx[j];
                
                wx[j] = fact * (fax + kfxi[j] * kqx[j]);
                
                wy[j] = fact * (fay + kfxi[j] * kqy[j]);
                
                wz[j] = fact * (faz + kfxi[j] * kqz[j]);
            }
            
            idx++;
        }
    }
    
    void
    compDistancesWA(      CMemBlock2D<double>& waDistances,
                    const CMemBlock2D<double>& wCoordinates,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const int32_t              iContrGto)
    {
        // skip computation for zero angular momentum on bra side
        
        if (braGtoBlock.getAngularMomentum() == 0) return;
        
        // set up pointers to primitives data on bra side
        
        auto brx = braGtoBlock.getCoordinatesX();
        
        auto bry = braGtoBlock.getCoordinatesY();
        
        auto brz = braGtoBlock.getCoordinatesZ();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoPairsBlock.getNumberOfScreenedPrimPairs();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of W
            
            auto wx = wCoordinates.data(3 * idx);
            
            auto wy = wCoordinates.data(3 * idx + 1);
            
            auto wz = wCoordinates.data(3 * idx + 2);
            
            // set up pmitive GTO data on bra side
            
            auto ax = brx[i];
            
            auto ay = bry[i];
            
            auto az = brz[i];
            
            // set up pointers to distances R(WA)
            
            auto wax = waDistances.data(3 * idx);
            
            auto way = waDistances.data(3 * idx + 1);
            
            auto waz = waDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(wx, wy, wz, wax, way, waz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                wax[j] = wx[j] - ax;
                
                way[j] = wy[j] - ay;
                
                waz[j] = wz[j] - az;
            }
            
            idx++;
        }
    }
    
    void
    compDistancesWD(      CMemBlock2D<double>& wdDistances,
                    const CMemBlock2D<double>& wCoordinates,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const int32_t              iContrGto)
    {
        // skip computation for zero angular momentum on ket side
        
        if ((ketGtoPairsBlock.getBraAngularMomentum() == 0) &&
            (ketGtoPairsBlock.getKetAngularMomentum() == 0)) return;
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto rdx = ketGtoPairsBlock.getCoordinatesBX();
        
        auto rdy = ketGtoPairsBlock.getCoordinatesBY();
        
        auto rdz = ketGtoPairsBlock.getCoordinatesBZ();
        
        auto nprim = ketGtoPairsBlock.getNumberOfScreenedPrimPairs();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of W
            
            auto wx = wCoordinates.data(3 * idx);
            
            auto wy = wCoordinates.data(3 * idx + 1);
            
            auto wz = wCoordinates.data(3 * idx + 2);
            
            // set up pointers to distances R(WD)
            
            auto wdx = wdDistances.data(3 * idx);
            
            auto wdy = wdDistances.data(3 * idx + 1);
            
            auto wdz = wdDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(wx, wy, wz, wdx, wdy, wdz, rdx, rdy,\
                                     rdz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                wdx[j] = wx[j] - rdx[j];
                
                wdy[j] = wy[j] - rdy[j];
                
                wdz[j] = wz[j] - rdz[j];
            }
            
            idx++;
        }
    }
    
    void
    compDistancesWQ(      CMemBlock2D<double>& wqDistances,
                    const CMemBlock2D<double>& wCoordinates,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const int32_t              iContrGto)
    {
        // skip computation for zero angular momentum on ket side
        
        if ((ketGtoPairsBlock.getBraAngularMomentum() == 0) &&
            (ketGtoPairsBlock.getKetAngularMomentum() == 0)) return;
    
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto rpx = ketGtoPairsBlock.getCoordinatesPX();
        
        auto rpy = ketGtoPairsBlock.getCoordinatesPY();
        
        auto rpz = ketGtoPairsBlock.getCoordinatesPZ();
        
        auto nprim = ketGtoPairsBlock.getNumberOfScreenedPrimPairs();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of W
            
            auto wx = wCoordinates.data(3 * idx);
            
            auto wy = wCoordinates.data(3 * idx + 1);
            
            auto wz = wCoordinates.data(3 * idx + 2);
            
            // set up pointers to distances R(WQ)
            
            auto wqx = wqDistances.data(3 * idx);
            
            auto wqy = wqDistances.data(3 * idx + 1);
            
            auto wqz = wqDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(wx, wy, wz, wqx, wqy, wqz, rpx, rpy,\
                                     rpz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                wqx[j] = wx[j] - rpx[j];
                
                wqy[j] = wy[j] - rpy[j];
                
                wqz[j] = wz[j] - rpz[j];
            }
            
            idx++;
        }
    }
    
    void
    compDistancesPQ(      CMemBlock2D<double>& pqDistances,
                    const CGtoPairsBlock&      braGtoPairsBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const int32_t              nKetPrimPairs,
                    const int32_t              iContrPair)
    {
        // set up pointers to primitive pairs data on bra side
        
        auto brx = braGtoPairsBlock.getCoordinatesPX();
        
        auto bry = braGtoPairsBlock.getCoordinatesPY();
        
        auto brz = braGtoPairsBlock.getCoordinatesPZ();
        
        auto spos = braGtoPairsBlock.getStartPositions();
        
        auto epos = braGtoPairsBlock.getEndPositions();
        
        // loop over componets of contracted GTOs pair
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
        {
            // compute distances
            
            mathfunc::distances(pqDistances.data(3 * idx),
                                pqDistances.data(3 * idx + 1),
                                pqDistances.data(3 * idx + 2),
                                brx[i], bry[i], brz[i],
                                ketGtoPairsBlock.getCoordinatesPX(),
                                ketGtoPairsBlock.getCoordinatesPY(),
                                ketGtoPairsBlock.getCoordinatesPZ(),
                                nKetPrimPairs);
            
            idx++; 
        }
    }
    
    void
    compFactorsForElectronRepulsion(      CMemBlock2D<double>& osFactors,
                                    const CGtoPairsBlock&      braGtoPairsBlock,
                                    const CGtoPairsBlock&      ketGtoPairsBlock,
                                    const int32_t              nKetPrimPairs,
                                    const int32_t              iContrPair)
    {
        // set up angular momentum data
        
        auto bang = braGtoPairsBlock.getBraAngularMomentum()
        
                  + braGtoPairsBlock.getKetAngularMomentum();
        
        auto kang = ketGtoPairsBlock.getBraAngularMomentum()
        
                  + ketGtoPairsBlock.getKetAngularMomentum();
        
        // set up pointers to primitive pairs data on bra side
        
        auto bfxi = braGtoPairsBlock.getFactorsXi();
        
        auto boxi = braGtoPairsBlock.getFactorsOneOverXi();
        
        auto spos = braGtoPairsBlock.getStartPositions();
        
        auto epos = braGtoPairsBlock.getEndPositions();
        
         // set up pointers to primitive pairs data on ket side
        
        auto kfxi = ketGtoPairsBlock.getFactorsXi();
        
        auto koxi = ketGtoPairsBlock.getFactorsOneOverXi();
        
        // loop over componets of contracted GTOs pair
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(4 * idx);
            
            auto fz = osFactors.data(4 * idx + 1);
            
            auto fb = bfxi[i];
            
            #pragma omp simd aligned(fx, fz, kfxi: VLX_ALIGN)
            for (int32_t j = 0; j < nKetPrimPairs; j++)
            {
                fx[j] = 1.0 / (fb + kfxi[j]);
                
                fz[j] = fb * kfxi[j] * fx[j];
            }
            
            if (bang > 1)
            {
                auto ta = osFactors.data(4 * idx + 2);
                
                auto fga = boxi[i];
                
                #pragma omp simd aligned(ta, fz: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    ta[j] = fga * fz[j];
                }
            }
            
            if (kang > 1)
            {
                auto td = osFactors.data(4 * idx + 3);
                
                #pragma omp simd aligned(td, koxi, fz: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    td[j] = koxi[j] * fz[j];
                }
            }
            
            idx++;
        }
    }
    
    void
    compCoordinatesForW(      CMemBlock2D<double>& wCoordinates,
                        const CMemBlock2D<double>& osFactors,
                        const int32_t              nFactors,
                        const CGtoPairsBlock&      braGtoPairsBlock,
                        const CGtoPairsBlock&      ketGtoPairsBlock,
                        const int32_t              nKetPrimPairs,
                        const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side
        
        auto bpx = braGtoPairsBlock.getCoordinatesPX();
        
        auto bpy = braGtoPairsBlock.getCoordinatesPY();
        
        auto bpz = braGtoPairsBlock.getCoordinatesPZ();
        
        auto bfxi = braGtoPairsBlock.getFactorsXi();
        
        auto spos = braGtoPairsBlock.getStartPositions();
        
        auto epos = braGtoPairsBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto kqx = ketGtoPairsBlock.getCoordinatesPX();
        
        auto kqy = ketGtoPairsBlock.getCoordinatesPY();
        
        auto kqz = ketGtoPairsBlock.getCoordinatesPZ();
        
        auto kfxi = ketGtoPairsBlock.getFactorsXi();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
        {
            // set up pointers to prefactors
            
            auto fx = osFactors.data(nFactors * idx);
            
            // set up primitive GTO pair data on bra side
            
            auto fax = bfxi[i] * bpx[i];
            
            auto fay = bfxi[i] * bpy[i];
            
            auto faz = bfxi[i] * bpz[i];
            
            // set up pointers to coordinates of W
            
            auto wx = wCoordinates.data(3 * idx);
            
            auto wy = wCoordinates.data(3 * idx + 1);
            
            auto wz = wCoordinates.data(3 * idx + 2);
            
            #pragma omp simd aligned(fx, wx, wy, wz, kfxi, kqx, kqy,\
                                     kqz: VLX_ALIGN)
            for (int32_t j = 0; j < nKetPrimPairs; j++)
            {
                double fact = fx[j];
                
                wx[j] = fact * (fax + kfxi[j] * kqx[j]);
                
                wy[j] = fact * (fay + kfxi[j] * kqy[j]);
                
                wz[j] = fact * (faz + kfxi[j] * kqz[j]);
            }
            
            idx++;
        }
    }
    
    void
    compDistancesWP(      CMemBlock2D<double>& wpDistances,
                    const CMemBlock2D<double>& wCoordinates,
                    const CGtoPairsBlock&      braGtoPairsBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const int32_t              nKetPrimPairs,
                    const int32_t              iContrPair)
    {
        // skip computation for zero angular momentum on bra side
        
        if ((braGtoPairsBlock.getBraAngularMomentum() == 0) &&
            (braGtoPairsBlock.getKetAngularMomentum() == 0)) return;
        
        // set up pointers to primitive pairs data on bra side
        
        auto rpx = braGtoPairsBlock.getCoordinatesPX();
        
        auto rpy = braGtoPairsBlock.getCoordinatesPY();
        
        auto rpz = braGtoPairsBlock.getCoordinatesPZ();
        
        auto spos = braGtoPairsBlock.getStartPositions();
        
        auto epos = braGtoPairsBlock.getEndPositions();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
        {
            // set up pointers to coordinates of W
            
            auto wx = wCoordinates.data(3 * idx);
            
            auto wy = wCoordinates.data(3 * idx + 1);
            
            auto wz = wCoordinates.data(3 * idx + 2);
            
            // set up pointers to distances R(WP)
            
            auto wpx = wpDistances.data(3 * idx);
            
            auto wpy = wpDistances.data(3 * idx + 1);
            
            auto wpz = wpDistances.data(3 * idx + 2);
            
            // set up distances
            
            auto cpx = rpx[i];
            
            auto cpy = rpy[i];
            
            auto cpz = rpz[i];
            
            #pragma omp simd aligned(wx, wy, wz, wpx, wpy, wpz: VLX_ALIGN)
            for (int32_t j = 0; j < nKetPrimPairs; j++)
            {
                wpx[j] = wx[j] - cpx;
                
                wpy[j] = wy[j] - cpy;
                
                wpz[j] = wz[j] - cpz;
            }
            
            idx++;
        }
    }
    
    void
    compDistancesWQ(      CMemBlock2D<double>& wqDistances,
                    const CMemBlock2D<double>& wCoordinates,
                    const CGtoPairsBlock&      braGtoPairsBlock,
                    const CGtoPairsBlock&      ketGtoPairsBlock,
                    const int32_t              nKetPrimPairs,
                    const int32_t              iContrPair)
    {
        // skip computation for zero angular momentum on ket side
        
        if ((ketGtoPairsBlock.getBraAngularMomentum() == 0) &&
            (ketGtoPairsBlock.getKetAngularMomentum() == 0)) return;
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoPairsBlock.getStartPositions();
        
        auto epos = braGtoPairsBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto rqx = ketGtoPairsBlock.getCoordinatesPX();
        
        auto rqy = ketGtoPairsBlock.getCoordinatesPY();
        
        auto rqz = ketGtoPairsBlock.getCoordinatesPZ();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
        {
            // set up pointers to coordinates of W
            
            auto wx = wCoordinates.data(3 * idx);
            
            auto wy = wCoordinates.data(3 * idx + 1);
            
            auto wz = wCoordinates.data(3 * idx + 2);
            
            // set up pointers to distances R(WQ)
            
            auto wqx = wqDistances.data(3 * idx);
            
            auto wqy = wqDistances.data(3 * idx + 1);
            
            auto wqz = wqDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(wx, wy, wz, wqx, wqy, wqz, rqx, rqy,\
                                     rqz: VLX_ALIGN)
            for (int32_t j = 0; j < nKetPrimPairs; j++)
            {
                wqx[j] = wx[j] - rqx[j];
                
                wqy[j] = wy[j] - rqy[j];
                
                wqz[j] = wz[j] - rqz[j];
            }
            
            idx++;
        }
    }
    
    void
    compEffectiveDistancesPQ(      CMemBlock<double>& pqDistances,
                             const CGtoPairsBlock&    braGtoPairsBlock,
                             const CGtoPairsBlock&    ketGtoPairsBlock,
                             const bool               isBraEqualKet,
                             const int32_t            iContrPair)
    {
        // set up dimensions on ket side
        
        auto kdim = ketGtoPairsBlock.getNumberOfScreenedContrPairs();
        
        if (isBraEqualKet) kdim = iContrPair + 1;
        
        // set up pointer to effective P coordinates on ket side
        
        auto krpx = ketGtoPairsBlock.getEffectiveCoordinatesPX();
        
        auto krpy = ketGtoPairsBlock.getEffectiveCoordinatesPY();
        
        auto krpz = ketGtoPairsBlock.getEffectiveCoordinatesPZ();
        
        // set up selected effective P coordinates on bra side
        
        auto bpx = (braGtoPairsBlock.getEffectiveCoordinatesPX())[iContrPair];
        
        auto bpy = (braGtoPairsBlock.getEffectiveCoordinatesPY())[iContrPair];
        
        auto bpz = (braGtoPairsBlock.getEffectiveCoordinatesPZ())[iContrPair];
        
        // set up pointer to effecttive PQ distances
        
        auto rpq = pqDistances.data();
        
        #pragma omp simd aligned(rpq, krpx, krpy, krpz)
        for (int32_t i = 0; i < kdim; i++)
        {
            double pqx = bpx - krpx[i];
            
            double pqy = bpy - krpy[i];
            
            double pqz = bpz - krpz[i];
            
            rpq[i] = std::sqrt(pqx * pqx + pqy * pqy + pqz * pqz); 
        }
    }
    
} // twointsfunc namespace
