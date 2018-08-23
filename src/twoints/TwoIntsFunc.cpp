//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

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
        
        // set up pointers to primitives data on ket side
        
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
    
    
} // twointsfunc namespace
