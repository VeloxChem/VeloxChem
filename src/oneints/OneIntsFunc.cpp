//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OneIntsFunc.hpp"

#include "MathFunc.hpp"

namespace intsfunc { // intsfunc namespace
    
    void
    compDistancesAB(      CMemBlock2D<double>& abDistances,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
                    const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto brx = braGtoBlock.getCoordinatesX();
        
        auto bry = braGtoBlock.getCoordinatesY();
        
        auto brz = braGtoBlock.getCoordinatesZ();
    
        auto spos = braGtoBlock.getStartPositions();
        
        // compute distances
        
        mathfunc::distances(abDistances.data(0),
                            abDistances.data(1),
                            abDistances.data(2),
                            brx[spos[iContrGto]],
                            bry[spos[iContrGto]],
                            brz[spos[iContrGto]],
                            ketGtoBlock.getCoordinatesX(),
                            ketGtoBlock.getCoordinatesY(),
                            ketGtoBlock.getCoordinatesZ(),
                            ketGtoBlock.getNumberOfPrimGtos());
    }

    void
    compFactorsForOverlap(      CMemBlock2D<double>& osFactors,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto bexp = braGtoBlock.getExponents();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto kexp = ketGtoBlock.getExponents();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(2 * idx);
            
            auto fz = osFactors.data(2 * idx + 1);
            
            auto fb = bexp[i];
            
            #pragma omp simd aligned(fx, fz, kexp: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fx[i] = 1.0 / (fb + kexp[j]);
                
                fz[i] = fb * kexp[i] * fx[i];
            }
            
            idx++;
        }
    }
    
    void
    compDistancesPA(      CMemBlock2D<double>& paDistances,
                    const CMemBlock2D<double>& abDistances,
                    const CMemBlock2D<double>& osFactors,
                    const int32_t              nFactors,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
                    const int32_t              iContrGto)
    {
        // skip computation for zero angular momentum on bra side
        
        if (braGtoBlock.getAngularMomentum() == 0) return;
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto kexp = ketGtoBlock.getExponents();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to prefactor
            
            auto fx = osFactors.data(nFactors * idx);
           
            // set up pointers to distances
            
            auto pax = paDistances.data(3 * idx);
            
            auto pay = paDistances.data(3 * idx + 1);
            
            auto paz = paDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(kexp, abx, aby, abz, fx, pax, pay, \
                                     paz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = -kexp[i] * fx[i];
                
                pax[i] = fact * abx[i];
                
                pay[i] = fact * aby[i];
                
                paz[i] = fact * abz[i];
            }
            
            idx++;
        }
    }
    
    void
    compDistancesPB(      CMemBlock2D<double>& pbDistances,
                    const CMemBlock2D<double>& abDistances,
                    const CMemBlock2D<double>& osFactors,
                    const int32_t              nFactors,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
                    const int32_t              iContrGto)
    {
        // skip computation for zero angular momentum on ket side
        
        if (ketGtoBlock.getAngularMomentum() == 0) return;
        
        // set up pointers to primitives data on bra side
        
        auto bexp = braGtoBlock.getExponents();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to prefactors
            
            auto fx = osFactors.data(nFactors * idx);
            
            auto fb = bexp[i];
            
            // set up pointers to distances
            
            auto pbx = pbDistances.data(3 * idx);
            
            auto pby = pbDistances.data(3 * idx + 1);
            
            auto pbz = pbDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(abx, aby, abz, fx, pbx, pby, pbz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = fb * fx[i];
                
                pbx[i] = fact * abx[i];
                
                pby[i] = fact * aby[i];
                
                pbz[i] = fact * abz[i];
            }
            
            idx++;
        }
    }
    
    void
    compXiAndZeta(      double* factorsXi,
                        double* factorsZeta,
                  const double  braExponent,
                  const double* ketExponents,
                  const int32_t nElements)
    {
        #pragma omp simd aligned(factorsXi, factorsZeta, ketExponents: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++)
        {
            factorsXi[i] = 1.0 / (braExponent + ketExponents[i]);
            
            factorsZeta[i] = braExponent * ketExponents[i] * factorsXi[i];
        }
    }
    
    void
    compDistancesPA(      double* xDistancesPA,
                          double* yDistancesPA,
                          double* zDistancesPA,
                    const double* ketExponents,
                    const double* factorsXi,
                    const double* xDistancesAB,
                    const double* yDistancesAB,
                    const double* zDistancesAB,
                    const int32_t nElements)
    {
        #pragma omp simd aligned(xDistancesPA, yDistancesPA, zDistancesPA,\
                                 ketExponents, factorsXi, xDistancesAB, \
                                 yDistancesAB, zDistancesAB: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++)
        {
            double fact = - ketExponents[i] * factorsXi[i];
            
            xDistancesPA[i] = fact * xDistancesAB[i];
            
            yDistancesPA[i] = fact * yDistancesAB[i];
            
            zDistancesPA[i] = fact * zDistancesAB[i];
        }
    }
    
    void
    compDistancesPB(      double* xDistancesPB,
                          double* yDistancesPB,
                          double* zDistancesPB,
                    const double  braExponent,
                    const double* factorsXi,
                    const double* xDistancesAB,
                    const double* yDistancesAB,
                    const double* zDistancesAB,
                    const int32_t nElements)
    {
       #pragma omp simd aligned(xDistancesPB, yDistancesPB, zDistancesPB,\
                                factorsXi, xDistancesAB, yDistancesAB,\
                                zDistancesAB: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++)
        {
            double fact = braExponent * factorsXi[i];
            
            xDistancesPB[i] = fact * xDistancesAB[i];
            
            yDistancesPB[i] = fact * yDistancesAB[i];
            
            zDistancesPB[i] = fact * zDistancesAB[i];
        }
    }
    
} // intsfunc namespace
