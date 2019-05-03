//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OneIntsFunc.hpp"

#include "MathFunc.hpp"
#include "AngularMomentum.hpp"
#include "GenFunc.hpp"

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
                fx[j] = 1.0 / (fb + kexp[j]);
                
                fz[j] = fb * kexp[j] * fx[j];
            }
            
            idx++;
        }
    }
    
    void
    compFactorsForKineticEnergy(      CMemBlock2D<double>& osFactors,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
    {
        // set up angular momentum for bra and ket sides
        
        auto bang = braGtoBlock.getAngularMomentum();
        
        auto kang = ketGtoBlock.getAngularMomentum();
        
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
            
            auto fx = osFactors.data(4 * idx);
            
            auto fz = osFactors.data(4 * idx + 1);
            
            auto fb = bexp[i];
            
            #pragma omp simd aligned(fx, fz, kexp: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fx[j] = 1.0 / (fb + kexp[j]);
                
                fz[j] = fb * kexp[j] * fx[j];
            }
            
            if (bang > 1)
            {
                auto ga = osFactors.data(4 * idx + 2);
                
                auto fga = 1.0 / fb;
                
                #pragma omp simd aligned(ga: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++) ga[j] = fga;
            }
            
            if (kang > 1)
            {
                auto gb = osFactors.data(4 * idx + 3);
             
                if (idx == 0)
                {
                    #pragma omp simd aligned(gb, kexp: VLX_ALIGN)
                    for (int32_t j = 0; j < nprim; j++) gb[j] = 1.0 / kexp[j];
                }
                else
                {
                    auto gb0 = osFactors.data(3);
                    
                    #pragma omp simd aligned(gb, gb0: VLX_ALIGN)
                    for (int32_t j = 0; j < nprim; j++) gb[j] = gb0[j];
                }
            }
            
            idx++;
        }
    }
    
    void
    compFactorsForNuclearPotential(      CMemBlock2D<double>& osFactors,
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
            
            auto fx = osFactors.data(3 * idx);
            
            auto fz = osFactors.data(3 * idx + 1);
            
            auto fg = osFactors.data(3 * idx + 2);
            
            auto fb = bexp[i];
            
            #pragma omp simd aligned(fx, fz, fg, kexp: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = fb + kexp[j];
                
                fg[j] = fact;
                
                fx[j] = 1.0 / fact;
                
                fz[j] = fb * kexp[j] * fx[j];
            }
            
            idx++;
        }
    }
    
    void
    compFactorsForElectronicPotential(      CMemBlock2D<double>& osFactors,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto)
    {
        // set up angular momentum for bra and ket sides
        
        auto bang = braGtoBlock.getAngularMomentum();
        
        auto kang = ketGtoBlock.getAngularMomentum();
        
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
            
            auto fx = osFactors.data(6 * idx);
            
            auto fz = osFactors.data(6 * idx + 1);
            
            auto fia = osFactors.data(6 * idx + 2);
            
            auto fib = osFactors.data(6 * idx + 4);
            
            auto fb = bexp[i];
            
            auto fact = 1.0 / fb;
            
            #pragma omp simd aligned(fx, fz, fia, fib, kexp: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fx[j] = 1.0 / (fb + kexp[j]);
                
                fz[j] = fb * kexp[j] * fx[j];
                
                fia[j] = fact;
                
                fib[j] = 1.0 / kexp[j];
            }
            
           if (bang > 1)
           {
                auto fra = osFactors.data(6 * idx + 3);
                
                #pragma omp simd aligned(fz, fra: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    fra[j] = fact * fz[j];
                }
            }
            
            if (kang > 1)
            {
                auto frb = osFactors.data(6 * idx + 5);
                
                if (idx == 0)
                {
                    #pragma omp simd aligned(fz, fib, frb: VLX_ALIGN)
                    for (int32_t j = 0; j < nprim; j++)
                    {
                        frb[j] = fib[j] * fz[j];
                    }
                }
                else
                {
                    #pragma omp simd aligned(fz, fib, frb: VLX_ALIGN)
                    for (int32_t j = 0; j < nprim; j++)
                    {
                        frb[j] = fib[j] * fz[j];
                    }
                }
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
           
            // set up pointers to distances R(PA)
            
            auto pax = paDistances.data(3 * idx);
            
            auto pay = paDistances.data(3 * idx + 1);
            
            auto paz = paDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(kexp, abx, aby, abz, fx, pax, pay, \
                                     paz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = -kexp[j] * fx[j];
                
                pax[j] = fact * abx[j];
                
                pay[j] = fact * aby[j];
                
                paz[j] = fact * abz[j];
            }
            
            idx++;
        }
    }
    
    void
    compTensorsPA(      CMemBlock2D<double>& paDistances,
                  const CMemBlock2D<double>& abDistances,
                  const CMemBlock2D<double>& osFactors,
                  const int32_t              nFactors,
                  const CGtoBlock&           braGtoBlock,
                  const CGtoBlock&           ketGtoBlock,
                  const int32_t              iContrGto)
    {
        // set up angular momentum of bra side
        
        auto bang = braGtoBlock.getAngularMomentum();
        
        if (bang > 0)
        {
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
            
            // determine number of tensor components
            
            auto tcomps = intsfunc::getNumberOfComponentsInDistancesTensor(bang);
        
            // loop over contracted GTO on bra side
        
            int32_t idx = 0;
        
            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to prefactor
            
                auto fx = osFactors.data(nFactors * idx);
            
                // set up pointers to distances R(PA)
            
                auto poff = tcomps * idx;
            
                auto pax = paDistances.data(poff);
            
                auto pay = paDistances.data(poff + 1);
            
                auto paz = paDistances.data(poff + 2);
            
                // compute R(PA) distances
            
                #pragma omp simd aligned(kexp, abx, aby, abz, fx, pax, pay, \
                                         paz: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fact = -kexp[j] * fx[j];
                
                    pax[j] = fact * abx[j];
                
                    pay[j] = fact * aby[j];
                
                    paz[j] = fact * abz[j];
                }
            
                // compute higher order reduced R(PA) tensors
            
                if (bang > 1)
                {
                    genfunc::compTensorTwoFromVector(paDistances, poff, poff + 3);
                
                    for (int32_t j = 3; j <= bang; j++)
                    {
                        auto t10ff = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 2);
                        
                        auto t2off = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 1);
                        
                        genfunc::compTensorFromVectorAndTensor(paDistances, poff,
                                                               t10ff,  t2off, j);
                    }
                }
            
                idx++;
            }
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
            
            // set up pointers to distances R(PB)
            
            auto pbx = pbDistances.data(3 * idx);
            
            auto pby = pbDistances.data(3 * idx + 1);
            
            auto pbz = pbDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(abx, aby, abz, fx, pbx, pby,\
                                     pbz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = fb * fx[j];
                
                pbx[j] = fact * abx[j];
                
                pby[j] = fact * aby[j];
                
                pbz[j] = fact * abz[j];
            }
            
            idx++;
        }
    }
    
    void
    compTensorsPB(      CMemBlock2D<double>& pbDistances,
                  const CMemBlock2D<double>& abDistances,
                  const CMemBlock2D<double>& osFactors,
                  const int32_t              nFactors,
                  const CGtoBlock&           braGtoBlock,
                  const CGtoBlock&           ketGtoBlock,
                  const int32_t              iContrGto)
    {
        // set up angular momentum of ket side
        
        auto kang = ketGtoBlock.getAngularMomentum();
        
        if (kang > 0)
        {
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
            
            // determine number of tensor components
            
            auto tcomps = intsfunc::getNumberOfComponentsInDistancesTensor(kang);
        
            // loop over contracted GTO on bra side
        
            int32_t idx = 0;
        
            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to prefactors
            
                auto fx = osFactors.data(nFactors * idx);
            
                auto fb = bexp[i];
            
                // set up pointers to distances R(PB)
            
                auto poff = tcomps * idx;
                
                auto pbx = pbDistances.data(poff);
            
                auto pby = pbDistances.data(poff + 1);
            
                auto pbz = pbDistances.data(poff + 2);
                
                // compute R(PB) distances
            
                #pragma omp simd aligned(abx, aby, abz, fx, pbx, pby,\
                                         pbz: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fact = fb * fx[j];
                
                    pbx[j] = fact * abx[j];
                
                    pby[j] = fact * aby[j];
                
                    pbz[j] = fact * abz[j];
                }
                
                // compute higher order reduced R(PB) tensors
                
                if (kang > 1)
                {
                    genfunc::compTensorTwoFromVector(pbDistances, poff, poff + 3);
                    
                    for (int32_t j = 3; j <= kang; j++)
                    {
                        auto t10ff = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 2);
                        
                        auto t2off = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 1);
                        
                        genfunc::compTensorFromVectorAndTensor(pbDistances, poff,
                                                               t10ff,  t2off, j);
                    }
                }
            
                idx++;
            }
        }
    }
    
    void
    compCoordinatesForP(      CMemBlock2D<double>& pCoordinates,
                        const CMemBlock2D<double>& osFactors,
                        const int32_t              nFactors,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
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
        
        auto krx = ketGtoBlock.getCoordinatesX();
        
        auto kry = ketGtoBlock.getCoordinatesY();
        
        auto krz = ketGtoBlock.getCoordinatesZ();
        
        auto kexp = ketGtoBlock.getExponents();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
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
            
            auto px = pCoordinates.data(3 * idx);
            
            auto py = pCoordinates.data(3 * idx + 1);
            
            auto pz = pCoordinates.data(3 * idx + 2);
            
            #pragma omp simd aligned(fx, px, py, pz, kexp, krx, kry,\
                                     krz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = fx[j];
                
                px[j] = fact * (fax + kexp[j] * krx[j]);
                
                py[j] = fact * (fay + kexp[j] * kry[j]);
                
                pz[j] = fact * (faz + kexp[j] * krz[j]);
            }
            
            idx++;
        }
    }
    
    void
    compTensorsPA(      CMemBlock2D<double>& paDistances,
                  const CMemBlock2D<double>& pCoordinates,
                  const CGtoBlock&           braGtoBlock,
                  const CGtoBlock&           ketGtoBlock,
                  const int32_t              iContrGto)
    {
        // skip computation for zero angular momentum on bra side
        
        auto bang = braGtoBlock.getAngularMomentum();
        
        if (bang == 0) return;
        
        // set up pointers to primitives data on bra side
        
        auto brx = braGtoBlock.getCoordinatesX();
        
        auto bry = braGtoBlock.getCoordinatesY();
        
        auto brz = braGtoBlock.getCoordinatesZ();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // determine number of tensor components
        
        auto tcomps = intsfunc::getNumberOfComponentsInDistancesTensor(bang);
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of P
            
            auto px = pCoordinates.data(3 * idx);
            
            auto py = pCoordinates.data(3 * idx + 1);
            
            auto pz = pCoordinates.data(3 * idx + 2);
            
            // set up pmitive GTO data on bra side
            
            auto ax = brx[i];
            
            auto ay = bry[i];
            
            auto az = brz[i];
            
            // set up pointers to distances R(PA)
            
            auto poff = tcomps * idx;
            
            auto pax = paDistances.data(poff);
            
            auto pay = paDistances.data(poff + 1);
            
            auto paz = paDistances.data(poff + 2);
            
            #pragma omp simd aligned(px, py, pz, pax, pay, paz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                pax[j] = px[j] - ax;
                
                pay[j] = py[j] - ay;
                
                paz[j] = pz[j] - az;
            }
            
            // compute higher order reduced R(PA) tensors
            
            if (bang > 1)
            {
                genfunc::compTensorTwoFromVector(paDistances, poff, poff + 3);
                
                for (int32_t j = 3; j <= bang; j++)
                {
                    auto t10ff = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 2);
                    
                    auto t2off = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 1);
                    
                    genfunc::compTensorFromVectorAndTensor(paDistances, poff,
                                                           t10ff,  t2off, j);
                }
            }
            
            idx++;
        }
    }
    
    void
    compDistancesPA(      CMemBlock2D<double>& paDistances,
                    const CMemBlock2D<double>& pCoordinates,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
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
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of P
            
            auto px = pCoordinates.data(3 * idx);
            
            auto py = pCoordinates.data(3 * idx + 1);
            
            auto pz = pCoordinates.data(3 * idx + 2);
            
            // set up pmitive GTO data on bra side
            
            auto ax = brx[i];
            
            auto ay = bry[i];
            
            auto az = brz[i];
            
            // set up pointers to distances R(PA)
            
            auto pax = paDistances.data(3 * idx);
            
            auto pay = paDistances.data(3 * idx + 1);
            
            auto paz = paDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(px, py, pz, pax, pay, paz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                pax[j] = px[j] - ax;
                
                pay[j] = py[j] - ay;
                
                paz[j] = pz[j] - az;
            }
            
            idx++;
        }
    }
    
    void compTensorsPB(      CMemBlock2D<double>& pbDistances,
                       const CMemBlock2D<double>& pCoordinates,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
    {
        
        // skip computation for zero angular momentum on ket side
        
        auto kang = ketGtoBlock.getAngularMomentum();
        
        if (kang == 0) return;
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto krx = ketGtoBlock.getCoordinatesX();
        
        auto kry = ketGtoBlock.getCoordinatesY();
        
        auto krz = ketGtoBlock.getCoordinatesZ();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // determine number of tensor components
        
        auto tcomps = intsfunc::getNumberOfComponentsInDistancesTensor(kang);
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of P
            
            auto px = pCoordinates.data(3 * idx);
            
            auto py = pCoordinates.data(3 * idx + 1);
            
            auto pz = pCoordinates.data(3 * idx + 2);
            
            // set up pointers to distances R(PB)
            
            auto poff = tcomps * idx;
            
            auto pbx = pbDistances.data(poff);
            
            auto pby = pbDistances.data(poff + 1);
            
            auto pbz = pbDistances.data(poff + 2);
            
            #pragma omp simd aligned(krx, kry, krz, px, py, pz, pbx, pby,\
                                     pbz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                pbx[j] = px[j] - krx[j];
                
                pby[j] = py[j] - kry[j];
                
                pbz[j] = pz[j] - krz[j];
            }
            
            // compute higher order reduced R(PB) tensors
            
            if (kang > 1)
            {
                genfunc::compTensorTwoFromVector(pbDistances, poff, poff + 3);
                
                for (int32_t j = 3; j <= kang; j++)
                {
                    auto t10ff = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 2);
                    
                    auto t2off = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 1);
                    
                    genfunc::compTensorFromVectorAndTensor(pbDistances, poff,
                                                           t10ff,  t2off, j);
                }
            }
            
            idx++;
        }
    }
    
    void
    compDistancesPB(      CMemBlock2D<double>& pbDistances,
                    const CMemBlock2D<double>& pCoordinates,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
                    const int32_t              iContrGto)
    {
        // skip computation for zero angular momentum on bra side
        
        if (ketGtoBlock.getAngularMomentum() == 0) return;
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto krx = ketGtoBlock.getCoordinatesX();
        
        auto kry = ketGtoBlock.getCoordinatesY();
        
        auto krz = ketGtoBlock.getCoordinatesZ();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of P
            
            auto px = pCoordinates.data(3 * idx);
            
            auto py = pCoordinates.data(3 * idx + 1);
            
            auto pz = pCoordinates.data(3 * idx + 2);
            
            // set up pointers to distances R(PB)
            
            auto pbx = pbDistances.data(3 * idx);
            
            auto pby = pbDistances.data(3 * idx + 1);
            
            auto pbz = pbDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(krx, kry, krz, px, py, pz, pbx, pby,\
                                     pbz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                pbx[j] = px[j] - krx[j];
                
                pby[j] = py[j] - kry[j];
                
                pbz[j] = pz[j] - krz[j];
            }
            
            idx++;
        }
    }
    
    void
    compTensorsPC(      CMemBlock2D<double>& pcDistances,
                  const CMemBlock2D<double>& pCoordinates,
                  const CMemBlock2D<double>& cCoordinates,
                  const int32_t              orderOfTensor,
                  const CGtoBlock&           braGtoBlock,
                  const CGtoBlock&           ketGtoBlock,
                  const int32_t              iContrGto,
                  const int32_t              iPointCharge)
    {
        // set up coordinates of point charges
        
        double crx = (cCoordinates.data(0))[iPointCharge];
        
        double cry = (cCoordinates.data(1))[iPointCharge];
        
        double crz = (cCoordinates.data(2))[iPointCharge];
        
        // compute PC distances : R(PC) = P - C
        
        intsfunc::compTensorsPC(pcDistances, pCoordinates, crx, cry, crz,
                                orderOfTensor, braGtoBlock, ketGtoBlock,
                                iContrGto); 
    }
    
    void
    compDistancesPC(      CMemBlock2D<double>& pcDistances,
                    const CMemBlock2D<double>& pCoordinates,
                    const CMemBlock2D<double>& cCoordinates,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
                    const int32_t              iContrGto,
                    const int32_t              iPointCharge)
    {
        // set up coordinates of point charges
        
        double crx = (cCoordinates.data(0))[iPointCharge];
        
        double cry = (cCoordinates.data(1))[iPointCharge];
        
        double crz = (cCoordinates.data(2))[iPointCharge];
        
        // compute PC distances : R(PC) = P - C
        
        intsfunc::compDistancesPC(pcDistances, pCoordinates, crx, cry, crz,
                                  braGtoBlock, ketGtoBlock, iContrGto);
    }
    
    void
    compDistancesPC(      CMemBlock2D<double>& pcDistances,
                    const CMemBlock2D<double>& pCoordinates,
                    const double               xCoordinateC,
                    const double               yCoordinateC,
                    const double               zCoordinateC,
                    const CGtoBlock&           braGtoBlock,
                    const CGtoBlock&           ketGtoBlock,
                    const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of P
            
            auto px = pCoordinates.data(3 * idx);
            
            auto py = pCoordinates.data(3 * idx + 1);
            
            auto pz = pCoordinates.data(3 * idx + 2);
            
            // set up pointers to distances R(PC)
            
            auto pcx = pcDistances.data(3 * idx);
            
            auto pcy = pcDistances.data(3 * idx + 1);
            
            auto pcz = pcDistances.data(3 * idx + 2);
            
            #pragma omp simd aligned(px, py, pz, pcx, pcy, pcz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                pcx[j] = px[j] - xCoordinateC;
                
                pcy[j] = py[j] - yCoordinateC;
                
                pcz[j] = pz[j] - zCoordinateC;
            }
            
            idx++;
        }
    }
    
    void
    compTensorsPC(      CMemBlock2D<double>& pcDistances,
                  const CMemBlock2D<double>& pCoordinates,
                  const double               xCoordinateC,
                  const double               yCoordinateC,
                  const double               zCoordinateC,
                  const int32_t              orderOfTensor,
                  const CGtoBlock&           braGtoBlock,
                  const CGtoBlock&           ketGtoBlock,
                  const int32_t              iContrGto)
    {
        if (orderOfTensor == 0) return;
        
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // determine number of tensor components
        
        auto tcomps = intsfunc::getNumberOfComponentsInDistancesTensor(orderOfTensor);
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to coordinates of P
            
            auto px = pCoordinates.data(3 * idx);
            
            auto py = pCoordinates.data(3 * idx + 1);
            
            auto pz = pCoordinates.data(3 * idx + 2);
            
            // set up pointers to distances R(PC)
            
            auto poff = tcomps * idx;
            
            auto pcx = pcDistances.data(poff);
            
            auto pcy = pcDistances.data(poff + 1);
            
            auto pcz = pcDistances.data(poff + 2);
            
            #pragma omp simd aligned(px, py, pz, pcx, pcy, pcz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                pcx[j] = px[j] - xCoordinateC;
                
                pcy[j] = py[j] - yCoordinateC;
                
                pcz[j] = pz[j] - zCoordinateC;
            }
            
            // compute higher order reduced R(PB) tensors
            
            if (orderOfTensor> 1)
            {
                genfunc::compTensorTwoFromVector(pcDistances, poff, poff + 3);
                
                for (int32_t j = 3; j <= orderOfTensor; j++)
                {
                    auto t10ff = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 2);
                    
                    auto t2off = poff + intsfunc::getNumberOfComponentsInDistancesTensor(j - 1);
                    
                    genfunc::compTensorFromVectorAndTensor(pcDistances, poff,
                                                           t10ff,  t2off, j);
                }
            }
            
            idx++;
        }
    }
    
    int32_t
    getNumberOfComponentsInDistancesTensor(const int32_t tensorOrder)
    {
        int32_t ncomps = 0;
        
        for (int32_t i = 1; i <= tensorOrder; i++)
        {
            ncomps += angmom::to_CartesianComponents(i); 
        }
        
        return ncomps;
    }
    
    void
    compTensorsProduct(      CMemBlock2D<double>& aTensor,
                       const CMemBlock2D<double>& bTensor,
                       const CMemBlock2D<double>& cTensor,
                       const int32_t              bDimensions,
                       const int32_t              cDimensions,
                       const int32_t              nBlocks)
    {
        // set up dimensions of tensors
        
        auto adim = bDimensions * cDimensions;
        
        // compute product of tensors
      
        auto ndim = aTensor.size(0);
        
        for (int32_t i = 0; i < nBlocks; i++)
        {
            for (int32_t j = 0; j < bDimensions; j++)
            {
                auto bvals = bTensor.data(bDimensions * i + j);
                
                for (int32_t k = 0; k < cDimensions; k++)
                {
                    auto cvals = cTensor.data(cDimensions * i + k);
                    
                    auto avals = aTensor.data(adim * i + j * cDimensions + k);
                    
                    #pragma omp simd aligned(avals, bvals, cvals: VLX_ALIGN)
                    for (int32_t l = 0; l < ndim; l++)
                    {
                        avals[l] = bvals[l] * cvals[l];
                    }
                }
            }
        }
    }
    
} // intsfunc namespace
