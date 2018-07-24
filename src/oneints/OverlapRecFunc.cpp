//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OverlapRecFunc.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace
    
    void
    compOverlapForSS(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& abDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto bnorm = braGtoBlock.getNormFactors();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto knorm = ketGtoBlock.getNormFactors();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // fetch up pi values
        
        auto fpi = mathconst::getPiValue();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(2 * idx);
            
            auto fz = osFactors.data(2 * idx + 1);
            
            auto fb = bnorm[i];
            
            // set up primitives buffer data
            
            auto fovl = primBuffer.data(idx);
            
            #pragma omp simd aligned(fovl, fx, fz, knorm, abx, aby,\
                                     abz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fovl[j] = fb * knorm[j] * std::pow(fpi * fx[j], 1.5)
                
                        * std::exp(-fz[j] * (abx[j] * abx[j] + aby[j] * aby[j] +
                                             
                                             abz[j] * abz[j]));
            }
            
            idx++;
        }
    }
    
    void
    compOverlapForSP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // compute number of primitives of bra side
        
        auto bdim = epos[iContrGto] - spos[iContrGto];
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to distances R(PB)
            
            auto pbx = pbDistances.data(3 * idx);
            
            auto pby = pbDistances.data(3 * idx + 1);
            
            auto pbz = pbDistances.data(3 * idx + 2);
            
            // set up pointers to (S|S) integrals
            
            auto fovl = primBuffer.data(idx);
            
            // set up pointers to (S|P) integrals
            
            auto s_0_x = primBuffer.data(bdim + 3 * idx);
            
            auto s_0_y = primBuffer.data(bdim + 3 * idx + 1);
            
            auto s_0_z = primBuffer.data(bdim + 3 * idx + 2);
            
            #pragma omp simd aligned(pbx, pby, pbz, fovl, s_0_x, s_0_y,\
                                     s_0_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = fovl[j];
                
                s_0_x[j] = pbx[j] * fact;
                
                s_0_y[j] = pby[j] * fact;
                
                s_0_z[j] = pbz[j] * fact;
            }
            
            idx++;
        }
    }
    
    void
    compOverlapForSD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // compute number of primitives of bra side
        
        auto bdim = epos[iContrGto] - spos[iContrGto];
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(2 * idx);
            
            // set up pointers to distances R(PB)
            
            auto pbx = pbDistances.data(3 * idx);
            
            auto pby = pbDistances.data(3 * idx + 1);
            
            auto pbz = pbDistances.data(3 * idx + 2);
            
            // set up pointers to (S|S) integrals
            
            auto s_0_0 = primBuffer.data(idx);
            
            // set up pointers to (S|P) integrals
            
            auto s_0_x = primBuffer.data(bdim + 3 * idx);
            
            auto s_0_y = primBuffer.data(bdim + 3 * idx + 1);
            
            auto s_0_z = primBuffer.data(bdim + 3 * idx + 2);
            
            // set up pointers to (S|D) integrals
            
            int32_t soff = 4 * bdim;
            
            auto s_0_xx = primBuffer.data(soff + 6 * idx);
            
            auto s_0_xy = primBuffer.data(soff + 6 * idx + 1);
            
            auto s_0_xz = primBuffer.data(soff + 6 * idx + 2);
            
            auto s_0_yy = primBuffer.data(soff + 6 * idx + 3);
            
            auto s_0_yz = primBuffer.data(soff + 6 * idx + 4);
            
            auto s_0_zz = primBuffer.data(soff + 6 * idx + 5);
            
            #pragma omp simd aligned(fx, pbx, pby, pbz, s_0_0, s_0_x, s_0_y,\
                                     s_0_z, s_0_xx, s_0_xy, s_0_xz, s_0_yy,\
                                     s_0_yz,s_0_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor
                
                double f2t = 0.50 * fx[j];
                
                // leading x component
                
                double fr = pbx[j];
                
                s_0_xx[j] = fr * s_0_x[j] + f2t * s_0_0[j];
                
                s_0_xy[j] = fr * s_0_y[j];
                
                s_0_xz[j] = fr * s_0_z[j];
                
                // leading y component
                
                fr = pby[j];
                
                s_0_yy[j] = fr * s_0_y[j] + f2t * s_0_0[j];
                
                s_0_yz[j] = fr * s_0_z[j];
                
                // leading z component
                
                s_0_zz[j] = pbz[j] * s_0_z[j] + f2t * s_0_0[j];
            }
            
            idx++;
        }
    }

    void
    compOverlapForSF(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute number of primitives of bra side

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|P) integrals

            auto s_0_x = primBuffer.data(bdim + 3 * idx);

            auto s_0_y = primBuffer.data(bdim + 3 * idx + 1);

            auto s_0_z = primBuffer.data(bdim + 3 * idx + 2);

            // set up pointers to (S|D) integrals

            int32_t t1off = 4 * bdim;

            auto s_0_xx = primBuffer.data(t1off + 6 * idx);

            auto s_0_xy = primBuffer.data(t1off + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(t1off + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(t1off + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(t1off + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(t1off + 6 * idx + 5);

            // set up pointers to (S|F) integrals

            int32_t soff = 10 * bdim;

            auto s_0_xxx = primBuffer.data(soff + 10 * idx);

            auto s_0_xxy = primBuffer.data(soff + 10 * idx + 1);

            auto s_0_xxz = primBuffer.data(soff + 10 * idx + 2);

            auto s_0_xyy = primBuffer.data(soff + 10 * idx + 3);

            auto s_0_xyz = primBuffer.data(soff + 10 * idx + 4);

            auto s_0_xzz = primBuffer.data(soff + 10 * idx + 5);

            auto s_0_yyy = primBuffer.data(soff + 10 * idx + 6);

            auto s_0_yyz = primBuffer.data(soff + 10 * idx + 7);

            auto s_0_yzz = primBuffer.data(soff + 10 * idx + 8);

            auto s_0_zzz = primBuffer.data(soff + 10 * idx + 9);

            #pragma omp simd aligned(fx, pbx, pby, pbz, s_0_x, s_0_y, s_0_z,\
                                     s_0_xx, s_0_xy, s_0_xz, s_0_yy, s_0_yz,\
                                     s_0_zz, s_0_xxx, s_0_xxy, s_0_xxz,\
                                     s_0_xyy, s_0_xyz, s_0_xzz, s_0_yyy,\
                                     s_0_yyz, s_0_yzz, s_0_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                // leading x component

                double fr = pbx[j];

                s_0_xxx[j] = fr * s_0_xx[j] + 2.0 * f2t * s_0_x[j];

                s_0_xxy[j] = fr * s_0_xy[j] + f2t * s_0_y[j];

                s_0_xxz[j] = fr * s_0_xz[j] + f2t * s_0_z[j];

                s_0_xyy[j] = fr * s_0_yy[j];

                s_0_xyz[j] = fr * s_0_yz[j];

                s_0_xzz[j] = fr * s_0_zz[j];

                // leading y component

                fr = pby[j];

                s_0_yyy[j] = fr * s_0_yy[j] + 2.0 * f2t * s_0_y[j];

                s_0_yyz[j] = fr * s_0_yz[j] + f2t * s_0_z[j];

                s_0_yzz[j] = fr * s_0_zz[j];

                // leading z component

                s_0_zzz[j] = pbz[j] * s_0_zz[j] + 2.0 * f2t * s_0_z[j];
            }

            idx++;
        }
    }

    void
    compOverlapForSG(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // compute number of primitives of bra side

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to distances R(PB)

            auto pbx = pbDistances.data(3 * idx);

            auto pby = pbDistances.data(3 * idx + 1);

            auto pbz = pbDistances.data(3 * idx + 2);

            // set up pointers to (S|D) integrals

            int32_t t2off = 4 * bdim;

            auto s_0_xx = primBuffer.data(t2off + 6 * idx);

            auto s_0_xy = primBuffer.data(t2off + 6 * idx + 1);

            auto s_0_xz = primBuffer.data(t2off + 6 * idx + 2);

            auto s_0_yy = primBuffer.data(t2off + 6 * idx + 3);

            auto s_0_yz = primBuffer.data(t2off + 6 * idx + 4);

            auto s_0_zz = primBuffer.data(t2off + 6 * idx + 5);

            // set up pointers to (S|F) integrals

            int32_t t1off = 10 * bdim;

            auto s_0_xxx = primBuffer.data(t1off + 10 * idx);

            auto s_0_xxy = primBuffer.data(t1off + 10 * idx + 1);

            auto s_0_xxz = primBuffer.data(t1off + 10 * idx + 2);

            auto s_0_xyy = primBuffer.data(t1off + 10 * idx + 3);

            auto s_0_xyz = primBuffer.data(t1off + 10 * idx + 4);

            auto s_0_xzz = primBuffer.data(t1off + 10 * idx + 5);

            auto s_0_yyy = primBuffer.data(t1off + 10 * idx + 6);

            auto s_0_yyz = primBuffer.data(t1off + 10 * idx + 7);

            auto s_0_yzz = primBuffer.data(t1off + 10 * idx + 8);

            auto s_0_zzz = primBuffer.data(t1off + 10 * idx + 9);

            // set up pointers to (S|G) integrals

            int32_t soff = 20 * bdim;

            auto s_0_xxxx = primBuffer.data(soff + 15 * idx);

            auto s_0_xxxy = primBuffer.data(soff + 15 * idx + 1);

            auto s_0_xxxz = primBuffer.data(soff + 15 * idx + 2);

            auto s_0_xxyy = primBuffer.data(soff + 15 * idx + 3);

            auto s_0_xxyz = primBuffer.data(soff + 15 * idx + 4);

            auto s_0_xxzz = primBuffer.data(soff + 15 * idx + 5);

            auto s_0_xyyy = primBuffer.data(soff + 15 * idx + 6);

            auto s_0_xyyz = primBuffer.data(soff + 15 * idx + 7);

            auto s_0_xyzz = primBuffer.data(soff + 15 * idx + 8);

            auto s_0_xzzz = primBuffer.data(soff + 15 * idx + 9);

            auto s_0_yyyy = primBuffer.data(soff + 15 * idx + 10);

            auto s_0_yyyz = primBuffer.data(soff + 15 * idx + 11);

            auto s_0_yyzz = primBuffer.data(soff + 15 * idx + 12);

            auto s_0_yzzz = primBuffer.data(soff + 15 * idx + 13);

            auto s_0_zzzz = primBuffer.data(soff + 15 * idx + 14);

            #pragma omp simd aligned(fx, pbx, pby, pbz, s_0_xx, s_0_xy, s_0_xz,\
                                     s_0_yy, s_0_yz, s_0_zz, s_0_xxx, s_0_xxy,\
                                     s_0_xxz, s_0_xyy, s_0_xyz, s_0_xzz, s_0_yyy,\
                                     s_0_yyz, s_0_yzz, s_0_zzz, s_0_xxxx,\
                                     s_0_xxxy, s_0_xxxz, s_0_xxyy, s_0_xxyz,\
                                     s_0_xxzz, s_0_xyyy, s_0_xyyz, s_0_xyzz,\
                                     s_0_xzzz, s_0_yyyy, s_0_yyyz, s_0_yyzz,\
                                     s_0_yzzz, s_0_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                // scaled prefactor

                double f2t = 0.50 * fx[j];

                // leading x component

                double fr = pbx[j];

                s_0_xxxx[j] = fr * s_0_xxx[j] + 3.0 * f2t * s_0_xx[j];

                s_0_xxxy[j] = fr * s_0_xxy[j] + 2.0 * f2t * s_0_xy[j];

                s_0_xxxz[j] = fr * s_0_xxz[j] + 2.0 * f2t * s_0_xz[j];

                s_0_xxyy[j] = fr * s_0_xyy[j] + f2t * s_0_yy[j];

                s_0_xxyz[j] = fr * s_0_xyz[j] + f2t * s_0_yz[j];

                s_0_xxzz[j] = fr * s_0_xzz[j] + f2t * s_0_zz[j];

                s_0_xyyy[j] = fr * s_0_yyy[j];

                s_0_xyyz[j] = fr * s_0_yyz[j];

                s_0_xyzz[j] = fr * s_0_yzz[j];

                s_0_xzzz[j] = fr * s_0_zzz[j];

                // leading y component

                fr = pby[j];

                s_0_yyyy[j] = fr * s_0_yyy[j] + 3.0 * f2t * s_0_yy[j];

                s_0_yyyz[j] = fr * s_0_yyz[j] + 2.0 * f2t * s_0_yz[j];

                s_0_yyzz[j] = fr * s_0_yzz[j] + f2t * s_0_zz[j];

                s_0_yzzz[j] = fr * s_0_zzz[j];

                // leading z component

                s_0_zzzz[j] = pbz[j] * s_0_zzz[j] + 3.0 * f2t * s_0_zz[j];
            }

            idx++;
        }
    }
    
    void
    compOverlapForPS(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& paDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // compute number of primitives of bra side
        
        auto bdim = epos[iContrGto] - spos[iContrGto];
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to distances R(PA)
            
            auto pax = paDistances.data(3 * idx);
            
            auto pay = paDistances.data(3 * idx + 1);
            
            auto paz = paDistances.data(3 * idx + 2);
            
            // set up pointers to (S|S) integrals
            
            auto fovl = primBuffer.data(idx);
            
            // set up pointers to (P|S) integrals
            
            auto s_x_0 = primBuffer.data(bdim + 3 * idx);
            
            auto s_y_0 = primBuffer.data(bdim + 3 * idx + 1);
            
            auto s_z_0 = primBuffer.data(bdim + 3 * idx + 2);
            
            #pragma omp simd aligned(pax, pay, paz, fovl, s_x_0, s_y_0,\
                                     s_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fact = fovl[j];
                
                s_x_0[j] = pax[j] * fact;
                
                s_y_0[j] = pay[j] * fact;
                
                s_z_0[j] = paz[j] * fact;
            }
            
            idx++;
        }
    }
    
   
    
} // ovlrecfunc namespace
