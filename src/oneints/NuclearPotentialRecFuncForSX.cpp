//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "NuclearPotentialRecFuncForSX.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace npotrecfunc { // npotrecfunc namespace

    void
    compNuclearPotentialForSS(      CMemBlock2D<double>&  primBuffer,
                                    CMemBlock2D<double>&  auxBuffer,
                              const CBoysFunction&        bfTable,
                                    CMemBlock<double>&    bfArguments,
                                    CMemBlock2D<double>&  bfValues,
                              const int32_t               bfOrder,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  abDistances,
                              const CMemBlock2D<double>&  pcDistances,
                              const int32_t               pcComponents,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto)
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
        
        // flag for (s|T|s) integrals generation
        
        bool doints = ((braGtoBlock.getAngularMomentum() == 0) &&
                       (ketGtoBlock.getAngularMomentum() == 0));
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up Obara-Saika prefactors
            
            auto fg = osFactors.data(3 * idx + 2);
            
            // set up pointers to ditances R(PC)
            
            auto pcx = pcDistances.data(pcComponents * idx);
            
            auto pcy = pcDistances.data(pcComponents * idx + 1);
            
            auto pcz = pcDistances.data(pcComponents * idx + 2);
            
            // compute Boys function argument
            
            auto fargs = bfArguments.data();
            
            #pragma omp simd aligned(fargs, fg, pcx, pcy, pcz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fargs[j] = fg[j] * (pcx[j] * pcx[j] + pcy[j] * pcy[j] +
                                    
                                    pcz[j] * pcz[j]);
            }
            
            // evaluate Boys function values
            
            bfTable.compute(bfValues, bfArguments, bfOrder);
            
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(3 * idx);
            
            auto fz = osFactors.data(3 * idx + 1);
            
            auto fb = bnorm[i];
            
            // fetch up pi values
            
            auto fpi = 2.0 * mathconst::getPiValue();
            
            // compute overlap scaling factor
            
            #pragma omp simd aligned(fx, fz, knorm, abx, aby, abz,\
                                     fargs: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fargs[j] = fb * knorm[j] * fpi * fx[j]
                
                         * std::exp(-fz[j] * (abx[j] * abx[j] + aby[j] * aby[j] +
                                     
                                     abz[j] * abz[j]));
            }
            
            // distribute (s|A(0)|s) integrals
            
            for (int32_t j = 0; j <= bfOrder; j++)
            {
                auto t_0_0 = (doints) ? primBuffer.data(idx) : auxBuffer.data((bfOrder + 1) * idx + j);
                
                auto bvals = bfValues.data(j);
                
                #pragma omp simd aligned(t_0_0, bvals, fargs: VLX_ALIGN)
                for (int32_t k = 0; k < nprim; k++)
                {
                    t_0_0[k] = bvals[k] * fargs[k];
                }
            }
            
            idx++;
        }
    }

    void
    compNuclearPotentialForSP(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        // Batch of Integrals (0,3)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(2 * idx);

            auto s_0_0_1 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_0_x = primBuffer.data(3 * idx);

            auto t_0_y = primBuffer.data(3 * idx + 1);

            auto t_0_z = primBuffer.data(3 * idx + 2);

            // Batch of Integrals (0,3)

            #pragma omp simd aligned(pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, s_0_0_0, s_0_0_1, t_0_x, t_0_y, t_0_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                t_0_x[j] = fl_s_0_0_0 * pb_x[j];

                t_0_x[j] += -fl_s_0_0_1 * pc_x[j];

                t_0_y[j] = fl_s_0_0_0 * pb_y[j];

                t_0_y[j] += -fl_s_0_0_1 * pc_y[j];

                t_0_z[j] = fl_s_0_0_0 * pb_z[j];

                t_0_z[j] += -fl_s_0_0_1 * pc_z[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPS(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        // Batch of Integrals (0,3)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(2 * idx);

            auto s_0_0_1 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_x_0 = primBuffer.data(3 * idx);

            auto t_y_0 = primBuffer.data(3 * idx + 1);

            auto t_z_0 = primBuffer.data(3 * idx + 2);

            // Batch of Integrals (0,3)

            #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, s_0_0_0, s_0_0_1, t_x_0, t_y_0, t_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                t_x_0[j] = fl_s_0_0_0 * pa_x[j];

                t_x_0[j] += -fl_s_0_0_1 * pc_x[j];

                t_y_0[j] = fl_s_0_0_0 * pa_y[j];

                t_y_0[j] += -fl_s_0_0_1 * pc_y[j];

                t_z_0[j] = fl_s_0_0_0 * pa_z[j];

                t_z_0[j] += -fl_s_0_0_1 * pc_z[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForSD(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        // Batch of Integrals (0,6)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(9 * idx);

            auto pc_y = pcDistances.data(9 * idx + 1);

            auto pc_z = pcDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(9 * idx + 3);

            auto pc_xy = pcDistances.data(9 * idx + 4);

            auto pc_xz = pcDistances.data(9 * idx + 5);

            auto pc_yy = pcDistances.data(9 * idx + 6);

            auto pc_yz = pcDistances.data(9 * idx + 7);

            auto pc_zz = pcDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(3 * idx);

            auto s_0_0_1 = auxBuffer.data(3 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(3 * idx + 2);

            // set up pointers to integrals

            auto t_0_xx = primBuffer.data(6 * idx);

            auto t_0_xy = primBuffer.data(6 * idx + 1);

            auto t_0_xz = primBuffer.data(6 * idx + 2);

            auto t_0_yy = primBuffer.data(6 * idx + 3);

            auto t_0_yz = primBuffer.data(6 * idx + 4);

            auto t_0_zz = primBuffer.data(6 * idx + 5);

            // Batch of Integrals (0,6)

            #pragma omp simd aligned(fx, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, pc_x, pc_xx, pc_xy, \
                                     pc_xz, pc_y, pc_yy, pc_yz, pc_z, pc_zz, s_0_0_0, s_0_0_1, s_0_0_2, t_0_xx, t_0_xy, \
                                     t_0_xz, t_0_yy, t_0_yz, t_0_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl1_fx = fx[j];

                t_0_xx[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pb_xx[j]);

                t_0_xx[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - 2.0 * pb_x[j] * pc_x[j]);

                t_0_xx[j] += fl_s_0_0_2 * pc_xx[j];

                t_0_xy[j] = fl_s_0_0_0 * pb_xy[j];

                t_0_xy[j] += fl_s_0_0_1 * (-pb_x[j] * pc_y[j] - pc_x[j] * pb_y[j]);

                t_0_xy[j] += fl_s_0_0_2 * pc_xy[j];

                t_0_xz[j] = fl_s_0_0_0 * pb_xz[j];

                t_0_xz[j] += fl_s_0_0_1 * (-pb_x[j] * pc_z[j] - pc_x[j] * pb_z[j]);

                t_0_xz[j] += fl_s_0_0_2 * pc_xz[j];

                t_0_yy[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pb_yy[j]);

                t_0_yy[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - 2.0 * pb_y[j] * pc_y[j]);

                t_0_yy[j] += fl_s_0_0_2 * pc_yy[j];

                t_0_yz[j] = fl_s_0_0_0 * pb_yz[j];

                t_0_yz[j] += fl_s_0_0_1 * (-pb_y[j] * pc_z[j] - pc_y[j] * pb_z[j]);

                t_0_yz[j] += fl_s_0_0_2 * pc_yz[j];

                t_0_zz[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pb_zz[j]);

                t_0_zz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - 2.0 * pb_z[j] * pc_z[j]);

                t_0_zz[j] += fl_s_0_0_2 * pc_zz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForDS(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        // Batch of Integrals (0,6)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(9 * idx);

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(9 * idx);

            auto pc_y = pcDistances.data(9 * idx + 1);

            auto pc_z = pcDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(9 * idx + 3);

            auto pc_xy = pcDistances.data(9 * idx + 4);

            auto pc_xz = pcDistances.data(9 * idx + 5);

            auto pc_yy = pcDistances.data(9 * idx + 6);

            auto pc_yz = pcDistances.data(9 * idx + 7);

            auto pc_zz = pcDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(3 * idx);

            auto s_0_0_1 = auxBuffer.data(3 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(3 * idx + 2);

            // set up pointers to integrals

            auto t_xx_0 = primBuffer.data(6 * idx);

            auto t_xy_0 = primBuffer.data(6 * idx + 1);

            auto t_xz_0 = primBuffer.data(6 * idx + 2);

            auto t_yy_0 = primBuffer.data(6 * idx + 3);

            auto t_yz_0 = primBuffer.data(6 * idx + 4);

            auto t_zz_0 = primBuffer.data(6 * idx + 5);

            // Batch of Integrals (0,6)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xy, pa_xz, pa_y, pa_yy, pa_yz, pa_z, pa_zz, pc_x, pc_xx, pc_xy, \
                                     pc_xz, pc_y, pc_yy, pc_yz, pc_z, pc_zz, s_0_0_0, s_0_0_1, s_0_0_2, t_xx_0, t_xy_0, \
                                     t_xz_0, t_yy_0, t_yz_0, t_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl1_fx = fx[j];

                t_xx_0[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pa_xx[j]);

                t_xx_0[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - 2.0 * pa_x[j] * pc_x[j]);

                t_xx_0[j] += fl_s_0_0_2 * pc_xx[j];

                t_xy_0[j] = fl_s_0_0_0 * pa_xy[j];

                t_xy_0[j] += fl_s_0_0_1 * (-pa_x[j] * pc_y[j] - pc_x[j] * pa_y[j]);

                t_xy_0[j] += fl_s_0_0_2 * pc_xy[j];

                t_xz_0[j] = fl_s_0_0_0 * pa_xz[j];

                t_xz_0[j] += fl_s_0_0_1 * (-pa_x[j] * pc_z[j] - pc_x[j] * pa_z[j]);

                t_xz_0[j] += fl_s_0_0_2 * pc_xz[j];

                t_yy_0[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pa_yy[j]);

                t_yy_0[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - 2.0 * pa_y[j] * pc_y[j]);

                t_yy_0[j] += fl_s_0_0_2 * pc_yy[j];

                t_yz_0[j] = fl_s_0_0_0 * pa_yz[j];

                t_yz_0[j] += fl_s_0_0_1 * (-pa_y[j] * pc_z[j] - pc_y[j] * pa_z[j]);

                t_yz_0[j] += fl_s_0_0_2 * pc_yz[j];

                t_zz_0[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pa_zz[j]);

                t_zz_0[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - 2.0 * pa_z[j] * pc_z[j]);

                t_zz_0[j] += fl_s_0_0_2 * pc_zz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForSF(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        // Batch of Integrals (0,10)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(19 * idx);

            auto pc_y = pcDistances.data(19 * idx + 1);

            auto pc_z = pcDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(19 * idx + 3);

            auto pc_xy = pcDistances.data(19 * idx + 4);

            auto pc_xz = pcDistances.data(19 * idx + 5);

            auto pc_yy = pcDistances.data(19 * idx + 6);

            auto pc_yz = pcDistances.data(19 * idx + 7);

            auto pc_zz = pcDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(19 * idx + 9);

            auto pc_xxy = pcDistances.data(19 * idx + 10);

            auto pc_xxz = pcDistances.data(19 * idx + 11);

            auto pc_xyy = pcDistances.data(19 * idx + 12);

            auto pc_xyz = pcDistances.data(19 * idx + 13);

            auto pc_xzz = pcDistances.data(19 * idx + 14);

            auto pc_yyy = pcDistances.data(19 * idx + 15);

            auto pc_yyz = pcDistances.data(19 * idx + 16);

            auto pc_yzz = pcDistances.data(19 * idx + 17);

            auto pc_zzz = pcDistances.data(19 * idx + 18);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(4 * idx);

            auto s_0_0_1 = auxBuffer.data(4 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(4 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(4 * idx + 3);

            // set up pointers to integrals

            auto t_0_xxx = primBuffer.data(10 * idx);

            auto t_0_xxy = primBuffer.data(10 * idx + 1);

            auto t_0_xxz = primBuffer.data(10 * idx + 2);

            auto t_0_xyy = primBuffer.data(10 * idx + 3);

            auto t_0_xyz = primBuffer.data(10 * idx + 4);

            auto t_0_xzz = primBuffer.data(10 * idx + 5);

            auto t_0_yyy = primBuffer.data(10 * idx + 6);

            auto t_0_yyz = primBuffer.data(10 * idx + 7);

            auto t_0_yzz = primBuffer.data(10 * idx + 8);

            auto t_0_zzz = primBuffer.data(10 * idx + 9);

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, \
                                     pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, pc_x, pc_xx, pc_xxx, \
                                     pc_xxy, pc_xxz, pc_xy, pc_xyy, pc_xyz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yyy, pc_yyz, pc_yz, \
                                     pc_yzz, pc_z, pc_zz, pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, t_0_xxx, t_0_xxy, \
                                     t_0_xxz, t_0_xyy, t_0_xyz, t_0_xzz, t_0_yyy, t_0_yyz, t_0_yzz, t_0_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl1_fx = fx[j];

                t_0_xxx[j] = fl_s_0_0_0 * (1.5 * pb_x[j] * fl1_fx + pb_xxx[j]);

                t_0_xxx[j] += fl_s_0_0_1 * (-1.5 * pb_x[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx - 3.0 * pb_xx[j] * pc_x[j]);

                t_0_xxx[j] += fl_s_0_0_2 * (1.5 * pc_x[j] * fl1_fx + 3.0 * pb_x[j] * pc_xx[j]);

                t_0_xxx[j] += -fl_s_0_0_3 * pc_xxx[j];

                t_0_xxy[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_y[j] + pb_xxy[j]);

                t_0_xxy[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_y[j] - 0.5 * fl1_fx * pb_y[j] - pb_xx[j] * pc_y[j] - 2.0 * pb_xy[j] * pc_x[j]);

                t_0_xxy[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_y[j] + 2.0 * pb_x[j] * pc_xy[j] + pc_xx[j] * pb_y[j]);

                t_0_xxy[j] += -fl_s_0_0_3 * pc_xxy[j];

                t_0_xxz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_z[j] + pb_xxz[j]);

                t_0_xxz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pb_z[j] - pb_xx[j] * pc_z[j] - 2.0 * pb_xz[j] * pc_x[j]);

                t_0_xxz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + 2.0 * pb_x[j] * pc_xz[j] + pc_xx[j] * pb_z[j]);

                t_0_xxz[j] += -fl_s_0_0_3 * pc_xxz[j];

                t_0_xyy[j] = fl_s_0_0_0 * (0.5 * pb_x[j] * fl1_fx + pb_xyy[j]);

                t_0_xyy[j] += fl_s_0_0_1 * (-0.5 * pb_x[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx - 2.0 * pb_xy[j] * pc_y[j] - pc_x[j] * pb_yy[j]);

                t_0_xyy[j] += fl_s_0_0_2 * (0.5 * pc_x[j] * fl1_fx + pb_x[j] * pc_yy[j] + 2.0 * pc_xy[j] * pb_y[j]);

                t_0_xyy[j] += -fl_s_0_0_3 * pc_xyy[j];

                t_0_xyz[j] = fl_s_0_0_0 * pb_xyz[j];

                t_0_xyz[j] += fl_s_0_0_1 * (-pb_xy[j] * pc_z[j] - pb_xz[j] * pc_y[j] - pc_x[j] * pb_yz[j]);

                t_0_xyz[j] += fl_s_0_0_2 * (pb_x[j] * pc_yz[j] + pc_xz[j] * pb_y[j] + pc_xy[j] * pb_z[j]);

                t_0_xyz[j] += -fl_s_0_0_3 * pc_xyz[j];

                t_0_xzz[j] = fl_s_0_0_0 * (0.5 * pb_x[j] * fl1_fx + pb_xzz[j]);

                t_0_xzz[j] += fl_s_0_0_1 * (-0.5 * pb_x[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx - 2.0 * pb_xz[j] * pc_z[j] - pc_x[j] * pb_zz[j]);

                t_0_xzz[j] += fl_s_0_0_2 * (0.5 * pc_x[j] * fl1_fx + pb_x[j] * pc_zz[j] + 2.0 * pc_xz[j] * pb_z[j]);

                t_0_xzz[j] += -fl_s_0_0_3 * pc_xzz[j];

                t_0_yyy[j] = fl_s_0_0_0 * (1.5 * pb_y[j] * fl1_fx + pb_yyy[j]);

                t_0_yyy[j] += fl_s_0_0_1 * (-1.5 * pb_y[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx - 3.0 * pb_yy[j] * pc_y[j]);

                t_0_yyy[j] += fl_s_0_0_2 * (1.5 * pc_y[j] * fl1_fx + 3.0 * pb_y[j] * pc_yy[j]);

                t_0_yyy[j] += -fl_s_0_0_3 * pc_yyy[j];

                t_0_yyz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_z[j] + pb_yyz[j]);

                t_0_yyz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pb_z[j] - pb_yy[j] * pc_z[j] - 2.0 * pb_yz[j] * pc_y[j]);

                t_0_yyz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + 2.0 * pb_y[j] * pc_yz[j] + pc_yy[j] * pb_z[j]);

                t_0_yyz[j] += -fl_s_0_0_3 * pc_yyz[j];

                t_0_yzz[j] = fl_s_0_0_0 * (0.5 * pb_y[j] * fl1_fx + pb_yzz[j]);

                t_0_yzz[j] += fl_s_0_0_1 * (-0.5 * pb_y[j] * fl1_fx - 0.5 * pc_y[j] * fl1_fx - 2.0 * pb_yz[j] * pc_z[j] - pc_y[j] * pb_zz[j]);

                t_0_yzz[j] += fl_s_0_0_2 * (0.5 * pc_y[j] * fl1_fx + pb_y[j] * pc_zz[j] + 2.0 * pc_yz[j] * pb_z[j]);

                t_0_yzz[j] += -fl_s_0_0_3 * pc_yzz[j];

                t_0_zzz[j] = fl_s_0_0_0 * (1.5 * pb_z[j] * fl1_fx + pb_zzz[j]);

                t_0_zzz[j] += fl_s_0_0_1 * (-1.5 * pb_z[j] * fl1_fx - 1.5 * pc_z[j] * fl1_fx - 3.0 * pb_zz[j] * pc_z[j]);

                t_0_zzz[j] += fl_s_0_0_2 * (1.5 * pc_z[j] * fl1_fx + 3.0 * pb_z[j] * pc_zz[j]);

                t_0_zzz[j] += -fl_s_0_0_3 * pc_zzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForFS(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        // Batch of Integrals (0,10)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(19 * idx);

            auto pc_y = pcDistances.data(19 * idx + 1);

            auto pc_z = pcDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(19 * idx + 3);

            auto pc_xy = pcDistances.data(19 * idx + 4);

            auto pc_xz = pcDistances.data(19 * idx + 5);

            auto pc_yy = pcDistances.data(19 * idx + 6);

            auto pc_yz = pcDistances.data(19 * idx + 7);

            auto pc_zz = pcDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(19 * idx + 9);

            auto pc_xxy = pcDistances.data(19 * idx + 10);

            auto pc_xxz = pcDistances.data(19 * idx + 11);

            auto pc_xyy = pcDistances.data(19 * idx + 12);

            auto pc_xyz = pcDistances.data(19 * idx + 13);

            auto pc_xzz = pcDistances.data(19 * idx + 14);

            auto pc_yyy = pcDistances.data(19 * idx + 15);

            auto pc_yyz = pcDistances.data(19 * idx + 16);

            auto pc_yzz = pcDistances.data(19 * idx + 17);

            auto pc_zzz = pcDistances.data(19 * idx + 18);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(4 * idx);

            auto s_0_0_1 = auxBuffer.data(4 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(4 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(4 * idx + 3);

            // set up pointers to integrals

            auto t_xxx_0 = primBuffer.data(10 * idx);

            auto t_xxy_0 = primBuffer.data(10 * idx + 1);

            auto t_xxz_0 = primBuffer.data(10 * idx + 2);

            auto t_xyy_0 = primBuffer.data(10 * idx + 3);

            auto t_xyz_0 = primBuffer.data(10 * idx + 4);

            auto t_xzz_0 = primBuffer.data(10 * idx + 5);

            auto t_yyy_0 = primBuffer.data(10 * idx + 6);

            auto t_yyz_0 = primBuffer.data(10 * idx + 7);

            auto t_yzz_0 = primBuffer.data(10 * idx + 8);

            auto t_zzz_0 = primBuffer.data(10 * idx + 9);

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxy, pa_xxz, pa_xy, pa_xyy, pa_xyz, pa_xz, pa_xzz, \
                                     pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pa_zzz, pc_x, pc_xx, pc_xxx, \
                                     pc_xxy, pc_xxz, pc_xy, pc_xyy, pc_xyz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yyy, pc_yyz, pc_yz, \
                                     pc_yzz, pc_z, pc_zz, pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, t_xxx_0, t_xxy_0, \
                                     t_xxz_0, t_xyy_0, t_xyz_0, t_xzz_0, t_yyy_0, t_yyz_0, t_yzz_0, t_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl1_fx = fx[j];

                t_xxx_0[j] = fl_s_0_0_0 * (1.5 * pa_x[j] * fl1_fx + pa_xxx[j]);

                t_xxx_0[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx - 3.0 * pa_xx[j] * pc_x[j]);

                t_xxx_0[j] += fl_s_0_0_2 * (1.5 * pc_x[j] * fl1_fx + 3.0 * pa_x[j] * pc_xx[j]);

                t_xxx_0[j] += -fl_s_0_0_3 * pc_xxx[j];

                t_xxy_0[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_y[j] + pa_xxy[j]);

                t_xxy_0[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_y[j] - 0.5 * fl1_fx * pa_y[j] - pa_xx[j] * pc_y[j] - 2.0 * pa_xy[j] * pc_x[j]);

                t_xxy_0[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_y[j] + 2.0 * pa_x[j] * pc_xy[j] + pc_xx[j] * pa_y[j]);

                t_xxy_0[j] += -fl_s_0_0_3 * pc_xxy[j];

                t_xxz_0[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_z[j] + pa_xxz[j]);

                t_xxz_0[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pa_z[j] - pa_xx[j] * pc_z[j] - 2.0 * pa_xz[j] * pc_x[j]);

                t_xxz_0[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + 2.0 * pa_x[j] * pc_xz[j] + pc_xx[j] * pa_z[j]);

                t_xxz_0[j] += -fl_s_0_0_3 * pc_xxz[j];

                t_xyy_0[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx + pa_xyy[j]);

                t_xyy_0[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx - 2.0 * pa_xy[j] * pc_y[j] - pc_x[j] * pa_yy[j]);

                t_xyy_0[j] += fl_s_0_0_2 * (0.5 * pc_x[j] * fl1_fx + pa_x[j] * pc_yy[j] + 2.0 * pc_xy[j] * pa_y[j]);

                t_xyy_0[j] += -fl_s_0_0_3 * pc_xyy[j];

                t_xyz_0[j] = fl_s_0_0_0 * pa_xyz[j];

                t_xyz_0[j] += fl_s_0_0_1 * (-pa_xy[j] * pc_z[j] - pa_xz[j] * pc_y[j] - pc_x[j] * pa_yz[j]);

                t_xyz_0[j] += fl_s_0_0_2 * (pa_x[j] * pc_yz[j] + pc_xz[j] * pa_y[j] + pc_xy[j] * pa_z[j]);

                t_xyz_0[j] += -fl_s_0_0_3 * pc_xyz[j];

                t_xzz_0[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx + pa_xzz[j]);

                t_xzz_0[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx - 2.0 * pa_xz[j] * pc_z[j] - pc_x[j] * pa_zz[j]);

                t_xzz_0[j] += fl_s_0_0_2 * (0.5 * pc_x[j] * fl1_fx + pa_x[j] * pc_zz[j] + 2.0 * pc_xz[j] * pa_z[j]);

                t_xzz_0[j] += -fl_s_0_0_3 * pc_xzz[j];

                t_yyy_0[j] = fl_s_0_0_0 * (1.5 * pa_y[j] * fl1_fx + pa_yyy[j]);

                t_yyy_0[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx - 3.0 * pa_yy[j] * pc_y[j]);

                t_yyy_0[j] += fl_s_0_0_2 * (1.5 * pc_y[j] * fl1_fx + 3.0 * pa_y[j] * pc_yy[j]);

                t_yyy_0[j] += -fl_s_0_0_3 * pc_yyy[j];

                t_yyz_0[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_z[j] + pa_yyz[j]);

                t_yyz_0[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pa_z[j] - pa_yy[j] * pc_z[j] - 2.0 * pa_yz[j] * pc_y[j]);

                t_yyz_0[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + 2.0 * pa_y[j] * pc_yz[j] + pc_yy[j] * pa_z[j]);

                t_yyz_0[j] += -fl_s_0_0_3 * pc_yyz[j];

                t_yzz_0[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl1_fx + pa_yzz[j]);

                t_yzz_0[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl1_fx - 0.5 * pc_y[j] * fl1_fx - 2.0 * pa_yz[j] * pc_z[j] - pc_y[j] * pa_zz[j]);

                t_yzz_0[j] += fl_s_0_0_2 * (0.5 * pc_y[j] * fl1_fx + pa_y[j] * pc_zz[j] + 2.0 * pc_yz[j] * pa_z[j]);

                t_yzz_0[j] += -fl_s_0_0_3 * pc_yzz[j];

                t_zzz_0[j] = fl_s_0_0_0 * (1.5 * pa_z[j] * fl1_fx + pa_zzz[j]);

                t_zzz_0[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * fl1_fx - 1.5 * pc_z[j] * fl1_fx - 3.0 * pa_zz[j] * pc_z[j]);

                t_zzz_0[j] += fl_s_0_0_2 * (1.5 * pc_z[j] * fl1_fx + 3.0 * pa_z[j] * pc_zz[j]);

                t_zzz_0[j] += -fl_s_0_0_3 * pc_zzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForSG(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForSG_0_5(primBuffer, auxBuffer, osFactors, pbDistances, pcDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForSG_5_10(primBuffer, auxBuffer, osFactors, pbDistances, pcDistances, 
                                                    braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForSG_10_15(primBuffer, auxBuffer, osFactors, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compNuclearPotentialForSG_0_5(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (0,5)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(34 * idx);

            auto pc_y = pcDistances.data(34 * idx + 1);

            auto pc_z = pcDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(34 * idx + 3);

            auto pc_xy = pcDistances.data(34 * idx + 4);

            auto pc_xz = pcDistances.data(34 * idx + 5);

            auto pc_yy = pcDistances.data(34 * idx + 6);

            auto pc_yz = pcDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(34 * idx + 9);

            auto pc_xxy = pcDistances.data(34 * idx + 10);

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(34 * idx + 19);

            auto pc_xxxy = pcDistances.data(34 * idx + 20);

            auto pc_xxxz = pcDistances.data(34 * idx + 21);

            auto pc_xxyy = pcDistances.data(34 * idx + 22);

            auto pc_xxyz = pcDistances.data(34 * idx + 23);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

            // set up pointers to integrals

            auto t_0_xxxx = primBuffer.data(15 * idx);

            auto t_0_xxxy = primBuffer.data(15 * idx + 1);

            auto t_0_xxxz = primBuffer.data(15 * idx + 2);

            auto t_0_xxyy = primBuffer.data(15 * idx + 3);

            auto t_0_xxyz = primBuffer.data(15 * idx + 4);

            // Batch of Integrals (0,5)

            #pragma omp simd aligned(fx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pc_x, pc_xx, pc_xxx, \
                                     pc_xxxx, pc_xxxy, pc_xxxz, pc_xxy, pc_xxyy, pc_xxyz, pc_xxz, pc_xy, pc_xyy, pc_xyz, \
                                     pc_xz, pc_y, pc_yy, pc_yz, pc_z, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, \
                                     t_0_xxxx, t_0_xxxy, t_0_xxxz, t_0_xxyy, t_0_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_0_xxxx[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 3.0 * pb_xx[j] * fl1_fx + pb_xxxx[j]);

                t_0_xxxx[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 3.0 * pb_xx[j] * fl1_fx - 6.0 * pb_x[j] * pc_x[j] * fl1_fx - 4.0 * pb_xxx[j] * pc_x[j]);

                t_0_xxxx[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 6.0 * pb_x[j] * pc_x[j] * fl1_fx + 3.0 * pc_xx[j] * fl1_fx + 6.0 * pb_xx[j] * pc_xx[j]);

                t_0_xxxx[j] += fl_s_0_0_3 * (-3.0 * pc_xx[j] * fl1_fx - 4.0 * pb_x[j] * pc_xxx[j]);

                t_0_xxxx[j] += fl_s_0_0_4 * pc_xxxx[j];

                t_0_xxxy[j] = fl_s_0_0_0 * (1.5 * pb_xy[j] * fl1_fx + pb_xxxy[j]);

                t_0_xxxy[j] += fl_s_0_0_1 * (-1.5 * pb_x[j] * fl1_fx * pc_y[j] - 1.5 * pb_xy[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pb_y[j] - pb_xxx[j] * pc_y[j] - 3.0 * pb_xxy[j] * pc_x[j]);

                t_0_xxxy[j] += fl_s_0_0_2 * (1.5 * pb_x[j] * fl1_fx * pc_y[j] + 1.5 * pc_xy[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pb_y[j] + 3.0 * pb_xx[j] * pc_xy[j] + 3.0 * pb_xy[j] * pc_xx[j]);

                t_0_xxxy[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - 3.0 * pb_x[j] * pc_xxy[j] - pc_xxx[j] * pb_y[j]);

                t_0_xxxy[j] += fl_s_0_0_4 * pc_xxxy[j];

                t_0_xxxz[j] = fl_s_0_0_0 * (1.5 * pb_xz[j] * fl1_fx + pb_xxxz[j]);

                t_0_xxxz[j] += fl_s_0_0_1 * (-1.5 * pb_x[j] * fl1_fx * pc_z[j] - 1.5 * pb_xz[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pb_z[j] - pb_xxx[j] * pc_z[j] - 3.0 * pb_xxz[j] * pc_x[j]);

                t_0_xxxz[j] += fl_s_0_0_2 * (1.5 * pb_x[j] * fl1_fx * pc_z[j] + 1.5 * pc_xz[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pb_z[j] + 3.0 * pb_xx[j] * pc_xz[j] + 3.0 * pb_xz[j] * pc_xx[j]);

                t_0_xxxz[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - 3.0 * pb_x[j] * pc_xxz[j] - pc_xxx[j] * pb_z[j]);

                t_0_xxxz[j] += fl_s_0_0_4 * pc_xxxz[j];

                t_0_xxyy[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pb_xx[j] * fl1_fx + 0.5 * fl1_fx * pb_yy[j] + pb_xxyy[j]);

                t_0_xxyy[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pb_xx[j] * fl1_fx - pb_x[j] * pc_x[j] * fl1_fx - fl1_fx * pb_y[j] * pc_y[j] - 0.5 * fl1_fx * pb_yy[j] - 2.0 * pb_xxy[j] * pc_y[j] - 2.0 * pb_xyy[j] * pc_x[j]);

                t_0_xxyy[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pb_x[j] * pc_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx + 0.5 * fl1_fx * pc_yy[j] + fl1_fx * pb_y[j] * pc_y[j] + pb_xx[j] * pc_yy[j] + 4.0 * pb_xy[j] * pc_xy[j] + pc_xx[j] * pb_yy[j]);

                t_0_xxyy[j] += fl_s_0_0_3 * (-0.5 * pc_xx[j] * fl1_fx - 0.5 * fl1_fx * pc_yy[j] - 2.0 * pb_x[j] * pc_xyy[j] - 2.0 * pc_xxy[j] * pb_y[j]);

                t_0_xxyy[j] += fl_s_0_0_4 * pc_xxyy[j];

                t_0_xxyz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_yz[j] + pb_xxyz[j]);

                t_0_xxyz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pb_y[j] * pc_z[j] - 0.5 * fl1_fx * pc_y[j] * pb_z[j] - 0.5 * fl1_fx * pb_yz[j] - pb_xxy[j] * pc_z[j] - pb_xxz[j] * pc_y[j] - 2.0 * pb_xyz[j] * pc_x[j]);

                t_0_xxyz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_yz[j] + 0.5 * fl1_fx * pb_y[j] * pc_z[j] + 0.5 * fl1_fx * pc_y[j] * pb_z[j] + pb_xx[j] * pc_yz[j] + 2.0 * pb_xy[j] * pc_xz[j] + 2.0 * pb_xz[j] * pc_xy[j] + pc_xx[j] * pb_yz[j]);

                t_0_xxyz[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_yz[j] - 2.0 * pb_x[j] * pc_xyz[j] - pc_xxz[j] * pb_y[j] - pc_xxy[j] * pb_z[j]);

                t_0_xxyz[j] += fl_s_0_0_4 * pc_xxyz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForSG_5_10(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (5,10)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(34 * idx);

            auto pc_y = pcDistances.data(34 * idx + 1);

            auto pc_z = pcDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(34 * idx + 3);

            auto pc_xy = pcDistances.data(34 * idx + 4);

            auto pc_xz = pcDistances.data(34 * idx + 5);

            auto pc_yy = pcDistances.data(34 * idx + 6);

            auto pc_yz = pcDistances.data(34 * idx + 7);

            auto pc_zz = pcDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            auto pc_xzz = pcDistances.data(34 * idx + 14);

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxzz = pcDistances.data(34 * idx + 24);

            auto pc_xyyy = pcDistances.data(34 * idx + 25);

            auto pc_xyyz = pcDistances.data(34 * idx + 26);

            auto pc_xyzz = pcDistances.data(34 * idx + 27);

            auto pc_xzzz = pcDistances.data(34 * idx + 28);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

            // set up pointers to integrals

            auto t_0_xxzz = primBuffer.data(15 * idx + 5);

            auto t_0_xyyy = primBuffer.data(15 * idx + 6);

            auto t_0_xyyz = primBuffer.data(15 * idx + 7);

            auto t_0_xyzz = primBuffer.data(15 * idx + 8);

            auto t_0_xzzz = primBuffer.data(15 * idx + 9);

            // Batch of Integrals (5,10)

            #pragma omp simd aligned(fx, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, pc_x, pc_xx, pc_xxz, pc_xxzz, pc_xy, pc_xyy, pc_xyyy, pc_xyyz, pc_xyz, pc_xyzz, \
                                     pc_xz, pc_xzz, pc_xzzz, pc_y, pc_yy, pc_yyy, pc_yyz, pc_yz, pc_yzz, pc_z, pc_zz, pc_zzz, \
                                     s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, t_0_xxzz, t_0_xyyy, t_0_xyyz, t_0_xyzz, \
                                     t_0_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_0_xxzz[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pb_xx[j] * fl1_fx + 0.5 * fl1_fx * pb_zz[j] + pb_xxzz[j]);

                t_0_xxzz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pb_xx[j] * fl1_fx - pb_x[j] * pc_x[j] * fl1_fx - fl1_fx * pb_z[j] * pc_z[j] - 0.5 * fl1_fx * pb_zz[j] - 2.0 * pb_xxz[j] * pc_z[j] - 2.0 * pb_xzz[j] * pc_x[j]);

                t_0_xxzz[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pb_x[j] * pc_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx + 0.5 * fl1_fx * pc_zz[j] + fl1_fx * pb_z[j] * pc_z[j] + pb_xx[j] * pc_zz[j] + 4.0 * pb_xz[j] * pc_xz[j] + pc_xx[j] * pb_zz[j]);

                t_0_xxzz[j] += fl_s_0_0_3 * (-0.5 * pc_xx[j] * fl1_fx - 0.5 * fl1_fx * pc_zz[j] - 2.0 * pb_x[j] * pc_xzz[j] - 2.0 * pc_xxz[j] * pb_z[j]);

                t_0_xxzz[j] += fl_s_0_0_4 * pc_xxzz[j];

                t_0_xyyy[j] = fl_s_0_0_0 * (1.5 * pb_xy[j] * fl1_fx + pb_xyyy[j]);

                t_0_xyyy[j] += fl_s_0_0_1 * (-1.5 * pb_xy[j] * fl1_fx - 1.5 * pb_x[j] * pc_y[j] * fl1_fx - 1.5 * pc_x[j] * pb_y[j] * fl1_fx - 3.0 * pb_xyy[j] * pc_y[j] - pc_x[j] * pb_yyy[j]);

                t_0_xyyy[j] += fl_s_0_0_2 * (1.5 * pb_x[j] * pc_y[j] * fl1_fx + 1.5 * pc_x[j] * pb_y[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx + 3.0 * pb_xy[j] * pc_yy[j] + 3.0 * pc_xy[j] * pb_yy[j]);

                t_0_xyyy[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - pb_x[j] * pc_yyy[j] - 3.0 * pc_xyy[j] * pb_y[j]);

                t_0_xyyy[j] += fl_s_0_0_4 * pc_xyyy[j];

                t_0_xyyz[j] = fl_s_0_0_0 * (0.5 * pb_xz[j] * fl1_fx + pb_xyyz[j]);

                t_0_xyyz[j] += fl_s_0_0_1 * (-0.5 * pb_x[j] * fl1_fx * pc_z[j] - 0.5 * pb_xz[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx * pb_z[j] - pb_xyy[j] * pc_z[j] - 2.0 * pb_xyz[j] * pc_y[j] - pc_x[j] * pb_yyz[j]);

                t_0_xyyz[j] += fl_s_0_0_2 * (0.5 * pb_x[j] * fl1_fx * pc_z[j] + 0.5 * pc_xz[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pb_z[j] + 2.0 * pb_xy[j] * pc_yz[j] + pb_xz[j] * pc_yy[j] + pc_xz[j] * pb_yy[j] + 2.0 * pc_xy[j] * pb_yz[j]);

                t_0_xyyz[j] += fl_s_0_0_3 * (-0.5 * pc_xz[j] * fl1_fx - pb_x[j] * pc_yyz[j] - 2.0 * pc_xyz[j] * pb_y[j] - pc_xyy[j] * pb_z[j]);

                t_0_xyyz[j] += fl_s_0_0_4 * pc_xyyz[j];

                t_0_xyzz[j] = fl_s_0_0_0 * (0.5 * pb_xy[j] * fl1_fx + pb_xyzz[j]);

                t_0_xyzz[j] += fl_s_0_0_1 * (-0.5 * pb_xy[j] * fl1_fx - 0.5 * pb_x[j] * pc_y[j] * fl1_fx - 0.5 * pc_x[j] * pb_y[j] * fl1_fx - 2.0 * pb_xyz[j] * pc_z[j] - pb_xzz[j] * pc_y[j] - pc_x[j] * pb_yzz[j]);

                t_0_xyzz[j] += fl_s_0_0_2 * (0.5 * pb_x[j] * pc_y[j] * fl1_fx + 0.5 * pc_x[j] * pb_y[j] * fl1_fx + 0.5 * pc_xy[j] * fl1_fx + pb_xy[j] * pc_zz[j] + 2.0 * pb_xz[j] * pc_yz[j] + 2.0 * pc_xz[j] * pb_yz[j] + pc_xy[j] * pb_zz[j]);

                t_0_xyzz[j] += fl_s_0_0_3 * (-0.5 * pc_xy[j] * fl1_fx - pb_x[j] * pc_yzz[j] - pc_xzz[j] * pb_y[j] - 2.0 * pc_xyz[j] * pb_z[j]);

                t_0_xyzz[j] += fl_s_0_0_4 * pc_xyzz[j];

                t_0_xzzz[j] = fl_s_0_0_0 * (1.5 * pb_xz[j] * fl1_fx + pb_xzzz[j]);

                t_0_xzzz[j] += fl_s_0_0_1 * (-1.5 * pb_xz[j] * fl1_fx - 1.5 * pb_x[j] * pc_z[j] * fl1_fx - 1.5 * pc_x[j] * pb_z[j] * fl1_fx - 3.0 * pb_xzz[j] * pc_z[j] - pc_x[j] * pb_zzz[j]);

                t_0_xzzz[j] += fl_s_0_0_2 * (1.5 * pb_x[j] * pc_z[j] * fl1_fx + 1.5 * pc_x[j] * pb_z[j] * fl1_fx + 1.5 * pc_xz[j] * fl1_fx + 3.0 * pb_xz[j] * pc_zz[j] + 3.0 * pc_xz[j] * pb_zz[j]);

                t_0_xzzz[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - pb_x[j] * pc_zzz[j] - 3.0 * pc_xzz[j] * pb_z[j]);

                t_0_xzzz[j] += fl_s_0_0_4 * pc_xzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForSG_10_15(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (10,15)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_y = pcDistances.data(34 * idx + 1);

            auto pc_z = pcDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_yy = pcDistances.data(34 * idx + 6);

            auto pc_yz = pcDistances.data(34 * idx + 7);

            auto pc_zz = pcDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_yyyy = pcDistances.data(34 * idx + 29);

            auto pc_yyyz = pcDistances.data(34 * idx + 30);

            auto pc_yyzz = pcDistances.data(34 * idx + 31);

            auto pc_yzzz = pcDistances.data(34 * idx + 32);

            auto pc_zzzz = pcDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

            // set up pointers to integrals

            auto t_0_yyyy = primBuffer.data(15 * idx + 10);

            auto t_0_yyyz = primBuffer.data(15 * idx + 11);

            auto t_0_yyzz = primBuffer.data(15 * idx + 12);

            auto t_0_yzzz = primBuffer.data(15 * idx + 13);

            auto t_0_zzzz = primBuffer.data(15 * idx + 14);

            // Batch of Integrals (10,15)

            #pragma omp simd aligned(fx, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, pc_y, pc_yy, pc_yyy, pc_yyyy, pc_yyyz, pc_yyz, \
                                     pc_yyzz, pc_yz, pc_yzz, pc_yzzz, pc_z, pc_zz, pc_zzz, pc_zzzz, s_0_0_0, s_0_0_1, \
                                     s_0_0_2, s_0_0_3, s_0_0_4, t_0_yyyy, t_0_yyyz, t_0_yyzz, t_0_yzzz, t_0_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_0_yyyy[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 3.0 * pb_yy[j] * fl1_fx + pb_yyyy[j]);

                t_0_yyyy[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 3.0 * pb_yy[j] * fl1_fx - 6.0 * pb_y[j] * pc_y[j] * fl1_fx - 4.0 * pb_yyy[j] * pc_y[j]);

                t_0_yyyy[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 6.0 * pb_y[j] * pc_y[j] * fl1_fx + 3.0 * pc_yy[j] * fl1_fx + 6.0 * pb_yy[j] * pc_yy[j]);

                t_0_yyyy[j] += fl_s_0_0_3 * (-3.0 * pc_yy[j] * fl1_fx - 4.0 * pb_y[j] * pc_yyy[j]);

                t_0_yyyy[j] += fl_s_0_0_4 * pc_yyyy[j];

                t_0_yyyz[j] = fl_s_0_0_0 * (1.5 * pb_yz[j] * fl1_fx + pb_yyyz[j]);

                t_0_yyyz[j] += fl_s_0_0_1 * (-1.5 * pb_y[j] * fl1_fx * pc_z[j] - 1.5 * pb_yz[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx * pb_z[j] - pb_yyy[j] * pc_z[j] - 3.0 * pb_yyz[j] * pc_y[j]);

                t_0_yyyz[j] += fl_s_0_0_2 * (1.5 * pb_y[j] * fl1_fx * pc_z[j] + 1.5 * pc_yz[j] * fl1_fx + 1.5 * pc_y[j] * fl1_fx * pb_z[j] + 3.0 * pb_yy[j] * pc_yz[j] + 3.0 * pb_yz[j] * pc_yy[j]);

                t_0_yyyz[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - 3.0 * pb_y[j] * pc_yyz[j] - pc_yyy[j] * pb_z[j]);

                t_0_yyyz[j] += fl_s_0_0_4 * pc_yyyz[j];

                t_0_yyzz[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pb_yy[j] * fl1_fx + 0.5 * fl1_fx * pb_zz[j] + pb_yyzz[j]);

                t_0_yyzz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pb_yy[j] * fl1_fx - pb_y[j] * pc_y[j] * fl1_fx - fl1_fx * pb_z[j] * pc_z[j] - 0.5 * fl1_fx * pb_zz[j] - 2.0 * pb_yyz[j] * pc_z[j] - 2.0 * pb_yzz[j] * pc_y[j]);

                t_0_yyzz[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pb_y[j] * pc_y[j] * fl1_fx + 0.5 * pc_yy[j] * fl1_fx + 0.5 * fl1_fx * pc_zz[j] + fl1_fx * pb_z[j] * pc_z[j] + pb_yy[j] * pc_zz[j] + 4.0 * pb_yz[j] * pc_yz[j] + pc_yy[j] * pb_zz[j]);

                t_0_yyzz[j] += fl_s_0_0_3 * (-0.5 * pc_yy[j] * fl1_fx - 0.5 * fl1_fx * pc_zz[j] - 2.0 * pb_y[j] * pc_yzz[j] - 2.0 * pc_yyz[j] * pb_z[j]);

                t_0_yyzz[j] += fl_s_0_0_4 * pc_yyzz[j];

                t_0_yzzz[j] = fl_s_0_0_0 * (1.5 * pb_yz[j] * fl1_fx + pb_yzzz[j]);

                t_0_yzzz[j] += fl_s_0_0_1 * (-1.5 * pb_yz[j] * fl1_fx - 1.5 * pb_y[j] * pc_z[j] * fl1_fx - 1.5 * pc_y[j] * pb_z[j] * fl1_fx - 3.0 * pb_yzz[j] * pc_z[j] - pc_y[j] * pb_zzz[j]);

                t_0_yzzz[j] += fl_s_0_0_2 * (1.5 * pb_y[j] * pc_z[j] * fl1_fx + 1.5 * pc_y[j] * pb_z[j] * fl1_fx + 1.5 * pc_yz[j] * fl1_fx + 3.0 * pb_yz[j] * pc_zz[j] + 3.0 * pc_yz[j] * pb_zz[j]);

                t_0_yzzz[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - pb_y[j] * pc_zzz[j] - 3.0 * pc_yzz[j] * pb_z[j]);

                t_0_yzzz[j] += fl_s_0_0_4 * pc_yzzz[j];

                t_0_zzzz[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 3.0 * pb_zz[j] * fl1_fx + pb_zzzz[j]);

                t_0_zzzz[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 3.0 * pb_zz[j] * fl1_fx - 6.0 * pb_z[j] * pc_z[j] * fl1_fx - 4.0 * pb_zzz[j] * pc_z[j]);

                t_0_zzzz[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 6.0 * pb_z[j] * pc_z[j] * fl1_fx + 3.0 * pc_zz[j] * fl1_fx + 6.0 * pb_zz[j] * pc_zz[j]);

                t_0_zzzz[j] += fl_s_0_0_3 * (-3.0 * pc_zz[j] * fl1_fx - 4.0 * pb_z[j] * pc_zzz[j]);

                t_0_zzzz[j] += fl_s_0_0_4 * pc_zzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGS(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForGS_0_5(primBuffer, auxBuffer, osFactors, paDistances, pcDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGS_5_10(primBuffer, auxBuffer, osFactors, paDistances, pcDistances, 
                                                    braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGS_10_15(primBuffer, auxBuffer, osFactors, paDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compNuclearPotentialForGS_0_5(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (0,5)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(34 * idx);

            auto pc_y = pcDistances.data(34 * idx + 1);

            auto pc_z = pcDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(34 * idx + 3);

            auto pc_xy = pcDistances.data(34 * idx + 4);

            auto pc_xz = pcDistances.data(34 * idx + 5);

            auto pc_yy = pcDistances.data(34 * idx + 6);

            auto pc_yz = pcDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(34 * idx + 9);

            auto pc_xxy = pcDistances.data(34 * idx + 10);

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(34 * idx + 19);

            auto pc_xxxy = pcDistances.data(34 * idx + 20);

            auto pc_xxxz = pcDistances.data(34 * idx + 21);

            auto pc_xxyy = pcDistances.data(34 * idx + 22);

            auto pc_xxyz = pcDistances.data(34 * idx + 23);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

            // set up pointers to integrals

            auto t_xxxx_0 = primBuffer.data(15 * idx);

            auto t_xxxy_0 = primBuffer.data(15 * idx + 1);

            auto t_xxxz_0 = primBuffer.data(15 * idx + 2);

            auto t_xxyy_0 = primBuffer.data(15 * idx + 3);

            auto t_xxyz_0 = primBuffer.data(15 * idx + 4);

            // Batch of Integrals (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pa_xxxy, pa_xxxz, pa_xxy, pa_xxyy, pa_xxyz, \
                                     pa_xxz, pa_xy, pa_xyy, pa_xyz, pa_xz, pa_y, pa_yy, pa_yz, pa_z, pc_x, pc_xx, pc_xxx, \
                                     pc_xxxx, pc_xxxy, pc_xxxz, pc_xxy, pc_xxyy, pc_xxyz, pc_xxz, pc_xy, pc_xyy, pc_xyz, \
                                     pc_xz, pc_y, pc_yy, pc_yz, pc_z, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, \
                                     t_xxxx_0, t_xxxy_0, t_xxxz_0, t_xxyy_0, t_xxyz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxxx_0[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 3.0 * pa_xx[j] * fl1_fx + pa_xxxx[j]);

                t_xxxx_0[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 3.0 * pa_xx[j] * fl1_fx - 6.0 * pa_x[j] * pc_x[j] * fl1_fx - 4.0 * pa_xxx[j] * pc_x[j]);

                t_xxxx_0[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 6.0 * pa_x[j] * pc_x[j] * fl1_fx + 3.0 * pc_xx[j] * fl1_fx + 6.0 * pa_xx[j] * pc_xx[j]);

                t_xxxx_0[j] += fl_s_0_0_3 * (-3.0 * pc_xx[j] * fl1_fx - 4.0 * pa_x[j] * pc_xxx[j]);

                t_xxxx_0[j] += fl_s_0_0_4 * pc_xxxx[j];

                t_xxxy_0[j] = fl_s_0_0_0 * (1.5 * pa_xy[j] * fl1_fx + pa_xxxy[j]);

                t_xxxy_0[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl1_fx * pc_y[j] - 1.5 * pa_xy[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pa_y[j] - pa_xxx[j] * pc_y[j] - 3.0 * pa_xxy[j] * pc_x[j]);

                t_xxxy_0[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * fl1_fx * pc_y[j] + 1.5 * pc_xy[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pa_y[j] + 3.0 * pa_xx[j] * pc_xy[j] + 3.0 * pa_xy[j] * pc_xx[j]);

                t_xxxy_0[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - 3.0 * pa_x[j] * pc_xxy[j] - pc_xxx[j] * pa_y[j]);

                t_xxxy_0[j] += fl_s_0_0_4 * pc_xxxy[j];

                t_xxxz_0[j] = fl_s_0_0_0 * (1.5 * pa_xz[j] * fl1_fx + pa_xxxz[j]);

                t_xxxz_0[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl1_fx * pc_z[j] - 1.5 * pa_xz[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pa_z[j] - pa_xxx[j] * pc_z[j] - 3.0 * pa_xxz[j] * pc_x[j]);

                t_xxxz_0[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * fl1_fx * pc_z[j] + 1.5 * pc_xz[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pa_z[j] + 3.0 * pa_xx[j] * pc_xz[j] + 3.0 * pa_xz[j] * pc_xx[j]);

                t_xxxz_0[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - 3.0 * pa_x[j] * pc_xxz[j] - pc_xxx[j] * pa_z[j]);

                t_xxxz_0[j] += fl_s_0_0_4 * pc_xxxz[j];

                t_xxyy_0[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * fl1_fx * pa_yy[j] + pa_xxyy[j]);

                t_xxyy_0[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_xx[j] * fl1_fx - pa_x[j] * pc_x[j] * fl1_fx - fl1_fx * pa_y[j] * pc_y[j] - 0.5 * fl1_fx * pa_yy[j] - 2.0 * pa_xxy[j] * pc_y[j] - 2.0 * pa_xyy[j] * pc_x[j]);

                t_xxyy_0[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pa_x[j] * pc_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx + 0.5 * fl1_fx * pc_yy[j] + fl1_fx * pa_y[j] * pc_y[j] + pa_xx[j] * pc_yy[j] + 4.0 * pa_xy[j] * pc_xy[j] + pc_xx[j] * pa_yy[j]);

                t_xxyy_0[j] += fl_s_0_0_3 * (-0.5 * pc_xx[j] * fl1_fx - 0.5 * fl1_fx * pc_yy[j] - 2.0 * pa_x[j] * pc_xyy[j] - 2.0 * pc_xxy[j] * pa_y[j]);

                t_xxyy_0[j] += fl_s_0_0_4 * pc_xxyy[j];

                t_xxyz_0[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_yz[j] + pa_xxyz[j]);

                t_xxyz_0[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pa_y[j] * pc_z[j] - 0.5 * fl1_fx * pc_y[j] * pa_z[j] - 0.5 * fl1_fx * pa_yz[j] - pa_xxy[j] * pc_z[j] - pa_xxz[j] * pc_y[j] - 2.0 * pa_xyz[j] * pc_x[j]);

                t_xxyz_0[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_yz[j] + 0.5 * fl1_fx * pa_y[j] * pc_z[j] + 0.5 * fl1_fx * pc_y[j] * pa_z[j] + pa_xx[j] * pc_yz[j] + 2.0 * pa_xy[j] * pc_xz[j] + 2.0 * pa_xz[j] * pc_xy[j] + pc_xx[j] * pa_yz[j]);

                t_xxyz_0[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_yz[j] - 2.0 * pa_x[j] * pc_xyz[j] - pc_xxz[j] * pa_y[j] - pc_xxy[j] * pa_z[j]);

                t_xxyz_0[j] += fl_s_0_0_4 * pc_xxyz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGS_5_10(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (5,10)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(34 * idx);

            auto pc_y = pcDistances.data(34 * idx + 1);

            auto pc_z = pcDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(34 * idx + 3);

            auto pc_xy = pcDistances.data(34 * idx + 4);

            auto pc_xz = pcDistances.data(34 * idx + 5);

            auto pc_yy = pcDistances.data(34 * idx + 6);

            auto pc_yz = pcDistances.data(34 * idx + 7);

            auto pc_zz = pcDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            auto pc_xzz = pcDistances.data(34 * idx + 14);

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxzz = pcDistances.data(34 * idx + 24);

            auto pc_xyyy = pcDistances.data(34 * idx + 25);

            auto pc_xyyz = pcDistances.data(34 * idx + 26);

            auto pc_xyzz = pcDistances.data(34 * idx + 27);

            auto pc_xzzz = pcDistances.data(34 * idx + 28);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

            // set up pointers to integrals

            auto t_xxzz_0 = primBuffer.data(15 * idx + 5);

            auto t_xyyy_0 = primBuffer.data(15 * idx + 6);

            auto t_xyyz_0 = primBuffer.data(15 * idx + 7);

            auto t_xyzz_0 = primBuffer.data(15 * idx + 8);

            auto t_xzzz_0 = primBuffer.data(15 * idx + 9);

            // Batch of Integrals (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xy, pa_xyy, pa_xyyy, pa_xyyz, pa_xyz, \
                                     pa_xyzz, pa_xz, pa_xzz, pa_xzzz, pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pa_zzz, pc_x, pc_xx, pc_xxz, pc_xxzz, pc_xy, pc_xyy, pc_xyyy, pc_xyyz, pc_xyz, pc_xyzz, \
                                     pc_xz, pc_xzz, pc_xzzz, pc_y, pc_yy, pc_yyy, pc_yyz, pc_yz, pc_yzz, pc_z, pc_zz, pc_zzz, \
                                     s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, t_xxzz_0, t_xyyy_0, t_xyyz_0, t_xyzz_0, \
                                     t_xzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxzz_0[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * fl1_fx * pa_zz[j] + pa_xxzz[j]);

                t_xxzz_0[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_xx[j] * fl1_fx - pa_x[j] * pc_x[j] * fl1_fx - fl1_fx * pa_z[j] * pc_z[j] - 0.5 * fl1_fx * pa_zz[j] - 2.0 * pa_xxz[j] * pc_z[j] - 2.0 * pa_xzz[j] * pc_x[j]);

                t_xxzz_0[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pa_x[j] * pc_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx + 0.5 * fl1_fx * pc_zz[j] + fl1_fx * pa_z[j] * pc_z[j] + pa_xx[j] * pc_zz[j] + 4.0 * pa_xz[j] * pc_xz[j] + pc_xx[j] * pa_zz[j]);

                t_xxzz_0[j] += fl_s_0_0_3 * (-0.5 * pc_xx[j] * fl1_fx - 0.5 * fl1_fx * pc_zz[j] - 2.0 * pa_x[j] * pc_xzz[j] - 2.0 * pc_xxz[j] * pa_z[j]);

                t_xxzz_0[j] += fl_s_0_0_4 * pc_xxzz[j];

                t_xyyy_0[j] = fl_s_0_0_0 * (1.5 * pa_xy[j] * fl1_fx + pa_xyyy[j]);

                t_xyyy_0[j] += fl_s_0_0_1 * (-1.5 * pa_xy[j] * fl1_fx - 1.5 * pa_x[j] * pc_y[j] * fl1_fx - 1.5 * pc_x[j] * pa_y[j] * fl1_fx - 3.0 * pa_xyy[j] * pc_y[j] - pc_x[j] * pa_yyy[j]);

                t_xyyy_0[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * pc_y[j] * fl1_fx + 1.5 * pc_x[j] * pa_y[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx + 3.0 * pa_xy[j] * pc_yy[j] + 3.0 * pc_xy[j] * pa_yy[j]);

                t_xyyy_0[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - pa_x[j] * pc_yyy[j] - 3.0 * pc_xyy[j] * pa_y[j]);

                t_xyyy_0[j] += fl_s_0_0_4 * pc_xyyy[j];

                t_xyyz_0[j] = fl_s_0_0_0 * (0.5 * pa_xz[j] * fl1_fx + pa_xyyz[j]);

                t_xyyz_0[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx * pc_z[j] - 0.5 * pa_xz[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx * pa_z[j] - pa_xyy[j] * pc_z[j] - 2.0 * pa_xyz[j] * pc_y[j] - pc_x[j] * pa_yyz[j]);

                t_xyyz_0[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl1_fx * pc_z[j] + 0.5 * pc_xz[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pa_z[j] + 2.0 * pa_xy[j] * pc_yz[j] + pa_xz[j] * pc_yy[j] + pc_xz[j] * pa_yy[j] + 2.0 * pc_xy[j] * pa_yz[j]);

                t_xyyz_0[j] += fl_s_0_0_3 * (-0.5 * pc_xz[j] * fl1_fx - pa_x[j] * pc_yyz[j] - 2.0 * pc_xyz[j] * pa_y[j] - pc_xyy[j] * pa_z[j]);

                t_xyyz_0[j] += fl_s_0_0_4 * pc_xyyz[j];

                t_xyzz_0[j] = fl_s_0_0_0 * (0.5 * pa_xy[j] * fl1_fx + pa_xyzz[j]);

                t_xyzz_0[j] += fl_s_0_0_1 * (-0.5 * pa_xy[j] * fl1_fx - 0.5 * pa_x[j] * pc_y[j] * fl1_fx - 0.5 * pc_x[j] * pa_y[j] * fl1_fx - 2.0 * pa_xyz[j] * pc_z[j] - pa_xzz[j] * pc_y[j] - pc_x[j] * pa_yzz[j]);

                t_xyzz_0[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * pc_y[j] * fl1_fx + 0.5 * pc_x[j] * pa_y[j] * fl1_fx + 0.5 * pc_xy[j] * fl1_fx + pa_xy[j] * pc_zz[j] + 2.0 * pa_xz[j] * pc_yz[j] + 2.0 * pc_xz[j] * pa_yz[j] + pc_xy[j] * pa_zz[j]);

                t_xyzz_0[j] += fl_s_0_0_3 * (-0.5 * pc_xy[j] * fl1_fx - pa_x[j] * pc_yzz[j] - pc_xzz[j] * pa_y[j] - 2.0 * pc_xyz[j] * pa_z[j]);

                t_xyzz_0[j] += fl_s_0_0_4 * pc_xyzz[j];

                t_xzzz_0[j] = fl_s_0_0_0 * (1.5 * pa_xz[j] * fl1_fx + pa_xzzz[j]);

                t_xzzz_0[j] += fl_s_0_0_1 * (-1.5 * pa_xz[j] * fl1_fx - 1.5 * pa_x[j] * pc_z[j] * fl1_fx - 1.5 * pc_x[j] * pa_z[j] * fl1_fx - 3.0 * pa_xzz[j] * pc_z[j] - pc_x[j] * pa_zzz[j]);

                t_xzzz_0[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * pc_z[j] * fl1_fx + 1.5 * pc_x[j] * pa_z[j] * fl1_fx + 1.5 * pc_xz[j] * fl1_fx + 3.0 * pa_xz[j] * pc_zz[j] + 3.0 * pc_xz[j] * pa_zz[j]);

                t_xzzz_0[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - pa_x[j] * pc_zzz[j] - 3.0 * pc_xzz[j] * pa_z[j]);

                t_xzzz_0[j] += fl_s_0_0_4 * pc_xzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGS_10_15(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (10,15)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_y = pcDistances.data(34 * idx + 1);

            auto pc_z = pcDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_yy = pcDistances.data(34 * idx + 6);

            auto pc_yz = pcDistances.data(34 * idx + 7);

            auto pc_zz = pcDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_yyyy = pcDistances.data(34 * idx + 29);

            auto pc_yyyz = pcDistances.data(34 * idx + 30);

            auto pc_yyzz = pcDistances.data(34 * idx + 31);

            auto pc_yzzz = pcDistances.data(34 * idx + 32);

            auto pc_zzzz = pcDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

            // set up pointers to integrals

            auto t_yyyy_0 = primBuffer.data(15 * idx + 10);

            auto t_yyyz_0 = primBuffer.data(15 * idx + 11);

            auto t_yyzz_0 = primBuffer.data(15 * idx + 12);

            auto t_yzzz_0 = primBuffer.data(15 * idx + 13);

            auto t_zzzz_0 = primBuffer.data(15 * idx + 14);

            // Batch of Integrals (10,15)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pa_yyyz, pa_yyz, pa_yyzz, pa_yz, pa_yzz, \
                                     pa_yzzz, pa_z, pa_zz, pa_zzz, pa_zzzz, pc_y, pc_yy, pc_yyy, pc_yyyy, pc_yyyz, pc_yyz, \
                                     pc_yyzz, pc_yz, pc_yzz, pc_yzzz, pc_z, pc_zz, pc_zzz, pc_zzzz, s_0_0_0, s_0_0_1, \
                                     s_0_0_2, s_0_0_3, s_0_0_4, t_yyyy_0, t_yyyz_0, t_yyzz_0, t_yzzz_0, t_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yyyy_0[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 3.0 * pa_yy[j] * fl1_fx + pa_yyyy[j]);

                t_yyyy_0[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 3.0 * pa_yy[j] * fl1_fx - 6.0 * pa_y[j] * pc_y[j] * fl1_fx - 4.0 * pa_yyy[j] * pc_y[j]);

                t_yyyy_0[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 6.0 * pa_y[j] * pc_y[j] * fl1_fx + 3.0 * pc_yy[j] * fl1_fx + 6.0 * pa_yy[j] * pc_yy[j]);

                t_yyyy_0[j] += fl_s_0_0_3 * (-3.0 * pc_yy[j] * fl1_fx - 4.0 * pa_y[j] * pc_yyy[j]);

                t_yyyy_0[j] += fl_s_0_0_4 * pc_yyyy[j];

                t_yyyz_0[j] = fl_s_0_0_0 * (1.5 * pa_yz[j] * fl1_fx + pa_yyyz[j]);

                t_yyyz_0[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl1_fx * pc_z[j] - 1.5 * pa_yz[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx * pa_z[j] - pa_yyy[j] * pc_z[j] - 3.0 * pa_yyz[j] * pc_y[j]);

                t_yyyz_0[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * fl1_fx * pc_z[j] + 1.5 * pc_yz[j] * fl1_fx + 1.5 * pc_y[j] * fl1_fx * pa_z[j] + 3.0 * pa_yy[j] * pc_yz[j] + 3.0 * pa_yz[j] * pc_yy[j]);

                t_yyyz_0[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - 3.0 * pa_y[j] * pc_yyz[j] - pc_yyy[j] * pa_z[j]);

                t_yyyz_0[j] += fl_s_0_0_4 * pc_yyyz[j];

                t_yyzz_0[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_yy[j] * fl1_fx + 0.5 * fl1_fx * pa_zz[j] + pa_yyzz[j]);

                t_yyzz_0[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_yy[j] * fl1_fx - pa_y[j] * pc_y[j] * fl1_fx - fl1_fx * pa_z[j] * pc_z[j] - 0.5 * fl1_fx * pa_zz[j] - 2.0 * pa_yyz[j] * pc_z[j] - 2.0 * pa_yzz[j] * pc_y[j]);

                t_yyzz_0[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pa_y[j] * pc_y[j] * fl1_fx + 0.5 * pc_yy[j] * fl1_fx + 0.5 * fl1_fx * pc_zz[j] + fl1_fx * pa_z[j] * pc_z[j] + pa_yy[j] * pc_zz[j] + 4.0 * pa_yz[j] * pc_yz[j] + pc_yy[j] * pa_zz[j]);

                t_yyzz_0[j] += fl_s_0_0_3 * (-0.5 * pc_yy[j] * fl1_fx - 0.5 * fl1_fx * pc_zz[j] - 2.0 * pa_y[j] * pc_yzz[j] - 2.0 * pc_yyz[j] * pa_z[j]);

                t_yyzz_0[j] += fl_s_0_0_4 * pc_yyzz[j];

                t_yzzz_0[j] = fl_s_0_0_0 * (1.5 * pa_yz[j] * fl1_fx + pa_yzzz[j]);

                t_yzzz_0[j] += fl_s_0_0_1 * (-1.5 * pa_yz[j] * fl1_fx - 1.5 * pa_y[j] * pc_z[j] * fl1_fx - 1.5 * pc_y[j] * pa_z[j] * fl1_fx - 3.0 * pa_yzz[j] * pc_z[j] - pc_y[j] * pa_zzz[j]);

                t_yzzz_0[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * pc_z[j] * fl1_fx + 1.5 * pc_y[j] * pa_z[j] * fl1_fx + 1.5 * pc_yz[j] * fl1_fx + 3.0 * pa_yz[j] * pc_zz[j] + 3.0 * pc_yz[j] * pa_zz[j]);

                t_yzzz_0[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - pa_y[j] * pc_zzz[j] - 3.0 * pc_yzz[j] * pa_z[j]);

                t_yzzz_0[j] += fl_s_0_0_4 * pc_yzzz[j];

                t_zzzz_0[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 3.0 * pa_zz[j] * fl1_fx + pa_zzzz[j]);

                t_zzzz_0[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 3.0 * pa_zz[j] * fl1_fx - 6.0 * pa_z[j] * pc_z[j] * fl1_fx - 4.0 * pa_zzz[j] * pc_z[j]);

                t_zzzz_0[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 6.0 * pa_z[j] * pc_z[j] * fl1_fx + 3.0 * pc_zz[j] * fl1_fx + 6.0 * pa_zz[j] * pc_zz[j]);

                t_zzzz_0[j] += fl_s_0_0_3 * (-3.0 * pc_zz[j] * fl1_fx - 4.0 * pa_z[j] * pc_zzz[j]);

                t_zzzz_0[j] += fl_s_0_0_4 * pc_zzzz[j];
            }

            idx++;
        }
    }


} // npotrecfunc namespace

