//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFuncForPX.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace

    void
    compOverlapForPP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // Batch of Integrals (0,9)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(9 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(9 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(9 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(9 * idx + 3);

            auto pa2pb_y_y = pa2pbDistances.data(9 * idx + 4);

            auto pa2pb_y_z = pa2pbDistances.data(9 * idx + 5);

            auto pa2pb_z_x = pa2pbDistances.data(9 * idx + 6);

            auto pa2pb_z_y = pa2pbDistances.data(9 * idx + 7);

            auto pa2pb_z_z = pa2pbDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_x = primBuffer.data(9 * idx);

            auto t_x_y = primBuffer.data(9 * idx + 1);

            auto t_x_z = primBuffer.data(9 * idx + 2);

            auto t_y_x = primBuffer.data(9 * idx + 3);

            auto t_y_y = primBuffer.data(9 * idx + 4);

            auto t_y_z = primBuffer.data(9 * idx + 5);

            auto t_z_x = primBuffer.data(9 * idx + 6);

            auto t_z_y = primBuffer.data(9 * idx + 7);

            auto t_z_z = primBuffer.data(9 * idx + 8);

            // Batch of Integrals (0,9)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, \
                                     pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, s_0_0, t_x_x, t_x_y, t_x_z, t_y_x, t_y_y, t_y_z, t_z_x, \
                                     t_z_y, t_z_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_x_x[j] = fl_s_0_0 * (0.5 * fl1_fx + pa2pb_x_x[j]);

                t_x_y[j] = fl_s_0_0 * pa2pb_x_y[j];

                t_x_z[j] = fl_s_0_0 * pa2pb_x_z[j];

                t_y_x[j] = fl_s_0_0 * pa2pb_y_x[j];

                t_y_y[j] = fl_s_0_0 * (0.5 * fl1_fx + pa2pb_y_y[j]);

                t_y_z[j] = fl_s_0_0 * pa2pb_y_z[j];

                t_z_x[j] = fl_s_0_0 * pa2pb_z_x[j];

                t_z_y[j] = fl_s_0_0 * pa2pb_z_y[j];

                t_z_z[j] = fl_s_0_0 * (0.5 * fl1_fx + pa2pb_z_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForPD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForPD_0_9(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForPD_9_18(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForPD_0_9(      CMemBlock2D<double>& primBuffer,
                         const CMemBlock2D<double>& auxBuffer,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CMemBlock2D<double>& pbDistances,
                         const CMemBlock2D<double>& pa2pbDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
    {
        // Batch of Integrals (0,9)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(27 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(27 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(27 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(27 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(27 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(27 * idx + 8);

            auto pa2pb_y_xx = pa2pbDistances.data(27 * idx + 12);

            auto pa2pb_y_xy = pa2pbDistances.data(27 * idx + 13);

            auto pa2pb_y_xz = pa2pbDistances.data(27 * idx + 14);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_xx = primBuffer.data(18 * idx);

            auto t_x_xy = primBuffer.data(18 * idx + 1);

            auto t_x_xz = primBuffer.data(18 * idx + 2);

            auto t_x_yy = primBuffer.data(18 * idx + 3);

            auto t_x_yz = primBuffer.data(18 * idx + 4);

            auto t_x_zz = primBuffer.data(18 * idx + 5);

            auto t_y_xx = primBuffer.data(18 * idx + 6);

            auto t_y_xy = primBuffer.data(18 * idx + 7);

            auto t_y_xz = primBuffer.data(18 * idx + 8);

            // Batch of Integrals (0,9)

            #pragma omp simd aligned(fx, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_x_zz, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_xz, pa_x, pa_y, pb_x, pb_y, pb_z, s_0_0, \
                                     t_x_xx, t_x_xy, t_x_xz, t_x_yy, t_x_yz, t_x_zz, t_y_xx, t_y_xy, t_y_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_x_xx[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl1_fx + pb_x[j] * fl1_fx + pa2pb_x_xx[j]);

                t_x_xy[j] = fl_s_0_0 * (0.5 * pb_y[j] * fl1_fx + pa2pb_x_xy[j]);

                t_x_xz[j] = fl_s_0_0 * (0.5 * pb_z[j] * fl1_fx + pa2pb_x_xz[j]);

                t_x_yy[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl1_fx + pa2pb_x_yy[j]);

                t_x_yz[j] = fl_s_0_0 * pa2pb_x_yz[j];

                t_x_zz[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl1_fx + pa2pb_x_zz[j]);

                t_y_xx[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl1_fx + pa2pb_y_xx[j]);

                t_y_xy[j] = fl_s_0_0 * (0.5 * pb_x[j] * fl1_fx + pa2pb_y_xy[j]);

                t_y_xz[j] = fl_s_0_0 * pa2pb_y_xz[j];
            }

            idx++;
        }
    }

    void
    compOverlapForPD_9_18(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
    {
        // Batch of Integrals (9,18)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_yy = pa2pbDistances.data(27 * idx + 15);

            auto pa2pb_y_yz = pa2pbDistances.data(27 * idx + 16);

            auto pa2pb_y_zz = pa2pbDistances.data(27 * idx + 17);

            auto pa2pb_z_xx = pa2pbDistances.data(27 * idx + 21);

            auto pa2pb_z_xy = pa2pbDistances.data(27 * idx + 22);

            auto pa2pb_z_xz = pa2pbDistances.data(27 * idx + 23);

            auto pa2pb_z_yy = pa2pbDistances.data(27 * idx + 24);

            auto pa2pb_z_yz = pa2pbDistances.data(27 * idx + 25);

            auto pa2pb_z_zz = pa2pbDistances.data(27 * idx + 26);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_y_yy = primBuffer.data(18 * idx + 9);

            auto t_y_yz = primBuffer.data(18 * idx + 10);

            auto t_y_zz = primBuffer.data(18 * idx + 11);

            auto t_z_xx = primBuffer.data(18 * idx + 12);

            auto t_z_xy = primBuffer.data(18 * idx + 13);

            auto t_z_xz = primBuffer.data(18 * idx + 14);

            auto t_z_yy = primBuffer.data(18 * idx + 15);

            auto t_z_yz = primBuffer.data(18 * idx + 16);

            auto t_z_zz = primBuffer.data(18 * idx + 17);

            // Batch of Integrals (9,18)

            #pragma omp simd aligned(fx, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_z_xx, pa2pb_z_xy, \
                                     pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, \
                                     t_y_yy, t_y_yz, t_y_zz, t_z_xx, t_z_xy, t_z_xz, t_z_yy, t_z_yz, t_z_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_y_yy[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl1_fx + pb_y[j] * fl1_fx + pa2pb_y_yy[j]);

                t_y_yz[j] = fl_s_0_0 * (0.5 * pb_z[j] * fl1_fx + pa2pb_y_yz[j]);

                t_y_zz[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl1_fx + pa2pb_y_zz[j]);

                t_z_xx[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl1_fx + pa2pb_z_xx[j]);

                t_z_xy[j] = fl_s_0_0 * pa2pb_z_xy[j];

                t_z_xz[j] = fl_s_0_0 * (0.5 * pb_x[j] * fl1_fx + pa2pb_z_xz[j]);

                t_z_yy[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl1_fx + pa2pb_z_yy[j]);

                t_z_yz[j] = fl_s_0_0 * (0.5 * pb_y[j] * fl1_fx + pa2pb_z_yz[j]);

                t_z_zz[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl1_fx + pb_z[j] * fl1_fx + pa2pb_z_zz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForDP_0_9(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDP_9_18(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForDP_0_9(      CMemBlock2D<double>& primBuffer,
                         const CMemBlock2D<double>& auxBuffer,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CMemBlock2D<double>& pbDistances,
                         const CMemBlock2D<double>& pa2pbDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
    {
        // Batch of Integrals (0,9)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(9 * idx);

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_xx_x = pa2pbDistances.data(27 * idx + 9);

            auto pa2pb_xx_y = pa2pbDistances.data(27 * idx + 10);

            auto pa2pb_xx_z = pa2pbDistances.data(27 * idx + 11);

            auto pa2pb_xy_x = pa2pbDistances.data(27 * idx + 12);

            auto pa2pb_xy_y = pa2pbDistances.data(27 * idx + 13);

            auto pa2pb_xy_z = pa2pbDistances.data(27 * idx + 14);

            auto pa2pb_xz_x = pa2pbDistances.data(27 * idx + 15);

            auto pa2pb_xz_y = pa2pbDistances.data(27 * idx + 16);

            auto pa2pb_xz_z = pa2pbDistances.data(27 * idx + 17);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_x = primBuffer.data(18 * idx);

            auto t_xx_y = primBuffer.data(18 * idx + 1);

            auto t_xx_z = primBuffer.data(18 * idx + 2);

            auto t_xy_x = primBuffer.data(18 * idx + 3);

            auto t_xy_y = primBuffer.data(18 * idx + 4);

            auto t_xy_z = primBuffer.data(18 * idx + 5);

            auto t_xz_x = primBuffer.data(18 * idx + 6);

            auto t_xz_y = primBuffer.data(18 * idx + 7);

            auto t_xz_z = primBuffer.data(18 * idx + 8);

            // Batch of Integrals (0,9)

            #pragma omp simd aligned(fx, pa2pb_xx_x, pa2pb_xx_y, pa2pb_xx_z, pa2pb_xy_x, pa2pb_xy_y, \
                                     pa2pb_xy_z, pa2pb_xz_x, pa2pb_xz_y, pa2pb_xz_z, pa_x, pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, \
                                     t_xx_x, t_xx_y, t_xx_z, t_xy_x, t_xy_y, t_xy_z, t_xz_x, t_xz_y, t_xz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_xx_x[j] = fl_s_0_0 * (pa_x[j] * fl1_fx + 0.5 * pb_x[j] * fl1_fx + pa2pb_xx_x[j]);

                t_xx_y[j] = fl_s_0_0 * (0.5 * pb_y[j] * fl1_fx + pa2pb_xx_y[j]);

                t_xx_z[j] = fl_s_0_0 * (0.5 * pb_z[j] * fl1_fx + pa2pb_xx_z[j]);

                t_xy_x[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl1_fx + pa2pb_xy_x[j]);

                t_xy_y[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl1_fx + pa2pb_xy_y[j]);

                t_xy_z[j] = fl_s_0_0 * pa2pb_xy_z[j];

                t_xz_x[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl1_fx + pa2pb_xz_x[j]);

                t_xz_y[j] = fl_s_0_0 * pa2pb_xz_y[j];

                t_xz_z[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl1_fx + pa2pb_xz_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDP_9_18(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
    {
        // Batch of Integrals (9,18)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_yy_x = pa2pbDistances.data(27 * idx + 18);

            auto pa2pb_yy_y = pa2pbDistances.data(27 * idx + 19);

            auto pa2pb_yy_z = pa2pbDistances.data(27 * idx + 20);

            auto pa2pb_yz_x = pa2pbDistances.data(27 * idx + 21);

            auto pa2pb_yz_y = pa2pbDistances.data(27 * idx + 22);

            auto pa2pb_yz_z = pa2pbDistances.data(27 * idx + 23);

            auto pa2pb_zz_x = pa2pbDistances.data(27 * idx + 24);

            auto pa2pb_zz_y = pa2pbDistances.data(27 * idx + 25);

            auto pa2pb_zz_z = pa2pbDistances.data(27 * idx + 26);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yy_x = primBuffer.data(18 * idx + 9);

            auto t_yy_y = primBuffer.data(18 * idx + 10);

            auto t_yy_z = primBuffer.data(18 * idx + 11);

            auto t_yz_x = primBuffer.data(18 * idx + 12);

            auto t_yz_y = primBuffer.data(18 * idx + 13);

            auto t_yz_z = primBuffer.data(18 * idx + 14);

            auto t_zz_x = primBuffer.data(18 * idx + 15);

            auto t_zz_y = primBuffer.data(18 * idx + 16);

            auto t_zz_z = primBuffer.data(18 * idx + 17);

            // Batch of Integrals (9,18)

            #pragma omp simd aligned(fx, pa2pb_yy_x, pa2pb_yy_y, pa2pb_yy_z, pa2pb_yz_x, pa2pb_yz_y, \
                                     pa2pb_yz_z, pa2pb_zz_x, pa2pb_zz_y, pa2pb_zz_z, pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, \
                                     t_yy_x, t_yy_y, t_yy_z, t_yz_x, t_yz_y, t_yz_z, t_zz_x, t_zz_y, t_zz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_yy_x[j] = fl_s_0_0 * (0.5 * pb_x[j] * fl1_fx + pa2pb_yy_x[j]);

                t_yy_y[j] = fl_s_0_0 * (pa_y[j] * fl1_fx + 0.5 * pb_y[j] * fl1_fx + pa2pb_yy_y[j]);

                t_yy_z[j] = fl_s_0_0 * (0.5 * pb_z[j] * fl1_fx + pa2pb_yy_z[j]);

                t_yz_x[j] = fl_s_0_0 * pa2pb_yz_x[j];

                t_yz_y[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl1_fx + pa2pb_yz_y[j]);

                t_yz_z[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl1_fx + pa2pb_yz_z[j]);

                t_zz_x[j] = fl_s_0_0 * (0.5 * pb_x[j] * fl1_fx + pa2pb_zz_x[j]);

                t_zz_y[j] = fl_s_0_0 * (0.5 * pb_y[j] * fl1_fx + pa2pb_zz_y[j]);

                t_zz_z[j] = fl_s_0_0 * (pa_z[j] * fl1_fx + 0.5 * pb_z[j] * fl1_fx + pa2pb_zz_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForPF(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForPF_0_15(primBuffer, auxBuffer, osFactors, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForPF_15_30(primBuffer, auxBuffer, osFactors, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForPF_0_15(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
    {
        // Batch of Integrals (0,15)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(57 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(57 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(57 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(57 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(57 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(57 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(57 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(57 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(57 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(57 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(57 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(57 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(57 * idx + 18);

            auto pa2pb_y_x = pa2pbDistances.data(57 * idx + 19);

            auto pa2pb_y_y = pa2pbDistances.data(57 * idx + 20);

            auto pa2pb_y_z = pa2pbDistances.data(57 * idx + 21);

            auto pa2pb_y_xxx = pa2pbDistances.data(57 * idx + 28);

            auto pa2pb_y_xxy = pa2pbDistances.data(57 * idx + 29);

            auto pa2pb_y_xxz = pa2pbDistances.data(57 * idx + 30);

            auto pa2pb_y_xyy = pa2pbDistances.data(57 * idx + 31);

            auto pa2pb_y_xyz = pa2pbDistances.data(57 * idx + 32);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_xxx = primBuffer.data(30 * idx);

            auto t_x_xxy = primBuffer.data(30 * idx + 1);

            auto t_x_xxz = primBuffer.data(30 * idx + 2);

            auto t_x_xyy = primBuffer.data(30 * idx + 3);

            auto t_x_xyz = primBuffer.data(30 * idx + 4);

            auto t_x_xzz = primBuffer.data(30 * idx + 5);

            auto t_x_yyy = primBuffer.data(30 * idx + 6);

            auto t_x_yyz = primBuffer.data(30 * idx + 7);

            auto t_x_yzz = primBuffer.data(30 * idx + 8);

            auto t_x_zzz = primBuffer.data(30 * idx + 9);

            auto t_y_xxx = primBuffer.data(30 * idx + 10);

            auto t_y_xxy = primBuffer.data(30 * idx + 11);

            auto t_y_xxz = primBuffer.data(30 * idx + 12);

            auto t_y_xyy = primBuffer.data(30 * idx + 13);

            auto t_y_xyz = primBuffer.data(30 * idx + 14);

            // Batch of Integrals (0,15)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, pa2pb_x_xyy, \
                                     pa2pb_x_xyz, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, \
                                     pa2pb_x_z, pa2pb_x_zzz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, \
                                     pa2pb_y_xyy, pa2pb_y_xyz, pa2pb_y_y, pa2pb_y_z, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, \
                                     s_0_0, t_x_xxx, t_x_xxy, t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy, t_x_yyz, \
                                     t_x_yzz, t_x_zzz, t_y_xxx, t_y_xxy, t_y_xxz, t_y_xyy, t_y_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_x_xxx[j] = fl_s_0_0 * (0.75 * fl2_fx + 1.5 * pa2pb_x_x[j] * fl1_fx + 1.5 * pb_xx[j] * fl1_fx + pa2pb_x_xxx[j]);

                t_x_xxy[j] = fl_s_0_0 * (0.5 * pa2pb_x_y[j] * fl1_fx + pb_xy[j] * fl1_fx + pa2pb_x_xxy[j]);

                t_x_xxz[j] = fl_s_0_0 * (0.5 * pa2pb_x_z[j] * fl1_fx + pb_xz[j] * fl1_fx + pa2pb_x_xxz[j]);

                t_x_xyy[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_x_x[j] * fl1_fx + 0.5 * pb_yy[j] * fl1_fx + pa2pb_x_xyy[j]);

                t_x_xyz[j] = fl_s_0_0 * (0.5 * pb_yz[j] * fl1_fx + pa2pb_x_xyz[j]);

                t_x_xzz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_x_x[j] * fl1_fx + 0.5 * pb_zz[j] * fl1_fx + pa2pb_x_xzz[j]);

                t_x_yyy[j] = fl_s_0_0 * (1.5 * pa2pb_x_y[j] * fl1_fx + pa2pb_x_yyy[j]);

                t_x_yyz[j] = fl_s_0_0 * (0.5 * pa2pb_x_z[j] * fl1_fx + pa2pb_x_yyz[j]);

                t_x_yzz[j] = fl_s_0_0 * (0.5 * pa2pb_x_y[j] * fl1_fx + pa2pb_x_yzz[j]);

                t_x_zzz[j] = fl_s_0_0 * (1.5 * pa2pb_x_z[j] * fl1_fx + pa2pb_x_zzz[j]);

                t_y_xxx[j] = fl_s_0_0 * (1.5 * pa2pb_y_x[j] * fl1_fx + pa2pb_y_xxx[j]);

                t_y_xxy[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_y_y[j] * fl1_fx + 0.5 * pb_xx[j] * fl1_fx + pa2pb_y_xxy[j]);

                t_y_xxz[j] = fl_s_0_0 * (0.5 * pa2pb_y_z[j] * fl1_fx + pa2pb_y_xxz[j]);

                t_y_xyy[j] = fl_s_0_0 * (0.5 * pa2pb_y_x[j] * fl1_fx + pb_xy[j] * fl1_fx + pa2pb_y_xyy[j]);

                t_y_xyz[j] = fl_s_0_0 * (0.5 * pb_xz[j] * fl1_fx + pa2pb_y_xyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForPF_15_30(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (15,30)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(57 * idx + 19);

            auto pa2pb_y_y = pa2pbDistances.data(57 * idx + 20);

            auto pa2pb_y_z = pa2pbDistances.data(57 * idx + 21);

            auto pa2pb_y_xzz = pa2pbDistances.data(57 * idx + 33);

            auto pa2pb_y_yyy = pa2pbDistances.data(57 * idx + 34);

            auto pa2pb_y_yyz = pa2pbDistances.data(57 * idx + 35);

            auto pa2pb_y_yzz = pa2pbDistances.data(57 * idx + 36);

            auto pa2pb_y_zzz = pa2pbDistances.data(57 * idx + 37);

            auto pa2pb_z_x = pa2pbDistances.data(57 * idx + 38);

            auto pa2pb_z_y = pa2pbDistances.data(57 * idx + 39);

            auto pa2pb_z_z = pa2pbDistances.data(57 * idx + 40);

            auto pa2pb_z_xxx = pa2pbDistances.data(57 * idx + 47);

            auto pa2pb_z_xxy = pa2pbDistances.data(57 * idx + 48);

            auto pa2pb_z_xxz = pa2pbDistances.data(57 * idx + 49);

            auto pa2pb_z_xyy = pa2pbDistances.data(57 * idx + 50);

            auto pa2pb_z_xyz = pa2pbDistances.data(57 * idx + 51);

            auto pa2pb_z_xzz = pa2pbDistances.data(57 * idx + 52);

            auto pa2pb_z_yyy = pa2pbDistances.data(57 * idx + 53);

            auto pa2pb_z_yyz = pa2pbDistances.data(57 * idx + 54);

            auto pa2pb_z_yzz = pa2pbDistances.data(57 * idx + 55);

            auto pa2pb_z_zzz = pa2pbDistances.data(57 * idx + 56);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_y_xzz = primBuffer.data(30 * idx + 15);

            auto t_y_yyy = primBuffer.data(30 * idx + 16);

            auto t_y_yyz = primBuffer.data(30 * idx + 17);

            auto t_y_yzz = primBuffer.data(30 * idx + 18);

            auto t_y_zzz = primBuffer.data(30 * idx + 19);

            auto t_z_xxx = primBuffer.data(30 * idx + 20);

            auto t_z_xxy = primBuffer.data(30 * idx + 21);

            auto t_z_xxz = primBuffer.data(30 * idx + 22);

            auto t_z_xyy = primBuffer.data(30 * idx + 23);

            auto t_z_xyz = primBuffer.data(30 * idx + 24);

            auto t_z_xzz = primBuffer.data(30 * idx + 25);

            auto t_z_yyy = primBuffer.data(30 * idx + 26);

            auto t_z_yyz = primBuffer.data(30 * idx + 27);

            auto t_z_yzz = primBuffer.data(30 * idx + 28);

            auto t_z_zzz = primBuffer.data(30 * idx + 29);

            // Batch of Integrals (15,30)

            #pragma omp simd aligned(fx, pa2pb_y_x, pa2pb_y_xzz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, \
                                     pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, \
                                     pa2pb_z_xxz, pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_z_y, pa2pb_z_yyy, \
                                     pa2pb_z_yyz, pa2pb_z_yzz, pa2pb_z_z, pa2pb_z_zzz, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, \
                                     pb_zz, s_0_0, t_y_xzz, t_y_yyy, t_y_yyz, t_y_yzz, t_y_zzz, t_z_xxx, t_z_xxy, \
                                     t_z_xxz, t_z_xyy, t_z_xyz, t_z_xzz, t_z_yyy, t_z_yyz, t_z_yzz, t_z_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_y_xzz[j] = fl_s_0_0 * (0.5 * pa2pb_y_x[j] * fl1_fx + pa2pb_y_xzz[j]);

                t_y_yyy[j] = fl_s_0_0 * (0.75 * fl2_fx + 1.5 * pa2pb_y_y[j] * fl1_fx + 1.5 * pb_yy[j] * fl1_fx + pa2pb_y_yyy[j]);

                t_y_yyz[j] = fl_s_0_0 * (0.5 * pa2pb_y_z[j] * fl1_fx + pb_yz[j] * fl1_fx + pa2pb_y_yyz[j]);

                t_y_yzz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_y_y[j] * fl1_fx + 0.5 * pb_zz[j] * fl1_fx + pa2pb_y_yzz[j]);

                t_y_zzz[j] = fl_s_0_0 * (1.5 * pa2pb_y_z[j] * fl1_fx + pa2pb_y_zzz[j]);

                t_z_xxx[j] = fl_s_0_0 * (1.5 * pa2pb_z_x[j] * fl1_fx + pa2pb_z_xxx[j]);

                t_z_xxy[j] = fl_s_0_0 * (0.5 * pa2pb_z_y[j] * fl1_fx + pa2pb_z_xxy[j]);

                t_z_xxz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_z_z[j] * fl1_fx + 0.5 * pb_xx[j] * fl1_fx + pa2pb_z_xxz[j]);

                t_z_xyy[j] = fl_s_0_0 * (0.5 * pa2pb_z_x[j] * fl1_fx + pa2pb_z_xyy[j]);

                t_z_xyz[j] = fl_s_0_0 * (0.5 * pb_xy[j] * fl1_fx + pa2pb_z_xyz[j]);

                t_z_xzz[j] = fl_s_0_0 * (0.5 * pa2pb_z_x[j] * fl1_fx + pb_xz[j] * fl1_fx + pa2pb_z_xzz[j]);

                t_z_yyy[j] = fl_s_0_0 * (1.5 * pa2pb_z_y[j] * fl1_fx + pa2pb_z_yyy[j]);

                t_z_yyz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_z_z[j] * fl1_fx + 0.5 * pb_yy[j] * fl1_fx + pa2pb_z_yyz[j]);

                t_z_yzz[j] = fl_s_0_0 * (0.5 * pa2pb_z_y[j] * fl1_fx + pb_yz[j] * fl1_fx + pa2pb_z_yzz[j]);

                t_z_zzz[j] = fl_s_0_0 * (0.75 * fl2_fx + 1.5 * pa2pb_z_z[j] * fl1_fx + 1.5 * pb_zz[j] * fl1_fx + pa2pb_z_zzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForFP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForFP_0_15(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForFP_15_30(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForFP_0_15(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
    {
        // Batch of Integrals (0,15)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(57 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(57 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(57 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(57 * idx + 3);

            auto pa2pb_y_y = pa2pbDistances.data(57 * idx + 4);

            auto pa2pb_y_z = pa2pbDistances.data(57 * idx + 5);

            auto pa2pb_z_x = pa2pbDistances.data(57 * idx + 6);

            auto pa2pb_z_y = pa2pbDistances.data(57 * idx + 7);

            auto pa2pb_z_z = pa2pbDistances.data(57 * idx + 8);

            auto pa2pb_xxx_x = pa2pbDistances.data(57 * idx + 27);

            auto pa2pb_xxx_y = pa2pbDistances.data(57 * idx + 28);

            auto pa2pb_xxx_z = pa2pbDistances.data(57 * idx + 29);

            auto pa2pb_xxy_x = pa2pbDistances.data(57 * idx + 30);

            auto pa2pb_xxy_y = pa2pbDistances.data(57 * idx + 31);

            auto pa2pb_xxy_z = pa2pbDistances.data(57 * idx + 32);

            auto pa2pb_xxz_x = pa2pbDistances.data(57 * idx + 33);

            auto pa2pb_xxz_y = pa2pbDistances.data(57 * idx + 34);

            auto pa2pb_xxz_z = pa2pbDistances.data(57 * idx + 35);

            auto pa2pb_xyy_x = pa2pbDistances.data(57 * idx + 36);

            auto pa2pb_xyy_y = pa2pbDistances.data(57 * idx + 37);

            auto pa2pb_xyy_z = pa2pbDistances.data(57 * idx + 38);

            auto pa2pb_xyz_x = pa2pbDistances.data(57 * idx + 39);

            auto pa2pb_xyz_y = pa2pbDistances.data(57 * idx + 40);

            auto pa2pb_xyz_z = pa2pbDistances.data(57 * idx + 41);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxx_x = primBuffer.data(30 * idx);

            auto t_xxx_y = primBuffer.data(30 * idx + 1);

            auto t_xxx_z = primBuffer.data(30 * idx + 2);

            auto t_xxy_x = primBuffer.data(30 * idx + 3);

            auto t_xxy_y = primBuffer.data(30 * idx + 4);

            auto t_xxy_z = primBuffer.data(30 * idx + 5);

            auto t_xxz_x = primBuffer.data(30 * idx + 6);

            auto t_xxz_y = primBuffer.data(30 * idx + 7);

            auto t_xxz_z = primBuffer.data(30 * idx + 8);

            auto t_xyy_x = primBuffer.data(30 * idx + 9);

            auto t_xyy_y = primBuffer.data(30 * idx + 10);

            auto t_xyy_z = primBuffer.data(30 * idx + 11);

            auto t_xyz_x = primBuffer.data(30 * idx + 12);

            auto t_xyz_y = primBuffer.data(30 * idx + 13);

            auto t_xyz_z = primBuffer.data(30 * idx + 14);

            // Batch of Integrals (0,15)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xxx_x, pa2pb_xxx_y, \
                                     pa2pb_xxx_z, pa2pb_xxy_x, pa2pb_xxy_y, pa2pb_xxy_z, pa2pb_xxz_x, pa2pb_xxz_y, \
                                     pa2pb_xxz_z, pa2pb_xyy_x, pa2pb_xyy_y, pa2pb_xyy_z, pa2pb_xyz_x, pa2pb_xyz_y, \
                                     pa2pb_xyz_z, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, pa_xx, \
                                     pa_xy, pa_xz, pa_yy, pa_yz, s_0_0, t_xxx_x, t_xxx_y, t_xxx_z, t_xxy_x, t_xxy_y, \
                                     t_xxy_z, t_xxz_x, t_xxz_y, t_xxz_z, t_xyy_x, t_xyy_y, t_xyy_z, t_xyz_x, t_xyz_y, \
                                     t_xyz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxx_x[j] = fl_s_0_0 * (0.75 * fl2_fx + 1.5 * pa_xx[j] * fl1_fx + 1.5 * pa2pb_x_x[j] * fl1_fx + pa2pb_xxx_x[j]);

                t_xxx_y[j] = fl_s_0_0 * (1.5 * pa2pb_x_y[j] * fl1_fx + pa2pb_xxx_y[j]);

                t_xxx_z[j] = fl_s_0_0 * (1.5 * pa2pb_x_z[j] * fl1_fx + pa2pb_xxx_z[j]);

                t_xxy_x[j] = fl_s_0_0 * (pa_xy[j] * fl1_fx + 0.5 * pa2pb_y_x[j] * fl1_fx + pa2pb_xxy_x[j]);

                t_xxy_y[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * pa2pb_y_y[j] * fl1_fx + pa2pb_xxy_y[j]);

                t_xxy_z[j] = fl_s_0_0 * (0.5 * pa2pb_y_z[j] * fl1_fx + pa2pb_xxy_z[j]);

                t_xxz_x[j] = fl_s_0_0 * (pa_xz[j] * fl1_fx + 0.5 * pa2pb_z_x[j] * fl1_fx + pa2pb_xxz_x[j]);

                t_xxz_y[j] = fl_s_0_0 * (0.5 * pa2pb_z_y[j] * fl1_fx + pa2pb_xxz_y[j]);

                t_xxz_z[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * pa2pb_z_z[j] * fl1_fx + pa2pb_xxz_z[j]);

                t_xyy_x[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_yy[j] * fl1_fx + 0.5 * pa2pb_x_x[j] * fl1_fx + pa2pb_xyy_x[j]);

                t_xyy_y[j] = fl_s_0_0 * (pa_xy[j] * fl1_fx + 0.5 * pa2pb_x_y[j] * fl1_fx + pa2pb_xyy_y[j]);

                t_xyy_z[j] = fl_s_0_0 * (0.5 * pa2pb_x_z[j] * fl1_fx + pa2pb_xyy_z[j]);

                t_xyz_x[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl1_fx + pa2pb_xyz_x[j]);

                t_xyz_y[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl1_fx + pa2pb_xyz_y[j]);

                t_xyz_z[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl1_fx + pa2pb_xyz_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForFP_15_30(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (15,30)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(57 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(57 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(57 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(57 * idx + 3);

            auto pa2pb_y_y = pa2pbDistances.data(57 * idx + 4);

            auto pa2pb_y_z = pa2pbDistances.data(57 * idx + 5);

            auto pa2pb_z_x = pa2pbDistances.data(57 * idx + 6);

            auto pa2pb_z_y = pa2pbDistances.data(57 * idx + 7);

            auto pa2pb_z_z = pa2pbDistances.data(57 * idx + 8);

            auto pa2pb_xzz_x = pa2pbDistances.data(57 * idx + 42);

            auto pa2pb_xzz_y = pa2pbDistances.data(57 * idx + 43);

            auto pa2pb_xzz_z = pa2pbDistances.data(57 * idx + 44);

            auto pa2pb_yyy_x = pa2pbDistances.data(57 * idx + 45);

            auto pa2pb_yyy_y = pa2pbDistances.data(57 * idx + 46);

            auto pa2pb_yyy_z = pa2pbDistances.data(57 * idx + 47);

            auto pa2pb_yyz_x = pa2pbDistances.data(57 * idx + 48);

            auto pa2pb_yyz_y = pa2pbDistances.data(57 * idx + 49);

            auto pa2pb_yyz_z = pa2pbDistances.data(57 * idx + 50);

            auto pa2pb_yzz_x = pa2pbDistances.data(57 * idx + 51);

            auto pa2pb_yzz_y = pa2pbDistances.data(57 * idx + 52);

            auto pa2pb_yzz_z = pa2pbDistances.data(57 * idx + 53);

            auto pa2pb_zzz_x = pa2pbDistances.data(57 * idx + 54);

            auto pa2pb_zzz_y = pa2pbDistances.data(57 * idx + 55);

            auto pa2pb_zzz_z = pa2pbDistances.data(57 * idx + 56);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xzz_x = primBuffer.data(30 * idx + 15);

            auto t_xzz_y = primBuffer.data(30 * idx + 16);

            auto t_xzz_z = primBuffer.data(30 * idx + 17);

            auto t_yyy_x = primBuffer.data(30 * idx + 18);

            auto t_yyy_y = primBuffer.data(30 * idx + 19);

            auto t_yyy_z = primBuffer.data(30 * idx + 20);

            auto t_yyz_x = primBuffer.data(30 * idx + 21);

            auto t_yyz_y = primBuffer.data(30 * idx + 22);

            auto t_yyz_z = primBuffer.data(30 * idx + 23);

            auto t_yzz_x = primBuffer.data(30 * idx + 24);

            auto t_yzz_y = primBuffer.data(30 * idx + 25);

            auto t_yzz_z = primBuffer.data(30 * idx + 26);

            auto t_zzz_x = primBuffer.data(30 * idx + 27);

            auto t_zzz_y = primBuffer.data(30 * idx + 28);

            auto t_zzz_z = primBuffer.data(30 * idx + 29);

            // Batch of Integrals (15,30)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xzz_x, pa2pb_xzz_y, \
                                     pa2pb_xzz_z, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa2pb_yyy_x, pa2pb_yyy_y, pa2pb_yyy_z, \
                                     pa2pb_yyz_x, pa2pb_yyz_y, pa2pb_yyz_z, pa2pb_yzz_x, pa2pb_yzz_y, pa2pb_yzz_z, \
                                     pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, pa2pb_zzz_x, pa2pb_zzz_y, pa2pb_zzz_z, pa_xz, \
                                     pa_yy, pa_yz, pa_zz, s_0_0, t_xzz_x, t_xzz_y, t_xzz_z, t_yyy_x, t_yyy_y, t_yyy_z, \
                                     t_yyz_x, t_yyz_y, t_yyz_z, t_yzz_x, t_yzz_y, t_yzz_z, t_zzz_x, t_zzz_y, t_zzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xzz_x[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_zz[j] * fl1_fx + 0.5 * pa2pb_x_x[j] * fl1_fx + pa2pb_xzz_x[j]);

                t_xzz_y[j] = fl_s_0_0 * (0.5 * pa2pb_x_y[j] * fl1_fx + pa2pb_xzz_y[j]);

                t_xzz_z[j] = fl_s_0_0 * (pa_xz[j] * fl1_fx + 0.5 * pa2pb_x_z[j] * fl1_fx + pa2pb_xzz_z[j]);

                t_yyy_x[j] = fl_s_0_0 * (1.5 * pa2pb_y_x[j] * fl1_fx + pa2pb_yyy_x[j]);

                t_yyy_y[j] = fl_s_0_0 * (0.75 * fl2_fx + 1.5 * pa_yy[j] * fl1_fx + 1.5 * pa2pb_y_y[j] * fl1_fx + pa2pb_yyy_y[j]);

                t_yyy_z[j] = fl_s_0_0 * (1.5 * pa2pb_y_z[j] * fl1_fx + pa2pb_yyy_z[j]);

                t_yyz_x[j] = fl_s_0_0 * (0.5 * pa2pb_z_x[j] * fl1_fx + pa2pb_yyz_x[j]);

                t_yyz_y[j] = fl_s_0_0 * (pa_yz[j] * fl1_fx + 0.5 * pa2pb_z_y[j] * fl1_fx + pa2pb_yyz_y[j]);

                t_yyz_z[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_yy[j] * fl1_fx + 0.5 * pa2pb_z_z[j] * fl1_fx + pa2pb_yyz_z[j]);

                t_yzz_x[j] = fl_s_0_0 * (0.5 * pa2pb_y_x[j] * fl1_fx + pa2pb_yzz_x[j]);

                t_yzz_y[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_zz[j] * fl1_fx + 0.5 * pa2pb_y_y[j] * fl1_fx + pa2pb_yzz_y[j]);

                t_yzz_z[j] = fl_s_0_0 * (pa_yz[j] * fl1_fx + 0.5 * pa2pb_y_z[j] * fl1_fx + pa2pb_yzz_z[j]);

                t_zzz_x[j] = fl_s_0_0 * (1.5 * pa2pb_z_x[j] * fl1_fx + pa2pb_zzz_x[j]);

                t_zzz_y[j] = fl_s_0_0 * (1.5 * pa2pb_z_y[j] * fl1_fx + pa2pb_zzz_y[j]);

                t_zzz_z[j] = fl_s_0_0 * (0.75 * fl2_fx + 1.5 * pa_zz[j] * fl1_fx + 1.5 * pa2pb_z_z[j] * fl1_fx + pa2pb_zzz_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForPG(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForPG_0_9(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForPG_9_18(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForPG_18_27(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForPG_27_36(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForPG_36_45(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForPG_0_9(      CMemBlock2D<double>& primBuffer,
                         const CMemBlock2D<double>& auxBuffer,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CMemBlock2D<double>& pbDistances,
                         const CMemBlock2D<double>& pa2pbDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
    {
        // Batch of Integrals (0,9)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(102 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(102 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(102 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(102 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(102 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(102 * idx + 8);

            auto pa2pb_x_xxxx = pa2pbDistances.data(102 * idx + 19);

            auto pa2pb_x_xxxy = pa2pbDistances.data(102 * idx + 20);

            auto pa2pb_x_xxxz = pa2pbDistances.data(102 * idx + 21);

            auto pa2pb_x_xxyy = pa2pbDistances.data(102 * idx + 22);

            auto pa2pb_x_xxyz = pa2pbDistances.data(102 * idx + 23);

            auto pa2pb_x_xxzz = pa2pbDistances.data(102 * idx + 24);

            auto pa2pb_x_xyyy = pa2pbDistances.data(102 * idx + 25);

            auto pa2pb_x_xyyz = pa2pbDistances.data(102 * idx + 26);

            auto pa2pb_x_xyzz = pa2pbDistances.data(102 * idx + 27);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_xxxx = primBuffer.data(45 * idx);

            auto t_x_xxxy = primBuffer.data(45 * idx + 1);

            auto t_x_xxxz = primBuffer.data(45 * idx + 2);

            auto t_x_xxyy = primBuffer.data(45 * idx + 3);

            auto t_x_xxyz = primBuffer.data(45 * idx + 4);

            auto t_x_xxzz = primBuffer.data(45 * idx + 5);

            auto t_x_xyyy = primBuffer.data(45 * idx + 6);

            auto t_x_xyyz = primBuffer.data(45 * idx + 7);

            auto t_x_xyzz = primBuffer.data(45 * idx + 8);

            // Batch of Integrals (0,9)

            #pragma omp simd aligned(fx, pa2pb_x_xx, pa2pb_x_xxxx, pa2pb_x_xxxy, pa2pb_x_xxxz, pa2pb_x_xxyy, \
                                     pa2pb_x_xxyz, pa2pb_x_xxzz, pa2pb_x_xy, pa2pb_x_xyyy, pa2pb_x_xyyz, pa2pb_x_xyzz, \
                                     pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa_x, pb_x, pb_xxx, pb_xxy, pb_xxz, \
                                     pb_xyy, pb_xyz, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, s_0_0, t_x_xxxx, t_x_xxxy, \
                                     t_x_xxxz, t_x_xxyy, t_x_xxyz, t_x_xxzz, t_x_xyyy, t_x_xyyz, t_x_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_x_xxxx[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 3.0 * pb_x[j] * fl2_fx + 3.0 * pa2pb_x_xx[j] * fl1_fx + 2.0 * pb_xxx[j] * fl1_fx + pa2pb_x_xxxx[j]);

                t_x_xxxy[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 1.5 * pa2pb_x_xy[j] * fl1_fx + 1.5 * pb_xxy[j] * fl1_fx + pa2pb_x_xxxy[j]);

                t_x_xxxz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 1.5 * pa2pb_x_xz[j] * fl1_fx + 1.5 * pb_xxz[j] * fl1_fx + pa2pb_x_xxxz[j]);

                t_x_xxyy[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pb_x[j] * fl2_fx + 0.5 * pa2pb_x_xx[j] * fl1_fx + 0.5 * pa2pb_x_yy[j] * fl1_fx + pb_xyy[j] * fl1_fx + pa2pb_x_xxyy[j]);

                t_x_xxyz[j] = fl_s_0_0 * (0.5 * pa2pb_x_yz[j] * fl1_fx + pb_xyz[j] * fl1_fx + pa2pb_x_xxyz[j]);

                t_x_xxzz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pb_x[j] * fl2_fx + 0.5 * pa2pb_x_xx[j] * fl1_fx + 0.5 * pa2pb_x_zz[j] * fl1_fx + pb_xzz[j] * fl1_fx + pa2pb_x_xxzz[j]);

                t_x_xyyy[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 1.5 * pa2pb_x_xy[j] * fl1_fx + 0.5 * pb_yyy[j] * fl1_fx + pa2pb_x_xyyy[j]);

                t_x_xyyz[j] = fl_s_0_0 * (0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_x_xz[j] * fl1_fx + 0.5 * pb_yyz[j] * fl1_fx + pa2pb_x_xyyz[j]);

                t_x_xyzz[j] = fl_s_0_0 * (0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_x_xy[j] * fl1_fx + 0.5 * pb_yzz[j] * fl1_fx + pa2pb_x_xyzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForPG_9_18(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
    {
        // Batch of Integrals (9,18)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xz = pa2pbDistances.data(102 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(102 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(102 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(102 * idx + 8);

            auto pa2pb_x_xzzz = pa2pbDistances.data(102 * idx + 28);

            auto pa2pb_x_yyyy = pa2pbDistances.data(102 * idx + 29);

            auto pa2pb_x_yyyz = pa2pbDistances.data(102 * idx + 30);

            auto pa2pb_x_yyzz = pa2pbDistances.data(102 * idx + 31);

            auto pa2pb_x_yzzz = pa2pbDistances.data(102 * idx + 32);

            auto pa2pb_x_zzzz = pa2pbDistances.data(102 * idx + 33);

            auto pa2pb_y_xx = pa2pbDistances.data(102 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(102 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(102 * idx + 39);

            auto pa2pb_y_xxxx = pa2pbDistances.data(102 * idx + 53);

            auto pa2pb_y_xxxy = pa2pbDistances.data(102 * idx + 54);

            auto pa2pb_y_xxxz = pa2pbDistances.data(102 * idx + 55);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_xzzz = primBuffer.data(45 * idx + 9);

            auto t_x_yyyy = primBuffer.data(45 * idx + 10);

            auto t_x_yyyz = primBuffer.data(45 * idx + 11);

            auto t_x_yyzz = primBuffer.data(45 * idx + 12);

            auto t_x_yzzz = primBuffer.data(45 * idx + 13);

            auto t_x_zzzz = primBuffer.data(45 * idx + 14);

            auto t_y_xxxx = primBuffer.data(45 * idx + 15);

            auto t_y_xxxy = primBuffer.data(45 * idx + 16);

            auto t_y_xxxz = primBuffer.data(45 * idx + 17);

            // Batch of Integrals (9,18)

            #pragma omp simd aligned(fx, pa2pb_x_xz, pa2pb_x_xzzz, pa2pb_x_yy, pa2pb_x_yyyy, pa2pb_x_yyyz, \
                                     pa2pb_x_yyzz, pa2pb_x_yz, pa2pb_x_yzzz, pa2pb_x_zz, pa2pb_x_zzzz, pa2pb_y_xx, \
                                     pa2pb_y_xxxx, pa2pb_y_xxxy, pa2pb_y_xxxz, pa2pb_y_xy, pa2pb_y_xz, pa_x, pa_y, pb_x, \
                                     pb_xxx, pb_z, pb_zzz, s_0_0, t_x_xzzz, t_x_yyyy, t_x_yyyz, t_x_yyzz, t_x_yzzz, \
                                     t_x_zzzz, t_y_xxxx, t_y_xxxy, t_y_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_x_xzzz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 1.5 * pa2pb_x_xz[j] * fl1_fx + 0.5 * pb_zzz[j] * fl1_fx + pa2pb_x_xzzz[j]);

                t_x_yyyy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 3.0 * pa2pb_x_yy[j] * fl1_fx + pa2pb_x_yyyy[j]);

                t_x_yyyz[j] = fl_s_0_0 * (1.5 * pa2pb_x_yz[j] * fl1_fx + pa2pb_x_yyyz[j]);

                t_x_yyzz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa2pb_x_yy[j] * fl1_fx + 0.5 * pa2pb_x_zz[j] * fl1_fx + pa2pb_x_yyzz[j]);

                t_x_yzzz[j] = fl_s_0_0 * (1.5 * pa2pb_x_yz[j] * fl1_fx + pa2pb_x_yzzz[j]);

                t_x_zzzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 3.0 * pa2pb_x_zz[j] * fl1_fx + pa2pb_x_zzzz[j]);

                t_y_xxxx[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 3.0 * pa2pb_y_xx[j] * fl1_fx + pa2pb_y_xxxx[j]);

                t_y_xxxy[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 1.5 * pa2pb_y_xy[j] * fl1_fx + 0.5 * pb_xxx[j] * fl1_fx + pa2pb_y_xxxy[j]);

                t_y_xxxz[j] = fl_s_0_0 * (1.5 * pa2pb_y_xz[j] * fl1_fx + pa2pb_y_xxxz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForPG_18_27(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (18,27)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(102 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(102 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(102 * idx + 39);

            auto pa2pb_y_yy = pa2pbDistances.data(102 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(102 * idx + 41);

            auto pa2pb_y_zz = pa2pbDistances.data(102 * idx + 42);

            auto pa2pb_y_xxyy = pa2pbDistances.data(102 * idx + 56);

            auto pa2pb_y_xxyz = pa2pbDistances.data(102 * idx + 57);

            auto pa2pb_y_xxzz = pa2pbDistances.data(102 * idx + 58);

            auto pa2pb_y_xyyy = pa2pbDistances.data(102 * idx + 59);

            auto pa2pb_y_xyyz = pa2pbDistances.data(102 * idx + 60);

            auto pa2pb_y_xyzz = pa2pbDistances.data(102 * idx + 61);

            auto pa2pb_y_xzzz = pa2pbDistances.data(102 * idx + 62);

            auto pa2pb_y_yyyy = pa2pbDistances.data(102 * idx + 63);

            auto pa2pb_y_yyyz = pa2pbDistances.data(102 * idx + 64);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_y_xxyy = primBuffer.data(45 * idx + 18);

            auto t_y_xxyz = primBuffer.data(45 * idx + 19);

            auto t_y_xxzz = primBuffer.data(45 * idx + 20);

            auto t_y_xyyy = primBuffer.data(45 * idx + 21);

            auto t_y_xyyz = primBuffer.data(45 * idx + 22);

            auto t_y_xyzz = primBuffer.data(45 * idx + 23);

            auto t_y_xzzz = primBuffer.data(45 * idx + 24);

            auto t_y_yyyy = primBuffer.data(45 * idx + 25);

            auto t_y_yyyz = primBuffer.data(45 * idx + 26);

            // Batch of Integrals (18,27)

            #pragma omp simd aligned(fx, pa2pb_y_xx, pa2pb_y_xxyy, pa2pb_y_xxyz, pa2pb_y_xxzz, pa2pb_y_xy, \
                                     pa2pb_y_xyyy, pa2pb_y_xyyz, pa2pb_y_xyzz, pa2pb_y_xz, pa2pb_y_xzzz, pa2pb_y_yy, \
                                     pa2pb_y_yyyy, pa2pb_y_yyyz, pa2pb_y_yz, pa2pb_y_zz, pa_y, pb_x, pb_xxy, pb_xxz, pb_xyy, \
                                     pb_xyz, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_z, s_0_0, t_y_xxyy, t_y_xxyz, t_y_xxzz, \
                                     t_y_xyyy, t_y_xyyz, t_y_xyzz, t_y_xzzz, t_y_yyyy, t_y_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_y_xxyy[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pb_y[j] * fl2_fx + 0.5 * pa2pb_y_xx[j] * fl1_fx + 0.5 * pa2pb_y_yy[j] * fl1_fx + pb_xxy[j] * fl1_fx + pa2pb_y_xxyy[j]);

                t_y_xxyz[j] = fl_s_0_0 * (0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_y_yz[j] * fl1_fx + 0.5 * pb_xxz[j] * fl1_fx + pa2pb_y_xxyz[j]);

                t_y_xxzz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa2pb_y_xx[j] * fl1_fx + 0.5 * pa2pb_y_zz[j] * fl1_fx + pa2pb_y_xxzz[j]);

                t_y_xyyy[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 1.5 * pa2pb_y_xy[j] * fl1_fx + 1.5 * pb_xyy[j] * fl1_fx + pa2pb_y_xyyy[j]);

                t_y_xyyz[j] = fl_s_0_0 * (0.5 * pa2pb_y_xz[j] * fl1_fx + pb_xyz[j] * fl1_fx + pa2pb_y_xyyz[j]);

                t_y_xyzz[j] = fl_s_0_0 * (0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_y_xy[j] * fl1_fx + 0.5 * pb_xzz[j] * fl1_fx + pa2pb_y_xyzz[j]);

                t_y_xzzz[j] = fl_s_0_0 * (1.5 * pa2pb_y_xz[j] * fl1_fx + pa2pb_y_xzzz[j]);

                t_y_yyyy[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 3.0 * pb_y[j] * fl2_fx + 3.0 * pa2pb_y_yy[j] * fl1_fx + 2.0 * pb_yyy[j] * fl1_fx + pa2pb_y_yyyy[j]);

                t_y_yyyz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 1.5 * pa2pb_y_yz[j] * fl1_fx + 1.5 * pb_yyz[j] * fl1_fx + pa2pb_y_yyyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForPG_27_36(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (27,36)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_yy = pa2pbDistances.data(102 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(102 * idx + 41);

            auto pa2pb_y_zz = pa2pbDistances.data(102 * idx + 42);

            auto pa2pb_y_yyzz = pa2pbDistances.data(102 * idx + 65);

            auto pa2pb_y_yzzz = pa2pbDistances.data(102 * idx + 66);

            auto pa2pb_y_zzzz = pa2pbDistances.data(102 * idx + 67);

            auto pa2pb_z_xx = pa2pbDistances.data(102 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(102 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(102 * idx + 73);

            auto pa2pb_z_yy = pa2pbDistances.data(102 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(102 * idx + 75);

            auto pa2pb_z_zz = pa2pbDistances.data(102 * idx + 76);

            auto pa2pb_z_xxxx = pa2pbDistances.data(102 * idx + 87);

            auto pa2pb_z_xxxy = pa2pbDistances.data(102 * idx + 88);

            auto pa2pb_z_xxxz = pa2pbDistances.data(102 * idx + 89);

            auto pa2pb_z_xxyy = pa2pbDistances.data(102 * idx + 90);

            auto pa2pb_z_xxyz = pa2pbDistances.data(102 * idx + 91);

            auto pa2pb_z_xxzz = pa2pbDistances.data(102 * idx + 92);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_y_yyzz = primBuffer.data(45 * idx + 27);

            auto t_y_yzzz = primBuffer.data(45 * idx + 28);

            auto t_y_zzzz = primBuffer.data(45 * idx + 29);

            auto t_z_xxxx = primBuffer.data(45 * idx + 30);

            auto t_z_xxxy = primBuffer.data(45 * idx + 31);

            auto t_z_xxxz = primBuffer.data(45 * idx + 32);

            auto t_z_xxyy = primBuffer.data(45 * idx + 33);

            auto t_z_xxyz = primBuffer.data(45 * idx + 34);

            auto t_z_xxzz = primBuffer.data(45 * idx + 35);

            // Batch of Integrals (27,36)

            #pragma omp simd aligned(fx, pa2pb_y_yy, pa2pb_y_yyzz, pa2pb_y_yz, pa2pb_y_yzzz, pa2pb_y_zz, \
                                     pa2pb_y_zzzz, pa2pb_z_xx, pa2pb_z_xxxx, pa2pb_z_xxxy, pa2pb_z_xxxz, pa2pb_z_xxyy, \
                                     pa2pb_z_xxyz, pa2pb_z_xxzz, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, \
                                     pa2pb_z_zz, pa_y, pa_z, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_y, pb_yzz, pb_z, pb_zzz, s_0_0, \
                                     t_y_yyzz, t_y_yzzz, t_y_zzzz, t_z_xxxx, t_z_xxxy, t_z_xxxz, t_z_xxyy, t_z_xxyz, \
                                     t_z_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_y_yyzz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pb_y[j] * fl2_fx + 0.5 * pa2pb_y_yy[j] * fl1_fx + 0.5 * pa2pb_y_zz[j] * fl1_fx + pb_yzz[j] * fl1_fx + pa2pb_y_yyzz[j]);

                t_y_yzzz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 1.5 * pa2pb_y_yz[j] * fl1_fx + 0.5 * pb_zzz[j] * fl1_fx + pa2pb_y_yzzz[j]);

                t_y_zzzz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 3.0 * pa2pb_y_zz[j] * fl1_fx + pa2pb_y_zzzz[j]);

                t_z_xxxx[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 3.0 * pa2pb_z_xx[j] * fl1_fx + pa2pb_z_xxxx[j]);

                t_z_xxxy[j] = fl_s_0_0 * (1.5 * pa2pb_z_xy[j] * fl1_fx + pa2pb_z_xxxy[j]);

                t_z_xxxz[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 1.5 * pa2pb_z_xz[j] * fl1_fx + 0.5 * pb_xxx[j] * fl1_fx + pa2pb_z_xxxz[j]);

                t_z_xxyy[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa2pb_z_xx[j] * fl1_fx + 0.5 * pa2pb_z_yy[j] * fl1_fx + pa2pb_z_xxyy[j]);

                t_z_xxyz[j] = fl_s_0_0 * (0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_z_yz[j] * fl1_fx + 0.5 * pb_xxy[j] * fl1_fx + pa2pb_z_xxyz[j]);

                t_z_xxzz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pb_z[j] * fl2_fx + 0.5 * pa2pb_z_xx[j] * fl1_fx + 0.5 * pa2pb_z_zz[j] * fl1_fx + pb_xxz[j] * fl1_fx + pa2pb_z_xxzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForPG_36_45(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (36,45)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_xy = pa2pbDistances.data(102 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(102 * idx + 73);

            auto pa2pb_z_yy = pa2pbDistances.data(102 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(102 * idx + 75);

            auto pa2pb_z_zz = pa2pbDistances.data(102 * idx + 76);

            auto pa2pb_z_xyyy = pa2pbDistances.data(102 * idx + 93);

            auto pa2pb_z_xyyz = pa2pbDistances.data(102 * idx + 94);

            auto pa2pb_z_xyzz = pa2pbDistances.data(102 * idx + 95);

            auto pa2pb_z_xzzz = pa2pbDistances.data(102 * idx + 96);

            auto pa2pb_z_yyyy = pa2pbDistances.data(102 * idx + 97);

            auto pa2pb_z_yyyz = pa2pbDistances.data(102 * idx + 98);

            auto pa2pb_z_yyzz = pa2pbDistances.data(102 * idx + 99);

            auto pa2pb_z_yzzz = pa2pbDistances.data(102 * idx + 100);

            auto pa2pb_z_zzzz = pa2pbDistances.data(102 * idx + 101);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_z_xyyy = primBuffer.data(45 * idx + 36);

            auto t_z_xyyz = primBuffer.data(45 * idx + 37);

            auto t_z_xyzz = primBuffer.data(45 * idx + 38);

            auto t_z_xzzz = primBuffer.data(45 * idx + 39);

            auto t_z_yyyy = primBuffer.data(45 * idx + 40);

            auto t_z_yyyz = primBuffer.data(45 * idx + 41);

            auto t_z_yyzz = primBuffer.data(45 * idx + 42);

            auto t_z_yzzz = primBuffer.data(45 * idx + 43);

            auto t_z_zzzz = primBuffer.data(45 * idx + 44);

            // Batch of Integrals (36,45)

            #pragma omp simd aligned(fx, pa2pb_z_xy, pa2pb_z_xyyy, pa2pb_z_xyyz, pa2pb_z_xyzz, pa2pb_z_xz, \
                                     pa2pb_z_xzzz, pa2pb_z_yy, pa2pb_z_yyyy, pa2pb_z_yyyz, pa2pb_z_yyzz, pa2pb_z_yz, \
                                     pa2pb_z_yzzz, pa2pb_z_zz, pa2pb_z_zzzz, pa_z, pb_x, pb_xyy, pb_xyz, pb_xzz, pb_y, pb_yyy, \
                                     pb_yyz, pb_yzz, pb_z, pb_zzz, s_0_0, t_z_xyyy, t_z_xyyz, t_z_xyzz, t_z_xzzz, \
                                     t_z_yyyy, t_z_yyyz, t_z_yyzz, t_z_yzzz, t_z_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_z_xyyy[j] = fl_s_0_0 * (1.5 * pa2pb_z_xy[j] * fl1_fx + pa2pb_z_xyyy[j]);

                t_z_xyyz[j] = fl_s_0_0 * (0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_z_xz[j] * fl1_fx + 0.5 * pb_xyy[j] * fl1_fx + pa2pb_z_xyyz[j]);

                t_z_xyzz[j] = fl_s_0_0 * (0.5 * pa2pb_z_xy[j] * fl1_fx + pb_xyz[j] * fl1_fx + pa2pb_z_xyzz[j]);

                t_z_xzzz[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 1.5 * pa2pb_z_xz[j] * fl1_fx + 1.5 * pb_xzz[j] * fl1_fx + pa2pb_z_xzzz[j]);

                t_z_yyyy[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 3.0 * pa2pb_z_yy[j] * fl1_fx + pa2pb_z_yyyy[j]);

                t_z_yyyz[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 1.5 * pa2pb_z_yz[j] * fl1_fx + 0.5 * pb_yyy[j] * fl1_fx + pa2pb_z_yyyz[j]);

                t_z_yyzz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pb_z[j] * fl2_fx + 0.5 * pa2pb_z_yy[j] * fl1_fx + 0.5 * pa2pb_z_zz[j] * fl1_fx + pb_yyz[j] * fl1_fx + pa2pb_z_yyzz[j]);

                t_z_yzzz[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 1.5 * pa2pb_z_yz[j] * fl1_fx + 1.5 * pb_yzz[j] * fl1_fx + pa2pb_z_yzzz[j]);

                t_z_zzzz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 3.0 * pb_z[j] * fl2_fx + 3.0 * pa2pb_z_zz[j] * fl1_fx + 2.0 * pb_zzz[j] * fl1_fx + pa2pb_z_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForGP_0_9(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                         braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGP_9_18(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGP_18_27(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGP_27_36(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGP_36_45(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForGP_0_9(      CMemBlock2D<double>& primBuffer,
                         const CMemBlock2D<double>& auxBuffer,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CMemBlock2D<double>& pbDistances,
                         const CMemBlock2D<double>& pa2pbDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
    {
        // Batch of Integrals (0,9)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_xx_x = pa2pbDistances.data(102 * idx + 9);

            auto pa2pb_xx_y = pa2pbDistances.data(102 * idx + 10);

            auto pa2pb_xx_z = pa2pbDistances.data(102 * idx + 11);

            auto pa2pb_xy_x = pa2pbDistances.data(102 * idx + 12);

            auto pa2pb_xy_y = pa2pbDistances.data(102 * idx + 13);

            auto pa2pb_xy_z = pa2pbDistances.data(102 * idx + 14);

            auto pa2pb_xz_x = pa2pbDistances.data(102 * idx + 15);

            auto pa2pb_xz_y = pa2pbDistances.data(102 * idx + 16);

            auto pa2pb_xz_z = pa2pbDistances.data(102 * idx + 17);

            auto pa2pb_xxxx_x = pa2pbDistances.data(102 * idx + 57);

            auto pa2pb_xxxx_y = pa2pbDistances.data(102 * idx + 58);

            auto pa2pb_xxxx_z = pa2pbDistances.data(102 * idx + 59);

            auto pa2pb_xxxy_x = pa2pbDistances.data(102 * idx + 60);

            auto pa2pb_xxxy_y = pa2pbDistances.data(102 * idx + 61);

            auto pa2pb_xxxy_z = pa2pbDistances.data(102 * idx + 62);

            auto pa2pb_xxxz_x = pa2pbDistances.data(102 * idx + 63);

            auto pa2pb_xxxz_y = pa2pbDistances.data(102 * idx + 64);

            auto pa2pb_xxxz_z = pa2pbDistances.data(102 * idx + 65);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxxx_x = primBuffer.data(45 * idx);

            auto t_xxxx_y = primBuffer.data(45 * idx + 1);

            auto t_xxxx_z = primBuffer.data(45 * idx + 2);

            auto t_xxxy_x = primBuffer.data(45 * idx + 3);

            auto t_xxxy_y = primBuffer.data(45 * idx + 4);

            auto t_xxxy_z = primBuffer.data(45 * idx + 5);

            auto t_xxxz_x = primBuffer.data(45 * idx + 6);

            auto t_xxxz_y = primBuffer.data(45 * idx + 7);

            auto t_xxxz_z = primBuffer.data(45 * idx + 8);

            // Batch of Integrals (0,9)

            #pragma omp simd aligned(fx, pa2pb_xx_x, pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxxx_x, pa2pb_xxxx_y, \
                                     pa2pb_xxxx_z, pa2pb_xxxy_x, pa2pb_xxxy_y, pa2pb_xxxy_z, pa2pb_xxxz_x, pa2pb_xxxz_y, \
                                     pa2pb_xxxz_z, pa2pb_xy_x, pa2pb_xy_y, pa2pb_xy_z, pa2pb_xz_x, pa2pb_xz_y, pa2pb_xz_z, \
                                     pa_x, pa_xxx, pa_xxy, pa_xxz, pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xxxx_x, t_xxxx_y, \
                                     t_xxxx_z, t_xxxy_x, t_xxxy_y, t_xxxy_z, t_xxxz_x, t_xxxz_y, t_xxxz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxxx_x[j] = fl_s_0_0 * (3.0 * pa_x[j] * fl2_fx + 2.0 * pa_xxx[j] * fl1_fx + 0.75 * pb_x[j] * fl2_fx + 3.0 * pa2pb_xx_x[j] * fl1_fx + pa2pb_xxxx_x[j]);

                t_xxxx_y[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 3.0 * pa2pb_xx_y[j] * fl1_fx + pa2pb_xxxx_y[j]);

                t_xxxx_z[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 3.0 * pa2pb_xx_z[j] * fl1_fx + pa2pb_xxxx_z[j]);

                t_xxxy_x[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pa_xxy[j] * fl1_fx + 1.5 * pa2pb_xy_x[j] * fl1_fx + pa2pb_xxxy_x[j]);

                t_xxxy_y[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 0.5 * pa_xxx[j] * fl1_fx + 1.5 * pa2pb_xy_y[j] * fl1_fx + pa2pb_xxxy_y[j]);

                t_xxxy_z[j] = fl_s_0_0 * (1.5 * pa2pb_xy_z[j] * fl1_fx + pa2pb_xxxy_z[j]);

                t_xxxz_x[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 1.5 * pa_xxz[j] * fl1_fx + 1.5 * pa2pb_xz_x[j] * fl1_fx + pa2pb_xxxz_x[j]);

                t_xxxz_y[j] = fl_s_0_0 * (1.5 * pa2pb_xz_y[j] * fl1_fx + pa2pb_xxxz_y[j]);

                t_xxxz_z[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 0.5 * pa_xxx[j] * fl1_fx + 1.5 * pa2pb_xz_z[j] * fl1_fx + pa2pb_xxxz_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGP_9_18(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
    {
        // Batch of Integrals (9,18)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_xx_x = pa2pbDistances.data(102 * idx + 9);

            auto pa2pb_xx_y = pa2pbDistances.data(102 * idx + 10);

            auto pa2pb_xx_z = pa2pbDistances.data(102 * idx + 11);

            auto pa2pb_yy_x = pa2pbDistances.data(102 * idx + 18);

            auto pa2pb_yy_y = pa2pbDistances.data(102 * idx + 19);

            auto pa2pb_yy_z = pa2pbDistances.data(102 * idx + 20);

            auto pa2pb_yz_x = pa2pbDistances.data(102 * idx + 21);

            auto pa2pb_yz_y = pa2pbDistances.data(102 * idx + 22);

            auto pa2pb_yz_z = pa2pbDistances.data(102 * idx + 23);

            auto pa2pb_zz_x = pa2pbDistances.data(102 * idx + 24);

            auto pa2pb_zz_y = pa2pbDistances.data(102 * idx + 25);

            auto pa2pb_zz_z = pa2pbDistances.data(102 * idx + 26);

            auto pa2pb_xxyy_x = pa2pbDistances.data(102 * idx + 66);

            auto pa2pb_xxyy_y = pa2pbDistances.data(102 * idx + 67);

            auto pa2pb_xxyy_z = pa2pbDistances.data(102 * idx + 68);

            auto pa2pb_xxyz_x = pa2pbDistances.data(102 * idx + 69);

            auto pa2pb_xxyz_y = pa2pbDistances.data(102 * idx + 70);

            auto pa2pb_xxyz_z = pa2pbDistances.data(102 * idx + 71);

            auto pa2pb_xxzz_x = pa2pbDistances.data(102 * idx + 72);

            auto pa2pb_xxzz_y = pa2pbDistances.data(102 * idx + 73);

            auto pa2pb_xxzz_z = pa2pbDistances.data(102 * idx + 74);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxyy_x = primBuffer.data(45 * idx + 9);

            auto t_xxyy_y = primBuffer.data(45 * idx + 10);

            auto t_xxyy_z = primBuffer.data(45 * idx + 11);

            auto t_xxyz_x = primBuffer.data(45 * idx + 12);

            auto t_xxyz_y = primBuffer.data(45 * idx + 13);

            auto t_xxyz_z = primBuffer.data(45 * idx + 14);

            auto t_xxzz_x = primBuffer.data(45 * idx + 15);

            auto t_xxzz_y = primBuffer.data(45 * idx + 16);

            auto t_xxzz_z = primBuffer.data(45 * idx + 17);

            // Batch of Integrals (9,18)

            #pragma omp simd aligned(fx, pa2pb_xx_x, pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxyy_x, pa2pb_xxyy_y, \
                                     pa2pb_xxyy_z, pa2pb_xxyz_x, pa2pb_xxyz_y, pa2pb_xxyz_z, pa2pb_xxzz_x, pa2pb_xxzz_y, \
                                     pa2pb_xxzz_z, pa2pb_yy_x, pa2pb_yy_y, pa2pb_yy_z, pa2pb_yz_x, pa2pb_yz_y, pa2pb_yz_z, \
                                     pa2pb_zz_x, pa2pb_zz_y, pa2pb_zz_z, pa_x, pa_xxy, pa_xxz, pa_xyy, pa_xyz, pa_xzz, pa_y, \
                                     pa_z, pb_x, pb_y, pb_z, s_0_0, t_xxyy_x, t_xxyy_y, t_xxyy_z, t_xxyz_x, t_xxyz_y, \
                                     t_xxyz_z, t_xxzz_x, t_xxzz_y, t_xxzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxyy_x[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl2_fx + pa_xyy[j] * fl1_fx + 0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_xx_x[j] * fl1_fx + 0.5 * pa2pb_yy_x[j] * fl1_fx + pa2pb_xxyy_x[j]);

                t_xxyy_y[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl2_fx + pa_xxy[j] * fl1_fx + 0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_xx_y[j] * fl1_fx + 0.5 * pa2pb_yy_y[j] * fl1_fx + pa2pb_xxyy_y[j]);

                t_xxyy_z[j] = fl_s_0_0 * (0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_xx_z[j] * fl1_fx + 0.5 * pa2pb_yy_z[j] * fl1_fx + pa2pb_xxyy_z[j]);

                t_xxyz_x[j] = fl_s_0_0 * (pa_xyz[j] * fl1_fx + 0.5 * pa2pb_yz_x[j] * fl1_fx + pa2pb_xxyz_x[j]);

                t_xxyz_y[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa_xxz[j] * fl1_fx + 0.5 * pa2pb_yz_y[j] * fl1_fx + pa2pb_xxyz_y[j]);

                t_xxyz_z[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa_xxy[j] * fl1_fx + 0.5 * pa2pb_yz_z[j] * fl1_fx + pa2pb_xxyz_z[j]);

                t_xxzz_x[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl2_fx + pa_xzz[j] * fl1_fx + 0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_xx_x[j] * fl1_fx + 0.5 * pa2pb_zz_x[j] * fl1_fx + pa2pb_xxzz_x[j]);

                t_xxzz_y[j] = fl_s_0_0 * (0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_xx_y[j] * fl1_fx + 0.5 * pa2pb_zz_y[j] * fl1_fx + pa2pb_xxzz_y[j]);

                t_xxzz_z[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl2_fx + pa_xxz[j] * fl1_fx + 0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_xx_z[j] * fl1_fx + 0.5 * pa2pb_zz_z[j] * fl1_fx + pa2pb_xxzz_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGP_18_27(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (18,27)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_xy_x = pa2pbDistances.data(102 * idx + 12);

            auto pa2pb_xy_y = pa2pbDistances.data(102 * idx + 13);

            auto pa2pb_xy_z = pa2pbDistances.data(102 * idx + 14);

            auto pa2pb_xz_x = pa2pbDistances.data(102 * idx + 15);

            auto pa2pb_xz_y = pa2pbDistances.data(102 * idx + 16);

            auto pa2pb_xz_z = pa2pbDistances.data(102 * idx + 17);

            auto pa2pb_xyyy_x = pa2pbDistances.data(102 * idx + 75);

            auto pa2pb_xyyy_y = pa2pbDistances.data(102 * idx + 76);

            auto pa2pb_xyyy_z = pa2pbDistances.data(102 * idx + 77);

            auto pa2pb_xyyz_x = pa2pbDistances.data(102 * idx + 78);

            auto pa2pb_xyyz_y = pa2pbDistances.data(102 * idx + 79);

            auto pa2pb_xyyz_z = pa2pbDistances.data(102 * idx + 80);

            auto pa2pb_xyzz_x = pa2pbDistances.data(102 * idx + 81);

            auto pa2pb_xyzz_y = pa2pbDistances.data(102 * idx + 82);

            auto pa2pb_xyzz_z = pa2pbDistances.data(102 * idx + 83);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xyyy_x = primBuffer.data(45 * idx + 18);

            auto t_xyyy_y = primBuffer.data(45 * idx + 19);

            auto t_xyyy_z = primBuffer.data(45 * idx + 20);

            auto t_xyyz_x = primBuffer.data(45 * idx + 21);

            auto t_xyyz_y = primBuffer.data(45 * idx + 22);

            auto t_xyyz_z = primBuffer.data(45 * idx + 23);

            auto t_xyzz_x = primBuffer.data(45 * idx + 24);

            auto t_xyzz_y = primBuffer.data(45 * idx + 25);

            auto t_xyzz_z = primBuffer.data(45 * idx + 26);

            // Batch of Integrals (18,27)

            #pragma omp simd aligned(fx, pa2pb_xy_x, pa2pb_xy_y, pa2pb_xy_z, pa2pb_xyyy_x, pa2pb_xyyy_y, \
                                     pa2pb_xyyy_z, pa2pb_xyyz_x, pa2pb_xyyz_y, pa2pb_xyyz_z, pa2pb_xyzz_x, pa2pb_xyzz_y, \
                                     pa2pb_xyzz_z, pa2pb_xz_x, pa2pb_xz_y, pa2pb_xz_z, pa_x, pa_xyy, pa_xyz, pa_xzz, pa_y, \
                                     pa_yyy, pa_yyz, pa_yzz, pa_z, s_0_0, t_xyyy_x, t_xyyy_y, t_xyyy_z, t_xyyz_x, \
                                     t_xyyz_y, t_xyyz_z, t_xyzz_x, t_xyzz_y, t_xyzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xyyy_x[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 0.5 * pa_yyy[j] * fl1_fx + 1.5 * pa2pb_xy_x[j] * fl1_fx + pa2pb_xyyy_x[j]);

                t_xyyy_y[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pa_xyy[j] * fl1_fx + 1.5 * pa2pb_xy_y[j] * fl1_fx + pa2pb_xyyy_y[j]);

                t_xyyy_z[j] = fl_s_0_0 * (1.5 * pa2pb_xy_z[j] * fl1_fx + pa2pb_xyyy_z[j]);

                t_xyyz_x[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa_yyz[j] * fl1_fx + 0.5 * pa2pb_xz_x[j] * fl1_fx + pa2pb_xyyz_x[j]);

                t_xyyz_y[j] = fl_s_0_0 * (pa_xyz[j] * fl1_fx + 0.5 * pa2pb_xz_y[j] * fl1_fx + pa2pb_xyyz_y[j]);

                t_xyyz_z[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa_xyy[j] * fl1_fx + 0.5 * pa2pb_xz_z[j] * fl1_fx + pa2pb_xyyz_z[j]);

                t_xyzz_x[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa_yzz[j] * fl1_fx + 0.5 * pa2pb_xy_x[j] * fl1_fx + pa2pb_xyzz_x[j]);

                t_xyzz_y[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa_xzz[j] * fl1_fx + 0.5 * pa2pb_xy_y[j] * fl1_fx + pa2pb_xyzz_y[j]);

                t_xyzz_z[j] = fl_s_0_0 * (pa_xyz[j] * fl1_fx + 0.5 * pa2pb_xy_z[j] * fl1_fx + pa2pb_xyzz_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGP_27_36(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (27,36)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_xz_x = pa2pbDistances.data(102 * idx + 15);

            auto pa2pb_xz_y = pa2pbDistances.data(102 * idx + 16);

            auto pa2pb_xz_z = pa2pbDistances.data(102 * idx + 17);

            auto pa2pb_yy_x = pa2pbDistances.data(102 * idx + 18);

            auto pa2pb_yy_y = pa2pbDistances.data(102 * idx + 19);

            auto pa2pb_yy_z = pa2pbDistances.data(102 * idx + 20);

            auto pa2pb_yz_x = pa2pbDistances.data(102 * idx + 21);

            auto pa2pb_yz_y = pa2pbDistances.data(102 * idx + 22);

            auto pa2pb_yz_z = pa2pbDistances.data(102 * idx + 23);

            auto pa2pb_xzzz_x = pa2pbDistances.data(102 * idx + 84);

            auto pa2pb_xzzz_y = pa2pbDistances.data(102 * idx + 85);

            auto pa2pb_xzzz_z = pa2pbDistances.data(102 * idx + 86);

            auto pa2pb_yyyy_x = pa2pbDistances.data(102 * idx + 87);

            auto pa2pb_yyyy_y = pa2pbDistances.data(102 * idx + 88);

            auto pa2pb_yyyy_z = pa2pbDistances.data(102 * idx + 89);

            auto pa2pb_yyyz_x = pa2pbDistances.data(102 * idx + 90);

            auto pa2pb_yyyz_y = pa2pbDistances.data(102 * idx + 91);

            auto pa2pb_yyyz_z = pa2pbDistances.data(102 * idx + 92);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xzzz_x = primBuffer.data(45 * idx + 27);

            auto t_xzzz_y = primBuffer.data(45 * idx + 28);

            auto t_xzzz_z = primBuffer.data(45 * idx + 29);

            auto t_yyyy_x = primBuffer.data(45 * idx + 30);

            auto t_yyyy_y = primBuffer.data(45 * idx + 31);

            auto t_yyyy_z = primBuffer.data(45 * idx + 32);

            auto t_yyyz_x = primBuffer.data(45 * idx + 33);

            auto t_yyyz_y = primBuffer.data(45 * idx + 34);

            auto t_yyyz_z = primBuffer.data(45 * idx + 35);

            // Batch of Integrals (27,36)

            #pragma omp simd aligned(fx, pa2pb_xz_x, pa2pb_xz_y, pa2pb_xz_z, pa2pb_xzzz_x, pa2pb_xzzz_y, \
                                     pa2pb_xzzz_z, pa2pb_yy_x, pa2pb_yy_y, pa2pb_yy_z, pa2pb_yyyy_x, pa2pb_yyyy_y, \
                                     pa2pb_yyyy_z, pa2pb_yyyz_x, pa2pb_yyyz_y, pa2pb_yyyz_z, pa2pb_yz_x, pa2pb_yz_y, \
                                     pa2pb_yz_z, pa_x, pa_xzz, pa_y, pa_yyy, pa_yyz, pa_z, pa_zzz, pb_x, pb_y, pb_z, s_0_0, \
                                     t_xzzz_x, t_xzzz_y, t_xzzz_z, t_yyyy_x, t_yyyy_y, t_yyyy_z, t_yyyz_x, t_yyyz_y, \
                                     t_yyyz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xzzz_x[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 0.5 * pa_zzz[j] * fl1_fx + 1.5 * pa2pb_xz_x[j] * fl1_fx + pa2pb_xzzz_x[j]);

                t_xzzz_y[j] = fl_s_0_0 * (1.5 * pa2pb_xz_y[j] * fl1_fx + pa2pb_xzzz_y[j]);

                t_xzzz_z[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pa_xzz[j] * fl1_fx + 1.5 * pa2pb_xz_z[j] * fl1_fx + pa2pb_xzzz_z[j]);

                t_yyyy_x[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 3.0 * pa2pb_yy_x[j] * fl1_fx + pa2pb_yyyy_x[j]);

                t_yyyy_y[j] = fl_s_0_0 * (3.0 * pa_y[j] * fl2_fx + 2.0 * pa_yyy[j] * fl1_fx + 0.75 * pb_y[j] * fl2_fx + 3.0 * pa2pb_yy_y[j] * fl1_fx + pa2pb_yyyy_y[j]);

                t_yyyy_z[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 3.0 * pa2pb_yy_z[j] * fl1_fx + pa2pb_yyyy_z[j]);

                t_yyyz_x[j] = fl_s_0_0 * (1.5 * pa2pb_yz_x[j] * fl1_fx + pa2pb_yyyz_x[j]);

                t_yyyz_y[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 1.5 * pa_yyz[j] * fl1_fx + 1.5 * pa2pb_yz_y[j] * fl1_fx + pa2pb_yyyz_y[j]);

                t_yyyz_z[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 0.5 * pa_yyy[j] * fl1_fx + 1.5 * pa2pb_yz_z[j] * fl1_fx + pa2pb_yyyz_z[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGP_36_45(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (36,45)

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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_yy_x = pa2pbDistances.data(102 * idx + 18);

            auto pa2pb_yy_y = pa2pbDistances.data(102 * idx + 19);

            auto pa2pb_yy_z = pa2pbDistances.data(102 * idx + 20);

            auto pa2pb_yz_x = pa2pbDistances.data(102 * idx + 21);

            auto pa2pb_yz_y = pa2pbDistances.data(102 * idx + 22);

            auto pa2pb_yz_z = pa2pbDistances.data(102 * idx + 23);

            auto pa2pb_zz_x = pa2pbDistances.data(102 * idx + 24);

            auto pa2pb_zz_y = pa2pbDistances.data(102 * idx + 25);

            auto pa2pb_zz_z = pa2pbDistances.data(102 * idx + 26);

            auto pa2pb_yyzz_x = pa2pbDistances.data(102 * idx + 93);

            auto pa2pb_yyzz_y = pa2pbDistances.data(102 * idx + 94);

            auto pa2pb_yyzz_z = pa2pbDistances.data(102 * idx + 95);

            auto pa2pb_yzzz_x = pa2pbDistances.data(102 * idx + 96);

            auto pa2pb_yzzz_y = pa2pbDistances.data(102 * idx + 97);

            auto pa2pb_yzzz_z = pa2pbDistances.data(102 * idx + 98);

            auto pa2pb_zzzz_x = pa2pbDistances.data(102 * idx + 99);

            auto pa2pb_zzzz_y = pa2pbDistances.data(102 * idx + 100);

            auto pa2pb_zzzz_z = pa2pbDistances.data(102 * idx + 101);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yyzz_x = primBuffer.data(45 * idx + 36);

            auto t_yyzz_y = primBuffer.data(45 * idx + 37);

            auto t_yyzz_z = primBuffer.data(45 * idx + 38);

            auto t_yzzz_x = primBuffer.data(45 * idx + 39);

            auto t_yzzz_y = primBuffer.data(45 * idx + 40);

            auto t_yzzz_z = primBuffer.data(45 * idx + 41);

            auto t_zzzz_x = primBuffer.data(45 * idx + 42);

            auto t_zzzz_y = primBuffer.data(45 * idx + 43);

            auto t_zzzz_z = primBuffer.data(45 * idx + 44);

            // Batch of Integrals (36,45)

            #pragma omp simd aligned(fx, pa2pb_yy_x, pa2pb_yy_y, pa2pb_yy_z, pa2pb_yyzz_x, pa2pb_yyzz_y, \
                                     pa2pb_yyzz_z, pa2pb_yz_x, pa2pb_yz_y, pa2pb_yz_z, pa2pb_yzzz_x, pa2pb_yzzz_y, \
                                     pa2pb_yzzz_z, pa2pb_zz_x, pa2pb_zz_y, pa2pb_zz_z, pa2pb_zzzz_x, pa2pb_zzzz_y, \
                                     pa2pb_zzzz_z, pa_y, pa_yyz, pa_yzz, pa_z, pa_zzz, pb_x, pb_y, pb_z, s_0_0, t_yyzz_x, t_yyzz_y, \
                                     t_yyzz_z, t_yzzz_x, t_yzzz_y, t_yzzz_z, t_zzzz_x, t_zzzz_y, t_zzzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yyzz_x[j] = fl_s_0_0 * (0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_yy_x[j] * fl1_fx + 0.5 * pa2pb_zz_x[j] * fl1_fx + pa2pb_yyzz_x[j]);

                t_yyzz_y[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl2_fx + pa_yzz[j] * fl1_fx + 0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_yy_y[j] * fl1_fx + 0.5 * pa2pb_zz_y[j] * fl1_fx + pa2pb_yyzz_y[j]);

                t_yyzz_z[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl2_fx + pa_yyz[j] * fl1_fx + 0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_yy_z[j] * fl1_fx + 0.5 * pa2pb_zz_z[j] * fl1_fx + pa2pb_yyzz_z[j]);

                t_yzzz_x[j] = fl_s_0_0 * (1.5 * pa2pb_yz_x[j] * fl1_fx + pa2pb_yzzz_x[j]);

                t_yzzz_y[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 0.5 * pa_zzz[j] * fl1_fx + 1.5 * pa2pb_yz_y[j] * fl1_fx + pa2pb_yzzz_y[j]);

                t_yzzz_z[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pa_yzz[j] * fl1_fx + 1.5 * pa2pb_yz_z[j] * fl1_fx + pa2pb_yzzz_z[j]);

                t_zzzz_x[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 3.0 * pa2pb_zz_x[j] * fl1_fx + pa2pb_zzzz_x[j]);

                t_zzzz_y[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 3.0 * pa2pb_zz_y[j] * fl1_fx + pa2pb_zzzz_y[j]);

                t_zzzz_z[j] = fl_s_0_0 * (3.0 * pa_z[j] * fl2_fx + 2.0 * pa_zzz[j] * fl1_fx + 0.75 * pb_z[j] * fl2_fx + 3.0 * pa2pb_zz_z[j] * fl1_fx + pa2pb_zzzz_z[j]);
            }

            idx++;
        }
    }


} // ovlrecfunc namespace

