//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFuncForDX.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace

    void
    compOverlapForDD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForDD_0_12(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDD_12_24(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDD_24_36(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForDD_0_12(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
    {
        // Batch of Integrals (0,12)

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

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(81 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(81 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(81 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(81 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(81 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(81 * idx + 11);

            auto pa2pb_xx_xx = pa2pbDistances.data(81 * idx + 30);

            auto pa2pb_xx_xy = pa2pbDistances.data(81 * idx + 31);

            auto pa2pb_xx_xz = pa2pbDistances.data(81 * idx + 32);

            auto pa2pb_xx_yy = pa2pbDistances.data(81 * idx + 33);

            auto pa2pb_xx_yz = pa2pbDistances.data(81 * idx + 34);

            auto pa2pb_xx_zz = pa2pbDistances.data(81 * idx + 35);

            auto pa2pb_xy_xx = pa2pbDistances.data(81 * idx + 39);

            auto pa2pb_xy_xy = pa2pbDistances.data(81 * idx + 40);

            auto pa2pb_xy_xz = pa2pbDistances.data(81 * idx + 41);

            auto pa2pb_xy_yy = pa2pbDistances.data(81 * idx + 42);

            auto pa2pb_xy_yz = pa2pbDistances.data(81 * idx + 43);

            auto pa2pb_xy_zz = pa2pbDistances.data(81 * idx + 44);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_xx = primBuffer.data(36 * idx);

            auto t_xx_xy = primBuffer.data(36 * idx + 1);

            auto t_xx_xz = primBuffer.data(36 * idx + 2);

            auto t_xx_yy = primBuffer.data(36 * idx + 3);

            auto t_xx_yz = primBuffer.data(36 * idx + 4);

            auto t_xx_zz = primBuffer.data(36 * idx + 5);

            auto t_xy_xx = primBuffer.data(36 * idx + 6);

            auto t_xy_xy = primBuffer.data(36 * idx + 7);

            auto t_xy_xz = primBuffer.data(36 * idx + 8);

            auto t_xy_yy = primBuffer.data(36 * idx + 9);

            auto t_xy_yz = primBuffer.data(36 * idx + 10);

            auto t_xy_zz = primBuffer.data(36 * idx + 11);

            // Batch of Integrals (0,12)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xy, \
                                     pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa2pb_xy_xx, pa2pb_xy_xy, \
                                     pa2pb_xy_xz, pa2pb_xy_yy, pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, \
                                     pa_xx, pa_xy, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, t_xx_xx, t_xx_xy, \
                                     t_xx_xz, t_xx_yy, t_xx_yz, t_xx_zz, t_xy_xx, t_xy_xy, t_xy_xz, t_xy_yy, t_xy_yz, \
                                     t_xy_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xx_xx[j] = fl_s_0_0 * (0.75 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 2.0 * pa2pb_x_x[j] * fl1_fx + 0.5 * pb_xx[j] * fl1_fx + pa2pb_xx_xx[j]);

                t_xx_xy[j] = fl_s_0_0 * (pa2pb_x_y[j] * fl1_fx + 0.5 * pb_xy[j] * fl1_fx + pa2pb_xx_xy[j]);

                t_xx_xz[j] = fl_s_0_0 * (pa2pb_x_z[j] * fl1_fx + 0.5 * pb_xz[j] * fl1_fx + pa2pb_xx_xz[j]);

                t_xx_yy[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * pb_yy[j] * fl1_fx + pa2pb_xx_yy[j]);

                t_xx_yz[j] = fl_s_0_0 * (0.5 * pb_yz[j] * fl1_fx + pa2pb_xx_yz[j]);

                t_xx_zz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * pb_zz[j] * fl1_fx + pa2pb_xx_zz[j]);

                t_xy_xx[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl1_fx + pa2pb_y_x[j] * fl1_fx + pa2pb_xy_xx[j]);

                t_xy_xy[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_x_x[j] * fl1_fx + 0.5 * pa2pb_y_y[j] * fl1_fx + pa2pb_xy_xy[j]);

                t_xy_xz[j] = fl_s_0_0 * (0.5 * pa2pb_y_z[j] * fl1_fx + pa2pb_xy_xz[j]);

                t_xy_yy[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl1_fx + pa2pb_x_y[j] * fl1_fx + pa2pb_xy_yy[j]);

                t_xy_yz[j] = fl_s_0_0 * (0.5 * pa2pb_x_z[j] * fl1_fx + pa2pb_xy_yz[j]);

                t_xy_zz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl1_fx + pa2pb_xy_zz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDD_12_24(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (12,24)

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

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(81 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(81 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(81 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(81 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(81 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(81 * idx + 11);

            auto pa2pb_z_x = pa2pbDistances.data(81 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(81 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(81 * idx + 20);

            auto pa2pb_xz_xx = pa2pbDistances.data(81 * idx + 48);

            auto pa2pb_xz_xy = pa2pbDistances.data(81 * idx + 49);

            auto pa2pb_xz_xz = pa2pbDistances.data(81 * idx + 50);

            auto pa2pb_xz_yy = pa2pbDistances.data(81 * idx + 51);

            auto pa2pb_xz_yz = pa2pbDistances.data(81 * idx + 52);

            auto pa2pb_xz_zz = pa2pbDistances.data(81 * idx + 53);

            auto pa2pb_yy_xx = pa2pbDistances.data(81 * idx + 57);

            auto pa2pb_yy_xy = pa2pbDistances.data(81 * idx + 58);

            auto pa2pb_yy_xz = pa2pbDistances.data(81 * idx + 59);

            auto pa2pb_yy_yy = pa2pbDistances.data(81 * idx + 60);

            auto pa2pb_yy_yz = pa2pbDistances.data(81 * idx + 61);

            auto pa2pb_yy_zz = pa2pbDistances.data(81 * idx + 62);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xz_xx = primBuffer.data(36 * idx + 12);

            auto t_xz_xy = primBuffer.data(36 * idx + 13);

            auto t_xz_xz = primBuffer.data(36 * idx + 14);

            auto t_xz_yy = primBuffer.data(36 * idx + 15);

            auto t_xz_yz = primBuffer.data(36 * idx + 16);

            auto t_xz_zz = primBuffer.data(36 * idx + 17);

            auto t_yy_xx = primBuffer.data(36 * idx + 18);

            auto t_yy_xy = primBuffer.data(36 * idx + 19);

            auto t_yy_xz = primBuffer.data(36 * idx + 20);

            auto t_yy_yy = primBuffer.data(36 * idx + 21);

            auto t_yy_yz = primBuffer.data(36 * idx + 22);

            auto t_yy_zz = primBuffer.data(36 * idx + 23);

            // Batch of Integrals (12,24)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xz_xx, pa2pb_xz_xy, \
                                     pa2pb_xz_xz, pa2pb_xz_yy, pa2pb_xz_yz, pa2pb_xz_zz, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, \
                                     pa2pb_yy_xx, pa2pb_yy_xy, pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, pa2pb_yy_zz, \
                                     pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, pa_xz, pa_yy, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, \
                                     s_0_0, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy, t_xz_yz, t_xz_zz, t_yy_xx, t_yy_xy, \
                                     t_yy_xz, t_yy_yy, t_yy_yz, t_yy_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xz_xx[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl1_fx + pa2pb_z_x[j] * fl1_fx + pa2pb_xz_xx[j]);

                t_xz_xy[j] = fl_s_0_0 * (0.5 * pa2pb_z_y[j] * fl1_fx + pa2pb_xz_xy[j]);

                t_xz_xz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_x_x[j] * fl1_fx + 0.5 * pa2pb_z_z[j] * fl1_fx + pa2pb_xz_xz[j]);

                t_xz_yy[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl1_fx + pa2pb_xz_yy[j]);

                t_xz_yz[j] = fl_s_0_0 * (0.5 * pa2pb_x_y[j] * fl1_fx + pa2pb_xz_yz[j]);

                t_xz_zz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl1_fx + pa2pb_x_z[j] * fl1_fx + pa2pb_xz_zz[j]);

                t_yy_xx[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_yy[j] * fl1_fx + 0.5 * pb_xx[j] * fl1_fx + pa2pb_yy_xx[j]);

                t_yy_xy[j] = fl_s_0_0 * (pa2pb_y_x[j] * fl1_fx + 0.5 * pb_xy[j] * fl1_fx + pa2pb_yy_xy[j]);

                t_yy_xz[j] = fl_s_0_0 * (0.5 * pb_xz[j] * fl1_fx + pa2pb_yy_xz[j]);

                t_yy_yy[j] = fl_s_0_0 * (0.75 * fl2_fx + 0.5 * pa_yy[j] * fl1_fx + 2.0 * pa2pb_y_y[j] * fl1_fx + 0.5 * pb_yy[j] * fl1_fx + pa2pb_yy_yy[j]);

                t_yy_yz[j] = fl_s_0_0 * (pa2pb_y_z[j] * fl1_fx + 0.5 * pb_yz[j] * fl1_fx + pa2pb_yy_yz[j]);

                t_yy_zz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_yy[j] * fl1_fx + 0.5 * pb_zz[j] * fl1_fx + pa2pb_yy_zz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDD_24_36(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (24,36)

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

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(81 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(81 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(81 * idx + 11);

            auto pa2pb_z_x = pa2pbDistances.data(81 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(81 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(81 * idx + 20);

            auto pa2pb_yz_xx = pa2pbDistances.data(81 * idx + 66);

            auto pa2pb_yz_xy = pa2pbDistances.data(81 * idx + 67);

            auto pa2pb_yz_xz = pa2pbDistances.data(81 * idx + 68);

            auto pa2pb_yz_yy = pa2pbDistances.data(81 * idx + 69);

            auto pa2pb_yz_yz = pa2pbDistances.data(81 * idx + 70);

            auto pa2pb_yz_zz = pa2pbDistances.data(81 * idx + 71);

            auto pa2pb_zz_xx = pa2pbDistances.data(81 * idx + 75);

            auto pa2pb_zz_xy = pa2pbDistances.data(81 * idx + 76);

            auto pa2pb_zz_xz = pa2pbDistances.data(81 * idx + 77);

            auto pa2pb_zz_yy = pa2pbDistances.data(81 * idx + 78);

            auto pa2pb_zz_yz = pa2pbDistances.data(81 * idx + 79);

            auto pa2pb_zz_zz = pa2pbDistances.data(81 * idx + 80);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yz_xx = primBuffer.data(36 * idx + 24);

            auto t_yz_xy = primBuffer.data(36 * idx + 25);

            auto t_yz_xz = primBuffer.data(36 * idx + 26);

            auto t_yz_yy = primBuffer.data(36 * idx + 27);

            auto t_yz_yz = primBuffer.data(36 * idx + 28);

            auto t_yz_zz = primBuffer.data(36 * idx + 29);

            auto t_zz_xx = primBuffer.data(36 * idx + 30);

            auto t_zz_xy = primBuffer.data(36 * idx + 31);

            auto t_zz_xz = primBuffer.data(36 * idx + 32);

            auto t_zz_yy = primBuffer.data(36 * idx + 33);

            auto t_zz_yz = primBuffer.data(36 * idx + 34);

            auto t_zz_zz = primBuffer.data(36 * idx + 35);

            // Batch of Integrals (24,36)

            #pragma omp simd aligned(fx, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa2pb_yz_xx, pa2pb_yz_xy, \
                                     pa2pb_yz_xz, pa2pb_yz_yy, pa2pb_yz_yz, pa2pb_yz_zz, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, \
                                     pa2pb_zz_xx, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zz_zz, pa_yz, \
                                     pa_zz, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, t_yz_xx, t_yz_xy, t_yz_xz, \
                                     t_yz_yy, t_yz_yz, t_yz_zz, t_zz_xx, t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz, t_zz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yz_xx[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl1_fx + pa2pb_yz_xx[j]);

                t_yz_xy[j] = fl_s_0_0 * (0.5 * pa2pb_z_x[j] * fl1_fx + pa2pb_yz_xy[j]);

                t_yz_xz[j] = fl_s_0_0 * (0.5 * pa2pb_y_x[j] * fl1_fx + pa2pb_yz_xz[j]);

                t_yz_yy[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl1_fx + pa2pb_z_y[j] * fl1_fx + pa2pb_yz_yy[j]);

                t_yz_yz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa2pb_y_y[j] * fl1_fx + 0.5 * pa2pb_z_z[j] * fl1_fx + pa2pb_yz_yz[j]);

                t_yz_zz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl1_fx + pa2pb_y_z[j] * fl1_fx + pa2pb_yz_zz[j]);

                t_zz_xx[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_zz[j] * fl1_fx + 0.5 * pb_xx[j] * fl1_fx + pa2pb_zz_xx[j]);

                t_zz_xy[j] = fl_s_0_0 * (0.5 * pb_xy[j] * fl1_fx + pa2pb_zz_xy[j]);

                t_zz_xz[j] = fl_s_0_0 * (pa2pb_z_x[j] * fl1_fx + 0.5 * pb_xz[j] * fl1_fx + pa2pb_zz_xz[j]);

                t_zz_yy[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_zz[j] * fl1_fx + 0.5 * pb_yy[j] * fl1_fx + pa2pb_zz_yy[j]);

                t_zz_yz[j] = fl_s_0_0 * (pa2pb_z_y[j] * fl1_fx + 0.5 * pb_yz[j] * fl1_fx + pa2pb_zz_yz[j]);

                t_zz_zz[j] = fl_s_0_0 * (0.75 * fl2_fx + 0.5 * pa_zz[j] * fl1_fx + 2.0 * pa2pb_z_z[j] * fl1_fx + 0.5 * pb_zz[j] * fl1_fx + pa2pb_zz_zz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDF(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForDF_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDF_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDF_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDF_30_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDF_40_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDF_50_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForDF_0_10(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(9 * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

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

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(171 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(171 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(171 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(171 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(171 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(171 * idx + 8);

            auto pa2pb_xx_x = pa2pbDistances.data(171 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(171 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(171 * idx + 59);

            auto pa2pb_xx_xxx = pa2pbDistances.data(171 * idx + 66);

            auto pa2pb_xx_xxy = pa2pbDistances.data(171 * idx + 67);

            auto pa2pb_xx_xxz = pa2pbDistances.data(171 * idx + 68);

            auto pa2pb_xx_xyy = pa2pbDistances.data(171 * idx + 69);

            auto pa2pb_xx_xyz = pa2pbDistances.data(171 * idx + 70);

            auto pa2pb_xx_xzz = pa2pbDistances.data(171 * idx + 71);

            auto pa2pb_xx_yyy = pa2pbDistances.data(171 * idx + 72);

            auto pa2pb_xx_yyz = pa2pbDistances.data(171 * idx + 73);

            auto pa2pb_xx_yzz = pa2pbDistances.data(171 * idx + 74);

            auto pa2pb_xx_zzz = pa2pbDistances.data(171 * idx + 75);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_xxx = primBuffer.data(60 * idx);

            auto t_xx_xxy = primBuffer.data(60 * idx + 1);

            auto t_xx_xxz = primBuffer.data(60 * idx + 2);

            auto t_xx_xyy = primBuffer.data(60 * idx + 3);

            auto t_xx_xyz = primBuffer.data(60 * idx + 4);

            auto t_xx_xzz = primBuffer.data(60 * idx + 5);

            auto t_xx_yyy = primBuffer.data(60 * idx + 6);

            auto t_xx_yyz = primBuffer.data(60 * idx + 7);

            auto t_xx_yzz = primBuffer.data(60 * idx + 8);

            auto t_xx_zzz = primBuffer.data(60 * idx + 9);

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_x_zz, pa2pb_xx_x, pa2pb_xx_xxx, pa2pb_xx_xxy, pa2pb_xx_xxz, pa2pb_xx_xyy, \
                                     pa2pb_xx_xyz, pa2pb_xx_xzz, pa2pb_xx_y, pa2pb_xx_yyy, pa2pb_xx_yyz, pa2pb_xx_yzz, \
                                     pa2pb_xx_z, pa2pb_xx_zzz, pa_x, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_xzz, pb_y, \
                                     pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, s_0_0, t_xx_xxx, t_xx_xxy, t_xx_xxz, t_xx_xyy, \
                                     t_xx_xyz, t_xx_xzz, t_xx_yyy, t_xx_yyz, t_xx_yzz, t_xx_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xx_xxx[j] = fl_s_0_0 * (1.5 * pa_x[j] * fl2_fx + 2.25 * pb_x[j] * fl2_fx + 1.5 * pa2pb_xx_x[j] * fl1_fx + 3.0 * pa2pb_x_xx[j] * fl1_fx + 0.5 * pb_xxx[j] * fl1_fx + pa2pb_xx_xxx[j]);

                t_xx_xxy[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 0.5 * pa2pb_xx_y[j] * fl1_fx + 2.0 * pa2pb_x_xy[j] * fl1_fx + 0.5 * pb_xxy[j] * fl1_fx + pa2pb_xx_xxy[j]);

                t_xx_xxz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 0.5 * pa2pb_xx_z[j] * fl1_fx + 2.0 * pa2pb_x_xz[j] * fl1_fx + 0.5 * pb_xxz[j] * fl1_fx + pa2pb_xx_xxz[j]);

                t_xx_xyy[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl2_fx + 0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_xx_x[j] * fl1_fx + pa2pb_x_yy[j] * fl1_fx + 0.5 * pb_xyy[j] * fl1_fx + pa2pb_xx_xyy[j]);

                t_xx_xyz[j] = fl_s_0_0 * (pa2pb_x_yz[j] * fl1_fx + 0.5 * pb_xyz[j] * fl1_fx + pa2pb_xx_xyz[j]);

                t_xx_xzz[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl2_fx + 0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_xx_x[j] * fl1_fx + pa2pb_x_zz[j] * fl1_fx + 0.5 * pb_xzz[j] * fl1_fx + pa2pb_xx_xzz[j]);

                t_xx_yyy[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 1.5 * pa2pb_xx_y[j] * fl1_fx + 0.5 * pb_yyy[j] * fl1_fx + pa2pb_xx_yyy[j]);

                t_xx_yyz[j] = fl_s_0_0 * (0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_xx_z[j] * fl1_fx + 0.5 * pb_yyz[j] * fl1_fx + pa2pb_xx_yyz[j]);

                t_xx_yzz[j] = fl_s_0_0 * (0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_xx_y[j] * fl1_fx + 0.5 * pb_yzz[j] * fl1_fx + pa2pb_xx_yzz[j]);

                t_xx_zzz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 1.5 * pa2pb_xx_z[j] * fl1_fx + 0.5 * pb_zzz[j] * fl1_fx + pa2pb_xx_zzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDF_10_20(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (10,20)

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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(171 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(171 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(171 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(171 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(171 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(171 * idx + 8);

            auto pa2pb_y_xx = pa2pbDistances.data(171 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(171 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(171 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(171 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(171 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(171 * idx + 27);

            auto pa2pb_xy_x = pa2pbDistances.data(171 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(171 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(171 * idx + 78);

            auto pa2pb_xy_xxx = pa2pbDistances.data(171 * idx + 85);

            auto pa2pb_xy_xxy = pa2pbDistances.data(171 * idx + 86);

            auto pa2pb_xy_xxz = pa2pbDistances.data(171 * idx + 87);

            auto pa2pb_xy_xyy = pa2pbDistances.data(171 * idx + 88);

            auto pa2pb_xy_xyz = pa2pbDistances.data(171 * idx + 89);

            auto pa2pb_xy_xzz = pa2pbDistances.data(171 * idx + 90);

            auto pa2pb_xy_yyy = pa2pbDistances.data(171 * idx + 91);

            auto pa2pb_xy_yyz = pa2pbDistances.data(171 * idx + 92);

            auto pa2pb_xy_yzz = pa2pbDistances.data(171 * idx + 93);

            auto pa2pb_xy_zzz = pa2pbDistances.data(171 * idx + 94);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xy_xxx = primBuffer.data(60 * idx + 10);

            auto t_xy_xxy = primBuffer.data(60 * idx + 11);

            auto t_xy_xxz = primBuffer.data(60 * idx + 12);

            auto t_xy_xyy = primBuffer.data(60 * idx + 13);

            auto t_xy_xyz = primBuffer.data(60 * idx + 14);

            auto t_xy_xzz = primBuffer.data(60 * idx + 15);

            auto t_xy_yyy = primBuffer.data(60 * idx + 16);

            auto t_xy_yyz = primBuffer.data(60 * idx + 17);

            auto t_xy_yzz = primBuffer.data(60 * idx + 18);

            auto t_xy_zzz = primBuffer.data(60 * idx + 19);

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fx, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_x_zz, pa2pb_xy_x, pa2pb_xy_xxx, pa2pb_xy_xxy, pa2pb_xy_xxz, pa2pb_xy_xyy, \
                                     pa2pb_xy_xyz, pa2pb_xy_xzz, pa2pb_xy_y, pa2pb_xy_yyy, pa2pb_xy_yyz, pa2pb_xy_yzz, \
                                     pa2pb_xy_z, pa2pb_xy_zzz, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, \
                                     pa2pb_y_yz, pa2pb_y_zz, pa_x, pa_y, pb_x, pb_y, pb_z, s_0_0, t_xy_xxx, t_xy_xxy, t_xy_xxz, \
                                     t_xy_xyy, t_xy_xyz, t_xy_xzz, t_xy_yyy, t_xy_yyz, t_xy_yzz, t_xy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xy_xxx[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pa2pb_xy_x[j] * fl1_fx + 1.5 * pa2pb_y_xx[j] * fl1_fx + pa2pb_xy_xxx[j]);

                t_xy_xxy[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pb_x[j] * fl2_fx + 0.5 * pa2pb_xy_y[j] * fl1_fx + 0.5 * pa2pb_x_xx[j] * fl1_fx + pa2pb_y_xy[j] * fl1_fx + pa2pb_xy_xxy[j]);

                t_xy_xxz[j] = fl_s_0_0 * (0.5 * pa2pb_xy_z[j] * fl1_fx + pa2pb_y_xz[j] * fl1_fx + pa2pb_xy_xxz[j]);

                t_xy_xyy[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pb_y[j] * fl2_fx + 0.5 * pa2pb_xy_x[j] * fl1_fx + pa2pb_x_xy[j] * fl1_fx + 0.5 * pa2pb_y_yy[j] * fl1_fx + pa2pb_xy_xyy[j]);

                t_xy_xyz[j] = fl_s_0_0 * (0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_x_xz[j] * fl1_fx + 0.5 * pa2pb_y_yz[j] * fl1_fx + pa2pb_xy_xyz[j]);

                t_xy_xzz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa2pb_xy_x[j] * fl1_fx + 0.5 * pa2pb_y_zz[j] * fl1_fx + pa2pb_xy_xzz[j]);

                t_xy_yyy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pa2pb_xy_y[j] * fl1_fx + 1.5 * pa2pb_x_yy[j] * fl1_fx + pa2pb_xy_yyy[j]);

                t_xy_yyz[j] = fl_s_0_0 * (0.5 * pa2pb_xy_z[j] * fl1_fx + pa2pb_x_yz[j] * fl1_fx + pa2pb_xy_yyz[j]);

                t_xy_yzz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa2pb_xy_y[j] * fl1_fx + 0.5 * pa2pb_x_zz[j] * fl1_fx + pa2pb_xy_yzz[j]);

                t_xy_zzz[j] = fl_s_0_0 * (1.5 * pa2pb_xy_z[j] * fl1_fx + pa2pb_xy_zzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDF_20_30(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (20,30)

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

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(171 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(171 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(171 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(171 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(171 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(171 * idx + 8);

            auto pa2pb_z_xx = pa2pbDistances.data(171 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(171 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(171 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(171 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(171 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(171 * idx + 46);

            auto pa2pb_xz_x = pa2pbDistances.data(171 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(171 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(171 * idx + 97);

            auto pa2pb_xz_xxx = pa2pbDistances.data(171 * idx + 104);

            auto pa2pb_xz_xxy = pa2pbDistances.data(171 * idx + 105);

            auto pa2pb_xz_xxz = pa2pbDistances.data(171 * idx + 106);

            auto pa2pb_xz_xyy = pa2pbDistances.data(171 * idx + 107);

            auto pa2pb_xz_xyz = pa2pbDistances.data(171 * idx + 108);

            auto pa2pb_xz_xzz = pa2pbDistances.data(171 * idx + 109);

            auto pa2pb_xz_yyy = pa2pbDistances.data(171 * idx + 110);

            auto pa2pb_xz_yyz = pa2pbDistances.data(171 * idx + 111);

            auto pa2pb_xz_yzz = pa2pbDistances.data(171 * idx + 112);

            auto pa2pb_xz_zzz = pa2pbDistances.data(171 * idx + 113);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xz_xxx = primBuffer.data(60 * idx + 20);

            auto t_xz_xxy = primBuffer.data(60 * idx + 21);

            auto t_xz_xxz = primBuffer.data(60 * idx + 22);

            auto t_xz_xyy = primBuffer.data(60 * idx + 23);

            auto t_xz_xyz = primBuffer.data(60 * idx + 24);

            auto t_xz_xzz = primBuffer.data(60 * idx + 25);

            auto t_xz_yyy = primBuffer.data(60 * idx + 26);

            auto t_xz_yyz = primBuffer.data(60 * idx + 27);

            auto t_xz_yzz = primBuffer.data(60 * idx + 28);

            auto t_xz_zzz = primBuffer.data(60 * idx + 29);

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fx, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_x_zz, pa2pb_xz_x, pa2pb_xz_xxx, pa2pb_xz_xxy, pa2pb_xz_xxz, pa2pb_xz_xyy, \
                                     pa2pb_xz_xyz, pa2pb_xz_xzz, pa2pb_xz_y, pa2pb_xz_yyy, pa2pb_xz_yyz, pa2pb_xz_yzz, \
                                     pa2pb_xz_z, pa2pb_xz_zzz, pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, \
                                     pa2pb_z_yz, pa2pb_z_zz, pa_x, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xz_xxx, t_xz_xxy, t_xz_xxz, \
                                     t_xz_xyy, t_xz_xyz, t_xz_xzz, t_xz_yyy, t_xz_yyz, t_xz_yzz, t_xz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xz_xxx[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 1.5 * pa2pb_xz_x[j] * fl1_fx + 1.5 * pa2pb_z_xx[j] * fl1_fx + pa2pb_xz_xxx[j]);

                t_xz_xxy[j] = fl_s_0_0 * (0.5 * pa2pb_xz_y[j] * fl1_fx + pa2pb_z_xy[j] * fl1_fx + pa2pb_xz_xxy[j]);

                t_xz_xxz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pb_x[j] * fl2_fx + 0.5 * pa2pb_xz_z[j] * fl1_fx + 0.5 * pa2pb_x_xx[j] * fl1_fx + pa2pb_z_xz[j] * fl1_fx + pa2pb_xz_xxz[j]);

                t_xz_xyy[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa2pb_xz_x[j] * fl1_fx + 0.5 * pa2pb_z_yy[j] * fl1_fx + pa2pb_xz_xyy[j]);

                t_xz_xyz[j] = fl_s_0_0 * (0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_x_xy[j] * fl1_fx + 0.5 * pa2pb_z_yz[j] * fl1_fx + pa2pb_xz_xyz[j]);

                t_xz_xzz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pb_z[j] * fl2_fx + 0.5 * pa2pb_xz_x[j] * fl1_fx + pa2pb_x_xz[j] * fl1_fx + 0.5 * pa2pb_z_zz[j] * fl1_fx + pa2pb_xz_xzz[j]);

                t_xz_yyy[j] = fl_s_0_0 * (1.5 * pa2pb_xz_y[j] * fl1_fx + pa2pb_xz_yyy[j]);

                t_xz_yyz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa2pb_xz_z[j] * fl1_fx + 0.5 * pa2pb_x_yy[j] * fl1_fx + pa2pb_xz_yyz[j]);

                t_xz_yzz[j] = fl_s_0_0 * (0.5 * pa2pb_xz_y[j] * fl1_fx + pa2pb_x_yz[j] * fl1_fx + pa2pb_xz_yzz[j]);

                t_xz_zzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pa2pb_xz_z[j] * fl1_fx + 1.5 * pa2pb_x_zz[j] * fl1_fx + pa2pb_xz_zzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDF_30_40(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (30,40)

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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

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

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(171 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(171 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(171 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(171 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(171 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(171 * idx + 27);

            auto pa2pb_yy_x = pa2pbDistances.data(171 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(171 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(171 * idx + 116);

            auto pa2pb_yy_xxx = pa2pbDistances.data(171 * idx + 123);

            auto pa2pb_yy_xxy = pa2pbDistances.data(171 * idx + 124);

            auto pa2pb_yy_xxz = pa2pbDistances.data(171 * idx + 125);

            auto pa2pb_yy_xyy = pa2pbDistances.data(171 * idx + 126);

            auto pa2pb_yy_xyz = pa2pbDistances.data(171 * idx + 127);

            auto pa2pb_yy_xzz = pa2pbDistances.data(171 * idx + 128);

            auto pa2pb_yy_yyy = pa2pbDistances.data(171 * idx + 129);

            auto pa2pb_yy_yyz = pa2pbDistances.data(171 * idx + 130);

            auto pa2pb_yy_yzz = pa2pbDistances.data(171 * idx + 131);

            auto pa2pb_yy_zzz = pa2pbDistances.data(171 * idx + 132);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yy_xxx = primBuffer.data(60 * idx + 30);

            auto t_yy_xxy = primBuffer.data(60 * idx + 31);

            auto t_yy_xxz = primBuffer.data(60 * idx + 32);

            auto t_yy_xyy = primBuffer.data(60 * idx + 33);

            auto t_yy_xyz = primBuffer.data(60 * idx + 34);

            auto t_yy_xzz = primBuffer.data(60 * idx + 35);

            auto t_yy_yyy = primBuffer.data(60 * idx + 36);

            auto t_yy_yyz = primBuffer.data(60 * idx + 37);

            auto t_yy_yzz = primBuffer.data(60 * idx + 38);

            auto t_yy_zzz = primBuffer.data(60 * idx + 39);

            // Batch of Integrals (30,40)

            #pragma omp simd aligned(fx, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, \
                                     pa2pb_y_zz, pa2pb_yy_x, pa2pb_yy_xxx, pa2pb_yy_xxy, pa2pb_yy_xxz, pa2pb_yy_xyy, \
                                     pa2pb_yy_xyz, pa2pb_yy_xzz, pa2pb_yy_y, pa2pb_yy_yyy, pa2pb_yy_yyz, pa2pb_yy_yzz, \
                                     pa2pb_yy_z, pa2pb_yy_zzz, pa_y, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_xzz, pb_y, \
                                     pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, s_0_0, t_yy_xxx, t_yy_xxy, t_yy_xxz, t_yy_xyy, \
                                     t_yy_xyz, t_yy_xzz, t_yy_yyy, t_yy_yyz, t_yy_yzz, t_yy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yy_xxx[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 1.5 * pa2pb_yy_x[j] * fl1_fx + 0.5 * pb_xxx[j] * fl1_fx + pa2pb_yy_xxx[j]);

                t_yy_xxy[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl2_fx + 0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_yy_y[j] * fl1_fx + pa2pb_y_xx[j] * fl1_fx + 0.5 * pb_xxy[j] * fl1_fx + pa2pb_yy_xxy[j]);

                t_yy_xxz[j] = fl_s_0_0 * (0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_yy_z[j] * fl1_fx + 0.5 * pb_xxz[j] * fl1_fx + pa2pb_yy_xxz[j]);

                t_yy_xyy[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 0.5 * pa2pb_yy_x[j] * fl1_fx + 2.0 * pa2pb_y_xy[j] * fl1_fx + 0.5 * pb_xyy[j] * fl1_fx + pa2pb_yy_xyy[j]);

                t_yy_xyz[j] = fl_s_0_0 * (pa2pb_y_xz[j] * fl1_fx + 0.5 * pb_xyz[j] * fl1_fx + pa2pb_yy_xyz[j]);

                t_yy_xzz[j] = fl_s_0_0 * (0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_yy_x[j] * fl1_fx + 0.5 * pb_xzz[j] * fl1_fx + pa2pb_yy_xzz[j]);

                t_yy_yyy[j] = fl_s_0_0 * (1.5 * pa_y[j] * fl2_fx + 2.25 * pb_y[j] * fl2_fx + 1.5 * pa2pb_yy_y[j] * fl1_fx + 3.0 * pa2pb_y_yy[j] * fl1_fx + 0.5 * pb_yyy[j] * fl1_fx + pa2pb_yy_yyy[j]);

                t_yy_yyz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 0.5 * pa2pb_yy_z[j] * fl1_fx + 2.0 * pa2pb_y_yz[j] * fl1_fx + 0.5 * pb_yyz[j] * fl1_fx + pa2pb_yy_yyz[j]);

                t_yy_yzz[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl2_fx + 0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_yy_y[j] * fl1_fx + pa2pb_y_zz[j] * fl1_fx + 0.5 * pb_yzz[j] * fl1_fx + pa2pb_yy_yzz[j]);

                t_yy_zzz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 1.5 * pa2pb_yy_z[j] * fl1_fx + 0.5 * pb_zzz[j] * fl1_fx + pa2pb_yy_zzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDF_40_50(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (40,50)

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

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(171 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(171 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(171 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(171 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(171 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(171 * idx + 27);

            auto pa2pb_z_xx = pa2pbDistances.data(171 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(171 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(171 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(171 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(171 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(171 * idx + 46);

            auto pa2pb_yz_x = pa2pbDistances.data(171 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(171 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(171 * idx + 135);

            auto pa2pb_yz_xxx = pa2pbDistances.data(171 * idx + 142);

            auto pa2pb_yz_xxy = pa2pbDistances.data(171 * idx + 143);

            auto pa2pb_yz_xxz = pa2pbDistances.data(171 * idx + 144);

            auto pa2pb_yz_xyy = pa2pbDistances.data(171 * idx + 145);

            auto pa2pb_yz_xyz = pa2pbDistances.data(171 * idx + 146);

            auto pa2pb_yz_xzz = pa2pbDistances.data(171 * idx + 147);

            auto pa2pb_yz_yyy = pa2pbDistances.data(171 * idx + 148);

            auto pa2pb_yz_yyz = pa2pbDistances.data(171 * idx + 149);

            auto pa2pb_yz_yzz = pa2pbDistances.data(171 * idx + 150);

            auto pa2pb_yz_zzz = pa2pbDistances.data(171 * idx + 151);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yz_xxx = primBuffer.data(60 * idx + 40);

            auto t_yz_xxy = primBuffer.data(60 * idx + 41);

            auto t_yz_xxz = primBuffer.data(60 * idx + 42);

            auto t_yz_xyy = primBuffer.data(60 * idx + 43);

            auto t_yz_xyz = primBuffer.data(60 * idx + 44);

            auto t_yz_xzz = primBuffer.data(60 * idx + 45);

            auto t_yz_yyy = primBuffer.data(60 * idx + 46);

            auto t_yz_yyz = primBuffer.data(60 * idx + 47);

            auto t_yz_yzz = primBuffer.data(60 * idx + 48);

            auto t_yz_zzz = primBuffer.data(60 * idx + 49);

            // Batch of Integrals (40,50)

            #pragma omp simd aligned(fx, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, \
                                     pa2pb_y_zz, pa2pb_yz_x, pa2pb_yz_xxx, pa2pb_yz_xxy, pa2pb_yz_xxz, pa2pb_yz_xyy, \
                                     pa2pb_yz_xyz, pa2pb_yz_xzz, pa2pb_yz_y, pa2pb_yz_yyy, pa2pb_yz_yyz, pa2pb_yz_yzz, \
                                     pa2pb_yz_z, pa2pb_yz_zzz, pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, \
                                     pa2pb_z_yz, pa2pb_z_zz, pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, t_yz_xxx, t_yz_xxy, t_yz_xxz, \
                                     t_yz_xyy, t_yz_xyz, t_yz_xzz, t_yz_yyy, t_yz_yyz, t_yz_yzz, t_yz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yz_xxx[j] = fl_s_0_0 * (1.5 * pa2pb_yz_x[j] * fl1_fx + pa2pb_yz_xxx[j]);

                t_yz_xxy[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa2pb_yz_y[j] * fl1_fx + 0.5 * pa2pb_z_xx[j] * fl1_fx + pa2pb_yz_xxy[j]);

                t_yz_xxz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa2pb_yz_z[j] * fl1_fx + 0.5 * pa2pb_y_xx[j] * fl1_fx + pa2pb_yz_xxz[j]);

                t_yz_xyy[j] = fl_s_0_0 * (0.5 * pa2pb_yz_x[j] * fl1_fx + pa2pb_z_xy[j] * fl1_fx + pa2pb_yz_xyy[j]);

                t_yz_xyz[j] = fl_s_0_0 * (0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_y_xy[j] * fl1_fx + 0.5 * pa2pb_z_xz[j] * fl1_fx + pa2pb_yz_xyz[j]);

                t_yz_xzz[j] = fl_s_0_0 * (0.5 * pa2pb_yz_x[j] * fl1_fx + pa2pb_y_xz[j] * fl1_fx + pa2pb_yz_xzz[j]);

                t_yz_yyy[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 1.5 * pa2pb_yz_y[j] * fl1_fx + 1.5 * pa2pb_z_yy[j] * fl1_fx + pa2pb_yz_yyy[j]);

                t_yz_yyz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pb_y[j] * fl2_fx + 0.5 * pa2pb_yz_z[j] * fl1_fx + 0.5 * pa2pb_y_yy[j] * fl1_fx + pa2pb_z_yz[j] * fl1_fx + pa2pb_yz_yyz[j]);

                t_yz_yzz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pb_z[j] * fl2_fx + 0.5 * pa2pb_yz_y[j] * fl1_fx + pa2pb_y_yz[j] * fl1_fx + 0.5 * pa2pb_z_zz[j] * fl1_fx + pa2pb_yz_yzz[j]);

                t_yz_zzz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pa2pb_yz_z[j] * fl1_fx + 1.5 * pa2pb_y_zz[j] * fl1_fx + pa2pb_yz_zzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDF_50_60(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (50,60)

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

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

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

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_xx = pa2pbDistances.data(171 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(171 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(171 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(171 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(171 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(171 * idx + 46);

            auto pa2pb_zz_x = pa2pbDistances.data(171 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(171 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(171 * idx + 154);

            auto pa2pb_zz_xxx = pa2pbDistances.data(171 * idx + 161);

            auto pa2pb_zz_xxy = pa2pbDistances.data(171 * idx + 162);

            auto pa2pb_zz_xxz = pa2pbDistances.data(171 * idx + 163);

            auto pa2pb_zz_xyy = pa2pbDistances.data(171 * idx + 164);

            auto pa2pb_zz_xyz = pa2pbDistances.data(171 * idx + 165);

            auto pa2pb_zz_xzz = pa2pbDistances.data(171 * idx + 166);

            auto pa2pb_zz_yyy = pa2pbDistances.data(171 * idx + 167);

            auto pa2pb_zz_yyz = pa2pbDistances.data(171 * idx + 168);

            auto pa2pb_zz_yzz = pa2pbDistances.data(171 * idx + 169);

            auto pa2pb_zz_zzz = pa2pbDistances.data(171 * idx + 170);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_zz_xxx = primBuffer.data(60 * idx + 50);

            auto t_zz_xxy = primBuffer.data(60 * idx + 51);

            auto t_zz_xxz = primBuffer.data(60 * idx + 52);

            auto t_zz_xyy = primBuffer.data(60 * idx + 53);

            auto t_zz_xyz = primBuffer.data(60 * idx + 54);

            auto t_zz_xzz = primBuffer.data(60 * idx + 55);

            auto t_zz_yyy = primBuffer.data(60 * idx + 56);

            auto t_zz_yyz = primBuffer.data(60 * idx + 57);

            auto t_zz_yzz = primBuffer.data(60 * idx + 58);

            auto t_zz_zzz = primBuffer.data(60 * idx + 59);

            // Batch of Integrals (50,60)

            #pragma omp simd aligned(fx, pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, \
                                     pa2pb_z_zz, pa2pb_zz_x, pa2pb_zz_xxx, pa2pb_zz_xxy, pa2pb_zz_xxz, pa2pb_zz_xyy, \
                                     pa2pb_zz_xyz, pa2pb_zz_xzz, pa2pb_zz_y, pa2pb_zz_yyy, pa2pb_zz_yyz, pa2pb_zz_yzz, \
                                     pa2pb_zz_z, pa2pb_zz_zzz, pa_z, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_xzz, pb_y, \
                                     pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, s_0_0, t_zz_xxx, t_zz_xxy, t_zz_xxz, t_zz_xyy, \
                                     t_zz_xyz, t_zz_xzz, t_zz_yyy, t_zz_yyz, t_zz_yzz, t_zz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_zz_xxx[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 1.5 * pa2pb_zz_x[j] * fl1_fx + 0.5 * pb_xxx[j] * fl1_fx + pa2pb_zz_xxx[j]);

                t_zz_xxy[j] = fl_s_0_0 * (0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_zz_y[j] * fl1_fx + 0.5 * pb_xxy[j] * fl1_fx + pa2pb_zz_xxy[j]);

                t_zz_xxz[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl2_fx + 0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_zz_z[j] * fl1_fx + pa2pb_z_xx[j] * fl1_fx + 0.5 * pb_xxz[j] * fl1_fx + pa2pb_zz_xxz[j]);

                t_zz_xyy[j] = fl_s_0_0 * (0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_zz_x[j] * fl1_fx + 0.5 * pb_xyy[j] * fl1_fx + pa2pb_zz_xyy[j]);

                t_zz_xyz[j] = fl_s_0_0 * (pa2pb_z_xy[j] * fl1_fx + 0.5 * pb_xyz[j] * fl1_fx + pa2pb_zz_xyz[j]);

                t_zz_xzz[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 0.5 * pa2pb_zz_x[j] * fl1_fx + 2.0 * pa2pb_z_xz[j] * fl1_fx + 0.5 * pb_xzz[j] * fl1_fx + pa2pb_zz_xzz[j]);

                t_zz_yyy[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 1.5 * pa2pb_zz_y[j] * fl1_fx + 0.5 * pb_yyy[j] * fl1_fx + pa2pb_zz_yyy[j]);

                t_zz_yyz[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl2_fx + 0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_zz_z[j] * fl1_fx + pa2pb_z_yy[j] * fl1_fx + 0.5 * pb_yyz[j] * fl1_fx + pa2pb_zz_yyz[j]);

                t_zz_yzz[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 0.5 * pa2pb_zz_y[j] * fl1_fx + 2.0 * pa2pb_z_yz[j] * fl1_fx + 0.5 * pb_yzz[j] * fl1_fx + pa2pb_zz_yzz[j]);

                t_zz_zzz[j] = fl_s_0_0 * (1.5 * pa_z[j] * fl2_fx + 2.25 * pb_z[j] * fl2_fx + 1.5 * pa2pb_zz_z[j] * fl1_fx + 3.0 * pa2pb_z_zz[j] * fl1_fx + 0.5 * pb_zzz[j] * fl1_fx + pa2pb_zz_zzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForFD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForFD_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForFD_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForFD_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForFD_30_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForFD_40_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForFD_50_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForFD_0_10(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(171 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(171 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(171 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(171 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(171 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(171 * idx + 8);

            auto pa2pb_y_xx = pa2pbDistances.data(171 * idx + 12);

            auto pa2pb_y_xy = pa2pbDistances.data(171 * idx + 13);

            auto pa2pb_y_xz = pa2pbDistances.data(171 * idx + 14);

            auto pa2pb_y_yy = pa2pbDistances.data(171 * idx + 15);

            auto pa2pb_xx_x = pa2pbDistances.data(171 * idx + 27);

            auto pa2pb_xx_y = pa2pbDistances.data(171 * idx + 28);

            auto pa2pb_xx_z = pa2pbDistances.data(171 * idx + 29);

            auto pa2pb_xy_x = pa2pbDistances.data(171 * idx + 36);

            auto pa2pb_xy_y = pa2pbDistances.data(171 * idx + 37);

            auto pa2pb_xy_z = pa2pbDistances.data(171 * idx + 38);

            auto pa2pb_xxx_xx = pa2pbDistances.data(171 * idx + 84);

            auto pa2pb_xxx_xy = pa2pbDistances.data(171 * idx + 85);

            auto pa2pb_xxx_xz = pa2pbDistances.data(171 * idx + 86);

            auto pa2pb_xxx_yy = pa2pbDistances.data(171 * idx + 87);

            auto pa2pb_xxx_yz = pa2pbDistances.data(171 * idx + 88);

            auto pa2pb_xxx_zz = pa2pbDistances.data(171 * idx + 89);

            auto pa2pb_xxy_xx = pa2pbDistances.data(171 * idx + 93);

            auto pa2pb_xxy_xy = pa2pbDistances.data(171 * idx + 94);

            auto pa2pb_xxy_xz = pa2pbDistances.data(171 * idx + 95);

            auto pa2pb_xxy_yy = pa2pbDistances.data(171 * idx + 96);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxx_xx = primBuffer.data(60 * idx);

            auto t_xxx_xy = primBuffer.data(60 * idx + 1);

            auto t_xxx_xz = primBuffer.data(60 * idx + 2);

            auto t_xxx_yy = primBuffer.data(60 * idx + 3);

            auto t_xxx_yz = primBuffer.data(60 * idx + 4);

            auto t_xxx_zz = primBuffer.data(60 * idx + 5);

            auto t_xxy_xx = primBuffer.data(60 * idx + 6);

            auto t_xxy_xy = primBuffer.data(60 * idx + 7);

            auto t_xxy_xz = primBuffer.data(60 * idx + 8);

            auto t_xxy_yy = primBuffer.data(60 * idx + 9);

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_x_zz, pa2pb_xx_x, pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxx_xx, pa2pb_xxx_xy, \
                                     pa2pb_xxx_xz, pa2pb_xxx_yy, pa2pb_xxx_yz, pa2pb_xxx_zz, pa2pb_xxy_xx, pa2pb_xxy_xy, \
                                     pa2pb_xxy_xz, pa2pb_xxy_yy, pa2pb_xy_x, pa2pb_xy_y, pa2pb_xy_z, pa2pb_y_xx, \
                                     pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, pa_x, pa_xxx, pa_xxy, pa_y, pb_x, pb_y, pb_z, s_0_0, \
                                     t_xxx_xx, t_xxx_xy, t_xxx_xz, t_xxx_yy, t_xxx_yz, t_xxx_zz, t_xxy_xx, t_xxy_xy, \
                                     t_xxy_xz, t_xxy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxx_xx[j] = fl_s_0_0 * (2.25 * pa_x[j] * fl2_fx + 1.5 * pb_x[j] * fl2_fx + 0.5 * pa_xxx[j] * fl1_fx + 3.0 * pa2pb_xx_x[j] * fl1_fx + 1.5 * pa2pb_x_xx[j] * fl1_fx + pa2pb_xxx_xx[j]);

                t_xxx_xy[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 1.5 * pa2pb_xx_y[j] * fl1_fx + 1.5 * pa2pb_x_xy[j] * fl1_fx + pa2pb_xxx_xy[j]);

                t_xxx_xz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 1.5 * pa2pb_xx_z[j] * fl1_fx + 1.5 * pa2pb_x_xz[j] * fl1_fx + pa2pb_xxx_xz[j]);

                t_xxx_yy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 0.5 * pa_xxx[j] * fl1_fx + 1.5 * pa2pb_x_yy[j] * fl1_fx + pa2pb_xxx_yy[j]);

                t_xxx_yz[j] = fl_s_0_0 * (1.5 * pa2pb_x_yz[j] * fl1_fx + pa2pb_xxx_yz[j]);

                t_xxx_zz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 0.5 * pa_xxx[j] * fl1_fx + 1.5 * pa2pb_x_zz[j] * fl1_fx + pa2pb_xxx_zz[j]);

                t_xxy_xx[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 0.5 * pa_xxy[j] * fl1_fx + 2.0 * pa2pb_xy_x[j] * fl1_fx + 0.5 * pa2pb_y_xx[j] * fl1_fx + pa2pb_xxy_xx[j]);

                t_xxy_xy[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl2_fx + 0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_xx_x[j] * fl1_fx + pa2pb_xy_y[j] * fl1_fx + 0.5 * pa2pb_y_xy[j] * fl1_fx + pa2pb_xxy_xy[j]);

                t_xxy_xz[j] = fl_s_0_0 * (pa2pb_xy_z[j] * fl1_fx + 0.5 * pa2pb_y_xz[j] * fl1_fx + pa2pb_xxy_xz[j]);

                t_xxy_yy[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pb_y[j] * fl2_fx + 0.5 * pa_xxy[j] * fl1_fx + pa2pb_xx_y[j] * fl1_fx + 0.5 * pa2pb_y_yy[j] * fl1_fx + pa2pb_xxy_yy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForFD_10_20(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (10,20)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(171 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(171 * idx + 4);

            auto pa2pb_y_yz = pa2pbDistances.data(171 * idx + 16);

            auto pa2pb_y_zz = pa2pbDistances.data(171 * idx + 17);

            auto pa2pb_z_xx = pa2pbDistances.data(171 * idx + 21);

            auto pa2pb_z_xy = pa2pbDistances.data(171 * idx + 22);

            auto pa2pb_z_xz = pa2pbDistances.data(171 * idx + 23);

            auto pa2pb_z_yy = pa2pbDistances.data(171 * idx + 24);

            auto pa2pb_z_yz = pa2pbDistances.data(171 * idx + 25);

            auto pa2pb_z_zz = pa2pbDistances.data(171 * idx + 26);

            auto pa2pb_xx_x = pa2pbDistances.data(171 * idx + 27);

            auto pa2pb_xx_y = pa2pbDistances.data(171 * idx + 28);

            auto pa2pb_xx_z = pa2pbDistances.data(171 * idx + 29);

            auto pa2pb_xy_x = pa2pbDistances.data(171 * idx + 36);

            auto pa2pb_xz_x = pa2pbDistances.data(171 * idx + 45);

            auto pa2pb_xz_y = pa2pbDistances.data(171 * idx + 46);

            auto pa2pb_xz_z = pa2pbDistances.data(171 * idx + 47);

            auto pa2pb_yy_x = pa2pbDistances.data(171 * idx + 54);

            auto pa2pb_yy_y = pa2pbDistances.data(171 * idx + 55);

            auto pa2pb_xxy_yz = pa2pbDistances.data(171 * idx + 97);

            auto pa2pb_xxy_zz = pa2pbDistances.data(171 * idx + 98);

            auto pa2pb_xxz_xx = pa2pbDistances.data(171 * idx + 102);

            auto pa2pb_xxz_xy = pa2pbDistances.data(171 * idx + 103);

            auto pa2pb_xxz_xz = pa2pbDistances.data(171 * idx + 104);

            auto pa2pb_xxz_yy = pa2pbDistances.data(171 * idx + 105);

            auto pa2pb_xxz_yz = pa2pbDistances.data(171 * idx + 106);

            auto pa2pb_xxz_zz = pa2pbDistances.data(171 * idx + 107);

            auto pa2pb_xyy_xx = pa2pbDistances.data(171 * idx + 111);

            auto pa2pb_xyy_xy = pa2pbDistances.data(171 * idx + 112);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxy_yz = primBuffer.data(60 * idx + 10);

            auto t_xxy_zz = primBuffer.data(60 * idx + 11);

            auto t_xxz_xx = primBuffer.data(60 * idx + 12);

            auto t_xxz_xy = primBuffer.data(60 * idx + 13);

            auto t_xxz_xz = primBuffer.data(60 * idx + 14);

            auto t_xxz_yy = primBuffer.data(60 * idx + 15);

            auto t_xxz_yz = primBuffer.data(60 * idx + 16);

            auto t_xxz_zz = primBuffer.data(60 * idx + 17);

            auto t_xyy_xx = primBuffer.data(60 * idx + 18);

            auto t_xyy_xy = primBuffer.data(60 * idx + 19);

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fx, pa2pb_x_xx, pa2pb_x_xy, pa2pb_xx_x, pa2pb_xx_y, pa2pb_xx_z, \
                                     pa2pb_xxy_yz, pa2pb_xxy_zz, pa2pb_xxz_xx, pa2pb_xxz_xy, pa2pb_xxz_xz, pa2pb_xxz_yy, \
                                     pa2pb_xxz_yz, pa2pb_xxz_zz, pa2pb_xy_x, pa2pb_xyy_xx, pa2pb_xyy_xy, pa2pb_xz_x, \
                                     pa2pb_xz_y, pa2pb_xz_z, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yy_x, pa2pb_yy_y, pa2pb_z_xx, \
                                     pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa_x, pa_xxy, pa_xxz, \
                                     pa_xyy, pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xxy_yz, t_xxy_zz, t_xxz_xx, t_xxz_xy, \
                                     t_xxz_xz, t_xxz_yy, t_xxz_yz, t_xxz_zz, t_xyy_xx, t_xyy_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxy_yz[j] = fl_s_0_0 * (0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_xx_z[j] * fl1_fx + 0.5 * pa2pb_y_yz[j] * fl1_fx + pa2pb_xxy_yz[j]);

                t_xxy_zz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa_xxy[j] * fl1_fx + 0.5 * pa2pb_y_zz[j] * fl1_fx + pa2pb_xxy_zz[j]);

                t_xxz_xx[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 0.5 * pa_xxz[j] * fl1_fx + 2.0 * pa2pb_xz_x[j] * fl1_fx + 0.5 * pa2pb_z_xx[j] * fl1_fx + pa2pb_xxz_xx[j]);

                t_xxz_xy[j] = fl_s_0_0 * (pa2pb_xz_y[j] * fl1_fx + 0.5 * pa2pb_z_xy[j] * fl1_fx + pa2pb_xxz_xy[j]);

                t_xxz_xz[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl2_fx + 0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_xx_x[j] * fl1_fx + pa2pb_xz_z[j] * fl1_fx + 0.5 * pa2pb_z_xz[j] * fl1_fx + pa2pb_xxz_xz[j]);

                t_xxz_yy[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa_xxz[j] * fl1_fx + 0.5 * pa2pb_z_yy[j] * fl1_fx + pa2pb_xxz_yy[j]);

                t_xxz_yz[j] = fl_s_0_0 * (0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_xx_y[j] * fl1_fx + 0.5 * pa2pb_z_yz[j] * fl1_fx + pa2pb_xxz_yz[j]);

                t_xxz_zz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pb_z[j] * fl2_fx + 0.5 * pa_xxz[j] * fl1_fx + pa2pb_xx_z[j] * fl1_fx + 0.5 * pa2pb_z_zz[j] * fl1_fx + pa2pb_xxz_zz[j]);

                t_xyy_xx[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pb_x[j] * fl2_fx + 0.5 * pa_xyy[j] * fl1_fx + pa2pb_yy_x[j] * fl1_fx + 0.5 * pa2pb_x_xx[j] * fl1_fx + pa2pb_xyy_xx[j]);

                t_xyy_xy[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl2_fx + 0.25 * pb_y[j] * fl2_fx + pa2pb_xy_x[j] * fl1_fx + 0.5 * pa2pb_yy_y[j] * fl1_fx + 0.5 * pa2pb_x_xy[j] * fl1_fx + pa2pb_xyy_xy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForFD_20_30(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (20,30)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xz = pa2pbDistances.data(171 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(171 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(171 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(171 * idx + 8);

            auto pa2pb_xy_x = pa2pbDistances.data(171 * idx + 36);

            auto pa2pb_xy_y = pa2pbDistances.data(171 * idx + 37);

            auto pa2pb_xy_z = pa2pbDistances.data(171 * idx + 38);

            auto pa2pb_xz_x = pa2pbDistances.data(171 * idx + 45);

            auto pa2pb_xz_y = pa2pbDistances.data(171 * idx + 46);

            auto pa2pb_xz_z = pa2pbDistances.data(171 * idx + 47);

            auto pa2pb_yy_z = pa2pbDistances.data(171 * idx + 56);

            auto pa2pb_yz_x = pa2pbDistances.data(171 * idx + 63);

            auto pa2pb_yz_y = pa2pbDistances.data(171 * idx + 64);

            auto pa2pb_yz_z = pa2pbDistances.data(171 * idx + 65);

            auto pa2pb_xyy_xz = pa2pbDistances.data(171 * idx + 113);

            auto pa2pb_xyy_yy = pa2pbDistances.data(171 * idx + 114);

            auto pa2pb_xyy_yz = pa2pbDistances.data(171 * idx + 115);

            auto pa2pb_xyy_zz = pa2pbDistances.data(171 * idx + 116);

            auto pa2pb_xyz_xx = pa2pbDistances.data(171 * idx + 120);

            auto pa2pb_xyz_xy = pa2pbDistances.data(171 * idx + 121);

            auto pa2pb_xyz_xz = pa2pbDistances.data(171 * idx + 122);

            auto pa2pb_xyz_yy = pa2pbDistances.data(171 * idx + 123);

            auto pa2pb_xyz_yz = pa2pbDistances.data(171 * idx + 124);

            auto pa2pb_xyz_zz = pa2pbDistances.data(171 * idx + 125);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xyy_xz = primBuffer.data(60 * idx + 20);

            auto t_xyy_yy = primBuffer.data(60 * idx + 21);

            auto t_xyy_yz = primBuffer.data(60 * idx + 22);

            auto t_xyy_zz = primBuffer.data(60 * idx + 23);

            auto t_xyz_xx = primBuffer.data(60 * idx + 24);

            auto t_xyz_xy = primBuffer.data(60 * idx + 25);

            auto t_xyz_xz = primBuffer.data(60 * idx + 26);

            auto t_xyz_yy = primBuffer.data(60 * idx + 27);

            auto t_xyz_yz = primBuffer.data(60 * idx + 28);

            auto t_xyz_zz = primBuffer.data(60 * idx + 29);

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fx, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xy_x, \
                                     pa2pb_xy_y, pa2pb_xy_z, pa2pb_xyy_xz, pa2pb_xyy_yy, pa2pb_xyy_yz, pa2pb_xyy_zz, \
                                     pa2pb_xyz_xx, pa2pb_xyz_xy, pa2pb_xyz_xz, pa2pb_xyz_yy, pa2pb_xyz_yz, pa2pb_xyz_zz, \
                                     pa2pb_xz_x, pa2pb_xz_y, pa2pb_xz_z, pa2pb_yy_z, pa2pb_yz_x, pa2pb_yz_y, pa2pb_yz_z, \
                                     pa_x, pa_xyy, pa_xyz, pa_y, pa_z, pb_z, s_0_0, t_xyy_xz, t_xyy_yy, t_xyy_yz, \
                                     t_xyy_zz, t_xyz_xx, t_xyz_xy, t_xyz_xz, t_xyz_yy, t_xyz_yz, t_xyz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xyy_xz[j] = fl_s_0_0 * (0.25 * pb_z[j] * fl2_fx + 0.5 * pa2pb_yy_z[j] * fl1_fx + 0.5 * pa2pb_x_xz[j] * fl1_fx + pa2pb_xyy_xz[j]);

                t_xyy_yy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 0.5 * pa_xyy[j] * fl1_fx + 2.0 * pa2pb_xy_y[j] * fl1_fx + 0.5 * pa2pb_x_yy[j] * fl1_fx + pa2pb_xyy_yy[j]);

                t_xyy_yz[j] = fl_s_0_0 * (pa2pb_xy_z[j] * fl1_fx + 0.5 * pa2pb_x_yz[j] * fl1_fx + pa2pb_xyy_yz[j]);

                t_xyy_zz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa_xyy[j] * fl1_fx + 0.5 * pa2pb_x_zz[j] * fl1_fx + pa2pb_xyy_zz[j]);

                t_xyz_xx[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl1_fx + pa2pb_yz_x[j] * fl1_fx + pa2pb_xyz_xx[j]);

                t_xyz_xy[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa2pb_xz_x[j] * fl1_fx + 0.5 * pa2pb_yz_y[j] * fl1_fx + pa2pb_xyz_xy[j]);

                t_xyz_xz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa2pb_xy_x[j] * fl1_fx + 0.5 * pa2pb_yz_z[j] * fl1_fx + pa2pb_xyz_xz[j]);

                t_xyz_yy[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl1_fx + pa2pb_xz_y[j] * fl1_fx + pa2pb_xyz_yy[j]);

                t_xyz_yz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa2pb_xy_y[j] * fl1_fx + 0.5 * pa2pb_xz_z[j] * fl1_fx + pa2pb_xyz_yz[j]);

                t_xyz_zz[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl1_fx + pa2pb_xy_z[j] * fl1_fx + pa2pb_xyz_zz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForFD_30_40(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (30,40)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(171 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(171 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(171 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(171 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(171 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(171 * idx + 8);

            auto pa2pb_y_xx = pa2pbDistances.data(171 * idx + 12);

            auto pa2pb_y_xy = pa2pbDistances.data(171 * idx + 13);

            auto pa2pb_y_xz = pa2pbDistances.data(171 * idx + 14);

            auto pa2pb_y_yy = pa2pbDistances.data(171 * idx + 15);

            auto pa2pb_xz_x = pa2pbDistances.data(171 * idx + 45);

            auto pa2pb_xz_y = pa2pbDistances.data(171 * idx + 46);

            auto pa2pb_xz_z = pa2pbDistances.data(171 * idx + 47);

            auto pa2pb_yy_x = pa2pbDistances.data(171 * idx + 54);

            auto pa2pb_yy_y = pa2pbDistances.data(171 * idx + 55);

            auto pa2pb_zz_x = pa2pbDistances.data(171 * idx + 72);

            auto pa2pb_zz_y = pa2pbDistances.data(171 * idx + 73);

            auto pa2pb_zz_z = pa2pbDistances.data(171 * idx + 74);

            auto pa2pb_xzz_xx = pa2pbDistances.data(171 * idx + 129);

            auto pa2pb_xzz_xy = pa2pbDistances.data(171 * idx + 130);

            auto pa2pb_xzz_xz = pa2pbDistances.data(171 * idx + 131);

            auto pa2pb_xzz_yy = pa2pbDistances.data(171 * idx + 132);

            auto pa2pb_xzz_yz = pa2pbDistances.data(171 * idx + 133);

            auto pa2pb_xzz_zz = pa2pbDistances.data(171 * idx + 134);

            auto pa2pb_yyy_xx = pa2pbDistances.data(171 * idx + 138);

            auto pa2pb_yyy_xy = pa2pbDistances.data(171 * idx + 139);

            auto pa2pb_yyy_xz = pa2pbDistances.data(171 * idx + 140);

            auto pa2pb_yyy_yy = pa2pbDistances.data(171 * idx + 141);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xzz_xx = primBuffer.data(60 * idx + 30);

            auto t_xzz_xy = primBuffer.data(60 * idx + 31);

            auto t_xzz_xz = primBuffer.data(60 * idx + 32);

            auto t_xzz_yy = primBuffer.data(60 * idx + 33);

            auto t_xzz_yz = primBuffer.data(60 * idx + 34);

            auto t_xzz_zz = primBuffer.data(60 * idx + 35);

            auto t_yyy_xx = primBuffer.data(60 * idx + 36);

            auto t_yyy_xy = primBuffer.data(60 * idx + 37);

            auto t_yyy_xz = primBuffer.data(60 * idx + 38);

            auto t_yyy_yy = primBuffer.data(60 * idx + 39);

            // Batch of Integrals (30,40)

            #pragma omp simd aligned(fx, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_x_zz, pa2pb_xz_x, pa2pb_xz_y, pa2pb_xz_z, pa2pb_xzz_xx, pa2pb_xzz_xy, \
                                     pa2pb_xzz_xz, pa2pb_xzz_yy, pa2pb_xzz_yz, pa2pb_xzz_zz, pa2pb_y_xx, pa2pb_y_xy, \
                                     pa2pb_y_xz, pa2pb_y_yy, pa2pb_yy_x, pa2pb_yy_y, pa2pb_yyy_xx, pa2pb_yyy_xy, \
                                     pa2pb_yyy_xz, pa2pb_yyy_yy, pa2pb_zz_x, pa2pb_zz_y, pa2pb_zz_z, pa_x, pa_xzz, pa_y, \
                                     pa_yyy, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xzz_xx, t_xzz_xy, t_xzz_xz, t_xzz_yy, \
                                     t_xzz_yz, t_xzz_zz, t_yyy_xx, t_yyy_xy, t_yyy_xz, t_yyy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xzz_xx[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pb_x[j] * fl2_fx + 0.5 * pa_xzz[j] * fl1_fx + pa2pb_zz_x[j] * fl1_fx + 0.5 * pa2pb_x_xx[j] * fl1_fx + pa2pb_xzz_xx[j]);

                t_xzz_xy[j] = fl_s_0_0 * (0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_zz_y[j] * fl1_fx + 0.5 * pa2pb_x_xy[j] * fl1_fx + pa2pb_xzz_xy[j]);

                t_xzz_xz[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl2_fx + 0.25 * pb_z[j] * fl2_fx + pa2pb_xz_x[j] * fl1_fx + 0.5 * pa2pb_zz_z[j] * fl1_fx + 0.5 * pa2pb_x_xz[j] * fl1_fx + pa2pb_xzz_xz[j]);

                t_xzz_yy[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa_xzz[j] * fl1_fx + 0.5 * pa2pb_x_yy[j] * fl1_fx + pa2pb_xzz_yy[j]);

                t_xzz_yz[j] = fl_s_0_0 * (pa2pb_xz_y[j] * fl1_fx + 0.5 * pa2pb_x_yz[j] * fl1_fx + pa2pb_xzz_yz[j]);

                t_xzz_zz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx + 0.5 * pa_xzz[j] * fl1_fx + 2.0 * pa2pb_xz_z[j] * fl1_fx + 0.5 * pa2pb_x_zz[j] * fl1_fx + pa2pb_xzz_zz[j]);

                t_yyy_xx[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 0.5 * pa_yyy[j] * fl1_fx + 1.5 * pa2pb_y_xx[j] * fl1_fx + pa2pb_yyy_xx[j]);

                t_yyy_xy[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 1.5 * pa2pb_yy_x[j] * fl1_fx + 1.5 * pa2pb_y_xy[j] * fl1_fx + pa2pb_yyy_xy[j]);

                t_yyy_xz[j] = fl_s_0_0 * (1.5 * pa2pb_y_xz[j] * fl1_fx + pa2pb_yyy_xz[j]);

                t_yyy_yy[j] = fl_s_0_0 * (2.25 * pa_y[j] * fl2_fx + 1.5 * pb_y[j] * fl2_fx + 0.5 * pa_yyy[j] * fl1_fx + 3.0 * pa2pb_yy_y[j] * fl1_fx + 1.5 * pa2pb_y_yy[j] * fl1_fx + pa2pb_yyy_yy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForFD_40_50(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (40,50)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(171 * idx + 12);

            auto pa2pb_y_xy = pa2pbDistances.data(171 * idx + 13);

            auto pa2pb_y_yz = pa2pbDistances.data(171 * idx + 16);

            auto pa2pb_y_zz = pa2pbDistances.data(171 * idx + 17);

            auto pa2pb_z_xx = pa2pbDistances.data(171 * idx + 21);

            auto pa2pb_z_xy = pa2pbDistances.data(171 * idx + 22);

            auto pa2pb_z_xz = pa2pbDistances.data(171 * idx + 23);

            auto pa2pb_z_yy = pa2pbDistances.data(171 * idx + 24);

            auto pa2pb_z_yz = pa2pbDistances.data(171 * idx + 25);

            auto pa2pb_z_zz = pa2pbDistances.data(171 * idx + 26);

            auto pa2pb_yy_x = pa2pbDistances.data(171 * idx + 54);

            auto pa2pb_yy_y = pa2pbDistances.data(171 * idx + 55);

            auto pa2pb_yy_z = pa2pbDistances.data(171 * idx + 56);

            auto pa2pb_yz_x = pa2pbDistances.data(171 * idx + 63);

            auto pa2pb_yz_y = pa2pbDistances.data(171 * idx + 64);

            auto pa2pb_yz_z = pa2pbDistances.data(171 * idx + 65);

            auto pa2pb_zz_x = pa2pbDistances.data(171 * idx + 72);

            auto pa2pb_yyy_yz = pa2pbDistances.data(171 * idx + 142);

            auto pa2pb_yyy_zz = pa2pbDistances.data(171 * idx + 143);

            auto pa2pb_yyz_xx = pa2pbDistances.data(171 * idx + 147);

            auto pa2pb_yyz_xy = pa2pbDistances.data(171 * idx + 148);

            auto pa2pb_yyz_xz = pa2pbDistances.data(171 * idx + 149);

            auto pa2pb_yyz_yy = pa2pbDistances.data(171 * idx + 150);

            auto pa2pb_yyz_yz = pa2pbDistances.data(171 * idx + 151);

            auto pa2pb_yyz_zz = pa2pbDistances.data(171 * idx + 152);

            auto pa2pb_yzz_xx = pa2pbDistances.data(171 * idx + 156);

            auto pa2pb_yzz_xy = pa2pbDistances.data(171 * idx + 157);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yyy_yz = primBuffer.data(60 * idx + 40);

            auto t_yyy_zz = primBuffer.data(60 * idx + 41);

            auto t_yyz_xx = primBuffer.data(60 * idx + 42);

            auto t_yyz_xy = primBuffer.data(60 * idx + 43);

            auto t_yyz_xz = primBuffer.data(60 * idx + 44);

            auto t_yyz_yy = primBuffer.data(60 * idx + 45);

            auto t_yyz_yz = primBuffer.data(60 * idx + 46);

            auto t_yyz_zz = primBuffer.data(60 * idx + 47);

            auto t_yzz_xx = primBuffer.data(60 * idx + 48);

            auto t_yzz_xy = primBuffer.data(60 * idx + 49);

            // Batch of Integrals (40,50)

            #pragma omp simd aligned(fx, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yy_x, \
                                     pa2pb_yy_y, pa2pb_yy_z, pa2pb_yyy_yz, pa2pb_yyy_zz, pa2pb_yyz_xx, pa2pb_yyz_xy, \
                                     pa2pb_yyz_xz, pa2pb_yyz_yy, pa2pb_yyz_yz, pa2pb_yyz_zz, pa2pb_yz_x, pa2pb_yz_y, \
                                     pa2pb_yz_z, pa2pb_yzz_xx, pa2pb_yzz_xy, pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, \
                                     pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa2pb_zz_x, pa_y, pa_yyy, pa_yyz, pa_yzz, pa_z, pb_x, \
                                     pb_y, pb_z, s_0_0, t_yyy_yz, t_yyy_zz, t_yyz_xx, t_yyz_xy, t_yyz_xz, t_yyz_yy, \
                                     t_yyz_yz, t_yyz_zz, t_yzz_xx, t_yzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yyy_yz[j] = fl_s_0_0 * (0.75 * pb_z[j] * fl2_fx + 1.5 * pa2pb_yy_z[j] * fl1_fx + 1.5 * pa2pb_y_yz[j] * fl1_fx + pa2pb_yyy_yz[j]);

                t_yyy_zz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 0.5 * pa_yyy[j] * fl1_fx + 1.5 * pa2pb_y_zz[j] * fl1_fx + pa2pb_yyy_zz[j]);

                t_yyz_xx[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa_yyz[j] * fl1_fx + 0.5 * pa2pb_z_xx[j] * fl1_fx + pa2pb_yyz_xx[j]);

                t_yyz_xy[j] = fl_s_0_0 * (pa2pb_yz_x[j] * fl1_fx + 0.5 * pa2pb_z_xy[j] * fl1_fx + pa2pb_yyz_xy[j]);

                t_yyz_xz[j] = fl_s_0_0 * (0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_yy_x[j] * fl1_fx + 0.5 * pa2pb_z_xz[j] * fl1_fx + pa2pb_yyz_xz[j]);

                t_yyz_yy[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 0.5 * pa_yyz[j] * fl1_fx + 2.0 * pa2pb_yz_y[j] * fl1_fx + 0.5 * pa2pb_z_yy[j] * fl1_fx + pa2pb_yyz_yy[j]);

                t_yyz_yz[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl2_fx + 0.25 * pb_y[j] * fl2_fx + 0.5 * pa2pb_yy_y[j] * fl1_fx + pa2pb_yz_z[j] * fl1_fx + 0.5 * pa2pb_z_yz[j] * fl1_fx + pa2pb_yyz_yz[j]);

                t_yyz_zz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pb_z[j] * fl2_fx + 0.5 * pa_yyz[j] * fl1_fx + pa2pb_yy_z[j] * fl1_fx + 0.5 * pa2pb_z_zz[j] * fl1_fx + pa2pb_yyz_zz[j]);

                t_yzz_xx[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa_yzz[j] * fl1_fx + 0.5 * pa2pb_y_xx[j] * fl1_fx + pa2pb_yzz_xx[j]);

                t_yzz_xy[j] = fl_s_0_0 * (0.25 * pb_x[j] * fl2_fx + 0.5 * pa2pb_zz_x[j] * fl1_fx + 0.5 * pa2pb_y_xy[j] * fl1_fx + pa2pb_yzz_xy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForFD_50_60(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (50,60)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xz = pa2pbDistances.data(171 * idx + 14);

            auto pa2pb_y_yy = pa2pbDistances.data(171 * idx + 15);

            auto pa2pb_y_yz = pa2pbDistances.data(171 * idx + 16);

            auto pa2pb_y_zz = pa2pbDistances.data(171 * idx + 17);

            auto pa2pb_z_xx = pa2pbDistances.data(171 * idx + 21);

            auto pa2pb_z_xy = pa2pbDistances.data(171 * idx + 22);

            auto pa2pb_z_xz = pa2pbDistances.data(171 * idx + 23);

            auto pa2pb_z_yy = pa2pbDistances.data(171 * idx + 24);

            auto pa2pb_z_yz = pa2pbDistances.data(171 * idx + 25);

            auto pa2pb_z_zz = pa2pbDistances.data(171 * idx + 26);

            auto pa2pb_yz_x = pa2pbDistances.data(171 * idx + 63);

            auto pa2pb_yz_y = pa2pbDistances.data(171 * idx + 64);

            auto pa2pb_yz_z = pa2pbDistances.data(171 * idx + 65);

            auto pa2pb_zz_x = pa2pbDistances.data(171 * idx + 72);

            auto pa2pb_zz_y = pa2pbDistances.data(171 * idx + 73);

            auto pa2pb_zz_z = pa2pbDistances.data(171 * idx + 74);

            auto pa2pb_yzz_xz = pa2pbDistances.data(171 * idx + 158);

            auto pa2pb_yzz_yy = pa2pbDistances.data(171 * idx + 159);

            auto pa2pb_yzz_yz = pa2pbDistances.data(171 * idx + 160);

            auto pa2pb_yzz_zz = pa2pbDistances.data(171 * idx + 161);

            auto pa2pb_zzz_xx = pa2pbDistances.data(171 * idx + 165);

            auto pa2pb_zzz_xy = pa2pbDistances.data(171 * idx + 166);

            auto pa2pb_zzz_xz = pa2pbDistances.data(171 * idx + 167);

            auto pa2pb_zzz_yy = pa2pbDistances.data(171 * idx + 168);

            auto pa2pb_zzz_yz = pa2pbDistances.data(171 * idx + 169);

            auto pa2pb_zzz_zz = pa2pbDistances.data(171 * idx + 170);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yzz_xz = primBuffer.data(60 * idx + 50);

            auto t_yzz_yy = primBuffer.data(60 * idx + 51);

            auto t_yzz_yz = primBuffer.data(60 * idx + 52);

            auto t_yzz_zz = primBuffer.data(60 * idx + 53);

            auto t_zzz_xx = primBuffer.data(60 * idx + 54);

            auto t_zzz_xy = primBuffer.data(60 * idx + 55);

            auto t_zzz_xz = primBuffer.data(60 * idx + 56);

            auto t_zzz_yy = primBuffer.data(60 * idx + 57);

            auto t_zzz_yz = primBuffer.data(60 * idx + 58);

            auto t_zzz_zz = primBuffer.data(60 * idx + 59);

            // Batch of Integrals (50,60)

            #pragma omp simd aligned(fx, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yz_x, \
                                     pa2pb_yz_y, pa2pb_yz_z, pa2pb_yzz_xz, pa2pb_yzz_yy, pa2pb_yzz_yz, pa2pb_yzz_zz, \
                                     pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa2pb_zz_x, \
                                     pa2pb_zz_y, pa2pb_zz_z, pa2pb_zzz_xx, pa2pb_zzz_xy, pa2pb_zzz_xz, pa2pb_zzz_yy, \
                                     pa2pb_zzz_yz, pa2pb_zzz_zz, pa_y, pa_yzz, pa_z, pa_zzz, pb_x, pb_y, pb_z, s_0_0, t_yzz_xz, \
                                     t_yzz_yy, t_yzz_yz, t_yzz_zz, t_zzz_xx, t_zzz_xy, t_zzz_xz, t_zzz_yy, t_zzz_yz, \
                                     t_zzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yzz_xz[j] = fl_s_0_0 * (pa2pb_yz_x[j] * fl1_fx + 0.5 * pa2pb_y_xz[j] * fl1_fx + pa2pb_yzz_xz[j]);

                t_yzz_yy[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pb_y[j] * fl2_fx + 0.5 * pa_yzz[j] * fl1_fx + pa2pb_zz_y[j] * fl1_fx + 0.5 * pa2pb_y_yy[j] * fl1_fx + pa2pb_yzz_yy[j]);

                t_yzz_yz[j] = fl_s_0_0 * (0.5 * pa_z[j] * fl2_fx + 0.25 * pb_z[j] * fl2_fx + pa2pb_yz_y[j] * fl1_fx + 0.5 * pa2pb_zz_z[j] * fl1_fx + 0.5 * pa2pb_y_yz[j] * fl1_fx + pa2pb_yzz_yz[j]);

                t_yzz_zz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx + 0.5 * pa_yzz[j] * fl1_fx + 2.0 * pa2pb_yz_z[j] * fl1_fx + 0.5 * pa2pb_y_zz[j] * fl1_fx + pa2pb_yzz_zz[j]);

                t_zzz_xx[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 0.5 * pa_zzz[j] * fl1_fx + 1.5 * pa2pb_z_xx[j] * fl1_fx + pa2pb_zzz_xx[j]);

                t_zzz_xy[j] = fl_s_0_0 * (1.5 * pa2pb_z_xy[j] * fl1_fx + pa2pb_zzz_xy[j]);

                t_zzz_xz[j] = fl_s_0_0 * (0.75 * pb_x[j] * fl2_fx + 1.5 * pa2pb_zz_x[j] * fl1_fx + 1.5 * pa2pb_z_xz[j] * fl1_fx + pa2pb_zzz_xz[j]);

                t_zzz_yy[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx + 0.5 * pa_zzz[j] * fl1_fx + 1.5 * pa2pb_z_yy[j] * fl1_fx + pa2pb_zzz_yy[j]);

                t_zzz_yz[j] = fl_s_0_0 * (0.75 * pb_y[j] * fl2_fx + 1.5 * pa2pb_zz_y[j] * fl1_fx + 1.5 * pa2pb_z_yz[j] * fl1_fx + pa2pb_zzz_yz[j]);

                t_zzz_zz[j] = fl_s_0_0 * (2.25 * pa_z[j] * fl2_fx + 1.5 * pb_z[j] * fl2_fx + 0.5 * pa_zzz[j] * fl1_fx + 3.0 * pa2pb_zz_z[j] * fl1_fx + 1.5 * pa2pb_z_zz[j] * fl1_fx + pa2pb_zzz_zz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForDG_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDG_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDG_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDG_30_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDG_40_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDG_50_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDG_60_70(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDG_70_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForDG_80_90(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForDG_0_10(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xx = paDistances.data(9 * idx + 3);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(306 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(306 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(306 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(306 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(306 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(306 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_xx_xx = pa2pbDistances.data(306 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(306 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(306 * idx + 107);

            auto pa2pb_xx_yy = pa2pbDistances.data(306 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(306 * idx + 109);

            auto pa2pb_xx_zz = pa2pbDistances.data(306 * idx + 110);

            auto pa2pb_xx_xxxx = pa2pbDistances.data(306 * idx + 121);

            auto pa2pb_xx_xxxy = pa2pbDistances.data(306 * idx + 122);

            auto pa2pb_xx_xxxz = pa2pbDistances.data(306 * idx + 123);

            auto pa2pb_xx_xxyy = pa2pbDistances.data(306 * idx + 124);

            auto pa2pb_xx_xxyz = pa2pbDistances.data(306 * idx + 125);

            auto pa2pb_xx_xxzz = pa2pbDistances.data(306 * idx + 126);

            auto pa2pb_xx_xyyy = pa2pbDistances.data(306 * idx + 127);

            auto pa2pb_xx_xyyz = pa2pbDistances.data(306 * idx + 128);

            auto pa2pb_xx_xyzz = pa2pbDistances.data(306 * idx + 129);

            auto pa2pb_xx_xzzz = pa2pbDistances.data(306 * idx + 130);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_xxxx = primBuffer.data(90 * idx);

            auto t_xx_xxxy = primBuffer.data(90 * idx + 1);

            auto t_xx_xxxz = primBuffer.data(90 * idx + 2);

            auto t_xx_xxyy = primBuffer.data(90 * idx + 3);

            auto t_xx_xxyz = primBuffer.data(90 * idx + 4);

            auto t_xx_xxzz = primBuffer.data(90 * idx + 5);

            auto t_xx_xyyy = primBuffer.data(90 * idx + 6);

            auto t_xx_xyyz = primBuffer.data(90 * idx + 7);

            auto t_xx_xyzz = primBuffer.data(90 * idx + 8);

            auto t_xx_xzzz = primBuffer.data(90 * idx + 9);

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, pa2pb_x_xyy, \
                                     pa2pb_x_xyz, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, \
                                     pa2pb_x_z, pa2pb_x_zzz, pa2pb_xx_xx, pa2pb_xx_xxxx, pa2pb_xx_xxxy, \
                                     pa2pb_xx_xxxz, pa2pb_xx_xxyy, pa2pb_xx_xxyz, pa2pb_xx_xxzz, pa2pb_xx_xy, \
                                     pa2pb_xx_xyyy, pa2pb_xx_xyyz, pa2pb_xx_xyzz, pa2pb_xx_xz, pa2pb_xx_xzzz, \
                                     pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa_xx, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, \
                                     pb_xxyz, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, pb_xzzz, pb_yy, pb_yz, pb_zz, \
                                     s_0_0, t_xx_xxxx, t_xx_xxxy, t_xx_xxxz, t_xx_xxyy, t_xx_xxyz, t_xx_xxzz, \
                                     t_xx_xyyy, t_xx_xyyz, t_xx_xyzz, t_xx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xx_xxxx[j] = fl_s_0_0 * (1.875 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 6.0 * pa2pb_x_x[j] * fl2_fx + 4.5 * pb_xx[j] * fl2_fx + 3.0 * pa2pb_xx_xx[j] * fl1_fx + 4.0 * pa2pb_x_xxx[j] * fl1_fx + 0.5 * pb_xxxx[j] * fl1_fx + pa2pb_xx_xxxx[j]);

                t_xx_xxxy[j] = fl_s_0_0 * (1.5 * pa2pb_x_y[j] * fl2_fx + 2.25 * pb_xy[j] * fl2_fx + 1.5 * pa2pb_xx_xy[j] * fl1_fx + 3.0 * pa2pb_x_xxy[j] * fl1_fx + 0.5 * pb_xxxy[j] * fl1_fx + pa2pb_xx_xxxy[j]);

                t_xx_xxxz[j] = fl_s_0_0 * (1.5 * pa2pb_x_z[j] * fl2_fx + 2.25 * pb_xz[j] * fl2_fx + 1.5 * pa2pb_xx_xz[j] * fl1_fx + 3.0 * pa2pb_x_xxz[j] * fl1_fx + 0.5 * pb_xxxz[j] * fl1_fx + pa2pb_xx_xxxz[j]);

                t_xx_xxyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + pa2pb_x_x[j] * fl2_fx + 0.75 * pb_yy[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_xx_xx[j] * fl1_fx + 0.5 * pa2pb_xx_yy[j] * fl1_fx + 2.0 * pa2pb_x_xyy[j] * fl1_fx + 0.5 * pb_xxyy[j] * fl1_fx + pa2pb_xx_xxyy[j]);

                t_xx_xxyz[j] = fl_s_0_0 * (0.75 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_xx_yz[j] * fl1_fx + 2.0 * pa2pb_x_xyz[j] * fl1_fx + 0.5 * pb_xxyz[j] * fl1_fx + pa2pb_xx_xxyz[j]);

                t_xx_xxzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + pa2pb_x_x[j] * fl2_fx + 0.75 * pb_zz[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_xx_xx[j] * fl1_fx + 0.5 * pa2pb_xx_zz[j] * fl1_fx + 2.0 * pa2pb_x_xzz[j] * fl1_fx + 0.5 * pb_xxzz[j] * fl1_fx + pa2pb_xx_xxzz[j]);

                t_xx_xyyy[j] = fl_s_0_0 * (1.5 * pa2pb_x_y[j] * fl2_fx + 0.75 * pb_xy[j] * fl2_fx + 1.5 * pa2pb_xx_xy[j] * fl1_fx + pa2pb_x_yyy[j] * fl1_fx + 0.5 * pb_xyyy[j] * fl1_fx + pa2pb_xx_xyyy[j]);

                t_xx_xyyz[j] = fl_s_0_0 * (0.5 * pa2pb_x_z[j] * fl2_fx + 0.25 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xx_xz[j] * fl1_fx + pa2pb_x_yyz[j] * fl1_fx + 0.5 * pb_xyyz[j] * fl1_fx + pa2pb_xx_xyyz[j]);

                t_xx_xyzz[j] = fl_s_0_0 * (0.5 * pa2pb_x_y[j] * fl2_fx + 0.25 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xx_xy[j] * fl1_fx + pa2pb_x_yzz[j] * fl1_fx + 0.5 * pb_xyzz[j] * fl1_fx + pa2pb_xx_xyzz[j]);

                t_xx_xzzz[j] = fl_s_0_0 * (1.5 * pa2pb_x_z[j] * fl2_fx + 0.75 * pb_xz[j] * fl2_fx + 1.5 * pa2pb_xx_xz[j] * fl1_fx + pa2pb_x_zzz[j] * fl1_fx + 0.5 * pb_xzzz[j] * fl1_fx + pa2pb_xx_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG_10_20(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (10,20)

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

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(306 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(306 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(306 * idx + 45);

            auto pa2pb_y_xyy = pa2pbDistances.data(306 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(306 * idx + 47);

            auto pa2pb_xx_yy = pa2pbDistances.data(306 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(306 * idx + 109);

            auto pa2pb_xx_zz = pa2pbDistances.data(306 * idx + 110);

            auto pa2pb_xx_yyyy = pa2pbDistances.data(306 * idx + 131);

            auto pa2pb_xx_yyyz = pa2pbDistances.data(306 * idx + 132);

            auto pa2pb_xx_yyzz = pa2pbDistances.data(306 * idx + 133);

            auto pa2pb_xx_yzzz = pa2pbDistances.data(306 * idx + 134);

            auto pa2pb_xx_zzzz = pa2pbDistances.data(306 * idx + 135);

            auto pa2pb_xy_xx = pa2pbDistances.data(306 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(306 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(306 * idx + 141);

            auto pa2pb_xy_yy = pa2pbDistances.data(306 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(306 * idx + 143);

            auto pa2pb_xy_xxxx = pa2pbDistances.data(306 * idx + 155);

            auto pa2pb_xy_xxxy = pa2pbDistances.data(306 * idx + 156);

            auto pa2pb_xy_xxxz = pa2pbDistances.data(306 * idx + 157);

            auto pa2pb_xy_xxyy = pa2pbDistances.data(306 * idx + 158);

            auto pa2pb_xy_xxyz = pa2pbDistances.data(306 * idx + 159);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_yyyy = primBuffer.data(90 * idx + 10);

            auto t_xx_yyyz = primBuffer.data(90 * idx + 11);

            auto t_xx_yyzz = primBuffer.data(90 * idx + 12);

            auto t_xx_yzzz = primBuffer.data(90 * idx + 13);

            auto t_xx_zzzz = primBuffer.data(90 * idx + 14);

            auto t_xy_xxxx = primBuffer.data(90 * idx + 15);

            auto t_xy_xxxy = primBuffer.data(90 * idx + 16);

            auto t_xy_xxxz = primBuffer.data(90 * idx + 17);

            auto t_xy_xxyy = primBuffer.data(90 * idx + 18);

            auto t_xy_xxyz = primBuffer.data(90 * idx + 19);

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, pa2pb_x_y, \
                                     pa2pb_x_z, pa2pb_xx_yy, pa2pb_xx_yyyy, pa2pb_xx_yyyz, pa2pb_xx_yyzz, \
                                     pa2pb_xx_yz, pa2pb_xx_yzzz, pa2pb_xx_zz, pa2pb_xx_zzzz, pa2pb_xy_xx, \
                                     pa2pb_xy_xxxx, pa2pb_xy_xxxy, pa2pb_xy_xxxz, pa2pb_xy_xxyy, pa2pb_xy_xxyz, \
                                     pa2pb_xy_xy, pa2pb_xy_xz, pa2pb_xy_yy, pa2pb_xy_yz, pa2pb_y_x, pa2pb_y_xxx, \
                                     pa2pb_y_xxy, pa2pb_y_xxz, pa2pb_y_xyy, pa2pb_y_xyz, pa2pb_y_y, pa2pb_y_z, pa_xx, \
                                     pa_xy, pb_xx, pb_xy, pb_xz, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_zz, \
                                     pb_zzzz, s_0_0, t_xx_yyyy, t_xx_yyyz, t_xx_yyzz, t_xx_yzzz, t_xx_zzzz, t_xy_xxxx, \
                                     t_xy_xxxy, t_xy_xxxz, t_xy_xxyy, t_xy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xx_yyyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 1.5 * pb_yy[j] * fl2_fx + 3.0 * pa2pb_xx_yy[j] * fl1_fx + 0.5 * pb_yyyy[j] * fl1_fx + pa2pb_xx_yyyy[j]);

                t_xx_yyyz[j] = fl_s_0_0 * (0.75 * pb_yz[j] * fl2_fx + 1.5 * pa2pb_xx_yz[j] * fl1_fx + 0.5 * pb_yyyz[j] * fl1_fx + pa2pb_xx_yyyz[j]);

                t_xx_yyzz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_xx_yy[j] * fl1_fx + 0.5 * pa2pb_xx_zz[j] * fl1_fx + 0.5 * pb_yyzz[j] * fl1_fx + pa2pb_xx_yyzz[j]);

                t_xx_yzzz[j] = fl_s_0_0 * (0.75 * pb_yz[j] * fl2_fx + 1.5 * pa2pb_xx_yz[j] * fl1_fx + 0.5 * pb_yzzz[j] * fl1_fx + pa2pb_xx_yzzz[j]);

                t_xx_zzzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 1.5 * pb_zz[j] * fl2_fx + 3.0 * pa2pb_xx_zz[j] * fl1_fx + 0.5 * pb_zzzz[j] * fl1_fx + pa2pb_xx_zzzz[j]);

                t_xy_xxxx[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 3.0 * pa2pb_y_x[j] * fl2_fx + 3.0 * pa2pb_xy_xx[j] * fl1_fx + 2.0 * pa2pb_y_xxx[j] * fl1_fx + pa2pb_xy_xxxx[j]);

                t_xy_xxxy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pb_xx[j] * fl2_fx + 1.5 * pa2pb_xy_xy[j] * fl1_fx + 0.5 * pa2pb_x_xxx[j] * fl1_fx + 1.5 * pa2pb_y_xxy[j] * fl1_fx + pa2pb_xy_xxxy[j]);

                t_xy_xxxz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl2_fx + 1.5 * pa2pb_xy_xz[j] * fl1_fx + 1.5 * pa2pb_y_xxz[j] * fl1_fx + pa2pb_xy_xxxz[j]);

                t_xy_xxyy[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_y_x[j] * fl2_fx + pb_xy[j] * fl2_fx + 0.5 * pa2pb_xy_xx[j] * fl1_fx + 0.5 * pa2pb_xy_yy[j] * fl1_fx + pa2pb_x_xxy[j] * fl1_fx + pa2pb_y_xyy[j] * fl1_fx + pa2pb_xy_xxyy[j]);

                t_xy_xxyz[j] = fl_s_0_0 * (0.25 * pa2pb_x_z[j] * fl2_fx + 0.5 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xy_yz[j] * fl1_fx + 0.5 * pa2pb_x_xxz[j] * fl1_fx + pa2pb_y_xyz[j] * fl1_fx + pa2pb_xy_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG_20_30(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (20,30)

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

            auto pa_xy = paDistances.data(9 * idx + 4);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_x_xyy = pa2pbDistances.data(306 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(306 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(306 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(306 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(306 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(306 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 36);

            auto pa2pb_y_xzz = pa2pbDistances.data(306 * idx + 48);

            auto pa2pb_y_yyy = pa2pbDistances.data(306 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(306 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(306 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(306 * idx + 52);

            auto pa2pb_xy_xx = pa2pbDistances.data(306 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(306 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(306 * idx + 141);

            auto pa2pb_xy_yy = pa2pbDistances.data(306 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(306 * idx + 143);

            auto pa2pb_xy_zz = pa2pbDistances.data(306 * idx + 144);

            auto pa2pb_xy_xxzz = pa2pbDistances.data(306 * idx + 160);

            auto pa2pb_xy_xyyy = pa2pbDistances.data(306 * idx + 161);

            auto pa2pb_xy_xyyz = pa2pbDistances.data(306 * idx + 162);

            auto pa2pb_xy_xyzz = pa2pbDistances.data(306 * idx + 163);

            auto pa2pb_xy_xzzz = pa2pbDistances.data(306 * idx + 164);

            auto pa2pb_xy_yyyy = pa2pbDistances.data(306 * idx + 165);

            auto pa2pb_xy_yyyz = pa2pbDistances.data(306 * idx + 166);

            auto pa2pb_xy_yyzz = pa2pbDistances.data(306 * idx + 167);

            auto pa2pb_xy_yzzz = pa2pbDistances.data(306 * idx + 168);

            auto pa2pb_xy_zzzz = pa2pbDistances.data(306 * idx + 169);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xy_xxzz = primBuffer.data(90 * idx + 20);

            auto t_xy_xyyy = primBuffer.data(90 * idx + 21);

            auto t_xy_xyyz = primBuffer.data(90 * idx + 22);

            auto t_xy_xyzz = primBuffer.data(90 * idx + 23);

            auto t_xy_xzzz = primBuffer.data(90 * idx + 24);

            auto t_xy_yyyy = primBuffer.data(90 * idx + 25);

            auto t_xy_yyyz = primBuffer.data(90 * idx + 26);

            auto t_xy_yyzz = primBuffer.data(90 * idx + 27);

            auto t_xy_yzzz = primBuffer.data(90 * idx + 28);

            auto t_xy_zzzz = primBuffer.data(90 * idx + 29);

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_xzz, pa2pb_x_y, \
                                     pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, pa2pb_x_z, pa2pb_x_zzz, pa2pb_xy_xx, \
                                     pa2pb_xy_xxzz, pa2pb_xy_xy, pa2pb_xy_xyyy, pa2pb_xy_xyyz, pa2pb_xy_xyzz, \
                                     pa2pb_xy_xz, pa2pb_xy_xzzz, pa2pb_xy_yy, pa2pb_xy_yyyy, pa2pb_xy_yyyz, \
                                     pa2pb_xy_yyzz, pa2pb_xy_yz, pa2pb_xy_yzzz, pa2pb_xy_zz, pa2pb_xy_zzzz, pa2pb_y_x, \
                                     pa2pb_y_xzz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, pa2pb_y_z, \
                                     pa2pb_y_zzz, pa_xy, pb_yy, pb_yz, pb_zz, s_0_0, t_xy_xxzz, t_xy_xyyy, t_xy_xyyz, \
                                     t_xy_xyzz, t_xy_xzzz, t_xy_yyyy, t_xy_yyyz, t_xy_yyzz, t_xy_yzzz, t_xy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xy_xxzz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa2pb_xy_xx[j] * fl1_fx + 0.5 * pa2pb_xy_zz[j] * fl1_fx + pa2pb_y_xzz[j] * fl1_fx + pa2pb_xy_xxzz[j]);

                t_xy_xyyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pb_yy[j] * fl2_fx + 1.5 * pa2pb_xy_xy[j] * fl1_fx + 1.5 * pa2pb_x_xyy[j] * fl1_fx + 0.5 * pa2pb_y_yyy[j] * fl1_fx + pa2pb_xy_xyyy[j]);

                t_xy_xyyz[j] = fl_s_0_0 * (0.25 * pa2pb_y_z[j] * fl2_fx + 0.5 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_xy_xz[j] * fl1_fx + pa2pb_x_xyz[j] * fl1_fx + 0.5 * pa2pb_y_yyz[j] * fl1_fx + pa2pb_xy_xyyz[j]);

                t_xy_xyzz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa2pb_x_x[j] * fl2_fx + 0.25 * pa2pb_y_y[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_xy_xy[j] * fl1_fx + 0.5 * pa2pb_x_xzz[j] * fl1_fx + 0.5 * pa2pb_y_yzz[j] * fl1_fx + pa2pb_xy_xyzz[j]);

                t_xy_xzzz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl2_fx + 1.5 * pa2pb_xy_xz[j] * fl1_fx + 0.5 * pa2pb_y_zzz[j] * fl1_fx + pa2pb_xy_xzzz[j]);

                t_xy_yyyy[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 3.0 * pa2pb_x_y[j] * fl2_fx + 3.0 * pa2pb_xy_yy[j] * fl1_fx + 2.0 * pa2pb_x_yyy[j] * fl1_fx + pa2pb_xy_yyyy[j]);

                t_xy_yyyz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl2_fx + 1.5 * pa2pb_xy_yz[j] * fl1_fx + 1.5 * pa2pb_x_yyz[j] * fl1_fx + pa2pb_xy_yyyz[j]);

                t_xy_yyzz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_xy_yy[j] * fl1_fx + 0.5 * pa2pb_xy_zz[j] * fl1_fx + pa2pb_x_yzz[j] * fl1_fx + pa2pb_xy_yyzz[j]);

                t_xy_yzzz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl2_fx + 1.5 * pa2pb_xy_yz[j] * fl1_fx + 0.5 * pa2pb_x_zzz[j] * fl1_fx + pa2pb_xy_yzzz[j]);

                t_xy_zzzz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 3.0 * pa2pb_xy_zz[j] * fl1_fx + pa2pb_xy_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG_30_40(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (30,40)

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

            auto pa_xz = paDistances.data(9 * idx + 5);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(306 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(306 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(306 * idx + 14);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(306 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(306 * idx + 78);

            auto pa2pb_z_xxz = pa2pbDistances.data(306 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(306 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(306 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(306 * idx + 82);

            auto pa2pb_z_yyy = pa2pbDistances.data(306 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(306 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(306 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(306 * idx + 86);

            auto pa2pb_xz_xx = pa2pbDistances.data(306 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(306 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(306 * idx + 175);

            auto pa2pb_xz_yy = pa2pbDistances.data(306 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(306 * idx + 177);

            auto pa2pb_xz_zz = pa2pbDistances.data(306 * idx + 178);

            auto pa2pb_xz_xxxx = pa2pbDistances.data(306 * idx + 189);

            auto pa2pb_xz_xxxy = pa2pbDistances.data(306 * idx + 190);

            auto pa2pb_xz_xxxz = pa2pbDistances.data(306 * idx + 191);

            auto pa2pb_xz_xxyy = pa2pbDistances.data(306 * idx + 192);

            auto pa2pb_xz_xxyz = pa2pbDistances.data(306 * idx + 193);

            auto pa2pb_xz_xxzz = pa2pbDistances.data(306 * idx + 194);

            auto pa2pb_xz_xyyy = pa2pbDistances.data(306 * idx + 195);

            auto pa2pb_xz_xyyz = pa2pbDistances.data(306 * idx + 196);

            auto pa2pb_xz_xyzz = pa2pbDistances.data(306 * idx + 197);

            auto pa2pb_xz_xzzz = pa2pbDistances.data(306 * idx + 198);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xz_xxxx = primBuffer.data(90 * idx + 30);

            auto t_xz_xxxy = primBuffer.data(90 * idx + 31);

            auto t_xz_xxxz = primBuffer.data(90 * idx + 32);

            auto t_xz_xxyy = primBuffer.data(90 * idx + 33);

            auto t_xz_xxyz = primBuffer.data(90 * idx + 34);

            auto t_xz_xxzz = primBuffer.data(90 * idx + 35);

            auto t_xz_xyyy = primBuffer.data(90 * idx + 36);

            auto t_xz_xyyz = primBuffer.data(90 * idx + 37);

            auto t_xz_xyzz = primBuffer.data(90 * idx + 38);

            auto t_xz_xzzz = primBuffer.data(90 * idx + 39);

            // Batch of Integrals (30,40)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, pa2pb_x_xyy, \
                                     pa2pb_x_xyz, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_z, pa2pb_xz_xx, pa2pb_xz_xxxx, \
                                     pa2pb_xz_xxxy, pa2pb_xz_xxxz, pa2pb_xz_xxyy, pa2pb_xz_xxyz, pa2pb_xz_xxzz, \
                                     pa2pb_xz_xy, pa2pb_xz_xyyy, pa2pb_xz_xyyz, pa2pb_xz_xyzz, pa2pb_xz_xz, \
                                     pa2pb_xz_xzzz, pa2pb_xz_yy, pa2pb_xz_yz, pa2pb_xz_zz, pa2pb_z_x, pa2pb_z_xxx, \
                                     pa2pb_z_xxy, pa2pb_z_xxz, pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_z_y, \
                                     pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, pa2pb_z_z, pa2pb_z_zzz, pa_xz, pb_xx, pb_xy, \
                                     pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, t_xz_xxxx, t_xz_xxxy, t_xz_xxxz, t_xz_xxyy, \
                                     t_xz_xxyz, t_xz_xxzz, t_xz_xyyy, t_xz_xyyz, t_xz_xyzz, t_xz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xz_xxxx[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 3.0 * pa2pb_z_x[j] * fl2_fx + 3.0 * pa2pb_xz_xx[j] * fl1_fx + 2.0 * pa2pb_z_xxx[j] * fl1_fx + pa2pb_xz_xxxx[j]);

                t_xz_xxxy[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl2_fx + 1.5 * pa2pb_xz_xy[j] * fl1_fx + 1.5 * pa2pb_z_xxy[j] * fl1_fx + pa2pb_xz_xxxy[j]);

                t_xz_xxxz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.75 * pb_xx[j] * fl2_fx + 1.5 * pa2pb_xz_xz[j] * fl1_fx + 0.5 * pa2pb_x_xxx[j] * fl1_fx + 1.5 * pa2pb_z_xxz[j] * fl1_fx + pa2pb_xz_xxxz[j]);

                t_xz_xxyy[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa2pb_xz_xx[j] * fl1_fx + 0.5 * pa2pb_xz_yy[j] * fl1_fx + pa2pb_z_xyy[j] * fl1_fx + pa2pb_xz_xxyy[j]);

                t_xz_xxyz[j] = fl_s_0_0 * (0.25 * pa2pb_x_y[j] * fl2_fx + 0.5 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xz_yz[j] * fl1_fx + 0.5 * pa2pb_x_xxy[j] * fl1_fx + pa2pb_z_xyz[j] * fl1_fx + pa2pb_xz_xxyz[j]);

                t_xz_xxzz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa2pb_z_x[j] * fl2_fx + pb_xz[j] * fl2_fx + 0.5 * pa2pb_xz_xx[j] * fl1_fx + 0.5 * pa2pb_xz_zz[j] * fl1_fx + pa2pb_x_xxz[j] * fl1_fx + pa2pb_z_xzz[j] * fl1_fx + pa2pb_xz_xxzz[j]);

                t_xz_xyyy[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl2_fx + 1.5 * pa2pb_xz_xy[j] * fl1_fx + 0.5 * pa2pb_z_yyy[j] * fl1_fx + pa2pb_xz_xyyy[j]);

                t_xz_xyyz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa2pb_x_x[j] * fl2_fx + 0.25 * pa2pb_z_z[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_xz_xz[j] * fl1_fx + 0.5 * pa2pb_x_xyy[j] * fl1_fx + 0.5 * pa2pb_z_yyz[j] * fl1_fx + pa2pb_xz_xyyz[j]);

                t_xz_xyzz[j] = fl_s_0_0 * (0.25 * pa2pb_z_y[j] * fl2_fx + 0.5 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_xz_xy[j] * fl1_fx + pa2pb_x_xyz[j] * fl1_fx + 0.5 * pa2pb_z_yzz[j] * fl1_fx + pa2pb_xz_xyzz[j]);

                t_xz_xzzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.75 * pb_zz[j] * fl2_fx + 1.5 * pa2pb_xz_xz[j] * fl1_fx + 1.5 * pa2pb_x_xzz[j] * fl1_fx + 0.5 * pa2pb_z_zzz[j] * fl1_fx + pa2pb_xz_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG_40_50(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (40,50)

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

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_x_yyy = pa2pbDistances.data(306 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(306 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(306 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(306 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(306 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(306 * idx + 45);

            auto pa2pb_xz_yy = pa2pbDistances.data(306 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(306 * idx + 177);

            auto pa2pb_xz_zz = pa2pbDistances.data(306 * idx + 178);

            auto pa2pb_xz_yyyy = pa2pbDistances.data(306 * idx + 199);

            auto pa2pb_xz_yyyz = pa2pbDistances.data(306 * idx + 200);

            auto pa2pb_xz_yyzz = pa2pbDistances.data(306 * idx + 201);

            auto pa2pb_xz_yzzz = pa2pbDistances.data(306 * idx + 202);

            auto pa2pb_xz_zzzz = pa2pbDistances.data(306 * idx + 203);

            auto pa2pb_yy_xx = pa2pbDistances.data(306 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(306 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(306 * idx + 209);

            auto pa2pb_yy_yy = pa2pbDistances.data(306 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(306 * idx + 211);

            auto pa2pb_yy_xxxx = pa2pbDistances.data(306 * idx + 223);

            auto pa2pb_yy_xxxy = pa2pbDistances.data(306 * idx + 224);

            auto pa2pb_yy_xxxz = pa2pbDistances.data(306 * idx + 225);

            auto pa2pb_yy_xxyy = pa2pbDistances.data(306 * idx + 226);

            auto pa2pb_yy_xxyz = pa2pbDistances.data(306 * idx + 227);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xz_yyyy = primBuffer.data(90 * idx + 40);

            auto t_xz_yyyz = primBuffer.data(90 * idx + 41);

            auto t_xz_yyzz = primBuffer.data(90 * idx + 42);

            auto t_xz_yzzz = primBuffer.data(90 * idx + 43);

            auto t_xz_zzzz = primBuffer.data(90 * idx + 44);

            auto t_yy_xxxx = primBuffer.data(90 * idx + 45);

            auto t_yy_xxxy = primBuffer.data(90 * idx + 46);

            auto t_yy_xxxz = primBuffer.data(90 * idx + 47);

            auto t_yy_xxyy = primBuffer.data(90 * idx + 48);

            auto t_yy_xxyz = primBuffer.data(90 * idx + 49);

            // Batch of Integrals (40,50)

            #pragma omp simd aligned(fx, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, pa2pb_x_z, \
                                     pa2pb_x_zzz, pa2pb_xz_yy, pa2pb_xz_yyyy, pa2pb_xz_yyyz, pa2pb_xz_yyzz, \
                                     pa2pb_xz_yz, pa2pb_xz_yzzz, pa2pb_xz_zz, pa2pb_xz_zzzz, pa2pb_y_x, pa2pb_y_xxx, \
                                     pa2pb_y_xxy, pa2pb_y_xxz, pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xxxx, \
                                     pa2pb_yy_xxxy, pa2pb_yy_xxxz, pa2pb_yy_xxyy, pa2pb_yy_xxyz, pa2pb_yy_xy, \
                                     pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, pa_xz, pa_yy, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_yy, pb_yz, s_0_0, t_xz_yyyy, t_xz_yyyz, t_xz_yyzz, \
                                     t_xz_yzzz, t_xz_zzzz, t_yy_xxxx, t_yy_xxxy, t_yy_xxxz, t_yy_xxyy, t_yy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xz_yyyy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 3.0 * pa2pb_xz_yy[j] * fl1_fx + pa2pb_xz_yyyy[j]);

                t_xz_yyyz[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl2_fx + 1.5 * pa2pb_xz_yz[j] * fl1_fx + 0.5 * pa2pb_x_yyy[j] * fl1_fx + pa2pb_xz_yyyz[j]);

                t_xz_yyzz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa2pb_xz_yy[j] * fl1_fx + 0.5 * pa2pb_xz_zz[j] * fl1_fx + pa2pb_x_yyz[j] * fl1_fx + pa2pb_xz_yyzz[j]);

                t_xz_yzzz[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl2_fx + 1.5 * pa2pb_xz_yz[j] * fl1_fx + 1.5 * pa2pb_x_yzz[j] * fl1_fx + pa2pb_xz_yzzz[j]);

                t_xz_zzzz[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 3.0 * pa2pb_x_z[j] * fl2_fx + 3.0 * pa2pb_xz_zz[j] * fl1_fx + 2.0 * pa2pb_x_zzz[j] * fl1_fx + pa2pb_xz_zzzz[j]);

                t_yy_xxxx[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 1.5 * pb_xx[j] * fl2_fx + 3.0 * pa2pb_yy_xx[j] * fl1_fx + 0.5 * pb_xxxx[j] * fl1_fx + pa2pb_yy_xxxx[j]);

                t_yy_xxxy[j] = fl_s_0_0 * (1.5 * pa2pb_y_x[j] * fl2_fx + 0.75 * pb_xy[j] * fl2_fx + 1.5 * pa2pb_yy_xy[j] * fl1_fx + pa2pb_y_xxx[j] * fl1_fx + 0.5 * pb_xxxy[j] * fl1_fx + pa2pb_yy_xxxy[j]);

                t_yy_xxxz[j] = fl_s_0_0 * (0.75 * pb_xz[j] * fl2_fx + 1.5 * pa2pb_yy_xz[j] * fl1_fx + 0.5 * pb_xxxz[j] * fl1_fx + pa2pb_yy_xxxz[j]);

                t_yy_xxyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + pa2pb_y_y[j] * fl2_fx + 0.75 * pb_xx[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_yy_xx[j] * fl1_fx + 0.5 * pa2pb_yy_yy[j] * fl1_fx + 2.0 * pa2pb_y_xxy[j] * fl1_fx + 0.5 * pb_xxyy[j] * fl1_fx + pa2pb_yy_xxyy[j]);

                t_yy_xxyz[j] = fl_s_0_0 * (0.5 * pa2pb_y_z[j] * fl2_fx + 0.25 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_yy_yz[j] * fl1_fx + pa2pb_y_xxz[j] * fl1_fx + 0.5 * pb_xxyz[j] * fl1_fx + pa2pb_yy_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG_50_60(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (50,60)

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

            auto pa_yy = paDistances.data(9 * idx + 6);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 36);

            auto pa2pb_y_xyy = pa2pbDistances.data(306 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(306 * idx + 47);

            auto pa2pb_y_xzz = pa2pbDistances.data(306 * idx + 48);

            auto pa2pb_y_yyy = pa2pbDistances.data(306 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(306 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(306 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(306 * idx + 52);

            auto pa2pb_yy_xx = pa2pbDistances.data(306 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(306 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(306 * idx + 209);

            auto pa2pb_yy_yy = pa2pbDistances.data(306 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(306 * idx + 211);

            auto pa2pb_yy_zz = pa2pbDistances.data(306 * idx + 212);

            auto pa2pb_yy_xxzz = pa2pbDistances.data(306 * idx + 228);

            auto pa2pb_yy_xyyy = pa2pbDistances.data(306 * idx + 229);

            auto pa2pb_yy_xyyz = pa2pbDistances.data(306 * idx + 230);

            auto pa2pb_yy_xyzz = pa2pbDistances.data(306 * idx + 231);

            auto pa2pb_yy_xzzz = pa2pbDistances.data(306 * idx + 232);

            auto pa2pb_yy_yyyy = pa2pbDistances.data(306 * idx + 233);

            auto pa2pb_yy_yyyz = pa2pbDistances.data(306 * idx + 234);

            auto pa2pb_yy_yyzz = pa2pbDistances.data(306 * idx + 235);

            auto pa2pb_yy_yzzz = pa2pbDistances.data(306 * idx + 236);

            auto pa2pb_yy_zzzz = pa2pbDistances.data(306 * idx + 237);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yy_xxzz = primBuffer.data(90 * idx + 50);

            auto t_yy_xyyy = primBuffer.data(90 * idx + 51);

            auto t_yy_xyyz = primBuffer.data(90 * idx + 52);

            auto t_yy_xyzz = primBuffer.data(90 * idx + 53);

            auto t_yy_xzzz = primBuffer.data(90 * idx + 54);

            auto t_yy_yyyy = primBuffer.data(90 * idx + 55);

            auto t_yy_yyyz = primBuffer.data(90 * idx + 56);

            auto t_yy_yyzz = primBuffer.data(90 * idx + 57);

            auto t_yy_yzzz = primBuffer.data(90 * idx + 58);

            auto t_yy_zzzz = primBuffer.data(90 * idx + 59);

            // Batch of Integrals (50,60)

            #pragma omp simd aligned(fx, pa2pb_y_x, pa2pb_y_xyy, pa2pb_y_xyz, pa2pb_y_xzz, pa2pb_y_y, \
                                     pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, pa2pb_yy_xx, \
                                     pa2pb_yy_xxzz, pa2pb_yy_xy, pa2pb_yy_xyyy, pa2pb_yy_xyyz, pa2pb_yy_xyzz, \
                                     pa2pb_yy_xz, pa2pb_yy_xzzz, pa2pb_yy_yy, pa2pb_yy_yyyy, pa2pb_yy_yyyz, \
                                     pa2pb_yy_yyzz, pa2pb_yy_yz, pa2pb_yy_yzzz, pa2pb_yy_zz, pa2pb_yy_zzzz, pa_yy, pb_xx, \
                                     pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, pb_xzzz, pb_yy, pb_yyyy, pb_yyyz, \
                                     pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, s_0_0, t_yy_xxzz, t_yy_xyyy, t_yy_xyyz, \
                                     t_yy_xyzz, t_yy_xzzz, t_yy_yyyy, t_yy_yyyz, t_yy_yyzz, t_yy_yzzz, t_yy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yy_xxzz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_yy_xx[j] * fl1_fx + 0.5 * pa2pb_yy_zz[j] * fl1_fx + 0.5 * pb_xxzz[j] * fl1_fx + pa2pb_yy_xxzz[j]);

                t_yy_xyyy[j] = fl_s_0_0 * (1.5 * pa2pb_y_x[j] * fl2_fx + 2.25 * pb_xy[j] * fl2_fx + 1.5 * pa2pb_yy_xy[j] * fl1_fx + 3.0 * pa2pb_y_xyy[j] * fl1_fx + 0.5 * pb_xyyy[j] * fl1_fx + pa2pb_yy_xyyy[j]);

                t_yy_xyyz[j] = fl_s_0_0 * (0.75 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_yy_xz[j] * fl1_fx + 2.0 * pa2pb_y_xyz[j] * fl1_fx + 0.5 * pb_xyyz[j] * fl1_fx + pa2pb_yy_xyyz[j]);

                t_yy_xyzz[j] = fl_s_0_0 * (0.5 * pa2pb_y_x[j] * fl2_fx + 0.25 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_yy_xy[j] * fl1_fx + pa2pb_y_xzz[j] * fl1_fx + 0.5 * pb_xyzz[j] * fl1_fx + pa2pb_yy_xyzz[j]);

                t_yy_xzzz[j] = fl_s_0_0 * (0.75 * pb_xz[j] * fl2_fx + 1.5 * pa2pb_yy_xz[j] * fl1_fx + 0.5 * pb_xzzz[j] * fl1_fx + pa2pb_yy_xzzz[j]);

                t_yy_yyyy[j] = fl_s_0_0 * (1.875 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 6.0 * pa2pb_y_y[j] * fl2_fx + 4.5 * pb_yy[j] * fl2_fx + 3.0 * pa2pb_yy_yy[j] * fl1_fx + 4.0 * pa2pb_y_yyy[j] * fl1_fx + 0.5 * pb_yyyy[j] * fl1_fx + pa2pb_yy_yyyy[j]);

                t_yy_yyyz[j] = fl_s_0_0 * (1.5 * pa2pb_y_z[j] * fl2_fx + 2.25 * pb_yz[j] * fl2_fx + 1.5 * pa2pb_yy_yz[j] * fl1_fx + 3.0 * pa2pb_y_yyz[j] * fl1_fx + 0.5 * pb_yyyz[j] * fl1_fx + pa2pb_yy_yyyz[j]);

                t_yy_yyzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + pa2pb_y_y[j] * fl2_fx + 0.75 * pb_zz[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_yy_yy[j] * fl1_fx + 0.5 * pa2pb_yy_zz[j] * fl1_fx + 2.0 * pa2pb_y_yzz[j] * fl1_fx + 0.5 * pb_yyzz[j] * fl1_fx + pa2pb_yy_yyzz[j]);

                t_yy_yzzz[j] = fl_s_0_0 * (1.5 * pa2pb_y_z[j] * fl2_fx + 0.75 * pb_yz[j] * fl2_fx + 1.5 * pa2pb_yy_yz[j] * fl1_fx + pa2pb_y_zzz[j] * fl1_fx + 0.5 * pb_yzzz[j] * fl1_fx + pa2pb_yy_yzzz[j]);

                t_yy_zzzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 1.5 * pb_zz[j] * fl2_fx + 3.0 * pa2pb_yy_zz[j] * fl1_fx + 0.5 * pb_zzzz[j] * fl1_fx + pa2pb_yy_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG_60_70(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (60,70)

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

            auto pa_yz = paDistances.data(9 * idx + 7);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(306 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(306 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(306 * idx + 45);

            auto pa2pb_y_xyy = pa2pbDistances.data(306 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(306 * idx + 47);

            auto pa2pb_y_xzz = pa2pbDistances.data(306 * idx + 48);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(306 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(306 * idx + 78);

            auto pa2pb_z_xxz = pa2pbDistances.data(306 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(306 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(306 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(306 * idx + 82);

            auto pa2pb_yz_xx = pa2pbDistances.data(306 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(306 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(306 * idx + 243);

            auto pa2pb_yz_yy = pa2pbDistances.data(306 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(306 * idx + 245);

            auto pa2pb_yz_zz = pa2pbDistances.data(306 * idx + 246);

            auto pa2pb_yz_xxxx = pa2pbDistances.data(306 * idx + 257);

            auto pa2pb_yz_xxxy = pa2pbDistances.data(306 * idx + 258);

            auto pa2pb_yz_xxxz = pa2pbDistances.data(306 * idx + 259);

            auto pa2pb_yz_xxyy = pa2pbDistances.data(306 * idx + 260);

            auto pa2pb_yz_xxyz = pa2pbDistances.data(306 * idx + 261);

            auto pa2pb_yz_xxzz = pa2pbDistances.data(306 * idx + 262);

            auto pa2pb_yz_xyyy = pa2pbDistances.data(306 * idx + 263);

            auto pa2pb_yz_xyyz = pa2pbDistances.data(306 * idx + 264);

            auto pa2pb_yz_xyzz = pa2pbDistances.data(306 * idx + 265);

            auto pa2pb_yz_xzzz = pa2pbDistances.data(306 * idx + 266);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yz_xxxx = primBuffer.data(90 * idx + 60);

            auto t_yz_xxxy = primBuffer.data(90 * idx + 61);

            auto t_yz_xxxz = primBuffer.data(90 * idx + 62);

            auto t_yz_xxyy = primBuffer.data(90 * idx + 63);

            auto t_yz_xxyz = primBuffer.data(90 * idx + 64);

            auto t_yz_xxzz = primBuffer.data(90 * idx + 65);

            auto t_yz_xyyy = primBuffer.data(90 * idx + 66);

            auto t_yz_xyyz = primBuffer.data(90 * idx + 67);

            auto t_yz_xyzz = primBuffer.data(90 * idx + 68);

            auto t_yz_xzzz = primBuffer.data(90 * idx + 69);

            // Batch of Integrals (60,70)

            #pragma omp simd aligned(fx, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, pa2pb_y_xyy, \
                                     pa2pb_y_xyz, pa2pb_y_xzz, pa2pb_y_y, pa2pb_y_z, pa2pb_yz_xx, pa2pb_yz_xxxx, \
                                     pa2pb_yz_xxxy, pa2pb_yz_xxxz, pa2pb_yz_xxyy, pa2pb_yz_xxyz, pa2pb_yz_xxzz, \
                                     pa2pb_yz_xy, pa2pb_yz_xyyy, pa2pb_yz_xyyz, pa2pb_yz_xyzz, pa2pb_yz_xz, \
                                     pa2pb_yz_xzzz, pa2pb_yz_yy, pa2pb_yz_yz, pa2pb_yz_zz, pa2pb_z_x, pa2pb_z_xxx, \
                                     pa2pb_z_xxy, pa2pb_z_xxz, pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_z_y, \
                                     pa2pb_z_z, pa_yz, pb_xx, pb_xy, pb_xz, s_0_0, t_yz_xxxx, t_yz_xxxy, t_yz_xxxz, \
                                     t_yz_xxyy, t_yz_xxyz, t_yz_xxzz, t_yz_xyyy, t_yz_xyyz, t_yz_xyzz, t_yz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yz_xxxx[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 3.0 * pa2pb_yz_xx[j] * fl1_fx + pa2pb_yz_xxxx[j]);

                t_yz_xxxy[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl2_fx + 1.5 * pa2pb_yz_xy[j] * fl1_fx + 0.5 * pa2pb_z_xxx[j] * fl1_fx + pa2pb_yz_xxxy[j]);

                t_yz_xxxz[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl2_fx + 1.5 * pa2pb_yz_xz[j] * fl1_fx + 0.5 * pa2pb_y_xxx[j] * fl1_fx + pa2pb_yz_xxxz[j]);

                t_yz_xxyy[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl2_fx + 0.5 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa2pb_yz_xx[j] * fl1_fx + 0.5 * pa2pb_yz_yy[j] * fl1_fx + pa2pb_z_xxy[j] * fl1_fx + pa2pb_yz_xxyy[j]);

                t_yz_xxyz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa2pb_y_y[j] * fl2_fx + 0.25 * pa2pb_z_z[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_yz_yz[j] * fl1_fx + 0.5 * pa2pb_y_xxy[j] * fl1_fx + 0.5 * pa2pb_z_xxz[j] * fl1_fx + pa2pb_yz_xxyz[j]);

                t_yz_xxzz[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl2_fx + 0.5 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa2pb_yz_xx[j] * fl1_fx + 0.5 * pa2pb_yz_zz[j] * fl1_fx + pa2pb_y_xxz[j] * fl1_fx + pa2pb_yz_xxzz[j]);

                t_yz_xyyy[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl2_fx + 1.5 * pa2pb_yz_xy[j] * fl1_fx + 1.5 * pa2pb_z_xyy[j] * fl1_fx + pa2pb_yz_xyyy[j]);

                t_yz_xyyz[j] = fl_s_0_0 * (0.25 * pa2pb_y_x[j] * fl2_fx + 0.5 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_yz_xz[j] * fl1_fx + 0.5 * pa2pb_y_xyy[j] * fl1_fx + pa2pb_z_xyz[j] * fl1_fx + pa2pb_yz_xyyz[j]);

                t_yz_xyzz[j] = fl_s_0_0 * (0.25 * pa2pb_z_x[j] * fl2_fx + 0.5 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_yz_xy[j] * fl1_fx + pa2pb_y_xyz[j] * fl1_fx + 0.5 * pa2pb_z_xzz[j] * fl1_fx + pa2pb_yz_xyzz[j]);

                t_yz_xzzz[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl2_fx + 1.5 * pa2pb_yz_xz[j] * fl1_fx + 1.5 * pa2pb_y_xzz[j] * fl1_fx + pa2pb_yz_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG_70_80(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (70,80)

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

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 36);

            auto pa2pb_y_yyy = pa2pbDistances.data(306 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(306 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(306 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(306 * idx + 52);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(306 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(306 * idx + 78);

            auto pa2pb_z_yyy = pa2pbDistances.data(306 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(306 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(306 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(306 * idx + 86);

            auto pa2pb_yz_yy = pa2pbDistances.data(306 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(306 * idx + 245);

            auto pa2pb_yz_zz = pa2pbDistances.data(306 * idx + 246);

            auto pa2pb_yz_yyyy = pa2pbDistances.data(306 * idx + 267);

            auto pa2pb_yz_yyyz = pa2pbDistances.data(306 * idx + 268);

            auto pa2pb_yz_yyzz = pa2pbDistances.data(306 * idx + 269);

            auto pa2pb_yz_yzzz = pa2pbDistances.data(306 * idx + 270);

            auto pa2pb_yz_zzzz = pa2pbDistances.data(306 * idx + 271);

            auto pa2pb_zz_xx = pa2pbDistances.data(306 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(306 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(306 * idx + 277);

            auto pa2pb_zz_yy = pa2pbDistances.data(306 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(306 * idx + 279);

            auto pa2pb_zz_xxxx = pa2pbDistances.data(306 * idx + 291);

            auto pa2pb_zz_xxxy = pa2pbDistances.data(306 * idx + 292);

            auto pa2pb_zz_xxxz = pa2pbDistances.data(306 * idx + 293);

            auto pa2pb_zz_xxyy = pa2pbDistances.data(306 * idx + 294);

            auto pa2pb_zz_xxyz = pa2pbDistances.data(306 * idx + 295);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yz_yyyy = primBuffer.data(90 * idx + 70);

            auto t_yz_yyyz = primBuffer.data(90 * idx + 71);

            auto t_yz_yyzz = primBuffer.data(90 * idx + 72);

            auto t_yz_yzzz = primBuffer.data(90 * idx + 73);

            auto t_yz_zzzz = primBuffer.data(90 * idx + 74);

            auto t_zz_xxxx = primBuffer.data(90 * idx + 75);

            auto t_zz_xxxy = primBuffer.data(90 * idx + 76);

            auto t_zz_xxxz = primBuffer.data(90 * idx + 77);

            auto t_zz_xxyy = primBuffer.data(90 * idx + 78);

            auto t_zz_xxyz = primBuffer.data(90 * idx + 79);

            // Batch of Integrals (70,80)

            #pragma omp simd aligned(fx, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, pa2pb_y_z, \
                                     pa2pb_y_zzz, pa2pb_yz_yy, pa2pb_yz_yyyy, pa2pb_yz_yyyz, pa2pb_yz_yyzz, \
                                     pa2pb_yz_yz, pa2pb_yz_yzzz, pa2pb_yz_zz, pa2pb_yz_zzzz, pa2pb_z_x, pa2pb_z_xxx, \
                                     pa2pb_z_xxy, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, pa2pb_z_z, \
                                     pa2pb_z_zzz, pa2pb_zz_xx, pa2pb_zz_xxxx, pa2pb_zz_xxxy, pa2pb_zz_xxxz, \
                                     pa2pb_zz_xxyy, pa2pb_zz_xxyz, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, pa2pb_zz_yz, \
                                     pa_yz, pa_zz, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_yy, \
                                     pb_yz, pb_zz, s_0_0, t_yz_yyyy, t_yz_yyyz, t_yz_yyzz, t_yz_yzzz, t_yz_zzzz, \
                                     t_zz_xxxx, t_zz_xxxy, t_zz_xxxz, t_zz_xxyy, t_zz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yz_yyyy[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 3.0 * pa2pb_z_y[j] * fl2_fx + 3.0 * pa2pb_yz_yy[j] * fl1_fx + 2.0 * pa2pb_z_yyy[j] * fl1_fx + pa2pb_yz_yyyy[j]);

                t_yz_yyyz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.75 * pb_yy[j] * fl2_fx + 1.5 * pa2pb_yz_yz[j] * fl1_fx + 0.5 * pa2pb_y_yyy[j] * fl1_fx + 1.5 * pa2pb_z_yyz[j] * fl1_fx + pa2pb_yz_yyyz[j]);

                t_yz_yyzz[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl2_fx + 0.5 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa2pb_z_y[j] * fl2_fx + pb_yz[j] * fl2_fx + 0.5 * pa2pb_yz_yy[j] * fl1_fx + 0.5 * pa2pb_yz_zz[j] * fl1_fx + pa2pb_y_yyz[j] * fl1_fx + pa2pb_z_yzz[j] * fl1_fx + pa2pb_yz_yyzz[j]);

                t_yz_yzzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.75 * pb_zz[j] * fl2_fx + 1.5 * pa2pb_yz_yz[j] * fl1_fx + 1.5 * pa2pb_y_yzz[j] * fl1_fx + 0.5 * pa2pb_z_zzz[j] * fl1_fx + pa2pb_yz_yzzz[j]);

                t_yz_zzzz[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 3.0 * pa2pb_y_z[j] * fl2_fx + 3.0 * pa2pb_yz_zz[j] * fl1_fx + 2.0 * pa2pb_y_zzz[j] * fl1_fx + pa2pb_yz_zzzz[j]);

                t_zz_xxxx[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 1.5 * pb_xx[j] * fl2_fx + 3.0 * pa2pb_zz_xx[j] * fl1_fx + 0.5 * pb_xxxx[j] * fl1_fx + pa2pb_zz_xxxx[j]);

                t_zz_xxxy[j] = fl_s_0_0 * (0.75 * pb_xy[j] * fl2_fx + 1.5 * pa2pb_zz_xy[j] * fl1_fx + 0.5 * pb_xxxy[j] * fl1_fx + pa2pb_zz_xxxy[j]);

                t_zz_xxxz[j] = fl_s_0_0 * (1.5 * pa2pb_z_x[j] * fl2_fx + 0.75 * pb_xz[j] * fl2_fx + 1.5 * pa2pb_zz_xz[j] * fl1_fx + pa2pb_z_xxx[j] * fl1_fx + 0.5 * pb_xxxz[j] * fl1_fx + pa2pb_zz_xxxz[j]);

                t_zz_xxyy[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_zz[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_zz_xx[j] * fl1_fx + 0.5 * pa2pb_zz_yy[j] * fl1_fx + 0.5 * pb_xxyy[j] * fl1_fx + pa2pb_zz_xxyy[j]);

                t_zz_xxyz[j] = fl_s_0_0 * (0.5 * pa2pb_z_y[j] * fl2_fx + 0.25 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_zz_yz[j] * fl1_fx + pa2pb_z_xxy[j] * fl1_fx + 0.5 * pb_xxyz[j] * fl1_fx + pa2pb_zz_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForDG_80_90(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (80,90)

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

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 70);

            auto pa2pb_z_xxz = pa2pbDistances.data(306 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(306 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(306 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(306 * idx + 82);

            auto pa2pb_z_yyy = pa2pbDistances.data(306 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(306 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(306 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(306 * idx + 86);

            auto pa2pb_zz_xx = pa2pbDistances.data(306 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(306 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(306 * idx + 277);

            auto pa2pb_zz_yy = pa2pbDistances.data(306 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(306 * idx + 279);

            auto pa2pb_zz_zz = pa2pbDistances.data(306 * idx + 280);

            auto pa2pb_zz_xxzz = pa2pbDistances.data(306 * idx + 296);

            auto pa2pb_zz_xyyy = pa2pbDistances.data(306 * idx + 297);

            auto pa2pb_zz_xyyz = pa2pbDistances.data(306 * idx + 298);

            auto pa2pb_zz_xyzz = pa2pbDistances.data(306 * idx + 299);

            auto pa2pb_zz_xzzz = pa2pbDistances.data(306 * idx + 300);

            auto pa2pb_zz_yyyy = pa2pbDistances.data(306 * idx + 301);

            auto pa2pb_zz_yyyz = pa2pbDistances.data(306 * idx + 302);

            auto pa2pb_zz_yyzz = pa2pbDistances.data(306 * idx + 303);

            auto pa2pb_zz_yzzz = pa2pbDistances.data(306 * idx + 304);

            auto pa2pb_zz_zzzz = pa2pbDistances.data(306 * idx + 305);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_zz_xxzz = primBuffer.data(90 * idx + 80);

            auto t_zz_xyyy = primBuffer.data(90 * idx + 81);

            auto t_zz_xyyz = primBuffer.data(90 * idx + 82);

            auto t_zz_xyzz = primBuffer.data(90 * idx + 83);

            auto t_zz_xzzz = primBuffer.data(90 * idx + 84);

            auto t_zz_yyyy = primBuffer.data(90 * idx + 85);

            auto t_zz_yyyz = primBuffer.data(90 * idx + 86);

            auto t_zz_yyzz = primBuffer.data(90 * idx + 87);

            auto t_zz_yzzz = primBuffer.data(90 * idx + 88);

            auto t_zz_zzzz = primBuffer.data(90 * idx + 89);

            // Batch of Integrals (80,90)

            #pragma omp simd aligned(fx, pa2pb_z_x, pa2pb_z_xxz, pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, \
                                     pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, pa2pb_z_z, pa2pb_z_zzz, \
                                     pa2pb_zz_xx, pa2pb_zz_xxzz, pa2pb_zz_xy, pa2pb_zz_xyyy, pa2pb_zz_xyyz, \
                                     pa2pb_zz_xyzz, pa2pb_zz_xz, pa2pb_zz_xzzz, pa2pb_zz_yy, pa2pb_zz_yyyy, \
                                     pa2pb_zz_yyyz, pa2pb_zz_yyzz, pa2pb_zz_yz, pa2pb_zz_yzzz, pa2pb_zz_zz, \
                                     pa2pb_zz_zzzz, pa_zz, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, pb_xzzz, pb_yy, \
                                     pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, s_0_0, t_zz_xxzz, \
                                     t_zz_xyyy, t_zz_xyyz, t_zz_xyzz, t_zz_xzzz, t_zz_yyyy, t_zz_yyyz, t_zz_yyzz, \
                                     t_zz_yzzz, t_zz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_zz_xxzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_zz[j] * fl2_fx + pa2pb_z_z[j] * fl2_fx + 0.75 * pb_xx[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_zz_xx[j] * fl1_fx + 0.5 * pa2pb_zz_zz[j] * fl1_fx + 2.0 * pa2pb_z_xxz[j] * fl1_fx + 0.5 * pb_xxzz[j] * fl1_fx + pa2pb_zz_xxzz[j]);

                t_zz_xyyy[j] = fl_s_0_0 * (0.75 * pb_xy[j] * fl2_fx + 1.5 * pa2pb_zz_xy[j] * fl1_fx + 0.5 * pb_xyyy[j] * fl1_fx + pa2pb_zz_xyyy[j]);

                t_zz_xyyz[j] = fl_s_0_0 * (0.5 * pa2pb_z_x[j] * fl2_fx + 0.25 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_zz_xz[j] * fl1_fx + pa2pb_z_xyy[j] * fl1_fx + 0.5 * pb_xyyz[j] * fl1_fx + pa2pb_zz_xyyz[j]);

                t_zz_xyzz[j] = fl_s_0_0 * (0.75 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_zz_xy[j] * fl1_fx + 2.0 * pa2pb_z_xyz[j] * fl1_fx + 0.5 * pb_xyzz[j] * fl1_fx + pa2pb_zz_xyzz[j]);

                t_zz_xzzz[j] = fl_s_0_0 * (1.5 * pa2pb_z_x[j] * fl2_fx + 2.25 * pb_xz[j] * fl2_fx + 1.5 * pa2pb_zz_xz[j] * fl1_fx + 3.0 * pa2pb_z_xzz[j] * fl1_fx + 0.5 * pb_xzzz[j] * fl1_fx + pa2pb_zz_xzzz[j]);

                t_zz_yyyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 1.5 * pb_yy[j] * fl2_fx + 3.0 * pa2pb_zz_yy[j] * fl1_fx + 0.5 * pb_yyyy[j] * fl1_fx + pa2pb_zz_yyyy[j]);

                t_zz_yyyz[j] = fl_s_0_0 * (1.5 * pa2pb_z_y[j] * fl2_fx + 0.75 * pb_yz[j] * fl2_fx + 1.5 * pa2pb_zz_yz[j] * fl1_fx + pa2pb_z_yyy[j] * fl1_fx + 0.5 * pb_yyyz[j] * fl1_fx + pa2pb_zz_yyyz[j]);

                t_zz_yyzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_zz[j] * fl2_fx + pa2pb_z_z[j] * fl2_fx + 0.75 * pb_yy[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_zz_yy[j] * fl1_fx + 0.5 * pa2pb_zz_zz[j] * fl1_fx + 2.0 * pa2pb_z_yyz[j] * fl1_fx + 0.5 * pb_yyzz[j] * fl1_fx + pa2pb_zz_yyzz[j]);

                t_zz_yzzz[j] = fl_s_0_0 * (1.5 * pa2pb_z_y[j] * fl2_fx + 2.25 * pb_yz[j] * fl2_fx + 1.5 * pa2pb_zz_yz[j] * fl1_fx + 3.0 * pa2pb_z_yzz[j] * fl1_fx + 0.5 * pb_yzzz[j] * fl1_fx + pa2pb_zz_yzzz[j]);

                t_zz_zzzz[j] = fl_s_0_0 * (1.875 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 6.0 * pa2pb_z_z[j] * fl2_fx + 4.5 * pb_zz[j] * fl2_fx + 3.0 * pa2pb_zz_zz[j] * fl1_fx + 4.0 * pa2pb_z_zzz[j] * fl1_fx + 0.5 * pb_zzzz[j] * fl1_fx + pa2pb_zz_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
                     const CMemBlock2D<double>& pbDistances,
                     const CMemBlock2D<double>& pa2pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        ovlrecfunc::compOverlapForGD_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGD_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGD_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGD_30_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGD_40_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGD_50_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGD_60_70(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGD_70_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGD_80_90(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances,
                                           braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForGD_0_10(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pbDistances,
                          const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_xx_xx = pa2pbDistances.data(306 * idx + 30);

            auto pa2pb_xx_xy = pa2pbDistances.data(306 * idx + 31);

            auto pa2pb_xx_xz = pa2pbDistances.data(306 * idx + 32);

            auto pa2pb_xx_yy = pa2pbDistances.data(306 * idx + 33);

            auto pa2pb_xx_yz = pa2pbDistances.data(306 * idx + 34);

            auto pa2pb_xx_zz = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_xy_xx = pa2pbDistances.data(306 * idx + 39);

            auto pa2pb_xy_xy = pa2pbDistances.data(306 * idx + 40);

            auto pa2pb_xy_xz = pa2pbDistances.data(306 * idx + 41);

            auto pa2pb_xy_yy = pa2pbDistances.data(306 * idx + 42);

            auto pa2pb_xxx_x = pa2pbDistances.data(306 * idx + 81);

            auto pa2pb_xxx_y = pa2pbDistances.data(306 * idx + 82);

            auto pa2pb_xxx_z = pa2pbDistances.data(306 * idx + 83);

            auto pa2pb_xxy_x = pa2pbDistances.data(306 * idx + 90);

            auto pa2pb_xxy_y = pa2pbDistances.data(306 * idx + 91);

            auto pa2pb_xxy_z = pa2pbDistances.data(306 * idx + 92);

            auto pa2pb_xxxx_xx = pa2pbDistances.data(306 * idx + 174);

            auto pa2pb_xxxx_xy = pa2pbDistances.data(306 * idx + 175);

            auto pa2pb_xxxx_xz = pa2pbDistances.data(306 * idx + 176);

            auto pa2pb_xxxx_yy = pa2pbDistances.data(306 * idx + 177);

            auto pa2pb_xxxx_yz = pa2pbDistances.data(306 * idx + 178);

            auto pa2pb_xxxx_zz = pa2pbDistances.data(306 * idx + 179);

            auto pa2pb_xxxy_xx = pa2pbDistances.data(306 * idx + 183);

            auto pa2pb_xxxy_xy = pa2pbDistances.data(306 * idx + 184);

            auto pa2pb_xxxy_xz = pa2pbDistances.data(306 * idx + 185);

            auto pa2pb_xxxy_yy = pa2pbDistances.data(306 * idx + 186);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxxx_xx = primBuffer.data(90 * idx);

            auto t_xxxx_xy = primBuffer.data(90 * idx + 1);

            auto t_xxxx_xz = primBuffer.data(90 * idx + 2);

            auto t_xxxx_yy = primBuffer.data(90 * idx + 3);

            auto t_xxxx_yz = primBuffer.data(90 * idx + 4);

            auto t_xxxx_zz = primBuffer.data(90 * idx + 5);

            auto t_xxxy_xx = primBuffer.data(90 * idx + 6);

            auto t_xxxy_xy = primBuffer.data(90 * idx + 7);

            auto t_xxxy_xz = primBuffer.data(90 * idx + 8);

            auto t_xxxy_yy = primBuffer.data(90 * idx + 9);

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xy, \
                                     pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa2pb_xxx_x, pa2pb_xxx_y, \
                                     pa2pb_xxx_z, pa2pb_xxxx_xx, pa2pb_xxxx_xy, pa2pb_xxxx_xz, pa2pb_xxxx_yy, \
                                     pa2pb_xxxx_yz, pa2pb_xxxx_zz, pa2pb_xxxy_xx, pa2pb_xxxy_xy, pa2pb_xxxy_xz, \
                                     pa2pb_xxxy_yy, pa2pb_xxy_x, pa2pb_xxy_y, pa2pb_xxy_z, pa2pb_xy_xx, pa2pb_xy_xy, \
                                     pa2pb_xy_xz, pa2pb_xy_yy, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa_xx, pa_xxxx, pa_xxxy, \
                                     pa_xy, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, t_xxxx_xx, t_xxxx_xy, \
                                     t_xxxx_xz, t_xxxx_yy, t_xxxx_yz, t_xxxx_zz, t_xxxy_xx, t_xxxy_xy, t_xxxy_xz, \
                                     t_xxxy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxxx_xx[j] = fl_s_0_0 * (1.875 * fl3_fx + 4.5 * pa_xx[j] * fl2_fx + 6.0 * pa2pb_x_x[j] * fl2_fx + 0.5 * pa_xxxx[j] * fl1_fx + 4.0 * pa2pb_xxx_x[j] * fl1_fx + 0.75 * pb_xx[j] * fl2_fx + 3.0 * pa2pb_xx_xx[j] * fl1_fx + pa2pb_xxxx_xx[j]);

                t_xxxx_xy[j] = fl_s_0_0 * (3.0 * pa2pb_x_y[j] * fl2_fx + 2.0 * pa2pb_xxx_y[j] * fl1_fx + 0.75 * pb_xy[j] * fl2_fx + 3.0 * pa2pb_xx_xy[j] * fl1_fx + pa2pb_xxxx_xy[j]);

                t_xxxx_xz[j] = fl_s_0_0 * (3.0 * pa2pb_x_z[j] * fl2_fx + 2.0 * pa2pb_xxx_z[j] * fl1_fx + 0.75 * pb_xz[j] * fl2_fx + 3.0 * pa2pb_xx_xz[j] * fl1_fx + pa2pb_xxxx_xz[j]);

                t_xxxx_yy[j] = fl_s_0_0 * (0.375 * fl3_fx + 1.5 * pa_xx[j] * fl2_fx + 0.5 * pa_xxxx[j] * fl1_fx + 0.75 * pb_yy[j] * fl2_fx + 3.0 * pa2pb_xx_yy[j] * fl1_fx + pa2pb_xxxx_yy[j]);

                t_xxxx_yz[j] = fl_s_0_0 * (0.75 * pb_yz[j] * fl2_fx + 3.0 * pa2pb_xx_yz[j] * fl1_fx + pa2pb_xxxx_yz[j]);

                t_xxxx_zz[j] = fl_s_0_0 * (0.375 * fl3_fx + 1.5 * pa_xx[j] * fl2_fx + 0.5 * pa_xxxx[j] * fl1_fx + 0.75 * pb_zz[j] * fl2_fx + 3.0 * pa2pb_xx_zz[j] * fl1_fx + pa2pb_xxxx_zz[j]);

                t_xxxy_xx[j] = fl_s_0_0 * (2.25 * pa_xy[j] * fl2_fx + 1.5 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa_xxxy[j] * fl1_fx + 3.0 * pa2pb_xxy_x[j] * fl1_fx + 1.5 * pa2pb_xy_xx[j] * fl1_fx + pa2pb_xxxy_xx[j]);

                t_xxxy_xy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.5 * pa2pb_xxx_x[j] * fl1_fx + 1.5 * pa2pb_xxy_y[j] * fl1_fx + 1.5 * pa2pb_xy_xy[j] * fl1_fx + pa2pb_xxxy_xy[j]);

                t_xxxy_xz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl2_fx + 1.5 * pa2pb_xxy_z[j] * fl1_fx + 1.5 * pa2pb_xy_xz[j] * fl1_fx + pa2pb_xxxy_xz[j]);

                t_xxxy_yy[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 1.5 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa_xxxy[j] * fl1_fx + pa2pb_xxx_y[j] * fl1_fx + 1.5 * pa2pb_xy_yy[j] * fl1_fx + pa2pb_xxxy_yy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD_10_20(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (10,20)

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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 20);

            auto pa2pb_xx_xx = pa2pbDistances.data(306 * idx + 30);

            auto pa2pb_xx_xy = pa2pbDistances.data(306 * idx + 31);

            auto pa2pb_xy_yz = pa2pbDistances.data(306 * idx + 43);

            auto pa2pb_xy_zz = pa2pbDistances.data(306 * idx + 44);

            auto pa2pb_xz_xx = pa2pbDistances.data(306 * idx + 48);

            auto pa2pb_xz_xy = pa2pbDistances.data(306 * idx + 49);

            auto pa2pb_xz_xz = pa2pbDistances.data(306 * idx + 50);

            auto pa2pb_xz_yy = pa2pbDistances.data(306 * idx + 51);

            auto pa2pb_xz_yz = pa2pbDistances.data(306 * idx + 52);

            auto pa2pb_xz_zz = pa2pbDistances.data(306 * idx + 53);

            auto pa2pb_yy_xx = pa2pbDistances.data(306 * idx + 57);

            auto pa2pb_yy_xy = pa2pbDistances.data(306 * idx + 58);

            auto pa2pb_xxx_x = pa2pbDistances.data(306 * idx + 81);

            auto pa2pb_xxx_y = pa2pbDistances.data(306 * idx + 82);

            auto pa2pb_xxx_z = pa2pbDistances.data(306 * idx + 83);

            auto pa2pb_xxy_x = pa2pbDistances.data(306 * idx + 90);

            auto pa2pb_xxz_x = pa2pbDistances.data(306 * idx + 99);

            auto pa2pb_xxz_y = pa2pbDistances.data(306 * idx + 100);

            auto pa2pb_xxz_z = pa2pbDistances.data(306 * idx + 101);

            auto pa2pb_xyy_x = pa2pbDistances.data(306 * idx + 108);

            auto pa2pb_xyy_y = pa2pbDistances.data(306 * idx + 109);

            auto pa2pb_xxxy_yz = pa2pbDistances.data(306 * idx + 187);

            auto pa2pb_xxxy_zz = pa2pbDistances.data(306 * idx + 188);

            auto pa2pb_xxxz_xx = pa2pbDistances.data(306 * idx + 192);

            auto pa2pb_xxxz_xy = pa2pbDistances.data(306 * idx + 193);

            auto pa2pb_xxxz_xz = pa2pbDistances.data(306 * idx + 194);

            auto pa2pb_xxxz_yy = pa2pbDistances.data(306 * idx + 195);

            auto pa2pb_xxxz_yz = pa2pbDistances.data(306 * idx + 196);

            auto pa2pb_xxxz_zz = pa2pbDistances.data(306 * idx + 197);

            auto pa2pb_xxyy_xx = pa2pbDistances.data(306 * idx + 201);

            auto pa2pb_xxyy_xy = pa2pbDistances.data(306 * idx + 202);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxxy_yz = primBuffer.data(90 * idx + 10);

            auto t_xxxy_zz = primBuffer.data(90 * idx + 11);

            auto t_xxxz_xx = primBuffer.data(90 * idx + 12);

            auto t_xxxz_xy = primBuffer.data(90 * idx + 13);

            auto t_xxxz_xz = primBuffer.data(90 * idx + 14);

            auto t_xxxz_yy = primBuffer.data(90 * idx + 15);

            auto t_xxxz_yz = primBuffer.data(90 * idx + 16);

            auto t_xxxz_zz = primBuffer.data(90 * idx + 17);

            auto t_xxyy_xx = primBuffer.data(90 * idx + 18);

            auto t_xxyy_xy = primBuffer.data(90 * idx + 19);

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xy, \
                                     pa2pb_xxx_x, pa2pb_xxx_y, pa2pb_xxx_z, pa2pb_xxxy_yz, pa2pb_xxxy_zz, \
                                     pa2pb_xxxz_xx, pa2pb_xxxz_xy, pa2pb_xxxz_xz, pa2pb_xxxz_yy, pa2pb_xxxz_yz, \
                                     pa2pb_xxxz_zz, pa2pb_xxy_x, pa2pb_xxyy_xx, pa2pb_xxyy_xy, pa2pb_xxz_x, pa2pb_xxz_y, \
                                     pa2pb_xxz_z, pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_xyy_x, pa2pb_xyy_y, pa2pb_xz_xx, \
                                     pa2pb_xz_xy, pa2pb_xz_xz, pa2pb_xz_yy, pa2pb_xz_yz, pa2pb_xz_zz, pa2pb_y_x, \
                                     pa2pb_yy_xx, pa2pb_yy_xy, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, pa_xx, pa_xxxy, pa_xxxz, \
                                     pa_xxyy, pa_xy, pa_xz, pa_yy, pb_xx, pb_xy, s_0_0, t_xxxy_yz, t_xxxy_zz, t_xxxz_xx, \
                                     t_xxxz_xy, t_xxxz_xz, t_xxxz_yy, t_xxxz_yz, t_xxxz_zz, t_xxyy_xx, t_xxyy_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxxy_yz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa2pb_xxx_z[j] * fl1_fx + 1.5 * pa2pb_xy_yz[j] * fl1_fx + pa2pb_xxxy_yz[j]);

                t_xxxy_zz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 0.5 * pa_xxxy[j] * fl1_fx + 1.5 * pa2pb_xy_zz[j] * fl1_fx + pa2pb_xxxy_zz[j]);

                t_xxxz_xx[j] = fl_s_0_0 * (2.25 * pa_xz[j] * fl2_fx + 1.5 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa_xxxz[j] * fl1_fx + 3.0 * pa2pb_xxz_x[j] * fl1_fx + 1.5 * pa2pb_xz_xx[j] * fl1_fx + pa2pb_xxxz_xx[j]);

                t_xxxz_xy[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl2_fx + 1.5 * pa2pb_xxz_y[j] * fl1_fx + 1.5 * pa2pb_xz_xy[j] * fl1_fx + pa2pb_xxxz_xy[j]);

                t_xxxz_xz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.5 * pa2pb_xxx_x[j] * fl1_fx + 1.5 * pa2pb_xxz_z[j] * fl1_fx + 1.5 * pa2pb_xz_xz[j] * fl1_fx + pa2pb_xxxz_xz[j]);

                t_xxxz_yy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 0.5 * pa_xxxz[j] * fl1_fx + 1.5 * pa2pb_xz_yy[j] * fl1_fx + pa2pb_xxxz_yy[j]);

                t_xxxz_yz[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_xxx_y[j] * fl1_fx + 1.5 * pa2pb_xz_yz[j] * fl1_fx + pa2pb_xxxz_yz[j]);

                t_xxxz_zz[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 1.5 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa_xxxz[j] * fl1_fx + pa2pb_xxx_z[j] * fl1_fx + 1.5 * pa2pb_xz_zz[j] * fl1_fx + pa2pb_xxxz_zz[j]);

                t_xxyy_xx[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.25 * pa_xx[j] * fl2_fx + pa2pb_x_x[j] * fl2_fx + 0.5 * pa_xxyy[j] * fl1_fx + 2.0 * pa2pb_xyy_x[j] * fl1_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_xx_xx[j] * fl1_fx + 0.5 * pa2pb_yy_xx[j] * fl1_fx + pa2pb_xxyy_xx[j]);

                t_xxyy_xy[j] = fl_s_0_0 * (pa_xy[j] * fl2_fx + 0.5 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_y_x[j] * fl2_fx + pa2pb_xxy_x[j] * fl1_fx + pa2pb_xyy_y[j] * fl1_fx + 0.25 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xx_xy[j] * fl1_fx + 0.5 * pa2pb_yy_xy[j] * fl1_fx + pa2pb_xxyy_xy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD_20_30(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (20,30)

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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 20);

            auto pa2pb_xx_xz = pa2pbDistances.data(306 * idx + 32);

            auto pa2pb_xx_yy = pa2pbDistances.data(306 * idx + 33);

            auto pa2pb_xx_yz = pa2pbDistances.data(306 * idx + 34);

            auto pa2pb_xx_zz = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_yy_xz = pa2pbDistances.data(306 * idx + 59);

            auto pa2pb_yy_yy = pa2pbDistances.data(306 * idx + 60);

            auto pa2pb_yy_yz = pa2pbDistances.data(306 * idx + 61);

            auto pa2pb_yy_zz = pa2pbDistances.data(306 * idx + 62);

            auto pa2pb_yz_xx = pa2pbDistances.data(306 * idx + 66);

            auto pa2pb_yz_xy = pa2pbDistances.data(306 * idx + 67);

            auto pa2pb_yz_xz = pa2pbDistances.data(306 * idx + 68);

            auto pa2pb_yz_yy = pa2pbDistances.data(306 * idx + 69);

            auto pa2pb_yz_yz = pa2pbDistances.data(306 * idx + 70);

            auto pa2pb_yz_zz = pa2pbDistances.data(306 * idx + 71);

            auto pa2pb_xxy_x = pa2pbDistances.data(306 * idx + 90);

            auto pa2pb_xxy_y = pa2pbDistances.data(306 * idx + 91);

            auto pa2pb_xxy_z = pa2pbDistances.data(306 * idx + 92);

            auto pa2pb_xxz_x = pa2pbDistances.data(306 * idx + 99);

            auto pa2pb_xxz_y = pa2pbDistances.data(306 * idx + 100);

            auto pa2pb_xxz_z = pa2pbDistances.data(306 * idx + 101);

            auto pa2pb_xyy_z = pa2pbDistances.data(306 * idx + 110);

            auto pa2pb_xyz_x = pa2pbDistances.data(306 * idx + 117);

            auto pa2pb_xyz_y = pa2pbDistances.data(306 * idx + 118);

            auto pa2pb_xyz_z = pa2pbDistances.data(306 * idx + 119);

            auto pa2pb_xxyy_xz = pa2pbDistances.data(306 * idx + 203);

            auto pa2pb_xxyy_yy = pa2pbDistances.data(306 * idx + 204);

            auto pa2pb_xxyy_yz = pa2pbDistances.data(306 * idx + 205);

            auto pa2pb_xxyy_zz = pa2pbDistances.data(306 * idx + 206);

            auto pa2pb_xxyz_xx = pa2pbDistances.data(306 * idx + 210);

            auto pa2pb_xxyz_xy = pa2pbDistances.data(306 * idx + 211);

            auto pa2pb_xxyz_xz = pa2pbDistances.data(306 * idx + 212);

            auto pa2pb_xxyz_yy = pa2pbDistances.data(306 * idx + 213);

            auto pa2pb_xxyz_yz = pa2pbDistances.data(306 * idx + 214);

            auto pa2pb_xxyz_zz = pa2pbDistances.data(306 * idx + 215);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxyy_xz = primBuffer.data(90 * idx + 20);

            auto t_xxyy_yy = primBuffer.data(90 * idx + 21);

            auto t_xxyy_yz = primBuffer.data(90 * idx + 22);

            auto t_xxyy_zz = primBuffer.data(90 * idx + 23);

            auto t_xxyz_xx = primBuffer.data(90 * idx + 24);

            auto t_xxyz_xy = primBuffer.data(90 * idx + 25);

            auto t_xxyz_xz = primBuffer.data(90 * idx + 26);

            auto t_xxyz_yy = primBuffer.data(90 * idx + 27);

            auto t_xxyz_yz = primBuffer.data(90 * idx + 28);

            auto t_xxyz_zz = primBuffer.data(90 * idx + 29);

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fx, pa2pb_x_z, pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, \
                                     pa2pb_xxy_x, pa2pb_xxy_y, pa2pb_xxy_z, pa2pb_xxyy_xz, pa2pb_xxyy_yy, \
                                     pa2pb_xxyy_yz, pa2pb_xxyy_zz, pa2pb_xxyz_xx, pa2pb_xxyz_xy, pa2pb_xxyz_xz, \
                                     pa2pb_xxyz_yy, pa2pb_xxyz_yz, pa2pb_xxyz_zz, pa2pb_xxz_x, pa2pb_xxz_y, pa2pb_xxz_z, \
                                     pa2pb_xyy_z, pa2pb_xyz_x, pa2pb_xyz_y, pa2pb_xyz_z, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, \
                                     pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, pa2pb_yy_zz, pa2pb_yz_xx, pa2pb_yz_xy, \
                                     pa2pb_yz_xz, pa2pb_yz_yy, pa2pb_yz_yz, pa2pb_yz_zz, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, \
                                     pa_xx, pa_xxyy, pa_xxyz, pa_xy, pa_xz, pa_yy, pa_yz, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, \
                                     t_xxyy_xz, t_xxyy_yy, t_xxyy_yz, t_xxyy_zz, t_xxyz_xx, t_xxyz_xy, t_xxyz_xz, \
                                     t_xxyz_yy, t_xxyz_yz, t_xxyz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxyy_xz[j] = fl_s_0_0 * (0.5 * pa2pb_x_z[j] * fl2_fx + pa2pb_xyy_z[j] * fl1_fx + 0.25 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xx_xz[j] * fl1_fx + 0.5 * pa2pb_yy_xz[j] * fl1_fx + pa2pb_xxyy_xz[j]);

                t_xxyy_yy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.25 * pa_yy[j] * fl2_fx + pa2pb_y_y[j] * fl2_fx + 0.5 * pa_xxyy[j] * fl1_fx + 2.0 * pa2pb_xxy_y[j] * fl1_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_xx_yy[j] * fl1_fx + 0.5 * pa2pb_yy_yy[j] * fl1_fx + pa2pb_xxyy_yy[j]);

                t_xxyy_yz[j] = fl_s_0_0 * (0.5 * pa2pb_y_z[j] * fl2_fx + pa2pb_xxy_z[j] * fl1_fx + 0.25 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_xx_yz[j] * fl1_fx + 0.5 * pa2pb_yy_yz[j] * fl1_fx + pa2pb_xxyy_yz[j]);

                t_xxyy_zz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + 0.25 * pa_yy[j] * fl2_fx + 0.5 * pa_xxyy[j] * fl1_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_xx_zz[j] * fl1_fx + 0.5 * pa2pb_yy_zz[j] * fl1_fx + pa2pb_xxyy_zz[j]);

                t_xxyz_xx[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 0.5 * pa_xxyz[j] * fl1_fx + 2.0 * pa2pb_xyz_x[j] * fl1_fx + 0.5 * pa2pb_yz_xx[j] * fl1_fx + pa2pb_xxyz_xx[j]);

                t_xxyz_xy[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + 0.25 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa2pb_xxz_x[j] * fl1_fx + pa2pb_xyz_y[j] * fl1_fx + 0.5 * pa2pb_yz_xy[j] * fl1_fx + pa2pb_xxyz_xy[j]);

                t_xxyz_xz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + 0.25 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa2pb_xxy_x[j] * fl1_fx + pa2pb_xyz_z[j] * fl1_fx + 0.5 * pa2pb_yz_xz[j] * fl1_fx + pa2pb_xxyz_xz[j]);

                t_xxyz_yy[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl2_fx + 0.5 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa_xxyz[j] * fl1_fx + pa2pb_xxz_y[j] * fl1_fx + 0.5 * pa2pb_yz_yy[j] * fl1_fx + pa2pb_xxyz_yy[j]);

                t_xxyz_yz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + 0.25 * pa2pb_y_y[j] * fl2_fx + 0.25 * pa2pb_z_z[j] * fl2_fx + 0.5 * pa2pb_xxy_y[j] * fl1_fx + 0.5 * pa2pb_xxz_z[j] * fl1_fx + 0.5 * pa2pb_yz_yz[j] * fl1_fx + pa2pb_xxyz_yz[j]);

                t_xxyz_zz[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl2_fx + 0.5 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa_xxyz[j] * fl1_fx + pa2pb_xxy_z[j] * fl1_fx + 0.5 * pa2pb_yz_zz[j] * fl1_fx + pa2pb_xxyz_zz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD_30_40(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (30,40)

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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 20);

            auto pa2pb_xx_xx = pa2pbDistances.data(306 * idx + 30);

            auto pa2pb_xx_xy = pa2pbDistances.data(306 * idx + 31);

            auto pa2pb_xx_xz = pa2pbDistances.data(306 * idx + 32);

            auto pa2pb_xx_yy = pa2pbDistances.data(306 * idx + 33);

            auto pa2pb_xx_yz = pa2pbDistances.data(306 * idx + 34);

            auto pa2pb_xx_zz = pa2pbDistances.data(306 * idx + 35);

            auto pa2pb_xy_xx = pa2pbDistances.data(306 * idx + 39);

            auto pa2pb_xy_xy = pa2pbDistances.data(306 * idx + 40);

            auto pa2pb_xy_xz = pa2pbDistances.data(306 * idx + 41);

            auto pa2pb_xy_yy = pa2pbDistances.data(306 * idx + 42);

            auto pa2pb_zz_xx = pa2pbDistances.data(306 * idx + 75);

            auto pa2pb_zz_xy = pa2pbDistances.data(306 * idx + 76);

            auto pa2pb_zz_xz = pa2pbDistances.data(306 * idx + 77);

            auto pa2pb_zz_yy = pa2pbDistances.data(306 * idx + 78);

            auto pa2pb_zz_yz = pa2pbDistances.data(306 * idx + 79);

            auto pa2pb_zz_zz = pa2pbDistances.data(306 * idx + 80);

            auto pa2pb_xxz_x = pa2pbDistances.data(306 * idx + 99);

            auto pa2pb_xxz_y = pa2pbDistances.data(306 * idx + 100);

            auto pa2pb_xxz_z = pa2pbDistances.data(306 * idx + 101);

            auto pa2pb_xyy_x = pa2pbDistances.data(306 * idx + 108);

            auto pa2pb_xyy_y = pa2pbDistances.data(306 * idx + 109);

            auto pa2pb_xzz_x = pa2pbDistances.data(306 * idx + 126);

            auto pa2pb_xzz_y = pa2pbDistances.data(306 * idx + 127);

            auto pa2pb_xzz_z = pa2pbDistances.data(306 * idx + 128);

            auto pa2pb_yyy_x = pa2pbDistances.data(306 * idx + 135);

            auto pa2pb_yyy_y = pa2pbDistances.data(306 * idx + 136);

            auto pa2pb_yyy_z = pa2pbDistances.data(306 * idx + 137);

            auto pa2pb_xxzz_xx = pa2pbDistances.data(306 * idx + 219);

            auto pa2pb_xxzz_xy = pa2pbDistances.data(306 * idx + 220);

            auto pa2pb_xxzz_xz = pa2pbDistances.data(306 * idx + 221);

            auto pa2pb_xxzz_yy = pa2pbDistances.data(306 * idx + 222);

            auto pa2pb_xxzz_yz = pa2pbDistances.data(306 * idx + 223);

            auto pa2pb_xxzz_zz = pa2pbDistances.data(306 * idx + 224);

            auto pa2pb_xyyy_xx = pa2pbDistances.data(306 * idx + 228);

            auto pa2pb_xyyy_xy = pa2pbDistances.data(306 * idx + 229);

            auto pa2pb_xyyy_xz = pa2pbDistances.data(306 * idx + 230);

            auto pa2pb_xyyy_yy = pa2pbDistances.data(306 * idx + 231);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxzz_xx = primBuffer.data(90 * idx + 30);

            auto t_xxzz_xy = primBuffer.data(90 * idx + 31);

            auto t_xxzz_xz = primBuffer.data(90 * idx + 32);

            auto t_xxzz_yy = primBuffer.data(90 * idx + 33);

            auto t_xxzz_yz = primBuffer.data(90 * idx + 34);

            auto t_xxzz_zz = primBuffer.data(90 * idx + 35);

            auto t_xyyy_xx = primBuffer.data(90 * idx + 36);

            auto t_xyyy_xy = primBuffer.data(90 * idx + 37);

            auto t_xyyy_xz = primBuffer.data(90 * idx + 38);

            auto t_xyyy_yy = primBuffer.data(90 * idx + 39);

            // Batch of Integrals (30,40)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xy, \
                                     pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa2pb_xxz_x, pa2pb_xxz_y, \
                                     pa2pb_xxz_z, pa2pb_xxzz_xx, pa2pb_xxzz_xy, pa2pb_xxzz_xz, pa2pb_xxzz_yy, \
                                     pa2pb_xxzz_yz, pa2pb_xxzz_zz, pa2pb_xy_xx, pa2pb_xy_xy, pa2pb_xy_xz, pa2pb_xy_yy, \
                                     pa2pb_xyy_x, pa2pb_xyy_y, pa2pb_xyyy_xx, pa2pb_xyyy_xy, pa2pb_xyyy_xz, \
                                     pa2pb_xyyy_yy, pa2pb_xzz_x, pa2pb_xzz_y, pa2pb_xzz_z, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, \
                                     pa2pb_yyy_x, pa2pb_yyy_y, pa2pb_yyy_z, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, pa2pb_zz_xx, \
                                     pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zz_zz, pa_xx, pa_xxzz, \
                                     pa_xy, pa_xyyy, pa_xz, pa_yy, pa_zz, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, \
                                     t_xxzz_xx, t_xxzz_xy, t_xxzz_xz, t_xxzz_yy, t_xxzz_yz, t_xxzz_zz, t_xyyy_xx, \
                                     t_xyyy_xy, t_xyyy_xz, t_xyyy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxzz_xx[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.25 * pa_xx[j] * fl2_fx + pa2pb_x_x[j] * fl2_fx + 0.5 * pa_xxzz[j] * fl1_fx + 2.0 * pa2pb_xzz_x[j] * fl1_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_xx_xx[j] * fl1_fx + 0.5 * pa2pb_zz_xx[j] * fl1_fx + pa2pb_xxzz_xx[j]);

                t_xxzz_xy[j] = fl_s_0_0 * (0.5 * pa2pb_x_y[j] * fl2_fx + pa2pb_xzz_y[j] * fl1_fx + 0.25 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xx_xy[j] * fl1_fx + 0.5 * pa2pb_zz_xy[j] * fl1_fx + pa2pb_xxzz_xy[j]);

                t_xxzz_xz[j] = fl_s_0_0 * (pa_xz[j] * fl2_fx + 0.5 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa2pb_z_x[j] * fl2_fx + pa2pb_xxz_x[j] * fl1_fx + pa2pb_xzz_z[j] * fl1_fx + 0.25 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xx_xz[j] * fl1_fx + 0.5 * pa2pb_zz_xz[j] * fl1_fx + pa2pb_xxzz_xz[j]);

                t_xxzz_yy[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + 0.25 * pa_zz[j] * fl2_fx + 0.5 * pa_xxzz[j] * fl1_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_xx_yy[j] * fl1_fx + 0.5 * pa2pb_zz_yy[j] * fl1_fx + pa2pb_xxzz_yy[j]);

                t_xxzz_yz[j] = fl_s_0_0 * (0.5 * pa2pb_z_y[j] * fl2_fx + pa2pb_xxz_y[j] * fl1_fx + 0.25 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_xx_yz[j] * fl1_fx + 0.5 * pa2pb_zz_yz[j] * fl1_fx + pa2pb_xxzz_yz[j]);

                t_xxzz_zz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.25 * pa_zz[j] * fl2_fx + pa2pb_z_z[j] * fl2_fx + 0.5 * pa_xxzz[j] * fl1_fx + 2.0 * pa2pb_xxz_z[j] * fl1_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_xx_zz[j] * fl1_fx + 0.5 * pa2pb_zz_zz[j] * fl1_fx + pa2pb_xxzz_zz[j]);

                t_xyyy_xx[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 1.5 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa_xyyy[j] * fl1_fx + pa2pb_yyy_x[j] * fl1_fx + 1.5 * pa2pb_xy_xx[j] * fl1_fx + pa2pb_xyyy_xx[j]);

                t_xyyy_xy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 1.5 * pa2pb_xyy_x[j] * fl1_fx + 0.5 * pa2pb_yyy_y[j] * fl1_fx + 1.5 * pa2pb_xy_xy[j] * fl1_fx + pa2pb_xyyy_xy[j]);

                t_xyyy_xz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa2pb_yyy_z[j] * fl1_fx + 1.5 * pa2pb_xy_xz[j] * fl1_fx + pa2pb_xyyy_xz[j]);

                t_xyyy_yy[j] = fl_s_0_0 * (2.25 * pa_xy[j] * fl2_fx + 1.5 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa_xyyy[j] * fl1_fx + 3.0 * pa2pb_xyy_y[j] * fl1_fx + 1.5 * pa2pb_xy_yy[j] * fl1_fx + pa2pb_xyyy_yy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD_40_50(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (40,50)

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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 20);

            auto pa2pb_xy_xx = pa2pbDistances.data(306 * idx + 39);

            auto pa2pb_xy_xy = pa2pbDistances.data(306 * idx + 40);

            auto pa2pb_xy_yz = pa2pbDistances.data(306 * idx + 43);

            auto pa2pb_xy_zz = pa2pbDistances.data(306 * idx + 44);

            auto pa2pb_xz_xx = pa2pbDistances.data(306 * idx + 48);

            auto pa2pb_xz_xy = pa2pbDistances.data(306 * idx + 49);

            auto pa2pb_xz_xz = pa2pbDistances.data(306 * idx + 50);

            auto pa2pb_xz_yy = pa2pbDistances.data(306 * idx + 51);

            auto pa2pb_xz_yz = pa2pbDistances.data(306 * idx + 52);

            auto pa2pb_xz_zz = pa2pbDistances.data(306 * idx + 53);

            auto pa2pb_xyy_x = pa2pbDistances.data(306 * idx + 108);

            auto pa2pb_xyy_y = pa2pbDistances.data(306 * idx + 109);

            auto pa2pb_xyy_z = pa2pbDistances.data(306 * idx + 110);

            auto pa2pb_xyz_x = pa2pbDistances.data(306 * idx + 117);

            auto pa2pb_xyz_y = pa2pbDistances.data(306 * idx + 118);

            auto pa2pb_xyz_z = pa2pbDistances.data(306 * idx + 119);

            auto pa2pb_xzz_x = pa2pbDistances.data(306 * idx + 126);

            auto pa2pb_yyz_x = pa2pbDistances.data(306 * idx + 144);

            auto pa2pb_yyz_y = pa2pbDistances.data(306 * idx + 145);

            auto pa2pb_yyz_z = pa2pbDistances.data(306 * idx + 146);

            auto pa2pb_yzz_x = pa2pbDistances.data(306 * idx + 153);

            auto pa2pb_yzz_y = pa2pbDistances.data(306 * idx + 154);

            auto pa2pb_xyyy_yz = pa2pbDistances.data(306 * idx + 232);

            auto pa2pb_xyyy_zz = pa2pbDistances.data(306 * idx + 233);

            auto pa2pb_xyyz_xx = pa2pbDistances.data(306 * idx + 237);

            auto pa2pb_xyyz_xy = pa2pbDistances.data(306 * idx + 238);

            auto pa2pb_xyyz_xz = pa2pbDistances.data(306 * idx + 239);

            auto pa2pb_xyyz_yy = pa2pbDistances.data(306 * idx + 240);

            auto pa2pb_xyyz_yz = pa2pbDistances.data(306 * idx + 241);

            auto pa2pb_xyyz_zz = pa2pbDistances.data(306 * idx + 242);

            auto pa2pb_xyzz_xx = pa2pbDistances.data(306 * idx + 246);

            auto pa2pb_xyzz_xy = pa2pbDistances.data(306 * idx + 247);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xyyy_yz = primBuffer.data(90 * idx + 40);

            auto t_xyyy_zz = primBuffer.data(90 * idx + 41);

            auto t_xyyz_xx = primBuffer.data(90 * idx + 42);

            auto t_xyyz_xy = primBuffer.data(90 * idx + 43);

            auto t_xyyz_xz = primBuffer.data(90 * idx + 44);

            auto t_xyyz_yy = primBuffer.data(90 * idx + 45);

            auto t_xyyz_yz = primBuffer.data(90 * idx + 46);

            auto t_xyyz_zz = primBuffer.data(90 * idx + 47);

            auto t_xyzz_xx = primBuffer.data(90 * idx + 48);

            auto t_xyzz_xy = primBuffer.data(90 * idx + 49);

            // Batch of Integrals (40,50)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xy_xx, pa2pb_xy_xy, \
                                     pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_xyy_x, pa2pb_xyy_y, pa2pb_xyy_z, pa2pb_xyyy_yz, \
                                     pa2pb_xyyy_zz, pa2pb_xyyz_xx, pa2pb_xyyz_xy, pa2pb_xyyz_xz, pa2pb_xyyz_yy, \
                                     pa2pb_xyyz_yz, pa2pb_xyyz_zz, pa2pb_xyz_x, pa2pb_xyz_y, pa2pb_xyz_z, pa2pb_xyzz_xx, \
                                     pa2pb_xyzz_xy, pa2pb_xz_xx, pa2pb_xz_xy, pa2pb_xz_xz, pa2pb_xz_yy, pa2pb_xz_yz, \
                                     pa2pb_xz_zz, pa2pb_xzz_x, pa2pb_y_x, pa2pb_y_y, pa2pb_yyz_x, pa2pb_yyz_y, \
                                     pa2pb_yyz_z, pa2pb_yzz_x, pa2pb_yzz_y, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, pa_xy, \
                                     pa_xyyy, pa_xyyz, pa_xyzz, pa_xz, pa_yy, pa_yz, pa_zz, s_0_0, t_xyyy_yz, t_xyyy_zz, \
                                     t_xyyz_xx, t_xyyz_xy, t_xyyz_xz, t_xyyz_yy, t_xyyz_yz, t_xyyz_zz, t_xyzz_xx, \
                                     t_xyzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyyy_yz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl2_fx + 1.5 * pa2pb_xyy_z[j] * fl1_fx + 1.5 * pa2pb_xy_yz[j] * fl1_fx + pa2pb_xyyy_yz[j]);

                t_xyyy_zz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 0.5 * pa_xyyy[j] * fl1_fx + 1.5 * pa2pb_xy_zz[j] * fl1_fx + pa2pb_xyyy_zz[j]);

                t_xyyz_xx[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa_xyyz[j] * fl1_fx + pa2pb_yyz_x[j] * fl1_fx + 0.5 * pa2pb_xz_xx[j] * fl1_fx + pa2pb_xyyz_xx[j]);

                t_xyyz_xy[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + 0.25 * pa2pb_z_y[j] * fl2_fx + pa2pb_xyz_x[j] * fl1_fx + 0.5 * pa2pb_yyz_y[j] * fl1_fx + 0.5 * pa2pb_xz_xy[j] * fl1_fx + pa2pb_xyyz_xy[j]);

                t_xyyz_xz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + 0.25 * pa2pb_x_x[j] * fl2_fx + 0.25 * pa2pb_z_z[j] * fl2_fx + 0.5 * pa2pb_xyy_x[j] * fl1_fx + 0.5 * pa2pb_yyz_z[j] * fl1_fx + 0.5 * pa2pb_xz_xz[j] * fl1_fx + pa2pb_xyyz_xz[j]);

                t_xyyz_yy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 0.5 * pa_xyyz[j] * fl1_fx + 2.0 * pa2pb_xyz_y[j] * fl1_fx + 0.5 * pa2pb_xz_yy[j] * fl1_fx + pa2pb_xyyz_yy[j]);

                t_xyyz_yz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + 0.25 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_xyy_y[j] * fl1_fx + pa2pb_xyz_z[j] * fl1_fx + 0.5 * pa2pb_xz_yz[j] * fl1_fx + pa2pb_xyyz_yz[j]);

                t_xyyz_zz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa_xyyz[j] * fl1_fx + pa2pb_xyy_z[j] * fl1_fx + 0.5 * pa2pb_xz_zz[j] * fl1_fx + pa2pb_xyyz_zz[j]);

                t_xyzz_xx[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa_xyzz[j] * fl1_fx + pa2pb_yzz_x[j] * fl1_fx + 0.5 * pa2pb_xy_xx[j] * fl1_fx + pa2pb_xyzz_xx[j]);

                t_xyzz_xy[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_zz[j] * fl2_fx + 0.25 * pa2pb_x_x[j] * fl2_fx + 0.25 * pa2pb_y_y[j] * fl2_fx + 0.5 * pa2pb_xzz_x[j] * fl1_fx + 0.5 * pa2pb_yzz_y[j] * fl1_fx + 0.5 * pa2pb_xy_xy[j] * fl1_fx + pa2pb_xyzz_xy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD_50_60(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (50,60)

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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(306 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(306 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(306 * idx + 2);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 20);

            auto pa2pb_xy_xz = pa2pbDistances.data(306 * idx + 41);

            auto pa2pb_xy_yy = pa2pbDistances.data(306 * idx + 42);

            auto pa2pb_xy_yz = pa2pbDistances.data(306 * idx + 43);

            auto pa2pb_xy_zz = pa2pbDistances.data(306 * idx + 44);

            auto pa2pb_xz_xx = pa2pbDistances.data(306 * idx + 48);

            auto pa2pb_xz_xy = pa2pbDistances.data(306 * idx + 49);

            auto pa2pb_xz_xz = pa2pbDistances.data(306 * idx + 50);

            auto pa2pb_xz_yy = pa2pbDistances.data(306 * idx + 51);

            auto pa2pb_xz_yz = pa2pbDistances.data(306 * idx + 52);

            auto pa2pb_xz_zz = pa2pbDistances.data(306 * idx + 53);

            auto pa2pb_xyz_x = pa2pbDistances.data(306 * idx + 117);

            auto pa2pb_xyz_y = pa2pbDistances.data(306 * idx + 118);

            auto pa2pb_xyz_z = pa2pbDistances.data(306 * idx + 119);

            auto pa2pb_xzz_x = pa2pbDistances.data(306 * idx + 126);

            auto pa2pb_xzz_y = pa2pbDistances.data(306 * idx + 127);

            auto pa2pb_xzz_z = pa2pbDistances.data(306 * idx + 128);

            auto pa2pb_yzz_z = pa2pbDistances.data(306 * idx + 155);

            auto pa2pb_zzz_x = pa2pbDistances.data(306 * idx + 162);

            auto pa2pb_zzz_y = pa2pbDistances.data(306 * idx + 163);

            auto pa2pb_zzz_z = pa2pbDistances.data(306 * idx + 164);

            auto pa2pb_xyzz_xz = pa2pbDistances.data(306 * idx + 248);

            auto pa2pb_xyzz_yy = pa2pbDistances.data(306 * idx + 249);

            auto pa2pb_xyzz_yz = pa2pbDistances.data(306 * idx + 250);

            auto pa2pb_xyzz_zz = pa2pbDistances.data(306 * idx + 251);

            auto pa2pb_xzzz_xx = pa2pbDistances.data(306 * idx + 255);

            auto pa2pb_xzzz_xy = pa2pbDistances.data(306 * idx + 256);

            auto pa2pb_xzzz_xz = pa2pbDistances.data(306 * idx + 257);

            auto pa2pb_xzzz_yy = pa2pbDistances.data(306 * idx + 258);

            auto pa2pb_xzzz_yz = pa2pbDistances.data(306 * idx + 259);

            auto pa2pb_xzzz_zz = pa2pbDistances.data(306 * idx + 260);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xyzz_xz = primBuffer.data(90 * idx + 50);

            auto t_xyzz_yy = primBuffer.data(90 * idx + 51);

            auto t_xyzz_yz = primBuffer.data(90 * idx + 52);

            auto t_xyzz_zz = primBuffer.data(90 * idx + 53);

            auto t_xzzz_xx = primBuffer.data(90 * idx + 54);

            auto t_xzzz_xy = primBuffer.data(90 * idx + 55);

            auto t_xzzz_xz = primBuffer.data(90 * idx + 56);

            auto t_xzzz_yy = primBuffer.data(90 * idx + 57);

            auto t_xzzz_yz = primBuffer.data(90 * idx + 58);

            auto t_xzzz_zz = primBuffer.data(90 * idx + 59);

            // Batch of Integrals (50,60)

            #pragma omp simd aligned(fx, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xy_xz, pa2pb_xy_yy, \
                                     pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_xyz_x, pa2pb_xyz_y, pa2pb_xyz_z, pa2pb_xyzz_xz, \
                                     pa2pb_xyzz_yy, pa2pb_xyzz_yz, pa2pb_xyzz_zz, pa2pb_xz_xx, pa2pb_xz_xy, pa2pb_xz_xz, \
                                     pa2pb_xz_yy, pa2pb_xz_yz, pa2pb_xz_zz, pa2pb_xzz_x, pa2pb_xzz_y, pa2pb_xzz_z, \
                                     pa2pb_xzzz_xx, pa2pb_xzzz_xy, pa2pb_xzzz_xz, pa2pb_xzzz_yy, pa2pb_xzzz_yz, \
                                     pa2pb_xzzz_zz, pa2pb_y_z, pa2pb_yzz_z, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, pa2pb_zzz_x, \
                                     pa2pb_zzz_y, pa2pb_zzz_z, pa_xy, pa_xyzz, pa_xz, pa_xzzz, pa_yz, pa_zz, s_0_0, t_xyzz_xz, \
                                     t_xyzz_yy, t_xyzz_yz, t_xyzz_zz, t_xzzz_xx, t_xzzz_xy, t_xzzz_xz, t_xzzz_yy, \
                                     t_xzzz_yz, t_xzzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyzz_xz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + 0.25 * pa2pb_y_z[j] * fl2_fx + pa2pb_xyz_x[j] * fl1_fx + 0.5 * pa2pb_yzz_z[j] * fl1_fx + 0.5 * pa2pb_xy_xz[j] * fl1_fx + pa2pb_xyzz_xz[j]);

                t_xyzz_yy[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa_xyzz[j] * fl1_fx + pa2pb_xzz_y[j] * fl1_fx + 0.5 * pa2pb_xy_yy[j] * fl1_fx + pa2pb_xyzz_yy[j]);

                t_xyzz_yz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + 0.25 * pa2pb_x_z[j] * fl2_fx + pa2pb_xyz_y[j] * fl1_fx + 0.5 * pa2pb_xzz_z[j] * fl1_fx + 0.5 * pa2pb_xy_yz[j] * fl1_fx + pa2pb_xyzz_yz[j]);

                t_xyzz_zz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 0.5 * pa_xyzz[j] * fl1_fx + 2.0 * pa2pb_xyz_z[j] * fl1_fx + 0.5 * pa2pb_xy_zz[j] * fl1_fx + pa2pb_xyzz_zz[j]);

                t_xzzz_xx[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 1.5 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa_xzzz[j] * fl1_fx + pa2pb_zzz_x[j] * fl1_fx + 1.5 * pa2pb_xz_xx[j] * fl1_fx + pa2pb_xzzz_xx[j]);

                t_xzzz_xy[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa2pb_zzz_y[j] * fl1_fx + 1.5 * pa2pb_xz_xy[j] * fl1_fx + pa2pb_xzzz_xy[j]);

                t_xzzz_xz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 1.5 * pa2pb_xzz_x[j] * fl1_fx + 0.5 * pa2pb_zzz_z[j] * fl1_fx + 1.5 * pa2pb_xz_xz[j] * fl1_fx + pa2pb_xzzz_xz[j]);

                t_xzzz_yy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 0.5 * pa_xzzz[j] * fl1_fx + 1.5 * pa2pb_xz_yy[j] * fl1_fx + pa2pb_xzzz_yy[j]);

                t_xzzz_yz[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl2_fx + 1.5 * pa2pb_xzz_y[j] * fl1_fx + 1.5 * pa2pb_xz_yz[j] * fl1_fx + pa2pb_xzzz_yz[j]);

                t_xzzz_zz[j] = fl_s_0_0 * (2.25 * pa_xz[j] * fl2_fx + 1.5 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa_xzzz[j] * fl1_fx + 3.0 * pa2pb_xzz_z[j] * fl1_fx + 1.5 * pa2pb_xz_zz[j] * fl1_fx + pa2pb_xzzz_zz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD_60_70(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (60,70)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 19);

            auto pa2pb_yy_xx = pa2pbDistances.data(306 * idx + 57);

            auto pa2pb_yy_xy = pa2pbDistances.data(306 * idx + 58);

            auto pa2pb_yy_xz = pa2pbDistances.data(306 * idx + 59);

            auto pa2pb_yy_yy = pa2pbDistances.data(306 * idx + 60);

            auto pa2pb_yy_yz = pa2pbDistances.data(306 * idx + 61);

            auto pa2pb_yy_zz = pa2pbDistances.data(306 * idx + 62);

            auto pa2pb_yz_xx = pa2pbDistances.data(306 * idx + 66);

            auto pa2pb_yz_xy = pa2pbDistances.data(306 * idx + 67);

            auto pa2pb_yz_xz = pa2pbDistances.data(306 * idx + 68);

            auto pa2pb_yz_yy = pa2pbDistances.data(306 * idx + 69);

            auto pa2pb_yyy_x = pa2pbDistances.data(306 * idx + 135);

            auto pa2pb_yyy_y = pa2pbDistances.data(306 * idx + 136);

            auto pa2pb_yyy_z = pa2pbDistances.data(306 * idx + 137);

            auto pa2pb_yyz_x = pa2pbDistances.data(306 * idx + 144);

            auto pa2pb_yyz_y = pa2pbDistances.data(306 * idx + 145);

            auto pa2pb_yyyy_xx = pa2pbDistances.data(306 * idx + 264);

            auto pa2pb_yyyy_xy = pa2pbDistances.data(306 * idx + 265);

            auto pa2pb_yyyy_xz = pa2pbDistances.data(306 * idx + 266);

            auto pa2pb_yyyy_yy = pa2pbDistances.data(306 * idx + 267);

            auto pa2pb_yyyy_yz = pa2pbDistances.data(306 * idx + 268);

            auto pa2pb_yyyy_zz = pa2pbDistances.data(306 * idx + 269);

            auto pa2pb_yyyz_xx = pa2pbDistances.data(306 * idx + 273);

            auto pa2pb_yyyz_xy = pa2pbDistances.data(306 * idx + 274);

            auto pa2pb_yyyz_xz = pa2pbDistances.data(306 * idx + 275);

            auto pa2pb_yyyz_yy = pa2pbDistances.data(306 * idx + 276);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yyyy_xx = primBuffer.data(90 * idx + 60);

            auto t_yyyy_xy = primBuffer.data(90 * idx + 61);

            auto t_yyyy_xz = primBuffer.data(90 * idx + 62);

            auto t_yyyy_yy = primBuffer.data(90 * idx + 63);

            auto t_yyyy_yz = primBuffer.data(90 * idx + 64);

            auto t_yyyy_zz = primBuffer.data(90 * idx + 65);

            auto t_yyyz_xx = primBuffer.data(90 * idx + 66);

            auto t_yyyz_xy = primBuffer.data(90 * idx + 67);

            auto t_yyyz_xz = primBuffer.data(90 * idx + 68);

            auto t_yyyz_yy = primBuffer.data(90 * idx + 69);

            // Batch of Integrals (60,70)

            #pragma omp simd aligned(fx, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xy, \
                                     pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, pa2pb_yy_zz, pa2pb_yyy_x, pa2pb_yyy_y, \
                                     pa2pb_yyy_z, pa2pb_yyyy_xx, pa2pb_yyyy_xy, pa2pb_yyyy_xz, pa2pb_yyyy_yy, \
                                     pa2pb_yyyy_yz, pa2pb_yyyy_zz, pa2pb_yyyz_xx, pa2pb_yyyz_xy, pa2pb_yyyz_xz, \
                                     pa2pb_yyyz_yy, pa2pb_yyz_x, pa2pb_yyz_y, pa2pb_yz_xx, pa2pb_yz_xy, pa2pb_yz_xz, \
                                     pa2pb_yz_yy, pa2pb_z_x, pa2pb_z_y, pa_yy, pa_yyyy, pa_yyyz, pa_yz, pb_xx, pb_xy, pb_xz, \
                                     pb_yy, pb_yz, pb_zz, s_0_0, t_yyyy_xx, t_yyyy_xy, t_yyyy_xz, t_yyyy_yy, t_yyyy_yz, \
                                     t_yyyy_zz, t_yyyz_xx, t_yyyz_xy, t_yyyz_xz, t_yyyz_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yyyy_xx[j] = fl_s_0_0 * (0.375 * fl3_fx + 1.5 * pa_yy[j] * fl2_fx + 0.5 * pa_yyyy[j] * fl1_fx + 0.75 * pb_xx[j] * fl2_fx + 3.0 * pa2pb_yy_xx[j] * fl1_fx + pa2pb_yyyy_xx[j]);

                t_yyyy_xy[j] = fl_s_0_0 * (3.0 * pa2pb_y_x[j] * fl2_fx + 2.0 * pa2pb_yyy_x[j] * fl1_fx + 0.75 * pb_xy[j] * fl2_fx + 3.0 * pa2pb_yy_xy[j] * fl1_fx + pa2pb_yyyy_xy[j]);

                t_yyyy_xz[j] = fl_s_0_0 * (0.75 * pb_xz[j] * fl2_fx + 3.0 * pa2pb_yy_xz[j] * fl1_fx + pa2pb_yyyy_xz[j]);

                t_yyyy_yy[j] = fl_s_0_0 * (1.875 * fl3_fx + 4.5 * pa_yy[j] * fl2_fx + 6.0 * pa2pb_y_y[j] * fl2_fx + 0.5 * pa_yyyy[j] * fl1_fx + 4.0 * pa2pb_yyy_y[j] * fl1_fx + 0.75 * pb_yy[j] * fl2_fx + 3.0 * pa2pb_yy_yy[j] * fl1_fx + pa2pb_yyyy_yy[j]);

                t_yyyy_yz[j] = fl_s_0_0 * (3.0 * pa2pb_y_z[j] * fl2_fx + 2.0 * pa2pb_yyy_z[j] * fl1_fx + 0.75 * pb_yz[j] * fl2_fx + 3.0 * pa2pb_yy_yz[j] * fl1_fx + pa2pb_yyyy_yz[j]);

                t_yyyy_zz[j] = fl_s_0_0 * (0.375 * fl3_fx + 1.5 * pa_yy[j] * fl2_fx + 0.5 * pa_yyyy[j] * fl1_fx + 0.75 * pb_zz[j] * fl2_fx + 3.0 * pa2pb_yy_zz[j] * fl1_fx + pa2pb_yyyy_zz[j]);

                t_yyyz_xx[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 0.5 * pa_yyyz[j] * fl1_fx + 1.5 * pa2pb_yz_xx[j] * fl1_fx + pa2pb_yyyz_xx[j]);

                t_yyyz_xy[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl2_fx + 1.5 * pa2pb_yyz_x[j] * fl1_fx + 1.5 * pa2pb_yz_xy[j] * fl1_fx + pa2pb_yyyz_xy[j]);

                t_yyyz_xz[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa2pb_yyy_x[j] * fl1_fx + 1.5 * pa2pb_yz_xz[j] * fl1_fx + pa2pb_yyyz_xz[j]);

                t_yyyz_yy[j] = fl_s_0_0 * (2.25 * pa_yz[j] * fl2_fx + 1.5 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa_yyyz[j] * fl1_fx + 3.0 * pa2pb_yyz_y[j] * fl1_fx + 1.5 * pa2pb_yz_yy[j] * fl1_fx + pa2pb_yyyz_yy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD_70_80(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (70,80)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 20);

            auto pa2pb_yy_xx = pa2pbDistances.data(306 * idx + 57);

            auto pa2pb_yy_xy = pa2pbDistances.data(306 * idx + 58);

            auto pa2pb_yy_xz = pa2pbDistances.data(306 * idx + 59);

            auto pa2pb_yy_yy = pa2pbDistances.data(306 * idx + 60);

            auto pa2pb_yy_yz = pa2pbDistances.data(306 * idx + 61);

            auto pa2pb_yy_zz = pa2pbDistances.data(306 * idx + 62);

            auto pa2pb_yz_xx = pa2pbDistances.data(306 * idx + 66);

            auto pa2pb_yz_xy = pa2pbDistances.data(306 * idx + 67);

            auto pa2pb_yz_yz = pa2pbDistances.data(306 * idx + 70);

            auto pa2pb_yz_zz = pa2pbDistances.data(306 * idx + 71);

            auto pa2pb_zz_xx = pa2pbDistances.data(306 * idx + 75);

            auto pa2pb_zz_xy = pa2pbDistances.data(306 * idx + 76);

            auto pa2pb_zz_xz = pa2pbDistances.data(306 * idx + 77);

            auto pa2pb_zz_yy = pa2pbDistances.data(306 * idx + 78);

            auto pa2pb_zz_yz = pa2pbDistances.data(306 * idx + 79);

            auto pa2pb_zz_zz = pa2pbDistances.data(306 * idx + 80);

            auto pa2pb_yyy_y = pa2pbDistances.data(306 * idx + 136);

            auto pa2pb_yyy_z = pa2pbDistances.data(306 * idx + 137);

            auto pa2pb_yyz_x = pa2pbDistances.data(306 * idx + 144);

            auto pa2pb_yyz_y = pa2pbDistances.data(306 * idx + 145);

            auto pa2pb_yyz_z = pa2pbDistances.data(306 * idx + 146);

            auto pa2pb_yzz_x = pa2pbDistances.data(306 * idx + 153);

            auto pa2pb_yzz_y = pa2pbDistances.data(306 * idx + 154);

            auto pa2pb_yzz_z = pa2pbDistances.data(306 * idx + 155);

            auto pa2pb_zzz_x = pa2pbDistances.data(306 * idx + 162);

            auto pa2pb_yyyz_yz = pa2pbDistances.data(306 * idx + 277);

            auto pa2pb_yyyz_zz = pa2pbDistances.data(306 * idx + 278);

            auto pa2pb_yyzz_xx = pa2pbDistances.data(306 * idx + 282);

            auto pa2pb_yyzz_xy = pa2pbDistances.data(306 * idx + 283);

            auto pa2pb_yyzz_xz = pa2pbDistances.data(306 * idx + 284);

            auto pa2pb_yyzz_yy = pa2pbDistances.data(306 * idx + 285);

            auto pa2pb_yyzz_yz = pa2pbDistances.data(306 * idx + 286);

            auto pa2pb_yyzz_zz = pa2pbDistances.data(306 * idx + 287);

            auto pa2pb_yzzz_xx = pa2pbDistances.data(306 * idx + 291);

            auto pa2pb_yzzz_xy = pa2pbDistances.data(306 * idx + 292);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yyyz_yz = primBuffer.data(90 * idx + 70);

            auto t_yyyz_zz = primBuffer.data(90 * idx + 71);

            auto t_yyzz_xx = primBuffer.data(90 * idx + 72);

            auto t_yyzz_xy = primBuffer.data(90 * idx + 73);

            auto t_yyzz_xz = primBuffer.data(90 * idx + 74);

            auto t_yyzz_yy = primBuffer.data(90 * idx + 75);

            auto t_yyzz_yz = primBuffer.data(90 * idx + 76);

            auto t_yyzz_zz = primBuffer.data(90 * idx + 77);

            auto t_yzzz_xx = primBuffer.data(90 * idx + 78);

            auto t_yzzz_xy = primBuffer.data(90 * idx + 79);

            // Batch of Integrals (70,80)

            #pragma omp simd aligned(fx, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xy, \
                                     pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, pa2pb_yy_zz, pa2pb_yyy_y, pa2pb_yyy_z, \
                                     pa2pb_yyyz_yz, pa2pb_yyyz_zz, pa2pb_yyz_x, pa2pb_yyz_y, pa2pb_yyz_z, pa2pb_yyzz_xx, \
                                     pa2pb_yyzz_xy, pa2pb_yyzz_xz, pa2pb_yyzz_yy, pa2pb_yyzz_yz, pa2pb_yyzz_zz, \
                                     pa2pb_yz_xx, pa2pb_yz_xy, pa2pb_yz_yz, pa2pb_yz_zz, pa2pb_yzz_x, pa2pb_yzz_y, \
                                     pa2pb_yzz_z, pa2pb_yzzz_xx, pa2pb_yzzz_xy, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, \
                                     pa2pb_zz_xx, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zz_zz, \
                                     pa2pb_zzz_x, pa_yy, pa_yyyz, pa_yyzz, pa_yz, pa_yzzz, pa_zz, pb_xx, pb_xy, pb_xz, pb_yy, \
                                     pb_yz, pb_zz, s_0_0, t_yyyz_yz, t_yyyz_zz, t_yyzz_xx, t_yyzz_xy, t_yyzz_xz, \
                                     t_yyzz_yy, t_yyzz_yz, t_yyzz_zz, t_yzzz_xx, t_yzzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yyyz_yz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.5 * pa2pb_yyy_y[j] * fl1_fx + 1.5 * pa2pb_yyz_z[j] * fl1_fx + 1.5 * pa2pb_yz_yz[j] * fl1_fx + pa2pb_yyyz_yz[j]);

                t_yyyz_zz[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 1.5 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa_yyyz[j] * fl1_fx + pa2pb_yyy_z[j] * fl1_fx + 1.5 * pa2pb_yz_zz[j] * fl1_fx + pa2pb_yyyz_zz[j]);

                t_yyzz_xx[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + 0.25 * pa_zz[j] * fl2_fx + 0.5 * pa_yyzz[j] * fl1_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_yy_xx[j] * fl1_fx + 0.5 * pa2pb_zz_xx[j] * fl1_fx + pa2pb_yyzz_xx[j]);

                t_yyzz_xy[j] = fl_s_0_0 * (0.5 * pa2pb_y_x[j] * fl2_fx + pa2pb_yzz_x[j] * fl1_fx + 0.25 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_yy_xy[j] * fl1_fx + 0.5 * pa2pb_zz_xy[j] * fl1_fx + pa2pb_yyzz_xy[j]);

                t_yyzz_xz[j] = fl_s_0_0 * (0.5 * pa2pb_z_x[j] * fl2_fx + pa2pb_yyz_x[j] * fl1_fx + 0.25 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_yy_xz[j] * fl1_fx + 0.5 * pa2pb_zz_xz[j] * fl1_fx + pa2pb_yyzz_xz[j]);

                t_yyzz_yy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.25 * pa_yy[j] * fl2_fx + pa2pb_y_y[j] * fl2_fx + 0.5 * pa_yyzz[j] * fl1_fx + 2.0 * pa2pb_yzz_y[j] * fl1_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_yy_yy[j] * fl1_fx + 0.5 * pa2pb_zz_yy[j] * fl1_fx + pa2pb_yyzz_yy[j]);

                t_yyzz_yz[j] = fl_s_0_0 * (pa_yz[j] * fl2_fx + 0.5 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa2pb_z_y[j] * fl2_fx + pa2pb_yyz_y[j] * fl1_fx + pa2pb_yzz_z[j] * fl1_fx + 0.25 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_yy_yz[j] * fl1_fx + 0.5 * pa2pb_zz_yz[j] * fl1_fx + pa2pb_yyzz_yz[j]);

                t_yyzz_zz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.25 * pa_zz[j] * fl2_fx + pa2pb_z_z[j] * fl2_fx + 0.5 * pa_yyzz[j] * fl1_fx + 2.0 * pa2pb_yyz_z[j] * fl1_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_yy_zz[j] * fl1_fx + 0.5 * pa2pb_zz_zz[j] * fl1_fx + pa2pb_yyzz_zz[j]);

                t_yzzz_xx[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 0.5 * pa_yzzz[j] * fl1_fx + 1.5 * pa2pb_yz_xx[j] * fl1_fx + pa2pb_yzzz_xx[j]);

                t_yzzz_xy[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa2pb_zzz_x[j] * fl1_fx + 1.5 * pa2pb_yz_xy[j] * fl1_fx + pa2pb_yzzz_xy[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGD_80_90(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (80,90)

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

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(306 * idx + 9);

            auto pa2pb_y_y = pa2pbDistances.data(306 * idx + 10);

            auto pa2pb_y_z = pa2pbDistances.data(306 * idx + 11);

            auto pa2pb_z_x = pa2pbDistances.data(306 * idx + 18);

            auto pa2pb_z_y = pa2pbDistances.data(306 * idx + 19);

            auto pa2pb_z_z = pa2pbDistances.data(306 * idx + 20);

            auto pa2pb_yz_xz = pa2pbDistances.data(306 * idx + 68);

            auto pa2pb_yz_yy = pa2pbDistances.data(306 * idx + 69);

            auto pa2pb_yz_yz = pa2pbDistances.data(306 * idx + 70);

            auto pa2pb_yz_zz = pa2pbDistances.data(306 * idx + 71);

            auto pa2pb_zz_xx = pa2pbDistances.data(306 * idx + 75);

            auto pa2pb_zz_xy = pa2pbDistances.data(306 * idx + 76);

            auto pa2pb_zz_xz = pa2pbDistances.data(306 * idx + 77);

            auto pa2pb_zz_yy = pa2pbDistances.data(306 * idx + 78);

            auto pa2pb_zz_yz = pa2pbDistances.data(306 * idx + 79);

            auto pa2pb_zz_zz = pa2pbDistances.data(306 * idx + 80);

            auto pa2pb_yzz_x = pa2pbDistances.data(306 * idx + 153);

            auto pa2pb_yzz_y = pa2pbDistances.data(306 * idx + 154);

            auto pa2pb_yzz_z = pa2pbDistances.data(306 * idx + 155);

            auto pa2pb_zzz_x = pa2pbDistances.data(306 * idx + 162);

            auto pa2pb_zzz_y = pa2pbDistances.data(306 * idx + 163);

            auto pa2pb_zzz_z = pa2pbDistances.data(306 * idx + 164);

            auto pa2pb_yzzz_xz = pa2pbDistances.data(306 * idx + 293);

            auto pa2pb_yzzz_yy = pa2pbDistances.data(306 * idx + 294);

            auto pa2pb_yzzz_yz = pa2pbDistances.data(306 * idx + 295);

            auto pa2pb_yzzz_zz = pa2pbDistances.data(306 * idx + 296);

            auto pa2pb_zzzz_xx = pa2pbDistances.data(306 * idx + 300);

            auto pa2pb_zzzz_xy = pa2pbDistances.data(306 * idx + 301);

            auto pa2pb_zzzz_xz = pa2pbDistances.data(306 * idx + 302);

            auto pa2pb_zzzz_yy = pa2pbDistances.data(306 * idx + 303);

            auto pa2pb_zzzz_yz = pa2pbDistances.data(306 * idx + 304);

            auto pa2pb_zzzz_zz = pa2pbDistances.data(306 * idx + 305);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yzzz_xz = primBuffer.data(90 * idx + 80);

            auto t_yzzz_yy = primBuffer.data(90 * idx + 81);

            auto t_yzzz_yz = primBuffer.data(90 * idx + 82);

            auto t_yzzz_zz = primBuffer.data(90 * idx + 83);

            auto t_zzzz_xx = primBuffer.data(90 * idx + 84);

            auto t_zzzz_xy = primBuffer.data(90 * idx + 85);

            auto t_zzzz_xz = primBuffer.data(90 * idx + 86);

            auto t_zzzz_yy = primBuffer.data(90 * idx + 87);

            auto t_zzzz_yz = primBuffer.data(90 * idx + 88);

            auto t_zzzz_zz = primBuffer.data(90 * idx + 89);

            // Batch of Integrals (80,90)

            #pragma omp simd aligned(fx, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa2pb_yz_xz, pa2pb_yz_yy, \
                                     pa2pb_yz_yz, pa2pb_yz_zz, pa2pb_yzz_x, pa2pb_yzz_y, pa2pb_yzz_z, pa2pb_yzzz_xz, \
                                     pa2pb_yzzz_yy, pa2pb_yzzz_yz, pa2pb_yzzz_zz, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, \
                                     pa2pb_zz_xx, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zz_zz, \
                                     pa2pb_zzz_x, pa2pb_zzz_y, pa2pb_zzz_z, pa2pb_zzzz_xx, pa2pb_zzzz_xy, \
                                     pa2pb_zzzz_xz, pa2pb_zzzz_yy, pa2pb_zzzz_yz, pa2pb_zzzz_zz, pa_yz, pa_yzzz, pa_zz, \
                                     pa_zzzz, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, t_yzzz_xz, t_yzzz_yy, \
                                     t_yzzz_yz, t_yzzz_zz, t_zzzz_xx, t_zzzz_xy, t_zzzz_xz, t_zzzz_yy, t_zzzz_yz, \
                                     t_zzzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yzzz_xz[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl2_fx + 1.5 * pa2pb_yzz_x[j] * fl1_fx + 1.5 * pa2pb_yz_xz[j] * fl1_fx + pa2pb_yzzz_xz[j]);

                t_yzzz_yy[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 1.5 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa_yzzz[j] * fl1_fx + pa2pb_zzz_y[j] * fl1_fx + 1.5 * pa2pb_yz_yy[j] * fl1_fx + pa2pb_yzzz_yy[j]);

                t_yzzz_yz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 1.5 * pa2pb_yzz_y[j] * fl1_fx + 0.5 * pa2pb_zzz_z[j] * fl1_fx + 1.5 * pa2pb_yz_yz[j] * fl1_fx + pa2pb_yzzz_yz[j]);

                t_yzzz_zz[j] = fl_s_0_0 * (2.25 * pa_yz[j] * fl2_fx + 1.5 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa_yzzz[j] * fl1_fx + 3.0 * pa2pb_yzz_z[j] * fl1_fx + 1.5 * pa2pb_yz_zz[j] * fl1_fx + pa2pb_yzzz_zz[j]);

                t_zzzz_xx[j] = fl_s_0_0 * (0.375 * fl3_fx + 1.5 * pa_zz[j] * fl2_fx + 0.5 * pa_zzzz[j] * fl1_fx + 0.75 * pb_xx[j] * fl2_fx + 3.0 * pa2pb_zz_xx[j] * fl1_fx + pa2pb_zzzz_xx[j]);

                t_zzzz_xy[j] = fl_s_0_0 * (0.75 * pb_xy[j] * fl2_fx + 3.0 * pa2pb_zz_xy[j] * fl1_fx + pa2pb_zzzz_xy[j]);

                t_zzzz_xz[j] = fl_s_0_0 * (3.0 * pa2pb_z_x[j] * fl2_fx + 2.0 * pa2pb_zzz_x[j] * fl1_fx + 0.75 * pb_xz[j] * fl2_fx + 3.0 * pa2pb_zz_xz[j] * fl1_fx + pa2pb_zzzz_xz[j]);

                t_zzzz_yy[j] = fl_s_0_0 * (0.375 * fl3_fx + 1.5 * pa_zz[j] * fl2_fx + 0.5 * pa_zzzz[j] * fl1_fx + 0.75 * pb_yy[j] * fl2_fx + 3.0 * pa2pb_zz_yy[j] * fl1_fx + pa2pb_zzzz_yy[j]);

                t_zzzz_yz[j] = fl_s_0_0 * (3.0 * pa2pb_z_y[j] * fl2_fx + 2.0 * pa2pb_zzz_y[j] * fl1_fx + 0.75 * pb_yz[j] * fl2_fx + 3.0 * pa2pb_zz_yz[j] * fl1_fx + pa2pb_zzzz_yz[j]);

                t_zzzz_zz[j] = fl_s_0_0 * (1.875 * fl3_fx + 4.5 * pa_zz[j] * fl2_fx + 6.0 * pa2pb_z_z[j] * fl2_fx + 0.5 * pa_zzzz[j] * fl1_fx + 4.0 * pa2pb_zzz_z[j] * fl1_fx + 0.75 * pb_zz[j] * fl2_fx + 3.0 * pa2pb_zz_zz[j] * fl1_fx + pa2pb_zzzz_zz[j]);
            }

            idx++;
        }
    }


} // ovlrecfunc namespace

