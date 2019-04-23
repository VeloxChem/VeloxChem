//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "NuclearPotentialRecFuncForPX.hpp"

namespace npotrecfunc { // npotrecfunc namespace

    void
    compNuclearPotentialForPP(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xy, pc_xz, pc_y, pc_yy, pc_yz, \
                                     pc_z, pc_zz, s_0_0_0, s_0_0_1, s_0_0_2, t_x_x, t_x_y, t_x_z, t_y_x, t_y_y, t_y_z, \
                                     t_z_x, t_z_y, t_z_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl1_fx = fx[j];

                t_x_x[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pa_x[j] * pb_x[j]);

                t_x_x[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - pa_x[j] * pc_x[j] - pc_x[j] * pb_x[j]);

                t_x_x[j] += fl_s_0_0_2 * pc_xx[j];

                t_x_y[j] = fl_s_0_0_0 * pa_x[j] * pb_y[j];

                t_x_y[j] += fl_s_0_0_1 * (-pa_x[j] * pc_y[j] - pc_x[j] * pb_y[j]);

                t_x_y[j] += fl_s_0_0_2 * pc_xy[j];

                t_x_z[j] = fl_s_0_0_0 * pa_x[j] * pb_z[j];

                t_x_z[j] += fl_s_0_0_1 * (-pa_x[j] * pc_z[j] - pc_x[j] * pb_z[j]);

                t_x_z[j] += fl_s_0_0_2 * pc_xz[j];

                t_y_x[j] = fl_s_0_0_0 * pa_y[j] * pb_x[j];

                t_y_x[j] += fl_s_0_0_1 * (-pa_y[j] * pc_x[j] - pc_y[j] * pb_x[j]);

                t_y_x[j] += fl_s_0_0_2 * pc_xy[j];

                t_y_y[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pa_y[j] * pb_y[j]);

                t_y_y[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - pa_y[j] * pc_y[j] - pc_y[j] * pb_y[j]);

                t_y_y[j] += fl_s_0_0_2 * pc_yy[j];

                t_y_z[j] = fl_s_0_0_0 * pa_y[j] * pb_z[j];

                t_y_z[j] += fl_s_0_0_1 * (-pa_y[j] * pc_z[j] - pc_y[j] * pb_z[j]);

                t_y_z[j] += fl_s_0_0_2 * pc_yz[j];

                t_z_x[j] = fl_s_0_0_0 * pa_z[j] * pb_x[j];

                t_z_x[j] += fl_s_0_0_1 * (-pa_z[j] * pc_x[j] - pc_z[j] * pb_x[j]);

                t_z_x[j] += fl_s_0_0_2 * pc_xz[j];

                t_z_y[j] = fl_s_0_0_0 * pa_z[j] * pb_y[j];

                t_z_y[j] += fl_s_0_0_1 * (-pa_z[j] * pc_y[j] - pc_z[j] * pb_y[j]);

                t_z_y[j] += fl_s_0_0_2 * pc_yz[j];

                t_z_z[j] = fl_s_0_0_0 * (0.5 * fl1_fx + pa_z[j] * pb_z[j]);

                t_z_z[j] += fl_s_0_0_1 * (-0.5 * fl1_fx - pa_z[j] * pc_z[j] - pc_z[j] * pb_z[j]);

                t_z_z[j] += fl_s_0_0_2 * pc_zz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPD(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForPD_0_9(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPD_9_18(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                    braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compNuclearPotentialForPD_0_9(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

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

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(4 * idx);

            auto s_0_0_1 = auxBuffer.data(4 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(4 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(4 * idx + 3);

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, pc_x, \
                                     pc_xx, pc_xxx, pc_xxy, pc_xxz, pc_xy, pc_xyy, pc_xyz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yz, \
                                     pc_z, pc_zz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, t_x_xx, t_x_xy, t_x_xz, t_x_yy, \
                                     t_x_yz, t_x_zz, t_y_xx, t_y_xy, t_y_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl1_fx = fx[j];

                t_x_xx[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx + fl1_fx * pb_x[j] + pa_x[j] * pb_xx[j]);

                t_x_xx[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx - fl1_fx * pb_x[j] - 2.0 * pa_x[j] * pb_x[j] * pc_x[j] - pc_x[j] * pb_xx[j]);

                t_x_xx[j] += fl_s_0_0_2 * (1.5 * pc_x[j] * fl1_fx + pa_x[j] * pc_xx[j] + 2.0 * pc_xx[j] * pb_x[j]);

                t_x_xx[j] += -fl_s_0_0_3 * pc_xxx[j];

                t_x_xy[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_y[j] + pa_x[j] * pb_xy[j]);

                t_x_xy[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_y[j] - 0.5 * fl1_fx * pb_y[j] - pa_x[j] * pb_x[j] * pc_y[j] - pa_x[j] * pc_x[j] * pb_y[j] - pc_x[j] * pb_xy[j]);

                t_x_xy[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_y[j] + pa_x[j] * pc_xy[j] + pc_xy[j] * pb_x[j] + pc_xx[j] * pb_y[j]);

                t_x_xy[j] += -fl_s_0_0_3 * pc_xxy[j];

                t_x_xz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_z[j] + pa_x[j] * pb_xz[j]);

                t_x_xz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pb_z[j] - pa_x[j] * pb_x[j] * pc_z[j] - pa_x[j] * pc_x[j] * pb_z[j] - pc_x[j] * pb_xz[j]);

                t_x_xz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + pa_x[j] * pc_xz[j] + pc_xz[j] * pb_x[j] + pc_xx[j] * pb_z[j]);

                t_x_xz[j] += -fl_s_0_0_3 * pc_xxz[j];

                t_x_yy[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx + pa_x[j] * pb_yy[j]);

                t_x_yy[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx - 2.0 * pa_x[j] * pb_y[j] * pc_y[j] - pc_x[j] * pb_yy[j]);

                t_x_yy[j] += fl_s_0_0_2 * (0.5 * pc_x[j] * fl1_fx + pa_x[j] * pc_yy[j] + 2.0 * pc_xy[j] * pb_y[j]);

                t_x_yy[j] += -fl_s_0_0_3 * pc_xyy[j];

                t_x_yz[j] = fl_s_0_0_0 * pa_x[j] * pb_yz[j];

                t_x_yz[j] += fl_s_0_0_1 * (-pa_x[j] * pb_y[j] * pc_z[j] - pa_x[j] * pc_y[j] * pb_z[j] - pc_x[j] * pb_yz[j]);

                t_x_yz[j] += fl_s_0_0_2 * (pa_x[j] * pc_yz[j] + pc_xz[j] * pb_y[j] + pc_xy[j] * pb_z[j]);

                t_x_yz[j] += -fl_s_0_0_3 * pc_xyz[j];

                t_x_zz[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx + pa_x[j] * pb_zz[j]);

                t_x_zz[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx - 2.0 * pa_x[j] * pb_z[j] * pc_z[j] - pc_x[j] * pb_zz[j]);

                t_x_zz[j] += fl_s_0_0_2 * (0.5 * pc_x[j] * fl1_fx + pa_x[j] * pc_zz[j] + 2.0 * pc_xz[j] * pb_z[j]);

                t_x_zz[j] += -fl_s_0_0_3 * pc_xzz[j];

                t_y_xx[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl1_fx + pa_y[j] * pb_xx[j]);

                t_y_xx[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl1_fx - 0.5 * pc_y[j] * fl1_fx - 2.0 * pa_y[j] * pb_x[j] * pc_x[j] - pc_y[j] * pb_xx[j]);

                t_y_xx[j] += fl_s_0_0_2 * (0.5 * pc_y[j] * fl1_fx + pa_y[j] * pc_xx[j] + 2.0 * pc_xy[j] * pb_x[j]);

                t_y_xx[j] += -fl_s_0_0_3 * pc_xxy[j];

                t_y_xy[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_x[j] + pa_y[j] * pb_xy[j]);

                t_y_xy[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_x[j] - 0.5 * fl1_fx * pb_x[j] - pa_y[j] * pb_x[j] * pc_y[j] - pa_y[j] * pc_x[j] * pb_y[j] - pc_y[j] * pb_xy[j]);

                t_y_xy[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_x[j] + pa_y[j] * pc_xy[j] + pc_yy[j] * pb_x[j] + pc_xy[j] * pb_y[j]);

                t_y_xy[j] += -fl_s_0_0_3 * pc_xyy[j];

                t_y_xz[j] = fl_s_0_0_0 * pa_y[j] * pb_xz[j];

                t_y_xz[j] += fl_s_0_0_1 * (-pa_y[j] * pb_x[j] * pc_z[j] - pa_y[j] * pc_x[j] * pb_z[j] - pc_y[j] * pb_xz[j]);

                t_y_xz[j] += fl_s_0_0_2 * (pa_y[j] * pc_xz[j] + pc_yz[j] * pb_x[j] + pc_xy[j] * pb_z[j]);

                t_y_xz[j] += -fl_s_0_0_3 * pc_xyz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPD_9_18(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

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

            auto pc_xxz = pcDistances.data(19 * idx + 11);

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

            #pragma omp simd aligned(fx, pa_y, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, pc_x, \
                                     pc_xx, pc_xxz, pc_xy, pc_xyz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yyy, pc_yyz, pc_yz, pc_yzz, \
                                     pc_z, pc_zz, pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, t_y_yy, t_y_yz, t_y_zz, \
                                     t_z_xx, t_z_xy, t_z_xz, t_z_yy, t_z_yz, t_z_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl1_fx = fx[j];

                t_y_yy[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl1_fx + fl1_fx * pb_y[j] + pa_y[j] * pb_yy[j]);

                t_y_yy[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx - fl1_fx * pb_y[j] - 2.0 * pa_y[j] * pb_y[j] * pc_y[j] - pc_y[j] * pb_yy[j]);

                t_y_yy[j] += fl_s_0_0_2 * (1.5 * pc_y[j] * fl1_fx + pa_y[j] * pc_yy[j] + 2.0 * pc_yy[j] * pb_y[j]);

                t_y_yy[j] += -fl_s_0_0_3 * pc_yyy[j];

                t_y_yz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_z[j] + pa_y[j] * pb_yz[j]);

                t_y_yz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pb_z[j] - pa_y[j] * pb_y[j] * pc_z[j] - pa_y[j] * pc_y[j] * pb_z[j] - pc_y[j] * pb_yz[j]);

                t_y_yz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + pa_y[j] * pc_yz[j] + pc_yz[j] * pb_y[j] + pc_yy[j] * pb_z[j]);

                t_y_yz[j] += -fl_s_0_0_3 * pc_yyz[j];

                t_y_zz[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl1_fx + pa_y[j] * pb_zz[j]);

                t_y_zz[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl1_fx - 0.5 * pc_y[j] * fl1_fx - 2.0 * pa_y[j] * pb_z[j] * pc_z[j] - pc_y[j] * pb_zz[j]);

                t_y_zz[j] += fl_s_0_0_2 * (0.5 * pc_y[j] * fl1_fx + pa_y[j] * pc_zz[j] + 2.0 * pc_yz[j] * pb_z[j]);

                t_y_zz[j] += -fl_s_0_0_3 * pc_yzz[j];

                t_z_xx[j] = fl_s_0_0_0 * (0.5 * pa_z[j] * fl1_fx + pa_z[j] * pb_xx[j]);

                t_z_xx[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * fl1_fx - 0.5 * pc_z[j] * fl1_fx - 2.0 * pa_z[j] * pb_x[j] * pc_x[j] - pc_z[j] * pb_xx[j]);

                t_z_xx[j] += fl_s_0_0_2 * (0.5 * pc_z[j] * fl1_fx + pa_z[j] * pc_xx[j] + 2.0 * pc_xz[j] * pb_x[j]);

                t_z_xx[j] += -fl_s_0_0_3 * pc_xxz[j];

                t_z_xy[j] = fl_s_0_0_0 * pa_z[j] * pb_xy[j];

                t_z_xy[j] += fl_s_0_0_1 * (-pa_z[j] * pb_x[j] * pc_y[j] - pa_z[j] * pc_x[j] * pb_y[j] - pc_z[j] * pb_xy[j]);

                t_z_xy[j] += fl_s_0_0_2 * (pa_z[j] * pc_xy[j] + pc_yz[j] * pb_x[j] + pc_xz[j] * pb_y[j]);

                t_z_xy[j] += -fl_s_0_0_3 * pc_xyz[j];

                t_z_xz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_x[j] + pa_z[j] * pb_xz[j]);

                t_z_xz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_x[j] - 0.5 * fl1_fx * pb_x[j] - pa_z[j] * pb_x[j] * pc_z[j] - pa_z[j] * pc_x[j] * pb_z[j] - pc_z[j] * pb_xz[j]);

                t_z_xz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_x[j] + pa_z[j] * pc_xz[j] + pc_zz[j] * pb_x[j] + pc_xz[j] * pb_z[j]);

                t_z_xz[j] += -fl_s_0_0_3 * pc_xzz[j];

                t_z_yy[j] = fl_s_0_0_0 * (0.5 * pa_z[j] * fl1_fx + pa_z[j] * pb_yy[j]);

                t_z_yy[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * fl1_fx - 0.5 * pc_z[j] * fl1_fx - 2.0 * pa_z[j] * pb_y[j] * pc_y[j] - pc_z[j] * pb_yy[j]);

                t_z_yy[j] += fl_s_0_0_2 * (0.5 * pc_z[j] * fl1_fx + pa_z[j] * pc_yy[j] + 2.0 * pc_yz[j] * pb_y[j]);

                t_z_yy[j] += -fl_s_0_0_3 * pc_yyz[j];

                t_z_yz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_y[j] + pa_z[j] * pb_yz[j]);

                t_z_yz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_y[j] - 0.5 * fl1_fx * pb_y[j] - pa_z[j] * pb_y[j] * pc_z[j] - pa_z[j] * pc_y[j] * pb_z[j] - pc_z[j] * pb_yz[j]);

                t_z_yz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_y[j] + pa_z[j] * pc_yz[j] + pc_zz[j] * pb_y[j] + pc_yz[j] * pb_z[j]);

                t_z_yz[j] += -fl_s_0_0_3 * pc_yzz[j];

                t_z_zz[j] = fl_s_0_0_0 * (0.5 * pa_z[j] * fl1_fx + fl1_fx * pb_z[j] + pa_z[j] * pb_zz[j]);

                t_z_zz[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * fl1_fx - 1.5 * pc_z[j] * fl1_fx - fl1_fx * pb_z[j] - 2.0 * pa_z[j] * pb_z[j] * pc_z[j] - pc_z[j] * pb_zz[j]);

                t_z_zz[j] += fl_s_0_0_2 * (1.5 * pc_z[j] * fl1_fx + pa_z[j] * pc_zz[j] + 2.0 * pc_zz[j] * pb_z[j]);

                t_z_zz[j] += -fl_s_0_0_3 * pc_zzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForDP(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForDP_0_9(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForDP_9_18(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                    braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compNuclearPotentialForDP_0_9(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(9 * idx);

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            auto pa_xz = paDistances.data(9 * idx + 5);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(4 * idx);

            auto s_0_0_1 = auxBuffer.data(4 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(4 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(4 * idx + 3);

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

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xy, pa_xz, pa_y, pa_z, pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xxx, \
                                     pc_xxy, pc_xxz, pc_xy, pc_xyy, pc_xyz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yz, pc_z, pc_zz, \
                                     s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, t_xx_x, t_xx_y, t_xx_z, t_xy_x, t_xy_y, t_xy_z, \
                                     t_xz_x, t_xz_y, t_xz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl1_fx = fx[j];

                t_xx_x[j] = fl_s_0_0_0 * (pa_x[j] * fl1_fx + 0.5 * fl1_fx * pb_x[j] + pa_xx[j] * pb_x[j]);

                t_xx_x[j] += fl_s_0_0_1 * (-pa_x[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx - 0.5 * fl1_fx * pb_x[j] - pa_xx[j] * pc_x[j] - 2.0 * pa_x[j] * pc_x[j] * pb_x[j]);

                t_xx_x[j] += fl_s_0_0_2 * (1.5 * pc_x[j] * fl1_fx + 2.0 * pa_x[j] * pc_xx[j] + pc_xx[j] * pb_x[j]);

                t_xx_x[j] += -fl_s_0_0_3 * pc_xxx[j];

                t_xx_y[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_y[j] + pa_xx[j] * pb_y[j]);

                t_xx_y[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_y[j] - 0.5 * fl1_fx * pb_y[j] - pa_xx[j] * pc_y[j] - 2.0 * pa_x[j] * pc_x[j] * pb_y[j]);

                t_xx_y[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_y[j] + 2.0 * pa_x[j] * pc_xy[j] + pc_xx[j] * pb_y[j]);

                t_xx_y[j] += -fl_s_0_0_3 * pc_xxy[j];

                t_xx_z[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_z[j] + pa_xx[j] * pb_z[j]);

                t_xx_z[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pb_z[j] - pa_xx[j] * pc_z[j] - 2.0 * pa_x[j] * pc_x[j] * pb_z[j]);

                t_xx_z[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + 2.0 * pa_x[j] * pc_xz[j] + pc_xx[j] * pb_z[j]);

                t_xx_z[j] += -fl_s_0_0_3 * pc_xxz[j];

                t_xy_x[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_y[j] + pa_xy[j] * pb_x[j]);

                t_xy_x[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_y[j] - 0.5 * fl1_fx * pa_y[j] - pa_xy[j] * pc_x[j] - pa_x[j] * pc_y[j] * pb_x[j] - pc_x[j] * pa_y[j] * pb_x[j]);

                t_xy_x[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_y[j] + pa_x[j] * pc_xy[j] + pc_xx[j] * pa_y[j] + pc_xy[j] * pb_x[j]);

                t_xy_x[j] += -fl_s_0_0_3 * pc_xxy[j];

                t_xy_y[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx + pa_xy[j] * pb_y[j]);

                t_xy_y[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx - pa_xy[j] * pc_y[j] - pa_x[j] * pc_y[j] * pb_y[j] - pc_x[j] * pa_y[j] * pb_y[j]);

                t_xy_y[j] += fl_s_0_0_2 * (0.5 * pc_x[j] * fl1_fx + pa_x[j] * pc_yy[j] + pc_xy[j] * pa_y[j] + pc_xy[j] * pb_y[j]);

                t_xy_y[j] += -fl_s_0_0_3 * pc_xyy[j];

                t_xy_z[j] = fl_s_0_0_0 * pa_xy[j] * pb_z[j];

                t_xy_z[j] += fl_s_0_0_1 * (-pa_xy[j] * pc_z[j] - pa_x[j] * pc_y[j] * pb_z[j] - pc_x[j] * pa_y[j] * pb_z[j]);

                t_xy_z[j] += fl_s_0_0_2 * (pa_x[j] * pc_yz[j] + pc_xz[j] * pa_y[j] + pc_xy[j] * pb_z[j]);

                t_xy_z[j] += -fl_s_0_0_3 * pc_xyz[j];

                t_xz_x[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_z[j] + pa_xz[j] * pb_x[j]);

                t_xz_x[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pa_z[j] - pa_xz[j] * pc_x[j] - pa_x[j] * pc_z[j] * pb_x[j] - pc_x[j] * pa_z[j] * pb_x[j]);

                t_xz_x[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + pa_x[j] * pc_xz[j] + pc_xx[j] * pa_z[j] + pc_xz[j] * pb_x[j]);

                t_xz_x[j] += -fl_s_0_0_3 * pc_xxz[j];

                t_xz_y[j] = fl_s_0_0_0 * pa_xz[j] * pb_y[j];

                t_xz_y[j] += fl_s_0_0_1 * (-pa_xz[j] * pc_y[j] - pa_x[j] * pc_z[j] * pb_y[j] - pc_x[j] * pa_z[j] * pb_y[j]);

                t_xz_y[j] += fl_s_0_0_2 * (pa_x[j] * pc_yz[j] + pc_xy[j] * pa_z[j] + pc_xz[j] * pb_y[j]);

                t_xz_y[j] += -fl_s_0_0_3 * pc_xyz[j];

                t_xz_z[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx + pa_xz[j] * pb_z[j]);

                t_xz_z[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx - pa_xz[j] * pc_z[j] - pa_x[j] * pc_z[j] * pb_z[j] - pc_x[j] * pa_z[j] * pb_z[j]);

                t_xz_z[j] += fl_s_0_0_2 * (0.5 * pc_x[j] * fl1_fx + pa_x[j] * pc_zz[j] + pc_xz[j] * pa_z[j] + pc_xz[j] * pb_z[j]);

                t_xz_z[j] += -fl_s_0_0_3 * pc_xzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForDP_9_18(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(9 * idx + 6);

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(19 * idx);

            auto pc_y = pcDistances.data(19 * idx + 1);

            auto pc_z = pcDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xy = pcDistances.data(19 * idx + 4);

            auto pc_xz = pcDistances.data(19 * idx + 5);

            auto pc_yy = pcDistances.data(19 * idx + 6);

            auto pc_yz = pcDistances.data(19 * idx + 7);

            auto pc_zz = pcDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

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

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yz, pa_z, pa_zz, pb_x, pb_y, pb_z, pc_x, pc_xy, pc_xyy, pc_xyz, \
                                     pc_xz, pc_xzz, pc_y, pc_yy, pc_yyy, pc_yyz, pc_yz, pc_yzz, pc_z, pc_zz, pc_zzz, s_0_0_0, \
                                     s_0_0_1, s_0_0_2, s_0_0_3, t_yy_x, t_yy_y, t_yy_z, t_yz_x, t_yz_y, t_yz_z, t_zz_x, \
                                     t_zz_y, t_zz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl1_fx = fx[j];

                t_yy_x[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_x[j] + pa_yy[j] * pb_x[j]);

                t_yy_x[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_x[j] - 0.5 * fl1_fx * pb_x[j] - pa_yy[j] * pc_x[j] - 2.0 * pa_y[j] * pc_y[j] * pb_x[j]);

                t_yy_x[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_x[j] + 2.0 * pa_y[j] * pc_xy[j] + pc_yy[j] * pb_x[j]);

                t_yy_x[j] += -fl_s_0_0_3 * pc_xyy[j];

                t_yy_y[j] = fl_s_0_0_0 * (pa_y[j] * fl1_fx + 0.5 * fl1_fx * pb_y[j] + pa_yy[j] * pb_y[j]);

                t_yy_y[j] += fl_s_0_0_1 * (-pa_y[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx - 0.5 * fl1_fx * pb_y[j] - pa_yy[j] * pc_y[j] - 2.0 * pa_y[j] * pc_y[j] * pb_y[j]);

                t_yy_y[j] += fl_s_0_0_2 * (1.5 * pc_y[j] * fl1_fx + 2.0 * pa_y[j] * pc_yy[j] + pc_yy[j] * pb_y[j]);

                t_yy_y[j] += -fl_s_0_0_3 * pc_yyy[j];

                t_yy_z[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_z[j] + pa_yy[j] * pb_z[j]);

                t_yy_z[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pb_z[j] - pa_yy[j] * pc_z[j] - 2.0 * pa_y[j] * pc_y[j] * pb_z[j]);

                t_yy_z[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + 2.0 * pa_y[j] * pc_yz[j] + pc_yy[j] * pb_z[j]);

                t_yy_z[j] += -fl_s_0_0_3 * pc_yyz[j];

                t_yz_x[j] = fl_s_0_0_0 * pa_yz[j] * pb_x[j];

                t_yz_x[j] += fl_s_0_0_1 * (-pa_yz[j] * pc_x[j] - pa_y[j] * pc_z[j] * pb_x[j] - pc_y[j] * pa_z[j] * pb_x[j]);

                t_yz_x[j] += fl_s_0_0_2 * (pa_y[j] * pc_xz[j] + pc_xy[j] * pa_z[j] + pc_yz[j] * pb_x[j]);

                t_yz_x[j] += -fl_s_0_0_3 * pc_xyz[j];

                t_yz_y[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_z[j] + pa_yz[j] * pb_y[j]);

                t_yz_y[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_z[j] - 0.5 * fl1_fx * pa_z[j] - pa_yz[j] * pc_y[j] - pa_y[j] * pc_z[j] * pb_y[j] - pc_y[j] * pa_z[j] * pb_y[j]);

                t_yz_y[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_z[j] + pa_y[j] * pc_yz[j] + pc_yy[j] * pa_z[j] + pc_yz[j] * pb_y[j]);

                t_yz_y[j] += -fl_s_0_0_3 * pc_yyz[j];

                t_yz_z[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl1_fx + pa_yz[j] * pb_z[j]);

                t_yz_z[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl1_fx - 0.5 * pc_y[j] * fl1_fx - pa_yz[j] * pc_z[j] - pa_y[j] * pc_z[j] * pb_z[j] - pc_y[j] * pa_z[j] * pb_z[j]);

                t_yz_z[j] += fl_s_0_0_2 * (0.5 * pc_y[j] * fl1_fx + pa_y[j] * pc_zz[j] + pc_yz[j] * pa_z[j] + pc_yz[j] * pb_z[j]);

                t_yz_z[j] += -fl_s_0_0_3 * pc_yzz[j];

                t_zz_x[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_x[j] + pa_zz[j] * pb_x[j]);

                t_zz_x[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_x[j] - 0.5 * fl1_fx * pb_x[j] - pa_zz[j] * pc_x[j] - 2.0 * pa_z[j] * pc_z[j] * pb_x[j]);

                t_zz_x[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_x[j] + 2.0 * pa_z[j] * pc_xz[j] + pc_zz[j] * pb_x[j]);

                t_zz_x[j] += -fl_s_0_0_3 * pc_xzz[j];

                t_zz_y[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_y[j] + pa_zz[j] * pb_y[j]);

                t_zz_y[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pc_y[j] - 0.5 * fl1_fx * pb_y[j] - pa_zz[j] * pc_y[j] - 2.0 * pa_z[j] * pc_z[j] * pb_y[j]);

                t_zz_y[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_y[j] + 2.0 * pa_z[j] * pc_yz[j] + pc_zz[j] * pb_y[j]);

                t_zz_y[j] += -fl_s_0_0_3 * pc_yzz[j];

                t_zz_z[j] = fl_s_0_0_0 * (pa_z[j] * fl1_fx + 0.5 * fl1_fx * pb_z[j] + pa_zz[j] * pb_z[j]);

                t_zz_z[j] += fl_s_0_0_1 * (-pa_z[j] * fl1_fx - 1.5 * pc_z[j] * fl1_fx - 0.5 * fl1_fx * pb_z[j] - pa_zz[j] * pc_z[j] - 2.0 * pa_z[j] * pc_z[j] * pb_z[j]);

                t_zz_z[j] += fl_s_0_0_2 * (1.5 * pc_z[j] * fl1_fx + 2.0 * pa_z[j] * pc_zz[j] + pc_zz[j] * pb_z[j]);

                t_zz_z[j] += -fl_s_0_0_3 * pc_zzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPF(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForPF_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                    braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPF_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPF_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compNuclearPotentialForPF_0_10(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

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

            auto pc_xxx = pcDistances.data(34 * idx + 9);

            auto pc_xxy = pcDistances.data(34 * idx + 10);

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            auto pc_xzz = pcDistances.data(34 * idx + 14);

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(34 * idx + 19);

            auto pc_xxxy = pcDistances.data(34 * idx + 20);

            auto pc_xxxz = pcDistances.data(34 * idx + 21);

            auto pc_xxyy = pcDistances.data(34 * idx + 22);

            auto pc_xxyz = pcDistances.data(34 * idx + 23);

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

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa_x, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, pc_x, pc_xx, \
                                     pc_xxx, pc_xxxx, pc_xxxy, pc_xxxz, pc_xxy, pc_xxyy, pc_xxyz, pc_xxz, pc_xxzz, pc_xy, \
                                     pc_xyy, pc_xyyy, pc_xyyz, pc_xyz, pc_xyzz, pc_xz, pc_xzz, pc_xzzz, pc_y, pc_yy, pc_yyy, \
                                     pc_yyz, pc_yz, pc_yzz, pc_z, pc_zz, pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, \
                                     s_0_0_4, t_x_xxx, t_x_xxy, t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy, t_x_yyz, \
                                     t_x_yzz, t_x_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_x_xxx[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 1.5 * pa_x[j] * pb_x[j] * fl1_fx + 1.5 * fl1_fx * pb_xx[j] + pa_x[j] * pb_xxx[j]);

                t_x_xxx[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 1.5 * pa_x[j] * pb_x[j] * fl1_fx - 1.5 * pa_x[j] * pc_x[j] * fl1_fx - 4.5 * pc_x[j] * pb_x[j] * fl1_fx - 1.5 * fl1_fx * pb_xx[j] - 3.0 * pa_x[j] * pb_xx[j] * pc_x[j] - pc_x[j] * pb_xxx[j]);

                t_x_xxx[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 1.5 * pa_x[j] * pc_x[j] * fl1_fx + 4.5 * pc_x[j] * pb_x[j] * fl1_fx + 3.0 * pc_xx[j] * fl1_fx + 3.0 * pa_x[j] * pb_x[j] * pc_xx[j] + 3.0 * pc_xx[j] * pb_xx[j]);

                t_x_xxx[j] += fl_s_0_0_3 * (-3.0 * pc_xx[j] * fl1_fx - pa_x[j] * pc_xxx[j] - 3.0 * pc_xxx[j] * pb_x[j]);

                t_x_xxx[j] += fl_s_0_0_4 * pc_xxxx[j];

                t_x_xxy[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx * pb_y[j] + fl1_fx * pb_xy[j] + pa_x[j] * pb_xxy[j]);

                t_x_xxy[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx * pc_y[j] - 0.5 * pa_x[j] * fl1_fx * pb_y[j] - 1.5 * pc_x[j] * fl1_fx * pb_y[j] - fl1_fx * pb_x[j] * pc_y[j] - fl1_fx * pb_xy[j] - pa_x[j] * pb_xx[j] * pc_y[j] - 2.0 * pa_x[j] * pb_xy[j] * pc_x[j] - pc_x[j] * pb_xxy[j]);

                t_x_xxy[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl1_fx * pc_y[j] + 1.5 * pc_xy[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pb_y[j] + fl1_fx * pb_x[j] * pc_y[j] + 2.0 * pa_x[j] * pb_x[j] * pc_xy[j] + pa_x[j] * pc_xx[j] * pb_y[j] + pc_xy[j] * pb_xx[j] + 2.0 * pc_xx[j] * pb_xy[j]);

                t_x_xxy[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - pa_x[j] * pc_xxy[j] - 2.0 * pc_xxy[j] * pb_x[j] - pc_xxx[j] * pb_y[j]);

                t_x_xxy[j] += fl_s_0_0_4 * pc_xxxy[j];

                t_x_xxz[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx * pb_z[j] + fl1_fx * pb_xz[j] + pa_x[j] * pb_xxz[j]);

                t_x_xxz[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx * pc_z[j] - 0.5 * pa_x[j] * fl1_fx * pb_z[j] - 1.5 * pc_x[j] * fl1_fx * pb_z[j] - fl1_fx * pb_x[j] * pc_z[j] - fl1_fx * pb_xz[j] - pa_x[j] * pb_xx[j] * pc_z[j] - 2.0 * pa_x[j] * pb_xz[j] * pc_x[j] - pc_x[j] * pb_xxz[j]);

                t_x_xxz[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl1_fx * pc_z[j] + 1.5 * pc_xz[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pb_z[j] + fl1_fx * pb_x[j] * pc_z[j] + 2.0 * pa_x[j] * pb_x[j] * pc_xz[j] + pa_x[j] * pc_xx[j] * pb_z[j] + pc_xz[j] * pb_xx[j] + 2.0 * pc_xx[j] * pb_xz[j]);

                t_x_xxz[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - pa_x[j] * pc_xxz[j] - 2.0 * pc_xxz[j] * pb_x[j] - pc_xxx[j] * pb_z[j]);

                t_x_xxz[j] += fl_s_0_0_4 * pc_xxxz[j];

                t_x_xyy[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_x[j] * pb_x[j] * fl1_fx + 0.5 * fl1_fx * pb_yy[j] + pa_x[j] * pb_xyy[j]);

                t_x_xyy[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_x[j] * pb_x[j] * fl1_fx - 0.5 * pa_x[j] * pc_x[j] * fl1_fx - 0.5 * pc_x[j] * pb_x[j] * fl1_fx - fl1_fx * pb_y[j] * pc_y[j] - 0.5 * fl1_fx * pb_yy[j] - 2.0 * pa_x[j] * pb_xy[j] * pc_y[j] - pa_x[j] * pc_x[j] * pb_yy[j] - pc_x[j] * pb_xyy[j]);

                t_x_xyy[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * pa_x[j] * pc_x[j] * fl1_fx + 0.5 * pc_x[j] * pb_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx + 0.5 * fl1_fx * pc_yy[j] + fl1_fx * pb_y[j] * pc_y[j] + pa_x[j] * pb_x[j] * pc_yy[j] + 2.0 * pa_x[j] * pc_xy[j] * pb_y[j] + 2.0 * pc_xy[j] * pb_xy[j] + pc_xx[j] * pb_yy[j]);

                t_x_xyy[j] += fl_s_0_0_3 * (-0.5 * pc_xx[j] * fl1_fx - 0.5 * fl1_fx * pc_yy[j] - pa_x[j] * pc_xyy[j] - pc_xyy[j] * pb_x[j] - 2.0 * pc_xxy[j] * pb_y[j]);

                t_x_xyy[j] += fl_s_0_0_4 * pc_xxyy[j];

                t_x_xyz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_yz[j] + pa_x[j] * pb_xyz[j]);

                t_x_xyz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pb_y[j] * pc_z[j] - 0.5 * fl1_fx * pc_y[j] * pb_z[j] - 0.5 * fl1_fx * pb_yz[j] - pa_x[j] * pb_xy[j] * pc_z[j] - pa_x[j] * pb_xz[j] * pc_y[j] - pa_x[j] * pc_x[j] * pb_yz[j] - pc_x[j] * pb_xyz[j]);

                t_x_xyz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_yz[j] + 0.5 * fl1_fx * pb_y[j] * pc_z[j] + 0.5 * fl1_fx * pc_y[j] * pb_z[j] + pa_x[j] * pb_x[j] * pc_yz[j] + pa_x[j] * pc_xz[j] * pb_y[j] + pa_x[j] * pc_xy[j] * pb_z[j] + pc_xz[j] * pb_xy[j] + pc_xy[j] * pb_xz[j] + pc_xx[j] * pb_yz[j]);

                t_x_xyz[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_yz[j] - pa_x[j] * pc_xyz[j] - pc_xyz[j] * pb_x[j] - pc_xxz[j] * pb_y[j] - pc_xxy[j] * pb_z[j]);

                t_x_xyz[j] += fl_s_0_0_4 * pc_xxyz[j];

                t_x_xzz[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_x[j] * pb_x[j] * fl1_fx + 0.5 * fl1_fx * pb_zz[j] + pa_x[j] * pb_xzz[j]);

                t_x_xzz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_x[j] * pb_x[j] * fl1_fx - 0.5 * pa_x[j] * pc_x[j] * fl1_fx - 0.5 * pc_x[j] * pb_x[j] * fl1_fx - fl1_fx * pb_z[j] * pc_z[j] - 0.5 * fl1_fx * pb_zz[j] - 2.0 * pa_x[j] * pb_xz[j] * pc_z[j] - pa_x[j] * pc_x[j] * pb_zz[j] - pc_x[j] * pb_xzz[j]);

                t_x_xzz[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * pa_x[j] * pc_x[j] * fl1_fx + 0.5 * pc_x[j] * pb_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx + 0.5 * fl1_fx * pc_zz[j] + fl1_fx * pb_z[j] * pc_z[j] + pa_x[j] * pb_x[j] * pc_zz[j] + 2.0 * pa_x[j] * pc_xz[j] * pb_z[j] + 2.0 * pc_xz[j] * pb_xz[j] + pc_xx[j] * pb_zz[j]);

                t_x_xzz[j] += fl_s_0_0_3 * (-0.5 * pc_xx[j] * fl1_fx - 0.5 * fl1_fx * pc_zz[j] - pa_x[j] * pc_xzz[j] - pc_xzz[j] * pb_x[j] - 2.0 * pc_xxz[j] * pb_z[j]);

                t_x_xzz[j] += fl_s_0_0_4 * pc_xxzz[j];

                t_x_yyy[j] = fl_s_0_0_0 * (1.5 * pa_x[j] * pb_y[j] * fl1_fx + pa_x[j] * pb_yyy[j]);

                t_x_yyy[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * pb_y[j] * fl1_fx - 1.5 * pa_x[j] * pc_y[j] * fl1_fx - 1.5 * pc_x[j] * pb_y[j] * fl1_fx - 3.0 * pa_x[j] * pb_yy[j] * pc_y[j] - pc_x[j] * pb_yyy[j]);

                t_x_yyy[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * pc_y[j] * fl1_fx + 1.5 * pc_x[j] * pb_y[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx + 3.0 * pa_x[j] * pb_y[j] * pc_yy[j] + 3.0 * pc_xy[j] * pb_yy[j]);

                t_x_yyy[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - pa_x[j] * pc_yyy[j] - 3.0 * pc_xyy[j] * pb_y[j]);

                t_x_yyy[j] += fl_s_0_0_4 * pc_xyyy[j];

                t_x_yyz[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx * pb_z[j] + pa_x[j] * pb_yyz[j]);

                t_x_yyz[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx * pc_z[j] - 0.5 * pa_x[j] * fl1_fx * pb_z[j] - 0.5 * pc_x[j] * fl1_fx * pb_z[j] - pa_x[j] * pb_yy[j] * pc_z[j] - 2.0 * pa_x[j] * pb_yz[j] * pc_y[j] - pc_x[j] * pb_yyz[j]);

                t_x_yyz[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl1_fx * pc_z[j] + 0.5 * pc_xz[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pb_z[j] + 2.0 * pa_x[j] * pb_y[j] * pc_yz[j] + pa_x[j] * pc_yy[j] * pb_z[j] + pc_xz[j] * pb_yy[j] + 2.0 * pc_xy[j] * pb_yz[j]);

                t_x_yyz[j] += fl_s_0_0_3 * (-0.5 * pc_xz[j] * fl1_fx - pa_x[j] * pc_yyz[j] - 2.0 * pc_xyz[j] * pb_y[j] - pc_xyy[j] * pb_z[j]);

                t_x_yyz[j] += fl_s_0_0_4 * pc_xyyz[j];

                t_x_yzz[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * pb_y[j] * fl1_fx + pa_x[j] * pb_yzz[j]);

                t_x_yzz[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * pb_y[j] * fl1_fx - 0.5 * pa_x[j] * pc_y[j] * fl1_fx - 0.5 * pc_x[j] * pb_y[j] * fl1_fx - 2.0 * pa_x[j] * pb_yz[j] * pc_z[j] - pa_x[j] * pc_y[j] * pb_zz[j] - pc_x[j] * pb_yzz[j]);

                t_x_yzz[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * pc_y[j] * fl1_fx + 0.5 * pc_x[j] * pb_y[j] * fl1_fx + 0.5 * pc_xy[j] * fl1_fx + pa_x[j] * pb_y[j] * pc_zz[j] + 2.0 * pa_x[j] * pc_yz[j] * pb_z[j] + 2.0 * pc_xz[j] * pb_yz[j] + pc_xy[j] * pb_zz[j]);

                t_x_yzz[j] += fl_s_0_0_3 * (-0.5 * pc_xy[j] * fl1_fx - pa_x[j] * pc_yzz[j] - pc_xzz[j] * pb_y[j] - 2.0 * pc_xyz[j] * pb_z[j]);

                t_x_yzz[j] += fl_s_0_0_4 * pc_xyzz[j];

                t_x_zzz[j] = fl_s_0_0_0 * (1.5 * pa_x[j] * pb_z[j] * fl1_fx + pa_x[j] * pb_zzz[j]);

                t_x_zzz[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * pb_z[j] * fl1_fx - 1.5 * pa_x[j] * pc_z[j] * fl1_fx - 1.5 * pc_x[j] * pb_z[j] * fl1_fx - 3.0 * pa_x[j] * pb_zz[j] * pc_z[j] - pc_x[j] * pb_zzz[j]);

                t_x_zzz[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * pc_z[j] * fl1_fx + 1.5 * pc_x[j] * pb_z[j] * fl1_fx + 1.5 * pc_xz[j] * fl1_fx + 3.0 * pa_x[j] * pb_z[j] * pc_zz[j] + 3.0 * pc_xz[j] * pb_zz[j]);

                t_x_zzz[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - pa_x[j] * pc_zzz[j] - 3.0 * pc_xzz[j] * pb_z[j]);

                t_x_zzz[j] += fl_s_0_0_4 * pc_xzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPF_10_20(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(3 * idx + 1);

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

            auto pc_xxx = pcDistances.data(34 * idx + 9);

            auto pc_xxy = pcDistances.data(34 * idx + 10);

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            auto pc_xzz = pcDistances.data(34 * idx + 14);

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxy = pcDistances.data(34 * idx + 20);

            auto pc_xxyy = pcDistances.data(34 * idx + 22);

            auto pc_xxyz = pcDistances.data(34 * idx + 23);

            auto pc_xyyy = pcDistances.data(34 * idx + 25);

            auto pc_xyyz = pcDistances.data(34 * idx + 26);

            auto pc_xyzz = pcDistances.data(34 * idx + 27);

            auto pc_yyyy = pcDistances.data(34 * idx + 29);

            auto pc_yyyz = pcDistances.data(34 * idx + 30);

            auto pc_yyzz = pcDistances.data(34 * idx + 31);

            auto pc_yzzz = pcDistances.data(34 * idx + 32);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

            // set up pointers to integrals

            auto t_y_xxx = primBuffer.data(30 * idx + 10);

            auto t_y_xxy = primBuffer.data(30 * idx + 11);

            auto t_y_xxz = primBuffer.data(30 * idx + 12);

            auto t_y_xyy = primBuffer.data(30 * idx + 13);

            auto t_y_xyz = primBuffer.data(30 * idx + 14);

            auto t_y_xzz = primBuffer.data(30 * idx + 15);

            auto t_y_yyy = primBuffer.data(30 * idx + 16);

            auto t_y_yyz = primBuffer.data(30 * idx + 17);

            auto t_y_yzz = primBuffer.data(30 * idx + 18);

            auto t_y_zzz = primBuffer.data(30 * idx + 19);

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fx, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, pc_x, pc_xx, \
                                     pc_xxx, pc_xxxy, pc_xxy, pc_xxyy, pc_xxyz, pc_xxz, pc_xy, pc_xyy, pc_xyyy, pc_xyyz, \
                                     pc_xyz, pc_xyzz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yyy, pc_yyyy, pc_yyyz, pc_yyz, pc_yyzz, \
                                     pc_yz, pc_yzz, pc_yzzz, pc_z, pc_zz, pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, \
                                     s_0_0_4, t_y_xxx, t_y_xxy, t_y_xxz, t_y_xyy, t_y_xyz, t_y_xzz, t_y_yyy, t_y_yyz, \
                                     t_y_yzz, t_y_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_y_xxx[j] = fl_s_0_0_0 * (1.5 * pa_y[j] * pb_x[j] * fl1_fx + pa_y[j] * pb_xxx[j]);

                t_y_xxx[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * pb_x[j] * fl1_fx - 1.5 * pa_y[j] * pc_x[j] * fl1_fx - 1.5 * pc_y[j] * pb_x[j] * fl1_fx - 3.0 * pa_y[j] * pb_xx[j] * pc_x[j] - pc_y[j] * pb_xxx[j]);

                t_y_xxx[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * pc_x[j] * fl1_fx + 1.5 * pc_y[j] * pb_x[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx + 3.0 * pa_y[j] * pb_x[j] * pc_xx[j] + 3.0 * pc_xy[j] * pb_xx[j]);

                t_y_xxx[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - pa_y[j] * pc_xxx[j] - 3.0 * pc_xxy[j] * pb_x[j]);

                t_y_xxx[j] += fl_s_0_0_4 * pc_xxxy[j];

                t_y_xxy[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_y[j] * fl1_fx * pb_y[j] + 0.5 * fl1_fx * pb_xx[j] + pa_y[j] * pb_xxy[j]);

                t_y_xxy[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_y[j] * fl1_fx * pc_y[j] - 0.5 * pa_y[j] * fl1_fx * pb_y[j] - 0.5 * pc_y[j] * fl1_fx * pb_y[j] - fl1_fx * pb_x[j] * pc_x[j] - 0.5 * fl1_fx * pb_xx[j] - pa_y[j] * pb_xx[j] * pc_y[j] - 2.0 * pa_y[j] * pb_xy[j] * pc_x[j] - pc_y[j] * pb_xxy[j]);

                t_y_xxy[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * pa_y[j] * fl1_fx * pc_y[j] + 0.5 * pc_yy[j] * fl1_fx + 0.5 * pc_y[j] * fl1_fx * pb_y[j] + 0.5 * fl1_fx * pc_xx[j] + fl1_fx * pb_x[j] * pc_x[j] + 2.0 * pa_y[j] * pb_x[j] * pc_xy[j] + pa_y[j] * pc_xx[j] * pb_y[j] + pc_yy[j] * pb_xx[j] + 2.0 * pc_xy[j] * pb_xy[j]);

                t_y_xxy[j] += fl_s_0_0_3 * (-0.5 * pc_yy[j] * fl1_fx - 0.5 * fl1_fx * pc_xx[j] - pa_y[j] * pc_xxy[j] - 2.0 * pc_xyy[j] * pb_x[j] - pc_xxy[j] * pb_y[j]);

                t_y_xxy[j] += fl_s_0_0_4 * pc_xxyy[j];

                t_y_xxz[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl1_fx * pb_z[j] + pa_y[j] * pb_xxz[j]);

                t_y_xxz[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl1_fx * pc_z[j] - 0.5 * pa_y[j] * fl1_fx * pb_z[j] - 0.5 * pc_y[j] * fl1_fx * pb_z[j] - pa_y[j] * pb_xx[j] * pc_z[j] - 2.0 * pa_y[j] * pb_xz[j] * pc_x[j] - pc_y[j] * pb_xxz[j]);

                t_y_xxz[j] += fl_s_0_0_2 * (0.5 * pa_y[j] * fl1_fx * pc_z[j] + 0.5 * pc_yz[j] * fl1_fx + 0.5 * pc_y[j] * fl1_fx * pb_z[j] + 2.0 * pa_y[j] * pb_x[j] * pc_xz[j] + pa_y[j] * pc_xx[j] * pb_z[j] + pc_yz[j] * pb_xx[j] + 2.0 * pc_xy[j] * pb_xz[j]);

                t_y_xxz[j] += fl_s_0_0_3 * (-0.5 * pc_yz[j] * fl1_fx - pa_y[j] * pc_xxz[j] - 2.0 * pc_xyz[j] * pb_x[j] - pc_xxy[j] * pb_z[j]);

                t_y_xxz[j] += fl_s_0_0_4 * pc_xxyz[j];

                t_y_xyy[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * pb_x[j] * fl1_fx + fl1_fx * pb_xy[j] + pa_y[j] * pb_xyy[j]);

                t_y_xyy[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * pb_x[j] * fl1_fx - 0.5 * pa_y[j] * pc_x[j] * fl1_fx - 1.5 * pc_y[j] * pb_x[j] * fl1_fx - fl1_fx * pc_x[j] * pb_y[j] - fl1_fx * pb_xy[j] - 2.0 * pa_y[j] * pb_xy[j] * pc_y[j] - pa_y[j] * pc_x[j] * pb_yy[j] - pc_y[j] * pb_xyy[j]);

                t_y_xyy[j] += fl_s_0_0_2 * (0.5 * pa_y[j] * pc_x[j] * fl1_fx + 1.5 * pc_y[j] * pb_x[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx + fl1_fx * pc_x[j] * pb_y[j] + pa_y[j] * pb_x[j] * pc_yy[j] + 2.0 * pa_y[j] * pc_xy[j] * pb_y[j] + 2.0 * pc_yy[j] * pb_xy[j] + pc_xy[j] * pb_yy[j]);

                t_y_xyy[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - pa_y[j] * pc_xyy[j] - pc_yyy[j] * pb_x[j] - 2.0 * pc_xyy[j] * pb_y[j]);

                t_y_xyy[j] += fl_s_0_0_4 * pc_xyyy[j];

                t_y_xyz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_xz[j] + pa_y[j] * pb_xyz[j]);

                t_y_xyz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pb_x[j] * pc_z[j] - 0.5 * fl1_fx * pc_x[j] * pb_z[j] - 0.5 * fl1_fx * pb_xz[j] - pa_y[j] * pb_xy[j] * pc_z[j] - pa_y[j] * pb_xz[j] * pc_y[j] - pa_y[j] * pc_x[j] * pb_yz[j] - pc_y[j] * pb_xyz[j]);

                t_y_xyz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_xz[j] + 0.5 * fl1_fx * pb_x[j] * pc_z[j] + 0.5 * fl1_fx * pc_x[j] * pb_z[j] + pa_y[j] * pb_x[j] * pc_yz[j] + pa_y[j] * pc_xz[j] * pb_y[j] + pa_y[j] * pc_xy[j] * pb_z[j] + pc_yz[j] * pb_xy[j] + pc_yy[j] * pb_xz[j] + pc_xy[j] * pb_yz[j]);

                t_y_xyz[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_xz[j] - pa_y[j] * pc_xyz[j] - pc_yyz[j] * pb_x[j] - pc_xyz[j] * pb_y[j] - pc_xyy[j] * pb_z[j]);

                t_y_xyz[j] += fl_s_0_0_4 * pc_xyyz[j];

                t_y_xzz[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * pb_x[j] * fl1_fx + pa_y[j] * pb_xzz[j]);

                t_y_xzz[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * pb_x[j] * fl1_fx - 0.5 * pa_y[j] * pc_x[j] * fl1_fx - 0.5 * pc_y[j] * pb_x[j] * fl1_fx - 2.0 * pa_y[j] * pb_xz[j] * pc_z[j] - pa_y[j] * pc_x[j] * pb_zz[j] - pc_y[j] * pb_xzz[j]);

                t_y_xzz[j] += fl_s_0_0_2 * (0.5 * pa_y[j] * pc_x[j] * fl1_fx + 0.5 * pc_y[j] * pb_x[j] * fl1_fx + 0.5 * pc_xy[j] * fl1_fx + pa_y[j] * pb_x[j] * pc_zz[j] + 2.0 * pa_y[j] * pc_xz[j] * pb_z[j] + 2.0 * pc_yz[j] * pb_xz[j] + pc_xy[j] * pb_zz[j]);

                t_y_xzz[j] += fl_s_0_0_3 * (-0.5 * pc_xy[j] * fl1_fx - pa_y[j] * pc_xzz[j] - pc_yzz[j] * pb_x[j] - 2.0 * pc_xyz[j] * pb_z[j]);

                t_y_xzz[j] += fl_s_0_0_4 * pc_xyzz[j];

                t_y_yyy[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 1.5 * pa_y[j] * pb_y[j] * fl1_fx + 1.5 * fl1_fx * pb_yy[j] + pa_y[j] * pb_yyy[j]);

                t_y_yyy[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 1.5 * pa_y[j] * pb_y[j] * fl1_fx - 1.5 * pa_y[j] * pc_y[j] * fl1_fx - 4.5 * pc_y[j] * pb_y[j] * fl1_fx - 1.5 * fl1_fx * pb_yy[j] - 3.0 * pa_y[j] * pb_yy[j] * pc_y[j] - pc_y[j] * pb_yyy[j]);

                t_y_yyy[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 1.5 * pa_y[j] * pc_y[j] * fl1_fx + 4.5 * pc_y[j] * pb_y[j] * fl1_fx + 3.0 * pc_yy[j] * fl1_fx + 3.0 * pa_y[j] * pb_y[j] * pc_yy[j] + 3.0 * pc_yy[j] * pb_yy[j]);

                t_y_yyy[j] += fl_s_0_0_3 * (-3.0 * pc_yy[j] * fl1_fx - pa_y[j] * pc_yyy[j] - 3.0 * pc_yyy[j] * pb_y[j]);

                t_y_yyy[j] += fl_s_0_0_4 * pc_yyyy[j];

                t_y_yyz[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl1_fx * pb_z[j] + fl1_fx * pb_yz[j] + pa_y[j] * pb_yyz[j]);

                t_y_yyz[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl1_fx * pc_z[j] - 0.5 * pa_y[j] * fl1_fx * pb_z[j] - 1.5 * pc_y[j] * fl1_fx * pb_z[j] - fl1_fx * pb_y[j] * pc_z[j] - fl1_fx * pb_yz[j] - pa_y[j] * pb_yy[j] * pc_z[j] - 2.0 * pa_y[j] * pb_yz[j] * pc_y[j] - pc_y[j] * pb_yyz[j]);

                t_y_yyz[j] += fl_s_0_0_2 * (0.5 * pa_y[j] * fl1_fx * pc_z[j] + 1.5 * pc_yz[j] * fl1_fx + 1.5 * pc_y[j] * fl1_fx * pb_z[j] + fl1_fx * pb_y[j] * pc_z[j] + 2.0 * pa_y[j] * pb_y[j] * pc_yz[j] + pa_y[j] * pc_yy[j] * pb_z[j] + pc_yz[j] * pb_yy[j] + 2.0 * pc_yy[j] * pb_yz[j]);

                t_y_yyz[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - pa_y[j] * pc_yyz[j] - 2.0 * pc_yyz[j] * pb_y[j] - pc_yyy[j] * pb_z[j]);

                t_y_yyz[j] += fl_s_0_0_4 * pc_yyyz[j];

                t_y_yzz[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_y[j] * pb_y[j] * fl1_fx + 0.5 * fl1_fx * pb_zz[j] + pa_y[j] * pb_yzz[j]);

                t_y_yzz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_y[j] * pb_y[j] * fl1_fx - 0.5 * pa_y[j] * pc_y[j] * fl1_fx - 0.5 * pc_y[j] * pb_y[j] * fl1_fx - fl1_fx * pb_z[j] * pc_z[j] - 0.5 * fl1_fx * pb_zz[j] - 2.0 * pa_y[j] * pb_yz[j] * pc_z[j] - pa_y[j] * pc_y[j] * pb_zz[j] - pc_y[j] * pb_yzz[j]);

                t_y_yzz[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * pa_y[j] * pc_y[j] * fl1_fx + 0.5 * pc_y[j] * pb_y[j] * fl1_fx + 0.5 * pc_yy[j] * fl1_fx + 0.5 * fl1_fx * pc_zz[j] + fl1_fx * pb_z[j] * pc_z[j] + pa_y[j] * pb_y[j] * pc_zz[j] + 2.0 * pa_y[j] * pc_yz[j] * pb_z[j] + 2.0 * pc_yz[j] * pb_yz[j] + pc_yy[j] * pb_zz[j]);

                t_y_yzz[j] += fl_s_0_0_3 * (-0.5 * pc_yy[j] * fl1_fx - 0.5 * fl1_fx * pc_zz[j] - pa_y[j] * pc_yzz[j] - pc_yzz[j] * pb_y[j] - 2.0 * pc_yyz[j] * pb_z[j]);

                t_y_yzz[j] += fl_s_0_0_4 * pc_yyzz[j];

                t_y_zzz[j] = fl_s_0_0_0 * (1.5 * pa_y[j] * pb_z[j] * fl1_fx + pa_y[j] * pb_zzz[j]);

                t_y_zzz[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * pb_z[j] * fl1_fx - 1.5 * pa_y[j] * pc_z[j] * fl1_fx - 1.5 * pc_y[j] * pb_z[j] * fl1_fx - 3.0 * pa_y[j] * pb_zz[j] * pc_z[j] - pc_y[j] * pb_zzz[j]);

                t_y_zzz[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * pc_z[j] * fl1_fx + 1.5 * pc_y[j] * pb_z[j] * fl1_fx + 1.5 * pc_yz[j] * fl1_fx + 3.0 * pa_y[j] * pb_z[j] * pc_zz[j] + 3.0 * pc_yz[j] * pb_zz[j]);

                t_y_zzz[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - pa_y[j] * pc_zzz[j] - 3.0 * pc_yzz[j] * pb_z[j]);

                t_y_zzz[j] += fl_s_0_0_4 * pc_yzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPF_20_30(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_z = paDistances.data(3 * idx + 2);

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

            auto pc_xxx = pcDistances.data(34 * idx + 9);

            auto pc_xxy = pcDistances.data(34 * idx + 10);

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            auto pc_xzz = pcDistances.data(34 * idx + 14);

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxz = pcDistances.data(34 * idx + 21);

            auto pc_xxyz = pcDistances.data(34 * idx + 23);

            auto pc_xxzz = pcDistances.data(34 * idx + 24);

            auto pc_xyyz = pcDistances.data(34 * idx + 26);

            auto pc_xyzz = pcDistances.data(34 * idx + 27);

            auto pc_xzzz = pcDistances.data(34 * idx + 28);

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

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fx, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, pc_x, pc_xx, \
                                     pc_xxx, pc_xxxz, pc_xxy, pc_xxyz, pc_xxz, pc_xxzz, pc_xy, pc_xyy, pc_xyyz, pc_xyz, \
                                     pc_xyzz, pc_xz, pc_xzz, pc_xzzz, pc_y, pc_yy, pc_yyy, pc_yyyz, pc_yyz, pc_yyzz, pc_yz, \
                                     pc_yzz, pc_yzzz, pc_z, pc_zz, pc_zzz, pc_zzzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, \
                                     s_0_0_4, t_z_xxx, t_z_xxy, t_z_xxz, t_z_xyy, t_z_xyz, t_z_xzz, t_z_yyy, t_z_yyz, \
                                     t_z_yzz, t_z_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_z_xxx[j] = fl_s_0_0_0 * (1.5 * pa_z[j] * pb_x[j] * fl1_fx + pa_z[j] * pb_xxx[j]);

                t_z_xxx[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * pb_x[j] * fl1_fx - 1.5 * pa_z[j] * pc_x[j] * fl1_fx - 1.5 * pc_z[j] * pb_x[j] * fl1_fx - 3.0 * pa_z[j] * pb_xx[j] * pc_x[j] - pc_z[j] * pb_xxx[j]);

                t_z_xxx[j] += fl_s_0_0_2 * (1.5 * pa_z[j] * pc_x[j] * fl1_fx + 1.5 * pc_z[j] * pb_x[j] * fl1_fx + 1.5 * pc_xz[j] * fl1_fx + 3.0 * pa_z[j] * pb_x[j] * pc_xx[j] + 3.0 * pc_xz[j] * pb_xx[j]);

                t_z_xxx[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - pa_z[j] * pc_xxx[j] - 3.0 * pc_xxz[j] * pb_x[j]);

                t_z_xxx[j] += fl_s_0_0_4 * pc_xxxz[j];

                t_z_xxy[j] = fl_s_0_0_0 * (0.5 * pa_z[j] * fl1_fx * pb_y[j] + pa_z[j] * pb_xxy[j]);

                t_z_xxy[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * fl1_fx * pc_y[j] - 0.5 * pa_z[j] * fl1_fx * pb_y[j] - 0.5 * pc_z[j] * fl1_fx * pb_y[j] - pa_z[j] * pb_xx[j] * pc_y[j] - 2.0 * pa_z[j] * pb_xy[j] * pc_x[j] - pc_z[j] * pb_xxy[j]);

                t_z_xxy[j] += fl_s_0_0_2 * (0.5 * pa_z[j] * fl1_fx * pc_y[j] + 0.5 * pc_yz[j] * fl1_fx + 0.5 * pc_z[j] * fl1_fx * pb_y[j] + 2.0 * pa_z[j] * pb_x[j] * pc_xy[j] + pa_z[j] * pc_xx[j] * pb_y[j] + pc_yz[j] * pb_xx[j] + 2.0 * pc_xz[j] * pb_xy[j]);

                t_z_xxy[j] += fl_s_0_0_3 * (-0.5 * pc_yz[j] * fl1_fx - pa_z[j] * pc_xxy[j] - 2.0 * pc_xyz[j] * pb_x[j] - pc_xxz[j] * pb_y[j]);

                t_z_xxy[j] += fl_s_0_0_4 * pc_xxyz[j];

                t_z_xxz[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_z[j] * fl1_fx * pb_z[j] + 0.5 * fl1_fx * pb_xx[j] + pa_z[j] * pb_xxz[j]);

                t_z_xxz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_z[j] * fl1_fx * pc_z[j] - 0.5 * pa_z[j] * fl1_fx * pb_z[j] - 0.5 * pc_z[j] * fl1_fx * pb_z[j] - fl1_fx * pb_x[j] * pc_x[j] - 0.5 * fl1_fx * pb_xx[j] - pa_z[j] * pb_xx[j] * pc_z[j] - 2.0 * pa_z[j] * pb_xz[j] * pc_x[j] - pc_z[j] * pb_xxz[j]);

                t_z_xxz[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * pa_z[j] * fl1_fx * pc_z[j] + 0.5 * pc_zz[j] * fl1_fx + 0.5 * pc_z[j] * fl1_fx * pb_z[j] + 0.5 * fl1_fx * pc_xx[j] + fl1_fx * pb_x[j] * pc_x[j] + 2.0 * pa_z[j] * pb_x[j] * pc_xz[j] + pa_z[j] * pc_xx[j] * pb_z[j] + pc_zz[j] * pb_xx[j] + 2.0 * pc_xz[j] * pb_xz[j]);

                t_z_xxz[j] += fl_s_0_0_3 * (-0.5 * pc_zz[j] * fl1_fx - 0.5 * fl1_fx * pc_xx[j] - pa_z[j] * pc_xxz[j] - 2.0 * pc_xzz[j] * pb_x[j] - pc_xxz[j] * pb_z[j]);

                t_z_xxz[j] += fl_s_0_0_4 * pc_xxzz[j];

                t_z_xyy[j] = fl_s_0_0_0 * (0.5 * pa_z[j] * pb_x[j] * fl1_fx + pa_z[j] * pb_xyy[j]);

                t_z_xyy[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * pb_x[j] * fl1_fx - 0.5 * pa_z[j] * pc_x[j] * fl1_fx - 0.5 * pc_z[j] * pb_x[j] * fl1_fx - 2.0 * pa_z[j] * pb_xy[j] * pc_y[j] - pa_z[j] * pc_x[j] * pb_yy[j] - pc_z[j] * pb_xyy[j]);

                t_z_xyy[j] += fl_s_0_0_2 * (0.5 * pa_z[j] * pc_x[j] * fl1_fx + 0.5 * pc_z[j] * pb_x[j] * fl1_fx + 0.5 * pc_xz[j] * fl1_fx + pa_z[j] * pb_x[j] * pc_yy[j] + 2.0 * pa_z[j] * pc_xy[j] * pb_y[j] + 2.0 * pc_yz[j] * pb_xy[j] + pc_xz[j] * pb_yy[j]);

                t_z_xyy[j] += fl_s_0_0_3 * (-0.5 * pc_xz[j] * fl1_fx - pa_z[j] * pc_xyy[j] - pc_yyz[j] * pb_x[j] - 2.0 * pc_xyz[j] * pb_y[j]);

                t_z_xyy[j] += fl_s_0_0_4 * pc_xyyz[j];

                t_z_xyz[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pb_xy[j] + pa_z[j] * pb_xyz[j]);

                t_z_xyz[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pb_x[j] * pc_y[j] - 0.5 * fl1_fx * pc_x[j] * pb_y[j] - 0.5 * fl1_fx * pb_xy[j] - pa_z[j] * pb_xy[j] * pc_z[j] - pa_z[j] * pb_xz[j] * pc_y[j] - pa_z[j] * pc_x[j] * pb_yz[j] - pc_z[j] * pb_xyz[j]);

                t_z_xyz[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_xy[j] + 0.5 * fl1_fx * pb_x[j] * pc_y[j] + 0.5 * fl1_fx * pc_x[j] * pb_y[j] + pa_z[j] * pb_x[j] * pc_yz[j] + pa_z[j] * pc_xz[j] * pb_y[j] + pa_z[j] * pc_xy[j] * pb_z[j] + pc_zz[j] * pb_xy[j] + pc_yz[j] * pb_xz[j] + pc_xz[j] * pb_yz[j]);

                t_z_xyz[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_xy[j] - pa_z[j] * pc_xyz[j] - pc_yzz[j] * pb_x[j] - pc_xzz[j] * pb_y[j] - pc_xyz[j] * pb_z[j]);

                t_z_xyz[j] += fl_s_0_0_4 * pc_xyzz[j];

                t_z_xzz[j] = fl_s_0_0_0 * (0.5 * pa_z[j] * pb_x[j] * fl1_fx + fl1_fx * pb_xz[j] + pa_z[j] * pb_xzz[j]);

                t_z_xzz[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * pb_x[j] * fl1_fx - 0.5 * pa_z[j] * pc_x[j] * fl1_fx - 1.5 * pc_z[j] * pb_x[j] * fl1_fx - fl1_fx * pc_x[j] * pb_z[j] - fl1_fx * pb_xz[j] - 2.0 * pa_z[j] * pb_xz[j] * pc_z[j] - pa_z[j] * pc_x[j] * pb_zz[j] - pc_z[j] * pb_xzz[j]);

                t_z_xzz[j] += fl_s_0_0_2 * (0.5 * pa_z[j] * pc_x[j] * fl1_fx + 1.5 * pc_z[j] * pb_x[j] * fl1_fx + 1.5 * pc_xz[j] * fl1_fx + fl1_fx * pc_x[j] * pb_z[j] + pa_z[j] * pb_x[j] * pc_zz[j] + 2.0 * pa_z[j] * pc_xz[j] * pb_z[j] + 2.0 * pc_zz[j] * pb_xz[j] + pc_xz[j] * pb_zz[j]);

                t_z_xzz[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - pa_z[j] * pc_xzz[j] - pc_zzz[j] * pb_x[j] - 2.0 * pc_xzz[j] * pb_z[j]);

                t_z_xzz[j] += fl_s_0_0_4 * pc_xzzz[j];

                t_z_yyy[j] = fl_s_0_0_0 * (1.5 * pa_z[j] * pb_y[j] * fl1_fx + pa_z[j] * pb_yyy[j]);

                t_z_yyy[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * pb_y[j] * fl1_fx - 1.5 * pa_z[j] * pc_y[j] * fl1_fx - 1.5 * pc_z[j] * pb_y[j] * fl1_fx - 3.0 * pa_z[j] * pb_yy[j] * pc_y[j] - pc_z[j] * pb_yyy[j]);

                t_z_yyy[j] += fl_s_0_0_2 * (1.5 * pa_z[j] * pc_y[j] * fl1_fx + 1.5 * pc_z[j] * pb_y[j] * fl1_fx + 1.5 * pc_yz[j] * fl1_fx + 3.0 * pa_z[j] * pb_y[j] * pc_yy[j] + 3.0 * pc_yz[j] * pb_yy[j]);

                t_z_yyy[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - pa_z[j] * pc_yyy[j] - 3.0 * pc_yyz[j] * pb_y[j]);

                t_z_yyy[j] += fl_s_0_0_4 * pc_yyyz[j];

                t_z_yyz[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_z[j] * fl1_fx * pb_z[j] + 0.5 * fl1_fx * pb_yy[j] + pa_z[j] * pb_yyz[j]);

                t_z_yyz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_z[j] * fl1_fx * pc_z[j] - 0.5 * pa_z[j] * fl1_fx * pb_z[j] - 0.5 * pc_z[j] * fl1_fx * pb_z[j] - fl1_fx * pb_y[j] * pc_y[j] - 0.5 * fl1_fx * pb_yy[j] - pa_z[j] * pb_yy[j] * pc_z[j] - 2.0 * pa_z[j] * pb_yz[j] * pc_y[j] - pc_z[j] * pb_yyz[j]);

                t_z_yyz[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * pa_z[j] * fl1_fx * pc_z[j] + 0.5 * pc_zz[j] * fl1_fx + 0.5 * pc_z[j] * fl1_fx * pb_z[j] + 0.5 * fl1_fx * pc_yy[j] + fl1_fx * pb_y[j] * pc_y[j] + 2.0 * pa_z[j] * pb_y[j] * pc_yz[j] + pa_z[j] * pc_yy[j] * pb_z[j] + pc_zz[j] * pb_yy[j] + 2.0 * pc_yz[j] * pb_yz[j]);

                t_z_yyz[j] += fl_s_0_0_3 * (-0.5 * pc_zz[j] * fl1_fx - 0.5 * fl1_fx * pc_yy[j] - pa_z[j] * pc_yyz[j] - 2.0 * pc_yzz[j] * pb_y[j] - pc_yyz[j] * pb_z[j]);

                t_z_yyz[j] += fl_s_0_0_4 * pc_yyzz[j];

                t_z_yzz[j] = fl_s_0_0_0 * (0.5 * pa_z[j] * pb_y[j] * fl1_fx + fl1_fx * pb_yz[j] + pa_z[j] * pb_yzz[j]);

                t_z_yzz[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * pb_y[j] * fl1_fx - 0.5 * pa_z[j] * pc_y[j] * fl1_fx - 1.5 * pc_z[j] * pb_y[j] * fl1_fx - fl1_fx * pc_y[j] * pb_z[j] - fl1_fx * pb_yz[j] - 2.0 * pa_z[j] * pb_yz[j] * pc_z[j] - pa_z[j] * pc_y[j] * pb_zz[j] - pc_z[j] * pb_yzz[j]);

                t_z_yzz[j] += fl_s_0_0_2 * (0.5 * pa_z[j] * pc_y[j] * fl1_fx + 1.5 * pc_z[j] * pb_y[j] * fl1_fx + 1.5 * pc_yz[j] * fl1_fx + fl1_fx * pc_y[j] * pb_z[j] + pa_z[j] * pb_y[j] * pc_zz[j] + 2.0 * pa_z[j] * pc_yz[j] * pb_z[j] + 2.0 * pc_zz[j] * pb_yz[j] + pc_yz[j] * pb_zz[j]);

                t_z_yzz[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - pa_z[j] * pc_yzz[j] - pc_zzz[j] * pb_y[j] - 2.0 * pc_yzz[j] * pb_z[j]);

                t_z_yzz[j] += fl_s_0_0_4 * pc_yzzz[j];

                t_z_zzz[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 1.5 * pa_z[j] * pb_z[j] * fl1_fx + 1.5 * fl1_fx * pb_zz[j] + pa_z[j] * pb_zzz[j]);

                t_z_zzz[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 1.5 * pa_z[j] * pb_z[j] * fl1_fx - 1.5 * pa_z[j] * pc_z[j] * fl1_fx - 4.5 * pc_z[j] * pb_z[j] * fl1_fx - 1.5 * fl1_fx * pb_zz[j] - 3.0 * pa_z[j] * pb_zz[j] * pc_z[j] - pc_z[j] * pb_zzz[j]);

                t_z_zzz[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 1.5 * pa_z[j] * pc_z[j] * fl1_fx + 4.5 * pc_z[j] * pb_z[j] * fl1_fx + 3.0 * pc_zz[j] * fl1_fx + 3.0 * pa_z[j] * pb_z[j] * pc_zz[j] + 3.0 * pc_zz[j] * pb_zz[j]);

                t_z_zzz[j] += fl_s_0_0_3 * (-3.0 * pc_zz[j] * fl1_fx - pa_z[j] * pc_zzz[j] - 3.0 * pc_zzz[j] * pb_z[j]);

                t_z_zzz[j] += fl_s_0_0_4 * pc_zzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForFP(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForFP_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                    braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForFP_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForFP_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compNuclearPotentialForFP_0_10(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            auto pc_xxx = pcDistances.data(34 * idx + 9);

            auto pc_xxy = pcDistances.data(34 * idx + 10);

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            auto pc_xzz = pcDistances.data(34 * idx + 14);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(34 * idx + 19);

            auto pc_xxxy = pcDistances.data(34 * idx + 20);

            auto pc_xxxz = pcDistances.data(34 * idx + 21);

            auto pc_xxyy = pcDistances.data(34 * idx + 22);

            auto pc_xxyz = pcDistances.data(34 * idx + 23);

            auto pc_xxzz = pcDistances.data(34 * idx + 24);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

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

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxy, pa_xxz, pa_xy, pa_xyy, pa_xz, pa_y, pa_yy, pa_z, \
                                     pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xxx, pc_xxxx, pc_xxxy, pc_xxxz, pc_xxy, pc_xxyy, \
                                     pc_xxyz, pc_xxz, pc_xxzz, pc_xy, pc_xyy, pc_xyz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yz, pc_z, \
                                     pc_zz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, t_xxx_x, t_xxx_y, t_xxx_z, \
                                     t_xxy_x, t_xxy_y, t_xxy_z, t_xxz_x, t_xxz_y, t_xxz_z, t_xyy_x: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxx_x[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 1.5 * pa_xx[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pb_x[j] + pa_xxx[j] * pb_x[j]);

                t_xxx_x[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 1.5 * pa_xx[j] * fl1_fx - 4.5 * pa_x[j] * pc_x[j] * fl1_fx - 1.5 * pa_x[j] * fl1_fx * pb_x[j] - 1.5 * pc_x[j] * fl1_fx * pb_x[j] - pa_xxx[j] * pc_x[j] - 3.0 * pa_xx[j] * pc_x[j] * pb_x[j]);

                t_xxx_x[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 4.5 * pa_x[j] * pc_x[j] * fl1_fx + 3.0 * pc_xx[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pb_x[j] + 3.0 * pa_xx[j] * pc_xx[j] + 3.0 * pa_x[j] * pc_xx[j] * pb_x[j]);

                t_xxx_x[j] += fl_s_0_0_3 * (-3.0 * pc_xx[j] * fl1_fx - 3.0 * pa_x[j] * pc_xxx[j] - pc_xxx[j] * pb_x[j]);

                t_xxx_x[j] += fl_s_0_0_4 * pc_xxxx[j];

                t_xxx_y[j] = fl_s_0_0_0 * (1.5 * pa_x[j] * fl1_fx * pb_y[j] + pa_xxx[j] * pb_y[j]);

                t_xxx_y[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl1_fx * pc_y[j] - 1.5 * pa_x[j] * fl1_fx * pb_y[j] - 1.5 * pc_x[j] * fl1_fx * pb_y[j] - pa_xxx[j] * pc_y[j] - 3.0 * pa_xx[j] * pc_x[j] * pb_y[j]);

                t_xxx_y[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * fl1_fx * pc_y[j] + 1.5 * pc_xy[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pb_y[j] + 3.0 * pa_xx[j] * pc_xy[j] + 3.0 * pa_x[j] * pc_xx[j] * pb_y[j]);

                t_xxx_y[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - 3.0 * pa_x[j] * pc_xxy[j] - pc_xxx[j] * pb_y[j]);

                t_xxx_y[j] += fl_s_0_0_4 * pc_xxxy[j];

                t_xxx_z[j] = fl_s_0_0_0 * (1.5 * pa_x[j] * fl1_fx * pb_z[j] + pa_xxx[j] * pb_z[j]);

                t_xxx_z[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl1_fx * pc_z[j] - 1.5 * pa_x[j] * fl1_fx * pb_z[j] - 1.5 * pc_x[j] * fl1_fx * pb_z[j] - pa_xxx[j] * pc_z[j] - 3.0 * pa_xx[j] * pc_x[j] * pb_z[j]);

                t_xxx_z[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * fl1_fx * pc_z[j] + 1.5 * pc_xz[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pb_z[j] + 3.0 * pa_xx[j] * pc_xz[j] + 3.0 * pa_x[j] * pc_xx[j] * pb_z[j]);

                t_xxx_z[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - 3.0 * pa_x[j] * pc_xxz[j] - pc_xxx[j] * pb_z[j]);

                t_xxx_z[j] += fl_s_0_0_4 * pc_xxxz[j];

                t_xxy_x[j] = fl_s_0_0_0 * (pa_xy[j] * fl1_fx + 0.5 * fl1_fx * pa_y[j] * pb_x[j] + pa_xxy[j] * pb_x[j]);

                t_xxy_x[j] += fl_s_0_0_1 * (-pa_x[j] * fl1_fx * pc_y[j] - pa_xy[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pa_y[j] - 0.5 * fl1_fx * pc_y[j] * pb_x[j] - 0.5 * fl1_fx * pa_y[j] * pb_x[j] - pa_xxy[j] * pc_x[j] - pa_xx[j] * pc_y[j] * pb_x[j] - 2.0 * pa_xy[j] * pc_x[j] * pb_x[j]);

                t_xxy_x[j] += fl_s_0_0_2 * (pa_x[j] * fl1_fx * pc_y[j] + 1.5 * pc_xy[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pa_y[j] + 0.5 * fl1_fx * pc_y[j] * pb_x[j] + pa_xx[j] * pc_xy[j] + 2.0 * pa_xy[j] * pc_xx[j] + 2.0 * pa_x[j] * pc_xy[j] * pb_x[j] + pc_xx[j] * pa_y[j] * pb_x[j]);

                t_xxy_x[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - 2.0 * pa_x[j] * pc_xxy[j] - pc_xxx[j] * pa_y[j] - pc_xxy[j] * pb_x[j]);

                t_xxy_x[j] += fl_s_0_0_4 * pc_xxxy[j];

                t_xxy_y[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * fl1_fx * pa_y[j] * pb_y[j] + pa_xxy[j] * pb_y[j]);

                t_xxy_y[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_xx[j] * fl1_fx - pa_x[j] * pc_x[j] * fl1_fx - 0.5 * fl1_fx * pa_y[j] * pc_y[j] - 0.5 * fl1_fx * pc_y[j] * pb_y[j] - 0.5 * fl1_fx * pa_y[j] * pb_y[j] - pa_xxy[j] * pc_y[j] - pa_xx[j] * pc_y[j] * pb_y[j] - 2.0 * pa_xy[j] * pc_x[j] * pb_y[j]);

                t_xxy_y[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pa_x[j] * pc_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx + 0.5 * fl1_fx * pc_yy[j] + 0.5 * fl1_fx * pa_y[j] * pc_y[j] + 0.5 * fl1_fx * pc_y[j] * pb_y[j] + pa_xx[j] * pc_yy[j] + 2.0 * pa_xy[j] * pc_xy[j] + 2.0 * pa_x[j] * pc_xy[j] * pb_y[j] + pc_xx[j] * pa_y[j] * pb_y[j]);

                t_xxy_y[j] += fl_s_0_0_3 * (-0.5 * pc_xx[j] * fl1_fx - 0.5 * fl1_fx * pc_yy[j] - 2.0 * pa_x[j] * pc_xyy[j] - pc_xxy[j] * pa_y[j] - pc_xxy[j] * pb_y[j]);

                t_xxy_y[j] += fl_s_0_0_4 * pc_xxyy[j];

                t_xxy_z[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_y[j] * pb_z[j] + pa_xxy[j] * pb_z[j]);

                t_xxy_z[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pa_y[j] * pc_z[j] - 0.5 * fl1_fx * pc_y[j] * pb_z[j] - 0.5 * fl1_fx * pa_y[j] * pb_z[j] - pa_xxy[j] * pc_z[j] - pa_xx[j] * pc_y[j] * pb_z[j] - 2.0 * pa_xy[j] * pc_x[j] * pb_z[j]);

                t_xxy_z[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_yz[j] + 0.5 * fl1_fx * pa_y[j] * pc_z[j] + 0.5 * fl1_fx * pc_y[j] * pb_z[j] + pa_xx[j] * pc_yz[j] + 2.0 * pa_xy[j] * pc_xz[j] + 2.0 * pa_x[j] * pc_xy[j] * pb_z[j] + pc_xx[j] * pa_y[j] * pb_z[j]);

                t_xxy_z[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_yz[j] - 2.0 * pa_x[j] * pc_xyz[j] - pc_xxz[j] * pa_y[j] - pc_xxy[j] * pb_z[j]);

                t_xxy_z[j] += fl_s_0_0_4 * pc_xxyz[j];

                t_xxz_x[j] = fl_s_0_0_0 * (pa_xz[j] * fl1_fx + 0.5 * fl1_fx * pa_z[j] * pb_x[j] + pa_xxz[j] * pb_x[j]);

                t_xxz_x[j] += fl_s_0_0_1 * (-pa_x[j] * fl1_fx * pc_z[j] - pa_xz[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pa_z[j] - 0.5 * fl1_fx * pc_z[j] * pb_x[j] - 0.5 * fl1_fx * pa_z[j] * pb_x[j] - pa_xxz[j] * pc_x[j] - pa_xx[j] * pc_z[j] * pb_x[j] - 2.0 * pa_xz[j] * pc_x[j] * pb_x[j]);

                t_xxz_x[j] += fl_s_0_0_2 * (pa_x[j] * fl1_fx * pc_z[j] + 1.5 * pc_xz[j] * fl1_fx + 1.5 * pc_x[j] * fl1_fx * pa_z[j] + 0.5 * fl1_fx * pc_z[j] * pb_x[j] + pa_xx[j] * pc_xz[j] + 2.0 * pa_xz[j] * pc_xx[j] + 2.0 * pa_x[j] * pc_xz[j] * pb_x[j] + pc_xx[j] * pa_z[j] * pb_x[j]);

                t_xxz_x[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - 2.0 * pa_x[j] * pc_xxz[j] - pc_xxx[j] * pa_z[j] - pc_xxz[j] * pb_x[j]);

                t_xxz_x[j] += fl_s_0_0_4 * pc_xxxz[j];

                t_xxz_y[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_z[j] * pb_y[j] + pa_xxz[j] * pb_y[j]);

                t_xxz_y[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pa_z[j] * pc_y[j] - 0.5 * fl1_fx * pc_z[j] * pb_y[j] - 0.5 * fl1_fx * pa_z[j] * pb_y[j] - pa_xxz[j] * pc_y[j] - pa_xx[j] * pc_z[j] * pb_y[j] - 2.0 * pa_xz[j] * pc_x[j] * pb_y[j]);

                t_xxz_y[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_yz[j] + 0.5 * fl1_fx * pa_z[j] * pc_y[j] + 0.5 * fl1_fx * pc_z[j] * pb_y[j] + pa_xx[j] * pc_yz[j] + 2.0 * pa_xz[j] * pc_xy[j] + 2.0 * pa_x[j] * pc_xz[j] * pb_y[j] + pc_xx[j] * pa_z[j] * pb_y[j]);

                t_xxz_y[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_yz[j] - 2.0 * pa_x[j] * pc_xyz[j] - pc_xxy[j] * pa_z[j] - pc_xxz[j] * pb_y[j]);

                t_xxz_y[j] += fl_s_0_0_4 * pc_xxyz[j];

                t_xxz_z[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * fl1_fx * pa_z[j] * pb_z[j] + pa_xxz[j] * pb_z[j]);

                t_xxz_z[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_xx[j] * fl1_fx - pa_x[j] * pc_x[j] * fl1_fx - 0.5 * fl1_fx * pa_z[j] * pc_z[j] - 0.5 * fl1_fx * pc_z[j] * pb_z[j] - 0.5 * fl1_fx * pa_z[j] * pb_z[j] - pa_xxz[j] * pc_z[j] - pa_xx[j] * pc_z[j] * pb_z[j] - 2.0 * pa_xz[j] * pc_x[j] * pb_z[j]);

                t_xxz_z[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pa_x[j] * pc_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx + 0.5 * fl1_fx * pc_zz[j] + 0.5 * fl1_fx * pa_z[j] * pc_z[j] + 0.5 * fl1_fx * pc_z[j] * pb_z[j] + pa_xx[j] * pc_zz[j] + 2.0 * pa_xz[j] * pc_xz[j] + 2.0 * pa_x[j] * pc_xz[j] * pb_z[j] + pc_xx[j] * pa_z[j] * pb_z[j]);

                t_xxz_z[j] += fl_s_0_0_3 * (-0.5 * pc_xx[j] * fl1_fx - 0.5 * fl1_fx * pc_zz[j] - 2.0 * pa_x[j] * pc_xzz[j] - pc_xxz[j] * pa_z[j] - pc_xxz[j] * pb_z[j]);

                t_xxz_z[j] += fl_s_0_0_4 * pc_xxzz[j];

                t_xyy_x[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * fl1_fx * pa_yy[j] + 0.5 * pa_x[j] * fl1_fx * pb_x[j] + pa_xyy[j] * pb_x[j]);

                t_xyy_x[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - fl1_fx * pa_y[j] * pc_y[j] - 0.5 * fl1_fx * pa_yy[j] - 0.5 * pa_x[j] * fl1_fx * pc_x[j] - 0.5 * pa_x[j] * fl1_fx * pb_x[j] - 0.5 * pc_x[j] * fl1_fx * pb_x[j] - pa_xyy[j] * pc_x[j] - 2.0 * pa_xy[j] * pc_y[j] * pb_x[j] - pc_x[j] * pa_yy[j] * pb_x[j]);

                t_xyy_x[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * fl1_fx * pc_yy[j] + fl1_fx * pa_y[j] * pc_y[j] + 0.5 * pa_x[j] * fl1_fx * pc_x[j] + 0.5 * pc_xx[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pb_x[j] + 2.0 * pa_xy[j] * pc_xy[j] + pa_x[j] * pc_yy[j] * pb_x[j] + pc_xx[j] * pa_yy[j] + 2.0 * pc_xy[j] * pa_y[j] * pb_x[j]);

                t_xyy_x[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_yy[j] - 0.5 * pc_xx[j] * fl1_fx - pa_x[j] * pc_xyy[j] - 2.0 * pc_xxy[j] * pa_y[j] - pc_xyy[j] * pb_x[j]);

                t_xyy_x[j] += fl_s_0_0_4 * pc_xxyy[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForFP_10_20(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            auto pc_xxy = pcDistances.data(34 * idx + 10);

            auto pc_xxz = pcDistances.data(34 * idx + 11);

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            auto pc_xzz = pcDistances.data(34 * idx + 14);

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxyz = pcDistances.data(34 * idx + 23);

            auto pc_xxzz = pcDistances.data(34 * idx + 24);

            auto pc_xyyy = pcDistances.data(34 * idx + 25);

            auto pc_xyyz = pcDistances.data(34 * idx + 26);

            auto pc_xyzz = pcDistances.data(34 * idx + 27);

            auto pc_xzzz = pcDistances.data(34 * idx + 28);

            auto pc_yyyy = pcDistances.data(34 * idx + 29);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(5 * idx);

            auto s_0_0_1 = auxBuffer.data(5 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(5 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(5 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(5 * idx + 4);

            // set up pointers to integrals

            auto t_xyy_y = primBuffer.data(30 * idx + 10);

            auto t_xyy_z = primBuffer.data(30 * idx + 11);

            auto t_xyz_x = primBuffer.data(30 * idx + 12);

            auto t_xyz_y = primBuffer.data(30 * idx + 13);

            auto t_xyz_z = primBuffer.data(30 * idx + 14);

            auto t_xzz_x = primBuffer.data(30 * idx + 15);

            auto t_xzz_y = primBuffer.data(30 * idx + 16);

            auto t_xzz_z = primBuffer.data(30 * idx + 17);

            auto t_yyy_x = primBuffer.data(30 * idx + 18);

            auto t_yyy_y = primBuffer.data(30 * idx + 19);

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyz, pa_xz, pa_xzz, pa_y, pa_yy, pa_yyy, pa_yz, pa_z, \
                                     pa_zz, pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xxy, pc_xxyz, pc_xxz, pc_xxzz, pc_xy, pc_xyy, \
                                     pc_xyyy, pc_xyyz, pc_xyz, pc_xyzz, pc_xz, pc_xzz, pc_xzzz, pc_y, pc_yy, pc_yyy, pc_yyyy, \
                                     pc_yyz, pc_yz, pc_yzz, pc_z, pc_zz, pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, \
                                     s_0_0_4, t_xyy_y, t_xyy_z, t_xyz_x, t_xyz_y, t_xyz_z, t_xzz_x, t_xzz_y, t_xzz_z, \
                                     t_yyy_x, t_yyy_y: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xyy_y[j] = fl_s_0_0_0 * (pa_xy[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_y[j] + pa_xyy[j] * pb_y[j]);

                t_xyy_y[j] += fl_s_0_0_1 * (-pa_xy[j] * fl1_fx - 1.5 * pa_x[j] * pc_y[j] * fl1_fx - pc_x[j] * pa_y[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fx * pb_y[j] - 0.5 * pc_x[j] * fl1_fx * pb_y[j] - pa_xyy[j] * pc_y[j] - 2.0 * pa_xy[j] * pc_y[j] * pb_y[j] - pc_x[j] * pa_yy[j] * pb_y[j]);

                t_xyy_y[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * pc_y[j] * fl1_fx + pc_x[j] * pa_y[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pb_y[j] + 2.0 * pa_xy[j] * pc_yy[j] + pa_x[j] * pc_yy[j] * pb_y[j] + pc_xy[j] * pa_yy[j] + 2.0 * pc_xy[j] * pa_y[j] * pb_y[j]);

                t_xyy_y[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - pa_x[j] * pc_yyy[j] - 2.0 * pc_xyy[j] * pa_y[j] - pc_xyy[j] * pb_y[j]);

                t_xyy_y[j] += fl_s_0_0_4 * pc_xyyy[j];

                t_xyy_z[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx * pb_z[j] + pa_xyy[j] * pb_z[j]);

                t_xyy_z[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx * pc_z[j] - 0.5 * pa_x[j] * fl1_fx * pb_z[j] - 0.5 * pc_x[j] * fl1_fx * pb_z[j] - pa_xyy[j] * pc_z[j] - 2.0 * pa_xy[j] * pc_y[j] * pb_z[j] - pc_x[j] * pa_yy[j] * pb_z[j]);

                t_xyy_z[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl1_fx * pc_z[j] + 0.5 * pc_xz[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pb_z[j] + 2.0 * pa_xy[j] * pc_yz[j] + pa_x[j] * pc_yy[j] * pb_z[j] + pc_xz[j] * pa_yy[j] + 2.0 * pc_xy[j] * pa_y[j] * pb_z[j]);

                t_xyy_z[j] += fl_s_0_0_3 * (-0.5 * pc_xz[j] * fl1_fx - pa_x[j] * pc_yyz[j] - 2.0 * pc_xyz[j] * pa_y[j] - pc_xyy[j] * pb_z[j]);

                t_xyy_z[j] += fl_s_0_0_4 * pc_xyyz[j];

                t_xyz_x[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_yz[j] + pa_xyz[j] * pb_x[j]);

                t_xyz_x[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pa_y[j] * pc_z[j] - 0.5 * fl1_fx * pc_y[j] * pa_z[j] - 0.5 * fl1_fx * pa_yz[j] - pa_xyz[j] * pc_x[j] - pa_xy[j] * pc_z[j] * pb_x[j] - pa_xz[j] * pc_y[j] * pb_x[j] - pc_x[j] * pa_yz[j] * pb_x[j]);

                t_xyz_x[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_yz[j] + 0.5 * fl1_fx * pa_y[j] * pc_z[j] + 0.5 * fl1_fx * pc_y[j] * pa_z[j] + pa_xy[j] * pc_xz[j] + pa_xz[j] * pc_xy[j] + pa_x[j] * pc_yz[j] * pb_x[j] + pc_xx[j] * pa_yz[j] + pc_xz[j] * pa_y[j] * pb_x[j] + pc_xy[j] * pa_z[j] * pb_x[j]);

                t_xyz_x[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_yz[j] - pa_x[j] * pc_xyz[j] - pc_xxz[j] * pa_y[j] - pc_xxy[j] * pa_z[j] - pc_xyz[j] * pb_x[j]);

                t_xyz_x[j] += fl_s_0_0_4 * pc_xxyz[j];

                t_xyz_y[j] = fl_s_0_0_0 * (0.5 * pa_xz[j] * fl1_fx + pa_xyz[j] * pb_y[j]);

                t_xyz_y[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx * pc_z[j] - 0.5 * pa_xz[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx * pa_z[j] - pa_xyz[j] * pc_y[j] - pa_xy[j] * pc_z[j] * pb_y[j] - pa_xz[j] * pc_y[j] * pb_y[j] - pc_x[j] * pa_yz[j] * pb_y[j]);

                t_xyz_y[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl1_fx * pc_z[j] + 0.5 * pc_xz[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pa_z[j] + pa_xy[j] * pc_yz[j] + pa_xz[j] * pc_yy[j] + pa_x[j] * pc_yz[j] * pb_y[j] + pc_xy[j] * pa_yz[j] + pc_xz[j] * pa_y[j] * pb_y[j] + pc_xy[j] * pa_z[j] * pb_y[j]);

                t_xyz_y[j] += fl_s_0_0_3 * (-0.5 * pc_xz[j] * fl1_fx - pa_x[j] * pc_yyz[j] - pc_xyz[j] * pa_y[j] - pc_xyy[j] * pa_z[j] - pc_xyz[j] * pb_y[j]);

                t_xyz_y[j] += fl_s_0_0_4 * pc_xyyz[j];

                t_xyz_z[j] = fl_s_0_0_0 * (0.5 * pa_xy[j] * fl1_fx + pa_xyz[j] * pb_z[j]);

                t_xyz_z[j] += fl_s_0_0_1 * (-0.5 * pa_xy[j] * fl1_fx - 0.5 * pa_x[j] * pc_y[j] * fl1_fx - 0.5 * pc_x[j] * pa_y[j] * fl1_fx - pa_xyz[j] * pc_z[j] - pa_xy[j] * pc_z[j] * pb_z[j] - pa_xz[j] * pc_y[j] * pb_z[j] - pc_x[j] * pa_yz[j] * pb_z[j]);

                t_xyz_z[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * pc_y[j] * fl1_fx + 0.5 * pc_x[j] * pa_y[j] * fl1_fx + 0.5 * pc_xy[j] * fl1_fx + pa_xy[j] * pc_zz[j] + pa_xz[j] * pc_yz[j] + pa_x[j] * pc_yz[j] * pb_z[j] + pc_xz[j] * pa_yz[j] + pc_xz[j] * pa_y[j] * pb_z[j] + pc_xy[j] * pa_z[j] * pb_z[j]);

                t_xyz_z[j] += fl_s_0_0_3 * (-0.5 * pc_xy[j] * fl1_fx - pa_x[j] * pc_yzz[j] - pc_xzz[j] * pa_y[j] - pc_xyz[j] * pa_z[j] - pc_xyz[j] * pb_z[j]);

                t_xyz_z[j] += fl_s_0_0_4 * pc_xyzz[j];

                t_xzz_x[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * fl1_fx * pa_zz[j] + 0.5 * pa_x[j] * fl1_fx * pb_x[j] + pa_xzz[j] * pb_x[j]);

                t_xzz_x[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - fl1_fx * pa_z[j] * pc_z[j] - 0.5 * fl1_fx * pa_zz[j] - 0.5 * pa_x[j] * fl1_fx * pc_x[j] - 0.5 * pa_x[j] * fl1_fx * pb_x[j] - 0.5 * pc_x[j] * fl1_fx * pb_x[j] - pa_xzz[j] * pc_x[j] - 2.0 * pa_xz[j] * pc_z[j] * pb_x[j] - pc_x[j] * pa_zz[j] * pb_x[j]);

                t_xzz_x[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * fl1_fx * pc_zz[j] + fl1_fx * pa_z[j] * pc_z[j] + 0.5 * pa_x[j] * fl1_fx * pc_x[j] + 0.5 * pc_xx[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pb_x[j] + 2.0 * pa_xz[j] * pc_xz[j] + pa_x[j] * pc_zz[j] * pb_x[j] + pc_xx[j] * pa_zz[j] + 2.0 * pc_xz[j] * pa_z[j] * pb_x[j]);

                t_xzz_x[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_zz[j] - 0.5 * pc_xx[j] * fl1_fx - pa_x[j] * pc_xzz[j] - 2.0 * pc_xxz[j] * pa_z[j] - pc_xzz[j] * pb_x[j]);

                t_xzz_x[j] += fl_s_0_0_4 * pc_xxzz[j];

                t_xzz_y[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx * pb_y[j] + pa_xzz[j] * pb_y[j]);

                t_xzz_y[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx * pc_y[j] - 0.5 * pa_x[j] * fl1_fx * pb_y[j] - 0.5 * pc_x[j] * fl1_fx * pb_y[j] - pa_xzz[j] * pc_y[j] - 2.0 * pa_xz[j] * pc_z[j] * pb_y[j] - pc_x[j] * pa_zz[j] * pb_y[j]);

                t_xzz_y[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl1_fx * pc_y[j] + 0.5 * pc_xy[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pb_y[j] + 2.0 * pa_xz[j] * pc_yz[j] + pa_x[j] * pc_zz[j] * pb_y[j] + pc_xy[j] * pa_zz[j] + 2.0 * pc_xz[j] * pa_z[j] * pb_y[j]);

                t_xzz_y[j] += fl_s_0_0_3 * (-0.5 * pc_xy[j] * fl1_fx - pa_x[j] * pc_yzz[j] - 2.0 * pc_xyz[j] * pa_z[j] - pc_xzz[j] * pb_y[j]);

                t_xzz_y[j] += fl_s_0_0_4 * pc_xyzz[j];

                t_xzz_z[j] = fl_s_0_0_0 * (pa_xz[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_z[j] + pa_xzz[j] * pb_z[j]);

                t_xzz_z[j] += fl_s_0_0_1 * (-pa_xz[j] * fl1_fx - 1.5 * pa_x[j] * pc_z[j] * fl1_fx - pc_x[j] * pa_z[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fx * pb_z[j] - 0.5 * pc_x[j] * fl1_fx * pb_z[j] - pa_xzz[j] * pc_z[j] - 2.0 * pa_xz[j] * pc_z[j] * pb_z[j] - pc_x[j] * pa_zz[j] * pb_z[j]);

                t_xzz_z[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * pc_z[j] * fl1_fx + pc_x[j] * pa_z[j] * fl1_fx + 1.5 * pc_xz[j] * fl1_fx + 0.5 * pc_x[j] * fl1_fx * pb_z[j] + 2.0 * pa_xz[j] * pc_zz[j] + pa_x[j] * pc_zz[j] * pb_z[j] + pc_xz[j] * pa_zz[j] + 2.0 * pc_xz[j] * pa_z[j] * pb_z[j]);

                t_xzz_z[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - pa_x[j] * pc_zzz[j] - 2.0 * pc_xzz[j] * pa_z[j] - pc_xzz[j] * pb_z[j]);

                t_xzz_z[j] += fl_s_0_0_4 * pc_xzzz[j];

                t_yyy_x[j] = fl_s_0_0_0 * (1.5 * pa_y[j] * fl1_fx * pb_x[j] + pa_yyy[j] * pb_x[j]);

                t_yyy_x[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl1_fx * pc_x[j] - 1.5 * pa_y[j] * fl1_fx * pb_x[j] - 1.5 * pc_y[j] * fl1_fx * pb_x[j] - pa_yyy[j] * pc_x[j] - 3.0 * pa_yy[j] * pc_y[j] * pb_x[j]);

                t_yyy_x[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * fl1_fx * pc_x[j] + 1.5 * pc_xy[j] * fl1_fx + 1.5 * pc_y[j] * fl1_fx * pb_x[j] + 3.0 * pa_yy[j] * pc_xy[j] + 3.0 * pa_y[j] * pc_yy[j] * pb_x[j]);

                t_yyy_x[j] += fl_s_0_0_3 * (-1.5 * pc_xy[j] * fl1_fx - 3.0 * pa_y[j] * pc_xyy[j] - pc_yyy[j] * pb_x[j]);

                t_yyy_x[j] += fl_s_0_0_4 * pc_xyyy[j];

                t_yyy_y[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 1.5 * pa_yy[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pb_y[j] + pa_yyy[j] * pb_y[j]);

                t_yyy_y[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 1.5 * pa_yy[j] * fl1_fx - 4.5 * pa_y[j] * pc_y[j] * fl1_fx - 1.5 * pa_y[j] * fl1_fx * pb_y[j] - 1.5 * pc_y[j] * fl1_fx * pb_y[j] - pa_yyy[j] * pc_y[j] - 3.0 * pa_yy[j] * pc_y[j] * pb_y[j]);

                t_yyy_y[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 4.5 * pa_y[j] * pc_y[j] * fl1_fx + 3.0 * pc_yy[j] * fl1_fx + 1.5 * pc_y[j] * fl1_fx * pb_y[j] + 3.0 * pa_yy[j] * pc_yy[j] + 3.0 * pa_y[j] * pc_yy[j] * pb_y[j]);

                t_yyy_y[j] += fl_s_0_0_3 * (-3.0 * pc_yy[j] * fl1_fx - 3.0 * pa_y[j] * pc_yyy[j] - pc_yyy[j] * pb_y[j]);

                t_yyy_y[j] += fl_s_0_0_4 * pc_yyyy[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForFP_20_30(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
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

            auto fx = osFactors.data(3 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(34 * idx);

            auto pc_y = pcDistances.data(34 * idx + 1);

            auto pc_z = pcDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xy = pcDistances.data(34 * idx + 4);

            auto pc_xz = pcDistances.data(34 * idx + 5);

            auto pc_yy = pcDistances.data(34 * idx + 6);

            auto pc_yz = pcDistances.data(34 * idx + 7);

            auto pc_zz = pcDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xyy = pcDistances.data(34 * idx + 12);

            auto pc_xyz = pcDistances.data(34 * idx + 13);

            auto pc_xzz = pcDistances.data(34 * idx + 14);

            auto pc_yyy = pcDistances.data(34 * idx + 15);

            auto pc_yyz = pcDistances.data(34 * idx + 16);

            auto pc_yzz = pcDistances.data(34 * idx + 17);

            auto pc_zzz = pcDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xyyz = pcDistances.data(34 * idx + 26);

            auto pc_xyzz = pcDistances.data(34 * idx + 27);

            auto pc_xzzz = pcDistances.data(34 * idx + 28);

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

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pa_zzz, pb_x, pb_y, \
                                     pb_z, pc_x, pc_xy, pc_xyy, pc_xyyz, pc_xyz, pc_xyzz, pc_xz, pc_xzz, pc_xzzz, pc_y, \
                                     pc_yy, pc_yyy, pc_yyyz, pc_yyz, pc_yyzz, pc_yz, pc_yzz, pc_yzzz, pc_z, pc_zz, pc_zzz, \
                                     pc_zzzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, t_yyy_z, t_yyz_x, t_yyz_y, \
                                     t_yyz_z, t_yzz_x, t_yzz_y, t_yzz_z, t_zzz_x, t_zzz_y, t_zzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yyy_z[j] = fl_s_0_0_0 * (1.5 * pa_y[j] * fl1_fx * pb_z[j] + pa_yyy[j] * pb_z[j]);

                t_yyy_z[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl1_fx * pc_z[j] - 1.5 * pa_y[j] * fl1_fx * pb_z[j] - 1.5 * pc_y[j] * fl1_fx * pb_z[j] - pa_yyy[j] * pc_z[j] - 3.0 * pa_yy[j] * pc_y[j] * pb_z[j]);

                t_yyy_z[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * fl1_fx * pc_z[j] + 1.5 * pc_yz[j] * fl1_fx + 1.5 * pc_y[j] * fl1_fx * pb_z[j] + 3.0 * pa_yy[j] * pc_yz[j] + 3.0 * pa_y[j] * pc_yy[j] * pb_z[j]);

                t_yyy_z[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - 3.0 * pa_y[j] * pc_yyz[j] - pc_yyy[j] * pb_z[j]);

                t_yyy_z[j] += fl_s_0_0_4 * pc_yyyz[j];

                t_yyz_x[j] = fl_s_0_0_0 * (0.5 * fl1_fx * pa_z[j] * pb_x[j] + pa_yyz[j] * pb_x[j]);

                t_yyz_x[j] += fl_s_0_0_1 * (-0.5 * fl1_fx * pa_z[j] * pc_x[j] - 0.5 * fl1_fx * pc_z[j] * pb_x[j] - 0.5 * fl1_fx * pa_z[j] * pb_x[j] - pa_yyz[j] * pc_x[j] - pa_yy[j] * pc_z[j] * pb_x[j] - 2.0 * pa_yz[j] * pc_y[j] * pb_x[j]);

                t_yyz_x[j] += fl_s_0_0_2 * (0.5 * fl1_fx * pc_xz[j] + 0.5 * fl1_fx * pa_z[j] * pc_x[j] + 0.5 * fl1_fx * pc_z[j] * pb_x[j] + pa_yy[j] * pc_xz[j] + 2.0 * pa_yz[j] * pc_xy[j] + 2.0 * pa_y[j] * pc_yz[j] * pb_x[j] + pc_yy[j] * pa_z[j] * pb_x[j]);

                t_yyz_x[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_xz[j] - 2.0 * pa_y[j] * pc_xyz[j] - pc_xyy[j] * pa_z[j] - pc_yyz[j] * pb_x[j]);

                t_yyz_x[j] += fl_s_0_0_4 * pc_xyyz[j];

                t_yyz_y[j] = fl_s_0_0_0 * (pa_yz[j] * fl1_fx + 0.5 * fl1_fx * pa_z[j] * pb_y[j] + pa_yyz[j] * pb_y[j]);

                t_yyz_y[j] += fl_s_0_0_1 * (-pa_y[j] * fl1_fx * pc_z[j] - pa_yz[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx * pa_z[j] - 0.5 * fl1_fx * pc_z[j] * pb_y[j] - 0.5 * fl1_fx * pa_z[j] * pb_y[j] - pa_yyz[j] * pc_y[j] - pa_yy[j] * pc_z[j] * pb_y[j] - 2.0 * pa_yz[j] * pc_y[j] * pb_y[j]);

                t_yyz_y[j] += fl_s_0_0_2 * (pa_y[j] * fl1_fx * pc_z[j] + 1.5 * pc_yz[j] * fl1_fx + 1.5 * pc_y[j] * fl1_fx * pa_z[j] + 0.5 * fl1_fx * pc_z[j] * pb_y[j] + pa_yy[j] * pc_yz[j] + 2.0 * pa_yz[j] * pc_yy[j] + 2.0 * pa_y[j] * pc_yz[j] * pb_y[j] + pc_yy[j] * pa_z[j] * pb_y[j]);

                t_yyz_y[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - 2.0 * pa_y[j] * pc_yyz[j] - pc_yyy[j] * pa_z[j] - pc_yyz[j] * pb_y[j]);

                t_yyz_y[j] += fl_s_0_0_4 * pc_yyyz[j];

                t_yyz_z[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * pa_yy[j] * fl1_fx + 0.5 * fl1_fx * pa_z[j] * pb_z[j] + pa_yyz[j] * pb_z[j]);

                t_yyz_z[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - 0.5 * pa_yy[j] * fl1_fx - pa_y[j] * pc_y[j] * fl1_fx - 0.5 * fl1_fx * pa_z[j] * pc_z[j] - 0.5 * fl1_fx * pc_z[j] * pb_z[j] - 0.5 * fl1_fx * pa_z[j] * pb_z[j] - pa_yyz[j] * pc_z[j] - pa_yy[j] * pc_z[j] * pb_z[j] - 2.0 * pa_yz[j] * pc_y[j] * pb_z[j]);

                t_yyz_z[j] += fl_s_0_0_2 * (0.25 * fl2_fx + pa_y[j] * pc_y[j] * fl1_fx + 0.5 * pc_yy[j] * fl1_fx + 0.5 * fl1_fx * pc_zz[j] + 0.5 * fl1_fx * pa_z[j] * pc_z[j] + 0.5 * fl1_fx * pc_z[j] * pb_z[j] + pa_yy[j] * pc_zz[j] + 2.0 * pa_yz[j] * pc_yz[j] + 2.0 * pa_y[j] * pc_yz[j] * pb_z[j] + pc_yy[j] * pa_z[j] * pb_z[j]);

                t_yyz_z[j] += fl_s_0_0_3 * (-0.5 * pc_yy[j] * fl1_fx - 0.5 * fl1_fx * pc_zz[j] - 2.0 * pa_y[j] * pc_yzz[j] - pc_yyz[j] * pa_z[j] - pc_yyz[j] * pb_z[j]);

                t_yyz_z[j] += fl_s_0_0_4 * pc_yyzz[j];

                t_yzz_x[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl1_fx * pb_x[j] + pa_yzz[j] * pb_x[j]);

                t_yzz_x[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl1_fx * pc_x[j] - 0.5 * pa_y[j] * fl1_fx * pb_x[j] - 0.5 * pc_y[j] * fl1_fx * pb_x[j] - pa_yzz[j] * pc_x[j] - 2.0 * pa_yz[j] * pc_z[j] * pb_x[j] - pc_y[j] * pa_zz[j] * pb_x[j]);

                t_yzz_x[j] += fl_s_0_0_2 * (0.5 * pa_y[j] * fl1_fx * pc_x[j] + 0.5 * pc_xy[j] * fl1_fx + 0.5 * pc_y[j] * fl1_fx * pb_x[j] + 2.0 * pa_yz[j] * pc_xz[j] + pa_y[j] * pc_zz[j] * pb_x[j] + pc_xy[j] * pa_zz[j] + 2.0 * pc_yz[j] * pa_z[j] * pb_x[j]);

                t_yzz_x[j] += fl_s_0_0_3 * (-0.5 * pc_xy[j] * fl1_fx - pa_y[j] * pc_xzz[j] - 2.0 * pc_xyz[j] * pa_z[j] - pc_yzz[j] * pb_x[j]);

                t_yzz_x[j] += fl_s_0_0_4 * pc_xyzz[j];

                t_yzz_y[j] = fl_s_0_0_0 * (0.25 * fl2_fx + 0.5 * fl1_fx * pa_zz[j] + 0.5 * pa_y[j] * fl1_fx * pb_y[j] + pa_yzz[j] * pb_y[j]);

                t_yzz_y[j] += fl_s_0_0_1 * (-0.5 * fl2_fx - fl1_fx * pa_z[j] * pc_z[j] - 0.5 * fl1_fx * pa_zz[j] - 0.5 * pa_y[j] * fl1_fx * pc_y[j] - 0.5 * pa_y[j] * fl1_fx * pb_y[j] - 0.5 * pc_y[j] * fl1_fx * pb_y[j] - pa_yzz[j] * pc_y[j] - 2.0 * pa_yz[j] * pc_z[j] * pb_y[j] - pc_y[j] * pa_zz[j] * pb_y[j]);

                t_yzz_y[j] += fl_s_0_0_2 * (0.25 * fl2_fx + 0.5 * fl1_fx * pc_zz[j] + fl1_fx * pa_z[j] * pc_z[j] + 0.5 * pa_y[j] * fl1_fx * pc_y[j] + 0.5 * pc_yy[j] * fl1_fx + 0.5 * pc_y[j] * fl1_fx * pb_y[j] + 2.0 * pa_yz[j] * pc_yz[j] + pa_y[j] * pc_zz[j] * pb_y[j] + pc_yy[j] * pa_zz[j] + 2.0 * pc_yz[j] * pa_z[j] * pb_y[j]);

                t_yzz_y[j] += fl_s_0_0_3 * (-0.5 * fl1_fx * pc_zz[j] - 0.5 * pc_yy[j] * fl1_fx - pa_y[j] * pc_yzz[j] - 2.0 * pc_yyz[j] * pa_z[j] - pc_yzz[j] * pb_y[j]);

                t_yzz_y[j] += fl_s_0_0_4 * pc_yyzz[j];

                t_yzz_z[j] = fl_s_0_0_0 * (pa_yz[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pb_z[j] + pa_yzz[j] * pb_z[j]);

                t_yzz_z[j] += fl_s_0_0_1 * (-pa_yz[j] * fl1_fx - 1.5 * pa_y[j] * pc_z[j] * fl1_fx - pc_y[j] * pa_z[j] * fl1_fx - 0.5 * pa_y[j] * fl1_fx * pb_z[j] - 0.5 * pc_y[j] * fl1_fx * pb_z[j] - pa_yzz[j] * pc_z[j] - 2.0 * pa_yz[j] * pc_z[j] * pb_z[j] - pc_y[j] * pa_zz[j] * pb_z[j]);

                t_yzz_z[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * pc_z[j] * fl1_fx + pc_y[j] * pa_z[j] * fl1_fx + 1.5 * pc_yz[j] * fl1_fx + 0.5 * pc_y[j] * fl1_fx * pb_z[j] + 2.0 * pa_yz[j] * pc_zz[j] + pa_y[j] * pc_zz[j] * pb_z[j] + pc_yz[j] * pa_zz[j] + 2.0 * pc_yz[j] * pa_z[j] * pb_z[j]);

                t_yzz_z[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - pa_y[j] * pc_zzz[j] - 2.0 * pc_yzz[j] * pa_z[j] - pc_yzz[j] * pb_z[j]);

                t_yzz_z[j] += fl_s_0_0_4 * pc_yzzz[j];

                t_zzz_x[j] = fl_s_0_0_0 * (1.5 * pa_z[j] * fl1_fx * pb_x[j] + pa_zzz[j] * pb_x[j]);

                t_zzz_x[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * fl1_fx * pc_x[j] - 1.5 * pa_z[j] * fl1_fx * pb_x[j] - 1.5 * pc_z[j] * fl1_fx * pb_x[j] - pa_zzz[j] * pc_x[j] - 3.0 * pa_zz[j] * pc_z[j] * pb_x[j]);

                t_zzz_x[j] += fl_s_0_0_2 * (1.5 * pa_z[j] * fl1_fx * pc_x[j] + 1.5 * pc_xz[j] * fl1_fx + 1.5 * pc_z[j] * fl1_fx * pb_x[j] + 3.0 * pa_zz[j] * pc_xz[j] + 3.0 * pa_z[j] * pc_zz[j] * pb_x[j]);

                t_zzz_x[j] += fl_s_0_0_3 * (-1.5 * pc_xz[j] * fl1_fx - 3.0 * pa_z[j] * pc_xzz[j] - pc_zzz[j] * pb_x[j]);

                t_zzz_x[j] += fl_s_0_0_4 * pc_xzzz[j];

                t_zzz_y[j] = fl_s_0_0_0 * (1.5 * pa_z[j] * fl1_fx * pb_y[j] + pa_zzz[j] * pb_y[j]);

                t_zzz_y[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * fl1_fx * pc_y[j] - 1.5 * pa_z[j] * fl1_fx * pb_y[j] - 1.5 * pc_z[j] * fl1_fx * pb_y[j] - pa_zzz[j] * pc_y[j] - 3.0 * pa_zz[j] * pc_z[j] * pb_y[j]);

                t_zzz_y[j] += fl_s_0_0_2 * (1.5 * pa_z[j] * fl1_fx * pc_y[j] + 1.5 * pc_yz[j] * fl1_fx + 1.5 * pc_z[j] * fl1_fx * pb_y[j] + 3.0 * pa_zz[j] * pc_yz[j] + 3.0 * pa_z[j] * pc_zz[j] * pb_y[j]);

                t_zzz_y[j] += fl_s_0_0_3 * (-1.5 * pc_yz[j] * fl1_fx - 3.0 * pa_z[j] * pc_yzz[j] - pc_zzz[j] * pb_y[j]);

                t_zzz_y[j] += fl_s_0_0_4 * pc_yzzz[j];

                t_zzz_z[j] = fl_s_0_0_0 * (0.75 * fl2_fx + 1.5 * pa_zz[j] * fl1_fx + 1.5 * pa_z[j] * fl1_fx * pb_z[j] + pa_zzz[j] * pb_z[j]);

                t_zzz_z[j] += fl_s_0_0_1 * (-1.5 * fl2_fx - 1.5 * pa_zz[j] * fl1_fx - 4.5 * pa_z[j] * pc_z[j] * fl1_fx - 1.5 * pa_z[j] * fl1_fx * pb_z[j] - 1.5 * pc_z[j] * fl1_fx * pb_z[j] - pa_zzz[j] * pc_z[j] - 3.0 * pa_zz[j] * pc_z[j] * pb_z[j]);

                t_zzz_z[j] += fl_s_0_0_2 * (0.75 * fl2_fx + 4.5 * pa_z[j] * pc_z[j] * fl1_fx + 3.0 * pc_zz[j] * fl1_fx + 1.5 * pc_z[j] * fl1_fx * pb_z[j] + 3.0 * pa_zz[j] * pc_zz[j] + 3.0 * pa_z[j] * pc_zz[j] * pb_z[j]);

                t_zzz_z[j] += fl_s_0_0_3 * (-3.0 * pc_zz[j] * fl1_fx - 3.0 * pa_z[j] * pc_zzz[j] - pc_zzz[j] * pb_z[j]);

                t_zzz_z[j] += fl_s_0_0_4 * pc_zzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForPG_0_5(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPG_5_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                    braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPG_10_15(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPG_15_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPG_20_25(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPG_25_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPG_30_35(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPG_35_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForPG_40_45(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compNuclearPotentialForPG_0_5(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
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

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

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

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(55 * idx + 9);

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(55 * idx + 19);

            auto pc_xxxy = pcDistances.data(55 * idx + 20);

            auto pc_xxxz = pcDistances.data(55 * idx + 21);

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxxxx = pcDistances.data(55 * idx + 34);

            auto pc_xxxxy = pcDistances.data(55 * idx + 35);

            auto pc_xxxxz = pcDistances.data(55 * idx + 36);

            auto pc_xxxyy = pcDistances.data(55 * idx + 37);

            auto pc_xxxyz = pcDistances.data(55 * idx + 38);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_x_xxxx = primBuffer.data(45 * idx);

            auto t_x_xxxy = primBuffer.data(45 * idx + 1);

            auto t_x_xxxz = primBuffer.data(45 * idx + 2);

            auto t_x_xxyy = primBuffer.data(45 * idx + 3);

            auto t_x_xxyz = primBuffer.data(45 * idx + 4);

            // Batch of Integrals (0,5)

            #pragma omp simd aligned(fx, pa_x, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pc_x, pc_xx, \
                                     pc_xxx, pc_xxxx, pc_xxxxx, pc_xxxxy, pc_xxxxz, pc_xxxy, pc_xxxyy, pc_xxxyz, \
                                     pc_xxxz, pc_xxy, pc_xxyy, pc_xxyz, pc_xxz, pc_xy, pc_xyy, pc_xyz, pc_xz, pc_y, pc_yy, \
                                     pc_yz, pc_z, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_x_xxxx, \
                                     t_x_xxxy, t_x_xxxz, t_x_xxyy, t_x_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_x_xxxx[j] = fl_s_0_0_0 * (0.75 * pa_x[j] * fl2_fx + 3.0 * fl2_fx * pb_x[j] + 3.0 * pa_x[j] * pb_xx[j] * fl1_fx + 2.0 * fl1_fx * pb_xxx[j] + pa_x[j] * pb_xxxx[j]);

                t_x_xxxx[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl2_fx - 3.75 * pc_x[j] * fl2_fx - 6.0 * fl2_fx * pb_x[j] - 3.0 * pa_x[j] * pb_xx[j] * fl1_fx - 6.0 * pa_x[j] * pb_x[j] * pc_x[j] * fl1_fx - 9.0 * pc_x[j] * pb_xx[j] * fl1_fx - 2.0 * fl1_fx * pb_xxx[j] - 4.0 * pa_x[j] * pb_xxx[j] * pc_x[j] - pc_x[j] * pb_xxxx[j]);

                t_x_xxxx[j] += fl_s_0_0_2 * (0.75 * pa_x[j] * fl2_fx + 7.5 * pc_x[j] * fl2_fx + 3.0 * fl2_fx * pb_x[j] + 6.0 * pa_x[j] * pb_x[j] * pc_x[j] * fl1_fx + 3.0 * pa_x[j] * pc_xx[j] * fl1_fx + 9.0 * pc_x[j] * pb_xx[j] * fl1_fx + 12.0 * pc_xx[j] * pb_x[j] * fl1_fx + 6.0 * pa_x[j] * pb_xx[j] * pc_xx[j] + 4.0 * pc_xx[j] * pb_xxx[j]);

                t_x_xxxx[j] += fl_s_0_0_3 * (-3.75 * pc_x[j] * fl2_fx - 3.0 * pa_x[j] * pc_xx[j] * fl1_fx - 12.0 * pc_xx[j] * pb_x[j] * fl1_fx - 5.0 * pc_xxx[j] * fl1_fx - 4.0 * pa_x[j] * pb_x[j] * pc_xxx[j] - 6.0 * pc_xxx[j] * pb_xx[j]);

                t_x_xxxx[j] += fl_s_0_0_4 * (5.0 * pc_xxx[j] * fl1_fx + pa_x[j] * pc_xxxx[j] + 4.0 * pc_xxxx[j] * pb_x[j]);

                t_x_xxxx[j] += -fl_s_0_0_5 * pc_xxxxx[j];

                t_x_xxxy[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_y[j] + 1.5 * pa_x[j] * pb_xy[j] * fl1_fx + 1.5 * fl1_fx * pb_xxy[j] + pa_x[j] * pb_xxxy[j]);

                t_x_xxxy[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_y[j] - 1.5 * fl2_fx * pb_y[j] - 1.5 * pa_x[j] * pb_x[j] * fl1_fx * pc_y[j] - 1.5 * pa_x[j] * pb_xy[j] * fl1_fx - 1.5 * pa_x[j] * pc_x[j] * fl1_fx * pb_y[j] - 4.5 * pc_x[j] * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * pb_xx[j] * pc_y[j] - 1.5 * fl1_fx * pb_xxy[j] - pa_x[j] * pb_xxx[j] * pc_y[j] - 3.0 * pa_x[j] * pb_xxy[j] * pc_x[j] - pc_x[j] * pb_xxxy[j]);

                t_x_xxxy[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.75 * fl2_fx * pb_y[j] + 1.5 * pa_x[j] * pb_x[j] * fl1_fx * pc_y[j] + 1.5 * pa_x[j] * pc_xy[j] * fl1_fx + 1.5 * pa_x[j] * pc_x[j] * fl1_fx * pb_y[j] + 4.5 * pc_xy[j] * pb_x[j] * fl1_fx + 4.5 * pc_x[j] * pb_xy[j] * fl1_fx + 3.0 * pc_xx[j] * fl1_fx * pb_y[j] + 1.5 * fl1_fx * pb_xx[j] * pc_y[j] + 3.0 * pa_x[j] * pb_xx[j] * pc_xy[j] + 3.0 * pa_x[j] * pb_xy[j] * pc_xx[j] + pc_xy[j] * pb_xxx[j] + 3.0 * pc_xx[j] * pb_xxy[j]);

                t_x_xxxy[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 1.5 * pa_x[j] * pc_xy[j] * fl1_fx - 4.5 * pc_xy[j] * pb_x[j] * fl1_fx - 3.0 * pc_xxy[j] * fl1_fx - 3.0 * pc_xx[j] * fl1_fx * pb_y[j] - 3.0 * pa_x[j] * pb_x[j] * pc_xxy[j] - pa_x[j] * pc_xxx[j] * pb_y[j] - 3.0 * pc_xxy[j] * pb_xx[j] - 3.0 * pc_xxx[j] * pb_xy[j]);

                t_x_xxxy[j] += fl_s_0_0_4 * (3.0 * pc_xxy[j] * fl1_fx + pa_x[j] * pc_xxxy[j] + 3.0 * pc_xxxy[j] * pb_x[j] + pc_xxxx[j] * pb_y[j]);

                t_x_xxxy[j] += -fl_s_0_0_5 * pc_xxxxy[j];

                t_x_xxxz[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_z[j] + 1.5 * pa_x[j] * pb_xz[j] * fl1_fx + 1.5 * fl1_fx * pb_xxz[j] + pa_x[j] * pb_xxxz[j]);

                t_x_xxxz[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * fl2_fx * pb_z[j] - 1.5 * pa_x[j] * pb_x[j] * fl1_fx * pc_z[j] - 1.5 * pa_x[j] * pb_xz[j] * fl1_fx - 1.5 * pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] - 4.5 * pc_x[j] * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * pb_xx[j] * pc_z[j] - 1.5 * fl1_fx * pb_xxz[j] - pa_x[j] * pb_xxx[j] * pc_z[j] - 3.0 * pa_x[j] * pb_xxz[j] * pc_x[j] - pc_x[j] * pb_xxxz[j]);

                t_x_xxxz[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pb_z[j] + 1.5 * pa_x[j] * pb_x[j] * fl1_fx * pc_z[j] + 1.5 * pa_x[j] * pc_xz[j] * fl1_fx + 1.5 * pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] + 4.5 * pc_xz[j] * pb_x[j] * fl1_fx + 4.5 * pc_x[j] * pb_xz[j] * fl1_fx + 3.0 * pc_xx[j] * fl1_fx * pb_z[j] + 1.5 * fl1_fx * pb_xx[j] * pc_z[j] + 3.0 * pa_x[j] * pb_xx[j] * pc_xz[j] + 3.0 * pa_x[j] * pb_xz[j] * pc_xx[j] + pc_xz[j] * pb_xxx[j] + 3.0 * pc_xx[j] * pb_xxz[j]);

                t_x_xxxz[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * pa_x[j] * pc_xz[j] * fl1_fx - 4.5 * pc_xz[j] * pb_x[j] * fl1_fx - 3.0 * pc_xxz[j] * fl1_fx - 3.0 * pc_xx[j] * fl1_fx * pb_z[j] - 3.0 * pa_x[j] * pb_x[j] * pc_xxz[j] - pa_x[j] * pc_xxx[j] * pb_z[j] - 3.0 * pc_xxz[j] * pb_xx[j] - 3.0 * pc_xxx[j] * pb_xz[j]);

                t_x_xxxz[j] += fl_s_0_0_4 * (3.0 * pc_xxz[j] * fl1_fx + pa_x[j] * pc_xxxz[j] + 3.0 * pc_xxxz[j] * pb_x[j] + pc_xxxx[j] * pb_z[j]);

                t_x_xxxz[j] += -fl_s_0_0_5 * pc_xxxxz[j];

                t_x_xxyy[j] = fl_s_0_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * fl2_fx * pb_x[j] + 0.5 * pa_x[j] * pb_xx[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_yy[j] + fl1_fx * pb_xyy[j] + pa_x[j] * pb_xxyy[j]);

                t_x_xxyy[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - fl2_fx * pb_x[j] - 0.5 * pa_x[j] * pb_xx[j] * fl1_fx - pa_x[j] * pb_x[j] * pc_x[j] * fl1_fx - pa_x[j] * fl1_fx * pb_y[j] * pc_y[j] - 0.5 * pa_x[j] * fl1_fx * pb_yy[j] - 0.5 * pc_x[j] * pb_xx[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pb_yy[j] - 2.0 * fl1_fx * pb_xy[j] * pc_y[j] - fl1_fx * pb_xyy[j] - 2.0 * pa_x[j] * pb_xxy[j] * pc_y[j] - 2.0 * pa_x[j] * pb_xyy[j] * pc_x[j] - pc_x[j] * pb_xxyy[j]);

                t_x_xxyy[j] += fl_s_0_0_2 * (0.25 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + 0.5 * fl2_fx * pb_x[j] + pa_x[j] * pb_x[j] * pc_x[j] * fl1_fx + 0.5 * pa_x[j] * pc_xx[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pc_yy[j] + pa_x[j] * fl1_fx * pb_y[j] * pc_y[j] + 0.5 * pc_x[j] * pb_xx[j] * fl1_fx + pc_xx[j] * pb_x[j] * fl1_fx + 3.0 * pc_xy[j] * fl1_fx * pb_y[j] + 1.5 * pc_x[j] * fl1_fx * pb_yy[j] + fl1_fx * pb_x[j] * pc_yy[j] + 2.0 * fl1_fx * pb_xy[j] * pc_y[j] + pa_x[j] * pb_xx[j] * pc_yy[j] + 4.0 * pa_x[j] * pb_xy[j] * pc_xy[j] + pa_x[j] * pc_xx[j] * pb_yy[j] + 2.0 * pc_xy[j] * pb_xxy[j] + 2.0 * pc_xx[j] * pb_xyy[j]);

                t_x_xxyy[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - 0.5 * pa_x[j] * pc_xx[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fx * pc_yy[j] - pc_xx[j] * pb_x[j] * fl1_fx - 0.5 * pc_xxx[j] * fl1_fx - 1.5 * pc_xyy[j] * fl1_fx - 3.0 * pc_xy[j] * fl1_fx * pb_y[j] - fl1_fx * pb_x[j] * pc_yy[j] - 2.0 * pa_x[j] * pb_x[j] * pc_xyy[j] - 2.0 * pa_x[j] * pc_xxy[j] * pb_y[j] - pc_xyy[j] * pb_xx[j] - 4.0 * pc_xxy[j] * pb_xy[j] - pc_xxx[j] * pb_yy[j]);

                t_x_xxyy[j] += fl_s_0_0_4 * (0.5 * pc_xxx[j] * fl1_fx + 1.5 * pc_xyy[j] * fl1_fx + pa_x[j] * pc_xxyy[j] + 2.0 * pc_xxyy[j] * pb_x[j] + 2.0 * pc_xxxy[j] * pb_y[j]);

                t_x_xxyy[j] += -fl_s_0_0_5 * pc_xxxyy[j];

                t_x_xxyz[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl1_fx * pb_yz[j] + fl1_fx * pb_xyz[j] + pa_x[j] * pb_xxyz[j]);

                t_x_xxyz[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl1_fx * pb_y[j] * pc_z[j] - 0.5 * pa_x[j] * fl1_fx * pc_y[j] * pb_z[j] - 0.5 * pa_x[j] * fl1_fx * pb_yz[j] - 1.5 * pc_x[j] * fl1_fx * pb_yz[j] - fl1_fx * pb_xy[j] * pc_z[j] - fl1_fx * pb_xz[j] * pc_y[j] - fl1_fx * pb_xyz[j] - pa_x[j] * pb_xxy[j] * pc_z[j] - pa_x[j] * pb_xxz[j] * pc_y[j] - 2.0 * pa_x[j] * pb_xyz[j] * pc_x[j] - pc_x[j] * pb_xxyz[j]);

                t_x_xxyz[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl1_fx * pc_yz[j] + 0.5 * pa_x[j] * fl1_fx * pb_y[j] * pc_z[j] + 0.5 * pa_x[j] * fl1_fx * pc_y[j] * pb_z[j] + 1.5 * pc_xz[j] * fl1_fx * pb_y[j] + 1.5 * pc_xy[j] * fl1_fx * pb_z[j] + 1.5 * pc_x[j] * fl1_fx * pb_yz[j] + fl1_fx * pb_x[j] * pc_yz[j] + fl1_fx * pb_xy[j] * pc_z[j] + fl1_fx * pb_xz[j] * pc_y[j] + pa_x[j] * pb_xx[j] * pc_yz[j] + 2.0 * pa_x[j] * pb_xy[j] * pc_xz[j] + 2.0 * pa_x[j] * pb_xz[j] * pc_xy[j] + pa_x[j] * pc_xx[j] * pb_yz[j] + pc_xz[j] * pb_xxy[j] + pc_xy[j] * pb_xxz[j] + 2.0 * pc_xx[j] * pb_xyz[j]);

                t_x_xxyz[j] += fl_s_0_0_3 * (-0.5 * pa_x[j] * fl1_fx * pc_yz[j] - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pb_y[j] - 1.5 * pc_xy[j] * fl1_fx * pb_z[j] - fl1_fx * pb_x[j] * pc_yz[j] - 2.0 * pa_x[j] * pb_x[j] * pc_xyz[j] - pa_x[j] * pc_xxz[j] * pb_y[j] - pa_x[j] * pc_xxy[j] * pb_z[j] - pc_xyz[j] * pb_xx[j] - 2.0 * pc_xxz[j] * pb_xy[j] - 2.0 * pc_xxy[j] * pb_xz[j] - pc_xxx[j] * pb_yz[j]);

                t_x_xxyz[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_x[j] * pc_xxyz[j] + 2.0 * pc_xxyz[j] * pb_x[j] + pc_xxxz[j] * pb_y[j] + pc_xxxy[j] * pb_z[j]);

                t_x_xxyz[j] += -fl_s_0_0_5 * pc_xxxyz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG_5_10(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

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

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(55 * idx + 9);

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxz = pcDistances.data(55 * idx + 21);

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            auto pc_xzzz = pcDistances.data(55 * idx + 28);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxxzz = pcDistances.data(55 * idx + 39);

            auto pc_xxyyy = pcDistances.data(55 * idx + 40);

            auto pc_xxyyz = pcDistances.data(55 * idx + 41);

            auto pc_xxyzz = pcDistances.data(55 * idx + 42);

            auto pc_xxzzz = pcDistances.data(55 * idx + 43);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_x_xxzz = primBuffer.data(45 * idx + 5);

            auto t_x_xyyy = primBuffer.data(45 * idx + 6);

            auto t_x_xyyz = primBuffer.data(45 * idx + 7);

            auto t_x_xyzz = primBuffer.data(45 * idx + 8);

            auto t_x_xzzz = primBuffer.data(45 * idx + 9);

            // Batch of Integrals (5,10)

            #pragma omp simd aligned(fx, pa_x, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, pc_x, pc_xx, pc_xxx, pc_xxxz, pc_xxxzz, pc_xxy, pc_xxyy, pc_xxyyy, pc_xxyyz, \
                                     pc_xxyz, pc_xxyzz, pc_xxz, pc_xxzz, pc_xxzzz, pc_xy, pc_xyy, pc_xyyy, pc_xyyz, pc_xyz, \
                                     pc_xyzz, pc_xz, pc_xzz, pc_xzzz, pc_y, pc_yy, pc_yyy, pc_yyz, pc_yz, pc_yzz, pc_z, pc_zz, \
                                     pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_x_xxzz, t_x_xyyy, \
                                     t_x_xyyz, t_x_xyzz, t_x_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_x_xxzz[j] = fl_s_0_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * fl2_fx * pb_x[j] + 0.5 * pa_x[j] * pb_xx[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_zz[j] + fl1_fx * pb_xzz[j] + pa_x[j] * pb_xxzz[j]);

                t_x_xxzz[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - fl2_fx * pb_x[j] - 0.5 * pa_x[j] * pb_xx[j] * fl1_fx - pa_x[j] * pb_x[j] * pc_x[j] * fl1_fx - pa_x[j] * fl1_fx * pb_z[j] * pc_z[j] - 0.5 * pa_x[j] * fl1_fx * pb_zz[j] - 0.5 * pc_x[j] * pb_xx[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pb_zz[j] - 2.0 * fl1_fx * pb_xz[j] * pc_z[j] - fl1_fx * pb_xzz[j] - 2.0 * pa_x[j] * pb_xxz[j] * pc_z[j] - 2.0 * pa_x[j] * pb_xzz[j] * pc_x[j] - pc_x[j] * pb_xxzz[j]);

                t_x_xxzz[j] += fl_s_0_0_2 * (0.25 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + 0.5 * fl2_fx * pb_x[j] + pa_x[j] * pb_x[j] * pc_x[j] * fl1_fx + 0.5 * pa_x[j] * pc_xx[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pc_zz[j] + pa_x[j] * fl1_fx * pb_z[j] * pc_z[j] + 0.5 * pc_x[j] * pb_xx[j] * fl1_fx + pc_xx[j] * pb_x[j] * fl1_fx + 3.0 * pc_xz[j] * fl1_fx * pb_z[j] + 1.5 * pc_x[j] * fl1_fx * pb_zz[j] + fl1_fx * pb_x[j] * pc_zz[j] + 2.0 * fl1_fx * pb_xz[j] * pc_z[j] + pa_x[j] * pb_xx[j] * pc_zz[j] + 4.0 * pa_x[j] * pb_xz[j] * pc_xz[j] + pa_x[j] * pc_xx[j] * pb_zz[j] + 2.0 * pc_xz[j] * pb_xxz[j] + 2.0 * pc_xx[j] * pb_xzz[j]);

                t_x_xxzz[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - 0.5 * pa_x[j] * pc_xx[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fx * pc_zz[j] - pc_xx[j] * pb_x[j] * fl1_fx - 0.5 * pc_xxx[j] * fl1_fx - 1.5 * pc_xzz[j] * fl1_fx - 3.0 * pc_xz[j] * fl1_fx * pb_z[j] - fl1_fx * pb_x[j] * pc_zz[j] - 2.0 * pa_x[j] * pb_x[j] * pc_xzz[j] - 2.0 * pa_x[j] * pc_xxz[j] * pb_z[j] - pc_xzz[j] * pb_xx[j] - 4.0 * pc_xxz[j] * pb_xz[j] - pc_xxx[j] * pb_zz[j]);

                t_x_xxzz[j] += fl_s_0_0_4 * (0.5 * pc_xxx[j] * fl1_fx + 1.5 * pc_xzz[j] * fl1_fx + pa_x[j] * pc_xxzz[j] + 2.0 * pc_xxzz[j] * pb_x[j] + 2.0 * pc_xxxz[j] * pb_z[j]);

                t_x_xxzz[j] += -fl_s_0_0_5 * pc_xxxzz[j];

                t_x_xyyy[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_y[j] + 1.5 * pa_x[j] * pb_xy[j] * fl1_fx + 0.5 * fl1_fx * pb_yyy[j] + pa_x[j] * pb_xyyy[j]);

                t_x_xyyy[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_y[j] - 0.75 * fl2_fx * pc_y[j] - 1.5 * pa_x[j] * pb_xy[j] * fl1_fx - 1.5 * pa_x[j] * pb_x[j] * pc_y[j] * fl1_fx - 1.5 * pa_x[j] * pc_x[j] * pb_y[j] * fl1_fx - 1.5 * pc_x[j] * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * pb_yy[j] * pc_y[j] - 0.5 * fl1_fx * pb_yyy[j] - 3.0 * pa_x[j] * pb_xyy[j] * pc_y[j] - pa_x[j] * pc_x[j] * pb_yyy[j] - pc_x[j] * pb_xyyy[j]);

                t_x_xyyy[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.75 * fl2_fx * pb_y[j] + 1.5 * pa_x[j] * pb_x[j] * pc_y[j] * fl1_fx + 1.5 * pa_x[j] * pc_x[j] * pb_y[j] * fl1_fx + 1.5 * pa_x[j] * pc_xy[j] * fl1_fx + 1.5 * pc_x[j] * pb_xy[j] * fl1_fx + 1.5 * pc_xy[j] * pb_x[j] * fl1_fx + 1.5 * pc_xx[j] * pb_y[j] * fl1_fx + 1.5 * fl1_fx * pb_y[j] * pc_yy[j] + 1.5 * fl1_fx * pb_yy[j] * pc_y[j] + 3.0 * pa_x[j] * pb_xy[j] * pc_yy[j] + 3.0 * pa_x[j] * pc_xy[j] * pb_yy[j] + 3.0 * pc_xy[j] * pb_xyy[j] + pc_xx[j] * pb_yyy[j]);

                t_x_xyyy[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 1.5 * pa_x[j] * pc_xy[j] * fl1_fx - 1.5 * pc_xy[j] * pb_x[j] * fl1_fx - 1.5 * pc_xx[j] * pb_y[j] * fl1_fx - 1.5 * pc_xxy[j] * fl1_fx - 0.5 * fl1_fx * pc_yyy[j] - 1.5 * fl1_fx * pb_y[j] * pc_yy[j] - pa_x[j] * pb_x[j] * pc_yyy[j] - 3.0 * pa_x[j] * pc_xyy[j] * pb_y[j] - 3.0 * pc_xyy[j] * pb_xy[j] - 3.0 * pc_xxy[j] * pb_yy[j]);

                t_x_xyyy[j] += fl_s_0_0_4 * (1.5 * pc_xxy[j] * fl1_fx + 0.5 * fl1_fx * pc_yyy[j] + pa_x[j] * pc_xyyy[j] + pc_xyyy[j] * pb_x[j] + 3.0 * pc_xxyy[j] * pb_y[j]);

                t_x_xyyy[j] += -fl_s_0_0_5 * pc_xxyyy[j];

                t_x_xyyz[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_z[j] + 0.5 * pa_x[j] * pb_xz[j] * fl1_fx + 0.5 * fl1_fx * pb_yyz[j] + pa_x[j] * pb_xyyz[j]);

                t_x_xyyz[j] += fl_s_0_0_1 * (-0.25 * fl2_fx * pc_z[j] - 0.5 * fl2_fx * pb_z[j] - 0.5 * pa_x[j] * pb_x[j] * fl1_fx * pc_z[j] - 0.5 * pa_x[j] * pb_xz[j] * fl1_fx - 0.5 * pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] - 0.5 * pc_x[j] * pb_xz[j] * fl1_fx - 0.5 * fl1_fx * pb_yy[j] * pc_z[j] - fl1_fx * pb_yz[j] * pc_y[j] - 0.5 * fl1_fx * pb_yyz[j] - pa_x[j] * pb_xyy[j] * pc_z[j] - 2.0 * pa_x[j] * pb_xyz[j] * pc_y[j] - pa_x[j] * pc_x[j] * pb_yyz[j] - pc_x[j] * pb_xyyz[j]);

                t_x_xyyz[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_z[j] + 0.25 * fl2_fx * pb_z[j] + 0.5 * pa_x[j] * pb_x[j] * fl1_fx * pc_z[j] + 0.5 * pa_x[j] * pc_xz[j] * fl1_fx + 0.5 * pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] + 0.5 * pc_xz[j] * pb_x[j] * fl1_fx + 0.5 * pc_x[j] * pb_xz[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx * pb_z[j] + fl1_fx * pb_y[j] * pc_yz[j] + 0.5 * fl1_fx * pc_yy[j] * pb_z[j] + 0.5 * fl1_fx * pb_yy[j] * pc_z[j] + fl1_fx * pb_yz[j] * pc_y[j] + 2.0 * pa_x[j] * pb_xy[j] * pc_yz[j] + pa_x[j] * pb_xz[j] * pc_yy[j] + pa_x[j] * pc_xz[j] * pb_yy[j] + 2.0 * pa_x[j] * pc_xy[j] * pb_yz[j] + pc_xz[j] * pb_xyy[j] + 2.0 * pc_xy[j] * pb_xyz[j] + pc_xx[j] * pb_yyz[j]);

                t_x_xyyz[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_z[j] - 0.5 * pa_x[j] * pc_xz[j] * fl1_fx - 0.5 * pc_xz[j] * pb_x[j] * fl1_fx - 0.5 * pc_xxz[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_yyz[j] - fl1_fx * pb_y[j] * pc_yz[j] - 0.5 * fl1_fx * pc_yy[j] * pb_z[j] - pa_x[j] * pb_x[j] * pc_yyz[j] - 2.0 * pa_x[j] * pc_xyz[j] * pb_y[j] - pa_x[j] * pc_xyy[j] * pb_z[j] - 2.0 * pc_xyz[j] * pb_xy[j] - pc_xyy[j] * pb_xz[j] - pc_xxz[j] * pb_yy[j] - 2.0 * pc_xxy[j] * pb_yz[j]);

                t_x_xyyz[j] += fl_s_0_0_4 * (0.5 * pc_xxz[j] * fl1_fx + 0.5 * fl1_fx * pc_yyz[j] + pa_x[j] * pc_xyyz[j] + pc_xyyz[j] * pb_x[j] + 2.0 * pc_xxyz[j] * pb_y[j] + pc_xxyy[j] * pb_z[j]);

                t_x_xyyz[j] += -fl_s_0_0_5 * pc_xxyyz[j];

                t_x_xyzz[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_y[j] + 0.5 * pa_x[j] * pb_xy[j] * fl1_fx + 0.5 * fl1_fx * pb_yzz[j] + pa_x[j] * pb_xyzz[j]);

                t_x_xyzz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx * pb_y[j] - 0.25 * fl2_fx * pc_y[j] - 0.5 * pa_x[j] * pb_xy[j] * fl1_fx - 0.5 * pa_x[j] * pb_x[j] * pc_y[j] * fl1_fx - 0.5 * pa_x[j] * pc_x[j] * pb_y[j] * fl1_fx - 0.5 * pc_x[j] * pb_xy[j] * fl1_fx - fl1_fx * pb_yz[j] * pc_z[j] - 0.5 * fl1_fx * pc_y[j] * pb_zz[j] - 0.5 * fl1_fx * pb_yzz[j] - 2.0 * pa_x[j] * pb_xyz[j] * pc_z[j] - pa_x[j] * pb_xzz[j] * pc_y[j] - pa_x[j] * pc_x[j] * pb_yzz[j] - pc_x[j] * pb_xyzz[j]);

                t_x_xyzz[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_y[j] + 0.25 * fl2_fx * pb_y[j] + 0.5 * pa_x[j] * pb_x[j] * pc_y[j] * fl1_fx + 0.5 * pa_x[j] * pc_x[j] * pb_y[j] * fl1_fx + 0.5 * pa_x[j] * pc_xy[j] * fl1_fx + 0.5 * pc_x[j] * pb_xy[j] * fl1_fx + 0.5 * pc_xy[j] * pb_x[j] * fl1_fx + 0.5 * pc_xx[j] * pb_y[j] * fl1_fx + 0.5 * fl1_fx * pb_y[j] * pc_zz[j] + fl1_fx * pc_yz[j] * pb_z[j] + fl1_fx * pb_yz[j] * pc_z[j] + 0.5 * fl1_fx * pc_y[j] * pb_zz[j] + pa_x[j] * pb_xy[j] * pc_zz[j] + 2.0 * pa_x[j] * pb_xz[j] * pc_yz[j] + 2.0 * pa_x[j] * pc_xz[j] * pb_yz[j] + pa_x[j] * pc_xy[j] * pb_zz[j] + 2.0 * pc_xz[j] * pb_xyz[j] + pc_xy[j] * pb_xzz[j] + pc_xx[j] * pb_yzz[j]);

                t_x_xyzz[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_y[j] - 0.5 * pa_x[j] * pc_xy[j] * fl1_fx - 0.5 * pc_xy[j] * pb_x[j] * fl1_fx - 0.5 * pc_xx[j] * pb_y[j] * fl1_fx - 0.5 * pc_xxy[j] * fl1_fx - 0.5 * fl1_fx * pc_yzz[j] - 0.5 * fl1_fx * pb_y[j] * pc_zz[j] - fl1_fx * pc_yz[j] * pb_z[j] - pa_x[j] * pb_x[j] * pc_yzz[j] - pa_x[j] * pc_xzz[j] * pb_y[j] - 2.0 * pa_x[j] * pc_xyz[j] * pb_z[j] - pc_xzz[j] * pb_xy[j] - 2.0 * pc_xyz[j] * pb_xz[j] - 2.0 * pc_xxz[j] * pb_yz[j] - pc_xxy[j] * pb_zz[j]);

                t_x_xyzz[j] += fl_s_0_0_4 * (0.5 * pc_xxy[j] * fl1_fx + 0.5 * fl1_fx * pc_yzz[j] + pa_x[j] * pc_xyzz[j] + pc_xyzz[j] * pb_x[j] + pc_xxzz[j] * pb_y[j] + 2.0 * pc_xxyz[j] * pb_z[j]);

                t_x_xyzz[j] += -fl_s_0_0_5 * pc_xxyzz[j];

                t_x_xzzz[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_z[j] + 1.5 * pa_x[j] * pb_xz[j] * fl1_fx + 0.5 * fl1_fx * pb_zzz[j] + pa_x[j] * pb_xzzz[j]);

                t_x_xzzz[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_z[j] - 0.75 * fl2_fx * pc_z[j] - 1.5 * pa_x[j] * pb_xz[j] * fl1_fx - 1.5 * pa_x[j] * pb_x[j] * pc_z[j] * fl1_fx - 1.5 * pa_x[j] * pc_x[j] * pb_z[j] * fl1_fx - 1.5 * pc_x[j] * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * pb_zz[j] * pc_z[j] - 0.5 * fl1_fx * pb_zzz[j] - 3.0 * pa_x[j] * pb_xzz[j] * pc_z[j] - pa_x[j] * pc_x[j] * pb_zzz[j] - pc_x[j] * pb_xzzz[j]);

                t_x_xzzz[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pb_z[j] + 1.5 * pa_x[j] * pb_x[j] * pc_z[j] * fl1_fx + 1.5 * pa_x[j] * pc_x[j] * pb_z[j] * fl1_fx + 1.5 * pa_x[j] * pc_xz[j] * fl1_fx + 1.5 * pc_x[j] * pb_xz[j] * fl1_fx + 1.5 * pc_xz[j] * pb_x[j] * fl1_fx + 1.5 * pc_xx[j] * pb_z[j] * fl1_fx + 1.5 * fl1_fx * pb_z[j] * pc_zz[j] + 1.5 * fl1_fx * pb_zz[j] * pc_z[j] + 3.0 * pa_x[j] * pb_xz[j] * pc_zz[j] + 3.0 * pa_x[j] * pc_xz[j] * pb_zz[j] + 3.0 * pc_xz[j] * pb_xzz[j] + pc_xx[j] * pb_zzz[j]);

                t_x_xzzz[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * pa_x[j] * pc_xz[j] * fl1_fx - 1.5 * pc_xz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xx[j] * pb_z[j] * fl1_fx - 1.5 * pc_xxz[j] * fl1_fx - 0.5 * fl1_fx * pc_zzz[j] - 1.5 * fl1_fx * pb_z[j] * pc_zz[j] - pa_x[j] * pb_x[j] * pc_zzz[j] - 3.0 * pa_x[j] * pc_xzz[j] * pb_z[j] - 3.0 * pc_xzz[j] * pb_xz[j] - 3.0 * pc_xxz[j] * pb_zz[j]);

                t_x_xzzz[j] += fl_s_0_0_4 * (1.5 * pc_xxz[j] * fl1_fx + 0.5 * fl1_fx * pc_zzz[j] + pa_x[j] * pc_xzzz[j] + pc_xzzz[j] * pb_x[j] + 3.0 * pc_xxzz[j] * pb_z[j]);

                t_x_xzzz[j] += -fl_s_0_0_5 * pc_xxzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG_10_15(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
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

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

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

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            auto pc_xzzz = pcDistances.data(55 * idx + 28);

            auto pc_yyyy = pcDistances.data(55 * idx + 29);

            auto pc_yyyz = pcDistances.data(55 * idx + 30);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            auto pc_yzzz = pcDistances.data(55 * idx + 32);

            auto pc_zzzz = pcDistances.data(55 * idx + 33);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xyyyy = pcDistances.data(55 * idx + 44);

            auto pc_xyyyz = pcDistances.data(55 * idx + 45);

            auto pc_xyyzz = pcDistances.data(55 * idx + 46);

            auto pc_xyzzz = pcDistances.data(55 * idx + 47);

            auto pc_xzzzz = pcDistances.data(55 * idx + 48);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_x_yyyy = primBuffer.data(45 * idx + 10);

            auto t_x_yyyz = primBuffer.data(45 * idx + 11);

            auto t_x_yyzz = primBuffer.data(45 * idx + 12);

            auto t_x_yzzz = primBuffer.data(45 * idx + 13);

            auto t_x_zzzz = primBuffer.data(45 * idx + 14);

            // Batch of Integrals (10,15)

            #pragma omp simd aligned(fx, pa_x, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, pc_x, pc_xy, pc_xyy, pc_xyyy, pc_xyyyy, pc_xyyyz, \
                                     pc_xyyz, pc_xyyzz, pc_xyz, pc_xyzz, pc_xyzzz, pc_xz, pc_xzz, pc_xzzz, pc_xzzzz, pc_y, \
                                     pc_yy, pc_yyy, pc_yyyy, pc_yyyz, pc_yyz, pc_yyzz, pc_yz, pc_yzz, pc_yzzz, pc_z, pc_zz, \
                                     pc_zzz, pc_zzzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_x_yyyy, \
                                     t_x_yyyz, t_x_yyzz, t_x_yzzz, t_x_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_x_yyyy[j] = fl_s_0_0_0 * (0.75 * pa_x[j] * fl2_fx + 3.0 * pa_x[j] * pb_yy[j] * fl1_fx + pa_x[j] * pb_yyyy[j]);

                t_x_yyyy[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - 3.0 * pa_x[j] * pb_yy[j] * fl1_fx - 6.0 * pa_x[j] * pb_y[j] * pc_y[j] * fl1_fx - 3.0 * pc_x[j] * pb_yy[j] * fl1_fx - 4.0 * pa_x[j] * pb_yyy[j] * pc_y[j] - pc_x[j] * pb_yyyy[j]);

                t_x_yyyy[j] += fl_s_0_0_2 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + 6.0 * pa_x[j] * pb_y[j] * pc_y[j] * fl1_fx + 3.0 * pa_x[j] * pc_yy[j] * fl1_fx + 3.0 * pc_x[j] * pb_yy[j] * fl1_fx + 6.0 * pc_xy[j] * pb_y[j] * fl1_fx + 6.0 * pa_x[j] * pb_yy[j] * pc_yy[j] + 4.0 * pc_xy[j] * pb_yyy[j]);

                t_x_yyyy[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - 3.0 * pa_x[j] * pc_yy[j] * fl1_fx - 6.0 * pc_xy[j] * pb_y[j] * fl1_fx - 3.0 * pc_xyy[j] * fl1_fx - 4.0 * pa_x[j] * pb_y[j] * pc_yyy[j] - 6.0 * pc_xyy[j] * pb_yy[j]);

                t_x_yyyy[j] += fl_s_0_0_4 * (3.0 * pc_xyy[j] * fl1_fx + pa_x[j] * pc_yyyy[j] + 4.0 * pc_xyyy[j] * pb_y[j]);

                t_x_yyyy[j] += -fl_s_0_0_5 * pc_xyyyy[j];

                t_x_yyyz[j] = fl_s_0_0_0 * (1.5 * pa_x[j] * pb_yz[j] * fl1_fx + pa_x[j] * pb_yyyz[j]);

                t_x_yyyz[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * pb_y[j] * fl1_fx * pc_z[j] - 1.5 * pa_x[j] * pb_yz[j] * fl1_fx - 1.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_z[j] - 1.5 * pc_x[j] * pb_yz[j] * fl1_fx - pa_x[j] * pb_yyy[j] * pc_z[j] - 3.0 * pa_x[j] * pb_yyz[j] * pc_y[j] - pc_x[j] * pb_yyyz[j]);

                t_x_yyyz[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * pb_y[j] * fl1_fx * pc_z[j] + 1.5 * pa_x[j] * pc_yz[j] * fl1_fx + 1.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_z[j] + 1.5 * pc_xz[j] * pb_y[j] * fl1_fx + 1.5 * pc_x[j] * pb_yz[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx * pb_z[j] + 3.0 * pa_x[j] * pb_yy[j] * pc_yz[j] + 3.0 * pa_x[j] * pb_yz[j] * pc_yy[j] + pc_xz[j] * pb_yyy[j] + 3.0 * pc_xy[j] * pb_yyz[j]);

                t_x_yyyz[j] += fl_s_0_0_3 * (-1.5 * pa_x[j] * pc_yz[j] * fl1_fx - 1.5 * pc_xz[j] * pb_y[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pb_z[j] - 3.0 * pa_x[j] * pb_y[j] * pc_yyz[j] - pa_x[j] * pc_yyy[j] * pb_z[j] - 3.0 * pc_xyz[j] * pb_yy[j] - 3.0 * pc_xyy[j] * pb_yz[j]);

                t_x_yyyz[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_x[j] * pc_yyyz[j] + 3.0 * pc_xyyz[j] * pb_y[j] + pc_xyyy[j] * pb_z[j]);

                t_x_yyyz[j] += -fl_s_0_0_5 * pc_xyyyz[j];

                t_x_yyzz[j] = fl_s_0_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa_x[j] * pb_yy[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_zz[j] + pa_x[j] * pb_yyzz[j]);

                t_x_yyzz[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl2_fx - 0.25 * pc_x[j] * fl2_fx - 0.5 * pa_x[j] * pb_yy[j] * fl1_fx - pa_x[j] * pb_y[j] * pc_y[j] * fl1_fx - pa_x[j] * fl1_fx * pb_z[j] * pc_z[j] - 0.5 * pa_x[j] * fl1_fx * pb_zz[j] - 0.5 * pc_x[j] * pb_yy[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx * pb_zz[j] - 2.0 * pa_x[j] * pb_yyz[j] * pc_z[j] - 2.0 * pa_x[j] * pb_yzz[j] * pc_y[j] - pc_x[j] * pb_yyzz[j]);

                t_x_yyzz[j] += fl_s_0_0_2 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pc_x[j] * fl2_fx + pa_x[j] * pb_y[j] * pc_y[j] * fl1_fx + 0.5 * pa_x[j] * pc_yy[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pc_zz[j] + pa_x[j] * fl1_fx * pb_z[j] * pc_z[j] + 0.5 * pc_x[j] * pb_yy[j] * fl1_fx + pc_xy[j] * pb_y[j] * fl1_fx + pc_xz[j] * fl1_fx * pb_z[j] + 0.5 * pc_x[j] * fl1_fx * pb_zz[j] + pa_x[j] * pb_yy[j] * pc_zz[j] + 4.0 * pa_x[j] * pb_yz[j] * pc_yz[j] + pa_x[j] * pc_yy[j] * pb_zz[j] + 2.0 * pc_xz[j] * pb_yyz[j] + 2.0 * pc_xy[j] * pb_yzz[j]);

                t_x_yyzz[j] += fl_s_0_0_3 * (-0.25 * pc_x[j] * fl2_fx - 0.5 * pa_x[j] * pc_yy[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fx * pc_zz[j] - pc_xy[j] * pb_y[j] * fl1_fx - 0.5 * pc_xyy[j] * fl1_fx - 0.5 * pc_xzz[j] * fl1_fx - pc_xz[j] * fl1_fx * pb_z[j] - 2.0 * pa_x[j] * pb_y[j] * pc_yzz[j] - 2.0 * pa_x[j] * pc_yyz[j] * pb_z[j] - pc_xzz[j] * pb_yy[j] - 4.0 * pc_xyz[j] * pb_yz[j] - pc_xyy[j] * pb_zz[j]);

                t_x_yyzz[j] += fl_s_0_0_4 * (0.5 * pc_xyy[j] * fl1_fx + 0.5 * pc_xzz[j] * fl1_fx + pa_x[j] * pc_yyzz[j] + 2.0 * pc_xyzz[j] * pb_y[j] + 2.0 * pc_xyyz[j] * pb_z[j]);

                t_x_yyzz[j] += -fl_s_0_0_5 * pc_xyyzz[j];

                t_x_yzzz[j] = fl_s_0_0_0 * (1.5 * pa_x[j] * pb_yz[j] * fl1_fx + pa_x[j] * pb_yzzz[j]);

                t_x_yzzz[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * pb_yz[j] * fl1_fx - 1.5 * pa_x[j] * pb_y[j] * pc_z[j] * fl1_fx - 1.5 * pa_x[j] * pc_y[j] * pb_z[j] * fl1_fx - 1.5 * pc_x[j] * pb_yz[j] * fl1_fx - 3.0 * pa_x[j] * pb_yzz[j] * pc_z[j] - pa_x[j] * pc_y[j] * pb_zzz[j] - pc_x[j] * pb_yzzz[j]);

                t_x_yzzz[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * pb_y[j] * pc_z[j] * fl1_fx + 1.5 * pa_x[j] * pc_y[j] * pb_z[j] * fl1_fx + 1.5 * pa_x[j] * pc_yz[j] * fl1_fx + 1.5 * pc_x[j] * pb_yz[j] * fl1_fx + 1.5 * pc_xz[j] * pb_y[j] * fl1_fx + 1.5 * pc_xy[j] * pb_z[j] * fl1_fx + 3.0 * pa_x[j] * pb_yz[j] * pc_zz[j] + 3.0 * pa_x[j] * pc_yz[j] * pb_zz[j] + 3.0 * pc_xz[j] * pb_yzz[j] + pc_xy[j] * pb_zzz[j]);

                t_x_yzzz[j] += fl_s_0_0_3 * (-1.5 * pa_x[j] * pc_yz[j] * fl1_fx - 1.5 * pc_xz[j] * pb_y[j] * fl1_fx - 1.5 * pc_xy[j] * pb_z[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - pa_x[j] * pb_y[j] * pc_zzz[j] - 3.0 * pa_x[j] * pc_yzz[j] * pb_z[j] - 3.0 * pc_xzz[j] * pb_yz[j] - 3.0 * pc_xyz[j] * pb_zz[j]);

                t_x_yzzz[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_x[j] * pc_yzzz[j] + pc_xzzz[j] * pb_y[j] + 3.0 * pc_xyzz[j] * pb_z[j]);

                t_x_yzzz[j] += -fl_s_0_0_5 * pc_xyzzz[j];

                t_x_zzzz[j] = fl_s_0_0_0 * (0.75 * pa_x[j] * fl2_fx + 3.0 * pa_x[j] * pb_zz[j] * fl1_fx + pa_x[j] * pb_zzzz[j]);

                t_x_zzzz[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - 3.0 * pa_x[j] * pb_zz[j] * fl1_fx - 6.0 * pa_x[j] * pb_z[j] * pc_z[j] * fl1_fx - 3.0 * pc_x[j] * pb_zz[j] * fl1_fx - 4.0 * pa_x[j] * pb_zzz[j] * pc_z[j] - pc_x[j] * pb_zzzz[j]);

                t_x_zzzz[j] += fl_s_0_0_2 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + 6.0 * pa_x[j] * pb_z[j] * pc_z[j] * fl1_fx + 3.0 * pa_x[j] * pc_zz[j] * fl1_fx + 3.0 * pc_x[j] * pb_zz[j] * fl1_fx + 6.0 * pc_xz[j] * pb_z[j] * fl1_fx + 6.0 * pa_x[j] * pb_zz[j] * pc_zz[j] + 4.0 * pc_xz[j] * pb_zzz[j]);

                t_x_zzzz[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - 3.0 * pa_x[j] * pc_zz[j] * fl1_fx - 6.0 * pc_xz[j] * pb_z[j] * fl1_fx - 3.0 * pc_xzz[j] * fl1_fx - 4.0 * pa_x[j] * pb_z[j] * pc_zzz[j] - 6.0 * pc_xzz[j] * pb_zz[j]);

                t_x_zzzz[j] += fl_s_0_0_4 * (3.0 * pc_xzz[j] * fl1_fx + pa_x[j] * pc_zzzz[j] + 4.0 * pc_xzzz[j] * pb_z[j]);

                t_x_zzzz[j] += -fl_s_0_0_5 * pc_xzzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG_15_20(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (15,20)

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

            auto pa_y = paDistances.data(3 * idx + 1);

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

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(55 * idx + 9);

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(55 * idx + 19);

            auto pc_xxxy = pcDistances.data(55 * idx + 20);

            auto pc_xxxz = pcDistances.data(55 * idx + 21);

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxxxy = pcDistances.data(55 * idx + 35);

            auto pc_xxxyy = pcDistances.data(55 * idx + 37);

            auto pc_xxxyz = pcDistances.data(55 * idx + 38);

            auto pc_xxyyy = pcDistances.data(55 * idx + 40);

            auto pc_xxyyz = pcDistances.data(55 * idx + 41);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_y_xxxx = primBuffer.data(45 * idx + 15);

            auto t_y_xxxy = primBuffer.data(45 * idx + 16);

            auto t_y_xxxz = primBuffer.data(45 * idx + 17);

            auto t_y_xxyy = primBuffer.data(45 * idx + 18);

            auto t_y_xxyz = primBuffer.data(45 * idx + 19);

            // Batch of Integrals (15,20)

            #pragma omp simd aligned(fx, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pc_x, pc_xx, \
                                     pc_xxx, pc_xxxx, pc_xxxxy, pc_xxxy, pc_xxxyy, pc_xxxyz, pc_xxxz, pc_xxy, pc_xxyy, \
                                     pc_xxyyy, pc_xxyyz, pc_xxyz, pc_xxz, pc_xy, pc_xyy, pc_xyyy, pc_xyyz, pc_xyz, pc_xz, \
                                     pc_y, pc_yy, pc_yyy, pc_yyz, pc_yz, pc_z, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, \
                                     s_0_0_4, s_0_0_5, t_y_xxxx, t_y_xxxy, t_y_xxxz, t_y_xxyy, t_y_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_y_xxxx[j] = fl_s_0_0_0 * (0.75 * pa_y[j] * fl2_fx + 3.0 * pa_y[j] * pb_xx[j] * fl1_fx + pa_y[j] * pb_xxxx[j]);

                t_y_xxxx[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl2_fx - 0.75 * pc_y[j] * fl2_fx - 3.0 * pa_y[j] * pb_xx[j] * fl1_fx - 6.0 * pa_y[j] * pb_x[j] * pc_x[j] * fl1_fx - 3.0 * pc_y[j] * pb_xx[j] * fl1_fx - 4.0 * pa_y[j] * pb_xxx[j] * pc_x[j] - pc_y[j] * pb_xxxx[j]);

                t_y_xxxx[j] += fl_s_0_0_2 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pc_y[j] * fl2_fx + 6.0 * pa_y[j] * pb_x[j] * pc_x[j] * fl1_fx + 3.0 * pa_y[j] * pc_xx[j] * fl1_fx + 3.0 * pc_y[j] * pb_xx[j] * fl1_fx + 6.0 * pc_xy[j] * pb_x[j] * fl1_fx + 6.0 * pa_y[j] * pb_xx[j] * pc_xx[j] + 4.0 * pc_xy[j] * pb_xxx[j]);

                t_y_xxxx[j] += fl_s_0_0_3 * (-0.75 * pc_y[j] * fl2_fx - 3.0 * pa_y[j] * pc_xx[j] * fl1_fx - 6.0 * pc_xy[j] * pb_x[j] * fl1_fx - 3.0 * pc_xxy[j] * fl1_fx - 4.0 * pa_y[j] * pb_x[j] * pc_xxx[j] - 6.0 * pc_xxy[j] * pb_xx[j]);

                t_y_xxxx[j] += fl_s_0_0_4 * (3.0 * pc_xxy[j] * fl1_fx + pa_y[j] * pc_xxxx[j] + 4.0 * pc_xxxy[j] * pb_x[j]);

                t_y_xxxx[j] += -fl_s_0_0_5 * pc_xxxxy[j];

                t_y_xxxy[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_x[j] + 1.5 * pa_y[j] * pb_xy[j] * fl1_fx + 0.5 * fl1_fx * pb_xxx[j] + pa_y[j] * pb_xxxy[j]);

                t_y_xxxy[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_x[j] - 0.75 * fl2_fx * pc_x[j] - 1.5 * pa_y[j] * pb_x[j] * fl1_fx * pc_y[j] - 1.5 * pa_y[j] * pb_xy[j] * fl1_fx - 1.5 * pa_y[j] * pc_x[j] * fl1_fx * pb_y[j] - 1.5 * pc_y[j] * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * pb_xx[j] * pc_x[j] - 0.5 * fl1_fx * pb_xxx[j] - pa_y[j] * pb_xxx[j] * pc_y[j] - 3.0 * pa_y[j] * pb_xxy[j] * pc_x[j] - pc_y[j] * pb_xxxy[j]);

                t_y_xxxy[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_x[j] + 0.75 * fl2_fx * pb_x[j] + 1.5 * pa_y[j] * pb_x[j] * fl1_fx * pc_y[j] + 1.5 * pa_y[j] * pc_xy[j] * fl1_fx + 1.5 * pa_y[j] * pc_x[j] * fl1_fx * pb_y[j] + 1.5 * pc_yy[j] * pb_x[j] * fl1_fx + 1.5 * pc_y[j] * pb_xy[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx * pb_y[j] + 1.5 * fl1_fx * pb_x[j] * pc_xx[j] + 1.5 * fl1_fx * pb_xx[j] * pc_x[j] + 3.0 * pa_y[j] * pb_xx[j] * pc_xy[j] + 3.0 * pa_y[j] * pb_xy[j] * pc_xx[j] + pc_yy[j] * pb_xxx[j] + 3.0 * pc_xy[j] * pb_xxy[j]);

                t_y_xxxy[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_x[j] - 1.5 * pa_y[j] * pc_xy[j] * fl1_fx - 1.5 * pc_yy[j] * pb_x[j] * fl1_fx - 1.5 * pc_xyy[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pb_y[j] - 0.5 * fl1_fx * pc_xxx[j] - 1.5 * fl1_fx * pb_x[j] * pc_xx[j] - 3.0 * pa_y[j] * pb_x[j] * pc_xxy[j] - pa_y[j] * pc_xxx[j] * pb_y[j] - 3.0 * pc_xyy[j] * pb_xx[j] - 3.0 * pc_xxy[j] * pb_xy[j]);

                t_y_xxxy[j] += fl_s_0_0_4 * (1.5 * pc_xyy[j] * fl1_fx + 0.5 * fl1_fx * pc_xxx[j] + pa_y[j] * pc_xxxy[j] + 3.0 * pc_xxyy[j] * pb_x[j] + pc_xxxy[j] * pb_y[j]);

                t_y_xxxy[j] += -fl_s_0_0_5 * pc_xxxyy[j];

                t_y_xxxz[j] = fl_s_0_0_0 * (1.5 * pa_y[j] * pb_xz[j] * fl1_fx + pa_y[j] * pb_xxxz[j]);

                t_y_xxxz[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * pb_x[j] * fl1_fx * pc_z[j] - 1.5 * pa_y[j] * pb_xz[j] * fl1_fx - 1.5 * pa_y[j] * pc_x[j] * fl1_fx * pb_z[j] - 1.5 * pc_y[j] * pb_xz[j] * fl1_fx - pa_y[j] * pb_xxx[j] * pc_z[j] - 3.0 * pa_y[j] * pb_xxz[j] * pc_x[j] - pc_y[j] * pb_xxxz[j]);

                t_y_xxxz[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * pb_x[j] * fl1_fx * pc_z[j] + 1.5 * pa_y[j] * pc_xz[j] * fl1_fx + 1.5 * pa_y[j] * pc_x[j] * fl1_fx * pb_z[j] + 1.5 * pc_yz[j] * pb_x[j] * fl1_fx + 1.5 * pc_y[j] * pb_xz[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx * pb_z[j] + 3.0 * pa_y[j] * pb_xx[j] * pc_xz[j] + 3.0 * pa_y[j] * pb_xz[j] * pc_xx[j] + pc_yz[j] * pb_xxx[j] + 3.0 * pc_xy[j] * pb_xxz[j]);

                t_y_xxxz[j] += fl_s_0_0_3 * (-1.5 * pa_y[j] * pc_xz[j] * fl1_fx - 1.5 * pc_yz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pb_z[j] - 3.0 * pa_y[j] * pb_x[j] * pc_xxz[j] - pa_y[j] * pc_xxx[j] * pb_z[j] - 3.0 * pc_xyz[j] * pb_xx[j] - 3.0 * pc_xxy[j] * pb_xz[j]);

                t_y_xxxz[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_y[j] * pc_xxxz[j] + 3.0 * pc_xxyz[j] * pb_x[j] + pc_xxxy[j] * pb_z[j]);

                t_y_xxxz[j] += -fl_s_0_0_5 * pc_xxxyz[j];

                t_y_xxyy[j] = fl_s_0_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * fl2_fx * pb_y[j] + 0.5 * pa_y[j] * pb_xx[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pb_yy[j] + fl1_fx * pb_xxy[j] + pa_y[j] * pb_xxyy[j]);

                t_y_xxyy[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl2_fx - 0.75 * pc_y[j] * fl2_fx - fl2_fx * pb_y[j] - 0.5 * pa_y[j] * pb_xx[j] * fl1_fx - pa_y[j] * pb_x[j] * pc_x[j] * fl1_fx - pa_y[j] * fl1_fx * pb_y[j] * pc_y[j] - 0.5 * pa_y[j] * fl1_fx * pb_yy[j] - 1.5 * pc_y[j] * pb_xx[j] * fl1_fx - 0.5 * pc_y[j] * fl1_fx * pb_yy[j] - 2.0 * fl1_fx * pb_xy[j] * pc_x[j] - fl1_fx * pb_xxy[j] - 2.0 * pa_y[j] * pb_xxy[j] * pc_y[j] - 2.0 * pa_y[j] * pb_xyy[j] * pc_x[j] - pc_y[j] * pb_xxyy[j]);

                t_y_xxyy[j] += fl_s_0_0_2 * (0.25 * pa_y[j] * fl2_fx + 1.5 * pc_y[j] * fl2_fx + 0.5 * fl2_fx * pb_y[j] + pa_y[j] * pb_x[j] * pc_x[j] * fl1_fx + 0.5 * pa_y[j] * pc_xx[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pc_yy[j] + pa_y[j] * fl1_fx * pb_y[j] * pc_y[j] + 1.5 * pc_y[j] * pb_xx[j] * fl1_fx + 3.0 * pc_xy[j] * pb_x[j] * fl1_fx + pc_yy[j] * fl1_fx * pb_y[j] + 0.5 * pc_y[j] * fl1_fx * pb_yy[j] + fl1_fx * pc_xx[j] * pb_y[j] + 2.0 * fl1_fx * pb_xy[j] * pc_x[j] + pa_y[j] * pb_xx[j] * pc_yy[j] + 4.0 * pa_y[j] * pb_xy[j] * pc_xy[j] + pa_y[j] * pc_xx[j] * pb_yy[j] + 2.0 * pc_yy[j] * pb_xxy[j] + 2.0 * pc_xy[j] * pb_xyy[j]);

                t_y_xxyy[j] += fl_s_0_0_3 * (-0.75 * pc_y[j] * fl2_fx - 0.5 * pa_y[j] * pc_xx[j] * fl1_fx - 0.5 * pa_y[j] * fl1_fx * pc_yy[j] - 3.0 * pc_xy[j] * pb_x[j] * fl1_fx - 1.5 * pc_xxy[j] * fl1_fx - 0.5 * pc_yyy[j] * fl1_fx - pc_yy[j] * fl1_fx * pb_y[j] - fl1_fx * pc_xx[j] * pb_y[j] - 2.0 * pa_y[j] * pb_x[j] * pc_xyy[j] - 2.0 * pa_y[j] * pc_xxy[j] * pb_y[j] - pc_yyy[j] * pb_xx[j] - 4.0 * pc_xyy[j] * pb_xy[j] - pc_xxy[j] * pb_yy[j]);

                t_y_xxyy[j] += fl_s_0_0_4 * (1.5 * pc_xxy[j] * fl1_fx + 0.5 * pc_yyy[j] * fl1_fx + pa_y[j] * pc_xxyy[j] + 2.0 * pc_xyyy[j] * pb_x[j] + 2.0 * pc_xxyy[j] * pb_y[j]);

                t_y_xxyy[j] += -fl_s_0_0_5 * pc_xxyyy[j];

                t_y_xxyz[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_z[j] + 0.5 * pa_y[j] * fl1_fx * pb_yz[j] + 0.5 * fl1_fx * pb_xxz[j] + pa_y[j] * pb_xxyz[j]);

                t_y_xxyz[j] += fl_s_0_0_1 * (-0.25 * fl2_fx * pc_z[j] - 0.5 * fl2_fx * pb_z[j] - 0.5 * pa_y[j] * fl1_fx * pb_y[j] * pc_z[j] - 0.5 * pa_y[j] * fl1_fx * pc_y[j] * pb_z[j] - 0.5 * pa_y[j] * fl1_fx * pb_yz[j] - 0.5 * pc_y[j] * fl1_fx * pb_yz[j] - 0.5 * fl1_fx * pb_xx[j] * pc_z[j] - fl1_fx * pb_xz[j] * pc_x[j] - 0.5 * fl1_fx * pb_xxz[j] - pa_y[j] * pb_xxy[j] * pc_z[j] - pa_y[j] * pb_xxz[j] * pc_y[j] - 2.0 * pa_y[j] * pb_xyz[j] * pc_x[j] - pc_y[j] * pb_xxyz[j]);

                t_y_xxyz[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_z[j] + 0.25 * fl2_fx * pb_z[j] + 0.5 * pa_y[j] * fl1_fx * pc_yz[j] + 0.5 * pa_y[j] * fl1_fx * pb_y[j] * pc_z[j] + 0.5 * pa_y[j] * fl1_fx * pc_y[j] * pb_z[j] + 0.5 * pc_yz[j] * fl1_fx * pb_y[j] + 0.5 * pc_yy[j] * fl1_fx * pb_z[j] + 0.5 * pc_y[j] * fl1_fx * pb_yz[j] + fl1_fx * pb_x[j] * pc_xz[j] + 0.5 * fl1_fx * pc_xx[j] * pb_z[j] + 0.5 * fl1_fx * pb_xx[j] * pc_z[j] + fl1_fx * pb_xz[j] * pc_x[j] + pa_y[j] * pb_xx[j] * pc_yz[j] + 2.0 * pa_y[j] * pb_xy[j] * pc_xz[j] + 2.0 * pa_y[j] * pb_xz[j] * pc_xy[j] + pa_y[j] * pc_xx[j] * pb_yz[j] + pc_yz[j] * pb_xxy[j] + pc_yy[j] * pb_xxz[j] + 2.0 * pc_xy[j] * pb_xyz[j]);

                t_y_xxyz[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_z[j] - 0.5 * pa_y[j] * fl1_fx * pc_yz[j] - 0.5 * pc_yyz[j] * fl1_fx - 0.5 * pc_yz[j] * fl1_fx * pb_y[j] - 0.5 * pc_yy[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_xxz[j] - fl1_fx * pb_x[j] * pc_xz[j] - 0.5 * fl1_fx * pc_xx[j] * pb_z[j] - 2.0 * pa_y[j] * pb_x[j] * pc_xyz[j] - pa_y[j] * pc_xxz[j] * pb_y[j] - pa_y[j] * pc_xxy[j] * pb_z[j] - pc_yyz[j] * pb_xx[j] - 2.0 * pc_xyz[j] * pb_xy[j] - 2.0 * pc_xyy[j] * pb_xz[j] - pc_xxy[j] * pb_yz[j]);

                t_y_xxyz[j] += fl_s_0_0_4 * (0.5 * pc_yyz[j] * fl1_fx + 0.5 * fl1_fx * pc_xxz[j] + pa_y[j] * pc_xxyz[j] + 2.0 * pc_xyyz[j] * pb_x[j] + pc_xxyz[j] * pb_y[j] + pc_xxyy[j] * pb_z[j]);

                t_y_xxyz[j] += -fl_s_0_0_5 * pc_xxyyz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG_20_25(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (20,25)

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

            auto pa_y = paDistances.data(3 * idx + 1);

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

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            auto pc_xzzz = pcDistances.data(55 * idx + 28);

            auto pc_yyyy = pcDistances.data(55 * idx + 29);

            auto pc_yyyz = pcDistances.data(55 * idx + 30);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            auto pc_yzzz = pcDistances.data(55 * idx + 32);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxyzz = pcDistances.data(55 * idx + 42);

            auto pc_xyyyy = pcDistances.data(55 * idx + 44);

            auto pc_xyyyz = pcDistances.data(55 * idx + 45);

            auto pc_xyyzz = pcDistances.data(55 * idx + 46);

            auto pc_xyzzz = pcDistances.data(55 * idx + 47);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_y_xxzz = primBuffer.data(45 * idx + 20);

            auto t_y_xyyy = primBuffer.data(45 * idx + 21);

            auto t_y_xyyz = primBuffer.data(45 * idx + 22);

            auto t_y_xyzz = primBuffer.data(45 * idx + 23);

            auto t_y_xzzz = primBuffer.data(45 * idx + 24);

            // Batch of Integrals (20,25)

            #pragma omp simd aligned(fx, pa_y, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, pc_x, pc_xx, pc_xxy, pc_xxyz, pc_xxyzz, pc_xxz, pc_xxzz, pc_xy, pc_xyy, \
                                     pc_xyyy, pc_xyyyy, pc_xyyyz, pc_xyyz, pc_xyyzz, pc_xyz, pc_xyzz, pc_xyzzz, pc_xz, \
                                     pc_xzz, pc_xzzz, pc_y, pc_yy, pc_yyy, pc_yyyy, pc_yyyz, pc_yyz, pc_yyzz, pc_yz, pc_yzz, \
                                     pc_yzzz, pc_z, pc_zz, pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, \
                                     t_y_xxzz, t_y_xyyy, t_y_xyyz, t_y_xyzz, t_y_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_y_xxzz[j] = fl_s_0_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pa_y[j] * pb_xx[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pb_zz[j] + pa_y[j] * pb_xxzz[j]);

                t_y_xxzz[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl2_fx - 0.25 * pc_y[j] * fl2_fx - 0.5 * pa_y[j] * pb_xx[j] * fl1_fx - pa_y[j] * pb_x[j] * pc_x[j] * fl1_fx - pa_y[j] * fl1_fx * pb_z[j] * pc_z[j] - 0.5 * pa_y[j] * fl1_fx * pb_zz[j] - 0.5 * pc_y[j] * pb_xx[j] * fl1_fx - 0.5 * pc_y[j] * fl1_fx * pb_zz[j] - 2.0 * pa_y[j] * pb_xxz[j] * pc_z[j] - 2.0 * pa_y[j] * pb_xzz[j] * pc_x[j] - pc_y[j] * pb_xxzz[j]);

                t_y_xxzz[j] += fl_s_0_0_2 * (0.25 * pa_y[j] * fl2_fx + 0.5 * pc_y[j] * fl2_fx + pa_y[j] * pb_x[j] * pc_x[j] * fl1_fx + 0.5 * pa_y[j] * pc_xx[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pc_zz[j] + pa_y[j] * fl1_fx * pb_z[j] * pc_z[j] + 0.5 * pc_y[j] * pb_xx[j] * fl1_fx + pc_xy[j] * pb_x[j] * fl1_fx + pc_yz[j] * fl1_fx * pb_z[j] + 0.5 * pc_y[j] * fl1_fx * pb_zz[j] + pa_y[j] * pb_xx[j] * pc_zz[j] + 4.0 * pa_y[j] * pb_xz[j] * pc_xz[j] + pa_y[j] * pc_xx[j] * pb_zz[j] + 2.0 * pc_yz[j] * pb_xxz[j] + 2.0 * pc_xy[j] * pb_xzz[j]);

                t_y_xxzz[j] += fl_s_0_0_3 * (-0.25 * pc_y[j] * fl2_fx - 0.5 * pa_y[j] * pc_xx[j] * fl1_fx - 0.5 * pa_y[j] * fl1_fx * pc_zz[j] - pc_xy[j] * pb_x[j] * fl1_fx - 0.5 * pc_xxy[j] * fl1_fx - 0.5 * pc_yzz[j] * fl1_fx - pc_yz[j] * fl1_fx * pb_z[j] - 2.0 * pa_y[j] * pb_x[j] * pc_xzz[j] - 2.0 * pa_y[j] * pc_xxz[j] * pb_z[j] - pc_yzz[j] * pb_xx[j] - 4.0 * pc_xyz[j] * pb_xz[j] - pc_xxy[j] * pb_zz[j]);

                t_y_xxzz[j] += fl_s_0_0_4 * (0.5 * pc_xxy[j] * fl1_fx + 0.5 * pc_yzz[j] * fl1_fx + pa_y[j] * pc_xxzz[j] + 2.0 * pc_xyzz[j] * pb_x[j] + 2.0 * pc_xxyz[j] * pb_z[j]);

                t_y_xxzz[j] += -fl_s_0_0_5 * pc_xxyzz[j];

                t_y_xyyy[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_x[j] + 1.5 * pa_y[j] * pb_xy[j] * fl1_fx + 1.5 * fl1_fx * pb_xyy[j] + pa_y[j] * pb_xyyy[j]);

                t_y_xyyy[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_x[j] - 0.75 * fl2_fx * pc_x[j] - 1.5 * pa_y[j] * pb_xy[j] * fl1_fx - 1.5 * pa_y[j] * pb_x[j] * pc_y[j] * fl1_fx - 1.5 * pa_y[j] * pc_x[j] * pb_y[j] * fl1_fx - 4.5 * pc_y[j] * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * pc_x[j] * pb_yy[j] - 1.5 * fl1_fx * pb_xyy[j] - 3.0 * pa_y[j] * pb_xyy[j] * pc_y[j] - pa_y[j] * pc_x[j] * pb_yyy[j] - pc_y[j] * pb_xyyy[j]);

                t_y_xyyy[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_x[j] + 0.75 * fl2_fx * pb_x[j] + 1.5 * pa_y[j] * pb_x[j] * pc_y[j] * fl1_fx + 1.5 * pa_y[j] * pc_x[j] * pb_y[j] * fl1_fx + 1.5 * pa_y[j] * pc_xy[j] * fl1_fx + 4.5 * pc_y[j] * pb_xy[j] * fl1_fx + 3.0 * pc_yy[j] * pb_x[j] * fl1_fx + 4.5 * pc_xy[j] * pb_y[j] * fl1_fx + 1.5 * fl1_fx * pc_x[j] * pb_yy[j] + 3.0 * pa_y[j] * pb_xy[j] * pc_yy[j] + 3.0 * pa_y[j] * pc_xy[j] * pb_yy[j] + 3.0 * pc_yy[j] * pb_xyy[j] + pc_xy[j] * pb_yyy[j]);

                t_y_xyyy[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_x[j] - 1.5 * pa_y[j] * pc_xy[j] * fl1_fx - 3.0 * pc_yy[j] * pb_x[j] * fl1_fx - 4.5 * pc_xy[j] * pb_y[j] * fl1_fx - 3.0 * pc_xyy[j] * fl1_fx - pa_y[j] * pb_x[j] * pc_yyy[j] - 3.0 * pa_y[j] * pc_xyy[j] * pb_y[j] - 3.0 * pc_yyy[j] * pb_xy[j] - 3.0 * pc_xyy[j] * pb_yy[j]);

                t_y_xyyy[j] += fl_s_0_0_4 * (3.0 * pc_xyy[j] * fl1_fx + pa_y[j] * pc_xyyy[j] + pc_yyyy[j] * pb_x[j] + 3.0 * pc_xyyy[j] * pb_y[j]);

                t_y_xyyy[j] += -fl_s_0_0_5 * pc_xyyyy[j];

                t_y_xyyz[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * pb_xz[j] * fl1_fx + fl1_fx * pb_xyz[j] + pa_y[j] * pb_xyyz[j]);

                t_y_xyyz[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * pb_x[j] * fl1_fx * pc_z[j] - 0.5 * pa_y[j] * pb_xz[j] * fl1_fx - 0.5 * pa_y[j] * pc_x[j] * fl1_fx * pb_z[j] - 1.5 * pc_y[j] * pb_xz[j] * fl1_fx - fl1_fx * pb_xy[j] * pc_z[j] - fl1_fx * pc_x[j] * pb_yz[j] - fl1_fx * pb_xyz[j] - pa_y[j] * pb_xyy[j] * pc_z[j] - 2.0 * pa_y[j] * pb_xyz[j] * pc_y[j] - pa_y[j] * pc_x[j] * pb_yyz[j] - pc_y[j] * pb_xyyz[j]);

                t_y_xyyz[j] += fl_s_0_0_2 * (0.5 * pa_y[j] * pb_x[j] * fl1_fx * pc_z[j] + 0.5 * pa_y[j] * pc_xz[j] * fl1_fx + 0.5 * pa_y[j] * pc_x[j] * fl1_fx * pb_z[j] + 1.5 * pc_yz[j] * pb_x[j] * fl1_fx + 1.5 * pc_y[j] * pb_xz[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx * pb_z[j] + fl1_fx * pc_xz[j] * pb_y[j] + fl1_fx * pb_xy[j] * pc_z[j] + fl1_fx * pc_x[j] * pb_yz[j] + 2.0 * pa_y[j] * pb_xy[j] * pc_yz[j] + pa_y[j] * pb_xz[j] * pc_yy[j] + pa_y[j] * pc_xz[j] * pb_yy[j] + 2.0 * pa_y[j] * pc_xy[j] * pb_yz[j] + pc_yz[j] * pb_xyy[j] + 2.0 * pc_yy[j] * pb_xyz[j] + pc_xy[j] * pb_yyz[j]);

                t_y_xyyz[j] += fl_s_0_0_3 * (-0.5 * pa_y[j] * pc_xz[j] * fl1_fx - 1.5 * pc_yz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pb_z[j] - fl1_fx * pc_xz[j] * pb_y[j] - pa_y[j] * pb_x[j] * pc_yyz[j] - 2.0 * pa_y[j] * pc_xyz[j] * pb_y[j] - pa_y[j] * pc_xyy[j] * pb_z[j] - 2.0 * pc_yyz[j] * pb_xy[j] - pc_yyy[j] * pb_xz[j] - pc_xyz[j] * pb_yy[j] - 2.0 * pc_xyy[j] * pb_yz[j]);

                t_y_xyyz[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_y[j] * pc_xyyz[j] + pc_yyyz[j] * pb_x[j] + 2.0 * pc_xyyz[j] * pb_y[j] + pc_xyyy[j] * pb_z[j]);

                t_y_xyyz[j] += -fl_s_0_0_5 * pc_xyyyz[j];

                t_y_xyzz[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_x[j] + 0.5 * pa_y[j] * pb_xy[j] * fl1_fx + 0.5 * fl1_fx * pb_xzz[j] + pa_y[j] * pb_xyzz[j]);

                t_y_xyzz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx * pb_x[j] - 0.25 * fl2_fx * pc_x[j] - 0.5 * pa_y[j] * pb_xy[j] * fl1_fx - 0.5 * pa_y[j] * pb_x[j] * pc_y[j] * fl1_fx - 0.5 * pa_y[j] * pc_x[j] * pb_y[j] * fl1_fx - 0.5 * pc_y[j] * pb_xy[j] * fl1_fx - fl1_fx * pb_xz[j] * pc_z[j] - 0.5 * fl1_fx * pc_x[j] * pb_zz[j] - 0.5 * fl1_fx * pb_xzz[j] - 2.0 * pa_y[j] * pb_xyz[j] * pc_z[j] - pa_y[j] * pb_xzz[j] * pc_y[j] - pa_y[j] * pc_x[j] * pb_yzz[j] - pc_y[j] * pb_xyzz[j]);

                t_y_xyzz[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_x[j] + 0.25 * fl2_fx * pb_x[j] + 0.5 * pa_y[j] * pb_x[j] * pc_y[j] * fl1_fx + 0.5 * pa_y[j] * pc_x[j] * pb_y[j] * fl1_fx + 0.5 * pa_y[j] * pc_xy[j] * fl1_fx + 0.5 * pc_y[j] * pb_xy[j] * fl1_fx + 0.5 * pc_yy[j] * pb_x[j] * fl1_fx + 0.5 * pc_xy[j] * pb_y[j] * fl1_fx + 0.5 * fl1_fx * pb_x[j] * pc_zz[j] + fl1_fx * pc_xz[j] * pb_z[j] + fl1_fx * pb_xz[j] * pc_z[j] + 0.5 * fl1_fx * pc_x[j] * pb_zz[j] + pa_y[j] * pb_xy[j] * pc_zz[j] + 2.0 * pa_y[j] * pb_xz[j] * pc_yz[j] + 2.0 * pa_y[j] * pc_xz[j] * pb_yz[j] + pa_y[j] * pc_xy[j] * pb_zz[j] + 2.0 * pc_yz[j] * pb_xyz[j] + pc_yy[j] * pb_xzz[j] + pc_xy[j] * pb_yzz[j]);

                t_y_xyzz[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_x[j] - 0.5 * pa_y[j] * pc_xy[j] * fl1_fx - 0.5 * pc_yy[j] * pb_x[j] * fl1_fx - 0.5 * pc_xy[j] * pb_y[j] * fl1_fx - 0.5 * pc_xyy[j] * fl1_fx - 0.5 * fl1_fx * pc_xzz[j] - 0.5 * fl1_fx * pb_x[j] * pc_zz[j] - fl1_fx * pc_xz[j] * pb_z[j] - pa_y[j] * pb_x[j] * pc_yzz[j] - pa_y[j] * pc_xzz[j] * pb_y[j] - 2.0 * pa_y[j] * pc_xyz[j] * pb_z[j] - pc_yzz[j] * pb_xy[j] - 2.0 * pc_yyz[j] * pb_xz[j] - 2.0 * pc_xyz[j] * pb_yz[j] - pc_xyy[j] * pb_zz[j]);

                t_y_xyzz[j] += fl_s_0_0_4 * (0.5 * pc_xyy[j] * fl1_fx + 0.5 * fl1_fx * pc_xzz[j] + pa_y[j] * pc_xyzz[j] + pc_yyzz[j] * pb_x[j] + pc_xyzz[j] * pb_y[j] + 2.0 * pc_xyyz[j] * pb_z[j]);

                t_y_xyzz[j] += -fl_s_0_0_5 * pc_xyyzz[j];

                t_y_xzzz[j] = fl_s_0_0_0 * (1.5 * pa_y[j] * pb_xz[j] * fl1_fx + pa_y[j] * pb_xzzz[j]);

                t_y_xzzz[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * pb_xz[j] * fl1_fx - 1.5 * pa_y[j] * pb_x[j] * pc_z[j] * fl1_fx - 1.5 * pa_y[j] * pc_x[j] * pb_z[j] * fl1_fx - 1.5 * pc_y[j] * pb_xz[j] * fl1_fx - 3.0 * pa_y[j] * pb_xzz[j] * pc_z[j] - pa_y[j] * pc_x[j] * pb_zzz[j] - pc_y[j] * pb_xzzz[j]);

                t_y_xzzz[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * pb_x[j] * pc_z[j] * fl1_fx + 1.5 * pa_y[j] * pc_x[j] * pb_z[j] * fl1_fx + 1.5 * pa_y[j] * pc_xz[j] * fl1_fx + 1.5 * pc_y[j] * pb_xz[j] * fl1_fx + 1.5 * pc_yz[j] * pb_x[j] * fl1_fx + 1.5 * pc_xy[j] * pb_z[j] * fl1_fx + 3.0 * pa_y[j] * pb_xz[j] * pc_zz[j] + 3.0 * pa_y[j] * pc_xz[j] * pb_zz[j] + 3.0 * pc_yz[j] * pb_xzz[j] + pc_xy[j] * pb_zzz[j]);

                t_y_xzzz[j] += fl_s_0_0_3 * (-1.5 * pa_y[j] * pc_xz[j] * fl1_fx - 1.5 * pc_yz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xy[j] * pb_z[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - pa_y[j] * pb_x[j] * pc_zzz[j] - 3.0 * pa_y[j] * pc_xzz[j] * pb_z[j] - 3.0 * pc_yzz[j] * pb_xz[j] - 3.0 * pc_xyz[j] * pb_zz[j]);

                t_y_xzzz[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_y[j] * pc_xzzz[j] + pc_yzzz[j] * pb_x[j] + 3.0 * pc_xyzz[j] * pb_z[j]);

                t_y_xzzz[j] += -fl_s_0_0_5 * pc_xyzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG_25_30(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (25,30)

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

            auto pa_y = paDistances.data(3 * idx + 1);

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

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_yyyy = pcDistances.data(55 * idx + 29);

            auto pc_yyyz = pcDistances.data(55 * idx + 30);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            auto pc_yzzz = pcDistances.data(55 * idx + 32);

            auto pc_zzzz = pcDistances.data(55 * idx + 33);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_yyyyy = pcDistances.data(55 * idx + 49);

            auto pc_yyyyz = pcDistances.data(55 * idx + 50);

            auto pc_yyyzz = pcDistances.data(55 * idx + 51);

            auto pc_yyzzz = pcDistances.data(55 * idx + 52);

            auto pc_yzzzz = pcDistances.data(55 * idx + 53);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_y_yyyy = primBuffer.data(45 * idx + 25);

            auto t_y_yyyz = primBuffer.data(45 * idx + 26);

            auto t_y_yyzz = primBuffer.data(45 * idx + 27);

            auto t_y_yzzz = primBuffer.data(45 * idx + 28);

            auto t_y_zzzz = primBuffer.data(45 * idx + 29);

            // Batch of Integrals (25,30)

            #pragma omp simd aligned(fx, pa_y, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, pc_y, pc_yy, pc_yyy, pc_yyyy, pc_yyyyy, pc_yyyyz, \
                                     pc_yyyz, pc_yyyzz, pc_yyz, pc_yyzz, pc_yyzzz, pc_yz, pc_yzz, pc_yzzz, pc_yzzzz, pc_z, \
                                     pc_zz, pc_zzz, pc_zzzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, \
                                     t_y_yyyy, t_y_yyyz, t_y_yyzz, t_y_yzzz, t_y_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_y_yyyy[j] = fl_s_0_0_0 * (0.75 * pa_y[j] * fl2_fx + 3.0 * fl2_fx * pb_y[j] + 3.0 * pa_y[j] * pb_yy[j] * fl1_fx + 2.0 * fl1_fx * pb_yyy[j] + pa_y[j] * pb_yyyy[j]);

                t_y_yyyy[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl2_fx - 3.75 * pc_y[j] * fl2_fx - 6.0 * fl2_fx * pb_y[j] - 3.0 * pa_y[j] * pb_yy[j] * fl1_fx - 6.0 * pa_y[j] * pb_y[j] * pc_y[j] * fl1_fx - 9.0 * pc_y[j] * pb_yy[j] * fl1_fx - 2.0 * fl1_fx * pb_yyy[j] - 4.0 * pa_y[j] * pb_yyy[j] * pc_y[j] - pc_y[j] * pb_yyyy[j]);

                t_y_yyyy[j] += fl_s_0_0_2 * (0.75 * pa_y[j] * fl2_fx + 7.5 * pc_y[j] * fl2_fx + 3.0 * fl2_fx * pb_y[j] + 6.0 * pa_y[j] * pb_y[j] * pc_y[j] * fl1_fx + 3.0 * pa_y[j] * pc_yy[j] * fl1_fx + 9.0 * pc_y[j] * pb_yy[j] * fl1_fx + 12.0 * pc_yy[j] * pb_y[j] * fl1_fx + 6.0 * pa_y[j] * pb_yy[j] * pc_yy[j] + 4.0 * pc_yy[j] * pb_yyy[j]);

                t_y_yyyy[j] += fl_s_0_0_3 * (-3.75 * pc_y[j] * fl2_fx - 3.0 * pa_y[j] * pc_yy[j] * fl1_fx - 12.0 * pc_yy[j] * pb_y[j] * fl1_fx - 5.0 * pc_yyy[j] * fl1_fx - 4.0 * pa_y[j] * pb_y[j] * pc_yyy[j] - 6.0 * pc_yyy[j] * pb_yy[j]);

                t_y_yyyy[j] += fl_s_0_0_4 * (5.0 * pc_yyy[j] * fl1_fx + pa_y[j] * pc_yyyy[j] + 4.0 * pc_yyyy[j] * pb_y[j]);

                t_y_yyyy[j] += -fl_s_0_0_5 * pc_yyyyy[j];

                t_y_yyyz[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_z[j] + 1.5 * pa_y[j] * pb_yz[j] * fl1_fx + 1.5 * fl1_fx * pb_yyz[j] + pa_y[j] * pb_yyyz[j]);

                t_y_yyyz[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * fl2_fx * pb_z[j] - 1.5 * pa_y[j] * pb_y[j] * fl1_fx * pc_z[j] - 1.5 * pa_y[j] * pb_yz[j] * fl1_fx - 1.5 * pa_y[j] * pc_y[j] * fl1_fx * pb_z[j] - 4.5 * pc_y[j] * pb_yz[j] * fl1_fx - 1.5 * fl1_fx * pb_yy[j] * pc_z[j] - 1.5 * fl1_fx * pb_yyz[j] - pa_y[j] * pb_yyy[j] * pc_z[j] - 3.0 * pa_y[j] * pb_yyz[j] * pc_y[j] - pc_y[j] * pb_yyyz[j]);

                t_y_yyyz[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pb_z[j] + 1.5 * pa_y[j] * pb_y[j] * fl1_fx * pc_z[j] + 1.5 * pa_y[j] * pc_yz[j] * fl1_fx + 1.5 * pa_y[j] * pc_y[j] * fl1_fx * pb_z[j] + 4.5 * pc_yz[j] * pb_y[j] * fl1_fx + 4.5 * pc_y[j] * pb_yz[j] * fl1_fx + 3.0 * pc_yy[j] * fl1_fx * pb_z[j] + 1.5 * fl1_fx * pb_yy[j] * pc_z[j] + 3.0 * pa_y[j] * pb_yy[j] * pc_yz[j] + 3.0 * pa_y[j] * pb_yz[j] * pc_yy[j] + pc_yz[j] * pb_yyy[j] + 3.0 * pc_yy[j] * pb_yyz[j]);

                t_y_yyyz[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * pa_y[j] * pc_yz[j] * fl1_fx - 4.5 * pc_yz[j] * pb_y[j] * fl1_fx - 3.0 * pc_yyz[j] * fl1_fx - 3.0 * pc_yy[j] * fl1_fx * pb_z[j] - 3.0 * pa_y[j] * pb_y[j] * pc_yyz[j] - pa_y[j] * pc_yyy[j] * pb_z[j] - 3.0 * pc_yyz[j] * pb_yy[j] - 3.0 * pc_yyy[j] * pb_yz[j]);

                t_y_yyyz[j] += fl_s_0_0_4 * (3.0 * pc_yyz[j] * fl1_fx + pa_y[j] * pc_yyyz[j] + 3.0 * pc_yyyz[j] * pb_y[j] + pc_yyyy[j] * pb_z[j]);

                t_y_yyyz[j] += -fl_s_0_0_5 * pc_yyyyz[j];

                t_y_yyzz[j] = fl_s_0_0_0 * (0.25 * pa_y[j] * fl2_fx + 0.5 * fl2_fx * pb_y[j] + 0.5 * pa_y[j] * pb_yy[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pb_zz[j] + fl1_fx * pb_yzz[j] + pa_y[j] * pb_yyzz[j]);

                t_y_yyzz[j] += fl_s_0_0_1 * (-0.5 * pa_y[j] * fl2_fx - 0.75 * pc_y[j] * fl2_fx - fl2_fx * pb_y[j] - 0.5 * pa_y[j] * pb_yy[j] * fl1_fx - pa_y[j] * pb_y[j] * pc_y[j] * fl1_fx - pa_y[j] * fl1_fx * pb_z[j] * pc_z[j] - 0.5 * pa_y[j] * fl1_fx * pb_zz[j] - 0.5 * pc_y[j] * pb_yy[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx * pb_zz[j] - 2.0 * fl1_fx * pb_yz[j] * pc_z[j] - fl1_fx * pb_yzz[j] - 2.0 * pa_y[j] * pb_yyz[j] * pc_z[j] - 2.0 * pa_y[j] * pb_yzz[j] * pc_y[j] - pc_y[j] * pb_yyzz[j]);

                t_y_yyzz[j] += fl_s_0_0_2 * (0.25 * pa_y[j] * fl2_fx + 1.5 * pc_y[j] * fl2_fx + 0.5 * fl2_fx * pb_y[j] + pa_y[j] * pb_y[j] * pc_y[j] * fl1_fx + 0.5 * pa_y[j] * pc_yy[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pc_zz[j] + pa_y[j] * fl1_fx * pb_z[j] * pc_z[j] + 0.5 * pc_y[j] * pb_yy[j] * fl1_fx + pc_yy[j] * pb_y[j] * fl1_fx + 3.0 * pc_yz[j] * fl1_fx * pb_z[j] + 1.5 * pc_y[j] * fl1_fx * pb_zz[j] + fl1_fx * pb_y[j] * pc_zz[j] + 2.0 * fl1_fx * pb_yz[j] * pc_z[j] + pa_y[j] * pb_yy[j] * pc_zz[j] + 4.0 * pa_y[j] * pb_yz[j] * pc_yz[j] + pa_y[j] * pc_yy[j] * pb_zz[j] + 2.0 * pc_yz[j] * pb_yyz[j] + 2.0 * pc_yy[j] * pb_yzz[j]);

                t_y_yyzz[j] += fl_s_0_0_3 * (-0.75 * pc_y[j] * fl2_fx - 0.5 * pa_y[j] * pc_yy[j] * fl1_fx - 0.5 * pa_y[j] * fl1_fx * pc_zz[j] - pc_yy[j] * pb_y[j] * fl1_fx - 0.5 * pc_yyy[j] * fl1_fx - 1.5 * pc_yzz[j] * fl1_fx - 3.0 * pc_yz[j] * fl1_fx * pb_z[j] - fl1_fx * pb_y[j] * pc_zz[j] - 2.0 * pa_y[j] * pb_y[j] * pc_yzz[j] - 2.0 * pa_y[j] * pc_yyz[j] * pb_z[j] - pc_yzz[j] * pb_yy[j] - 4.0 * pc_yyz[j] * pb_yz[j] - pc_yyy[j] * pb_zz[j]);

                t_y_yyzz[j] += fl_s_0_0_4 * (0.5 * pc_yyy[j] * fl1_fx + 1.5 * pc_yzz[j] * fl1_fx + pa_y[j] * pc_yyzz[j] + 2.0 * pc_yyzz[j] * pb_y[j] + 2.0 * pc_yyyz[j] * pb_z[j]);

                t_y_yyzz[j] += -fl_s_0_0_5 * pc_yyyzz[j];

                t_y_yzzz[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_z[j] + 1.5 * pa_y[j] * pb_yz[j] * fl1_fx + 0.5 * fl1_fx * pb_zzz[j] + pa_y[j] * pb_yzzz[j]);

                t_y_yzzz[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_z[j] - 0.75 * fl2_fx * pc_z[j] - 1.5 * pa_y[j] * pb_yz[j] * fl1_fx - 1.5 * pa_y[j] * pb_y[j] * pc_z[j] * fl1_fx - 1.5 * pa_y[j] * pc_y[j] * pb_z[j] * fl1_fx - 1.5 * pc_y[j] * pb_yz[j] * fl1_fx - 1.5 * fl1_fx * pb_zz[j] * pc_z[j] - 0.5 * fl1_fx * pb_zzz[j] - 3.0 * pa_y[j] * pb_yzz[j] * pc_z[j] - pa_y[j] * pc_y[j] * pb_zzz[j] - pc_y[j] * pb_yzzz[j]);

                t_y_yzzz[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pb_z[j] + 1.5 * pa_y[j] * pb_y[j] * pc_z[j] * fl1_fx + 1.5 * pa_y[j] * pc_y[j] * pb_z[j] * fl1_fx + 1.5 * pa_y[j] * pc_yz[j] * fl1_fx + 1.5 * pc_y[j] * pb_yz[j] * fl1_fx + 1.5 * pc_yz[j] * pb_y[j] * fl1_fx + 1.5 * pc_yy[j] * pb_z[j] * fl1_fx + 1.5 * fl1_fx * pb_z[j] * pc_zz[j] + 1.5 * fl1_fx * pb_zz[j] * pc_z[j] + 3.0 * pa_y[j] * pb_yz[j] * pc_zz[j] + 3.0 * pa_y[j] * pc_yz[j] * pb_zz[j] + 3.0 * pc_yz[j] * pb_yzz[j] + pc_yy[j] * pb_zzz[j]);

                t_y_yzzz[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * pa_y[j] * pc_yz[j] * fl1_fx - 1.5 * pc_yz[j] * pb_y[j] * fl1_fx - 1.5 * pc_yy[j] * pb_z[j] * fl1_fx - 1.5 * pc_yyz[j] * fl1_fx - 0.5 * fl1_fx * pc_zzz[j] - 1.5 * fl1_fx * pb_z[j] * pc_zz[j] - pa_y[j] * pb_y[j] * pc_zzz[j] - 3.0 * pa_y[j] * pc_yzz[j] * pb_z[j] - 3.0 * pc_yzz[j] * pb_yz[j] - 3.0 * pc_yyz[j] * pb_zz[j]);

                t_y_yzzz[j] += fl_s_0_0_4 * (1.5 * pc_yyz[j] * fl1_fx + 0.5 * fl1_fx * pc_zzz[j] + pa_y[j] * pc_yzzz[j] + pc_yzzz[j] * pb_y[j] + 3.0 * pc_yyzz[j] * pb_z[j]);

                t_y_yzzz[j] += -fl_s_0_0_5 * pc_yyzzz[j];

                t_y_zzzz[j] = fl_s_0_0_0 * (0.75 * pa_y[j] * fl2_fx + 3.0 * pa_y[j] * pb_zz[j] * fl1_fx + pa_y[j] * pb_zzzz[j]);

                t_y_zzzz[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl2_fx - 0.75 * pc_y[j] * fl2_fx - 3.0 * pa_y[j] * pb_zz[j] * fl1_fx - 6.0 * pa_y[j] * pb_z[j] * pc_z[j] * fl1_fx - 3.0 * pc_y[j] * pb_zz[j] * fl1_fx - 4.0 * pa_y[j] * pb_zzz[j] * pc_z[j] - pc_y[j] * pb_zzzz[j]);

                t_y_zzzz[j] += fl_s_0_0_2 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pc_y[j] * fl2_fx + 6.0 * pa_y[j] * pb_z[j] * pc_z[j] * fl1_fx + 3.0 * pa_y[j] * pc_zz[j] * fl1_fx + 3.0 * pc_y[j] * pb_zz[j] * fl1_fx + 6.0 * pc_yz[j] * pb_z[j] * fl1_fx + 6.0 * pa_y[j] * pb_zz[j] * pc_zz[j] + 4.0 * pc_yz[j] * pb_zzz[j]);

                t_y_zzzz[j] += fl_s_0_0_3 * (-0.75 * pc_y[j] * fl2_fx - 3.0 * pa_y[j] * pc_zz[j] * fl1_fx - 6.0 * pc_yz[j] * pb_z[j] * fl1_fx - 3.0 * pc_yzz[j] * fl1_fx - 4.0 * pa_y[j] * pb_z[j] * pc_zzz[j] - 6.0 * pc_yzz[j] * pb_zz[j]);

                t_y_zzzz[j] += fl_s_0_0_4 * (3.0 * pc_yzz[j] * fl1_fx + pa_y[j] * pc_zzzz[j] + 4.0 * pc_yzzz[j] * pb_z[j]);

                t_y_zzzz[j] += -fl_s_0_0_5 * pc_yzzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG_30_35(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (30,35)

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

            auto pa_z = paDistances.data(3 * idx + 2);

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

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(55 * idx + 9);

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(55 * idx + 19);

            auto pc_xxxy = pcDistances.data(55 * idx + 20);

            auto pc_xxxz = pcDistances.data(55 * idx + 21);

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxxxz = pcDistances.data(55 * idx + 36);

            auto pc_xxxyz = pcDistances.data(55 * idx + 38);

            auto pc_xxxzz = pcDistances.data(55 * idx + 39);

            auto pc_xxyyz = pcDistances.data(55 * idx + 41);

            auto pc_xxyzz = pcDistances.data(55 * idx + 42);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_z_xxxx = primBuffer.data(45 * idx + 30);

            auto t_z_xxxy = primBuffer.data(45 * idx + 31);

            auto t_z_xxxz = primBuffer.data(45 * idx + 32);

            auto t_z_xxyy = primBuffer.data(45 * idx + 33);

            auto t_z_xxyz = primBuffer.data(45 * idx + 34);

            // Batch of Integrals (30,35)

            #pragma omp simd aligned(fx, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pc_x, pc_xx, \
                                     pc_xxx, pc_xxxx, pc_xxxxz, pc_xxxy, pc_xxxyz, pc_xxxz, pc_xxxzz, pc_xxy, pc_xxyy, \
                                     pc_xxyyz, pc_xxyz, pc_xxyzz, pc_xxz, pc_xxzz, pc_xy, pc_xyy, pc_xyyz, pc_xyz, pc_xyzz, \
                                     pc_xz, pc_xzz, pc_y, pc_yy, pc_yyz, pc_yz, pc_yzz, pc_z, pc_zz, s_0_0_0, s_0_0_1, \
                                     s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_z_xxxx, t_z_xxxy, t_z_xxxz, t_z_xxyy, \
                                     t_z_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_z_xxxx[j] = fl_s_0_0_0 * (0.75 * pa_z[j] * fl2_fx + 3.0 * pa_z[j] * pb_xx[j] * fl1_fx + pa_z[j] * pb_xxxx[j]);

                t_z_xxxx[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * fl2_fx - 0.75 * pc_z[j] * fl2_fx - 3.0 * pa_z[j] * pb_xx[j] * fl1_fx - 6.0 * pa_z[j] * pb_x[j] * pc_x[j] * fl1_fx - 3.0 * pc_z[j] * pb_xx[j] * fl1_fx - 4.0 * pa_z[j] * pb_xxx[j] * pc_x[j] - pc_z[j] * pb_xxxx[j]);

                t_z_xxxx[j] += fl_s_0_0_2 * (0.75 * pa_z[j] * fl2_fx + 1.5 * pc_z[j] * fl2_fx + 6.0 * pa_z[j] * pb_x[j] * pc_x[j] * fl1_fx + 3.0 * pa_z[j] * pc_xx[j] * fl1_fx + 3.0 * pc_z[j] * pb_xx[j] * fl1_fx + 6.0 * pc_xz[j] * pb_x[j] * fl1_fx + 6.0 * pa_z[j] * pb_xx[j] * pc_xx[j] + 4.0 * pc_xz[j] * pb_xxx[j]);

                t_z_xxxx[j] += fl_s_0_0_3 * (-0.75 * pc_z[j] * fl2_fx - 3.0 * pa_z[j] * pc_xx[j] * fl1_fx - 6.0 * pc_xz[j] * pb_x[j] * fl1_fx - 3.0 * pc_xxz[j] * fl1_fx - 4.0 * pa_z[j] * pb_x[j] * pc_xxx[j] - 6.0 * pc_xxz[j] * pb_xx[j]);

                t_z_xxxx[j] += fl_s_0_0_4 * (3.0 * pc_xxz[j] * fl1_fx + pa_z[j] * pc_xxxx[j] + 4.0 * pc_xxxz[j] * pb_x[j]);

                t_z_xxxx[j] += -fl_s_0_0_5 * pc_xxxxz[j];

                t_z_xxxy[j] = fl_s_0_0_0 * (1.5 * pa_z[j] * pb_xy[j] * fl1_fx + pa_z[j] * pb_xxxy[j]);

                t_z_xxxy[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * pb_x[j] * fl1_fx * pc_y[j] - 1.5 * pa_z[j] * pb_xy[j] * fl1_fx - 1.5 * pa_z[j] * pc_x[j] * fl1_fx * pb_y[j] - 1.5 * pc_z[j] * pb_xy[j] * fl1_fx - pa_z[j] * pb_xxx[j] * pc_y[j] - 3.0 * pa_z[j] * pb_xxy[j] * pc_x[j] - pc_z[j] * pb_xxxy[j]);

                t_z_xxxy[j] += fl_s_0_0_2 * (1.5 * pa_z[j] * pb_x[j] * fl1_fx * pc_y[j] + 1.5 * pa_z[j] * pc_xy[j] * fl1_fx + 1.5 * pa_z[j] * pc_x[j] * fl1_fx * pb_y[j] + 1.5 * pc_yz[j] * pb_x[j] * fl1_fx + 1.5 * pc_z[j] * pb_xy[j] * fl1_fx + 1.5 * pc_xz[j] * fl1_fx * pb_y[j] + 3.0 * pa_z[j] * pb_xx[j] * pc_xy[j] + 3.0 * pa_z[j] * pb_xy[j] * pc_xx[j] + pc_yz[j] * pb_xxx[j] + 3.0 * pc_xz[j] * pb_xxy[j]);

                t_z_xxxy[j] += fl_s_0_0_3 * (-1.5 * pa_z[j] * pc_xy[j] * fl1_fx - 1.5 * pc_yz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pb_y[j] - 3.0 * pa_z[j] * pb_x[j] * pc_xxy[j] - pa_z[j] * pc_xxx[j] * pb_y[j] - 3.0 * pc_xyz[j] * pb_xx[j] - 3.0 * pc_xxz[j] * pb_xy[j]);

                t_z_xxxy[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_z[j] * pc_xxxy[j] + 3.0 * pc_xxyz[j] * pb_x[j] + pc_xxxz[j] * pb_y[j]);

                t_z_xxxy[j] += -fl_s_0_0_5 * pc_xxxyz[j];

                t_z_xxxz[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_x[j] + 1.5 * pa_z[j] * pb_xz[j] * fl1_fx + 0.5 * fl1_fx * pb_xxx[j] + pa_z[j] * pb_xxxz[j]);

                t_z_xxxz[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_x[j] - 0.75 * fl2_fx * pc_x[j] - 1.5 * pa_z[j] * pb_x[j] * fl1_fx * pc_z[j] - 1.5 * pa_z[j] * pb_xz[j] * fl1_fx - 1.5 * pa_z[j] * pc_x[j] * fl1_fx * pb_z[j] - 1.5 * pc_z[j] * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * pb_xx[j] * pc_x[j] - 0.5 * fl1_fx * pb_xxx[j] - pa_z[j] * pb_xxx[j] * pc_z[j] - 3.0 * pa_z[j] * pb_xxz[j] * pc_x[j] - pc_z[j] * pb_xxxz[j]);

                t_z_xxxz[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_x[j] + 0.75 * fl2_fx * pb_x[j] + 1.5 * pa_z[j] * pb_x[j] * fl1_fx * pc_z[j] + 1.5 * pa_z[j] * pc_xz[j] * fl1_fx + 1.5 * pa_z[j] * pc_x[j] * fl1_fx * pb_z[j] + 1.5 * pc_zz[j] * pb_x[j] * fl1_fx + 1.5 * pc_z[j] * pb_xz[j] * fl1_fx + 1.5 * pc_xz[j] * fl1_fx * pb_z[j] + 1.5 * fl1_fx * pb_x[j] * pc_xx[j] + 1.5 * fl1_fx * pb_xx[j] * pc_x[j] + 3.0 * pa_z[j] * pb_xx[j] * pc_xz[j] + 3.0 * pa_z[j] * pb_xz[j] * pc_xx[j] + pc_zz[j] * pb_xxx[j] + 3.0 * pc_xz[j] * pb_xxz[j]);

                t_z_xxxz[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_x[j] - 1.5 * pa_z[j] * pc_xz[j] * fl1_fx - 1.5 * pc_zz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xzz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_xxx[j] - 1.5 * fl1_fx * pb_x[j] * pc_xx[j] - 3.0 * pa_z[j] * pb_x[j] * pc_xxz[j] - pa_z[j] * pc_xxx[j] * pb_z[j] - 3.0 * pc_xzz[j] * pb_xx[j] - 3.0 * pc_xxz[j] * pb_xz[j]);

                t_z_xxxz[j] += fl_s_0_0_4 * (1.5 * pc_xzz[j] * fl1_fx + 0.5 * fl1_fx * pc_xxx[j] + pa_z[j] * pc_xxxz[j] + 3.0 * pc_xxzz[j] * pb_x[j] + pc_xxxz[j] * pb_z[j]);

                t_z_xxxz[j] += -fl_s_0_0_5 * pc_xxxzz[j];

                t_z_xxyy[j] = fl_s_0_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pa_z[j] * pb_xx[j] * fl1_fx + 0.5 * pa_z[j] * fl1_fx * pb_yy[j] + pa_z[j] * pb_xxyy[j]);

                t_z_xxyy[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * fl2_fx - 0.25 * pc_z[j] * fl2_fx - 0.5 * pa_z[j] * pb_xx[j] * fl1_fx - pa_z[j] * pb_x[j] * pc_x[j] * fl1_fx - pa_z[j] * fl1_fx * pb_y[j] * pc_y[j] - 0.5 * pa_z[j] * fl1_fx * pb_yy[j] - 0.5 * pc_z[j] * pb_xx[j] * fl1_fx - 0.5 * pc_z[j] * fl1_fx * pb_yy[j] - 2.0 * pa_z[j] * pb_xxy[j] * pc_y[j] - 2.0 * pa_z[j] * pb_xyy[j] * pc_x[j] - pc_z[j] * pb_xxyy[j]);

                t_z_xxyy[j] += fl_s_0_0_2 * (0.25 * pa_z[j] * fl2_fx + 0.5 * pc_z[j] * fl2_fx + pa_z[j] * pb_x[j] * pc_x[j] * fl1_fx + 0.5 * pa_z[j] * pc_xx[j] * fl1_fx + 0.5 * pa_z[j] * fl1_fx * pc_yy[j] + pa_z[j] * fl1_fx * pb_y[j] * pc_y[j] + 0.5 * pc_z[j] * pb_xx[j] * fl1_fx + pc_xz[j] * pb_x[j] * fl1_fx + pc_yz[j] * fl1_fx * pb_y[j] + 0.5 * pc_z[j] * fl1_fx * pb_yy[j] + pa_z[j] * pb_xx[j] * pc_yy[j] + 4.0 * pa_z[j] * pb_xy[j] * pc_xy[j] + pa_z[j] * pc_xx[j] * pb_yy[j] + 2.0 * pc_yz[j] * pb_xxy[j] + 2.0 * pc_xz[j] * pb_xyy[j]);

                t_z_xxyy[j] += fl_s_0_0_3 * (-0.25 * pc_z[j] * fl2_fx - 0.5 * pa_z[j] * pc_xx[j] * fl1_fx - 0.5 * pa_z[j] * fl1_fx * pc_yy[j] - pc_xz[j] * pb_x[j] * fl1_fx - 0.5 * pc_xxz[j] * fl1_fx - 0.5 * pc_yyz[j] * fl1_fx - pc_yz[j] * fl1_fx * pb_y[j] - 2.0 * pa_z[j] * pb_x[j] * pc_xyy[j] - 2.0 * pa_z[j] * pc_xxy[j] * pb_y[j] - pc_yyz[j] * pb_xx[j] - 4.0 * pc_xyz[j] * pb_xy[j] - pc_xxz[j] * pb_yy[j]);

                t_z_xxyy[j] += fl_s_0_0_4 * (0.5 * pc_xxz[j] * fl1_fx + 0.5 * pc_yyz[j] * fl1_fx + pa_z[j] * pc_xxyy[j] + 2.0 * pc_xyyz[j] * pb_x[j] + 2.0 * pc_xxyz[j] * pb_y[j]);

                t_z_xxyy[j] += -fl_s_0_0_5 * pc_xxyyz[j];

                t_z_xxyz[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_y[j] + 0.5 * pa_z[j] * fl1_fx * pb_yz[j] + 0.5 * fl1_fx * pb_xxy[j] + pa_z[j] * pb_xxyz[j]);

                t_z_xxyz[j] += fl_s_0_0_1 * (-0.25 * fl2_fx * pc_y[j] - 0.5 * fl2_fx * pb_y[j] - 0.5 * pa_z[j] * fl1_fx * pb_y[j] * pc_z[j] - 0.5 * pa_z[j] * fl1_fx * pc_y[j] * pb_z[j] - 0.5 * pa_z[j] * fl1_fx * pb_yz[j] - 0.5 * pc_z[j] * fl1_fx * pb_yz[j] - 0.5 * fl1_fx * pb_xx[j] * pc_y[j] - fl1_fx * pb_xy[j] * pc_x[j] - 0.5 * fl1_fx * pb_xxy[j] - pa_z[j] * pb_xxy[j] * pc_z[j] - pa_z[j] * pb_xxz[j] * pc_y[j] - 2.0 * pa_z[j] * pb_xyz[j] * pc_x[j] - pc_z[j] * pb_xxyz[j]);

                t_z_xxyz[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_y[j] + 0.25 * fl2_fx * pb_y[j] + 0.5 * pa_z[j] * fl1_fx * pc_yz[j] + 0.5 * pa_z[j] * fl1_fx * pb_y[j] * pc_z[j] + 0.5 * pa_z[j] * fl1_fx * pc_y[j] * pb_z[j] + 0.5 * pc_zz[j] * fl1_fx * pb_y[j] + 0.5 * pc_yz[j] * fl1_fx * pb_z[j] + 0.5 * pc_z[j] * fl1_fx * pb_yz[j] + fl1_fx * pb_x[j] * pc_xy[j] + 0.5 * fl1_fx * pc_xx[j] * pb_y[j] + 0.5 * fl1_fx * pb_xx[j] * pc_y[j] + fl1_fx * pb_xy[j] * pc_x[j] + pa_z[j] * pb_xx[j] * pc_yz[j] + 2.0 * pa_z[j] * pb_xy[j] * pc_xz[j] + 2.0 * pa_z[j] * pb_xz[j] * pc_xy[j] + pa_z[j] * pc_xx[j] * pb_yz[j] + pc_zz[j] * pb_xxy[j] + pc_yz[j] * pb_xxz[j] + 2.0 * pc_xz[j] * pb_xyz[j]);

                t_z_xxyz[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_y[j] - 0.5 * pa_z[j] * fl1_fx * pc_yz[j] - 0.5 * pc_yzz[j] * fl1_fx - 0.5 * pc_zz[j] * fl1_fx * pb_y[j] - 0.5 * pc_yz[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_xxy[j] - fl1_fx * pb_x[j] * pc_xy[j] - 0.5 * fl1_fx * pc_xx[j] * pb_y[j] - 2.0 * pa_z[j] * pb_x[j] * pc_xyz[j] - pa_z[j] * pc_xxz[j] * pb_y[j] - pa_z[j] * pc_xxy[j] * pb_z[j] - pc_yzz[j] * pb_xx[j] - 2.0 * pc_xzz[j] * pb_xy[j] - 2.0 * pc_xyz[j] * pb_xz[j] - pc_xxz[j] * pb_yz[j]);

                t_z_xxyz[j] += fl_s_0_0_4 * (0.5 * pc_yzz[j] * fl1_fx + 0.5 * fl1_fx * pc_xxy[j] + pa_z[j] * pc_xxyz[j] + 2.0 * pc_xyzz[j] * pb_x[j] + pc_xxzz[j] * pb_y[j] + pc_xxyz[j] * pb_z[j]);

                t_z_xxyz[j] += -fl_s_0_0_5 * pc_xxyzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG_35_40(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (35,40)

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

            auto pa_z = paDistances.data(3 * idx + 2);

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

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            auto pc_xzzz = pcDistances.data(55 * idx + 28);

            auto pc_yyyz = pcDistances.data(55 * idx + 30);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            auto pc_yzzz = pcDistances.data(55 * idx + 32);

            auto pc_zzzz = pcDistances.data(55 * idx + 33);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxzzz = pcDistances.data(55 * idx + 43);

            auto pc_xyyyz = pcDistances.data(55 * idx + 45);

            auto pc_xyyzz = pcDistances.data(55 * idx + 46);

            auto pc_xyzzz = pcDistances.data(55 * idx + 47);

            auto pc_xzzzz = pcDistances.data(55 * idx + 48);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_z_xxzz = primBuffer.data(45 * idx + 35);

            auto t_z_xyyy = primBuffer.data(45 * idx + 36);

            auto t_z_xyyz = primBuffer.data(45 * idx + 37);

            auto t_z_xyzz = primBuffer.data(45 * idx + 38);

            auto t_z_xzzz = primBuffer.data(45 * idx + 39);

            // Batch of Integrals (35,40)

            #pragma omp simd aligned(fx, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, pc_x, pc_xx, pc_xxz, pc_xxzz, pc_xxzzz, pc_xy, pc_xyy, pc_xyyy, pc_xyyyz, \
                                     pc_xyyz, pc_xyyzz, pc_xyz, pc_xyzz, pc_xyzzz, pc_xz, pc_xzz, pc_xzzz, pc_xzzzz, pc_y, \
                                     pc_yy, pc_yyy, pc_yyyz, pc_yyz, pc_yyzz, pc_yz, pc_yzz, pc_yzzz, pc_z, pc_zz, pc_zzz, \
                                     pc_zzzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_z_xxzz, t_z_xyyy, \
                                     t_z_xyyz, t_z_xyzz, t_z_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_z_xxzz[j] = fl_s_0_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * fl2_fx * pb_z[j] + 0.5 * pa_z[j] * pb_xx[j] * fl1_fx + 0.5 * pa_z[j] * fl1_fx * pb_zz[j] + fl1_fx * pb_xxz[j] + pa_z[j] * pb_xxzz[j]);

                t_z_xxzz[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * fl2_fx - 0.75 * pc_z[j] * fl2_fx - fl2_fx * pb_z[j] - 0.5 * pa_z[j] * pb_xx[j] * fl1_fx - pa_z[j] * pb_x[j] * pc_x[j] * fl1_fx - pa_z[j] * fl1_fx * pb_z[j] * pc_z[j] - 0.5 * pa_z[j] * fl1_fx * pb_zz[j] - 1.5 * pc_z[j] * pb_xx[j] * fl1_fx - 0.5 * pc_z[j] * fl1_fx * pb_zz[j] - 2.0 * fl1_fx * pb_xz[j] * pc_x[j] - fl1_fx * pb_xxz[j] - 2.0 * pa_z[j] * pb_xxz[j] * pc_z[j] - 2.0 * pa_z[j] * pb_xzz[j] * pc_x[j] - pc_z[j] * pb_xxzz[j]);

                t_z_xxzz[j] += fl_s_0_0_2 * (0.25 * pa_z[j] * fl2_fx + 1.5 * pc_z[j] * fl2_fx + 0.5 * fl2_fx * pb_z[j] + pa_z[j] * pb_x[j] * pc_x[j] * fl1_fx + 0.5 * pa_z[j] * pc_xx[j] * fl1_fx + 0.5 * pa_z[j] * fl1_fx * pc_zz[j] + pa_z[j] * fl1_fx * pb_z[j] * pc_z[j] + 1.5 * pc_z[j] * pb_xx[j] * fl1_fx + 3.0 * pc_xz[j] * pb_x[j] * fl1_fx + pc_zz[j] * fl1_fx * pb_z[j] + 0.5 * pc_z[j] * fl1_fx * pb_zz[j] + fl1_fx * pc_xx[j] * pb_z[j] + 2.0 * fl1_fx * pb_xz[j] * pc_x[j] + pa_z[j] * pb_xx[j] * pc_zz[j] + 4.0 * pa_z[j] * pb_xz[j] * pc_xz[j] + pa_z[j] * pc_xx[j] * pb_zz[j] + 2.0 * pc_zz[j] * pb_xxz[j] + 2.0 * pc_xz[j] * pb_xzz[j]);

                t_z_xxzz[j] += fl_s_0_0_3 * (-0.75 * pc_z[j] * fl2_fx - 0.5 * pa_z[j] * pc_xx[j] * fl1_fx - 0.5 * pa_z[j] * fl1_fx * pc_zz[j] - 3.0 * pc_xz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xxz[j] * fl1_fx - 0.5 * pc_zzz[j] * fl1_fx - pc_zz[j] * fl1_fx * pb_z[j] - fl1_fx * pc_xx[j] * pb_z[j] - 2.0 * pa_z[j] * pb_x[j] * pc_xzz[j] - 2.0 * pa_z[j] * pc_xxz[j] * pb_z[j] - pc_zzz[j] * pb_xx[j] - 4.0 * pc_xzz[j] * pb_xz[j] - pc_xxz[j] * pb_zz[j]);

                t_z_xxzz[j] += fl_s_0_0_4 * (1.5 * pc_xxz[j] * fl1_fx + 0.5 * pc_zzz[j] * fl1_fx + pa_z[j] * pc_xxzz[j] + 2.0 * pc_xzzz[j] * pb_x[j] + 2.0 * pc_xxzz[j] * pb_z[j]);

                t_z_xxzz[j] += -fl_s_0_0_5 * pc_xxzzz[j];

                t_z_xyyy[j] = fl_s_0_0_0 * (1.5 * pa_z[j] * pb_xy[j] * fl1_fx + pa_z[j] * pb_xyyy[j]);

                t_z_xyyy[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * pb_xy[j] * fl1_fx - 1.5 * pa_z[j] * pb_x[j] * pc_y[j] * fl1_fx - 1.5 * pa_z[j] * pc_x[j] * pb_y[j] * fl1_fx - 1.5 * pc_z[j] * pb_xy[j] * fl1_fx - 3.0 * pa_z[j] * pb_xyy[j] * pc_y[j] - pa_z[j] * pc_x[j] * pb_yyy[j] - pc_z[j] * pb_xyyy[j]);

                t_z_xyyy[j] += fl_s_0_0_2 * (1.5 * pa_z[j] * pb_x[j] * pc_y[j] * fl1_fx + 1.5 * pa_z[j] * pc_x[j] * pb_y[j] * fl1_fx + 1.5 * pa_z[j] * pc_xy[j] * fl1_fx + 1.5 * pc_z[j] * pb_xy[j] * fl1_fx + 1.5 * pc_yz[j] * pb_x[j] * fl1_fx + 1.5 * pc_xz[j] * pb_y[j] * fl1_fx + 3.0 * pa_z[j] * pb_xy[j] * pc_yy[j] + 3.0 * pa_z[j] * pc_xy[j] * pb_yy[j] + 3.0 * pc_yz[j] * pb_xyy[j] + pc_xz[j] * pb_yyy[j]);

                t_z_xyyy[j] += fl_s_0_0_3 * (-1.5 * pa_z[j] * pc_xy[j] * fl1_fx - 1.5 * pc_yz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xz[j] * pb_y[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - pa_z[j] * pb_x[j] * pc_yyy[j] - 3.0 * pa_z[j] * pc_xyy[j] * pb_y[j] - 3.0 * pc_yyz[j] * pb_xy[j] - 3.0 * pc_xyz[j] * pb_yy[j]);

                t_z_xyyy[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_z[j] * pc_xyyy[j] + pc_yyyz[j] * pb_x[j] + 3.0 * pc_xyyz[j] * pb_y[j]);

                t_z_xyyy[j] += -fl_s_0_0_5 * pc_xyyyz[j];

                t_z_xyyz[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_x[j] + 0.5 * pa_z[j] * pb_xz[j] * fl1_fx + 0.5 * fl1_fx * pb_xyy[j] + pa_z[j] * pb_xyyz[j]);

                t_z_xyyz[j] += fl_s_0_0_1 * (-0.5 * fl2_fx * pb_x[j] - 0.25 * fl2_fx * pc_x[j] - 0.5 * pa_z[j] * pb_x[j] * fl1_fx * pc_z[j] - 0.5 * pa_z[j] * pb_xz[j] * fl1_fx - 0.5 * pa_z[j] * pc_x[j] * fl1_fx * pb_z[j] - 0.5 * pc_z[j] * pb_xz[j] * fl1_fx - fl1_fx * pb_xy[j] * pc_y[j] - 0.5 * fl1_fx * pc_x[j] * pb_yy[j] - 0.5 * fl1_fx * pb_xyy[j] - pa_z[j] * pb_xyy[j] * pc_z[j] - 2.0 * pa_z[j] * pb_xyz[j] * pc_y[j] - pa_z[j] * pc_x[j] * pb_yyz[j] - pc_z[j] * pb_xyyz[j]);

                t_z_xyyz[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_x[j] + 0.25 * fl2_fx * pb_x[j] + 0.5 * pa_z[j] * pb_x[j] * fl1_fx * pc_z[j] + 0.5 * pa_z[j] * pc_xz[j] * fl1_fx + 0.5 * pa_z[j] * pc_x[j] * fl1_fx * pb_z[j] + 0.5 * pc_zz[j] * pb_x[j] * fl1_fx + 0.5 * pc_z[j] * pb_xz[j] * fl1_fx + 0.5 * pc_xz[j] * fl1_fx * pb_z[j] + 0.5 * fl1_fx * pb_x[j] * pc_yy[j] + fl1_fx * pc_xy[j] * pb_y[j] + fl1_fx * pb_xy[j] * pc_y[j] + 0.5 * fl1_fx * pc_x[j] * pb_yy[j] + 2.0 * pa_z[j] * pb_xy[j] * pc_yz[j] + pa_z[j] * pb_xz[j] * pc_yy[j] + pa_z[j] * pc_xz[j] * pb_yy[j] + 2.0 * pa_z[j] * pc_xy[j] * pb_yz[j] + pc_zz[j] * pb_xyy[j] + 2.0 * pc_yz[j] * pb_xyz[j] + pc_xz[j] * pb_yyz[j]);

                t_z_xyyz[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_x[j] - 0.5 * pa_z[j] * pc_xz[j] * fl1_fx - 0.5 * pc_zz[j] * pb_x[j] * fl1_fx - 0.5 * pc_xzz[j] * fl1_fx - 0.5 * pc_xz[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_xyy[j] - 0.5 * fl1_fx * pb_x[j] * pc_yy[j] - fl1_fx * pc_xy[j] * pb_y[j] - pa_z[j] * pb_x[j] * pc_yyz[j] - 2.0 * pa_z[j] * pc_xyz[j] * pb_y[j] - pa_z[j] * pc_xyy[j] * pb_z[j] - 2.0 * pc_yzz[j] * pb_xy[j] - pc_yyz[j] * pb_xz[j] - pc_xzz[j] * pb_yy[j] - 2.0 * pc_xyz[j] * pb_yz[j]);

                t_z_xyyz[j] += fl_s_0_0_4 * (0.5 * pc_xzz[j] * fl1_fx + 0.5 * fl1_fx * pc_xyy[j] + pa_z[j] * pc_xyyz[j] + pc_yyzz[j] * pb_x[j] + 2.0 * pc_xyzz[j] * pb_y[j] + pc_xyyz[j] * pb_z[j]);

                t_z_xyyz[j] += -fl_s_0_0_5 * pc_xyyzz[j];

                t_z_xyzz[j] = fl_s_0_0_0 * (0.5 * pa_z[j] * pb_xy[j] * fl1_fx + fl1_fx * pb_xyz[j] + pa_z[j] * pb_xyzz[j]);

                t_z_xyzz[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * pb_xy[j] * fl1_fx - 0.5 * pa_z[j] * pb_x[j] * pc_y[j] * fl1_fx - 0.5 * pa_z[j] * pc_x[j] * pb_y[j] * fl1_fx - 1.5 * pc_z[j] * pb_xy[j] * fl1_fx - fl1_fx * pb_xz[j] * pc_y[j] - fl1_fx * pc_x[j] * pb_yz[j] - fl1_fx * pb_xyz[j] - 2.0 * pa_z[j] * pb_xyz[j] * pc_z[j] - pa_z[j] * pb_xzz[j] * pc_y[j] - pa_z[j] * pc_x[j] * pb_yzz[j] - pc_z[j] * pb_xyzz[j]);

                t_z_xyzz[j] += fl_s_0_0_2 * (0.5 * pa_z[j] * pb_x[j] * pc_y[j] * fl1_fx + 0.5 * pa_z[j] * pc_x[j] * pb_y[j] * fl1_fx + 0.5 * pa_z[j] * pc_xy[j] * fl1_fx + 1.5 * pc_z[j] * pb_xy[j] * fl1_fx + 1.5 * pc_yz[j] * pb_x[j] * fl1_fx + 1.5 * pc_xz[j] * pb_y[j] * fl1_fx + fl1_fx * pc_xy[j] * pb_z[j] + fl1_fx * pb_xz[j] * pc_y[j] + fl1_fx * pc_x[j] * pb_yz[j] + pa_z[j] * pb_xy[j] * pc_zz[j] + 2.0 * pa_z[j] * pb_xz[j] * pc_yz[j] + 2.0 * pa_z[j] * pc_xz[j] * pb_yz[j] + pa_z[j] * pc_xy[j] * pb_zz[j] + 2.0 * pc_zz[j] * pb_xyz[j] + pc_yz[j] * pb_xzz[j] + pc_xz[j] * pb_yzz[j]);

                t_z_xyzz[j] += fl_s_0_0_3 * (-0.5 * pa_z[j] * pc_xy[j] * fl1_fx - 1.5 * pc_yz[j] * pb_x[j] * fl1_fx - 1.5 * pc_xz[j] * pb_y[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - fl1_fx * pc_xy[j] * pb_z[j] - pa_z[j] * pb_x[j] * pc_yzz[j] - pa_z[j] * pc_xzz[j] * pb_y[j] - 2.0 * pa_z[j] * pc_xyz[j] * pb_z[j] - pc_zzz[j] * pb_xy[j] - 2.0 * pc_yzz[j] * pb_xz[j] - 2.0 * pc_xzz[j] * pb_yz[j] - pc_xyz[j] * pb_zz[j]);

                t_z_xyzz[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_z[j] * pc_xyzz[j] + pc_yzzz[j] * pb_x[j] + pc_xzzz[j] * pb_y[j] + 2.0 * pc_xyzz[j] * pb_z[j]);

                t_z_xyzz[j] += -fl_s_0_0_5 * pc_xyzzz[j];

                t_z_xzzz[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_x[j] + 1.5 * pa_z[j] * pb_xz[j] * fl1_fx + 1.5 * fl1_fx * pb_xzz[j] + pa_z[j] * pb_xzzz[j]);

                t_z_xzzz[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_x[j] - 0.75 * fl2_fx * pc_x[j] - 1.5 * pa_z[j] * pb_xz[j] * fl1_fx - 1.5 * pa_z[j] * pb_x[j] * pc_z[j] * fl1_fx - 1.5 * pa_z[j] * pc_x[j] * pb_z[j] * fl1_fx - 4.5 * pc_z[j] * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * pc_x[j] * pb_zz[j] - 1.5 * fl1_fx * pb_xzz[j] - 3.0 * pa_z[j] * pb_xzz[j] * pc_z[j] - pa_z[j] * pc_x[j] * pb_zzz[j] - pc_z[j] * pb_xzzz[j]);

                t_z_xzzz[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_x[j] + 0.75 * fl2_fx * pb_x[j] + 1.5 * pa_z[j] * pb_x[j] * pc_z[j] * fl1_fx + 1.5 * pa_z[j] * pc_x[j] * pb_z[j] * fl1_fx + 1.5 * pa_z[j] * pc_xz[j] * fl1_fx + 4.5 * pc_z[j] * pb_xz[j] * fl1_fx + 3.0 * pc_zz[j] * pb_x[j] * fl1_fx + 4.5 * pc_xz[j] * pb_z[j] * fl1_fx + 1.5 * fl1_fx * pc_x[j] * pb_zz[j] + 3.0 * pa_z[j] * pb_xz[j] * pc_zz[j] + 3.0 * pa_z[j] * pc_xz[j] * pb_zz[j] + 3.0 * pc_zz[j] * pb_xzz[j] + pc_xz[j] * pb_zzz[j]);

                t_z_xzzz[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_x[j] - 1.5 * pa_z[j] * pc_xz[j] * fl1_fx - 3.0 * pc_zz[j] * pb_x[j] * fl1_fx - 4.5 * pc_xz[j] * pb_z[j] * fl1_fx - 3.0 * pc_xzz[j] * fl1_fx - pa_z[j] * pb_x[j] * pc_zzz[j] - 3.0 * pa_z[j] * pc_xzz[j] * pb_z[j] - 3.0 * pc_zzz[j] * pb_xz[j] - 3.0 * pc_xzz[j] * pb_zz[j]);

                t_z_xzzz[j] += fl_s_0_0_4 * (3.0 * pc_xzz[j] * fl1_fx + pa_z[j] * pc_xzzz[j] + pc_zzzz[j] * pb_x[j] + 3.0 * pc_xzzz[j] * pb_z[j]);

                t_z_xzzz[j] += -fl_s_0_0_5 * pc_xzzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForPG_40_45(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (40,45)

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

            auto pa_z = paDistances.data(3 * idx + 2);

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

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_yyyy = pcDistances.data(55 * idx + 29);

            auto pc_yyyz = pcDistances.data(55 * idx + 30);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            auto pc_yzzz = pcDistances.data(55 * idx + 32);

            auto pc_zzzz = pcDistances.data(55 * idx + 33);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_yyyyz = pcDistances.data(55 * idx + 50);

            auto pc_yyyzz = pcDistances.data(55 * idx + 51);

            auto pc_yyzzz = pcDistances.data(55 * idx + 52);

            auto pc_yzzzz = pcDistances.data(55 * idx + 53);

            auto pc_zzzzz = pcDistances.data(55 * idx + 54);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_z_yyyy = primBuffer.data(45 * idx + 40);

            auto t_z_yyyz = primBuffer.data(45 * idx + 41);

            auto t_z_yyzz = primBuffer.data(45 * idx + 42);

            auto t_z_yzzz = primBuffer.data(45 * idx + 43);

            auto t_z_zzzz = primBuffer.data(45 * idx + 44);

            // Batch of Integrals (40,45)

            #pragma omp simd aligned(fx, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, pc_y, pc_yy, pc_yyy, pc_yyyy, pc_yyyyz, pc_yyyz, \
                                     pc_yyyzz, pc_yyz, pc_yyzz, pc_yyzzz, pc_yz, pc_yzz, pc_yzzz, pc_yzzzz, pc_z, pc_zz, \
                                     pc_zzz, pc_zzzz, pc_zzzzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, \
                                     t_z_yyyy, t_z_yyyz, t_z_yyzz, t_z_yzzz, t_z_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_z_yyyy[j] = fl_s_0_0_0 * (0.75 * pa_z[j] * fl2_fx + 3.0 * pa_z[j] * pb_yy[j] * fl1_fx + pa_z[j] * pb_yyyy[j]);

                t_z_yyyy[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * fl2_fx - 0.75 * pc_z[j] * fl2_fx - 3.0 * pa_z[j] * pb_yy[j] * fl1_fx - 6.0 * pa_z[j] * pb_y[j] * pc_y[j] * fl1_fx - 3.0 * pc_z[j] * pb_yy[j] * fl1_fx - 4.0 * pa_z[j] * pb_yyy[j] * pc_y[j] - pc_z[j] * pb_yyyy[j]);

                t_z_yyyy[j] += fl_s_0_0_2 * (0.75 * pa_z[j] * fl2_fx + 1.5 * pc_z[j] * fl2_fx + 6.0 * pa_z[j] * pb_y[j] * pc_y[j] * fl1_fx + 3.0 * pa_z[j] * pc_yy[j] * fl1_fx + 3.0 * pc_z[j] * pb_yy[j] * fl1_fx + 6.0 * pc_yz[j] * pb_y[j] * fl1_fx + 6.0 * pa_z[j] * pb_yy[j] * pc_yy[j] + 4.0 * pc_yz[j] * pb_yyy[j]);

                t_z_yyyy[j] += fl_s_0_0_3 * (-0.75 * pc_z[j] * fl2_fx - 3.0 * pa_z[j] * pc_yy[j] * fl1_fx - 6.0 * pc_yz[j] * pb_y[j] * fl1_fx - 3.0 * pc_yyz[j] * fl1_fx - 4.0 * pa_z[j] * pb_y[j] * pc_yyy[j] - 6.0 * pc_yyz[j] * pb_yy[j]);

                t_z_yyyy[j] += fl_s_0_0_4 * (3.0 * pc_yyz[j] * fl1_fx + pa_z[j] * pc_yyyy[j] + 4.0 * pc_yyyz[j] * pb_y[j]);

                t_z_yyyy[j] += -fl_s_0_0_5 * pc_yyyyz[j];

                t_z_yyyz[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_y[j] + 1.5 * pa_z[j] * pb_yz[j] * fl1_fx + 0.5 * fl1_fx * pb_yyy[j] + pa_z[j] * pb_yyyz[j]);

                t_z_yyyz[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_y[j] - 0.75 * fl2_fx * pc_y[j] - 1.5 * pa_z[j] * pb_y[j] * fl1_fx * pc_z[j] - 1.5 * pa_z[j] * pb_yz[j] * fl1_fx - 1.5 * pa_z[j] * pc_y[j] * fl1_fx * pb_z[j] - 1.5 * pc_z[j] * pb_yz[j] * fl1_fx - 1.5 * fl1_fx * pb_yy[j] * pc_y[j] - 0.5 * fl1_fx * pb_yyy[j] - pa_z[j] * pb_yyy[j] * pc_z[j] - 3.0 * pa_z[j] * pb_yyz[j] * pc_y[j] - pc_z[j] * pb_yyyz[j]);

                t_z_yyyz[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.75 * fl2_fx * pb_y[j] + 1.5 * pa_z[j] * pb_y[j] * fl1_fx * pc_z[j] + 1.5 * pa_z[j] * pc_yz[j] * fl1_fx + 1.5 * pa_z[j] * pc_y[j] * fl1_fx * pb_z[j] + 1.5 * pc_zz[j] * pb_y[j] * fl1_fx + 1.5 * pc_z[j] * pb_yz[j] * fl1_fx + 1.5 * pc_yz[j] * fl1_fx * pb_z[j] + 1.5 * fl1_fx * pb_y[j] * pc_yy[j] + 1.5 * fl1_fx * pb_yy[j] * pc_y[j] + 3.0 * pa_z[j] * pb_yy[j] * pc_yz[j] + 3.0 * pa_z[j] * pb_yz[j] * pc_yy[j] + pc_zz[j] * pb_yyy[j] + 3.0 * pc_yz[j] * pb_yyz[j]);

                t_z_yyyz[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 1.5 * pa_z[j] * pc_yz[j] * fl1_fx - 1.5 * pc_zz[j] * pb_y[j] * fl1_fx - 1.5 * pc_yzz[j] * fl1_fx - 1.5 * pc_yz[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_yyy[j] - 1.5 * fl1_fx * pb_y[j] * pc_yy[j] - 3.0 * pa_z[j] * pb_y[j] * pc_yyz[j] - pa_z[j] * pc_yyy[j] * pb_z[j] - 3.0 * pc_yzz[j] * pb_yy[j] - 3.0 * pc_yyz[j] * pb_yz[j]);

                t_z_yyyz[j] += fl_s_0_0_4 * (1.5 * pc_yzz[j] * fl1_fx + 0.5 * fl1_fx * pc_yyy[j] + pa_z[j] * pc_yyyz[j] + 3.0 * pc_yyzz[j] * pb_y[j] + pc_yyyz[j] * pb_z[j]);

                t_z_yyyz[j] += -fl_s_0_0_5 * pc_yyyzz[j];

                t_z_yyzz[j] = fl_s_0_0_0 * (0.25 * pa_z[j] * fl2_fx + 0.5 * fl2_fx * pb_z[j] + 0.5 * pa_z[j] * pb_yy[j] * fl1_fx + 0.5 * pa_z[j] * fl1_fx * pb_zz[j] + fl1_fx * pb_yyz[j] + pa_z[j] * pb_yyzz[j]);

                t_z_yyzz[j] += fl_s_0_0_1 * (-0.5 * pa_z[j] * fl2_fx - 0.75 * pc_z[j] * fl2_fx - fl2_fx * pb_z[j] - 0.5 * pa_z[j] * pb_yy[j] * fl1_fx - pa_z[j] * pb_y[j] * pc_y[j] * fl1_fx - pa_z[j] * fl1_fx * pb_z[j] * pc_z[j] - 0.5 * pa_z[j] * fl1_fx * pb_zz[j] - 1.5 * pc_z[j] * pb_yy[j] * fl1_fx - 0.5 * pc_z[j] * fl1_fx * pb_zz[j] - 2.0 * fl1_fx * pb_yz[j] * pc_y[j] - fl1_fx * pb_yyz[j] - 2.0 * pa_z[j] * pb_yyz[j] * pc_z[j] - 2.0 * pa_z[j] * pb_yzz[j] * pc_y[j] - pc_z[j] * pb_yyzz[j]);

                t_z_yyzz[j] += fl_s_0_0_2 * (0.25 * pa_z[j] * fl2_fx + 1.5 * pc_z[j] * fl2_fx + 0.5 * fl2_fx * pb_z[j] + pa_z[j] * pb_y[j] * pc_y[j] * fl1_fx + 0.5 * pa_z[j] * pc_yy[j] * fl1_fx + 0.5 * pa_z[j] * fl1_fx * pc_zz[j] + pa_z[j] * fl1_fx * pb_z[j] * pc_z[j] + 1.5 * pc_z[j] * pb_yy[j] * fl1_fx + 3.0 * pc_yz[j] * pb_y[j] * fl1_fx + pc_zz[j] * fl1_fx * pb_z[j] + 0.5 * pc_z[j] * fl1_fx * pb_zz[j] + fl1_fx * pc_yy[j] * pb_z[j] + 2.0 * fl1_fx * pb_yz[j] * pc_y[j] + pa_z[j] * pb_yy[j] * pc_zz[j] + 4.0 * pa_z[j] * pb_yz[j] * pc_yz[j] + pa_z[j] * pc_yy[j] * pb_zz[j] + 2.0 * pc_zz[j] * pb_yyz[j] + 2.0 * pc_yz[j] * pb_yzz[j]);

                t_z_yyzz[j] += fl_s_0_0_3 * (-0.75 * pc_z[j] * fl2_fx - 0.5 * pa_z[j] * pc_yy[j] * fl1_fx - 0.5 * pa_z[j] * fl1_fx * pc_zz[j] - 3.0 * pc_yz[j] * pb_y[j] * fl1_fx - 1.5 * pc_yyz[j] * fl1_fx - 0.5 * pc_zzz[j] * fl1_fx - pc_zz[j] * fl1_fx * pb_z[j] - fl1_fx * pc_yy[j] * pb_z[j] - 2.0 * pa_z[j] * pb_y[j] * pc_yzz[j] - 2.0 * pa_z[j] * pc_yyz[j] * pb_z[j] - pc_zzz[j] * pb_yy[j] - 4.0 * pc_yzz[j] * pb_yz[j] - pc_yyz[j] * pb_zz[j]);

                t_z_yyzz[j] += fl_s_0_0_4 * (1.5 * pc_yyz[j] * fl1_fx + 0.5 * pc_zzz[j] * fl1_fx + pa_z[j] * pc_yyzz[j] + 2.0 * pc_yzzz[j] * pb_y[j] + 2.0 * pc_yyzz[j] * pb_z[j]);

                t_z_yyzz[j] += -fl_s_0_0_5 * pc_yyzzz[j];

                t_z_yzzz[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_y[j] + 1.5 * pa_z[j] * pb_yz[j] * fl1_fx + 1.5 * fl1_fx * pb_yzz[j] + pa_z[j] * pb_yzzz[j]);

                t_z_yzzz[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pb_y[j] - 0.75 * fl2_fx * pc_y[j] - 1.5 * pa_z[j] * pb_yz[j] * fl1_fx - 1.5 * pa_z[j] * pb_y[j] * pc_z[j] * fl1_fx - 1.5 * pa_z[j] * pc_y[j] * pb_z[j] * fl1_fx - 4.5 * pc_z[j] * pb_yz[j] * fl1_fx - 1.5 * fl1_fx * pc_y[j] * pb_zz[j] - 1.5 * fl1_fx * pb_yzz[j] - 3.0 * pa_z[j] * pb_yzz[j] * pc_z[j] - pa_z[j] * pc_y[j] * pb_zzz[j] - pc_z[j] * pb_yzzz[j]);

                t_z_yzzz[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.75 * fl2_fx * pb_y[j] + 1.5 * pa_z[j] * pb_y[j] * pc_z[j] * fl1_fx + 1.5 * pa_z[j] * pc_y[j] * pb_z[j] * fl1_fx + 1.5 * pa_z[j] * pc_yz[j] * fl1_fx + 4.5 * pc_z[j] * pb_yz[j] * fl1_fx + 3.0 * pc_zz[j] * pb_y[j] * fl1_fx + 4.5 * pc_yz[j] * pb_z[j] * fl1_fx + 1.5 * fl1_fx * pc_y[j] * pb_zz[j] + 3.0 * pa_z[j] * pb_yz[j] * pc_zz[j] + 3.0 * pa_z[j] * pc_yz[j] * pb_zz[j] + 3.0 * pc_zz[j] * pb_yzz[j] + pc_yz[j] * pb_zzz[j]);

                t_z_yzzz[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 1.5 * pa_z[j] * pc_yz[j] * fl1_fx - 3.0 * pc_zz[j] * pb_y[j] * fl1_fx - 4.5 * pc_yz[j] * pb_z[j] * fl1_fx - 3.0 * pc_yzz[j] * fl1_fx - pa_z[j] * pb_y[j] * pc_zzz[j] - 3.0 * pa_z[j] * pc_yzz[j] * pb_z[j] - 3.0 * pc_zzz[j] * pb_yz[j] - 3.0 * pc_yzz[j] * pb_zz[j]);

                t_z_yzzz[j] += fl_s_0_0_4 * (3.0 * pc_yzz[j] * fl1_fx + pa_z[j] * pc_yzzz[j] + pc_zzzz[j] * pb_y[j] + 3.0 * pc_yzzz[j] * pb_z[j]);

                t_z_yzzz[j] += -fl_s_0_0_5 * pc_yzzzz[j];

                t_z_zzzz[j] = fl_s_0_0_0 * (0.75 * pa_z[j] * fl2_fx + 3.0 * fl2_fx * pb_z[j] + 3.0 * pa_z[j] * pb_zz[j] * fl1_fx + 2.0 * fl1_fx * pb_zzz[j] + pa_z[j] * pb_zzzz[j]);

                t_z_zzzz[j] += fl_s_0_0_1 * (-1.5 * pa_z[j] * fl2_fx - 3.75 * pc_z[j] * fl2_fx - 6.0 * fl2_fx * pb_z[j] - 3.0 * pa_z[j] * pb_zz[j] * fl1_fx - 6.0 * pa_z[j] * pb_z[j] * pc_z[j] * fl1_fx - 9.0 * pc_z[j] * pb_zz[j] * fl1_fx - 2.0 * fl1_fx * pb_zzz[j] - 4.0 * pa_z[j] * pb_zzz[j] * pc_z[j] - pc_z[j] * pb_zzzz[j]);

                t_z_zzzz[j] += fl_s_0_0_2 * (0.75 * pa_z[j] * fl2_fx + 7.5 * pc_z[j] * fl2_fx + 3.0 * fl2_fx * pb_z[j] + 6.0 * pa_z[j] * pb_z[j] * pc_z[j] * fl1_fx + 3.0 * pa_z[j] * pc_zz[j] * fl1_fx + 9.0 * pc_z[j] * pb_zz[j] * fl1_fx + 12.0 * pc_zz[j] * pb_z[j] * fl1_fx + 6.0 * pa_z[j] * pb_zz[j] * pc_zz[j] + 4.0 * pc_zz[j] * pb_zzz[j]);

                t_z_zzzz[j] += fl_s_0_0_3 * (-3.75 * pc_z[j] * fl2_fx - 3.0 * pa_z[j] * pc_zz[j] * fl1_fx - 12.0 * pc_zz[j] * pb_z[j] * fl1_fx - 5.0 * pc_zzz[j] * fl1_fx - 4.0 * pa_z[j] * pb_z[j] * pc_zzz[j] - 6.0 * pc_zzz[j] * pb_zz[j]);

                t_z_zzzz[j] += fl_s_0_0_4 * (5.0 * pc_zzz[j] * fl1_fx + pa_z[j] * pc_zzzz[j] + 4.0 * pc_zzzz[j] * pb_z[j]);

                t_z_zzzz[j] += -fl_s_0_0_5 * pc_zzzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP(      CMemBlock2D<double>& primBuffer,
                              const CMemBlock2D<double>& auxBuffer,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pbDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForGP_0_5(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGP_5_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                    braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGP_10_15(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGP_15_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGP_20_25(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGP_25_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGP_30_35(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGP_35_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 

        npotrecfunc::compNuclearPotentialForGP_40_45(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pcDistances, 
                                                     braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compNuclearPotentialForGP_0_5(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
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

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(55 * idx + 9);

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(55 * idx + 19);

            auto pc_xxxy = pcDistances.data(55 * idx + 20);

            auto pc_xxxz = pcDistances.data(55 * idx + 21);

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxxxx = pcDistances.data(55 * idx + 34);

            auto pc_xxxxy = pcDistances.data(55 * idx + 35);

            auto pc_xxxxz = pcDistances.data(55 * idx + 36);

            auto pc_xxxyy = pcDistances.data(55 * idx + 37);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_xxxx_x = primBuffer.data(45 * idx);

            auto t_xxxx_y = primBuffer.data(45 * idx + 1);

            auto t_xxxx_z = primBuffer.data(45 * idx + 2);

            auto t_xxxy_x = primBuffer.data(45 * idx + 3);

            auto t_xxxy_y = primBuffer.data(45 * idx + 4);

            // Batch of Integrals (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_y, pb_z, \
                                     pc_x, pc_xx, pc_xxx, pc_xxxx, pc_xxxxx, pc_xxxxy, pc_xxxxz, pc_xxxy, pc_xxxyy, \
                                     pc_xxxz, pc_xxy, pc_xxyy, pc_xxz, pc_xy, pc_xyy, pc_xz, pc_y, pc_yy, pc_z, s_0_0_0, \
                                     s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_xxxx_x, t_xxxx_y, t_xxxx_z, t_xxxy_x, \
                                     t_xxxy_y: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxxx_x[j] = fl_s_0_0_0 * (3.0 * pa_x[j] * fl2_fx + 2.0 * pa_xxx[j] * fl1_fx + 0.75 * fl2_fx * pb_x[j] + 3.0 * pa_xx[j] * fl1_fx * pb_x[j] + pa_xxxx[j] * pb_x[j]);

                t_xxxx_x[j] += fl_s_0_0_1 * (-6.0 * pa_x[j] * fl2_fx - 3.75 * pc_x[j] * fl2_fx - 2.0 * pa_xxx[j] * fl1_fx - 9.0 * pa_xx[j] * pc_x[j] * fl1_fx - 1.5 * fl2_fx * pb_x[j] - 3.0 * pa_xx[j] * fl1_fx * pb_x[j] - 6.0 * pa_x[j] * pc_x[j] * fl1_fx * pb_x[j] - pa_xxxx[j] * pc_x[j] - 4.0 * pa_xxx[j] * pc_x[j] * pb_x[j]);

                t_xxxx_x[j] += fl_s_0_0_2 * (3.0 * pa_x[j] * fl2_fx + 7.5 * pc_x[j] * fl2_fx + 9.0 * pa_xx[j] * pc_x[j] * fl1_fx + 12.0 * pa_x[j] * pc_xx[j] * fl1_fx + 0.75 * fl2_fx * pb_x[j] + 6.0 * pa_x[j] * pc_x[j] * fl1_fx * pb_x[j] + 3.0 * pc_xx[j] * fl1_fx * pb_x[j] + 4.0 * pa_xxx[j] * pc_xx[j] + 6.0 * pa_xx[j] * pc_xx[j] * pb_x[j]);

                t_xxxx_x[j] += fl_s_0_0_3 * (-3.75 * pc_x[j] * fl2_fx - 12.0 * pa_x[j] * pc_xx[j] * fl1_fx - 5.0 * pc_xxx[j] * fl1_fx - 3.0 * pc_xx[j] * fl1_fx * pb_x[j] - 6.0 * pa_xx[j] * pc_xxx[j] - 4.0 * pa_x[j] * pc_xxx[j] * pb_x[j]);

                t_xxxx_x[j] += fl_s_0_0_4 * (5.0 * pc_xxx[j] * fl1_fx + 4.0 * pa_x[j] * pc_xxxx[j] + pc_xxxx[j] * pb_x[j]);

                t_xxxx_x[j] += -fl_s_0_0_5 * pc_xxxxx[j];

                t_xxxx_y[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_y[j] + 3.0 * pa_xx[j] * fl1_fx * pb_y[j] + pa_xxxx[j] * pb_y[j]);

                t_xxxx_y[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_y[j] - 1.5 * fl2_fx * pb_y[j] - 3.0 * pa_xx[j] * fl1_fx * pc_y[j] - 3.0 * pa_xx[j] * fl1_fx * pb_y[j] - 6.0 * pa_x[j] * pc_x[j] * fl1_fx * pb_y[j] - pa_xxxx[j] * pc_y[j] - 4.0 * pa_xxx[j] * pc_x[j] * pb_y[j]);

                t_xxxx_y[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.75 * fl2_fx * pb_y[j] + 3.0 * pa_xx[j] * fl1_fx * pc_y[j] + 6.0 * pa_x[j] * pc_xy[j] * fl1_fx + 6.0 * pa_x[j] * pc_x[j] * fl1_fx * pb_y[j] + 3.0 * pc_xx[j] * fl1_fx * pb_y[j] + 4.0 * pa_xxx[j] * pc_xy[j] + 6.0 * pa_xx[j] * pc_xx[j] * pb_y[j]);

                t_xxxx_y[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 6.0 * pa_x[j] * pc_xy[j] * fl1_fx - 3.0 * pc_xxy[j] * fl1_fx - 3.0 * pc_xx[j] * fl1_fx * pb_y[j] - 6.0 * pa_xx[j] * pc_xxy[j] - 4.0 * pa_x[j] * pc_xxx[j] * pb_y[j]);

                t_xxxx_y[j] += fl_s_0_0_4 * (3.0 * pc_xxy[j] * fl1_fx + 4.0 * pa_x[j] * pc_xxxy[j] + pc_xxxx[j] * pb_y[j]);

                t_xxxx_y[j] += -fl_s_0_0_5 * pc_xxxxy[j];

                t_xxxx_z[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_z[j] + 3.0 * pa_xx[j] * fl1_fx * pb_z[j] + pa_xxxx[j] * pb_z[j]);

                t_xxxx_z[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * fl2_fx * pb_z[j] - 3.0 * pa_xx[j] * fl1_fx * pc_z[j] - 3.0 * pa_xx[j] * fl1_fx * pb_z[j] - 6.0 * pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] - pa_xxxx[j] * pc_z[j] - 4.0 * pa_xxx[j] * pc_x[j] * pb_z[j]);

                t_xxxx_z[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pb_z[j] + 3.0 * pa_xx[j] * fl1_fx * pc_z[j] + 6.0 * pa_x[j] * pc_xz[j] * fl1_fx + 6.0 * pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] + 3.0 * pc_xx[j] * fl1_fx * pb_z[j] + 4.0 * pa_xxx[j] * pc_xz[j] + 6.0 * pa_xx[j] * pc_xx[j] * pb_z[j]);

                t_xxxx_z[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 6.0 * pa_x[j] * pc_xz[j] * fl1_fx - 3.0 * pc_xxz[j] * fl1_fx - 3.0 * pc_xx[j] * fl1_fx * pb_z[j] - 6.0 * pa_xx[j] * pc_xxz[j] - 4.0 * pa_x[j] * pc_xxx[j] * pb_z[j]);

                t_xxxx_z[j] += fl_s_0_0_4 * (3.0 * pc_xxz[j] * fl1_fx + 4.0 * pa_x[j] * pc_xxxz[j] + pc_xxxx[j] * pb_z[j]);

                t_xxxx_z[j] += -fl_s_0_0_5 * pc_xxxxz[j];

                t_xxxy_x[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pa_y[j] + 1.5 * pa_xxy[j] * fl1_fx + 1.5 * pa_xy[j] * fl1_fx * pb_x[j] + pa_xxxy[j] * pb_x[j]);

                t_xxxy_x[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_y[j] - 1.5 * fl2_fx * pa_y[j] - 1.5 * pa_xx[j] * fl1_fx * pc_y[j] - 1.5 * pa_xxy[j] * fl1_fx - 4.5 * pa_xy[j] * pc_x[j] * fl1_fx - 1.5 * pa_x[j] * fl1_fx * pc_y[j] * pb_x[j] - 1.5 * pa_xy[j] * fl1_fx * pb_x[j] - 1.5 * pc_x[j] * fl1_fx * pa_y[j] * pb_x[j] - pa_xxxy[j] * pc_x[j] - pa_xxx[j] * pc_y[j] * pb_x[j] - 3.0 * pa_xxy[j] * pc_x[j] * pb_x[j]);

                t_xxxy_x[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.75 * fl2_fx * pa_y[j] + 1.5 * pa_xx[j] * fl1_fx * pc_y[j] + 4.5 * pa_x[j] * pc_xy[j] * fl1_fx + 4.5 * pa_xy[j] * pc_x[j] * fl1_fx + 3.0 * pc_xx[j] * fl1_fx * pa_y[j] + 1.5 * pa_x[j] * fl1_fx * pc_y[j] * pb_x[j] + 1.5 * pc_xy[j] * fl1_fx * pb_x[j] + 1.5 * pc_x[j] * fl1_fx * pa_y[j] * pb_x[j] + pa_xxx[j] * pc_xy[j] + 3.0 * pa_xxy[j] * pc_xx[j] + 3.0 * pa_xx[j] * pc_xy[j] * pb_x[j] + 3.0 * pa_xy[j] * pc_xx[j] * pb_x[j]);

                t_xxxy_x[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 4.5 * pa_x[j] * pc_xy[j] * fl1_fx - 3.0 * pc_xxy[j] * fl1_fx - 3.0 * pc_xx[j] * fl1_fx * pa_y[j] - 1.5 * pc_xy[j] * fl1_fx * pb_x[j] - 3.0 * pa_xx[j] * pc_xxy[j] - 3.0 * pa_xy[j] * pc_xxx[j] - 3.0 * pa_x[j] * pc_xxy[j] * pb_x[j] - pc_xxx[j] * pa_y[j] * pb_x[j]);

                t_xxxy_x[j] += fl_s_0_0_4 * (3.0 * pc_xxy[j] * fl1_fx + 3.0 * pa_x[j] * pc_xxxy[j] + pc_xxxx[j] * pa_y[j] + pc_xxxy[j] * pb_x[j]);

                t_xxxy_x[j] += -fl_s_0_0_5 * pc_xxxxy[j];

                t_xxxy_y[j] = fl_s_0_0_0 * (0.75 * pa_x[j] * fl2_fx + 0.5 * pa_xxx[j] * fl1_fx + 1.5 * pa_xy[j] * fl1_fx * pb_y[j] + pa_xxxy[j] * pb_y[j]);

                t_xxxy_y[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - 0.5 * pa_xxx[j] * fl1_fx - 1.5 * pa_xx[j] * pc_x[j] * fl1_fx - 1.5 * pa_xy[j] * fl1_fx * pc_y[j] - 1.5 * pa_x[j] * fl1_fx * pc_y[j] * pb_y[j] - 1.5 * pa_xy[j] * fl1_fx * pb_y[j] - 1.5 * pc_x[j] * fl1_fx * pa_y[j] * pb_y[j] - pa_xxxy[j] * pc_y[j] - pa_xxx[j] * pc_y[j] * pb_y[j] - 3.0 * pa_xxy[j] * pc_x[j] * pb_y[j]);

                t_xxxy_y[j] += fl_s_0_0_2 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + 1.5 * pa_xx[j] * pc_x[j] * fl1_fx + 1.5 * pa_x[j] * pc_xx[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pc_yy[j] + 1.5 * pa_xy[j] * fl1_fx * pc_y[j] + 1.5 * pa_x[j] * fl1_fx * pc_y[j] * pb_y[j] + 1.5 * pc_xy[j] * fl1_fx * pa_y[j] + 1.5 * pc_xy[j] * fl1_fx * pb_y[j] + 1.5 * pc_x[j] * fl1_fx * pa_y[j] * pb_y[j] + pa_xxx[j] * pc_yy[j] + 3.0 * pa_xxy[j] * pc_xy[j] + 3.0 * pa_xx[j] * pc_xy[j] * pb_y[j] + 3.0 * pa_xy[j] * pc_xx[j] * pb_y[j]);

                t_xxxy_y[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - 1.5 * pa_x[j] * pc_xx[j] * fl1_fx - 0.5 * pc_xxx[j] * fl1_fx - 1.5 * pa_x[j] * fl1_fx * pc_yy[j] - 1.5 * pc_xyy[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pa_y[j] - 1.5 * pc_xy[j] * fl1_fx * pb_y[j] - 3.0 * pa_xx[j] * pc_xyy[j] - 3.0 * pa_xy[j] * pc_xxy[j] - 3.0 * pa_x[j] * pc_xxy[j] * pb_y[j] - pc_xxx[j] * pa_y[j] * pb_y[j]);

                t_xxxy_y[j] += fl_s_0_0_4 * (0.5 * pc_xxx[j] * fl1_fx + 1.5 * pc_xyy[j] * fl1_fx + 3.0 * pa_x[j] * pc_xxyy[j] + pc_xxxy[j] * pa_y[j] + pc_xxxy[j] * pb_y[j]);

                t_xxxy_y[j] += -fl_s_0_0_5 * pc_xxxyy[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP_5_10(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(55 * idx + 9);

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxx = pcDistances.data(55 * idx + 19);

            auto pc_xxxy = pcDistances.data(55 * idx + 20);

            auto pc_xxxz = pcDistances.data(55 * idx + 21);

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxxxz = pcDistances.data(55 * idx + 36);

            auto pc_xxxyy = pcDistances.data(55 * idx + 37);

            auto pc_xxxyz = pcDistances.data(55 * idx + 38);

            auto pc_xxxzz = pcDistances.data(55 * idx + 39);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_xxxy_z = primBuffer.data(45 * idx + 5);

            auto t_xxxz_x = primBuffer.data(45 * idx + 6);

            auto t_xxxz_y = primBuffer.data(45 * idx + 7);

            auto t_xxxz_z = primBuffer.data(45 * idx + 8);

            auto t_xxyy_x = primBuffer.data(45 * idx + 9);

            // Batch of Integrals (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxxz, pa_xxy, pa_xxyy, pa_xxz, pa_xy, \
                                     pa_xyy, pa_xz, pa_y, pa_yy, pa_z, pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xxx, pc_xxxx, pc_xxxxz, \
                                     pc_xxxy, pc_xxxyy, pc_xxxyz, pc_xxxz, pc_xxxzz, pc_xxy, pc_xxyy, pc_xxyz, pc_xxz, \
                                     pc_xxzz, pc_xy, pc_xyy, pc_xyz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yz, pc_z, pc_zz, s_0_0_0, \
                                     s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_xxxy_z, t_xxxz_x, t_xxxz_y, t_xxxz_z, \
                                     t_xxyy_x: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxxy_z[j] = fl_s_0_0_0 * (1.5 * pa_xy[j] * fl1_fx * pb_z[j] + pa_xxxy[j] * pb_z[j]);

                t_xxxy_z[j] += fl_s_0_0_1 * (-1.5 * pa_xy[j] * fl1_fx * pc_z[j] - 1.5 * pa_x[j] * fl1_fx * pc_y[j] * pb_z[j] - 1.5 * pa_xy[j] * fl1_fx * pb_z[j] - 1.5 * pc_x[j] * fl1_fx * pa_y[j] * pb_z[j] - pa_xxxy[j] * pc_z[j] - pa_xxx[j] * pc_y[j] * pb_z[j] - 3.0 * pa_xxy[j] * pc_x[j] * pb_z[j]);

                t_xxxy_z[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * fl1_fx * pc_yz[j] + 1.5 * pa_xy[j] * fl1_fx * pc_z[j] + 1.5 * pa_x[j] * fl1_fx * pc_y[j] * pb_z[j] + 1.5 * pc_xz[j] * fl1_fx * pa_y[j] + 1.5 * pc_xy[j] * fl1_fx * pb_z[j] + 1.5 * pc_x[j] * fl1_fx * pa_y[j] * pb_z[j] + pa_xxx[j] * pc_yz[j] + 3.0 * pa_xxy[j] * pc_xz[j] + 3.0 * pa_xx[j] * pc_xy[j] * pb_z[j] + 3.0 * pa_xy[j] * pc_xx[j] * pb_z[j]);

                t_xxxy_z[j] += fl_s_0_0_3 * (-1.5 * pa_x[j] * fl1_fx * pc_yz[j] - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pa_y[j] - 1.5 * pc_xy[j] * fl1_fx * pb_z[j] - 3.0 * pa_xx[j] * pc_xyz[j] - 3.0 * pa_xy[j] * pc_xxz[j] - 3.0 * pa_x[j] * pc_xxy[j] * pb_z[j] - pc_xxx[j] * pa_y[j] * pb_z[j]);

                t_xxxy_z[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + 3.0 * pa_x[j] * pc_xxyz[j] + pc_xxxz[j] * pa_y[j] + pc_xxxy[j] * pb_z[j]);

                t_xxxy_z[j] += -fl_s_0_0_5 * pc_xxxyz[j];

                t_xxxz_x[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pa_z[j] + 1.5 * pa_xxz[j] * fl1_fx + 1.5 * pa_xz[j] * fl1_fx * pb_x[j] + pa_xxxz[j] * pb_x[j]);

                t_xxxz_x[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * fl2_fx * pa_z[j] - 1.5 * pa_xx[j] * fl1_fx * pc_z[j] - 1.5 * pa_xxz[j] * fl1_fx - 4.5 * pa_xz[j] * pc_x[j] * fl1_fx - 1.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_x[j] - 1.5 * pa_xz[j] * fl1_fx * pb_x[j] - 1.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_x[j] - pa_xxxz[j] * pc_x[j] - pa_xxx[j] * pc_z[j] * pb_x[j] - 3.0 * pa_xxz[j] * pc_x[j] * pb_x[j]);

                t_xxxz_x[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pa_z[j] + 1.5 * pa_xx[j] * fl1_fx * pc_z[j] + 4.5 * pa_x[j] * pc_xz[j] * fl1_fx + 4.5 * pa_xz[j] * pc_x[j] * fl1_fx + 3.0 * pc_xx[j] * fl1_fx * pa_z[j] + 1.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_x[j] + 1.5 * pc_xz[j] * fl1_fx * pb_x[j] + 1.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_x[j] + pa_xxx[j] * pc_xz[j] + 3.0 * pa_xxz[j] * pc_xx[j] + 3.0 * pa_xx[j] * pc_xz[j] * pb_x[j] + 3.0 * pa_xz[j] * pc_xx[j] * pb_x[j]);

                t_xxxz_x[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 4.5 * pa_x[j] * pc_xz[j] * fl1_fx - 3.0 * pc_xxz[j] * fl1_fx - 3.0 * pc_xx[j] * fl1_fx * pa_z[j] - 1.5 * pc_xz[j] * fl1_fx * pb_x[j] - 3.0 * pa_xx[j] * pc_xxz[j] - 3.0 * pa_xz[j] * pc_xxx[j] - 3.0 * pa_x[j] * pc_xxz[j] * pb_x[j] - pc_xxx[j] * pa_z[j] * pb_x[j]);

                t_xxxz_x[j] += fl_s_0_0_4 * (3.0 * pc_xxz[j] * fl1_fx + 3.0 * pa_x[j] * pc_xxxz[j] + pc_xxxx[j] * pa_z[j] + pc_xxxz[j] * pb_x[j]);

                t_xxxz_x[j] += -fl_s_0_0_5 * pc_xxxxz[j];

                t_xxxz_y[j] = fl_s_0_0_0 * (1.5 * pa_xz[j] * fl1_fx * pb_y[j] + pa_xxxz[j] * pb_y[j]);

                t_xxxz_y[j] += fl_s_0_0_1 * (-1.5 * pa_xz[j] * fl1_fx * pc_y[j] - 1.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_y[j] - 1.5 * pa_xz[j] * fl1_fx * pb_y[j] - 1.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_y[j] - pa_xxxz[j] * pc_y[j] - pa_xxx[j] * pc_z[j] * pb_y[j] - 3.0 * pa_xxz[j] * pc_x[j] * pb_y[j]);

                t_xxxz_y[j] += fl_s_0_0_2 * (1.5 * pa_x[j] * fl1_fx * pc_yz[j] + 1.5 * pa_xz[j] * fl1_fx * pc_y[j] + 1.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_y[j] + 1.5 * pc_xy[j] * fl1_fx * pa_z[j] + 1.5 * pc_xz[j] * fl1_fx * pb_y[j] + 1.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_y[j] + pa_xxx[j] * pc_yz[j] + 3.0 * pa_xxz[j] * pc_xy[j] + 3.0 * pa_xx[j] * pc_xz[j] * pb_y[j] + 3.0 * pa_xz[j] * pc_xx[j] * pb_y[j]);

                t_xxxz_y[j] += fl_s_0_0_3 * (-1.5 * pa_x[j] * fl1_fx * pc_yz[j] - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pa_z[j] - 1.5 * pc_xz[j] * fl1_fx * pb_y[j] - 3.0 * pa_xx[j] * pc_xyz[j] - 3.0 * pa_xz[j] * pc_xxy[j] - 3.0 * pa_x[j] * pc_xxz[j] * pb_y[j] - pc_xxx[j] * pa_z[j] * pb_y[j]);

                t_xxxz_y[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + 3.0 * pa_x[j] * pc_xxyz[j] + pc_xxxy[j] * pa_z[j] + pc_xxxz[j] * pb_y[j]);

                t_xxxz_y[j] += -fl_s_0_0_5 * pc_xxxyz[j];

                t_xxxz_z[j] = fl_s_0_0_0 * (0.75 * pa_x[j] * fl2_fx + 0.5 * pa_xxx[j] * fl1_fx + 1.5 * pa_xz[j] * fl1_fx * pb_z[j] + pa_xxxz[j] * pb_z[j]);

                t_xxxz_z[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - 0.5 * pa_xxx[j] * fl1_fx - 1.5 * pa_xx[j] * pc_x[j] * fl1_fx - 1.5 * pa_xz[j] * fl1_fx * pc_z[j] - 1.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_z[j] - 1.5 * pa_xz[j] * fl1_fx * pb_z[j] - 1.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_z[j] - pa_xxxz[j] * pc_z[j] - pa_xxx[j] * pc_z[j] * pb_z[j] - 3.0 * pa_xxz[j] * pc_x[j] * pb_z[j]);

                t_xxxz_z[j] += fl_s_0_0_2 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + 1.5 * pa_xx[j] * pc_x[j] * fl1_fx + 1.5 * pa_x[j] * pc_xx[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pc_zz[j] + 1.5 * pa_xz[j] * fl1_fx * pc_z[j] + 1.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_z[j] + 1.5 * pc_xz[j] * fl1_fx * pa_z[j] + 1.5 * pc_xz[j] * fl1_fx * pb_z[j] + 1.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_z[j] + pa_xxx[j] * pc_zz[j] + 3.0 * pa_xxz[j] * pc_xz[j] + 3.0 * pa_xx[j] * pc_xz[j] * pb_z[j] + 3.0 * pa_xz[j] * pc_xx[j] * pb_z[j]);

                t_xxxz_z[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - 1.5 * pa_x[j] * pc_xx[j] * fl1_fx - 0.5 * pc_xxx[j] * fl1_fx - 1.5 * pa_x[j] * fl1_fx * pc_zz[j] - 1.5 * pc_xzz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pa_z[j] - 1.5 * pc_xz[j] * fl1_fx * pb_z[j] - 3.0 * pa_xx[j] * pc_xzz[j] - 3.0 * pa_xz[j] * pc_xxz[j] - 3.0 * pa_x[j] * pc_xxz[j] * pb_z[j] - pc_xxx[j] * pa_z[j] * pb_z[j]);

                t_xxxz_z[j] += fl_s_0_0_4 * (0.5 * pc_xxx[j] * fl1_fx + 1.5 * pc_xzz[j] * fl1_fx + 3.0 * pa_x[j] * pc_xxzz[j] + pc_xxxz[j] * pa_z[j] + pc_xxxz[j] * pb_z[j]);

                t_xxxz_z[j] += -fl_s_0_0_5 * pc_xxxzz[j];

                t_xxyy_x[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl2_fx + pa_xyy[j] * fl1_fx + 0.25 * fl2_fx * pb_x[j] + 0.5 * pa_xx[j] * fl1_fx * pb_x[j] + 0.5 * fl1_fx * pa_yy[j] * pb_x[j] + pa_xxyy[j] * pb_x[j]);

                t_xxyy_x[j] += fl_s_0_0_1 * (-pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - 2.0 * pa_xy[j] * fl1_fx * pc_y[j] - pa_xyy[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pa_yy[j] - 0.5 * fl2_fx * pb_x[j] - 0.5 * pa_xx[j] * fl1_fx * pc_x[j] - 0.5 * pa_xx[j] * fl1_fx * pb_x[j] - pa_x[j] * pc_x[j] * fl1_fx * pb_x[j] - fl1_fx * pa_y[j] * pc_y[j] * pb_x[j] - 0.5 * fl1_fx * pa_yy[j] * pb_x[j] - pa_xxyy[j] * pc_x[j] - 2.0 * pa_xxy[j] * pc_y[j] * pb_x[j] - 2.0 * pa_xyy[j] * pc_x[j] * pb_x[j]);

                t_xxyy_x[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + pa_x[j] * fl1_fx * pc_yy[j] + 2.0 * pa_xy[j] * fl1_fx * pc_y[j] + 3.0 * pc_xy[j] * fl1_fx * pa_y[j] + 1.5 * pc_x[j] * fl1_fx * pa_yy[j] + 0.25 * fl2_fx * pb_x[j] + 0.5 * pa_xx[j] * fl1_fx * pc_x[j] + pa_x[j] * pc_xx[j] * fl1_fx + pa_x[j] * pc_x[j] * fl1_fx * pb_x[j] + 0.5 * pc_xx[j] * fl1_fx * pb_x[j] + 0.5 * fl1_fx * pc_yy[j] * pb_x[j] + fl1_fx * pa_y[j] * pc_y[j] * pb_x[j] + 2.0 * pa_xxy[j] * pc_xy[j] + pa_xx[j] * pc_yy[j] * pb_x[j] + 2.0 * pa_xyy[j] * pc_xx[j] + 4.0 * pa_xy[j] * pc_xy[j] * pb_x[j] + pc_xx[j] * pa_yy[j] * pb_x[j]);

                t_xxyy_x[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - pa_x[j] * fl1_fx * pc_yy[j] - 1.5 * pc_xyy[j] * fl1_fx - 3.0 * pc_xy[j] * fl1_fx * pa_y[j] - pa_x[j] * pc_xx[j] * fl1_fx - 0.5 * pc_xxx[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pb_x[j] - 0.5 * fl1_fx * pc_yy[j] * pb_x[j] - pa_xx[j] * pc_xyy[j] - 4.0 * pa_xy[j] * pc_xxy[j] - 2.0 * pa_x[j] * pc_xyy[j] * pb_x[j] - pc_xxx[j] * pa_yy[j] - 2.0 * pc_xxy[j] * pa_y[j] * pb_x[j]);

                t_xxyy_x[j] += fl_s_0_0_4 * (1.5 * pc_xyy[j] * fl1_fx + 0.5 * pc_xxx[j] * fl1_fx + 2.0 * pa_x[j] * pc_xxyy[j] + 2.0 * pc_xxxy[j] * pa_y[j] + pc_xxyy[j] * pb_x[j]);

                t_xxyy_x[j] += -fl_s_0_0_5 * pc_xxxyy[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP_10_15(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
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

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(55 * idx + 9);

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxy = pcDistances.data(55 * idx + 20);

            auto pc_xxxz = pcDistances.data(55 * idx + 21);

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxxyz = pcDistances.data(55 * idx + 38);

            auto pc_xxyyy = pcDistances.data(55 * idx + 40);

            auto pc_xxyyz = pcDistances.data(55 * idx + 41);

            auto pc_xxyzz = pcDistances.data(55 * idx + 42);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_xxyy_y = primBuffer.data(45 * idx + 10);

            auto t_xxyy_z = primBuffer.data(45 * idx + 11);

            auto t_xxyz_x = primBuffer.data(45 * idx + 12);

            auto t_xxyz_y = primBuffer.data(45 * idx + 13);

            auto t_xxyz_z = primBuffer.data(45 * idx + 14);

            // Batch of Integrals (10,15)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xxyz, pa_xxz, pa_xy, pa_xyy, pa_xyz, pa_xz, \
                                     pa_y, pa_yy, pa_yz, pa_z, pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xxx, pc_xxxy, pc_xxxyz, \
                                     pc_xxxz, pc_xxy, pc_xxyy, pc_xxyyy, pc_xxyyz, pc_xxyz, pc_xxyzz, pc_xxz, pc_xxzz, \
                                     pc_xy, pc_xyy, pc_xyyy, pc_xyyz, pc_xyz, pc_xyzz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yyy, \
                                     pc_yyz, pc_yz, pc_yzz, pc_z, pc_zz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, \
                                     s_0_0_5, t_xxyy_y, t_xxyy_z, t_xxyz_x, t_xxyz_y, t_xxyz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxyy_y[j] = fl_s_0_0_0 * (0.5 * fl2_fx * pa_y[j] + pa_xxy[j] * fl1_fx + 0.25 * fl2_fx * pb_y[j] + 0.5 * pa_xx[j] * fl1_fx * pb_y[j] + 0.5 * fl1_fx * pa_yy[j] * pb_y[j] + pa_xxyy[j] * pb_y[j]);

                t_xxyy_y[j] += fl_s_0_0_1 * (-fl2_fx * pa_y[j] - 0.75 * fl2_fx * pc_y[j] - pa_xxy[j] * fl1_fx - 1.5 * pa_xx[j] * pc_y[j] * fl1_fx - 2.0 * pa_xy[j] * pc_x[j] * fl1_fx - 0.5 * fl2_fx * pb_y[j] - 0.5 * pa_xx[j] * fl1_fx * pb_y[j] - pa_x[j] * pc_x[j] * fl1_fx * pb_y[j] - 0.5 * fl1_fx * pa_yy[j] * pc_y[j] - fl1_fx * pa_y[j] * pc_y[j] * pb_y[j] - 0.5 * fl1_fx * pa_yy[j] * pb_y[j] - pa_xxyy[j] * pc_y[j] - 2.0 * pa_xxy[j] * pc_y[j] * pb_y[j] - 2.0 * pa_xyy[j] * pc_x[j] * pb_y[j]);

                t_xxyy_y[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.5 * fl2_fx * pa_y[j] + 1.5 * pa_xx[j] * pc_y[j] * fl1_fx + 2.0 * pa_xy[j] * pc_x[j] * fl1_fx + 3.0 * pa_x[j] * pc_xy[j] * fl1_fx + pc_xx[j] * pa_y[j] * fl1_fx + 0.25 * fl2_fx * pb_y[j] + pa_x[j] * pc_x[j] * fl1_fx * pb_y[j] + 0.5 * pc_xx[j] * fl1_fx * pb_y[j] + fl1_fx * pa_y[j] * pc_yy[j] + 0.5 * fl1_fx * pc_yy[j] * pb_y[j] + 0.5 * fl1_fx * pa_yy[j] * pc_y[j] + fl1_fx * pa_y[j] * pc_y[j] * pb_y[j] + 2.0 * pa_xxy[j] * pc_yy[j] + pa_xx[j] * pc_yy[j] * pb_y[j] + 2.0 * pa_xyy[j] * pc_xy[j] + 4.0 * pa_xy[j] * pc_xy[j] * pb_y[j] + pc_xx[j] * pa_yy[j] * pb_y[j]);

                t_xxyy_y[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 3.0 * pa_x[j] * pc_xy[j] * fl1_fx - pc_xx[j] * pa_y[j] * fl1_fx - 1.5 * pc_xxy[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pb_y[j] - 0.5 * fl1_fx * pc_yyy[j] - fl1_fx * pa_y[j] * pc_yy[j] - 0.5 * fl1_fx * pc_yy[j] * pb_y[j] - pa_xx[j] * pc_yyy[j] - 4.0 * pa_xy[j] * pc_xyy[j] - 2.0 * pa_x[j] * pc_xyy[j] * pb_y[j] - pc_xxy[j] * pa_yy[j] - 2.0 * pc_xxy[j] * pa_y[j] * pb_y[j]);

                t_xxyy_y[j] += fl_s_0_0_4 * (1.5 * pc_xxy[j] * fl1_fx + 0.5 * fl1_fx * pc_yyy[j] + 2.0 * pa_x[j] * pc_xyyy[j] + 2.0 * pc_xxyy[j] * pa_y[j] + pc_xxyy[j] * pb_y[j]);

                t_xxyy_y[j] += -fl_s_0_0_5 * pc_xxyyy[j];

                t_xxyy_z[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_z[j] + 0.5 * pa_xx[j] * fl1_fx * pb_z[j] + 0.5 * fl1_fx * pa_yy[j] * pb_z[j] + pa_xxyy[j] * pb_z[j]);

                t_xxyy_z[j] += fl_s_0_0_1 * (-0.25 * fl2_fx * pc_z[j] - 0.5 * fl2_fx * pb_z[j] - 0.5 * pa_xx[j] * fl1_fx * pc_z[j] - 0.5 * pa_xx[j] * fl1_fx * pb_z[j] - pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pa_yy[j] * pc_z[j] - fl1_fx * pa_y[j] * pc_y[j] * pb_z[j] - 0.5 * fl1_fx * pa_yy[j] * pb_z[j] - pa_xxyy[j] * pc_z[j] - 2.0 * pa_xxy[j] * pc_y[j] * pb_z[j] - 2.0 * pa_xyy[j] * pc_x[j] * pb_z[j]);

                t_xxyy_z[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_z[j] + 0.25 * fl2_fx * pb_z[j] + 0.5 * pa_xx[j] * fl1_fx * pc_z[j] + pa_x[j] * pc_xz[j] * fl1_fx + pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] + 0.5 * pc_xx[j] * fl1_fx * pb_z[j] + fl1_fx * pa_y[j] * pc_yz[j] + 0.5 * fl1_fx * pc_yy[j] * pb_z[j] + 0.5 * fl1_fx * pa_yy[j] * pc_z[j] + fl1_fx * pa_y[j] * pc_y[j] * pb_z[j] + 2.0 * pa_xxy[j] * pc_yz[j] + pa_xx[j] * pc_yy[j] * pb_z[j] + 2.0 * pa_xyy[j] * pc_xz[j] + 4.0 * pa_xy[j] * pc_xy[j] * pb_z[j] + pc_xx[j] * pa_yy[j] * pb_z[j]);

                t_xxyy_z[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_z[j] - pa_x[j] * pc_xz[j] * fl1_fx - 0.5 * pc_xxz[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_yyz[j] - fl1_fx * pa_y[j] * pc_yz[j] - 0.5 * fl1_fx * pc_yy[j] * pb_z[j] - pa_xx[j] * pc_yyz[j] - 4.0 * pa_xy[j] * pc_xyz[j] - 2.0 * pa_x[j] * pc_xyy[j] * pb_z[j] - pc_xxz[j] * pa_yy[j] - 2.0 * pc_xxy[j] * pa_y[j] * pb_z[j]);

                t_xxyy_z[j] += fl_s_0_0_4 * (0.5 * pc_xxz[j] * fl1_fx + 0.5 * fl1_fx * pc_yyz[j] + 2.0 * pa_x[j] * pc_xyyz[j] + 2.0 * pc_xxyz[j] * pa_y[j] + pc_xxyy[j] * pb_z[j]);

                t_xxyy_z[j] += -fl_s_0_0_5 * pc_xxyyz[j];

                t_xxyz_x[j] = fl_s_0_0_0 * (pa_xyz[j] * fl1_fx + 0.5 * fl1_fx * pa_yz[j] * pb_x[j] + pa_xxyz[j] * pb_x[j]);

                t_xxyz_x[j] += fl_s_0_0_1 * (-pa_xy[j] * fl1_fx * pc_z[j] - pa_xz[j] * fl1_fx * pc_y[j] - pa_xyz[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pa_yz[j] - 0.5 * fl1_fx * pa_y[j] * pc_z[j] * pb_x[j] - 0.5 * fl1_fx * pc_y[j] * pa_z[j] * pb_x[j] - 0.5 * fl1_fx * pa_yz[j] * pb_x[j] - pa_xxyz[j] * pc_x[j] - pa_xxy[j] * pc_z[j] * pb_x[j] - pa_xxz[j] * pc_y[j] * pb_x[j] - 2.0 * pa_xyz[j] * pc_x[j] * pb_x[j]);

                t_xxyz_x[j] += fl_s_0_0_2 * (pa_x[j] * fl1_fx * pc_yz[j] + pa_xy[j] * fl1_fx * pc_z[j] + pa_xz[j] * fl1_fx * pc_y[j] + 1.5 * pc_xz[j] * fl1_fx * pa_y[j] + 1.5 * pc_xy[j] * fl1_fx * pa_z[j] + 1.5 * pc_x[j] * fl1_fx * pa_yz[j] + 0.5 * fl1_fx * pc_yz[j] * pb_x[j] + 0.5 * fl1_fx * pa_y[j] * pc_z[j] * pb_x[j] + 0.5 * fl1_fx * pc_y[j] * pa_z[j] * pb_x[j] + pa_xxy[j] * pc_xz[j] + pa_xxz[j] * pc_xy[j] + pa_xx[j] * pc_yz[j] * pb_x[j] + 2.0 * pa_xyz[j] * pc_xx[j] + 2.0 * pa_xy[j] * pc_xz[j] * pb_x[j] + 2.0 * pa_xz[j] * pc_xy[j] * pb_x[j] + pc_xx[j] * pa_yz[j] * pb_x[j]);

                t_xxyz_x[j] += fl_s_0_0_3 * (-pa_x[j] * fl1_fx * pc_yz[j] - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pa_y[j] - 1.5 * pc_xy[j] * fl1_fx * pa_z[j] - 0.5 * fl1_fx * pc_yz[j] * pb_x[j] - pa_xx[j] * pc_xyz[j] - 2.0 * pa_xy[j] * pc_xxz[j] - 2.0 * pa_xz[j] * pc_xxy[j] - 2.0 * pa_x[j] * pc_xyz[j] * pb_x[j] - pc_xxx[j] * pa_yz[j] - pc_xxz[j] * pa_y[j] * pb_x[j] - pc_xxy[j] * pa_z[j] * pb_x[j]);

                t_xxyz_x[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + 2.0 * pa_x[j] * pc_xxyz[j] + pc_xxxz[j] * pa_y[j] + pc_xxxy[j] * pa_z[j] + pc_xxyz[j] * pb_x[j]);

                t_xxyz_x[j] += -fl_s_0_0_5 * pc_xxxyz[j];

                t_xxyz_y[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pa_z[j] + 0.5 * pa_xxz[j] * fl1_fx + 0.5 * fl1_fx * pa_yz[j] * pb_y[j] + pa_xxyz[j] * pb_y[j]);

                t_xxyz_y[j] += fl_s_0_0_1 * (-0.25 * fl2_fx * pc_z[j] - 0.5 * fl2_fx * pa_z[j] - 0.5 * pa_xx[j] * fl1_fx * pc_z[j] - 0.5 * pa_xxz[j] * fl1_fx - pa_xz[j] * pc_x[j] * fl1_fx - 0.5 * fl1_fx * pa_yz[j] * pc_y[j] - 0.5 * fl1_fx * pa_y[j] * pc_z[j] * pb_y[j] - 0.5 * fl1_fx * pc_y[j] * pa_z[j] * pb_y[j] - 0.5 * fl1_fx * pa_yz[j] * pb_y[j] - pa_xxyz[j] * pc_y[j] - pa_xxy[j] * pc_z[j] * pb_y[j] - pa_xxz[j] * pc_y[j] * pb_y[j] - 2.0 * pa_xyz[j] * pc_x[j] * pb_y[j]);

                t_xxyz_y[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_z[j] + 0.25 * fl2_fx * pa_z[j] + 0.5 * pa_xx[j] * fl1_fx * pc_z[j] + pa_x[j] * pc_xz[j] * fl1_fx + pa_xz[j] * pc_x[j] * fl1_fx + 0.5 * pc_xx[j] * fl1_fx * pa_z[j] + 0.5 * fl1_fx * pa_y[j] * pc_yz[j] + 0.5 * fl1_fx * pc_yy[j] * pa_z[j] + 0.5 * fl1_fx * pc_yz[j] * pb_y[j] + 0.5 * fl1_fx * pa_yz[j] * pc_y[j] + 0.5 * fl1_fx * pa_y[j] * pc_z[j] * pb_y[j] + 0.5 * fl1_fx * pc_y[j] * pa_z[j] * pb_y[j] + pa_xxy[j] * pc_yz[j] + pa_xxz[j] * pc_yy[j] + pa_xx[j] * pc_yz[j] * pb_y[j] + 2.0 * pa_xyz[j] * pc_xy[j] + 2.0 * pa_xy[j] * pc_xz[j] * pb_y[j] + 2.0 * pa_xz[j] * pc_xy[j] * pb_y[j] + pc_xx[j] * pa_yz[j] * pb_y[j]);

                t_xxyz_y[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_z[j] - pa_x[j] * pc_xz[j] * fl1_fx - 0.5 * pc_xxz[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pa_z[j] - 0.5 * fl1_fx * pc_yyz[j] - 0.5 * fl1_fx * pa_y[j] * pc_yz[j] - 0.5 * fl1_fx * pc_yy[j] * pa_z[j] - 0.5 * fl1_fx * pc_yz[j] * pb_y[j] - pa_xx[j] * pc_yyz[j] - 2.0 * pa_xy[j] * pc_xyz[j] - 2.0 * pa_xz[j] * pc_xyy[j] - 2.0 * pa_x[j] * pc_xyz[j] * pb_y[j] - pc_xxy[j] * pa_yz[j] - pc_xxz[j] * pa_y[j] * pb_y[j] - pc_xxy[j] * pa_z[j] * pb_y[j]);

                t_xxyz_y[j] += fl_s_0_0_4 * (0.5 * pc_xxz[j] * fl1_fx + 0.5 * fl1_fx * pc_yyz[j] + 2.0 * pa_x[j] * pc_xyyz[j] + pc_xxyz[j] * pa_y[j] + pc_xxyy[j] * pa_z[j] + pc_xxyz[j] * pb_y[j]);

                t_xxyz_y[j] += -fl_s_0_0_5 * pc_xxyyz[j];

                t_xxyz_z[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pa_y[j] + 0.5 * pa_xxy[j] * fl1_fx + 0.5 * fl1_fx * pa_yz[j] * pb_z[j] + pa_xxyz[j] * pb_z[j]);

                t_xxyz_z[j] += fl_s_0_0_1 * (-0.5 * fl2_fx * pa_y[j] - 0.25 * fl2_fx * pc_y[j] - 0.5 * pa_xxy[j] * fl1_fx - 0.5 * pa_xx[j] * pc_y[j] * fl1_fx - pa_xy[j] * pc_x[j] * fl1_fx - 0.5 * fl1_fx * pa_yz[j] * pc_z[j] - 0.5 * fl1_fx * pa_y[j] * pc_z[j] * pb_z[j] - 0.5 * fl1_fx * pc_y[j] * pa_z[j] * pb_z[j] - 0.5 * fl1_fx * pa_yz[j] * pb_z[j] - pa_xxyz[j] * pc_z[j] - pa_xxy[j] * pc_z[j] * pb_z[j] - pa_xxz[j] * pc_y[j] * pb_z[j] - 2.0 * pa_xyz[j] * pc_x[j] * pb_z[j]);

                t_xxyz_z[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_y[j] + 0.25 * fl2_fx * pa_y[j] + 0.5 * pa_xx[j] * pc_y[j] * fl1_fx + pa_xy[j] * pc_x[j] * fl1_fx + pa_x[j] * pc_xy[j] * fl1_fx + 0.5 * pc_xx[j] * pa_y[j] * fl1_fx + 0.5 * fl1_fx * pa_y[j] * pc_zz[j] + 0.5 * fl1_fx * pc_yz[j] * pa_z[j] + 0.5 * fl1_fx * pc_yz[j] * pb_z[j] + 0.5 * fl1_fx * pa_yz[j] * pc_z[j] + 0.5 * fl1_fx * pa_y[j] * pc_z[j] * pb_z[j] + 0.5 * fl1_fx * pc_y[j] * pa_z[j] * pb_z[j] + pa_xxy[j] * pc_zz[j] + pa_xxz[j] * pc_yz[j] + pa_xx[j] * pc_yz[j] * pb_z[j] + 2.0 * pa_xyz[j] * pc_xz[j] + 2.0 * pa_xy[j] * pc_xz[j] * pb_z[j] + 2.0 * pa_xz[j] * pc_xy[j] * pb_z[j] + pc_xx[j] * pa_yz[j] * pb_z[j]);

                t_xxyz_z[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_y[j] - pa_x[j] * pc_xy[j] * fl1_fx - 0.5 * pc_xx[j] * pa_y[j] * fl1_fx - 0.5 * pc_xxy[j] * fl1_fx - 0.5 * fl1_fx * pc_yzz[j] - 0.5 * fl1_fx * pa_y[j] * pc_zz[j] - 0.5 * fl1_fx * pc_yz[j] * pa_z[j] - 0.5 * fl1_fx * pc_yz[j] * pb_z[j] - pa_xx[j] * pc_yzz[j] - 2.0 * pa_xy[j] * pc_xzz[j] - 2.0 * pa_xz[j] * pc_xyz[j] - 2.0 * pa_x[j] * pc_xyz[j] * pb_z[j] - pc_xxz[j] * pa_yz[j] - pc_xxz[j] * pa_y[j] * pb_z[j] - pc_xxy[j] * pa_z[j] * pb_z[j]);

                t_xxyz_z[j] += fl_s_0_0_4 * (0.5 * pc_xxy[j] * fl1_fx + 0.5 * fl1_fx * pc_yzz[j] + 2.0 * pa_x[j] * pc_xyzz[j] + pc_xxzz[j] * pa_y[j] + pc_xxyz[j] * pa_z[j] + pc_xxyz[j] * pb_z[j]);

                t_xxyz_z[j] += -fl_s_0_0_5 * pc_xxyzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP_15_20(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (15,20)

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

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxx = pcDistances.data(55 * idx + 9);

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxxz = pcDistances.data(55 * idx + 21);

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            auto pc_xzzz = pcDistances.data(55 * idx + 28);

            auto pc_yyyy = pcDistances.data(55 * idx + 29);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxxzz = pcDistances.data(55 * idx + 39);

            auto pc_xxyyy = pcDistances.data(55 * idx + 40);

            auto pc_xxyzz = pcDistances.data(55 * idx + 42);

            auto pc_xxzzz = pcDistances.data(55 * idx + 43);

            auto pc_xyyyy = pcDistances.data(55 * idx + 44);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_xxzz_x = primBuffer.data(45 * idx + 15);

            auto t_xxzz_y = primBuffer.data(45 * idx + 16);

            auto t_xxzz_z = primBuffer.data(45 * idx + 17);

            auto t_xyyy_x = primBuffer.data(45 * idx + 18);

            auto t_xyyy_y = primBuffer.data(45 * idx + 19);

            // Batch of Integrals (15,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xy, pa_xyy, pa_xyyy, pa_xz, pa_xzz, pa_y, \
                                     pa_yy, pa_yyy, pa_z, pa_zz, pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xxx, pc_xxxz, pc_xxxzz, \
                                     pc_xxy, pc_xxyy, pc_xxyyy, pc_xxyz, pc_xxyzz, pc_xxz, pc_xxzz, pc_xxzzz, pc_xy, \
                                     pc_xyy, pc_xyyy, pc_xyyyy, pc_xyz, pc_xyzz, pc_xz, pc_xzz, pc_xzzz, pc_y, pc_yy, \
                                     pc_yyy, pc_yyyy, pc_yz, pc_yzz, pc_z, pc_zz, pc_zzz, s_0_0_0, s_0_0_1, s_0_0_2, \
                                     s_0_0_3, s_0_0_4, s_0_0_5, t_xxzz_x, t_xxzz_y, t_xxzz_z, t_xyyy_x, t_xyyy_y: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxzz_x[j] = fl_s_0_0_0 * (0.5 * pa_x[j] * fl2_fx + pa_xzz[j] * fl1_fx + 0.25 * fl2_fx * pb_x[j] + 0.5 * pa_xx[j] * fl1_fx * pb_x[j] + 0.5 * fl1_fx * pa_zz[j] * pb_x[j] + pa_xxzz[j] * pb_x[j]);

                t_xxzz_x[j] += fl_s_0_0_1 * (-pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - 2.0 * pa_xz[j] * fl1_fx * pc_z[j] - pa_xzz[j] * fl1_fx - 1.5 * pc_x[j] * fl1_fx * pa_zz[j] - 0.5 * fl2_fx * pb_x[j] - 0.5 * pa_xx[j] * fl1_fx * pc_x[j] - 0.5 * pa_xx[j] * fl1_fx * pb_x[j] - pa_x[j] * pc_x[j] * fl1_fx * pb_x[j] - fl1_fx * pa_z[j] * pc_z[j] * pb_x[j] - 0.5 * fl1_fx * pa_zz[j] * pb_x[j] - pa_xxzz[j] * pc_x[j] - 2.0 * pa_xxz[j] * pc_z[j] * pb_x[j] - 2.0 * pa_xzz[j] * pc_x[j] * pb_x[j]);

                t_xxzz_x[j] += fl_s_0_0_2 * (0.5 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + pa_x[j] * fl1_fx * pc_zz[j] + 2.0 * pa_xz[j] * fl1_fx * pc_z[j] + 3.0 * pc_xz[j] * fl1_fx * pa_z[j] + 1.5 * pc_x[j] * fl1_fx * pa_zz[j] + 0.25 * fl2_fx * pb_x[j] + 0.5 * pa_xx[j] * fl1_fx * pc_x[j] + pa_x[j] * pc_xx[j] * fl1_fx + pa_x[j] * pc_x[j] * fl1_fx * pb_x[j] + 0.5 * pc_xx[j] * fl1_fx * pb_x[j] + 0.5 * fl1_fx * pc_zz[j] * pb_x[j] + fl1_fx * pa_z[j] * pc_z[j] * pb_x[j] + 2.0 * pa_xxz[j] * pc_xz[j] + pa_xx[j] * pc_zz[j] * pb_x[j] + 2.0 * pa_xzz[j] * pc_xx[j] + 4.0 * pa_xz[j] * pc_xz[j] * pb_x[j] + pc_xx[j] * pa_zz[j] * pb_x[j]);

                t_xxzz_x[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - pa_x[j] * fl1_fx * pc_zz[j] - 1.5 * pc_xzz[j] * fl1_fx - 3.0 * pc_xz[j] * fl1_fx * pa_z[j] - pa_x[j] * pc_xx[j] * fl1_fx - 0.5 * pc_xxx[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pb_x[j] - 0.5 * fl1_fx * pc_zz[j] * pb_x[j] - pa_xx[j] * pc_xzz[j] - 4.0 * pa_xz[j] * pc_xxz[j] - 2.0 * pa_x[j] * pc_xzz[j] * pb_x[j] - pc_xxx[j] * pa_zz[j] - 2.0 * pc_xxz[j] * pa_z[j] * pb_x[j]);

                t_xxzz_x[j] += fl_s_0_0_4 * (1.5 * pc_xzz[j] * fl1_fx + 0.5 * pc_xxx[j] * fl1_fx + 2.0 * pa_x[j] * pc_xxzz[j] + 2.0 * pc_xxxz[j] * pa_z[j] + pc_xxzz[j] * pb_x[j]);

                t_xxzz_x[j] += -fl_s_0_0_5 * pc_xxxzz[j];

                t_xxzz_y[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_y[j] + 0.5 * pa_xx[j] * fl1_fx * pb_y[j] + 0.5 * fl1_fx * pa_zz[j] * pb_y[j] + pa_xxzz[j] * pb_y[j]);

                t_xxzz_y[j] += fl_s_0_0_1 * (-0.25 * fl2_fx * pc_y[j] - 0.5 * fl2_fx * pb_y[j] - 0.5 * pa_xx[j] * fl1_fx * pc_y[j] - 0.5 * pa_xx[j] * fl1_fx * pb_y[j] - pa_x[j] * pc_x[j] * fl1_fx * pb_y[j] - 0.5 * fl1_fx * pa_zz[j] * pc_y[j] - fl1_fx * pa_z[j] * pc_z[j] * pb_y[j] - 0.5 * fl1_fx * pa_zz[j] * pb_y[j] - pa_xxzz[j] * pc_y[j] - 2.0 * pa_xxz[j] * pc_z[j] * pb_y[j] - 2.0 * pa_xzz[j] * pc_x[j] * pb_y[j]);

                t_xxzz_y[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_y[j] + 0.25 * fl2_fx * pb_y[j] + 0.5 * pa_xx[j] * fl1_fx * pc_y[j] + pa_x[j] * pc_xy[j] * fl1_fx + pa_x[j] * pc_x[j] * fl1_fx * pb_y[j] + 0.5 * pc_xx[j] * fl1_fx * pb_y[j] + fl1_fx * pa_z[j] * pc_yz[j] + 0.5 * fl1_fx * pc_zz[j] * pb_y[j] + 0.5 * fl1_fx * pa_zz[j] * pc_y[j] + fl1_fx * pa_z[j] * pc_z[j] * pb_y[j] + 2.0 * pa_xxz[j] * pc_yz[j] + pa_xx[j] * pc_zz[j] * pb_y[j] + 2.0 * pa_xzz[j] * pc_xy[j] + 4.0 * pa_xz[j] * pc_xz[j] * pb_y[j] + pc_xx[j] * pa_zz[j] * pb_y[j]);

                t_xxzz_y[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_y[j] - pa_x[j] * pc_xy[j] * fl1_fx - 0.5 * pc_xxy[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pb_y[j] - 0.5 * fl1_fx * pc_yzz[j] - fl1_fx * pa_z[j] * pc_yz[j] - 0.5 * fl1_fx * pc_zz[j] * pb_y[j] - pa_xx[j] * pc_yzz[j] - 4.0 * pa_xz[j] * pc_xyz[j] - 2.0 * pa_x[j] * pc_xzz[j] * pb_y[j] - pc_xxy[j] * pa_zz[j] - 2.0 * pc_xxz[j] * pa_z[j] * pb_y[j]);

                t_xxzz_y[j] += fl_s_0_0_4 * (0.5 * pc_xxy[j] * fl1_fx + 0.5 * fl1_fx * pc_yzz[j] + 2.0 * pa_x[j] * pc_xyzz[j] + 2.0 * pc_xxyz[j] * pa_z[j] + pc_xxzz[j] * pb_y[j]);

                t_xxzz_y[j] += -fl_s_0_0_5 * pc_xxyzz[j];

                t_xxzz_z[j] = fl_s_0_0_0 * (0.5 * fl2_fx * pa_z[j] + pa_xxz[j] * fl1_fx + 0.25 * fl2_fx * pb_z[j] + 0.5 * pa_xx[j] * fl1_fx * pb_z[j] + 0.5 * fl1_fx * pa_zz[j] * pb_z[j] + pa_xxzz[j] * pb_z[j]);

                t_xxzz_z[j] += fl_s_0_0_1 * (-fl2_fx * pa_z[j] - 0.75 * fl2_fx * pc_z[j] - pa_xxz[j] * fl1_fx - 1.5 * pa_xx[j] * pc_z[j] * fl1_fx - 2.0 * pa_xz[j] * pc_x[j] * fl1_fx - 0.5 * fl2_fx * pb_z[j] - 0.5 * pa_xx[j] * fl1_fx * pb_z[j] - pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pa_zz[j] * pc_z[j] - fl1_fx * pa_z[j] * pc_z[j] * pb_z[j] - 0.5 * fl1_fx * pa_zz[j] * pb_z[j] - pa_xxzz[j] * pc_z[j] - 2.0 * pa_xxz[j] * pc_z[j] * pb_z[j] - 2.0 * pa_xzz[j] * pc_x[j] * pb_z[j]);

                t_xxzz_z[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.5 * fl2_fx * pa_z[j] + 1.5 * pa_xx[j] * pc_z[j] * fl1_fx + 2.0 * pa_xz[j] * pc_x[j] * fl1_fx + 3.0 * pa_x[j] * pc_xz[j] * fl1_fx + pc_xx[j] * pa_z[j] * fl1_fx + 0.25 * fl2_fx * pb_z[j] + pa_x[j] * pc_x[j] * fl1_fx * pb_z[j] + 0.5 * pc_xx[j] * fl1_fx * pb_z[j] + fl1_fx * pa_z[j] * pc_zz[j] + 0.5 * fl1_fx * pc_zz[j] * pb_z[j] + 0.5 * fl1_fx * pa_zz[j] * pc_z[j] + fl1_fx * pa_z[j] * pc_z[j] * pb_z[j] + 2.0 * pa_xxz[j] * pc_zz[j] + pa_xx[j] * pc_zz[j] * pb_z[j] + 2.0 * pa_xzz[j] * pc_xz[j] + 4.0 * pa_xz[j] * pc_xz[j] * pb_z[j] + pc_xx[j] * pa_zz[j] * pb_z[j]);

                t_xxzz_z[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 3.0 * pa_x[j] * pc_xz[j] * fl1_fx - pc_xx[j] * pa_z[j] * fl1_fx - 1.5 * pc_xxz[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_zzz[j] - fl1_fx * pa_z[j] * pc_zz[j] - 0.5 * fl1_fx * pc_zz[j] * pb_z[j] - pa_xx[j] * pc_zzz[j] - 4.0 * pa_xz[j] * pc_xzz[j] - 2.0 * pa_x[j] * pc_xzz[j] * pb_z[j] - pc_xxz[j] * pa_zz[j] - 2.0 * pc_xxz[j] * pa_z[j] * pb_z[j]);

                t_xxzz_z[j] += fl_s_0_0_4 * (1.5 * pc_xxz[j] * fl1_fx + 0.5 * fl1_fx * pc_zzz[j] + 2.0 * pa_x[j] * pc_xzzz[j] + 2.0 * pc_xxzz[j] * pa_z[j] + pc_xxzz[j] * pb_z[j]);

                t_xxzz_z[j] += -fl_s_0_0_5 * pc_xxzzz[j];

                t_xyyy_x[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pa_y[j] + 0.5 * fl1_fx * pa_yyy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_x[j] + pa_xyyy[j] * pb_x[j]);

                t_xyyy_x[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pa_y[j] - 0.75 * fl2_fx * pc_y[j] - 1.5 * fl1_fx * pa_yy[j] * pc_y[j] - 0.5 * fl1_fx * pa_yyy[j] - 1.5 * pa_xy[j] * fl1_fx * pc_x[j] - 1.5 * pa_xy[j] * fl1_fx * pb_x[j] - 1.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_x[j] - 1.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_x[j] - pa_xyyy[j] * pc_x[j] - 3.0 * pa_xyy[j] * pc_y[j] * pb_x[j] - pc_x[j] * pa_yyy[j] * pb_x[j]);

                t_xyyy_x[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.75 * fl2_fx * pa_y[j] + 1.5 * fl1_fx * pa_y[j] * pc_yy[j] + 1.5 * fl1_fx * pa_yy[j] * pc_y[j] + 1.5 * pa_xy[j] * fl1_fx * pc_x[j] + 1.5 * pa_x[j] * pc_xy[j] * fl1_fx + 1.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_x[j] + 1.5 * pc_xx[j] * pa_y[j] * fl1_fx + 1.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_x[j] + 1.5 * pc_xy[j] * fl1_fx * pb_x[j] + 3.0 * pa_xyy[j] * pc_xy[j] + 3.0 * pa_xy[j] * pc_yy[j] * pb_x[j] + pc_xx[j] * pa_yyy[j] + 3.0 * pc_xy[j] * pa_yy[j] * pb_x[j]);

                t_xyyy_x[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 0.5 * fl1_fx * pc_yyy[j] - 1.5 * fl1_fx * pa_y[j] * pc_yy[j] - 1.5 * pa_x[j] * pc_xy[j] * fl1_fx - 1.5 * pc_xx[j] * pa_y[j] * fl1_fx - 1.5 * pc_xxy[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pb_x[j] - 3.0 * pa_xy[j] * pc_xyy[j] - pa_x[j] * pc_yyy[j] * pb_x[j] - 3.0 * pc_xxy[j] * pa_yy[j] - 3.0 * pc_xyy[j] * pa_y[j] * pb_x[j]);

                t_xyyy_x[j] += fl_s_0_0_4 * (0.5 * fl1_fx * pc_yyy[j] + 1.5 * pc_xxy[j] * fl1_fx + pa_x[j] * pc_xyyy[j] + 3.0 * pc_xxyy[j] * pa_y[j] + pc_xyyy[j] * pb_x[j]);

                t_xyyy_x[j] += -fl_s_0_0_5 * pc_xxyyy[j];

                t_xyyy_y[j] = fl_s_0_0_0 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pa_xyy[j] * fl1_fx + 1.5 * pa_xy[j] * fl1_fx * pb_y[j] + pa_xyyy[j] * pb_y[j]);

                t_xyyy_y[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - 1.5 * pa_xyy[j] * fl1_fx - 4.5 * pa_xy[j] * pc_y[j] * fl1_fx - 1.5 * pc_x[j] * pa_yy[j] * fl1_fx - 1.5 * pa_xy[j] * fl1_fx * pb_y[j] - 1.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_y[j] - 1.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_y[j] - pa_xyyy[j] * pc_y[j] - 3.0 * pa_xyy[j] * pc_y[j] * pb_y[j] - pc_x[j] * pa_yyy[j] * pb_y[j]);

                t_xyyy_y[j] += fl_s_0_0_2 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + 4.5 * pa_xy[j] * pc_y[j] * fl1_fx + 3.0 * pa_x[j] * pc_yy[j] * fl1_fx + 1.5 * pc_x[j] * pa_yy[j] * fl1_fx + 4.5 * pc_xy[j] * pa_y[j] * fl1_fx + 1.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_y[j] + 1.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_y[j] + 1.5 * pc_xy[j] * fl1_fx * pb_y[j] + 3.0 * pa_xyy[j] * pc_yy[j] + 3.0 * pa_xy[j] * pc_yy[j] * pb_y[j] + pc_xy[j] * pa_yyy[j] + 3.0 * pc_xy[j] * pa_yy[j] * pb_y[j]);

                t_xyyy_y[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - 3.0 * pa_x[j] * pc_yy[j] * fl1_fx - 4.5 * pc_xy[j] * pa_y[j] * fl1_fx - 3.0 * pc_xyy[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pb_y[j] - 3.0 * pa_xy[j] * pc_yyy[j] - pa_x[j] * pc_yyy[j] * pb_y[j] - 3.0 * pc_xyy[j] * pa_yy[j] - 3.0 * pc_xyy[j] * pa_y[j] * pb_y[j]);

                t_xyyy_y[j] += fl_s_0_0_4 * (3.0 * pc_xyy[j] * fl1_fx + pa_x[j] * pc_yyyy[j] + 3.0 * pc_xyyy[j] * pa_y[j] + pc_xyyy[j] * pb_y[j]);

                t_xyyy_y[j] += -fl_s_0_0_5 * pc_xyyyy[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP_20_25(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (20,25)

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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxy = pcDistances.data(55 * idx + 10);

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxyy = pcDistances.data(55 * idx + 22);

            auto pc_xxyz = pcDistances.data(55 * idx + 23);

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            auto pc_yyyz = pcDistances.data(55 * idx + 30);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxyyz = pcDistances.data(55 * idx + 41);

            auto pc_xxyzz = pcDistances.data(55 * idx + 42);

            auto pc_xyyyz = pcDistances.data(55 * idx + 45);

            auto pc_xyyzz = pcDistances.data(55 * idx + 46);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_xyyy_z = primBuffer.data(45 * idx + 20);

            auto t_xyyz_x = primBuffer.data(45 * idx + 21);

            auto t_xyyz_y = primBuffer.data(45 * idx + 22);

            auto t_xyyz_z = primBuffer.data(45 * idx + 23);

            auto t_xyzz_x = primBuffer.data(45 * idx + 24);

            // Batch of Integrals (20,25)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_xyyz, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, \
                                     pa_yy, pa_yyy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xxy, \
                                     pc_xxyy, pc_xxyyz, pc_xxyz, pc_xxyzz, pc_xxz, pc_xxzz, pc_xy, pc_xyy, pc_xyyy, \
                                     pc_xyyyz, pc_xyyz, pc_xyyzz, pc_xyz, pc_xyzz, pc_xz, pc_xzz, pc_y, pc_yy, pc_yyy, \
                                     pc_yyyz, pc_yyz, pc_yyzz, pc_yz, pc_yzz, pc_z, pc_zz, s_0_0_0, s_0_0_1, s_0_0_2, \
                                     s_0_0_3, s_0_0_4, s_0_0_5, t_xyyy_z, t_xyyz_x, t_xyyz_y, t_xyyz_z, t_xyzz_x: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xyyy_z[j] = fl_s_0_0_0 * (1.5 * pa_xy[j] * fl1_fx * pb_z[j] + pa_xyyy[j] * pb_z[j]);

                t_xyyy_z[j] += fl_s_0_0_1 * (-1.5 * pa_xy[j] * fl1_fx * pc_z[j] - 1.5 * pa_xy[j] * fl1_fx * pb_z[j] - 1.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_z[j] - 1.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_z[j] - pa_xyyy[j] * pc_z[j] - 3.0 * pa_xyy[j] * pc_y[j] * pb_z[j] - pc_x[j] * pa_yyy[j] * pb_z[j]);

                t_xyyy_z[j] += fl_s_0_0_2 * (1.5 * pa_xy[j] * fl1_fx * pc_z[j] + 1.5 * pa_x[j] * pc_yz[j] * fl1_fx + 1.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_z[j] + 1.5 * pc_xz[j] * pa_y[j] * fl1_fx + 1.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_z[j] + 1.5 * pc_xy[j] * fl1_fx * pb_z[j] + 3.0 * pa_xyy[j] * pc_yz[j] + 3.0 * pa_xy[j] * pc_yy[j] * pb_z[j] + pc_xz[j] * pa_yyy[j] + 3.0 * pc_xy[j] * pa_yy[j] * pb_z[j]);

                t_xyyy_z[j] += fl_s_0_0_3 * (-1.5 * pa_x[j] * pc_yz[j] * fl1_fx - 1.5 * pc_xz[j] * pa_y[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pb_z[j] - 3.0 * pa_xy[j] * pc_yyz[j] - pa_x[j] * pc_yyy[j] * pb_z[j] - 3.0 * pc_xyz[j] * pa_yy[j] - 3.0 * pc_xyy[j] * pa_y[j] * pb_z[j]);

                t_xyyy_z[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_x[j] * pc_yyyz[j] + 3.0 * pc_xyyz[j] * pa_y[j] + pc_xyyy[j] * pb_z[j]);

                t_xyyy_z[j] += -fl_s_0_0_5 * pc_xyyyz[j];

                t_xyyz_x[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pa_z[j] + 0.5 * fl1_fx * pa_yyz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_x[j] + pa_xyyz[j] * pb_x[j]);

                t_xyyz_x[j] += fl_s_0_0_1 * (-0.25 * fl2_fx * pc_z[j] - 0.5 * fl2_fx * pa_z[j] - 0.5 * fl1_fx * pa_yy[j] * pc_z[j] - fl1_fx * pa_yz[j] * pc_y[j] - 0.5 * fl1_fx * pa_yyz[j] - 0.5 * pa_xz[j] * fl1_fx * pc_x[j] - 0.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_x[j] - 0.5 * pa_xz[j] * fl1_fx * pb_x[j] - 0.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_x[j] - pa_xyyz[j] * pc_x[j] - pa_xyy[j] * pc_z[j] * pb_x[j] - 2.0 * pa_xyz[j] * pc_y[j] * pb_x[j] - pc_x[j] * pa_yyz[j] * pb_x[j]);

                t_xyyz_x[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_z[j] + 0.25 * fl2_fx * pa_z[j] + fl1_fx * pa_y[j] * pc_yz[j] + 0.5 * fl1_fx * pc_yy[j] * pa_z[j] + 0.5 * fl1_fx * pa_yy[j] * pc_z[j] + fl1_fx * pa_yz[j] * pc_y[j] + 0.5 * pa_x[j] * fl1_fx * pc_xz[j] + 0.5 * pa_xz[j] * fl1_fx * pc_x[j] + 0.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_x[j] + 0.5 * pc_xx[j] * fl1_fx * pa_z[j] + 0.5 * pc_xz[j] * fl1_fx * pb_x[j] + 0.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_x[j] + pa_xyy[j] * pc_xz[j] + 2.0 * pa_xyz[j] * pc_xy[j] + 2.0 * pa_xy[j] * pc_yz[j] * pb_x[j] + pa_xz[j] * pc_yy[j] * pb_x[j] + pc_xx[j] * pa_yyz[j] + pc_xz[j] * pa_yy[j] * pb_x[j] + 2.0 * pc_xy[j] * pa_yz[j] * pb_x[j]);

                t_xyyz_x[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_z[j] - 0.5 * fl1_fx * pc_yyz[j] - fl1_fx * pa_y[j] * pc_yz[j] - 0.5 * fl1_fx * pc_yy[j] * pa_z[j] - 0.5 * pa_x[j] * fl1_fx * pc_xz[j] - 0.5 * pc_xxz[j] * fl1_fx - 0.5 * pc_xx[j] * fl1_fx * pa_z[j] - 0.5 * pc_xz[j] * fl1_fx * pb_x[j] - 2.0 * pa_xy[j] * pc_xyz[j] - pa_xz[j] * pc_xyy[j] - pa_x[j] * pc_yyz[j] * pb_x[j] - pc_xxz[j] * pa_yy[j] - 2.0 * pc_xxy[j] * pa_yz[j] - 2.0 * pc_xyz[j] * pa_y[j] * pb_x[j] - pc_xyy[j] * pa_z[j] * pb_x[j]);

                t_xyyz_x[j] += fl_s_0_0_4 * (0.5 * fl1_fx * pc_yyz[j] + 0.5 * pc_xxz[j] * fl1_fx + pa_x[j] * pc_xyyz[j] + 2.0 * pc_xxyz[j] * pa_y[j] + pc_xxyy[j] * pa_z[j] + pc_xyyz[j] * pb_x[j]);

                t_xyyz_x[j] += -fl_s_0_0_5 * pc_xxyyz[j];

                t_xyyz_y[j] = fl_s_0_0_0 * (pa_xyz[j] * fl1_fx + 0.5 * pa_xz[j] * fl1_fx * pb_y[j] + pa_xyyz[j] * pb_y[j]);

                t_xyyz_y[j] += fl_s_0_0_1 * (-pa_xy[j] * fl1_fx * pc_z[j] - pa_xyz[j] * fl1_fx - 1.5 * pa_xz[j] * pc_y[j] * fl1_fx - pc_x[j] * pa_yz[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_y[j] - 0.5 * pa_xz[j] * fl1_fx * pb_y[j] - 0.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_y[j] - pa_xyyz[j] * pc_y[j] - pa_xyy[j] * pc_z[j] * pb_y[j] - 2.0 * pa_xyz[j] * pc_y[j] * pb_y[j] - pc_x[j] * pa_yyz[j] * pb_y[j]);

                t_xyyz_y[j] += fl_s_0_0_2 * (pa_xy[j] * fl1_fx * pc_z[j] + 1.5 * pa_x[j] * pc_yz[j] * fl1_fx + 1.5 * pa_xz[j] * pc_y[j] * fl1_fx + pc_xz[j] * pa_y[j] * fl1_fx + pc_x[j] * pa_yz[j] * fl1_fx + 1.5 * pc_xy[j] * fl1_fx * pa_z[j] + 0.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_y[j] + 0.5 * pc_xz[j] * fl1_fx * pb_y[j] + 0.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_y[j] + pa_xyy[j] * pc_yz[j] + 2.0 * pa_xyz[j] * pc_yy[j] + 2.0 * pa_xy[j] * pc_yz[j] * pb_y[j] + pa_xz[j] * pc_yy[j] * pb_y[j] + pc_xy[j] * pa_yyz[j] + pc_xz[j] * pa_yy[j] * pb_y[j] + 2.0 * pc_xy[j] * pa_yz[j] * pb_y[j]);

                t_xyyz_y[j] += fl_s_0_0_3 * (-1.5 * pa_x[j] * pc_yz[j] * fl1_fx - pc_xz[j] * pa_y[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pa_z[j] - 0.5 * pc_xz[j] * fl1_fx * pb_y[j] - 2.0 * pa_xy[j] * pc_yyz[j] - pa_xz[j] * pc_yyy[j] - pa_x[j] * pc_yyz[j] * pb_y[j] - pc_xyz[j] * pa_yy[j] - 2.0 * pc_xyy[j] * pa_yz[j] - 2.0 * pc_xyz[j] * pa_y[j] * pb_y[j] - pc_xyy[j] * pa_z[j] * pb_y[j]);

                t_xyyz_y[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_x[j] * pc_yyyz[j] + 2.0 * pc_xyyz[j] * pa_y[j] + pc_xyyy[j] * pa_z[j] + pc_xyyz[j] * pb_y[j]);

                t_xyyz_y[j] += -fl_s_0_0_5 * pc_xyyyz[j];

                t_xyyz_z[j] = fl_s_0_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa_xyy[j] * fl1_fx + 0.5 * pa_xz[j] * fl1_fx * pb_z[j] + pa_xyyz[j] * pb_z[j]);

                t_xyyz_z[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl2_fx - 0.25 * pc_x[j] * fl2_fx - 0.5 * pa_xyy[j] * fl1_fx - pa_xy[j] * pc_y[j] * fl1_fx - 0.5 * pc_x[j] * pa_yy[j] * fl1_fx - 0.5 * pa_xz[j] * fl1_fx * pc_z[j] - 0.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_z[j] - 0.5 * pa_xz[j] * fl1_fx * pb_z[j] - 0.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_z[j] - pa_xyyz[j] * pc_z[j] - pa_xyy[j] * pc_z[j] * pb_z[j] - 2.0 * pa_xyz[j] * pc_y[j] * pb_z[j] - pc_x[j] * pa_yyz[j] * pb_z[j]);

                t_xyyz_z[j] += fl_s_0_0_2 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pc_x[j] * fl2_fx + pa_xy[j] * pc_y[j] * fl1_fx + 0.5 * pa_x[j] * pc_yy[j] * fl1_fx + 0.5 * pc_x[j] * pa_yy[j] * fl1_fx + pc_xy[j] * pa_y[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pc_zz[j] + 0.5 * pa_xz[j] * fl1_fx * pc_z[j] + 0.5 * pa_x[j] * fl1_fx * pc_z[j] * pb_z[j] + 0.5 * pc_xz[j] * fl1_fx * pa_z[j] + 0.5 * pc_xz[j] * fl1_fx * pb_z[j] + 0.5 * pc_x[j] * fl1_fx * pa_z[j] * pb_z[j] + pa_xyy[j] * pc_zz[j] + 2.0 * pa_xyz[j] * pc_yz[j] + 2.0 * pa_xy[j] * pc_yz[j] * pb_z[j] + pa_xz[j] * pc_yy[j] * pb_z[j] + pc_xz[j] * pa_yyz[j] + pc_xz[j] * pa_yy[j] * pb_z[j] + 2.0 * pc_xy[j] * pa_yz[j] * pb_z[j]);

                t_xyyz_z[j] += fl_s_0_0_3 * (-0.25 * pc_x[j] * fl2_fx - 0.5 * pa_x[j] * pc_yy[j] * fl1_fx - pc_xy[j] * pa_y[j] * fl1_fx - 0.5 * pc_xyy[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fx * pc_zz[j] - 0.5 * pc_xzz[j] * fl1_fx - 0.5 * pc_xz[j] * fl1_fx * pa_z[j] - 0.5 * pc_xz[j] * fl1_fx * pb_z[j] - 2.0 * pa_xy[j] * pc_yzz[j] - pa_xz[j] * pc_yyz[j] - pa_x[j] * pc_yyz[j] * pb_z[j] - pc_xzz[j] * pa_yy[j] - 2.0 * pc_xyz[j] * pa_yz[j] - 2.0 * pc_xyz[j] * pa_y[j] * pb_z[j] - pc_xyy[j] * pa_z[j] * pb_z[j]);

                t_xyyz_z[j] += fl_s_0_0_4 * (0.5 * pc_xyy[j] * fl1_fx + 0.5 * pc_xzz[j] * fl1_fx + pa_x[j] * pc_yyzz[j] + 2.0 * pc_xyzz[j] * pa_y[j] + pc_xyyz[j] * pa_z[j] + pc_xyyz[j] * pb_z[j]);

                t_xyyz_z[j] += -fl_s_0_0_5 * pc_xyyzz[j];

                t_xyzz_x[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pa_y[j] + 0.5 * fl1_fx * pa_yzz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_x[j] + pa_xyzz[j] * pb_x[j]);

                t_xyzz_x[j] += fl_s_0_0_1 * (-0.5 * fl2_fx * pa_y[j] - 0.25 * fl2_fx * pc_y[j] - fl1_fx * pa_yz[j] * pc_z[j] - 0.5 * fl1_fx * pc_y[j] * pa_zz[j] - 0.5 * fl1_fx * pa_yzz[j] - 0.5 * pa_xy[j] * fl1_fx * pc_x[j] - 0.5 * pa_xy[j] * fl1_fx * pb_x[j] - 0.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_x[j] - 0.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_x[j] - pa_xyzz[j] * pc_x[j] - 2.0 * pa_xyz[j] * pc_z[j] * pb_x[j] - pa_xzz[j] * pc_y[j] * pb_x[j] - pc_x[j] * pa_yzz[j] * pb_x[j]);

                t_xyzz_x[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_y[j] + 0.25 * fl2_fx * pa_y[j] + 0.5 * fl1_fx * pa_y[j] * pc_zz[j] + fl1_fx * pc_yz[j] * pa_z[j] + fl1_fx * pa_yz[j] * pc_z[j] + 0.5 * fl1_fx * pc_y[j] * pa_zz[j] + 0.5 * pa_xy[j] * fl1_fx * pc_x[j] + 0.5 * pa_x[j] * pc_xy[j] * fl1_fx + 0.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_x[j] + 0.5 * pc_xx[j] * pa_y[j] * fl1_fx + 0.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_x[j] + 0.5 * pc_xy[j] * fl1_fx * pb_x[j] + 2.0 * pa_xyz[j] * pc_xz[j] + pa_xy[j] * pc_zz[j] * pb_x[j] + pa_xzz[j] * pc_xy[j] + 2.0 * pa_xz[j] * pc_yz[j] * pb_x[j] + pc_xx[j] * pa_yzz[j] + 2.0 * pc_xz[j] * pa_yz[j] * pb_x[j] + pc_xy[j] * pa_zz[j] * pb_x[j]);

                t_xyzz_x[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_y[j] - 0.5 * fl1_fx * pc_yzz[j] - 0.5 * fl1_fx * pa_y[j] * pc_zz[j] - fl1_fx * pc_yz[j] * pa_z[j] - 0.5 * pa_x[j] * pc_xy[j] * fl1_fx - 0.5 * pc_xx[j] * pa_y[j] * fl1_fx - 0.5 * pc_xxy[j] * fl1_fx - 0.5 * pc_xy[j] * fl1_fx * pb_x[j] - pa_xy[j] * pc_xzz[j] - 2.0 * pa_xz[j] * pc_xyz[j] - pa_x[j] * pc_yzz[j] * pb_x[j] - 2.0 * pc_xxz[j] * pa_yz[j] - pc_xzz[j] * pa_y[j] * pb_x[j] - pc_xxy[j] * pa_zz[j] - 2.0 * pc_xyz[j] * pa_z[j] * pb_x[j]);

                t_xyzz_x[j] += fl_s_0_0_4 * (0.5 * fl1_fx * pc_yzz[j] + 0.5 * pc_xxy[j] * fl1_fx + pa_x[j] * pc_xyzz[j] + pc_xxzz[j] * pa_y[j] + 2.0 * pc_xxyz[j] * pa_z[j] + pc_xyzz[j] * pb_x[j]);

                t_xyzz_x[j] += -fl_s_0_0_5 * pc_xxyzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP_25_30(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (25,30)

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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xx = pcDistances.data(55 * idx + 3);

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xxz = pcDistances.data(55 * idx + 11);

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xxzz = pcDistances.data(55 * idx + 24);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            auto pc_xzzz = pcDistances.data(55 * idx + 28);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            auto pc_yzzz = pcDistances.data(55 * idx + 32);

            auto pc_zzzz = pcDistances.data(55 * idx + 33);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xxzzz = pcDistances.data(55 * idx + 43);

            auto pc_xyyzz = pcDistances.data(55 * idx + 46);

            auto pc_xyzzz = pcDistances.data(55 * idx + 47);

            auto pc_xzzzz = pcDistances.data(55 * idx + 48);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_xyzz_y = primBuffer.data(45 * idx + 25);

            auto t_xyzz_z = primBuffer.data(45 * idx + 26);

            auto t_xzzz_x = primBuffer.data(45 * idx + 27);

            auto t_xzzz_y = primBuffer.data(45 * idx + 28);

            auto t_xzzz_z = primBuffer.data(45 * idx + 29);

            // Batch of Integrals (25,30)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_xzzz, pa_y, pa_yz, pa_yzz, \
                                     pa_z, pa_zz, pa_zzz, pb_x, pb_y, pb_z, pc_x, pc_xx, pc_xxz, pc_xxzz, pc_xxzzz, pc_xy, \
                                     pc_xyy, pc_xyyz, pc_xyyzz, pc_xyz, pc_xyzz, pc_xyzzz, pc_xz, pc_xzz, pc_xzzz, \
                                     pc_xzzzz, pc_y, pc_yy, pc_yyz, pc_yyzz, pc_yz, pc_yzz, pc_yzzz, pc_z, pc_zz, pc_zzz, \
                                     pc_zzzz, s_0_0_0, s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_xyzz_y, t_xyzz_z, \
                                     t_xzzz_x, t_xzzz_y, t_xzzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xyzz_y[j] = fl_s_0_0_0 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pa_xzz[j] * fl1_fx + 0.5 * pa_xy[j] * fl1_fx * pb_y[j] + pa_xyzz[j] * pb_y[j]);

                t_xyzz_y[j] += fl_s_0_0_1 * (-0.5 * pa_x[j] * fl2_fx - 0.25 * pc_x[j] * fl2_fx - pa_xz[j] * fl1_fx * pc_z[j] - 0.5 * pa_xzz[j] * fl1_fx - 0.5 * pc_x[j] * fl1_fx * pa_zz[j] - 0.5 * pa_xy[j] * fl1_fx * pc_y[j] - 0.5 * pa_xy[j] * fl1_fx * pb_y[j] - 0.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_y[j] - 0.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_y[j] - pa_xyzz[j] * pc_y[j] - 2.0 * pa_xyz[j] * pc_z[j] * pb_y[j] - pa_xzz[j] * pc_y[j] * pb_y[j] - pc_x[j] * pa_yzz[j] * pb_y[j]);

                t_xyzz_y[j] += fl_s_0_0_2 * (0.25 * pa_x[j] * fl2_fx + 0.5 * pc_x[j] * fl2_fx + 0.5 * pa_x[j] * fl1_fx * pc_zz[j] + pa_xz[j] * fl1_fx * pc_z[j] + pc_xz[j] * fl1_fx * pa_z[j] + 0.5 * pc_x[j] * fl1_fx * pa_zz[j] + 0.5 * pa_xy[j] * fl1_fx * pc_y[j] + 0.5 * pa_x[j] * pc_yy[j] * fl1_fx + 0.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_y[j] + 0.5 * pc_xy[j] * pa_y[j] * fl1_fx + 0.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_y[j] + 0.5 * pc_xy[j] * fl1_fx * pb_y[j] + 2.0 * pa_xyz[j] * pc_yz[j] + pa_xy[j] * pc_zz[j] * pb_y[j] + pa_xzz[j] * pc_yy[j] + 2.0 * pa_xz[j] * pc_yz[j] * pb_y[j] + pc_xy[j] * pa_yzz[j] + 2.0 * pc_xz[j] * pa_yz[j] * pb_y[j] + pc_xy[j] * pa_zz[j] * pb_y[j]);

                t_xyzz_y[j] += fl_s_0_0_3 * (-0.25 * pc_x[j] * fl2_fx - 0.5 * pa_x[j] * fl1_fx * pc_zz[j] - 0.5 * pc_xzz[j] * fl1_fx - pc_xz[j] * fl1_fx * pa_z[j] - 0.5 * pa_x[j] * pc_yy[j] * fl1_fx - 0.5 * pc_xy[j] * pa_y[j] * fl1_fx - 0.5 * pc_xyy[j] * fl1_fx - 0.5 * pc_xy[j] * fl1_fx * pb_y[j] - pa_xy[j] * pc_yzz[j] - 2.0 * pa_xz[j] * pc_yyz[j] - pa_x[j] * pc_yzz[j] * pb_y[j] - 2.0 * pc_xyz[j] * pa_yz[j] - pc_xzz[j] * pa_y[j] * pb_y[j] - pc_xyy[j] * pa_zz[j] - 2.0 * pc_xyz[j] * pa_z[j] * pb_y[j]);

                t_xyzz_y[j] += fl_s_0_0_4 * (0.5 * pc_xzz[j] * fl1_fx + 0.5 * pc_xyy[j] * fl1_fx + pa_x[j] * pc_yyzz[j] + pc_xyzz[j] * pa_y[j] + 2.0 * pc_xyyz[j] * pa_z[j] + pc_xyzz[j] * pb_y[j]);

                t_xyzz_y[j] += -fl_s_0_0_5 * pc_xyyzz[j];

                t_xyzz_z[j] = fl_s_0_0_0 * (pa_xyz[j] * fl1_fx + 0.5 * pa_xy[j] * fl1_fx * pb_z[j] + pa_xyzz[j] * pb_z[j]);

                t_xyzz_z[j] += fl_s_0_0_1 * (-pa_xyz[j] * fl1_fx - 1.5 * pa_xy[j] * pc_z[j] * fl1_fx - pa_xz[j] * pc_y[j] * fl1_fx - pc_x[j] * pa_yz[j] * fl1_fx - 0.5 * pa_xy[j] * fl1_fx * pb_z[j] - 0.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_z[j] - 0.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_z[j] - pa_xyzz[j] * pc_z[j] - 2.0 * pa_xyz[j] * pc_z[j] * pb_z[j] - pa_xzz[j] * pc_y[j] * pb_z[j] - pc_x[j] * pa_yzz[j] * pb_z[j]);

                t_xyzz_z[j] += fl_s_0_0_2 * (1.5 * pa_xy[j] * pc_z[j] * fl1_fx + pa_xz[j] * pc_y[j] * fl1_fx + 1.5 * pa_x[j] * pc_yz[j] * fl1_fx + pc_x[j] * pa_yz[j] * fl1_fx + 1.5 * pc_xz[j] * pa_y[j] * fl1_fx + pc_xy[j] * pa_z[j] * fl1_fx + 0.5 * pa_x[j] * pc_y[j] * fl1_fx * pb_z[j] + 0.5 * pc_x[j] * pa_y[j] * fl1_fx * pb_z[j] + 0.5 * pc_xy[j] * fl1_fx * pb_z[j] + 2.0 * pa_xyz[j] * pc_zz[j] + pa_xy[j] * pc_zz[j] * pb_z[j] + pa_xzz[j] * pc_yz[j] + 2.0 * pa_xz[j] * pc_yz[j] * pb_z[j] + pc_xz[j] * pa_yzz[j] + 2.0 * pc_xz[j] * pa_yz[j] * pb_z[j] + pc_xy[j] * pa_zz[j] * pb_z[j]);

                t_xyzz_z[j] += fl_s_0_0_3 * (-1.5 * pa_x[j] * pc_yz[j] * fl1_fx - 1.5 * pc_xz[j] * pa_y[j] * fl1_fx - pc_xy[j] * pa_z[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 0.5 * pc_xy[j] * fl1_fx * pb_z[j] - pa_xy[j] * pc_zzz[j] - 2.0 * pa_xz[j] * pc_yzz[j] - pa_x[j] * pc_yzz[j] * pb_z[j] - 2.0 * pc_xzz[j] * pa_yz[j] - pc_xzz[j] * pa_y[j] * pb_z[j] - pc_xyz[j] * pa_zz[j] - 2.0 * pc_xyz[j] * pa_z[j] * pb_z[j]);

                t_xyzz_z[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_x[j] * pc_yzzz[j] + pc_xzzz[j] * pa_y[j] + 2.0 * pc_xyzz[j] * pa_z[j] + pc_xyzz[j] * pb_z[j]);

                t_xyzz_z[j] += -fl_s_0_0_5 * pc_xyzzz[j];

                t_xzzz_x[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pa_z[j] + 0.5 * fl1_fx * pa_zzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_x[j] + pa_xzzz[j] * pb_x[j]);

                t_xzzz_x[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pa_z[j] - 0.75 * fl2_fx * pc_z[j] - 1.5 * fl1_fx * pa_zz[j] * pc_z[j] - 0.5 * fl1_fx * pa_zzz[j] - 1.5 * pa_xz[j] * fl1_fx * pc_x[j] - 1.5 * pa_xz[j] * fl1_fx * pb_x[j] - 1.5 * pa_x[j] * pc_z[j] * fl1_fx * pb_x[j] - 1.5 * pc_x[j] * pa_z[j] * fl1_fx * pb_x[j] - pa_xzzz[j] * pc_x[j] - 3.0 * pa_xzz[j] * pc_z[j] * pb_x[j] - pc_x[j] * pa_zzz[j] * pb_x[j]);

                t_xzzz_x[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pa_z[j] + 1.5 * fl1_fx * pa_z[j] * pc_zz[j] + 1.5 * fl1_fx * pa_zz[j] * pc_z[j] + 1.5 * pa_xz[j] * fl1_fx * pc_x[j] + 1.5 * pa_x[j] * pc_xz[j] * fl1_fx + 1.5 * pa_x[j] * pc_z[j] * fl1_fx * pb_x[j] + 1.5 * pc_xx[j] * pa_z[j] * fl1_fx + 1.5 * pc_x[j] * pa_z[j] * fl1_fx * pb_x[j] + 1.5 * pc_xz[j] * fl1_fx * pb_x[j] + 3.0 * pa_xzz[j] * pc_xz[j] + 3.0 * pa_xz[j] * pc_zz[j] * pb_x[j] + pc_xx[j] * pa_zzz[j] + 3.0 * pc_xz[j] * pa_zz[j] * pb_x[j]);

                t_xzzz_x[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 0.5 * fl1_fx * pc_zzz[j] - 1.5 * fl1_fx * pa_z[j] * pc_zz[j] - 1.5 * pa_x[j] * pc_xz[j] * fl1_fx - 1.5 * pc_xx[j] * pa_z[j] * fl1_fx - 1.5 * pc_xxz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pb_x[j] - 3.0 * pa_xz[j] * pc_xzz[j] - pa_x[j] * pc_zzz[j] * pb_x[j] - 3.0 * pc_xxz[j] * pa_zz[j] - 3.0 * pc_xzz[j] * pa_z[j] * pb_x[j]);

                t_xzzz_x[j] += fl_s_0_0_4 * (0.5 * fl1_fx * pc_zzz[j] + 1.5 * pc_xxz[j] * fl1_fx + pa_x[j] * pc_xzzz[j] + 3.0 * pc_xxzz[j] * pa_z[j] + pc_xzzz[j] * pb_x[j]);

                t_xzzz_x[j] += -fl_s_0_0_5 * pc_xxzzz[j];

                t_xzzz_y[j] = fl_s_0_0_0 * (1.5 * pa_xz[j] * fl1_fx * pb_y[j] + pa_xzzz[j] * pb_y[j]);

                t_xzzz_y[j] += fl_s_0_0_1 * (-1.5 * pa_xz[j] * fl1_fx * pc_y[j] - 1.5 * pa_xz[j] * fl1_fx * pb_y[j] - 1.5 * pa_x[j] * pc_z[j] * fl1_fx * pb_y[j] - 1.5 * pc_x[j] * pa_z[j] * fl1_fx * pb_y[j] - pa_xzzz[j] * pc_y[j] - 3.0 * pa_xzz[j] * pc_z[j] * pb_y[j] - pc_x[j] * pa_zzz[j] * pb_y[j]);

                t_xzzz_y[j] += fl_s_0_0_2 * (1.5 * pa_xz[j] * fl1_fx * pc_y[j] + 1.5 * pa_x[j] * pc_yz[j] * fl1_fx + 1.5 * pa_x[j] * pc_z[j] * fl1_fx * pb_y[j] + 1.5 * pc_xy[j] * pa_z[j] * fl1_fx + 1.5 * pc_x[j] * pa_z[j] * fl1_fx * pb_y[j] + 1.5 * pc_xz[j] * fl1_fx * pb_y[j] + 3.0 * pa_xzz[j] * pc_yz[j] + 3.0 * pa_xz[j] * pc_zz[j] * pb_y[j] + pc_xy[j] * pa_zzz[j] + 3.0 * pc_xz[j] * pa_zz[j] * pb_y[j]);

                t_xzzz_y[j] += fl_s_0_0_3 * (-1.5 * pa_x[j] * pc_yz[j] * fl1_fx - 1.5 * pc_xy[j] * pa_z[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pb_y[j] - 3.0 * pa_xz[j] * pc_yzz[j] - pa_x[j] * pc_zzz[j] * pb_y[j] - 3.0 * pc_xyz[j] * pa_zz[j] - 3.0 * pc_xzz[j] * pa_z[j] * pb_y[j]);

                t_xzzz_y[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_x[j] * pc_yzzz[j] + 3.0 * pc_xyzz[j] * pa_z[j] + pc_xzzz[j] * pb_y[j]);

                t_xzzz_y[j] += -fl_s_0_0_5 * pc_xyzzz[j];

                t_xzzz_z[j] = fl_s_0_0_0 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pa_xzz[j] * fl1_fx + 1.5 * pa_xz[j] * fl1_fx * pb_z[j] + pa_xzzz[j] * pb_z[j]);

                t_xzzz_z[j] += fl_s_0_0_1 * (-1.5 * pa_x[j] * fl2_fx - 0.75 * pc_x[j] * fl2_fx - 1.5 * pa_xzz[j] * fl1_fx - 4.5 * pa_xz[j] * pc_z[j] * fl1_fx - 1.5 * pc_x[j] * pa_zz[j] * fl1_fx - 1.5 * pa_xz[j] * fl1_fx * pb_z[j] - 1.5 * pa_x[j] * pc_z[j] * fl1_fx * pb_z[j] - 1.5 * pc_x[j] * pa_z[j] * fl1_fx * pb_z[j] - pa_xzzz[j] * pc_z[j] - 3.0 * pa_xzz[j] * pc_z[j] * pb_z[j] - pc_x[j] * pa_zzz[j] * pb_z[j]);

                t_xzzz_z[j] += fl_s_0_0_2 * (0.75 * pa_x[j] * fl2_fx + 1.5 * pc_x[j] * fl2_fx + 4.5 * pa_xz[j] * pc_z[j] * fl1_fx + 3.0 * pa_x[j] * pc_zz[j] * fl1_fx + 1.5 * pc_x[j] * pa_zz[j] * fl1_fx + 4.5 * pc_xz[j] * pa_z[j] * fl1_fx + 1.5 * pa_x[j] * pc_z[j] * fl1_fx * pb_z[j] + 1.5 * pc_x[j] * pa_z[j] * fl1_fx * pb_z[j] + 1.5 * pc_xz[j] * fl1_fx * pb_z[j] + 3.0 * pa_xzz[j] * pc_zz[j] + 3.0 * pa_xz[j] * pc_zz[j] * pb_z[j] + pc_xz[j] * pa_zzz[j] + 3.0 * pc_xz[j] * pa_zz[j] * pb_z[j]);

                t_xzzz_z[j] += fl_s_0_0_3 * (-0.75 * pc_x[j] * fl2_fx - 3.0 * pa_x[j] * pc_zz[j] * fl1_fx - 4.5 * pc_xz[j] * pa_z[j] * fl1_fx - 3.0 * pc_xzz[j] * fl1_fx - 1.5 * pc_xz[j] * fl1_fx * pb_z[j] - 3.0 * pa_xz[j] * pc_zzz[j] - pa_x[j] * pc_zzz[j] * pb_z[j] - 3.0 * pc_xzz[j] * pa_zz[j] - 3.0 * pc_xzz[j] * pa_z[j] * pb_z[j]);

                t_xzzz_z[j] += fl_s_0_0_4 * (3.0 * pc_xzz[j] * fl1_fx + pa_x[j] * pc_zzzz[j] + 3.0 * pc_xzzz[j] * pa_z[j] + pc_xzzz[j] * pb_z[j]);

                t_xzzz_z[j] += -fl_s_0_0_5 * pc_xzzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP_30_35(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (30,35)

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

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xyyy = pcDistances.data(55 * idx + 25);

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_yyyy = pcDistances.data(55 * idx + 29);

            auto pc_yyyz = pcDistances.data(55 * idx + 30);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xyyyy = pcDistances.data(55 * idx + 44);

            auto pc_xyyyz = pcDistances.data(55 * idx + 45);

            auto pc_yyyyy = pcDistances.data(55 * idx + 49);

            auto pc_yyyyz = pcDistances.data(55 * idx + 50);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_yyyy_x = primBuffer.data(45 * idx + 30);

            auto t_yyyy_y = primBuffer.data(45 * idx + 31);

            auto t_yyyy_z = primBuffer.data(45 * idx + 32);

            auto t_yyyz_x = primBuffer.data(45 * idx + 33);

            auto t_yyyz_y = primBuffer.data(45 * idx + 34);

            // Batch of Integrals (30,35)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_y, pb_z, \
                                     pc_x, pc_xy, pc_xyy, pc_xyyy, pc_xyyyy, pc_xyyyz, pc_xyyz, pc_xyz, pc_xz, pc_y, \
                                     pc_yy, pc_yyy, pc_yyyy, pc_yyyyy, pc_yyyyz, pc_yyyz, pc_yyz, pc_yz, pc_z, s_0_0_0, \
                                     s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_yyyy_x, t_yyyy_y, t_yyyy_z, t_yyyz_x, \
                                     t_yyyz_y: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yyyy_x[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_x[j] + 3.0 * pa_yy[j] * fl1_fx * pb_x[j] + pa_yyyy[j] * pb_x[j]);

                t_yyyy_x[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_x[j] - 1.5 * fl2_fx * pb_x[j] - 3.0 * pa_yy[j] * fl1_fx * pc_x[j] - 3.0 * pa_yy[j] * fl1_fx * pb_x[j] - 6.0 * pa_y[j] * pc_y[j] * fl1_fx * pb_x[j] - pa_yyyy[j] * pc_x[j] - 4.0 * pa_yyy[j] * pc_y[j] * pb_x[j]);

                t_yyyy_x[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_x[j] + 0.75 * fl2_fx * pb_x[j] + 3.0 * pa_yy[j] * fl1_fx * pc_x[j] + 6.0 * pa_y[j] * pc_xy[j] * fl1_fx + 6.0 * pa_y[j] * pc_y[j] * fl1_fx * pb_x[j] + 3.0 * pc_yy[j] * fl1_fx * pb_x[j] + 4.0 * pa_yyy[j] * pc_xy[j] + 6.0 * pa_yy[j] * pc_yy[j] * pb_x[j]);

                t_yyyy_x[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_x[j] - 6.0 * pa_y[j] * pc_xy[j] * fl1_fx - 3.0 * pc_xyy[j] * fl1_fx - 3.0 * pc_yy[j] * fl1_fx * pb_x[j] - 6.0 * pa_yy[j] * pc_xyy[j] - 4.0 * pa_y[j] * pc_yyy[j] * pb_x[j]);

                t_yyyy_x[j] += fl_s_0_0_4 * (3.0 * pc_xyy[j] * fl1_fx + 4.0 * pa_y[j] * pc_xyyy[j] + pc_yyyy[j] * pb_x[j]);

                t_yyyy_x[j] += -fl_s_0_0_5 * pc_xyyyy[j];

                t_yyyy_y[j] = fl_s_0_0_0 * (3.0 * pa_y[j] * fl2_fx + 2.0 * pa_yyy[j] * fl1_fx + 0.75 * fl2_fx * pb_y[j] + 3.0 * pa_yy[j] * fl1_fx * pb_y[j] + pa_yyyy[j] * pb_y[j]);

                t_yyyy_y[j] += fl_s_0_0_1 * (-6.0 * pa_y[j] * fl2_fx - 3.75 * pc_y[j] * fl2_fx - 2.0 * pa_yyy[j] * fl1_fx - 9.0 * pa_yy[j] * pc_y[j] * fl1_fx - 1.5 * fl2_fx * pb_y[j] - 3.0 * pa_yy[j] * fl1_fx * pb_y[j] - 6.0 * pa_y[j] * pc_y[j] * fl1_fx * pb_y[j] - pa_yyyy[j] * pc_y[j] - 4.0 * pa_yyy[j] * pc_y[j] * pb_y[j]);

                t_yyyy_y[j] += fl_s_0_0_2 * (3.0 * pa_y[j] * fl2_fx + 7.5 * pc_y[j] * fl2_fx + 9.0 * pa_yy[j] * pc_y[j] * fl1_fx + 12.0 * pa_y[j] * pc_yy[j] * fl1_fx + 0.75 * fl2_fx * pb_y[j] + 6.0 * pa_y[j] * pc_y[j] * fl1_fx * pb_y[j] + 3.0 * pc_yy[j] * fl1_fx * pb_y[j] + 4.0 * pa_yyy[j] * pc_yy[j] + 6.0 * pa_yy[j] * pc_yy[j] * pb_y[j]);

                t_yyyy_y[j] += fl_s_0_0_3 * (-3.75 * pc_y[j] * fl2_fx - 12.0 * pa_y[j] * pc_yy[j] * fl1_fx - 5.0 * pc_yyy[j] * fl1_fx - 3.0 * pc_yy[j] * fl1_fx * pb_y[j] - 6.0 * pa_yy[j] * pc_yyy[j] - 4.0 * pa_y[j] * pc_yyy[j] * pb_y[j]);

                t_yyyy_y[j] += fl_s_0_0_4 * (5.0 * pc_yyy[j] * fl1_fx + 4.0 * pa_y[j] * pc_yyyy[j] + pc_yyyy[j] * pb_y[j]);

                t_yyyy_y[j] += -fl_s_0_0_5 * pc_yyyyy[j];

                t_yyyy_z[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_z[j] + 3.0 * pa_yy[j] * fl1_fx * pb_z[j] + pa_yyyy[j] * pb_z[j]);

                t_yyyy_z[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * fl2_fx * pb_z[j] - 3.0 * pa_yy[j] * fl1_fx * pc_z[j] - 3.0 * pa_yy[j] * fl1_fx * pb_z[j] - 6.0 * pa_y[j] * pc_y[j] * fl1_fx * pb_z[j] - pa_yyyy[j] * pc_z[j] - 4.0 * pa_yyy[j] * pc_y[j] * pb_z[j]);

                t_yyyy_z[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pb_z[j] + 3.0 * pa_yy[j] * fl1_fx * pc_z[j] + 6.0 * pa_y[j] * pc_yz[j] * fl1_fx + 6.0 * pa_y[j] * pc_y[j] * fl1_fx * pb_z[j] + 3.0 * pc_yy[j] * fl1_fx * pb_z[j] + 4.0 * pa_yyy[j] * pc_yz[j] + 6.0 * pa_yy[j] * pc_yy[j] * pb_z[j]);

                t_yyyy_z[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 6.0 * pa_y[j] * pc_yz[j] * fl1_fx - 3.0 * pc_yyz[j] * fl1_fx - 3.0 * pc_yy[j] * fl1_fx * pb_z[j] - 6.0 * pa_yy[j] * pc_yyz[j] - 4.0 * pa_y[j] * pc_yyy[j] * pb_z[j]);

                t_yyyy_z[j] += fl_s_0_0_4 * (3.0 * pc_yyz[j] * fl1_fx + 4.0 * pa_y[j] * pc_yyyz[j] + pc_yyyy[j] * pb_z[j]);

                t_yyyy_z[j] += -fl_s_0_0_5 * pc_yyyyz[j];

                t_yyyz_x[j] = fl_s_0_0_0 * (1.5 * pa_yz[j] * fl1_fx * pb_x[j] + pa_yyyz[j] * pb_x[j]);

                t_yyyz_x[j] += fl_s_0_0_1 * (-1.5 * pa_yz[j] * fl1_fx * pc_x[j] - 1.5 * pa_y[j] * fl1_fx * pc_z[j] * pb_x[j] - 1.5 * pa_yz[j] * fl1_fx * pb_x[j] - 1.5 * pc_y[j] * fl1_fx * pa_z[j] * pb_x[j] - pa_yyyz[j] * pc_x[j] - pa_yyy[j] * pc_z[j] * pb_x[j] - 3.0 * pa_yyz[j] * pc_y[j] * pb_x[j]);

                t_yyyz_x[j] += fl_s_0_0_2 * (1.5 * pa_y[j] * fl1_fx * pc_xz[j] + 1.5 * pa_yz[j] * fl1_fx * pc_x[j] + 1.5 * pa_y[j] * fl1_fx * pc_z[j] * pb_x[j] + 1.5 * pc_xy[j] * fl1_fx * pa_z[j] + 1.5 * pc_yz[j] * fl1_fx * pb_x[j] + 1.5 * pc_y[j] * fl1_fx * pa_z[j] * pb_x[j] + pa_yyy[j] * pc_xz[j] + 3.0 * pa_yyz[j] * pc_xy[j] + 3.0 * pa_yy[j] * pc_yz[j] * pb_x[j] + 3.0 * pa_yz[j] * pc_yy[j] * pb_x[j]);

                t_yyyz_x[j] += fl_s_0_0_3 * (-1.5 * pa_y[j] * fl1_fx * pc_xz[j] - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_xy[j] * fl1_fx * pa_z[j] - 1.5 * pc_yz[j] * fl1_fx * pb_x[j] - 3.0 * pa_yy[j] * pc_xyz[j] - 3.0 * pa_yz[j] * pc_xyy[j] - 3.0 * pa_y[j] * pc_yyz[j] * pb_x[j] - pc_yyy[j] * pa_z[j] * pb_x[j]);

                t_yyyz_x[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + 3.0 * pa_y[j] * pc_xyyz[j] + pc_xyyy[j] * pa_z[j] + pc_yyyz[j] * pb_x[j]);

                t_yyyz_x[j] += -fl_s_0_0_5 * pc_xyyyz[j];

                t_yyyz_y[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pa_z[j] + 1.5 * pa_yyz[j] * fl1_fx + 1.5 * pa_yz[j] * fl1_fx * pb_y[j] + pa_yyyz[j] * pb_y[j]);

                t_yyyz_y[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_z[j] - 1.5 * fl2_fx * pa_z[j] - 1.5 * pa_yy[j] * fl1_fx * pc_z[j] - 1.5 * pa_yyz[j] * fl1_fx - 4.5 * pa_yz[j] * pc_y[j] * fl1_fx - 1.5 * pa_y[j] * fl1_fx * pc_z[j] * pb_y[j] - 1.5 * pa_yz[j] * fl1_fx * pb_y[j] - 1.5 * pc_y[j] * fl1_fx * pa_z[j] * pb_y[j] - pa_yyyz[j] * pc_y[j] - pa_yyy[j] * pc_z[j] * pb_y[j] - 3.0 * pa_yyz[j] * pc_y[j] * pb_y[j]);

                t_yyyz_y[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pa_z[j] + 1.5 * pa_yy[j] * fl1_fx * pc_z[j] + 4.5 * pa_y[j] * pc_yz[j] * fl1_fx + 4.5 * pa_yz[j] * pc_y[j] * fl1_fx + 3.0 * pc_yy[j] * fl1_fx * pa_z[j] + 1.5 * pa_y[j] * fl1_fx * pc_z[j] * pb_y[j] + 1.5 * pc_yz[j] * fl1_fx * pb_y[j] + 1.5 * pc_y[j] * fl1_fx * pa_z[j] * pb_y[j] + pa_yyy[j] * pc_yz[j] + 3.0 * pa_yyz[j] * pc_yy[j] + 3.0 * pa_yy[j] * pc_yz[j] * pb_y[j] + 3.0 * pa_yz[j] * pc_yy[j] * pb_y[j]);

                t_yyyz_y[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 4.5 * pa_y[j] * pc_yz[j] * fl1_fx - 3.0 * pc_yyz[j] * fl1_fx - 3.0 * pc_yy[j] * fl1_fx * pa_z[j] - 1.5 * pc_yz[j] * fl1_fx * pb_y[j] - 3.0 * pa_yy[j] * pc_yyz[j] - 3.0 * pa_yz[j] * pc_yyy[j] - 3.0 * pa_y[j] * pc_yyz[j] * pb_y[j] - pc_yyy[j] * pa_z[j] * pb_y[j]);

                t_yyyz_y[j] += fl_s_0_0_4 * (3.0 * pc_yyz[j] * fl1_fx + 3.0 * pa_y[j] * pc_yyyz[j] + pc_yyyy[j] * pa_z[j] + pc_yyyz[j] * pb_y[j]);

                t_yyyz_y[j] += -fl_s_0_0_5 * pc_yyyyz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP_35_40(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (35,40)

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

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xy = pcDistances.data(55 * idx + 4);

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xyy = pcDistances.data(55 * idx + 12);

            auto pc_xyz = pcDistances.data(55 * idx + 13);

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyy = pcDistances.data(55 * idx + 15);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xyyz = pcDistances.data(55 * idx + 26);

            auto pc_xyzz = pcDistances.data(55 * idx + 27);

            auto pc_xzzz = pcDistances.data(55 * idx + 28);

            auto pc_yyyz = pcDistances.data(55 * idx + 30);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            auto pc_yzzz = pcDistances.data(55 * idx + 32);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xyyzz = pcDistances.data(55 * idx + 46);

            auto pc_xyzzz = pcDistances.data(55 * idx + 47);

            auto pc_yyyzz = pcDistances.data(55 * idx + 51);

            auto pc_yyzzz = pcDistances.data(55 * idx + 52);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_yyyz_z = primBuffer.data(45 * idx + 35);

            auto t_yyzz_x = primBuffer.data(45 * idx + 36);

            auto t_yyzz_y = primBuffer.data(45 * idx + 37);

            auto t_yyzz_z = primBuffer.data(45 * idx + 38);

            auto t_yzzz_x = primBuffer.data(45 * idx + 39);

            // Batch of Integrals (35,40)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_yzzz, pa_z, \
                                     pa_zz, pa_zzz, pb_x, pb_y, pb_z, pc_x, pc_xy, pc_xyy, pc_xyyz, pc_xyyzz, pc_xyz, \
                                     pc_xyzz, pc_xyzzz, pc_xz, pc_xzz, pc_xzzz, pc_y, pc_yy, pc_yyy, pc_yyyz, pc_yyyzz, \
                                     pc_yyz, pc_yyzz, pc_yyzzz, pc_yz, pc_yzz, pc_yzzz, pc_z, pc_zz, pc_zzz, s_0_0_0, \
                                     s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_yyyz_z, t_yyzz_x, t_yyzz_y, t_yyzz_z, \
                                     t_yzzz_x: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yyyz_z[j] = fl_s_0_0_0 * (0.75 * pa_y[j] * fl2_fx + 0.5 * pa_yyy[j] * fl1_fx + 1.5 * pa_yz[j] * fl1_fx * pb_z[j] + pa_yyyz[j] * pb_z[j]);

                t_yyyz_z[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl2_fx - 0.75 * pc_y[j] * fl2_fx - 0.5 * pa_yyy[j] * fl1_fx - 1.5 * pa_yy[j] * pc_y[j] * fl1_fx - 1.5 * pa_yz[j] * fl1_fx * pc_z[j] - 1.5 * pa_y[j] * fl1_fx * pc_z[j] * pb_z[j] - 1.5 * pa_yz[j] * fl1_fx * pb_z[j] - 1.5 * pc_y[j] * fl1_fx * pa_z[j] * pb_z[j] - pa_yyyz[j] * pc_z[j] - pa_yyy[j] * pc_z[j] * pb_z[j] - 3.0 * pa_yyz[j] * pc_y[j] * pb_z[j]);

                t_yyyz_z[j] += fl_s_0_0_2 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pc_y[j] * fl2_fx + 1.5 * pa_yy[j] * pc_y[j] * fl1_fx + 1.5 * pa_y[j] * pc_yy[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pc_zz[j] + 1.5 * pa_yz[j] * fl1_fx * pc_z[j] + 1.5 * pa_y[j] * fl1_fx * pc_z[j] * pb_z[j] + 1.5 * pc_yz[j] * fl1_fx * pa_z[j] + 1.5 * pc_yz[j] * fl1_fx * pb_z[j] + 1.5 * pc_y[j] * fl1_fx * pa_z[j] * pb_z[j] + pa_yyy[j] * pc_zz[j] + 3.0 * pa_yyz[j] * pc_yz[j] + 3.0 * pa_yy[j] * pc_yz[j] * pb_z[j] + 3.0 * pa_yz[j] * pc_yy[j] * pb_z[j]);

                t_yyyz_z[j] += fl_s_0_0_3 * (-0.75 * pc_y[j] * fl2_fx - 1.5 * pa_y[j] * pc_yy[j] * fl1_fx - 0.5 * pc_yyy[j] * fl1_fx - 1.5 * pa_y[j] * fl1_fx * pc_zz[j] - 1.5 * pc_yzz[j] * fl1_fx - 1.5 * pc_yz[j] * fl1_fx * pa_z[j] - 1.5 * pc_yz[j] * fl1_fx * pb_z[j] - 3.0 * pa_yy[j] * pc_yzz[j] - 3.0 * pa_yz[j] * pc_yyz[j] - 3.0 * pa_y[j] * pc_yyz[j] * pb_z[j] - pc_yyy[j] * pa_z[j] * pb_z[j]);

                t_yyyz_z[j] += fl_s_0_0_4 * (0.5 * pc_yyy[j] * fl1_fx + 1.5 * pc_yzz[j] * fl1_fx + 3.0 * pa_y[j] * pc_yyzz[j] + pc_yyyz[j] * pa_z[j] + pc_yyyz[j] * pb_z[j]);

                t_yyyz_z[j] += -fl_s_0_0_5 * pc_yyyzz[j];

                t_yyzz_x[j] = fl_s_0_0_0 * (0.25 * fl2_fx * pb_x[j] + 0.5 * pa_yy[j] * fl1_fx * pb_x[j] + 0.5 * fl1_fx * pa_zz[j] * pb_x[j] + pa_yyzz[j] * pb_x[j]);

                t_yyzz_x[j] += fl_s_0_0_1 * (-0.25 * fl2_fx * pc_x[j] - 0.5 * fl2_fx * pb_x[j] - 0.5 * pa_yy[j] * fl1_fx * pc_x[j] - 0.5 * pa_yy[j] * fl1_fx * pb_x[j] - pa_y[j] * pc_y[j] * fl1_fx * pb_x[j] - 0.5 * fl1_fx * pa_zz[j] * pc_x[j] - fl1_fx * pa_z[j] * pc_z[j] * pb_x[j] - 0.5 * fl1_fx * pa_zz[j] * pb_x[j] - pa_yyzz[j] * pc_x[j] - 2.0 * pa_yyz[j] * pc_z[j] * pb_x[j] - 2.0 * pa_yzz[j] * pc_y[j] * pb_x[j]);

                t_yyzz_x[j] += fl_s_0_0_2 * (0.5 * fl2_fx * pc_x[j] + 0.25 * fl2_fx * pb_x[j] + 0.5 * pa_yy[j] * fl1_fx * pc_x[j] + pa_y[j] * pc_xy[j] * fl1_fx + pa_y[j] * pc_y[j] * fl1_fx * pb_x[j] + 0.5 * pc_yy[j] * fl1_fx * pb_x[j] + fl1_fx * pa_z[j] * pc_xz[j] + 0.5 * fl1_fx * pc_zz[j] * pb_x[j] + 0.5 * fl1_fx * pa_zz[j] * pc_x[j] + fl1_fx * pa_z[j] * pc_z[j] * pb_x[j] + 2.0 * pa_yyz[j] * pc_xz[j] + pa_yy[j] * pc_zz[j] * pb_x[j] + 2.0 * pa_yzz[j] * pc_xy[j] + 4.0 * pa_yz[j] * pc_yz[j] * pb_x[j] + pc_yy[j] * pa_zz[j] * pb_x[j]);

                t_yyzz_x[j] += fl_s_0_0_3 * (-0.25 * fl2_fx * pc_x[j] - pa_y[j] * pc_xy[j] * fl1_fx - 0.5 * pc_xyy[j] * fl1_fx - 0.5 * pc_yy[j] * fl1_fx * pb_x[j] - 0.5 * fl1_fx * pc_xzz[j] - fl1_fx * pa_z[j] * pc_xz[j] - 0.5 * fl1_fx * pc_zz[j] * pb_x[j] - pa_yy[j] * pc_xzz[j] - 4.0 * pa_yz[j] * pc_xyz[j] - 2.0 * pa_y[j] * pc_yzz[j] * pb_x[j] - pc_xyy[j] * pa_zz[j] - 2.0 * pc_yyz[j] * pa_z[j] * pb_x[j]);

                t_yyzz_x[j] += fl_s_0_0_4 * (0.5 * pc_xyy[j] * fl1_fx + 0.5 * fl1_fx * pc_xzz[j] + 2.0 * pa_y[j] * pc_xyzz[j] + 2.0 * pc_xyyz[j] * pa_z[j] + pc_yyzz[j] * pb_x[j]);

                t_yyzz_x[j] += -fl_s_0_0_5 * pc_xyyzz[j];

                t_yyzz_y[j] = fl_s_0_0_0 * (0.5 * pa_y[j] * fl2_fx + pa_yzz[j] * fl1_fx + 0.25 * fl2_fx * pb_y[j] + 0.5 * pa_yy[j] * fl1_fx * pb_y[j] + 0.5 * fl1_fx * pa_zz[j] * pb_y[j] + pa_yyzz[j] * pb_y[j]);

                t_yyzz_y[j] += fl_s_0_0_1 * (-pa_y[j] * fl2_fx - 0.75 * pc_y[j] * fl2_fx - 2.0 * pa_yz[j] * fl1_fx * pc_z[j] - pa_yzz[j] * fl1_fx - 1.5 * pc_y[j] * fl1_fx * pa_zz[j] - 0.5 * fl2_fx * pb_y[j] - 0.5 * pa_yy[j] * fl1_fx * pc_y[j] - 0.5 * pa_yy[j] * fl1_fx * pb_y[j] - pa_y[j] * pc_y[j] * fl1_fx * pb_y[j] - fl1_fx * pa_z[j] * pc_z[j] * pb_y[j] - 0.5 * fl1_fx * pa_zz[j] * pb_y[j] - pa_yyzz[j] * pc_y[j] - 2.0 * pa_yyz[j] * pc_z[j] * pb_y[j] - 2.0 * pa_yzz[j] * pc_y[j] * pb_y[j]);

                t_yyzz_y[j] += fl_s_0_0_2 * (0.5 * pa_y[j] * fl2_fx + 1.5 * pc_y[j] * fl2_fx + pa_y[j] * fl1_fx * pc_zz[j] + 2.0 * pa_yz[j] * fl1_fx * pc_z[j] + 3.0 * pc_yz[j] * fl1_fx * pa_z[j] + 1.5 * pc_y[j] * fl1_fx * pa_zz[j] + 0.25 * fl2_fx * pb_y[j] + 0.5 * pa_yy[j] * fl1_fx * pc_y[j] + pa_y[j] * pc_yy[j] * fl1_fx + pa_y[j] * pc_y[j] * fl1_fx * pb_y[j] + 0.5 * pc_yy[j] * fl1_fx * pb_y[j] + 0.5 * fl1_fx * pc_zz[j] * pb_y[j] + fl1_fx * pa_z[j] * pc_z[j] * pb_y[j] + 2.0 * pa_yyz[j] * pc_yz[j] + pa_yy[j] * pc_zz[j] * pb_y[j] + 2.0 * pa_yzz[j] * pc_yy[j] + 4.0 * pa_yz[j] * pc_yz[j] * pb_y[j] + pc_yy[j] * pa_zz[j] * pb_y[j]);

                t_yyzz_y[j] += fl_s_0_0_3 * (-0.75 * pc_y[j] * fl2_fx - pa_y[j] * fl1_fx * pc_zz[j] - 1.5 * pc_yzz[j] * fl1_fx - 3.0 * pc_yz[j] * fl1_fx * pa_z[j] - pa_y[j] * pc_yy[j] * fl1_fx - 0.5 * pc_yyy[j] * fl1_fx - 0.5 * pc_yy[j] * fl1_fx * pb_y[j] - 0.5 * fl1_fx * pc_zz[j] * pb_y[j] - pa_yy[j] * pc_yzz[j] - 4.0 * pa_yz[j] * pc_yyz[j] - 2.0 * pa_y[j] * pc_yzz[j] * pb_y[j] - pc_yyy[j] * pa_zz[j] - 2.0 * pc_yyz[j] * pa_z[j] * pb_y[j]);

                t_yyzz_y[j] += fl_s_0_0_4 * (1.5 * pc_yzz[j] * fl1_fx + 0.5 * pc_yyy[j] * fl1_fx + 2.0 * pa_y[j] * pc_yyzz[j] + 2.0 * pc_yyyz[j] * pa_z[j] + pc_yyzz[j] * pb_y[j]);

                t_yyzz_y[j] += -fl_s_0_0_5 * pc_yyyzz[j];

                t_yyzz_z[j] = fl_s_0_0_0 * (0.5 * fl2_fx * pa_z[j] + pa_yyz[j] * fl1_fx + 0.25 * fl2_fx * pb_z[j] + 0.5 * pa_yy[j] * fl1_fx * pb_z[j] + 0.5 * fl1_fx * pa_zz[j] * pb_z[j] + pa_yyzz[j] * pb_z[j]);

                t_yyzz_z[j] += fl_s_0_0_1 * (-fl2_fx * pa_z[j] - 0.75 * fl2_fx * pc_z[j] - pa_yyz[j] * fl1_fx - 1.5 * pa_yy[j] * pc_z[j] * fl1_fx - 2.0 * pa_yz[j] * pc_y[j] * fl1_fx - 0.5 * fl2_fx * pb_z[j] - 0.5 * pa_yy[j] * fl1_fx * pb_z[j] - pa_y[j] * pc_y[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pa_zz[j] * pc_z[j] - fl1_fx * pa_z[j] * pc_z[j] * pb_z[j] - 0.5 * fl1_fx * pa_zz[j] * pb_z[j] - pa_yyzz[j] * pc_z[j] - 2.0 * pa_yyz[j] * pc_z[j] * pb_z[j] - 2.0 * pa_yzz[j] * pc_y[j] * pb_z[j]);

                t_yyzz_z[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.5 * fl2_fx * pa_z[j] + 1.5 * pa_yy[j] * pc_z[j] * fl1_fx + 2.0 * pa_yz[j] * pc_y[j] * fl1_fx + 3.0 * pa_y[j] * pc_yz[j] * fl1_fx + pc_yy[j] * pa_z[j] * fl1_fx + 0.25 * fl2_fx * pb_z[j] + pa_y[j] * pc_y[j] * fl1_fx * pb_z[j] + 0.5 * pc_yy[j] * fl1_fx * pb_z[j] + fl1_fx * pa_z[j] * pc_zz[j] + 0.5 * fl1_fx * pc_zz[j] * pb_z[j] + 0.5 * fl1_fx * pa_zz[j] * pc_z[j] + fl1_fx * pa_z[j] * pc_z[j] * pb_z[j] + 2.0 * pa_yyz[j] * pc_zz[j] + pa_yy[j] * pc_zz[j] * pb_z[j] + 2.0 * pa_yzz[j] * pc_yz[j] + 4.0 * pa_yz[j] * pc_yz[j] * pb_z[j] + pc_yy[j] * pa_zz[j] * pb_z[j]);

                t_yyzz_z[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 3.0 * pa_y[j] * pc_yz[j] * fl1_fx - pc_yy[j] * pa_z[j] * fl1_fx - 1.5 * pc_yyz[j] * fl1_fx - 0.5 * pc_yy[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fx * pc_zzz[j] - fl1_fx * pa_z[j] * pc_zz[j] - 0.5 * fl1_fx * pc_zz[j] * pb_z[j] - pa_yy[j] * pc_zzz[j] - 4.0 * pa_yz[j] * pc_yzz[j] - 2.0 * pa_y[j] * pc_yzz[j] * pb_z[j] - pc_yyz[j] * pa_zz[j] - 2.0 * pc_yyz[j] * pa_z[j] * pb_z[j]);

                t_yyzz_z[j] += fl_s_0_0_4 * (1.5 * pc_yyz[j] * fl1_fx + 0.5 * fl1_fx * pc_zzz[j] + 2.0 * pa_y[j] * pc_yzzz[j] + 2.0 * pc_yyzz[j] * pa_z[j] + pc_yyzz[j] * pb_z[j]);

                t_yyzz_z[j] += -fl_s_0_0_5 * pc_yyzzz[j];

                t_yzzz_x[j] = fl_s_0_0_0 * (1.5 * pa_yz[j] * fl1_fx * pb_x[j] + pa_yzzz[j] * pb_x[j]);

                t_yzzz_x[j] += fl_s_0_0_1 * (-1.5 * pa_yz[j] * fl1_fx * pc_x[j] - 1.5 * pa_yz[j] * fl1_fx * pb_x[j] - 1.5 * pa_y[j] * pc_z[j] * fl1_fx * pb_x[j] - 1.5 * pc_y[j] * pa_z[j] * fl1_fx * pb_x[j] - pa_yzzz[j] * pc_x[j] - 3.0 * pa_yzz[j] * pc_z[j] * pb_x[j] - pc_y[j] * pa_zzz[j] * pb_x[j]);

                t_yzzz_x[j] += fl_s_0_0_2 * (1.5 * pa_yz[j] * fl1_fx * pc_x[j] + 1.5 * pa_y[j] * pc_xz[j] * fl1_fx + 1.5 * pa_y[j] * pc_z[j] * fl1_fx * pb_x[j] + 1.5 * pc_xy[j] * pa_z[j] * fl1_fx + 1.5 * pc_y[j] * pa_z[j] * fl1_fx * pb_x[j] + 1.5 * pc_yz[j] * fl1_fx * pb_x[j] + 3.0 * pa_yzz[j] * pc_xz[j] + 3.0 * pa_yz[j] * pc_zz[j] * pb_x[j] + pc_xy[j] * pa_zzz[j] + 3.0 * pc_yz[j] * pa_zz[j] * pb_x[j]);

                t_yzzz_x[j] += fl_s_0_0_3 * (-1.5 * pa_y[j] * pc_xz[j] * fl1_fx - 1.5 * pc_xy[j] * pa_z[j] * fl1_fx - 1.5 * pc_xyz[j] * fl1_fx - 1.5 * pc_yz[j] * fl1_fx * pb_x[j] - 3.0 * pa_yz[j] * pc_xzz[j] - pa_y[j] * pc_zzz[j] * pb_x[j] - 3.0 * pc_xyz[j] * pa_zz[j] - 3.0 * pc_yzz[j] * pa_z[j] * pb_x[j]);

                t_yzzz_x[j] += fl_s_0_0_4 * (1.5 * pc_xyz[j] * fl1_fx + pa_y[j] * pc_xzzz[j] + 3.0 * pc_xyzz[j] * pa_z[j] + pc_yzzz[j] * pb_x[j]);

                t_yzzz_x[j] += -fl_s_0_0_5 * pc_xyzzz[j];
            }

            idx++;
        }
    }

    void
    compNuclearPotentialForGP_40_45(      CMemBlock2D<double>& primBuffer,
                                    const CMemBlock2D<double>& auxBuffer,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pbDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (40,45)

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

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PC)

            auto pc_x = pcDistances.data(55 * idx);

            auto pc_y = pcDistances.data(55 * idx + 1);

            auto pc_z = pcDistances.data(55 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PC)

            auto pc_xz = pcDistances.data(55 * idx + 5);

            auto pc_yy = pcDistances.data(55 * idx + 6);

            auto pc_yz = pcDistances.data(55 * idx + 7);

            auto pc_zz = pcDistances.data(55 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PC)

            auto pc_xzz = pcDistances.data(55 * idx + 14);

            auto pc_yyz = pcDistances.data(55 * idx + 16);

            auto pc_yzz = pcDistances.data(55 * idx + 17);

            auto pc_zzz = pcDistances.data(55 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PC)

            auto pc_xzzz = pcDistances.data(55 * idx + 28);

            auto pc_yyzz = pcDistances.data(55 * idx + 31);

            auto pc_yzzz = pcDistances.data(55 * idx + 32);

            auto pc_zzzz = pcDistances.data(55 * idx + 33);

            // set up pointers to 5-th order tensor of distance R(PC)

            auto pc_xzzzz = pcDistances.data(55 * idx + 48);

            auto pc_yyzzz = pcDistances.data(55 * idx + 52);

            auto pc_yzzzz = pcDistances.data(55 * idx + 53);

            auto pc_zzzzz = pcDistances.data(55 * idx + 54);

            // set up pointers to auxilary integrals

            auto s_0_0_0 = auxBuffer.data(6 * idx);

            auto s_0_0_1 = auxBuffer.data(6 * idx + 1);

            auto s_0_0_2 = auxBuffer.data(6 * idx + 2);

            auto s_0_0_3 = auxBuffer.data(6 * idx + 3);

            auto s_0_0_4 = auxBuffer.data(6 * idx + 4);

            auto s_0_0_5 = auxBuffer.data(6 * idx + 5);

            // set up pointers to integrals

            auto t_yzzz_y = primBuffer.data(45 * idx + 40);

            auto t_yzzz_z = primBuffer.data(45 * idx + 41);

            auto t_zzzz_x = primBuffer.data(45 * idx + 42);

            auto t_zzzz_y = primBuffer.data(45 * idx + 43);

            auto t_zzzz_z = primBuffer.data(45 * idx + 44);

            // Batch of Integrals (40,45)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_y, pb_z, \
                                     pc_x, pc_xz, pc_xzz, pc_xzzz, pc_xzzzz, pc_y, pc_yy, pc_yyz, pc_yyzz, pc_yyzzz, \
                                     pc_yz, pc_yzz, pc_yzzz, pc_yzzzz, pc_z, pc_zz, pc_zzz, pc_zzzz, pc_zzzzz, s_0_0_0, \
                                     s_0_0_1, s_0_0_2, s_0_0_3, s_0_0_4, s_0_0_5, t_yzzz_y, t_yzzz_z, t_zzzz_x, t_zzzz_y, \
                                     t_zzzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0_0 = s_0_0_0[j];

                double fl_s_0_0_1 = s_0_0_1[j];

                double fl_s_0_0_2 = s_0_0_2[j];

                double fl_s_0_0_3 = s_0_0_3[j];

                double fl_s_0_0_4 = s_0_0_4[j];

                double fl_s_0_0_5 = s_0_0_5[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_yzzz_y[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pa_z[j] + 0.5 * fl1_fx * pa_zzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_y[j] + pa_yzzz[j] * pb_y[j]);

                t_yzzz_y[j] += fl_s_0_0_1 * (-1.5 * fl2_fx * pa_z[j] - 0.75 * fl2_fx * pc_z[j] - 1.5 * fl1_fx * pa_zz[j] * pc_z[j] - 0.5 * fl1_fx * pa_zzz[j] - 1.5 * pa_yz[j] * fl1_fx * pc_y[j] - 1.5 * pa_yz[j] * fl1_fx * pb_y[j] - 1.5 * pa_y[j] * pc_z[j] * fl1_fx * pb_y[j] - 1.5 * pc_y[j] * pa_z[j] * fl1_fx * pb_y[j] - pa_yzzz[j] * pc_y[j] - 3.0 * pa_yzz[j] * pc_z[j] * pb_y[j] - pc_y[j] * pa_zzz[j] * pb_y[j]);

                t_yzzz_y[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_z[j] + 0.75 * fl2_fx * pa_z[j] + 1.5 * fl1_fx * pa_z[j] * pc_zz[j] + 1.5 * fl1_fx * pa_zz[j] * pc_z[j] + 1.5 * pa_yz[j] * fl1_fx * pc_y[j] + 1.5 * pa_y[j] * pc_yz[j] * fl1_fx + 1.5 * pa_y[j] * pc_z[j] * fl1_fx * pb_y[j] + 1.5 * pc_yy[j] * pa_z[j] * fl1_fx + 1.5 * pc_y[j] * pa_z[j] * fl1_fx * pb_y[j] + 1.5 * pc_yz[j] * fl1_fx * pb_y[j] + 3.0 * pa_yzz[j] * pc_yz[j] + 3.0 * pa_yz[j] * pc_zz[j] * pb_y[j] + pc_yy[j] * pa_zzz[j] + 3.0 * pc_yz[j] * pa_zz[j] * pb_y[j]);

                t_yzzz_y[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_z[j] - 0.5 * fl1_fx * pc_zzz[j] - 1.5 * fl1_fx * pa_z[j] * pc_zz[j] - 1.5 * pa_y[j] * pc_yz[j] * fl1_fx - 1.5 * pc_yy[j] * pa_z[j] * fl1_fx - 1.5 * pc_yyz[j] * fl1_fx - 1.5 * pc_yz[j] * fl1_fx * pb_y[j] - 3.0 * pa_yz[j] * pc_yzz[j] - pa_y[j] * pc_zzz[j] * pb_y[j] - 3.0 * pc_yyz[j] * pa_zz[j] - 3.0 * pc_yzz[j] * pa_z[j] * pb_y[j]);

                t_yzzz_y[j] += fl_s_0_0_4 * (0.5 * fl1_fx * pc_zzz[j] + 1.5 * pc_yyz[j] * fl1_fx + pa_y[j] * pc_yzzz[j] + 3.0 * pc_yyzz[j] * pa_z[j] + pc_yzzz[j] * pb_y[j]);

                t_yzzz_y[j] += -fl_s_0_0_5 * pc_yyzzz[j];

                t_yzzz_z[j] = fl_s_0_0_0 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pa_yzz[j] * fl1_fx + 1.5 * pa_yz[j] * fl1_fx * pb_z[j] + pa_yzzz[j] * pb_z[j]);

                t_yzzz_z[j] += fl_s_0_0_1 * (-1.5 * pa_y[j] * fl2_fx - 0.75 * pc_y[j] * fl2_fx - 1.5 * pa_yzz[j] * fl1_fx - 4.5 * pa_yz[j] * pc_z[j] * fl1_fx - 1.5 * pc_y[j] * pa_zz[j] * fl1_fx - 1.5 * pa_yz[j] * fl1_fx * pb_z[j] - 1.5 * pa_y[j] * pc_z[j] * fl1_fx * pb_z[j] - 1.5 * pc_y[j] * pa_z[j] * fl1_fx * pb_z[j] - pa_yzzz[j] * pc_z[j] - 3.0 * pa_yzz[j] * pc_z[j] * pb_z[j] - pc_y[j] * pa_zzz[j] * pb_z[j]);

                t_yzzz_z[j] += fl_s_0_0_2 * (0.75 * pa_y[j] * fl2_fx + 1.5 * pc_y[j] * fl2_fx + 4.5 * pa_yz[j] * pc_z[j] * fl1_fx + 3.0 * pa_y[j] * pc_zz[j] * fl1_fx + 1.5 * pc_y[j] * pa_zz[j] * fl1_fx + 4.5 * pc_yz[j] * pa_z[j] * fl1_fx + 1.5 * pa_y[j] * pc_z[j] * fl1_fx * pb_z[j] + 1.5 * pc_y[j] * pa_z[j] * fl1_fx * pb_z[j] + 1.5 * pc_yz[j] * fl1_fx * pb_z[j] + 3.0 * pa_yzz[j] * pc_zz[j] + 3.0 * pa_yz[j] * pc_zz[j] * pb_z[j] + pc_yz[j] * pa_zzz[j] + 3.0 * pc_yz[j] * pa_zz[j] * pb_z[j]);

                t_yzzz_z[j] += fl_s_0_0_3 * (-0.75 * pc_y[j] * fl2_fx - 3.0 * pa_y[j] * pc_zz[j] * fl1_fx - 4.5 * pc_yz[j] * pa_z[j] * fl1_fx - 3.0 * pc_yzz[j] * fl1_fx - 1.5 * pc_yz[j] * fl1_fx * pb_z[j] - 3.0 * pa_yz[j] * pc_zzz[j] - pa_y[j] * pc_zzz[j] * pb_z[j] - 3.0 * pc_yzz[j] * pa_zz[j] - 3.0 * pc_yzz[j] * pa_z[j] * pb_z[j]);

                t_yzzz_z[j] += fl_s_0_0_4 * (3.0 * pc_yzz[j] * fl1_fx + pa_y[j] * pc_zzzz[j] + 3.0 * pc_yzzz[j] * pa_z[j] + pc_yzzz[j] * pb_z[j]);

                t_yzzz_z[j] += -fl_s_0_0_5 * pc_yzzzz[j];

                t_zzzz_x[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_x[j] + 3.0 * pa_zz[j] * fl1_fx * pb_x[j] + pa_zzzz[j] * pb_x[j]);

                t_zzzz_x[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_x[j] - 1.5 * fl2_fx * pb_x[j] - 3.0 * pa_zz[j] * fl1_fx * pc_x[j] - 3.0 * pa_zz[j] * fl1_fx * pb_x[j] - 6.0 * pa_z[j] * pc_z[j] * fl1_fx * pb_x[j] - pa_zzzz[j] * pc_x[j] - 4.0 * pa_zzz[j] * pc_z[j] * pb_x[j]);

                t_zzzz_x[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_x[j] + 0.75 * fl2_fx * pb_x[j] + 3.0 * pa_zz[j] * fl1_fx * pc_x[j] + 6.0 * pa_z[j] * pc_xz[j] * fl1_fx + 6.0 * pa_z[j] * pc_z[j] * fl1_fx * pb_x[j] + 3.0 * pc_zz[j] * fl1_fx * pb_x[j] + 4.0 * pa_zzz[j] * pc_xz[j] + 6.0 * pa_zz[j] * pc_zz[j] * pb_x[j]);

                t_zzzz_x[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_x[j] - 6.0 * pa_z[j] * pc_xz[j] * fl1_fx - 3.0 * pc_xzz[j] * fl1_fx - 3.0 * pc_zz[j] * fl1_fx * pb_x[j] - 6.0 * pa_zz[j] * pc_xzz[j] - 4.0 * pa_z[j] * pc_zzz[j] * pb_x[j]);

                t_zzzz_x[j] += fl_s_0_0_4 * (3.0 * pc_xzz[j] * fl1_fx + 4.0 * pa_z[j] * pc_xzzz[j] + pc_zzzz[j] * pb_x[j]);

                t_zzzz_x[j] += -fl_s_0_0_5 * pc_xzzzz[j];

                t_zzzz_y[j] = fl_s_0_0_0 * (0.75 * fl2_fx * pb_y[j] + 3.0 * pa_zz[j] * fl1_fx * pb_y[j] + pa_zzzz[j] * pb_y[j]);

                t_zzzz_y[j] += fl_s_0_0_1 * (-0.75 * fl2_fx * pc_y[j] - 1.5 * fl2_fx * pb_y[j] - 3.0 * pa_zz[j] * fl1_fx * pc_y[j] - 3.0 * pa_zz[j] * fl1_fx * pb_y[j] - 6.0 * pa_z[j] * pc_z[j] * fl1_fx * pb_y[j] - pa_zzzz[j] * pc_y[j] - 4.0 * pa_zzz[j] * pc_z[j] * pb_y[j]);

                t_zzzz_y[j] += fl_s_0_0_2 * (1.5 * fl2_fx * pc_y[j] + 0.75 * fl2_fx * pb_y[j] + 3.0 * pa_zz[j] * fl1_fx * pc_y[j] + 6.0 * pa_z[j] * pc_yz[j] * fl1_fx + 6.0 * pa_z[j] * pc_z[j] * fl1_fx * pb_y[j] + 3.0 * pc_zz[j] * fl1_fx * pb_y[j] + 4.0 * pa_zzz[j] * pc_yz[j] + 6.0 * pa_zz[j] * pc_zz[j] * pb_y[j]);

                t_zzzz_y[j] += fl_s_0_0_3 * (-0.75 * fl2_fx * pc_y[j] - 6.0 * pa_z[j] * pc_yz[j] * fl1_fx - 3.0 * pc_yzz[j] * fl1_fx - 3.0 * pc_zz[j] * fl1_fx * pb_y[j] - 6.0 * pa_zz[j] * pc_yzz[j] - 4.0 * pa_z[j] * pc_zzz[j] * pb_y[j]);

                t_zzzz_y[j] += fl_s_0_0_4 * (3.0 * pc_yzz[j] * fl1_fx + 4.0 * pa_z[j] * pc_yzzz[j] + pc_zzzz[j] * pb_y[j]);

                t_zzzz_y[j] += -fl_s_0_0_5 * pc_yzzzz[j];

                t_zzzz_z[j] = fl_s_0_0_0 * (3.0 * pa_z[j] * fl2_fx + 2.0 * pa_zzz[j] * fl1_fx + 0.75 * fl2_fx * pb_z[j] + 3.0 * pa_zz[j] * fl1_fx * pb_z[j] + pa_zzzz[j] * pb_z[j]);

                t_zzzz_z[j] += fl_s_0_0_1 * (-6.0 * pa_z[j] * fl2_fx - 3.75 * pc_z[j] * fl2_fx - 2.0 * pa_zzz[j] * fl1_fx - 9.0 * pa_zz[j] * pc_z[j] * fl1_fx - 1.5 * fl2_fx * pb_z[j] - 3.0 * pa_zz[j] * fl1_fx * pb_z[j] - 6.0 * pa_z[j] * pc_z[j] * fl1_fx * pb_z[j] - pa_zzzz[j] * pc_z[j] - 4.0 * pa_zzz[j] * pc_z[j] * pb_z[j]);

                t_zzzz_z[j] += fl_s_0_0_2 * (3.0 * pa_z[j] * fl2_fx + 7.5 * pc_z[j] * fl2_fx + 9.0 * pa_zz[j] * pc_z[j] * fl1_fx + 12.0 * pa_z[j] * pc_zz[j] * fl1_fx + 0.75 * fl2_fx * pb_z[j] + 6.0 * pa_z[j] * pc_z[j] * fl1_fx * pb_z[j] + 3.0 * pc_zz[j] * fl1_fx * pb_z[j] + 4.0 * pa_zzz[j] * pc_zz[j] + 6.0 * pa_zz[j] * pc_zz[j] * pb_z[j]);

                t_zzzz_z[j] += fl_s_0_0_3 * (-3.75 * pc_z[j] * fl2_fx - 12.0 * pa_z[j] * pc_zz[j] * fl1_fx - 5.0 * pc_zzz[j] * fl1_fx - 3.0 * pc_zz[j] * fl1_fx * pb_z[j] - 6.0 * pa_zz[j] * pc_zzz[j] - 4.0 * pa_z[j] * pc_zzz[j] * pb_z[j]);

                t_zzzz_z[j] += fl_s_0_0_4 * (5.0 * pc_zzz[j] * fl1_fx + 4.0 * pa_z[j] * pc_zzzz[j] + pc_zzzz[j] * pb_z[j]);

                t_zzzz_z[j] += -fl_s_0_0_5 * pc_zzzzz[j];
            }

            idx++;
        }
    }


} // npotrecfunc namespace

