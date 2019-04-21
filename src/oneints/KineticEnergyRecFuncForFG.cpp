//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForFG.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForFG(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForFG_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_30_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_40_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_50_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_60_70(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_70_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_80_90(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_90_100(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                  braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_100_110(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_110_120(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_120_130(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_130_140(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_140_150(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compKineticEnergyForFG_0_10(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxx_xxxx = primBuffer.data(150 * idx);

            auto t_xxx_xxxy = primBuffer.data(150 * idx + 1);

            auto t_xxx_xxxz = primBuffer.data(150 * idx + 2);

            auto t_xxx_xxyy = primBuffer.data(150 * idx + 3);

            auto t_xxx_xxyz = primBuffer.data(150 * idx + 4);

            auto t_xxx_xxzz = primBuffer.data(150 * idx + 5);

            auto t_xxx_xyyy = primBuffer.data(150 * idx + 6);

            auto t_xxx_xyyz = primBuffer.data(150 * idx + 7);

            auto t_xxx_xyzz = primBuffer.data(150 * idx + 8);

            auto t_xxx_xzzz = primBuffer.data(150 * idx + 9);

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxx_xxxx, t_xxx_xxxy, t_xxx_xxxz, \
                                     t_xxx_xxyy, t_xxx_xxyz, t_xxx_xxzz, t_xxx_xyyy, t_xxx_xyyz, t_xxx_xyzz, t_xxx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxx_xxxx[j] = fl_s_0_0 * (5.625 * pa_x[j] * fl3_fx + 7.5 * fl3_fx * pb_x[j] + 0.75 * pa_xxx[j] * fl2_fx + 9.0 * pa_xx[j] * fl2_fx * pb_x[j] + 13.5 * pa_x[j] * fl2_fx * pb_xx[j] + 3.0 * fl2_fx * pb_xxx[j] + 3.0 * pa_xxx[j] * pb_xx[j] * fl1_fx + 6.0 * pa_xx[j] * fl1_fx * pb_xxx[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxxx[j] + pa_xxx[j] * pb_xxxx[j]);

                t_xxx_xxxx[j] += fl_r_0_0 * (-13.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 45.0 * pa_x[j] * fl3_fx * fl1_fz + 60.0 * fl3_fx * fl1_fz * pb_x[j] - 2.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 9.0 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 3.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 7.5 * pa_xxx[j] * fl1_fz * fl2_fx + 90.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 135.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xx[j] - 9.0 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 9.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 6.0 * fl1_fx * fl1_fz * fl1_fga * pb_xxx[j] - 6.0 * pa_xxx[j] * pb_xx[j] * fl1_fz * fl1_fgb + 30.0 * fl2_fx * fl1_fz * pb_xxx[j] + 36.0 * pa_xxx[j] * fl1_fz * pb_xx[j] * fl1_fx + 72.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxx[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxxx[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxx[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xxxx[j]);

                t_xxx_xxxy[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_y[j] + 2.25 * pa_xx[j] * fl2_fx * pb_y[j] + 6.75 * pa_x[j] * fl2_fx * pb_xy[j] + 2.25 * fl2_fx * pb_xxy[j] + 1.5 * pa_xxx[j] * pb_xy[j] * fl1_fx + 4.5 * pa_xx[j] * fl1_fx * pb_xxy[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxxy[j] + pa_xxx[j] * pb_xxxy[j]);

                t_xxx_xxxy[j] += fl_r_0_0 * (15.0 * fl3_fx * fl1_fz * pb_y[j] - 2.25 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 4.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] + 22.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] + 67.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xy[j] - 4.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxy[j] - 3.0 * pa_xxx[j] * pb_xy[j] * fl1_fz * fl1_fgb + 22.5 * fl2_fx * fl1_fz * pb_xxy[j] + 18.0 * pa_xxx[j] * fl1_fz * pb_xy[j] * fl1_fx + 54.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxy[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxxy[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxy[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xxxy[j]);

                t_xxx_xxxz[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_z[j] + 2.25 * pa_xx[j] * fl2_fx * pb_z[j] + 6.75 * pa_x[j] * fl2_fx * pb_xz[j] + 2.25 * fl2_fx * pb_xxz[j] + 1.5 * pa_xxx[j] * pb_xz[j] * fl1_fx + 4.5 * pa_xx[j] * fl1_fx * pb_xxz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxxz[j] + pa_xxx[j] * pb_xxxz[j]);

                t_xxx_xxxz[j] += fl_r_0_0 * (15.0 * fl3_fx * fl1_fz * pb_z[j] - 2.25 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 4.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 22.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] + 67.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xz[j] - 4.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxz[j] - 3.0 * pa_xxx[j] * pb_xz[j] * fl1_fz * fl1_fgb + 22.5 * fl2_fx * fl1_fz * pb_xxz[j] + 18.0 * pa_xxx[j] * fl1_fz * pb_xz[j] * fl1_fx + 54.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxxz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xxxz[j]);

                t_xxx_xxyy[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * fl3_fx * pb_x[j] + 0.25 * pa_xxx[j] * fl2_fx + 1.5 * pa_xx[j] * fl2_fx * pb_x[j] + 2.25 * pa_x[j] * fl2_fx * pb_yy[j] + 0.75 * pa_x[j] * fl2_fx * pb_xx[j] + 1.5 * fl2_fx * pb_xyy[j] + 0.5 * pa_xxx[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_yy[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xyy[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxyy[j] + pa_xxx[j] * pb_xxyy[j]);

                t_xxx_xxyy[j] += fl_r_0_0 * (-3.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * fl1_fz * pb_x[j] + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx + 15.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 22.5 * pa_x[j] * fl2_fx * fl1_fz * pb_yy[j] - 1.5 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yy[j] - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 1.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xyy[j] - pa_xxx[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xxx[j] * fl1_fz * fl1_fgb * pb_yy[j] + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xx[j] + 15.0 * fl2_fx * fl1_fz * pb_xyy[j] + 6.0 * pa_xxx[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xxx[j] * fl1_fz * fl1_fx * pb_yy[j] + 36.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xyy[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxyy[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxyy[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xxyy[j]);

                t_xxx_xxyz[j] = fl_s_0_0 * (2.25 * pa_x[j] * fl2_fx * pb_yz[j] + 1.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_xxx[j] * fl1_fx * pb_yz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxyz[j] + pa_xxx[j] * pb_xxyz[j]);

                t_xxx_xxyz[j] += fl_r_0_0 * (22.5 * pa_x[j] * fl2_fx * fl1_fz * pb_yz[j] - 1.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yz[j] - 1.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xyz[j] - pa_xxx[j] * fl1_fz * fl1_fgb * pb_yz[j] + 15.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_xxx[j] * fl1_fz * fl1_fx * pb_yz[j] + 36.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xyz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxyz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxyz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xxyz[j]);

                t_xxx_xxzz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * fl3_fx * pb_x[j] + 0.25 * pa_xxx[j] * fl2_fx + 1.5 * pa_xx[j] * fl2_fx * pb_x[j] + 2.25 * pa_x[j] * fl2_fx * pb_zz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xx[j] + 1.5 * fl2_fx * pb_xzz[j] + 0.5 * pa_xxx[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_zz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxzz[j] + pa_xxx[j] * pb_xxzz[j]);

                t_xxx_xxzz[j] += fl_r_0_0 * (-3.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * fl1_fz * pb_x[j] + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx + 15.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 22.5 * pa_x[j] * fl2_fx * fl1_fz * pb_zz[j] - 1.5 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 1.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xzz[j] - pa_xxx[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xxx[j] * fl1_fz * fl1_fgb * pb_zz[j] + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xx[j] + 15.0 * fl2_fx * fl1_fz * pb_xzz[j] + 6.0 * pa_xxx[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xxx[j] * fl1_fz * fl1_fx * pb_zz[j] + 36.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xzz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxzz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxzz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xxzz[j]);

                t_xxx_xyyy[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_y[j] + 2.25 * pa_xx[j] * fl2_fx * pb_y[j] + 2.25 * pa_x[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pb_yyy[j] + 1.5 * pa_xxx[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_x[j] * fl1_fx * pb_xyyy[j] + pa_xxx[j] * pb_xyyy[j]);

                t_xxx_xyyy[j] += fl_r_0_0 * (-2.25 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 4.5 * pa_xx[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 9.0 * fl3_fx * fl1_fz * pb_y[j] + 22.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] - 4.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyy[j] - 3.0 * pa_xxx[j] * pb_xy[j] * fl1_fz * fl1_fgb + 22.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pb_yyy[j] + 18.0 * pa_xxx[j] * fl1_fz * pb_xy[j] * fl1_fx + 18.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yyy[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xyyy[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyyy[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xyyy[j]);

                t_xxx_xyyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_z[j] + 0.75 * pa_xx[j] * fl2_fx * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pb_yyz[j] + 0.5 * pa_xxx[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xyyz[j] + pa_xxx[j] * pb_xyyz[j]);

                t_xxx_xyyz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 3.0 * fl3_fx * fl1_fz * pb_z[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] - 1.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyz[j] - pa_xxx[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pb_yyz[j] + 6.0 * pa_xxx[j] * fl1_fz * pb_xz[j] * fl1_fx + 18.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yyz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xyyz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyyz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xyyz[j]);

                t_xxx_xyzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_y[j] + 0.75 * pa_xx[j] * fl2_fx * pb_y[j] + 0.75 * pa_x[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pb_yzz[j] + 0.5 * pa_xxx[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xyzz[j] + pa_xxx[j] * pb_xyzz[j]);

                t_xxx_xyzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 1.5 * pa_xx[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] - 1.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yzz[j] - pa_xxx[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pb_yzz[j] + 6.0 * pa_xxx[j] * fl1_fz * pb_xy[j] * fl1_fx + 18.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yzz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xyzz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyzz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xyzz[j]);

                t_xxx_xzzz[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_z[j] + 2.25 * pa_xx[j] * fl2_fx * pb_z[j] + 2.25 * pa_x[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pb_zzz[j] + 1.5 * pa_xxx[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xzzz[j] + pa_xxx[j] * pb_xzzz[j]);

                t_xxx_xzzz[j] += fl_r_0_0 * (-2.25 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 4.5 * pa_xx[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 9.0 * fl3_fx * fl1_fz * pb_z[j] + 22.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] - 4.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_zzz[j] - 3.0 * pa_xxx[j] * pb_xz[j] * fl1_fz * fl1_fgb + 22.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pb_zzz[j] + 18.0 * pa_xxx[j] * fl1_fz * pb_xz[j] * fl1_fx + 18.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_zzz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xzzz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xzzz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_10_20(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxx_yyyy = primBuffer.data(150 * idx + 10);

            auto t_xxx_yyyz = primBuffer.data(150 * idx + 11);

            auto t_xxx_yyzz = primBuffer.data(150 * idx + 12);

            auto t_xxx_yzzz = primBuffer.data(150 * idx + 13);

            auto t_xxx_zzzz = primBuffer.data(150 * idx + 14);

            auto t_xxy_xxxx = primBuffer.data(150 * idx + 15);

            auto t_xxy_xxxy = primBuffer.data(150 * idx + 16);

            auto t_xxy_xxxz = primBuffer.data(150 * idx + 17);

            auto t_xxy_xxyy = primBuffer.data(150 * idx + 18);

            auto t_xxy_xxyz = primBuffer.data(150 * idx + 19);

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_z, pb_zz, pb_zzzz, \
                                     r_0_0, s_0_0, t_xxx_yyyy, t_xxx_yyyz, t_xxx_yyzz, t_xxx_yzzz, t_xxx_zzzz, \
                                     t_xxy_xxxx, t_xxy_xxxy, t_xxy_xxxz, t_xxy_xxyy, t_xxy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxx_yyyy[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pa_xxx[j] * fl2_fx + 4.5 * pa_x[j] * fl2_fx * pb_yy[j] + 3.0 * pa_xxx[j] * pb_yy[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pb_yyyy[j] + pa_xxx[j] * pb_yyyy[j]);

                t_xxx_yyyy[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl1_fz * fl3_fx + 7.5 * pa_xxx[j] * fl1_fz * fl2_fx - 9.0 * pa_x[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 9.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 6.0 * pa_xxx[j] * pb_yy[j] * fl1_fz * fl1_fgb + 45.0 * pa_x[j] * fl1_fz * fl2_fx * pb_yy[j] + 36.0 * pa_xxx[j] * fl1_fz * pb_yy[j] * fl1_fx - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yyyy[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyyy[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_yyyy[j]);

                t_xxx_yyyz[j] = fl_s_0_0 * (2.25 * pa_x[j] * fl2_fx * pb_yz[j] + 1.5 * pa_xxx[j] * pb_yz[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pb_yyyz[j] + pa_xxx[j] * pb_yyyz[j]);

                t_xxx_yyyz[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 3.0 * pa_xxx[j] * pb_yz[j] * fl1_fz * fl1_fgb + 22.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yz[j] + 18.0 * pa_xxx[j] * fl1_fz * pb_yz[j] * fl1_fx - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yyyz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyyz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_yyyz[j]);

                t_xxx_yyzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xxx[j] * fl2_fx + 0.75 * pa_x[j] * fl2_fx * pb_yy[j] + 0.75 * pa_x[j] * fl2_fx * pb_zz[j] + 0.5 * pa_xxx[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_zz[j] + 1.5 * pa_x[j] * fl1_fx * pb_yyzz[j] + pa_xxx[j] * pb_yyzz[j]);

                t_xxx_yyzz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx - 1.5 * pa_x[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 1.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - pa_xxx[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_xxx[j] * fl1_fz * fl1_fgb * pb_zz[j] + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yy[j] + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_zz[j] + 6.0 * pa_xxx[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_xxx[j] * fl1_fz * fl1_fx * pb_zz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yyzz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyzz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_yyzz[j]);

                t_xxx_yzzz[j] = fl_s_0_0 * (2.25 * pa_x[j] * fl2_fx * pb_yz[j] + 1.5 * pa_xxx[j] * pb_yz[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pb_yzzz[j] + pa_xxx[j] * pb_yzzz[j]);

                t_xxx_yzzz[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 3.0 * pa_xxx[j] * pb_yz[j] * fl1_fz * fl1_fgb + 22.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yz[j] + 18.0 * pa_xxx[j] * fl1_fz * pb_yz[j] * fl1_fx - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yzzz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yzzz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_yzzz[j]);

                t_xxx_zzzz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pa_xxx[j] * fl2_fx + 4.5 * pa_x[j] * fl2_fx * pb_zz[j] + 3.0 * pa_xxx[j] * pb_zz[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pb_zzzz[j] + pa_xxx[j] * pb_zzzz[j]);

                t_xxx_zzzz[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl1_fz * fl3_fx + 7.5 * pa_xxx[j] * fl1_fz * fl2_fx - 9.0 * pa_x[j] * fl1_fx * pb_zz[j] * fl1_fz * fl1_fgb - 9.0 * pa_x[j] * fl1_fz * fl1_fga * pb_zz[j] * fl1_fx - 6.0 * pa_xxx[j] * pb_zz[j] * fl1_fz * fl1_fgb + 45.0 * pa_x[j] * fl1_fz * fl2_fx * pb_zz[j] + 36.0 * pa_xxx[j] * fl1_fz * pb_zz[j] * fl1_fx - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_zzzz[j] + 18.0 * pa_x[j] * fl1_fz * fl1_fx * pb_zzzz[j] + 14.0 * pa_xxx[j] * fl1_fz * pb_zzzz[j]);

                t_xxy_xxxx[j] = fl_s_0_0 * (1.875 * fl3_fx * pa_y[j] + 0.75 * pa_xxy[j] * fl2_fx + 6.0 * pa_xy[j] * fl2_fx * pb_x[j] + 4.5 * fl2_fx * pa_y[j] * pb_xx[j] + 3.0 * pa_xxy[j] * pb_xx[j] * fl1_fx + 4.0 * pa_xy[j] * fl1_fx * pb_xxx[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxxx[j] + pa_xxy[j] * pb_xxxx[j]);

                t_xxy_xxxx[j] += fl_r_0_0 * (-4.5 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + 15.0 * fl3_fx * fl1_fz * pa_y[j] - 0.75 * fl1_fz * fl1_fga * pa_y[j] * fl2_fx - 3.0 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 7.5 * pa_xxy[j] * fl1_fz * fl2_fx + 60.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_x[j] + 45.0 * fl2_fx * fl1_fz * pa_y[j] * pb_xx[j] - 3.0 * fl1_fx * pa_y[j] * pb_xx[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_y[j] * pb_xx[j] * fl1_fx - 6.0 * pa_xxy[j] * pb_xx[j] * fl1_fz * fl1_fgb + 36.0 * pa_xxy[j] * fl1_fz * pb_xx[j] * fl1_fx + 48.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xxx[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxxx[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxxx[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xxxx[j]);

                t_xxy_xxxy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 1.125 * fl3_fx * pb_x[j] + 0.75 * pa_xx[j] * fl2_fx * pb_x[j] + 1.5 * pa_xy[j] * fl2_fx * pb_y[j] + 1.5 * pa_x[j] * fl2_fx * pb_xx[j] + 2.25 * fl2_fx * pa_y[j] * pb_xy[j] + 0.25 * fl2_fx * pb_xxx[j] + 1.5 * pa_xxy[j] * pb_xy[j] * fl1_fx + 0.5 * pa_xx[j] * fl1_fx * pb_xxx[j] + 3.0 * pa_xy[j] * fl1_fx * pb_xxy[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxxy[j] + pa_xxy[j] * pb_xxxy[j]);

                t_xxy_xxxy[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz + 9.0 * fl3_fx * fl1_fz * pb_x[j] - 0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 1.5 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 3.0 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 15.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_y[j] + 15.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xx[j] + 22.5 * fl2_fx * fl1_fz * pa_y[j] * pb_xy[j] - 1.5 * fl1_fx * pa_y[j] * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_xy[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xxx[j] - 3.0 * pa_xxy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_xxx[j] + 18.0 * pa_xxy[j] * fl1_fz * pb_xy[j] * fl1_fx + 6.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxx[j] + 36.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xxy[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxxy[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxxy[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xxxy[j]);

                t_xxy_xxxz[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl2_fx * pb_z[j] + 2.25 * fl2_fx * pa_y[j] * pb_xz[j] + 1.5 * pa_xxy[j] * pb_xz[j] * fl1_fx + 3.0 * pa_xy[j] * fl1_fx * pb_xxz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxxz[j] + pa_xxy[j] * pb_xxxz[j]);

                t_xxy_xxxz[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 15.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_z[j] + 22.5 * fl2_fx * fl1_fz * pa_y[j] * pb_xz[j] - 1.5 * fl1_fx * pa_y[j] * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_xz[j] * fl1_fx - 3.0 * pa_xxy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa_xxy[j] * fl1_fz * pb_xz[j] * fl1_fx + 36.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xxz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxxz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxxz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xxxz[j]);

                t_xxy_xxyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_y[j] + 0.75 * fl3_fx * pb_y[j] + 0.25 * pa_xxy[j] * fl2_fx + 0.5 * pa_xx[j] * fl2_fx * pb_y[j] + pa_xy[j] * fl2_fx * pb_x[j] + 2.0 * pa_x[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_y[j] * pb_yy[j] + 0.25 * fl2_fx * pa_y[j] * pb_xx[j] + 0.5 * fl2_fx * pb_xxy[j] + 0.5 * pa_xxy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxy[j] * fl1_fx * pb_yy[j] + pa_xx[j] * fl1_fx * pb_xxy[j] + 2.0 * pa_xy[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxyy[j] + pa_xxy[j] * pb_xxyy[j]);

                t_xxy_xxyy[j] += fl_r_0_0 * (-fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_y[j] + 6.0 * fl3_fx * fl1_fz * pb_y[j] - 0.5 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.25 * fl1_fz * fl1_fga * pa_y[j] * fl2_fx - 0.5 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 2.0 * pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_xxy[j] * fl1_fz * fl2_fx + 5.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] + 10.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_x[j] + 20.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pa_y[j] * pb_yy[j] - 0.5 * fl1_fx * pa_y[j] * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_y[j] * fl1_fz * fl1_fgb * pb_yy[j] - 0.5 * fl1_fz * fl1_fga * pa_y[j] * pb_xx[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_y[j] * fl1_fx * pb_yy[j] - fl1_fz * fl1_fga * fl1_fx * pb_xxy[j] - pa_xxy[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xxy[j] * fl1_fz * fl1_fgb * pb_yy[j] + 2.5 * fl2_fx * fl1_fz * pa_y[j] * pb_xx[j] + 5.0 * fl2_fx * fl1_fz * pb_xxy[j] + 6.0 * pa_xxy[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xxy[j] * fl1_fz * fl1_fx * pb_yy[j] + 12.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxy[j] + 24.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xyy[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxyy[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxyy[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xxyy[j]);

                t_xxy_xxyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_z[j] + 0.25 * pa_xx[j] * fl2_fx * pb_z[j] + pa_x[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_y[j] * pb_yz[j] + 0.25 * fl2_fx * pb_xxz[j] + 0.5 * pa_xxy[j] * fl1_fx * pb_yz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxz[j] + 2.0 * pa_xy[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxyz[j] + pa_xxy[j] * pb_xxyz[j]);

                t_xxy_xxyz[j] += fl_r_0_0 * (3.0 * fl3_fx * fl1_fz * pb_z[j] - 0.25 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.25 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 2.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] + 10.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pa_y[j] * pb_yz[j] - 0.5 * fl1_fx * pa_y[j] * fl1_fz * fl1_fgb * pb_yz[j] - 0.5 * fl1_fz * fl1_fga * pa_y[j] * fl1_fx * pb_yz[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xxz[j] - pa_xxy[j] * fl1_fz * fl1_fgb * pb_yz[j] + 2.5 * fl2_fx * fl1_fz * pb_xxz[j] + 6.0 * pa_xxy[j] * fl1_fz * fl1_fx * pb_yz[j] + 6.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxz[j] + 24.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xyz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxyz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxyz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_20_30(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxy = paDistances.data(19 * idx + 10);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxy_xxzz = primBuffer.data(150 * idx + 20);

            auto t_xxy_xyyy = primBuffer.data(150 * idx + 21);

            auto t_xxy_xyyz = primBuffer.data(150 * idx + 22);

            auto t_xxy_xyzz = primBuffer.data(150 * idx + 23);

            auto t_xxy_xzzz = primBuffer.data(150 * idx + 24);

            auto t_xxy_yyyy = primBuffer.data(150 * idx + 25);

            auto t_xxy_yyyz = primBuffer.data(150 * idx + 26);

            auto t_xxy_yyzz = primBuffer.data(150 * idx + 27);

            auto t_xxy_yzzz = primBuffer.data(150 * idx + 28);

            auto t_xxy_zzzz = primBuffer.data(150 * idx + 29);

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     r_0_0, s_0_0, t_xxy_xxzz, t_xxy_xyyy, t_xxy_xyyz, t_xxy_xyzz, t_xxy_xzzz, \
                                     t_xxy_yyyy, t_xxy_yyyz, t_xxy_yyzz, t_xxy_yzzz, t_xxy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxy_xxzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_y[j] + 0.25 * pa_xxy[j] * fl2_fx + pa_xy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_y[j] * pb_zz[j] + 0.25 * fl2_fx * pa_y[j] * pb_xx[j] + 0.5 * pa_xxy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxy[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xy[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxzz[j] + pa_xxy[j] * pb_xxzz[j]);

                t_xxy_xxzz[j] += fl_r_0_0 * (-fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_y[j] - 0.25 * fl1_fz * fl1_fga * pa_y[j] * fl2_fx - pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_xxy[j] * fl1_fz * fl2_fx + 10.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_x[j] + 7.5 * fl2_fx * fl1_fz * pa_y[j] * pb_zz[j] - 0.5 * fl1_fx * pa_y[j] * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_y[j] * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * fl1_fz * fl1_fga * pa_y[j] * pb_xx[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_y[j] * fl1_fx * pb_zz[j] - pa_xxy[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xxy[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * fl2_fx * fl1_fz * pa_y[j] * pb_xx[j] + 6.0 * pa_xxy[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xxy[j] * fl1_fz * fl1_fx * pb_zz[j] + 24.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xzz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxzz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxzz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xxzz[j]);

                t_xxy_xyyy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 0.375 * fl3_fx * pb_x[j] + 0.75 * pa_xx[j] * fl2_fx * pb_x[j] + 1.5 * pa_xy[j] * fl2_fx * pb_y[j] + 1.5 * pa_x[j] * fl2_fx * pb_yy[j] + 0.75 * fl2_fx * pa_y[j] * pb_xy[j] + 0.75 * fl2_fx * pb_xyy[j] + 1.5 * pa_xxy[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_xyy[j] + pa_xy[j] * fl1_fx * pb_yyy[j] + 0.5 * fl1_fx * pa_y[j] * pb_xyyy[j] + pa_xxy[j] * pb_xyyy[j]);

                t_xxy_xyyy[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 1.5 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 3.0 * pa_xy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 15.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_y[j] + 15.0 * pa_x[j] * fl2_fx * fl1_fz * pb_yy[j] - 1.5 * fl1_fx * pa_y[j] * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_xy[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_xyy[j] - 3.0 * pa_xxy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_y[j] * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pb_xyy[j] + 18.0 * pa_xxy[j] * fl1_fz * pb_xy[j] * fl1_fx + 18.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_yyy[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xyyy[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xyyy[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xyyy[j]);

                t_xxy_xyyz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx * pb_z[j] + pa_x[j] * fl2_fx * pb_yz[j] + 0.25 * fl2_fx * pa_y[j] * pb_xz[j] + 0.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_xxy[j] * pb_xz[j] * fl1_fx + pa_xx[j] * fl1_fx * pb_xyz[j] + pa_xy[j] * fl1_fx * pb_yyz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xyyz[j] + pa_xxy[j] * pb_xyyz[j]);

                t_xxy_xyyz[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 5.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_z[j] + 10.0 * pa_x[j] * fl2_fx * fl1_fz * pb_yz[j] - 0.5 * fl1_fx * pa_y[j] * pb_xz[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_y[j] * pb_xz[j] * fl1_fx - fl1_fz * fl1_fga * fl1_fx * pb_xyz[j] - pa_xxy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pa_y[j] * pb_xz[j] + 5.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_xxy[j] * fl1_fz * pb_xz[j] * fl1_fx + 12.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_yyz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xyyz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xyyz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xyyz[j]);

                t_xxy_xyzz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl3_fx + 0.125 * fl3_fx * pb_x[j] + 0.25 * pa_xx[j] * fl2_fx * pb_x[j] + 0.5 * pa_xy[j] * fl2_fx * pb_y[j] + 0.5 * pa_x[j] * fl2_fx * pb_zz[j] + 0.25 * fl2_fx * pa_y[j] * pb_xy[j] + 0.25 * fl2_fx * pb_xzz[j] + 0.5 * pa_xxy[j] * pb_xy[j] * fl1_fx + 0.5 * pa_xx[j] * fl1_fx * pb_xzz[j] + pa_xy[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xyzz[j] + pa_xxy[j] * pb_xyzz[j]);

                t_xxy_xyzz[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 0.5 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - pa_xy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + fl3_fx * fl1_fz * pb_x[j] + 2.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 5.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_y[j] + 5.0 * pa_x[j] * fl2_fx * fl1_fz * pb_zz[j] - 0.5 * fl1_fx * pa_y[j] * pb_xy[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_y[j] * pb_xy[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xzz[j] - pa_xxy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pa_y[j] * pb_xy[j] + 2.5 * fl2_fx * fl1_fz * pb_xzz[j] + 6.0 * pa_xxy[j] * fl1_fz * pb_xy[j] * fl1_fx + 6.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xzz[j] + 12.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_yzz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xyzz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xyzz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xyzz[j]);

                t_xxy_xzzz[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl2_fx * pb_z[j] + 0.75 * fl2_fx * pa_y[j] * pb_xz[j] + 1.5 * pa_xxy[j] * pb_xz[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xzzz[j] + pa_xxy[j] * pb_xzzz[j]);

                t_xxy_xzzz[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 15.0 * pa_xy[j] * fl2_fx * fl1_fz * pb_z[j] - 1.5 * fl1_fx * pa_y[j] * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_xz[j] * fl1_fx - 3.0 * pa_xxy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_y[j] * pb_xz[j] + 18.0 * pa_xxy[j] * fl1_fz * pb_xz[j] * fl1_fx + 12.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_zzz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xzzz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_xzzz[j]);

                t_xxy_yyyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_y[j] + 1.5 * fl3_fx * pb_y[j] + 0.75 * pa_xxy[j] * fl2_fx + 3.0 * pa_xx[j] * fl2_fx * pb_y[j] + 1.5 * fl2_fx * pa_y[j] * pb_yy[j] + fl2_fx * pb_yyy[j] + 3.0 * pa_xxy[j] * pb_yy[j] * fl1_fx + 2.0 * pa_xx[j] * fl1_fx * pb_yyy[j] + 0.5 * fl1_fx * pa_y[j] * pb_yyyy[j] + pa_xxy[j] * pb_yyyy[j]);

                t_xxy_yyyy[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb - 3.0 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * pa_y[j] * fl2_fx - 3.0 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - 3.0 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa_xx[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_y[j] + 12.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * pa_xxy[j] * fl1_fz * fl2_fx + 30.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] - 3.0 * fl1_fx * pa_y[j] * pb_yy[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_y[j] * pb_yy[j] * fl1_fx - 2.0 * fl1_fz * fl1_fga * fl1_fx * pb_yyy[j] - 6.0 * pa_xxy[j] * pb_yy[j] * fl1_fz * fl1_fgb + 15.0 * fl2_fx * fl1_fz * pa_y[j] * pb_yy[j] + 10.0 * fl2_fx * fl1_fz * pb_yyy[j] + 36.0 * pa_xxy[j] * fl1_fz * pb_yy[j] * fl1_fx + 24.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yyy[j] - fl1_fz * fl1_fga * pa_y[j] * pb_yyyy[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_yyyy[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_yyyy[j]);

                t_xxy_yyyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_z[j] + 0.75 * pa_xx[j] * fl2_fx * pb_z[j] + 0.75 * fl2_fx * pa_y[j] * pb_yz[j] + 0.75 * fl2_fx * pb_yyz[j] + 1.5 * pa_xxy[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_yyz[j] + 0.5 * fl1_fx * pa_y[j] * pb_yyyz[j] + pa_xxy[j] * pb_yyyz[j]);

                t_xxy_yyyz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 3.0 * fl3_fx * fl1_fz * pb_z[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] - 1.5 * fl1_fx * pa_y[j] * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_yz[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_yyz[j] - 3.0 * pa_xxy[j] * pb_yz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_y[j] * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pb_yyz[j] + 18.0 * pa_xxy[j] * fl1_fz * pb_yz[j] * fl1_fx + 18.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yyz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_yyyz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_yyyz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_yyyz[j]);

                t_xxy_yyzz[j] = fl_s_0_0 * (0.125 * fl3_fx * pa_y[j] + 0.25 * fl3_fx * pb_y[j] + 0.25 * pa_xxy[j] * fl2_fx + 0.5 * pa_xx[j] * fl2_fx * pb_y[j] + 0.25 * fl2_fx * pa_y[j] * pb_yy[j] + 0.25 * fl2_fx * pa_y[j] * pb_zz[j] + 0.5 * fl2_fx * pb_yzz[j] + 0.5 * pa_xxy[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxy[j] * fl1_fx * pb_zz[j] + pa_xx[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_y[j] * pb_yyzz[j] + pa_xxy[j] * pb_yyzz[j]);

                t_xxy_yyzz[j] += fl_r_0_0 * (-0.5 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb - 0.5 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * pa_y[j] * fl2_fx - 0.5 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - pa_xx[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + fl3_fx * fl1_fz * pa_y[j] + 2.0 * fl3_fx * fl1_fz * pb_y[j] + 2.5 * pa_xxy[j] * fl1_fz * fl2_fx + 5.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] - 0.5 * fl1_fx * pa_y[j] * pb_yy[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_y[j] * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * fl1_fz * fl1_fga * pa_y[j] * pb_yy[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_y[j] * fl1_fx * pb_zz[j] - fl1_fz * fl1_fga * fl1_fx * pb_yzz[j] - pa_xxy[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_xxy[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * fl2_fx * fl1_fz * pa_y[j] * pb_yy[j] + 2.5 * fl2_fx * fl1_fz * pa_y[j] * pb_zz[j] + 5.0 * fl2_fx * fl1_fz * pb_yzz[j] + 6.0 * pa_xxy[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_xxy[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yzz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_yyzz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_yyzz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_yyzz[j]);

                t_xxy_yzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_z[j] + 0.75 * pa_xx[j] * fl2_fx * pb_z[j] + 0.75 * fl2_fx * pa_y[j] * pb_yz[j] + 0.25 * fl2_fx * pb_zzz[j] + 1.5 * pa_xxy[j] * pb_yz[j] * fl1_fx + 0.5 * pa_xx[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_y[j] * pb_yzzz[j] + pa_xxy[j] * pb_yzzz[j]);

                t_xxy_yzzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - 1.5 * pa_xx[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_z[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] - 1.5 * fl1_fx * pa_y[j] * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_yz[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_zzz[j] - 3.0 * pa_xxy[j] * pb_yz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_y[j] * pb_yz[j] + 2.5 * fl2_fx * fl1_fz * pb_zzz[j] + 18.0 * pa_xxy[j] * fl1_fz * pb_yz[j] * fl1_fx + 6.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_zzz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_yzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_yzzz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_yzzz[j]);

                t_xxy_zzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_y[j] + 0.75 * pa_xxy[j] * fl2_fx + 1.5 * fl2_fx * pa_y[j] * pb_zz[j] + 3.0 * pa_xxy[j] * pb_zz[j] * fl1_fx + 0.5 * fl1_fx * pa_y[j] * pb_zzzz[j] + pa_xxy[j] * pb_zzzz[j]);

                t_xxy_zzzz[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * pa_y[j] * fl2_fx - 3.0 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_y[j] + 7.5 * pa_xxy[j] * fl1_fz * fl2_fx - 3.0 * fl1_fx * pa_y[j] * pb_zz[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_y[j] * pb_zz[j] * fl1_fx - 6.0 * pa_xxy[j] * pb_zz[j] * fl1_fz * fl1_fgb + 15.0 * fl2_fx * fl1_fz * pa_y[j] * pb_zz[j] + 36.0 * pa_xxy[j] * fl1_fz * pb_zz[j] * fl1_fx - fl1_fz * fl1_fga * pa_y[j] * pb_zzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_y[j] * pb_zzzz[j] + 14.0 * pa_xxy[j] * fl1_fz * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_30_40(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xz = paDistances.data(19 * idx + 5);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxz = paDistances.data(19 * idx + 11);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxz_xxxx = primBuffer.data(150 * idx + 30);

            auto t_xxz_xxxy = primBuffer.data(150 * idx + 31);

            auto t_xxz_xxxz = primBuffer.data(150 * idx + 32);

            auto t_xxz_xxyy = primBuffer.data(150 * idx + 33);

            auto t_xxz_xxyz = primBuffer.data(150 * idx + 34);

            auto t_xxz_xxzz = primBuffer.data(150 * idx + 35);

            auto t_xxz_xyyy = primBuffer.data(150 * idx + 36);

            auto t_xxz_xyyz = primBuffer.data(150 * idx + 37);

            auto t_xxz_xyzz = primBuffer.data(150 * idx + 38);

            auto t_xxz_xzzz = primBuffer.data(150 * idx + 39);

            // Batch of Integrals (30,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxz_xxxx, t_xxz_xxxy, t_xxz_xxxz, \
                                     t_xxz_xxyy, t_xxz_xxyz, t_xxz_xxzz, t_xxz_xyyy, t_xxz_xyyz, t_xxz_xyzz, t_xxz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxz_xxxx[j] = fl_s_0_0 * (1.875 * fl3_fx * pa_z[j] + 0.75 * pa_xxz[j] * fl2_fx + 6.0 * pa_xz[j] * fl2_fx * pb_x[j] + 4.5 * fl2_fx * pa_z[j] * pb_xx[j] + 3.0 * pa_xxz[j] * pb_xx[j] * fl1_fx + 4.0 * pa_xz[j] * fl1_fx * pb_xxx[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxxx[j] + pa_xxz[j] * pb_xxxx[j]);

                t_xxz_xxxx[j] += fl_r_0_0 * (-4.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 15.0 * fl3_fx * fl1_fz * pa_z[j] - 0.75 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 3.0 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 7.5 * pa_xxz[j] * fl1_fz * fl2_fx + 60.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_x[j] + 45.0 * fl2_fx * fl1_fz * pa_z[j] * pb_xx[j] - 3.0 * fl1_fx * pa_z[j] * pb_xx[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_z[j] * pb_xx[j] * fl1_fx - 6.0 * pa_xxz[j] * pb_xx[j] * fl1_fz * fl1_fgb + 36.0 * pa_xxz[j] * fl1_fz * pb_xx[j] * fl1_fx + 48.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xxx[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxxx[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxxx[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xxxx[j]);

                t_xxz_xxxy[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl2_fx * pb_y[j] + 2.25 * fl2_fx * pa_z[j] * pb_xy[j] + 1.5 * pa_xxz[j] * pb_xy[j] * fl1_fx + 3.0 * pa_xz[j] * fl1_fx * pb_xxy[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxxy[j] + pa_xxz[j] * pb_xxxy[j]);

                t_xxz_xxxy[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] + 15.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_y[j] + 22.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xy[j] - 1.5 * fl1_fx * pa_z[j] * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xy[j] * fl1_fx - 3.0 * pa_xxz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa_xxz[j] * fl1_fz * pb_xy[j] * fl1_fx + 36.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xxy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxxy[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxxy[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xxxy[j]);

                t_xxz_xxxz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 1.125 * fl3_fx * pb_x[j] + 0.75 * pa_xx[j] * fl2_fx * pb_x[j] + 1.5 * pa_xz[j] * fl2_fx * pb_z[j] + 1.5 * pa_x[j] * fl2_fx * pb_xx[j] + 2.25 * fl2_fx * pa_z[j] * pb_xz[j] + 0.25 * fl2_fx * pb_xxx[j] + 1.5 * pa_xxz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xx[j] * fl1_fx * pb_xxx[j] + 3.0 * pa_xz[j] * fl1_fx * pb_xxz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxxz[j] + pa_xxz[j] * pb_xxxz[j]);

                t_xxz_xxxz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz + 9.0 * fl3_fx * fl1_fz * pb_x[j] - 0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 1.5 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 3.0 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 15.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_z[j] + 15.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xx[j] + 22.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xz[j] - 1.5 * fl1_fx * pa_z[j] * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xz[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xxx[j] - 3.0 * pa_xxz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_xxx[j] + 18.0 * pa_xxz[j] * fl1_fz * pb_xz[j] * fl1_fx + 6.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxx[j] + 36.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xxz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxxz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxxz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xxxz[j]);

                t_xxz_xxyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 0.25 * pa_xxz[j] * fl2_fx + pa_xz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_z[j] * pb_yy[j] + 0.25 * fl2_fx * pa_z[j] * pb_xx[j] + 0.5 * pa_xxz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxz[j] * fl1_fx * pb_yy[j] + 2.0 * pa_xz[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxyy[j] + pa_xxz[j] * pb_xxyy[j]);

                t_xxz_xxyy[j] += fl_r_0_0 * (-fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] - 0.25 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_xxz[j] * fl1_fz * fl2_fx + 10.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_x[j] + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yy[j] - 0.5 * fl1_fx * pa_z[j] * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_yy[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xx[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_yy[j] - pa_xxz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xxz[j] * fl1_fz * fl1_fgb * pb_yy[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xx[j] + 6.0 * pa_xxz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xxz[j] * fl1_fz * fl1_fx * pb_yy[j] + 24.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xyy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxyy[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxyy[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xxyy[j]);

                t_xxz_xxyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_y[j] + 0.25 * pa_xx[j] * fl2_fx * pb_y[j] + pa_x[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_z[j] * pb_yz[j] + 0.25 * fl2_fx * pb_xxy[j] + 0.5 * pa_xxz[j] * fl1_fx * pb_yz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxy[j] + 2.0 * pa_xz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxyz[j] + pa_xxz[j] * pb_xxyz[j]);

                t_xxz_xxyz[j] += fl_r_0_0 * (3.0 * fl3_fx * fl1_fz * pb_y[j] - 0.25 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.25 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] + 2.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] + 10.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yz[j] - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_yz[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_yz[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xxy[j] - pa_xxz[j] * fl1_fz * fl1_fgb * pb_yz[j] + 2.5 * fl2_fx * fl1_fz * pb_xxy[j] + 6.0 * pa_xxz[j] * fl1_fz * fl1_fx * pb_yz[j] + 6.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxy[j] + 24.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xyz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxyz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxyz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xxyz[j]);

                t_xxz_xxzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 0.75 * fl3_fx * pb_z[j] + 0.25 * pa_xxz[j] * fl2_fx + 0.5 * pa_xx[j] * fl2_fx * pb_z[j] + pa_xz[j] * fl2_fx * pb_x[j] + 2.0 * pa_x[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_z[j] * pb_zz[j] + 0.25 * fl2_fx * pa_z[j] * pb_xx[j] + 0.5 * fl2_fx * pb_xxz[j] + 0.5 * pa_xxz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxz[j] * fl1_fx * pb_zz[j] + pa_xx[j] * fl1_fx * pb_xxz[j] + 2.0 * pa_xz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxzz[j] + pa_xxz[j] * pb_xxzz[j]);

                t_xxz_xxzz[j] += fl_r_0_0 * (-fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] + 6.0 * fl3_fx * fl1_fz * pb_z[j] - 0.5 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.25 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 0.5 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 2.0 * pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_xxz[j] * fl1_fz * fl2_fx + 5.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] + 10.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_x[j] + 20.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_zz[j] - 0.5 * fl1_fx * pa_z[j] * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xx[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_zz[j] - fl1_fz * fl1_fga * fl1_fx * pb_xxz[j] - pa_xxz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xxz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xx[j] + 5.0 * fl2_fx * fl1_fz * pb_xxz[j] + 6.0 * pa_xxz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xxz[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xxz[j] + 24.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxzz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xxzz[j]);

                t_xxz_xyyy[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pa_z[j] * pb_xy[j] + 1.5 * pa_xxz[j] * pb_xy[j] * fl1_fx + pa_xz[j] * fl1_fx * pb_yyy[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyyy[j] + pa_xxz[j] * pb_xyyy[j]);

                t_xxz_xyyy[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 15.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_y[j] - 1.5 * fl1_fx * pa_z[j] * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xy[j] * fl1_fx - 3.0 * pa_xxz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xy[j] + 18.0 * pa_xxz[j] * fl1_fz * pb_xy[j] * fl1_fx + 12.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yyy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyyy[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyyy[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xyyy[j]);

                t_xxz_xyyz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl3_fx + 0.125 * fl3_fx * pb_x[j] + 0.25 * pa_xx[j] * fl2_fx * pb_x[j] + 0.5 * pa_xz[j] * fl2_fx * pb_z[j] + 0.5 * pa_x[j] * fl2_fx * pb_yy[j] + 0.25 * fl2_fx * pa_z[j] * pb_xz[j] + 0.25 * fl2_fx * pb_xyy[j] + 0.5 * pa_xxz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xx[j] * fl1_fx * pb_xyy[j] + pa_xz[j] * fl1_fx * pb_yyz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyyz[j] + pa_xxz[j] * pb_xyyz[j]);

                t_xxz_xyyz[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 0.5 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + fl3_fx * fl1_fz * pb_x[j] + 2.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 5.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_z[j] + 5.0 * pa_x[j] * fl2_fx * fl1_fz * pb_yy[j] - 0.5 * fl1_fx * pa_z[j] * pb_xz[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xz[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xyy[j] - pa_xxz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xz[j] + 2.5 * fl2_fx * fl1_fz * pb_xyy[j] + 6.0 * pa_xxz[j] * fl1_fz * pb_xz[j] * fl1_fx + 6.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yyz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyyz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyyz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xyyz[j]);

                t_xxz_xyzz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx * pb_y[j] + pa_x[j] * fl2_fx * pb_yz[j] + 0.25 * fl2_fx * pa_z[j] * pb_xy[j] + 0.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_xxz[j] * pb_xy[j] * fl1_fx + pa_xx[j] * fl1_fx * pb_xyz[j] + pa_xz[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyzz[j] + pa_xxz[j] * pb_xyzz[j]);

                t_xxz_xyzz[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 5.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_y[j] + 10.0 * pa_x[j] * fl2_fx * fl1_fz * pb_yz[j] - 0.5 * fl1_fx * pa_z[j] * pb_xy[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xy[j] * fl1_fx - fl1_fz * fl1_fga * fl1_fx * pb_xyz[j] - pa_xxz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xy[j] + 5.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_xxz[j] * fl1_fz * pb_xy[j] * fl1_fx + 12.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyzz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xyzz[j]);

                t_xxz_xzzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 0.375 * fl3_fx * pb_x[j] + 0.75 * pa_xx[j] * fl2_fx * pb_x[j] + 1.5 * pa_xz[j] * fl2_fx * pb_z[j] + 1.5 * pa_x[j] * fl2_fx * pb_zz[j] + 0.75 * fl2_fx * pa_z[j] * pb_xz[j] + 0.75 * fl2_fx * pb_xzz[j] + 1.5 * pa_xxz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_xzz[j] + pa_xz[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xzzz[j] + pa_xxz[j] * pb_xzzz[j]);

                t_xxz_xzzz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 1.5 * pa_xx[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 3.0 * pa_xz[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_x[j] + 15.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_z[j] + 15.0 * pa_x[j] * fl2_fx * fl1_fz * pb_zz[j] - 1.5 * fl1_fx * pa_z[j] * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xz[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_xzz[j] - 3.0 * pa_xxz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pb_xzz[j] + 18.0 * pa_xxz[j] * fl1_fz * pb_xz[j] * fl1_fx + 18.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xzz[j] + 12.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_zzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xzzz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_40_50(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_yy = paDistances.data(19 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxz_yyyy = primBuffer.data(150 * idx + 40);

            auto t_xxz_yyyz = primBuffer.data(150 * idx + 41);

            auto t_xxz_yyzz = primBuffer.data(150 * idx + 42);

            auto t_xxz_yzzz = primBuffer.data(150 * idx + 43);

            auto t_xxz_zzzz = primBuffer.data(150 * idx + 44);

            auto t_xyy_xxxx = primBuffer.data(150 * idx + 45);

            auto t_xyy_xxxy = primBuffer.data(150 * idx + 46);

            auto t_xyy_xxxz = primBuffer.data(150 * idx + 47);

            auto t_xyy_xxyy = primBuffer.data(150 * idx + 48);

            auto t_xyy_xxyz = primBuffer.data(150 * idx + 49);

            // Batch of Integrals (40,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xy, pa_xyy, pa_y, pa_yy, pa_z, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xxz_yyyy, t_xxz_yyyz, \
                                     t_xxz_yyzz, t_xxz_yzzz, t_xxz_zzzz, t_xyy_xxxx, t_xyy_xxxy, t_xyy_xxxz, t_xyy_xxyy, \
                                     t_xyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xxz_yyyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 0.75 * pa_xxz[j] * fl2_fx + 1.5 * fl2_fx * pa_z[j] * pb_yy[j] + 3.0 * pa_xxz[j] * pb_yy[j] * fl1_fx + 0.5 * fl1_fx * pa_z[j] * pb_yyyy[j] + pa_xxz[j] * pb_yyyy[j]);

                t_xxz_yyyy[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 3.0 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] + 7.5 * pa_xxz[j] * fl1_fz * fl2_fx - 3.0 * fl1_fx * pa_z[j] * pb_yy[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_z[j] * pb_yy[j] * fl1_fx - 6.0 * pa_xxz[j] * pb_yy[j] * fl1_fz * fl1_fgb + 15.0 * fl2_fx * fl1_fz * pa_z[j] * pb_yy[j] + 36.0 * pa_xxz[j] * fl1_fz * pb_yy[j] * fl1_fx - fl1_fz * fl1_fga * pa_z[j] * pb_yyyy[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyyy[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_yyyy[j]);

                t_xxz_yyyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_y[j] + 0.75 * pa_xx[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pa_z[j] * pb_yz[j] + 0.25 * fl2_fx * pb_yyy[j] + 1.5 * pa_xxz[j] * pb_yz[j] * fl1_fx + 0.5 * pa_xx[j] * fl1_fx * pb_yyy[j] + 0.5 * fl1_fx * pa_z[j] * pb_yyyz[j] + pa_xxz[j] * pb_yyyz[j]);

                t_xxz_yyyz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - 1.5 * pa_xx[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] - 1.5 * fl1_fx * pa_z[j] * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_yz[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_yyy[j] - 3.0 * pa_xxz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yz[j] + 2.5 * fl2_fx * fl1_fz * pb_yyy[j] + 18.0 * pa_xxz[j] * fl1_fz * pb_yz[j] * fl1_fx + 6.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yyy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yyyz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyyz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_yyyz[j]);

                t_xxz_yyzz[j] = fl_s_0_0 * (0.125 * fl3_fx * pa_z[j] + 0.25 * fl3_fx * pb_z[j] + 0.25 * pa_xxz[j] * fl2_fx + 0.5 * pa_xx[j] * fl2_fx * pb_z[j] + 0.25 * fl2_fx * pa_z[j] * pb_yy[j] + 0.25 * fl2_fx * pa_z[j] * pb_zz[j] + 0.5 * fl2_fx * pb_yyz[j] + 0.5 * pa_xxz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxz[j] * fl1_fx * pb_zz[j] + pa_xx[j] * fl1_fx * pb_yyz[j] + 0.5 * fl1_fx * pa_z[j] * pb_yyzz[j] + pa_xxz[j] * pb_yyzz[j]);

                t_xxz_yyzz[j] += fl_r_0_0 * (-0.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb - 0.5 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.25 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 0.5 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + fl3_fx * fl1_fz * pa_z[j] + 2.0 * fl3_fx * fl1_fz * pb_z[j] + 2.5 * pa_xxz[j] * fl1_fz * fl2_fx + 5.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] - 0.5 * fl1_fx * pa_z[j] * pb_yy[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_yy[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_zz[j] - fl1_fz * fl1_fga * fl1_fx * pb_yyz[j] - pa_xxz[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_xxz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yy[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_zz[j] + 5.0 * fl2_fx * fl1_fz * pb_yyz[j] + 6.0 * pa_xxz[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_xxz[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yyz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yyzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyzz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_yyzz[j]);

                t_xxz_yzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_y[j] + 0.75 * pa_xx[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pa_z[j] * pb_yz[j] + 0.75 * fl2_fx * pb_yzz[j] + 1.5 * pa_xxz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_yzzz[j] + pa_xxz[j] * pb_yzzz[j]);

                t_xxz_yzzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - 1.5 * pa_xx[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * pa_xx[j] * fl1_fz * fl2_fx * pb_y[j] - 1.5 * fl1_fx * pa_z[j] * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_yz[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_yzz[j] - 3.0 * pa_xxz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pb_yzz[j] + 18.0 * pa_xxz[j] * fl1_fz * pb_yz[j] * fl1_fx + 18.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yzzz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_yzzz[j]);

                t_xxz_zzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 1.5 * fl3_fx * pb_z[j] + 0.75 * pa_xxz[j] * fl2_fx + 3.0 * pa_xx[j] * fl2_fx * pb_z[j] + 1.5 * fl2_fx * pa_z[j] * pb_zz[j] + fl2_fx * pb_zzz[j] + 3.0 * pa_xxz[j] * pb_zz[j] * fl1_fx + 2.0 * pa_xx[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_zzzz[j] + pa_xxz[j] * pb_zzzz[j]);

                t_xxz_zzzz[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb - 3.0 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 3.0 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - 3.0 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa_xx[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] + 12.0 * fl3_fx * fl1_fz * pb_z[j] + 7.5 * pa_xxz[j] * fl1_fz * fl2_fx + 30.0 * pa_xx[j] * fl1_fz * fl2_fx * pb_z[j] - 3.0 * fl1_fx * pa_z[j] * pb_zz[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_z[j] * pb_zz[j] * fl1_fx - 2.0 * fl1_fz * fl1_fga * fl1_fx * pb_zzz[j] - 6.0 * pa_xxz[j] * pb_zz[j] * fl1_fz * fl1_fgb + 15.0 * fl2_fx * fl1_fz * pa_z[j] * pb_zz[j] + 10.0 * fl2_fx * fl1_fz * pb_zzz[j] + 36.0 * pa_xxz[j] * fl1_fz * pb_zz[j] * fl1_fx + 24.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_zzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_zzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_zzzz[j] + 14.0 * pa_xxz[j] * fl1_fz * pb_zzzz[j]);

                t_xyy_xxxx[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 1.5 * fl3_fx * pb_x[j] + 0.75 * pa_xyy[j] * fl2_fx + 3.0 * fl2_fx * pa_yy[j] * pb_x[j] + 1.5 * pa_x[j] * fl2_fx * pb_xx[j] + fl2_fx * pb_xxx[j] + 3.0 * pa_xyy[j] * pb_xx[j] * fl1_fx + 2.0 * fl1_fx * pa_yy[j] * pb_xxx[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxxx[j] + pa_xyy[j] * pb_xxxx[j]);

                t_xyy_xxxx[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 3.0 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 3.0 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx - 6.0 * fl1_fx * pa_yy[j] * pb_x[j] * fl1_fz * fl1_fgb + 12.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx + 30.0 * fl2_fx * pa_yy[j] * fl1_fz * pb_x[j] - 3.0 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 2.0 * fl1_fx * fl1_fz * fl1_fga * pb_xxx[j] - 6.0 * pa_xyy[j] * pb_xx[j] * fl1_fz * fl1_fgb + 15.0 * pa_x[j] * fl1_fz * fl2_fx * pb_xx[j] + 10.0 * fl2_fx * fl1_fz * pb_xxx[j] + 36.0 * pa_xyy[j] * fl1_fz * pb_xx[j] * fl1_fx + 24.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xxx[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxxx[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxx[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xxxx[j]);

                t_xyy_xxxy[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_y[j] + 0.375 * fl3_fx * pb_y[j] + 1.5 * pa_xy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yy[j] * pb_y[j] + 1.5 * fl2_fx * pa_y[j] * pb_xx[j] + 0.75 * pa_x[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pb_xxy[j] + 1.5 * pa_xyy[j] * pb_xy[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_yy[j] * pb_xxy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxxy[j] + pa_xyy[j] * pb_xxxy[j]);

                t_xyy_xxxy[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * pa_y[j] * fl1_fz - 0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 3.0 * pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_yy[j] * fl1_fz * fl1_fgb * pb_y[j] + 3.0 * fl3_fx * fl1_fz * pb_y[j] + 15.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_x[j] + 7.5 * fl2_fx * pa_yy[j] * fl1_fz * pb_y[j] + 15.0 * fl2_fx * pa_y[j] * fl1_fz * pb_xx[j] - 1.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxy[j] - 3.0 * pa_xyy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pb_xxy[j] + 18.0 * pa_xyy[j] * fl1_fz * pb_xy[j] * fl1_fx + 12.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xxx[j] + 18.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xxy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxxy[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxy[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xxxy[j]);

                t_xyy_xxxz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_z[j] + 0.75 * fl2_fx * pa_yy[j] * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pb_xxz[j] + 1.5 * pa_xyy[j] * pb_xz[j] * fl1_fx + 1.5 * fl1_fx * pa_yy[j] * pb_xxz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxxz[j] + pa_xyy[j] * pb_xxxz[j]);

                t_xyy_xxxz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 1.5 * fl1_fx * pa_yy[j] * fl1_fz * fl1_fgb * pb_z[j] + 3.0 * fl3_fx * fl1_fz * pb_z[j] + 7.5 * fl2_fx * pa_yy[j] * fl1_fz * pb_z[j] - 1.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxz[j] - 3.0 * pa_xyy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pb_xxz[j] + 18.0 * pa_xyy[j] * fl1_fz * pb_xz[j] * fl1_fx + 18.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xxz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxxz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xxxz[j]);

                t_xyy_xxyy[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * fl3_fx * pb_x[j] + 0.25 * pa_xyy[j] * fl2_fx + pa_xy[j] * fl2_fx * pb_y[j] + 0.75 * pa_x[j] * fl2_fx * pb_xx[j] + 0.5 * fl2_fx * pa_yy[j] * pb_x[j] + 2.0 * fl2_fx * pa_y[j] * pb_xy[j] + 0.25 * pa_x[j] * fl2_fx * pb_yy[j] + 0.5 * fl2_fx * pb_xyy[j] + 0.5 * pa_xyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyy[j] * fl1_fx * pb_yy[j] + 2.0 * pa_xy[j] * fl1_fx * pb_xxy[j] + fl1_fx * pa_yy[j] * pb_xyy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxyy[j] + pa_xyy[j] * pb_xxyy[j]);

                t_xyy_xxyy[j] += fl_r_0_0 * (-pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 6.0 * fl3_fx * fl1_fz * pb_x[j] - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - fl1_fx * pa_yy[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_xyy[j] * fl1_fz * fl2_fx + 10.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_y[j] + 7.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * fl2_fx * pa_yy[j] * fl1_fz * pb_x[j] + 20.0 * fl2_fx * pa_y[j] * fl1_fz * pb_xy[j] - 0.5 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yy[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - fl1_fx * fl1_fz * fl1_fga * pb_xyy[j] - pa_xyy[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xyy[j] * fl1_fz * fl1_fgb * pb_yy[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yy[j] + 5.0 * fl2_fx * fl1_fz * pb_xyy[j] + 6.0 * pa_xyy[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xyy[j] * fl1_fz * fl1_fx * pb_yy[j] + 24.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xyy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxyy[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxyy[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xxyy[j]);

                t_xyy_xxyz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx * pb_z[j] + fl2_fx * pa_y[j] * pb_xz[j] + 0.25 * pa_x[j] * fl2_fx * pb_yz[j] + 0.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_xyy[j] * fl1_fx * pb_yz[j] + pa_xy[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_yy[j] * pb_xyz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxyz[j] + pa_xyy[j] * pb_xxyz[j]);

                t_xyy_xxyz[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 5.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_z[j] + 10.0 * fl2_fx * pa_y[j] * fl1_fz * pb_xz[j] - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yz[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - fl1_fx * fl1_fz * fl1_fga * pb_xyz[j] - pa_xyy[j] * fl1_fz * fl1_fgb * pb_yz[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yz[j] + 5.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_xyy[j] * fl1_fz * fl1_fx * pb_yz[j] + 12.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xyz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxyz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxyz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_50_60(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_yy = paDistances.data(19 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyy = paDistances.data(19 * idx + 12);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyy_xxzz = primBuffer.data(150 * idx + 50);

            auto t_xyy_xyyy = primBuffer.data(150 * idx + 51);

            auto t_xyy_xyyz = primBuffer.data(150 * idx + 52);

            auto t_xyy_xyzz = primBuffer.data(150 * idx + 53);

            auto t_xyy_xzzz = primBuffer.data(150 * idx + 54);

            auto t_xyy_yyyy = primBuffer.data(150 * idx + 55);

            auto t_xyy_yyyz = primBuffer.data(150 * idx + 56);

            auto t_xyy_yyzz = primBuffer.data(150 * idx + 57);

            auto t_xyy_yzzz = primBuffer.data(150 * idx + 58);

            auto t_xyy_zzzz = primBuffer.data(150 * idx + 59);

            // Batch of Integrals (50,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     r_0_0, s_0_0, t_xyy_xxzz, t_xyy_xyyy, t_xyy_xyyz, t_xyy_xyzz, t_xyy_xzzz, \
                                     t_xyy_yyyy, t_xyy_yyyz, t_xyy_yyzz, t_xyy_yzzz, t_xyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyy_xxzz[j] = fl_s_0_0 * (0.125 * pa_x[j] * fl3_fx + 0.25 * fl3_fx * pb_x[j] + 0.25 * pa_xyy[j] * fl2_fx + 0.5 * fl2_fx * pa_yy[j] * pb_x[j] + 0.25 * pa_x[j] * fl2_fx * pb_xx[j] + 0.25 * pa_x[j] * fl2_fx * pb_zz[j] + 0.5 * fl2_fx * pb_xzz[j] + 0.5 * pa_xyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyy[j] * fl1_fx * pb_zz[j] + fl1_fx * pa_yy[j] * pb_xzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxzz[j] + pa_xyy[j] * pb_xxzz[j]);

                t_xyy_xxzz[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + pa_x[j] * fl1_fz * fl3_fx - fl1_fx * pa_yy[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl3_fx * fl1_fz * pb_x[j] + 2.5 * pa_xyy[j] * fl1_fz * fl2_fx + 5.0 * fl2_fx * pa_yy[j] * fl1_fz * pb_x[j] - 0.5 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - fl1_fx * fl1_fz * fl1_fga * pb_xzz[j] - pa_xyy[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xyy[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xx[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_zz[j] + 5.0 * fl2_fx * fl1_fz * pb_xzz[j] + 6.0 * pa_xyy[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xyy[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxzz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xxzz[j]);

                t_xyy_xyyy[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_y[j] + 1.125 * fl3_fx * pb_y[j] + 1.5 * pa_xy[j] * fl2_fx * pb_x[j] + 2.25 * pa_x[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_yy[j] * pb_y[j] + 1.5 * fl2_fx * pa_y[j] * pb_yy[j] + 0.25 * fl2_fx * pb_yyy[j] + 1.5 * pa_xyy[j] * pb_xy[j] * fl1_fx + 3.0 * pa_xy[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yyy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyyy[j] + pa_xyy[j] * pb_xyyy[j]);

                t_xyy_xyyy[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * pa_y[j] * fl1_fz + 9.0 * fl3_fx * fl1_fz * pb_y[j] - 0.75 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 3.0 * pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_yy[j] * pb_y[j] * fl1_fz * fl1_fgb + 15.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_x[j] + 22.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xy[j] + 7.5 * fl2_fx * pa_yy[j] * fl1_fz * pb_y[j] + 15.0 * fl2_fx * pa_y[j] * fl1_fz * pb_yy[j] - 1.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyy[j] - 3.0 * pa_xyy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_yyy[j] + 18.0 * pa_xyy[j] * fl1_fz * pb_xy[j] * fl1_fx + 36.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xyy[j] + 6.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_yyy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyyy[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyyy[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xyyy[j]);

                t_xyy_xyyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pa_yy[j] * pb_z[j] + fl2_fx * pa_y[j] * pb_yz[j] + 0.25 * fl2_fx * pb_yyz[j] + 0.5 * pa_xyy[j] * pb_xz[j] * fl1_fx + 2.0 * pa_xy[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yyz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyyz[j] + pa_xyy[j] * pb_xyyz[j]);

                t_xyy_xyyz[j] += fl_r_0_0 * (3.0 * fl3_fx * fl1_fz * pb_z[j] - 0.25 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.25 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 0.5 * fl1_fx * pa_yy[j] * fl1_fz * fl1_fgb * pb_z[j] + 7.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xz[j] + 2.5 * fl2_fx * pa_yy[j] * fl1_fz * pb_z[j] + 10.0 * fl2_fx * pa_y[j] * fl1_fz * pb_yz[j] - 0.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyz[j] - pa_xyy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_yyz[j] + 6.0 * pa_xyy[j] * fl1_fz * pb_xz[j] * fl1_fx + 24.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xyz[j] + 6.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_yyz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyyz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyyz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xyyz[j]);

                t_xyy_xyzz[j] = fl_s_0_0 * (0.25 * fl3_fx * pa_y[j] + 0.125 * fl3_fx * pb_y[j] + 0.5 * pa_xy[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pa_yy[j] * pb_y[j] + 0.5 * fl2_fx * pa_y[j] * pb_zz[j] + 0.25 * pa_x[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pb_yzz[j] + 0.5 * pa_xyy[j] * pb_xy[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyzz[j] + pa_xyy[j] * pb_xyzz[j]);

                t_xyy_xyzz[j] += fl_r_0_0 * (-0.5 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + 2.0 * fl3_fx * pa_y[j] * fl1_fz - 0.25 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_yy[j] * pb_y[j] * fl1_fz * fl1_fgb + fl3_fx * fl1_fz * pb_y[j] + 5.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_x[j] + 2.5 * fl2_fx * pa_yy[j] * fl1_fz * pb_y[j] + 5.0 * fl2_fx * pa_y[j] * fl1_fz * pb_zz[j] - 0.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yzz[j] - pa_xyy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xy[j] + 2.5 * fl2_fx * fl1_fz * pb_yzz[j] + 6.0 * pa_xyy[j] * fl1_fz * pb_xy[j] * fl1_fx + 12.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xzz[j] + 6.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_yzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyzz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xyzz[j]);

                t_xyy_xzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_z[j] + 0.75 * fl2_fx * pa_yy[j] * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pb_zzz[j] + 1.5 * pa_xyy[j] * pb_xz[j] * fl1_fx + 0.5 * fl1_fx * pa_yy[j] * pb_zzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xzzz[j] + pa_xyy[j] * pb_xzzz[j]);

                t_xyy_xzzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 1.5 * fl1_fx * pa_yy[j] * pb_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_z[j] + 7.5 * fl2_fx * pa_yy[j] * fl1_fz * pb_z[j] - 1.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_zzz[j] - 3.0 * pa_xyy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xz[j] + 2.5 * fl2_fx * fl1_fz * pb_zzz[j] + 18.0 * pa_xyy[j] * fl1_fz * pb_xz[j] * fl1_fx + 6.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_zzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xzzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xzzz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_xzzz[j]);

                t_xyy_yyyy[j] = fl_s_0_0 * (1.875 * pa_x[j] * fl3_fx + 0.75 * pa_xyy[j] * fl2_fx + 6.0 * pa_xy[j] * fl2_fx * pb_y[j] + 4.5 * pa_x[j] * fl2_fx * pb_yy[j] + 3.0 * pa_xyy[j] * pb_yy[j] * fl1_fx + 4.0 * pa_xy[j] * fl1_fx * pb_yyy[j] + 0.5 * pa_x[j] * fl1_fx * pb_yyyy[j] + pa_xyy[j] * pb_yyyy[j]);

                t_xyy_yyyy[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa_xy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx + 60.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_y[j] + 45.0 * pa_x[j] * fl2_fx * fl1_fz * pb_yy[j] - 3.0 * pa_x[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 6.0 * pa_xyy[j] * pb_yy[j] * fl1_fz * fl1_fgb + 36.0 * pa_xyy[j] * fl1_fz * pb_yy[j] * fl1_fx + 48.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yyy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yyyy[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyyy[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_yyyy[j]);

                t_xyy_yyyz[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl2_fx * pb_z[j] + 2.25 * pa_x[j] * fl2_fx * pb_yz[j] + 1.5 * pa_xyy[j] * pb_yz[j] * fl1_fx + 3.0 * pa_xy[j] * fl1_fx * pb_yyz[j] + 0.5 * pa_x[j] * fl1_fx * pb_yyyz[j] + pa_xyy[j] * pb_yyyz[j]);

                t_xyy_yyyz[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 15.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_z[j] + 22.5 * pa_x[j] * fl2_fx * fl1_fz * pb_yz[j] - 1.5 * pa_x[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 3.0 * pa_xyy[j] * pb_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa_xyy[j] * fl1_fz * pb_yz[j] * fl1_fx + 36.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yyz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yyyz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyyz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_yyyz[j]);

                t_xyy_yyzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xyy[j] * fl2_fx + pa_xy[j] * fl2_fx * pb_y[j] + 0.75 * pa_x[j] * fl2_fx * pb_zz[j] + 0.25 * pa_x[j] * fl2_fx * pb_yy[j] + 0.5 * pa_xyy[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xyy[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xy[j] * fl1_fx * pb_yzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_yyzz[j] + pa_xyy[j] * pb_yyzz[j]);

                t_xyy_yyzz[j] += fl_r_0_0 * (-pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_xy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 2.5 * pa_xyy[j] * fl1_fz * fl2_fx + 10.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_y[j] + 7.5 * pa_x[j] * fl2_fx * fl1_fz * pb_zz[j] - 0.5 * pa_x[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - pa_xyy[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_xyy[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yy[j] + 6.0 * pa_xyy[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_xyy[j] * fl1_fz * fl1_fx * pb_zz[j] + 24.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yyzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyzz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_yyzz[j]);

                t_xyy_yzzz[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl2_fx * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_yz[j] + 1.5 * pa_xyy[j] * pb_yz[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_zzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_yzzz[j] + pa_xyy[j] * pb_yzzz[j]);

                t_xyy_yzzz[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 15.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_z[j] - 1.5 * pa_x[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 3.0 * pa_xyy[j] * pb_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yz[j] + 18.0 * pa_xyy[j] * fl1_fz * pb_yz[j] * fl1_fx + 12.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_zzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yzzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yzzz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_yzzz[j]);

                t_xyy_zzzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xyy[j] * fl2_fx + 1.5 * pa_x[j] * fl2_fx * pb_zz[j] + 3.0 * pa_xyy[j] * pb_zz[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_zzzz[j] + pa_xyy[j] * pb_zzzz[j]);

                t_xyy_zzzz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx - 3.0 * pa_x[j] * fl1_fx * pb_zz[j] * fl1_fz * fl1_fgb - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_zz[j] * fl1_fx - 6.0 * pa_xyy[j] * pb_zz[j] * fl1_fz * fl1_fgb + 15.0 * pa_x[j] * fl1_fz * fl2_fx * pb_zz[j] + 36.0 * pa_xyy[j] * fl1_fz * pb_zz[j] * fl1_fx - pa_x[j] * fl1_fz * fl1_fga * pb_zzzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_zzzz[j] + 14.0 * pa_xyy[j] * fl1_fz * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_60_70(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yz = paDistances.data(19 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyz = paDistances.data(19 * idx + 13);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyz_xxxx = primBuffer.data(150 * idx + 60);

            auto t_xyz_xxxy = primBuffer.data(150 * idx + 61);

            auto t_xyz_xxxz = primBuffer.data(150 * idx + 62);

            auto t_xyz_xxyy = primBuffer.data(150 * idx + 63);

            auto t_xyz_xxyz = primBuffer.data(150 * idx + 64);

            auto t_xyz_xxzz = primBuffer.data(150 * idx + 65);

            auto t_xyz_xyyy = primBuffer.data(150 * idx + 66);

            auto t_xyz_xyyz = primBuffer.data(150 * idx + 67);

            auto t_xyz_xyzz = primBuffer.data(150 * idx + 68);

            auto t_xyz_xzzz = primBuffer.data(150 * idx + 69);

            // Batch of Integrals (60,70)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyz_xxxx, t_xyz_xxxy, \
                                     t_xyz_xxxz, t_xyz_xxyy, t_xyz_xxyz, t_xyz_xxzz, t_xyz_xyyy, t_xyz_xyyz, t_xyz_xyzz, \
                                     t_xyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyz_xxxx[j] = fl_s_0_0 * (0.75 * pa_xyz[j] * fl2_fx + 3.0 * fl2_fx * pa_yz[j] * pb_x[j] + 3.0 * pa_xyz[j] * pb_xx[j] * fl1_fx + 2.0 * fl1_fx * pa_yz[j] * pb_xxx[j] + pa_xyz[j] * pb_xxxx[j]);

                t_xyz_xxxx[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * fl1_fx * pa_yz[j] * pb_x[j] * fl1_fz * fl1_fgb + 7.5 * pa_xyz[j] * fl1_fz * fl2_fx + 30.0 * fl2_fx * pa_yz[j] * fl1_fz * pb_x[j] - 6.0 * pa_xyz[j] * pb_xx[j] * fl1_fz * fl1_fgb + 36.0 * pa_xyz[j] * fl1_fz * pb_xx[j] * fl1_fx + 24.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xxx[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xxxx[j]);

                t_xyz_xxxy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 0.75 * pa_xz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yz[j] * pb_y[j] + 0.75 * fl2_fx * pa_z[j] * pb_xx[j] + 1.5 * pa_xyz[j] * pb_xy[j] * fl1_fx + 0.5 * pa_xz[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_yz[j] * pb_xxy[j] + pa_xyz[j] * pb_xxxy[j]);

                t_xyz_xxxy[j] += fl_r_0_0 * (-0.75 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] - 1.5 * pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_yz[j] * fl1_fz * fl1_fgb * pb_y[j] + 7.5 * pa_xz[j] * fl2_fx * fl1_fz * pb_x[j] + 7.5 * fl2_fx * pa_yz[j] * fl1_fz * pb_y[j] + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xx[j] - 3.0 * pa_xyz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa_xyz[j] * fl1_fz * pb_xy[j] * fl1_fx + 6.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xxx[j] + 18.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xxy[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xxxy[j]);

                t_xyz_xxxz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_y[j] + 0.75 * pa_xy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yz[j] * pb_z[j] + 0.75 * fl2_fx * pa_y[j] * pb_xx[j] + 1.5 * pa_xyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xy[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_yz[j] * pb_xxz[j] + pa_xyz[j] * pb_xxxz[j]);

                t_xyz_xxxz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * pa_y[j] * fl1_fz - 1.5 * pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_yz[j] * fl1_fz * fl1_fgb * pb_z[j] + 7.5 * pa_xy[j] * fl1_fz * fl2_fx * pb_x[j] + 7.5 * fl2_fx * pa_yz[j] * fl1_fz * pb_z[j] + 7.5 * fl2_fx * pa_y[j] * fl1_fz * pb_xx[j] - 3.0 * pa_xyz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa_xyz[j] * fl1_fz * pb_xz[j] * fl1_fx + 6.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xxx[j] + 18.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xxz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xxxz[j]);

                t_xyz_xxyy[j] = fl_s_0_0 * (0.25 * pa_xyz[j] * fl2_fx + 0.5 * pa_xz[j] * fl2_fx * pb_y[j] + 0.5 * fl2_fx * pa_yz[j] * pb_x[j] + fl2_fx * pa_z[j] * pb_xy[j] + 0.5 * pa_xyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyz[j] * fl1_fx * pb_yy[j] + pa_xz[j] * fl1_fx * pb_xxy[j] + fl1_fx * pa_yz[j] * pb_xyy[j] + pa_xyz[j] * pb_xxyy[j]);

                t_xyz_xxyy[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - fl1_fx * pa_yz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_xyz[j] * fl1_fz * fl2_fx + 5.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_y[j] + 5.0 * fl2_fx * pa_yz[j] * fl1_fz * pb_x[j] + 10.0 * fl2_fx * fl1_fz * pa_z[j] * pb_xy[j] - pa_xyz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xyz[j] * fl1_fz * fl1_fgb * pb_yy[j] + 6.0 * pa_xyz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xyz[j] * fl1_fz * fl1_fx * pb_yy[j] + 12.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xxy[j] + 12.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xyy[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xxyy[j]);

                t_xyz_xxyz[j] = fl_s_0_0 * (0.125 * pa_x[j] * fl3_fx + 0.25 * fl3_fx * pb_x[j] + 0.25 * pa_xy[j] * fl2_fx * pb_y[j] + 0.25 * pa_xz[j] * fl2_fx * pb_z[j] + 0.25 * pa_x[j] * fl2_fx * pb_xx[j] + 0.5 * fl2_fx * pa_y[j] * pb_xy[j] + 0.5 * fl2_fx * pa_z[j] * pb_xz[j] + 0.5 * pa_xyz[j] * fl1_fx * pb_yz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xxy[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_yz[j] * pb_xyz[j] + pa_xyz[j] * pb_xxyz[j]);

                t_xyz_xxyz[j] += fl_r_0_0 * (-0.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + pa_x[j] * fl3_fx * fl1_fz + 2.0 * fl3_fx * fl1_fz * pb_x[j] - 0.5 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 2.5 * pa_xy[j] * fl1_fz * fl2_fx * pb_y[j] + 2.5 * pa_xz[j] * fl2_fx * fl1_fz * pb_z[j] + 2.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * fl2_fx * pa_y[j] * fl1_fz * pb_xy[j] + 5.0 * fl2_fx * fl1_fz * pa_z[j] * pb_xz[j] - pa_xyz[j] * fl1_fz * fl1_fgb * pb_yz[j] + 6.0 * pa_xyz[j] * fl1_fz * fl1_fx * pb_yz[j] + 6.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xxy[j] + 6.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xxz[j] + 12.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xyz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xxyz[j]);

                t_xyz_xxzz[j] = fl_s_0_0 * (0.25 * pa_xyz[j] * fl2_fx + 0.5 * pa_xy[j] * fl2_fx * pb_z[j] + 0.5 * fl2_fx * pa_yz[j] * pb_x[j] + fl2_fx * pa_y[j] * pb_xz[j] + 0.5 * pa_xyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyz[j] * fl1_fx * pb_zz[j] + pa_xy[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_yz[j] * pb_xzz[j] + pa_xyz[j] * pb_xxzz[j]);

                t_xyz_xxzz[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - fl1_fx * pa_yz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_xyz[j] * fl1_fz * fl2_fx + 5.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_z[j] + 5.0 * fl2_fx * pa_yz[j] * fl1_fz * pb_x[j] + 10.0 * fl2_fx * pa_y[j] * fl1_fz * pb_xz[j] - pa_xyz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xyz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 6.0 * pa_xyz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xyz[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xzz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xxzz[j]);

                t_xyz_xyyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 0.75 * pa_xz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yz[j] * pb_y[j] + 0.75 * fl2_fx * pa_z[j] * pb_yy[j] + 1.5 * pa_xyz[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xz[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yyy[j] + pa_xyz[j] * pb_xyyy[j]);

                t_xyz_xyyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] - 1.5 * pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_yz[j] * pb_y[j] * fl1_fz * fl1_fgb + 7.5 * pa_xz[j] * fl2_fx * fl1_fz * pb_x[j] + 7.5 * fl2_fx * pa_yz[j] * fl1_fz * pb_y[j] + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yy[j] - 3.0 * pa_xyz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa_xyz[j] * fl1_fz * pb_xy[j] * fl1_fx + 18.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xyy[j] + 6.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_yyy[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xyyy[j]);

                t_xyz_xyyz[j] = fl_s_0_0 * (0.125 * fl3_fx * pa_y[j] + 0.25 * fl3_fx * pb_y[j] + 0.25 * pa_xy[j] * fl2_fx * pb_x[j] + 0.5 * pa_x[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pa_yz[j] * pb_z[j] + 0.25 * fl2_fx * pa_y[j] * pb_yy[j] + 0.5 * fl2_fx * pa_z[j] * pb_yz[j] + 0.5 * pa_xyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xy[j] * fl1_fx * pb_xyy[j] + pa_xz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yyz[j] + pa_xyz[j] * pb_xyyz[j]);

                t_xyz_xyyz[j] += fl_r_0_0 * (-0.25 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + fl3_fx * pa_y[j] * fl1_fz + 2.0 * fl3_fx * fl1_fz * pb_y[j] - 0.5 * pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_yz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.5 * pa_xy[j] * fl1_fz * fl2_fx * pb_x[j] + 5.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xy[j] + 2.5 * fl2_fx * pa_yz[j] * fl1_fz * pb_z[j] + 2.5 * fl2_fx * pa_y[j] * fl1_fz * pb_yy[j] + 5.0 * fl2_fx * fl1_fz * pa_z[j] * pb_yz[j] - pa_xyz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 6.0 * pa_xyz[j] * fl1_fz * pb_xz[j] * fl1_fx + 6.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xyz[j] + 6.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_yyz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xyyz[j]);

                t_xyz_xyzz[j] = fl_s_0_0 * (0.125 * fl3_fx * pa_z[j] + 0.25 * fl3_fx * pb_z[j] + 0.25 * pa_xz[j] * fl2_fx * pb_x[j] + 0.5 * pa_x[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pa_yz[j] * pb_y[j] + 0.5 * fl2_fx * pa_y[j] * pb_yz[j] + 0.25 * fl2_fx * pa_z[j] * pb_zz[j] + 0.5 * pa_xyz[j] * pb_xy[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_xyz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yzz[j] + pa_xyz[j] * pb_xyzz[j]);

                t_xyz_xyzz[j] += fl_r_0_0 * (-0.25 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + fl3_fx * fl1_fz * pa_z[j] + 2.0 * fl3_fx * fl1_fz * pb_z[j] - 0.5 * pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_yz[j] * pb_y[j] * fl1_fz * fl1_fgb + 2.5 * pa_xz[j] * fl2_fx * fl1_fz * pb_x[j] + 5.0 * pa_x[j] * fl2_fx * fl1_fz * pb_xz[j] + 2.5 * fl2_fx * pa_yz[j] * fl1_fz * pb_y[j] + 5.0 * fl2_fx * pa_y[j] * fl1_fz * pb_yz[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_zz[j] - pa_xyz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 6.0 * pa_xyz[j] * fl1_fz * pb_xy[j] * fl1_fx + 12.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xyz[j] + 6.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xzz[j] + 6.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_yzz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xyzz[j]);

                t_xyz_xzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_y[j] + 0.75 * pa_xy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yz[j] * pb_z[j] + 0.75 * fl2_fx * pa_y[j] * pb_zz[j] + 1.5 * pa_xyz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xy[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_zzz[j] + pa_xyz[j] * pb_xzzz[j]);

                t_xyz_xzzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pa_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * pa_y[j] * fl1_fz - 1.5 * pa_xy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_yz[j] * pb_z[j] * fl1_fz * fl1_fgb + 7.5 * pa_xy[j] * fl1_fz * fl2_fx * pb_x[j] + 7.5 * fl2_fx * pa_yz[j] * fl1_fz * pb_z[j] + 7.5 * fl2_fx * pa_y[j] * fl1_fz * pb_zz[j] - 3.0 * pa_xyz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa_xyz[j] * fl1_fz * pb_xz[j] * fl1_fx + 18.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xzz[j] + 6.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_zzz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_70_80(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyz_yyyy = primBuffer.data(150 * idx + 70);

            auto t_xyz_yyyz = primBuffer.data(150 * idx + 71);

            auto t_xyz_yyzz = primBuffer.data(150 * idx + 72);

            auto t_xyz_yzzz = primBuffer.data(150 * idx + 73);

            auto t_xyz_zzzz = primBuffer.data(150 * idx + 74);

            auto t_xzz_xxxx = primBuffer.data(150 * idx + 75);

            auto t_xzz_xxxy = primBuffer.data(150 * idx + 76);

            auto t_xzz_xxxz = primBuffer.data(150 * idx + 77);

            auto t_xzz_xxyy = primBuffer.data(150 * idx + 78);

            auto t_xzz_xxyz = primBuffer.data(150 * idx + 79);

            // Batch of Integrals (70,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, \
                                     pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xyz_yyyy, t_xyz_yyyz, t_xyz_yyzz, \
                                     t_xyz_yzzz, t_xyz_zzzz, t_xzz_xxxx, t_xzz_xxxy, t_xzz_xxxz, t_xzz_xxyy, t_xzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyz_yyyy[j] = fl_s_0_0 * (0.75 * pa_xyz[j] * fl2_fx + 3.0 * pa_xz[j] * fl2_fx * pb_y[j] + 3.0 * pa_xyz[j] * pb_yy[j] * fl1_fx + 2.0 * pa_xz[j] * fl1_fx * pb_yyy[j] + pa_xyz[j] * pb_yyyy[j]);

                t_xyz_yyyy[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa_xz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 7.5 * pa_xyz[j] * fl1_fz * fl2_fx + 30.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_y[j] - 6.0 * pa_xyz[j] * pb_yy[j] * fl1_fz * fl1_fgb + 36.0 * pa_xyz[j] * fl1_fz * pb_yy[j] * fl1_fx + 24.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yyy[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_yyyy[j]);

                t_xyz_yyyz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xy[j] * fl2_fx * pb_y[j] + 0.75 * pa_xz[j] * fl2_fx * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_yy[j] + 1.5 * pa_xyz[j] * pb_yz[j] * fl1_fx + 0.5 * pa_xy[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_yyz[j] + pa_xyz[j] * pb_yyyz[j]);

                t_xyz_yyyz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz - 1.5 * pa_xy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 7.5 * pa_xy[j] * fl1_fz * fl2_fx * pb_y[j] + 7.5 * pa_xz[j] * fl2_fx * fl1_fz * pb_z[j] + 7.5 * pa_x[j] * fl2_fx * fl1_fz * pb_yy[j] - 3.0 * pa_xyz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa_xyz[j] * fl1_fz * pb_yz[j] * fl1_fx + 6.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yyy[j] + 18.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yyz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_yyyz[j]);

                t_xyz_yyzz[j] = fl_s_0_0 * (0.25 * pa_xyz[j] * fl2_fx + 0.5 * pa_xy[j] * fl2_fx * pb_z[j] + 0.5 * pa_xz[j] * fl2_fx * pb_y[j] + pa_x[j] * fl2_fx * pb_yz[j] + 0.5 * pa_xyz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xyz[j] * fl1_fx * pb_zz[j] + pa_xy[j] * fl1_fx * pb_yyz[j] + pa_xz[j] * fl1_fx * pb_yzz[j] + pa_xyz[j] * pb_yyzz[j]);

                t_xyz_yyzz[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - pa_xz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 2.5 * pa_xyz[j] * fl1_fz * fl2_fx + 5.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_z[j] + 5.0 * pa_xz[j] * fl2_fx * fl1_fz * pb_y[j] + 10.0 * pa_x[j] * fl2_fx * fl1_fz * pb_yz[j] - pa_xyz[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_xyz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 6.0 * pa_xyz[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_xyz[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yyz[j] + 12.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yzz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_yyzz[j]);

                t_xyz_yzzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xy[j] * fl2_fx * pb_y[j] + 0.75 * pa_xz[j] * fl2_fx * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_zz[j] + 1.5 * pa_xyz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xy[j] * fl1_fx * pb_yzz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_zzz[j] + pa_xyz[j] * pb_yzzz[j]);

                t_xyz_yzzz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz - 1.5 * pa_xy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * pa_xz[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 7.5 * pa_xy[j] * fl1_fz * fl2_fx * pb_y[j] + 7.5 * pa_xz[j] * fl2_fx * fl1_fz * pb_z[j] + 7.5 * pa_x[j] * fl2_fx * fl1_fz * pb_zz[j] - 3.0 * pa_xyz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa_xyz[j] * fl1_fz * pb_yz[j] * fl1_fx + 18.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yzz[j] + 6.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_zzz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_yzzz[j]);

                t_xyz_zzzz[j] = fl_s_0_0 * (0.75 * pa_xyz[j] * fl2_fx + 3.0 * pa_xy[j] * fl2_fx * pb_z[j] + 3.0 * pa_xyz[j] * pb_zz[j] * fl1_fx + 2.0 * pa_xy[j] * fl1_fx * pb_zzz[j] + pa_xyz[j] * pb_zzzz[j]);

                t_xyz_zzzz[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa_xy[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 7.5 * pa_xyz[j] * fl1_fz * fl2_fx + 30.0 * pa_xy[j] * fl1_fz * fl2_fx * pb_z[j] - 6.0 * pa_xyz[j] * pb_zz[j] * fl1_fz * fl1_fgb + 36.0 * pa_xyz[j] * fl1_fz * pb_zz[j] * fl1_fx + 24.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_zzz[j] + 14.0 * pa_xyz[j] * fl1_fz * pb_zzzz[j]);

                t_xzz_xxxx[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 1.5 * fl3_fx * pb_x[j] + 0.75 * pa_xzz[j] * fl2_fx + 3.0 * fl2_fx * pa_zz[j] * pb_x[j] + 1.5 * pa_x[j] * fl2_fx * pb_xx[j] + fl2_fx * pb_xxx[j] + 3.0 * pa_xzz[j] * pb_xx[j] * fl1_fx + 2.0 * fl1_fx * pa_zz[j] * pb_xxx[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxxx[j] + pa_xzz[j] * pb_xxxx[j]);

                t_xzz_xxxx[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 3.0 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 3.0 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx - 6.0 * fl1_fx * pa_zz[j] * pb_x[j] * fl1_fz * fl1_fgb + 12.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * pa_xzz[j] * fl1_fz * fl2_fx + 30.0 * fl2_fx * pa_zz[j] * fl1_fz * pb_x[j] - 3.0 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 2.0 * fl1_fx * fl1_fz * fl1_fga * pb_xxx[j] - 6.0 * pa_xzz[j] * pb_xx[j] * fl1_fz * fl1_fgb + 15.0 * pa_x[j] * fl1_fz * fl2_fx * pb_xx[j] + 10.0 * fl2_fx * fl1_fz * pb_xxx[j] + 36.0 * pa_xzz[j] * fl1_fz * pb_xx[j] * fl1_fx + 24.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xxx[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxxx[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxx[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xxxx[j]);

                t_xzz_xxxy[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_y[j] + 0.75 * fl2_fx * pa_zz[j] * pb_y[j] + 0.75 * pa_x[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pb_xxy[j] + 1.5 * pa_xzz[j] * pb_xy[j] * fl1_fx + 1.5 * fl1_fx * pa_zz[j] * pb_xxy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxxy[j] + pa_xzz[j] * pb_xxxy[j]);

                t_xzz_xxxy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 1.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb * pb_y[j] + 3.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_y[j] - 1.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxy[j] - 3.0 * pa_xzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pb_xxy[j] + 18.0 * pa_xzz[j] * fl1_fz * pb_xy[j] * fl1_fx + 18.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xxy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxxy[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxy[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xxxy[j]);

                t_xzz_xxxz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] + 0.375 * fl3_fx * pb_z[j] + 1.5 * pa_xz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_zz[j] * pb_z[j] + 1.5 * fl2_fx * pa_z[j] * pb_xx[j] + 0.75 * pa_x[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pb_xxz[j] + 1.5 * pa_xzz[j] * pb_xz[j] * fl1_fx + pa_xz[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_zz[j] * pb_xxz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxxz[j] + pa_xzz[j] * pb_xxxz[j]);

                t_xzz_xxxz[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * pa_z[j] * fl1_fz - 0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 3.0 * pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb * pb_z[j] + 3.0 * fl3_fx * fl1_fz * pb_z[j] + 15.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_x[j] + 7.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_z[j] + 15.0 * fl2_fx * pa_z[j] * fl1_fz * pb_xx[j] - 1.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxz[j] - 3.0 * pa_xzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pb_xxz[j] + 18.0 * pa_xzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 12.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xxx[j] + 18.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xxz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxxz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxxz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xxxz[j]);

                t_xzz_xxyy[j] = fl_s_0_0 * (0.125 * pa_x[j] * fl3_fx + 0.25 * fl3_fx * pb_x[j] + 0.25 * pa_xzz[j] * fl2_fx + 0.5 * fl2_fx * pa_zz[j] * pb_x[j] + 0.25 * pa_x[j] * fl2_fx * pb_xx[j] + 0.25 * pa_x[j] * fl2_fx * pb_yy[j] + 0.5 * fl2_fx * pb_xyy[j] + 0.5 * pa_xzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xzz[j] * fl1_fx * pb_yy[j] + fl1_fx * pa_zz[j] * pb_xyy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxyy[j] + pa_xzz[j] * pb_xxyy[j]);

                t_xzz_xxyy[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + pa_x[j] * fl1_fz * fl3_fx - fl1_fx * pa_zz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl3_fx * fl1_fz * pb_x[j] + 2.5 * pa_xzz[j] * fl1_fz * fl2_fx + 5.0 * fl2_fx * pa_zz[j] * fl1_fz * pb_x[j] - 0.5 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yy[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - fl1_fx * fl1_fz * fl1_fga * pb_xyy[j] - pa_xzz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xzz[j] * fl1_fz * fl1_fgb * pb_yy[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xx[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yy[j] + 5.0 * fl2_fx * fl1_fz * pb_xyy[j] + 6.0 * pa_xzz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xzz[j] * fl1_fz * fl1_fx * pb_yy[j] + 12.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xyy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxyy[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxyy[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xxyy[j]);

                t_xzz_xxyz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx * pb_y[j] + fl2_fx * pa_z[j] * pb_xy[j] + 0.25 * pa_x[j] * fl2_fx * pb_yz[j] + 0.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_xzz[j] * fl1_fx * pb_yz[j] + pa_xz[j] * fl1_fx * pb_xxy[j] + fl1_fx * pa_zz[j] * pb_xyz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxyz[j] + pa_xzz[j] * pb_xxyz[j]);

                t_xzz_xxyz[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] + 5.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_y[j] + 10.0 * fl2_fx * pa_z[j] * fl1_fz * pb_xy[j] - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yz[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - fl1_fx * fl1_fz * fl1_fga * pb_xyz[j] - pa_xzz[j] * fl1_fz * fl1_fgb * pb_yz[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yz[j] + 5.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_xzz[j] * fl1_fz * fl1_fx * pb_yz[j] + 12.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xyz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxyz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxyz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_80_90(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xzz = paDistances.data(19 * idx + 14);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzz_xxzz = primBuffer.data(150 * idx + 80);

            auto t_xzz_xyyy = primBuffer.data(150 * idx + 81);

            auto t_xzz_xyyz = primBuffer.data(150 * idx + 82);

            auto t_xzz_xyzz = primBuffer.data(150 * idx + 83);

            auto t_xzz_xzzz = primBuffer.data(150 * idx + 84);

            auto t_xzz_yyyy = primBuffer.data(150 * idx + 85);

            auto t_xzz_yyyz = primBuffer.data(150 * idx + 86);

            auto t_xzz_yyzz = primBuffer.data(150 * idx + 87);

            auto t_xzz_yzzz = primBuffer.data(150 * idx + 88);

            auto t_xzz_zzzz = primBuffer.data(150 * idx + 89);

            // Batch of Integrals (80,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, r_0_0, s_0_0, t_xzz_xxzz, t_xzz_xyyy, t_xzz_xyyz, t_xzz_xyzz, t_xzz_xzzz, \
                                     t_xzz_yyyy, t_xzz_yyyz, t_xzz_yyzz, t_xzz_yzzz, t_xzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xzz_xxzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * fl3_fx * pb_x[j] + 0.25 * pa_xzz[j] * fl2_fx + pa_xz[j] * fl2_fx * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_xx[j] + 0.5 * fl2_fx * pa_zz[j] * pb_x[j] + 2.0 * fl2_fx * pa_z[j] * pb_xz[j] + 0.25 * pa_x[j] * fl2_fx * pb_zz[j] + 0.5 * fl2_fx * pb_xzz[j] + 0.5 * pa_xzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xz[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_zz[j] * pb_xzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxzz[j] + pa_xzz[j] * pb_xxzz[j]);

                t_xzz_xxzz[j] += fl_r_0_0 * (-pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 6.0 * fl3_fx * fl1_fz * pb_x[j] - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - fl1_fx * pa_zz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_xzz[j] * fl1_fz * fl2_fx + 10.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_z[j] + 7.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * fl2_fx * pa_zz[j] * fl1_fz * pb_x[j] + 20.0 * fl2_fx * pa_z[j] * fl1_fz * pb_xz[j] - 0.5 * pa_x[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - fl1_fx * fl1_fz * fl1_fga * pb_xzz[j] - pa_xzz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_xzz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_zz[j] + 5.0 * fl2_fx * fl1_fz * pb_xzz[j] + 6.0 * pa_xzz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_xzz[j] * fl1_fz * fl1_fx * pb_zz[j] + 24.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxzz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xxzz[j]);

                t_xzz_xyyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_y[j] + 0.75 * fl2_fx * pa_zz[j] * pb_y[j] + 0.75 * pa_x[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pb_yyy[j] + 1.5 * pa_xzz[j] * pb_xy[j] * fl1_fx + 0.5 * fl1_fx * pa_zz[j] * pb_yyy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyyy[j] + pa_xzz[j] * pb_xyyy[j]);

                t_xzz_xyyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 1.5 * fl1_fx * pa_zz[j] * pb_y[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_y[j] - 1.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyy[j] - 3.0 * pa_xzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xy[j] + 2.5 * fl2_fx * fl1_fz * pb_yyy[j] + 18.0 * pa_xzz[j] * fl1_fz * pb_xy[j] * fl1_fx + 6.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yyy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyyy[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyyy[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xyyy[j]);

                t_xzz_xyyz[j] = fl_s_0_0 * (0.25 * fl3_fx * pa_z[j] + 0.125 * fl3_fx * pb_z[j] + 0.5 * pa_xz[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pa_zz[j] * pb_z[j] + 0.5 * fl2_fx * pa_z[j] * pb_yy[j] + 0.25 * pa_x[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pb_yyz[j] + 0.5 * pa_xzz[j] * pb_xz[j] * fl1_fx + pa_xz[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yyz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyyz[j] + pa_xzz[j] * pb_xyyz[j]);

                t_xzz_xyyz[j] += fl_r_0_0 * (-0.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 2.0 * fl3_fx * pa_z[j] * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.25 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb * pb_z[j] + fl3_fx * fl1_fz * pb_z[j] + 5.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_x[j] + 2.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_z[j] + 5.0 * fl2_fx * pa_z[j] * fl1_fz * pb_yy[j] - 0.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyz[j] - pa_xzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_xz[j] + 2.5 * fl2_fx * fl1_fz * pb_yyz[j] + 6.0 * pa_xzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 12.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xyy[j] + 6.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yyz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyyz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyyz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xyyz[j]);

                t_xzz_xyzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_y[j] + 0.75 * pa_x[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pa_zz[j] * pb_y[j] + fl2_fx * pa_z[j] * pb_yz[j] + 0.25 * fl2_fx * pb_yzz[j] + 0.5 * pa_xzz[j] * pb_xy[j] * fl1_fx + 2.0 * pa_xz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyzz[j] + pa_xzz[j] * pb_xyzz[j]);

                t_xzz_xyzz[j] += fl_r_0_0 * (3.0 * fl3_fx * fl1_fz * pb_y[j] - 0.25 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 0.5 * fl1_fx * pa_zz[j] * pb_y[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xy[j] + 2.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_y[j] + 10.0 * fl2_fx * pa_z[j] * fl1_fz * pb_yz[j] - 0.5 * pa_x[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yzz[j] - pa_xzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_yzz[j] + 6.0 * pa_xzz[j] * fl1_fz * pb_xy[j] * fl1_fx + 24.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xyz[j] + 6.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyzz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xyzz[j]);

                t_xzz_xzzz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] + 1.125 * fl3_fx * pb_z[j] + 1.5 * pa_xz[j] * fl2_fx * pb_x[j] + 2.25 * pa_x[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_z[j] + 1.5 * fl2_fx * pa_z[j] * pb_zz[j] + 0.25 * fl2_fx * pb_zzz[j] + 1.5 * pa_xzz[j] * pb_xz[j] * fl1_fx + 3.0 * pa_xz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_zzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xzzz[j] + pa_xzz[j] * pb_xzzz[j]);

                t_xzz_xzzz[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * pa_z[j] * fl1_fz + 9.0 * fl3_fx * fl1_fz * pb_z[j] - 0.75 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 3.0 * pa_xz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_zz[j] * pb_z[j] * fl1_fz * fl1_fgb + 15.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_x[j] + 22.5 * pa_x[j] * fl2_fx * fl1_fz * pb_xz[j] + 7.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_z[j] + 15.0 * fl2_fx * pa_z[j] * fl1_fz * pb_zz[j] - 1.5 * pa_x[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_zzz[j] - 3.0 * pa_xzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_zzz[j] + 18.0 * pa_xzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 36.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xzz[j] + 6.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_zzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xzzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xzzz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_xzzz[j]);

                t_xzz_yyyy[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xzz[j] * fl2_fx + 1.5 * pa_x[j] * fl2_fx * pb_yy[j] + 3.0 * pa_xzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_yyyy[j] + pa_xzz[j] * pb_yyyy[j]);

                t_xzz_yyyy[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx + 7.5 * pa_xzz[j] * fl1_fz * fl2_fx - 3.0 * pa_x[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 6.0 * pa_xzz[j] * pb_yy[j] * fl1_fz * fl1_fgb + 15.0 * pa_x[j] * fl1_fz * fl2_fx * pb_yy[j] + 36.0 * pa_xzz[j] * fl1_fz * pb_yy[j] * fl1_fx - pa_x[j] * fl1_fz * fl1_fga * pb_yyyy[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyyy[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_yyyy[j]);

                t_xzz_yyyz[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl2_fx * pb_y[j] + 0.75 * pa_x[j] * fl2_fx * pb_yz[j] + 1.5 * pa_xzz[j] * pb_yz[j] * fl1_fx + pa_xz[j] * fl1_fx * pb_yyy[j] + 0.5 * pa_x[j] * fl1_fx * pb_yyyz[j] + pa_xzz[j] * pb_yyyz[j]);

                t_xzz_yyyz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 15.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_y[j] - 1.5 * pa_x[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 3.0 * pa_xzz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa_x[j] * fl1_fz * fl2_fx * pb_yz[j] + 18.0 * pa_xzz[j] * fl1_fz * pb_yz[j] * fl1_fx + 12.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_yyy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yyyz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyyz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_yyyz[j]);

                t_xzz_yyzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xzz[j] * fl2_fx + pa_xz[j] * fl2_fx * pb_z[j] + 0.75 * pa_x[j] * fl2_fx * pb_yy[j] + 0.25 * pa_x[j] * fl2_fx * pb_zz[j] + 0.5 * pa_xzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xz[j] * fl1_fx * pb_yyz[j] + 0.5 * pa_x[j] * fl1_fx * pb_yyzz[j] + pa_xzz[j] * pb_yyzz[j]);

                t_xzz_yyzz[j] += fl_r_0_0 * (-pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 2.5 * pa_xzz[j] * fl1_fz * fl2_fx + 10.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_z[j] + 7.5 * pa_x[j] * fl2_fx * fl1_fz * pb_yy[j] - 0.5 * pa_x[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - pa_xzz[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_xzz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * pa_x[j] * fl1_fz * fl2_fx * pb_zz[j] + 6.0 * pa_xzz[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_xzz[j] * fl1_fz * fl1_fx * pb_zz[j] + 24.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_yyz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yyzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyzz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_yyzz[j]);

                t_xzz_yzzz[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl2_fx * pb_y[j] + 2.25 * pa_x[j] * fl2_fx * pb_yz[j] + 1.5 * pa_xzz[j] * pb_yz[j] * fl1_fx + 3.0 * pa_xz[j] * fl1_fx * pb_yzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_yzzz[j] + pa_xzz[j] * pb_yzzz[j]);

                t_xzz_yzzz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 15.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_y[j] + 22.5 * pa_x[j] * fl2_fx * fl1_fz * pb_yz[j] - 1.5 * pa_x[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 3.0 * pa_xzz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa_xzz[j] * fl1_fz * pb_yz[j] * fl1_fx + 36.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_yzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yzzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yzzz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_yzzz[j]);

                t_xzz_zzzz[j] = fl_s_0_0 * (1.875 * pa_x[j] * fl3_fx + 0.75 * pa_xzz[j] * fl2_fx + 6.0 * pa_xz[j] * fl2_fx * pb_z[j] + 4.5 * pa_x[j] * fl2_fx * pb_zz[j] + 3.0 * pa_xzz[j] * pb_zz[j] * fl1_fx + 4.0 * pa_xz[j] * fl1_fx * pb_zzz[j] + 0.5 * pa_x[j] * fl1_fx * pb_zzzz[j] + pa_xzz[j] * pb_zzzz[j]);

                t_xzz_zzzz[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa_xz[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 7.5 * pa_xzz[j] * fl1_fz * fl2_fx + 60.0 * pa_xz[j] * fl1_fz * fl2_fx * pb_z[j] + 45.0 * pa_x[j] * fl2_fx * fl1_fz * pb_zz[j] - 3.0 * pa_x[j] * fl1_fx * pb_zz[j] * fl1_fz * fl1_fgb - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_zz[j] * fl1_fx - 6.0 * pa_xzz[j] * pb_zz[j] * fl1_fz * fl1_fgb + 36.0 * pa_xzz[j] * fl1_fz * pb_zz[j] * fl1_fx + 48.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_zzz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_zzzz[j] + 6.0 * pa_x[j] * fl1_fz * fl1_fx * pb_zzzz[j] + 14.0 * pa_xzz[j] * fl1_fz * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_90_100(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (90,100)

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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(19 * idx + 1);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(19 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(19 * idx + 15);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyy_xxxx = primBuffer.data(150 * idx + 90);

            auto t_yyy_xxxy = primBuffer.data(150 * idx + 91);

            auto t_yyy_xxxz = primBuffer.data(150 * idx + 92);

            auto t_yyy_xxyy = primBuffer.data(150 * idx + 93);

            auto t_yyy_xxyz = primBuffer.data(150 * idx + 94);

            auto t_yyy_xxzz = primBuffer.data(150 * idx + 95);

            auto t_yyy_xyyy = primBuffer.data(150 * idx + 96);

            auto t_yyy_xyyz = primBuffer.data(150 * idx + 97);

            auto t_yyy_xyzz = primBuffer.data(150 * idx + 98);

            auto t_yyy_xzzz = primBuffer.data(150 * idx + 99);

            // Batch of Integrals (90,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, \
                                     t_yyy_xxxx, t_yyy_xxxy, t_yyy_xxxz, t_yyy_xxyy, t_yyy_xxyz, t_yyy_xxzz, t_yyy_xyyy, \
                                     t_yyy_xyyz, t_yyy_xyzz, t_yyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yyy_xxxx[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pa_yyy[j] * fl2_fx + 4.5 * pa_y[j] * fl2_fx * pb_xx[j] + 3.0 * pa_yyy[j] * pb_xx[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pb_xxxx[j] + pa_yyy[j] * pb_xxxx[j]);

                t_yyy_xxxx[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl1_fz * fl3_fx + 7.5 * pa_yyy[j] * fl1_fz * fl2_fx - 9.0 * pa_y[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 9.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 6.0 * pa_yyy[j] * pb_xx[j] * fl1_fz * fl1_fgb + 45.0 * pa_y[j] * fl1_fz * fl2_fx * pb_xx[j] + 36.0 * pa_yyy[j] * fl1_fz * pb_xx[j] * fl1_fx - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxxx[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxxx[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xxxx[j]);

                t_yyy_xxxy[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_x[j] + 2.25 * pa_yy[j] * fl2_fx * pb_x[j] + 2.25 * pa_y[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pb_xxx[j] + 1.5 * pa_yyy[j] * pb_xy[j] * fl1_fx + 1.5 * pa_yy[j] * fl1_fx * pb_xxx[j] + 1.5 * pa_y[j] * fl1_fx * pb_xxxy[j] + pa_yyy[j] * pb_xxxy[j]);

                t_yyy_xxxy[j] += fl_r_0_0 * (-2.25 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 4.5 * pa_yy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 9.0 * fl3_fx * fl1_fz * pb_x[j] + 22.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_x[j] - 4.5 * pa_y[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxx[j] - 3.0 * pa_yyy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 22.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pb_xxx[j] + 18.0 * pa_yyy[j] * fl1_fz * pb_xy[j] * fl1_fx + 18.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xxx[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxxy[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxxy[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xxxy[j]);

                t_yyy_xxxz[j] = fl_s_0_0 * (2.25 * pa_y[j] * fl2_fx * pb_xz[j] + 1.5 * pa_yyy[j] * pb_xz[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pb_xxxz[j] + pa_yyy[j] * pb_xxxz[j]);

                t_yyy_xxxz[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 3.0 * pa_yyy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 22.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xz[j] + 18.0 * pa_yyy[j] * fl1_fz * pb_xz[j] * fl1_fx - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxxz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxxz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xxxz[j]);

                t_yyy_xxyy[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * fl3_fx * pb_y[j] + 0.25 * pa_yyy[j] * fl2_fx + 1.5 * pa_yy[j] * fl2_fx * pb_y[j] + 2.25 * pa_y[j] * fl2_fx * pb_xx[j] + 0.75 * pa_y[j] * fl2_fx * pb_yy[j] + 1.5 * fl2_fx * pb_xxy[j] + 0.5 * pa_yyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyy[j] * fl1_fx * pb_yy[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xxy[j] + 1.5 * pa_y[j] * fl1_fx * pb_xxyy[j] + pa_yyy[j] * pb_xxyy[j]);

                t_yyy_xxyy[j] += fl_r_0_0 * (-3.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 1.5 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] + 6.0 * fl3_fx * fl1_fz * pb_y[j] + 2.5 * pa_yyy[j] * fl1_fz * fl2_fx + 15.0 * pa_yy[j] * fl1_fz * fl2_fx * pb_y[j] + 22.5 * pa_y[j] * fl2_fx * fl1_fz * pb_xx[j] - 1.5 * pa_y[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yy[j] - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 1.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xxy[j] - pa_yyy[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_yyy[j] * fl1_fz * fl1_fgb * pb_yy[j] + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_yy[j] + 15.0 * fl2_fx * fl1_fz * pb_xxy[j] + 6.0 * pa_yyy[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_yyy[j] * fl1_fz * fl1_fx * pb_yy[j] + 36.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xxy[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxyy[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxyy[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xxyy[j]);

                t_yyy_xxyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_z[j] + 0.75 * pa_yy[j] * fl2_fx * pb_z[j] + 0.75 * pa_y[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pb_xxz[j] + 0.5 * pa_yyy[j] * fl1_fx * pb_yz[j] + 1.5 * pa_yy[j] * fl1_fx * pb_xxz[j] + 1.5 * pa_y[j] * fl1_fx * pb_xxyz[j] + pa_yyy[j] * pb_xxyz[j]);

                t_yyy_xxyz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 1.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 3.0 * fl3_fx * fl1_fz * pb_z[j] + 7.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_z[j] - 1.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yz[j] - 1.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxz[j] - pa_yyy[j] * fl1_fz * fl1_fgb * pb_yz[j] + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pb_xxz[j] + 6.0 * pa_yyy[j] * fl1_fz * fl1_fx * pb_yz[j] + 18.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xxz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxyz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxyz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xxyz[j]);

                t_yyy_xxzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.25 * pa_yyy[j] * fl2_fx + 0.75 * pa_y[j] * fl2_fx * pb_xx[j] + 0.75 * pa_y[j] * fl2_fx * pb_zz[j] + 0.5 * pa_yyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyy[j] * fl1_fx * pb_zz[j] + 1.5 * pa_y[j] * fl1_fx * pb_xxzz[j] + pa_yyy[j] * pb_xxzz[j]);

                t_yyy_xxzz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl1_fz * fl3_fx + 2.5 * pa_yyy[j] * fl1_fz * fl2_fx - 1.5 * pa_y[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 1.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - pa_yyy[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_yyy[j] * fl1_fz * fl1_fgb * pb_zz[j] + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xx[j] + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_zz[j] + 6.0 * pa_yyy[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_yyy[j] * fl1_fz * fl1_fx * pb_zz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxzz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxzz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xxzz[j]);

                t_yyy_xyyy[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_x[j] + 2.25 * pa_yy[j] * fl2_fx * pb_x[j] + 6.75 * pa_y[j] * fl2_fx * pb_xy[j] + 2.25 * fl2_fx * pb_xyy[j] + 1.5 * pa_yyy[j] * pb_xy[j] * fl1_fx + 4.5 * pa_yy[j] * fl1_fx * pb_xyy[j] + 1.5 * pa_y[j] * fl1_fx * pb_xyyy[j] + pa_yyy[j] * pb_xyyy[j]);

                t_yyy_xyyy[j] += fl_r_0_0 * (15.0 * fl3_fx * fl1_fz * pb_x[j] - 2.25 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 4.5 * pa_yy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 22.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_x[j] + 67.5 * pa_y[j] * fl2_fx * fl1_fz * pb_xy[j] - 4.5 * pa_y[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_xyy[j] - 3.0 * pa_yyy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 22.5 * fl2_fx * fl1_fz * pb_xyy[j] + 18.0 * pa_yyy[j] * fl1_fz * pb_xy[j] * fl1_fx + 54.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xyy[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xyyy[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyyy[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xyyy[j]);

                t_yyy_xyyz[j] = fl_s_0_0 * (2.25 * pa_y[j] * fl2_fx * pb_xz[j] + 1.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_yyy[j] * pb_xz[j] * fl1_fx + 3.0 * pa_yy[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_y[j] * fl1_fx * pb_xyyz[j] + pa_yyy[j] * pb_xyyz[j]);

                t_yyy_xyyz[j] += fl_r_0_0 * (22.5 * pa_y[j] * fl2_fx * fl1_fz * pb_xz[j] - 1.5 * pa_y[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xyz[j] - pa_yyy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 15.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_yyy[j] * fl1_fz * pb_xz[j] * fl1_fx + 36.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xyz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xyyz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyyz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xyyz[j]);

                t_yyy_xyzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_x[j] + 0.75 * pa_yy[j] * fl2_fx * pb_x[j] + 0.75 * pa_y[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pb_xzz[j] + 0.5 * pa_yyy[j] * pb_xy[j] * fl1_fx + 1.5 * pa_yy[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_y[j] * fl1_fx * pb_xyzz[j] + pa_yyy[j] * pb_xyzz[j]);

                t_yyy_xyzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 1.5 * pa_yy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_x[j] - 1.5 * pa_y[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xzz[j] - pa_yyy[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pb_xzz[j] + 6.0 * pa_yyy[j] * fl1_fz * pb_xy[j] * fl1_fx + 18.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xzz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xyzz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyzz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xyzz[j]);

                t_yyy_xzzz[j] = fl_s_0_0 * (2.25 * pa_y[j] * fl2_fx * pb_xz[j] + 1.5 * pa_yyy[j] * pb_xz[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pb_xzzz[j] + pa_yyy[j] * pb_xzzz[j]);

                t_yyy_xzzz[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 3.0 * pa_yyy[j] * pb_xz[j] * fl1_fz * fl1_fgb + 22.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xz[j] + 18.0 * pa_yyy[j] * fl1_fz * pb_xz[j] * fl1_fx - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xzzz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xzzz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_100_110(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (100,110)

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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyy_yyyy = primBuffer.data(150 * idx + 100);

            auto t_yyy_yyyz = primBuffer.data(150 * idx + 101);

            auto t_yyy_yyzz = primBuffer.data(150 * idx + 102);

            auto t_yyy_yzzz = primBuffer.data(150 * idx + 103);

            auto t_yyy_zzzz = primBuffer.data(150 * idx + 104);

            auto t_yyz_xxxx = primBuffer.data(150 * idx + 105);

            auto t_yyz_xxxy = primBuffer.data(150 * idx + 106);

            auto t_yyz_xxxz = primBuffer.data(150 * idx + 107);

            auto t_yyz_xxyy = primBuffer.data(150 * idx + 108);

            auto t_yyz_xxyz = primBuffer.data(150 * idx + 109);

            // Batch of Integrals (100,110)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, r_0_0, s_0_0, t_yyy_yyyy, t_yyy_yyyz, t_yyy_yyzz, t_yyy_yzzz, t_yyy_zzzz, \
                                     t_yyz_xxxx, t_yyz_xxxy, t_yyz_xxxz, t_yyz_xxyy, t_yyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yyy_yyyy[j] = fl_s_0_0 * (5.625 * pa_y[j] * fl3_fx + 7.5 * fl3_fx * pb_y[j] + 0.75 * pa_yyy[j] * fl2_fx + 9.0 * pa_yy[j] * fl2_fx * pb_y[j] + 13.5 * pa_y[j] * fl2_fx * pb_yy[j] + 3.0 * fl2_fx * pb_yyy[j] + 3.0 * pa_yyy[j] * pb_yy[j] * fl1_fx + 6.0 * pa_yy[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_y[j] * fl1_fx * pb_yyyy[j] + pa_yyy[j] * pb_yyyy[j]);

                t_yyy_yyyy[j] += fl_r_0_0 * (-13.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 45.0 * pa_y[j] * fl3_fx * fl1_fz + 60.0 * fl3_fx * fl1_fz * pb_y[j] - 2.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 9.0 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 3.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa_yy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 7.5 * pa_yyy[j] * fl1_fz * fl2_fx + 90.0 * pa_yy[j] * fl1_fz * fl2_fx * pb_y[j] + 135.0 * pa_y[j] * fl2_fx * fl1_fz * pb_yy[j] - 9.0 * pa_y[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 9.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 6.0 * fl1_fx * fl1_fz * fl1_fga * pb_yyy[j] - 6.0 * pa_yyy[j] * pb_yy[j] * fl1_fz * fl1_fgb + 30.0 * fl2_fx * fl1_fz * pb_yyy[j] + 36.0 * pa_yyy[j] * fl1_fz * pb_yy[j] * fl1_fx + 72.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yyy[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yyyy[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyyy[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_yyyy[j]);

                t_yyy_yyyz[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_z[j] + 2.25 * pa_yy[j] * fl2_fx * pb_z[j] + 6.75 * pa_y[j] * fl2_fx * pb_yz[j] + 2.25 * fl2_fx * pb_yyz[j] + 1.5 * pa_yyy[j] * pb_yz[j] * fl1_fx + 4.5 * pa_yy[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_y[j] * fl1_fx * pb_yyyz[j] + pa_yyy[j] * pb_yyyz[j]);

                t_yyy_yyyz[j] += fl_r_0_0 * (15.0 * fl3_fx * fl1_fz * pb_z[j] - 2.25 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 4.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 22.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_z[j] + 67.5 * pa_y[j] * fl2_fx * fl1_fz * pb_yz[j] - 4.5 * pa_y[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyz[j] - 3.0 * pa_yyy[j] * pb_yz[j] * fl1_fz * fl1_fgb + 22.5 * fl2_fx * fl1_fz * pb_yyz[j] + 18.0 * pa_yyy[j] * fl1_fz * pb_yz[j] * fl1_fx + 54.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yyz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yyyz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyyz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_yyyz[j]);

                t_yyy_yyzz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * fl3_fx * pb_y[j] + 0.25 * pa_yyy[j] * fl2_fx + 1.5 * pa_yy[j] * fl2_fx * pb_y[j] + 2.25 * pa_y[j] * fl2_fx * pb_zz[j] + 0.75 * pa_y[j] * fl2_fx * pb_yy[j] + 1.5 * fl2_fx * pb_yzz[j] + 0.5 * pa_yyy[j] * pb_yy[j] * fl1_fx + 0.5 * pa_yyy[j] * fl1_fx * pb_zz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_y[j] * fl1_fx * pb_yyzz[j] + pa_yyy[j] * pb_yyzz[j]);

                t_yyy_yyzz[j] += fl_r_0_0 * (-3.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa_yy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * fl1_fz * pb_y[j] + 2.5 * pa_yyy[j] * fl1_fz * fl2_fx + 15.0 * pa_yy[j] * fl1_fz * fl2_fx * pb_y[j] + 22.5 * pa_y[j] * fl2_fx * fl1_fz * pb_zz[j] - 1.5 * pa_y[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 1.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_yzz[j] - pa_yyy[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_yyy[j] * fl1_fz * fl1_fgb * pb_zz[j] + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_yy[j] + 15.0 * fl2_fx * fl1_fz * pb_yzz[j] + 6.0 * pa_yyy[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_yyy[j] * fl1_fz * fl1_fx * pb_zz[j] + 36.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yzz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yyzz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyzz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_yyzz[j]);

                t_yyy_yzzz[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_z[j] + 2.25 * pa_yy[j] * fl2_fx * pb_z[j] + 2.25 * pa_y[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pb_zzz[j] + 1.5 * pa_yyy[j] * pb_yz[j] * fl1_fx + 1.5 * pa_yy[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_y[j] * fl1_fx * pb_yzzz[j] + pa_yyy[j] * pb_yzzz[j]);

                t_yyy_yzzz[j] += fl_r_0_0 * (-2.25 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 4.5 * pa_yy[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 9.0 * fl3_fx * fl1_fz * pb_z[j] + 22.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_z[j] - 4.5 * pa_y[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_zzz[j] - 3.0 * pa_yyy[j] * pb_yz[j] * fl1_fz * fl1_fgb + 22.5 * pa_y[j] * fl1_fz * fl2_fx * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pb_zzz[j] + 18.0 * pa_yyy[j] * fl1_fz * pb_yz[j] * fl1_fx + 18.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_zzz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yzzz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yzzz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_yzzz[j]);

                t_yyy_zzzz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pa_yyy[j] * fl2_fx + 4.5 * pa_y[j] * fl2_fx * pb_zz[j] + 3.0 * pa_yyy[j] * pb_zz[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pb_zzzz[j] + pa_yyy[j] * pb_zzzz[j]);

                t_yyy_zzzz[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl1_fz * fl3_fx + 7.5 * pa_yyy[j] * fl1_fz * fl2_fx - 9.0 * pa_y[j] * fl1_fx * pb_zz[j] * fl1_fz * fl1_fgb - 9.0 * pa_y[j] * fl1_fz * fl1_fga * pb_zz[j] * fl1_fx - 6.0 * pa_yyy[j] * pb_zz[j] * fl1_fz * fl1_fgb + 45.0 * pa_y[j] * fl1_fz * fl2_fx * pb_zz[j] + 36.0 * pa_yyy[j] * fl1_fz * pb_zz[j] * fl1_fx - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_zzzz[j] + 18.0 * pa_y[j] * fl1_fz * fl1_fx * pb_zzzz[j] + 14.0 * pa_yyy[j] * fl1_fz * pb_zzzz[j]);

                t_yyz_xxxx[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 0.75 * pa_yyz[j] * fl2_fx + 1.5 * fl2_fx * pa_z[j] * pb_xx[j] + 3.0 * pa_yyz[j] * pb_xx[j] * fl1_fx + 0.5 * fl1_fx * pa_z[j] * pb_xxxx[j] + pa_yyz[j] * pb_xxxx[j]);

                t_yyz_xxxx[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 3.0 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] + 7.5 * pa_yyz[j] * fl1_fz * fl2_fx - 3.0 * fl1_fx * pa_z[j] * pb_xx[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_z[j] * pb_xx[j] * fl1_fx - 6.0 * pa_yyz[j] * pb_xx[j] * fl1_fz * fl1_fgb + 15.0 * fl2_fx * fl1_fz * pa_z[j] * pb_xx[j] + 36.0 * pa_yyz[j] * fl1_fz * pb_xx[j] * fl1_fx - fl1_fz * fl1_fga * pa_z[j] * pb_xxxx[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxxx[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xxxx[j]);

                t_yyz_xxxy[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_z[j] * pb_xy[j] + 1.5 * pa_yyz[j] * pb_xy[j] * fl1_fx + pa_yz[j] * fl1_fx * pb_xxx[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxxy[j] + pa_yyz[j] * pb_xxxy[j]);

                t_yyz_xxxy[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 15.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_x[j] - 1.5 * fl1_fx * pa_z[j] * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xy[j] * fl1_fx - 3.0 * pa_yyz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xy[j] + 18.0 * pa_yyz[j] * fl1_fz * pb_xy[j] * fl1_fx + 12.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xxx[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxxy[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxxy[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xxxy[j]);

                t_yyz_xxxz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_x[j] + 0.75 * pa_yy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_z[j] * pb_xz[j] + 0.25 * fl2_fx * pb_xxx[j] + 1.5 * pa_yyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_yy[j] * fl1_fx * pb_xxx[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxxz[j] + pa_yyz[j] * pb_xxxz[j]);

                t_yyz_xxxz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 1.5 * pa_yy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_x[j] - 1.5 * fl1_fx * pa_z[j] * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xz[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xxx[j] - 3.0 * pa_yyz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xz[j] + 2.5 * fl2_fx * fl1_fz * pb_xxx[j] + 18.0 * pa_yyz[j] * fl1_fz * pb_xz[j] * fl1_fx + 6.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xxx[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxxz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxxz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xxxz[j]);

                t_yyz_xxyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 0.25 * pa_yyz[j] * fl2_fx + pa_yz[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pa_z[j] * pb_xx[j] + 0.25 * fl2_fx * pa_z[j] * pb_yy[j] + 0.5 * pa_yyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyz[j] * fl1_fx * pb_yy[j] + 2.0 * pa_yz[j] * fl1_fx * pb_xxy[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxyy[j] + pa_yyz[j] * pb_xxyy[j]);

                t_yyz_xxyy[j] += fl_r_0_0 * (-fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] - 0.25 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] + 2.5 * pa_yyz[j] * fl1_fz * fl2_fx + 10.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_y[j] + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xx[j] - 0.5 * fl1_fx * pa_z[j] * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_yy[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xx[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_yy[j] - pa_yyz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_yyz[j] * fl1_fz * fl1_fgb * pb_yy[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yy[j] + 6.0 * pa_yyz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_yyz[j] * fl1_fz * fl1_fx * pb_yy[j] + 24.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xxy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxyy[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxyy[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xxyy[j]);

                t_yyz_xxyz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl3_fx + 0.125 * fl3_fx * pb_y[j] + 0.25 * pa_yy[j] * fl2_fx * pb_y[j] + 0.5 * pa_yz[j] * fl2_fx * pb_z[j] + 0.5 * pa_y[j] * fl2_fx * pb_xx[j] + 0.25 * fl2_fx * pa_z[j] * pb_yz[j] + 0.25 * fl2_fx * pb_xxy[j] + 0.5 * pa_yyz[j] * fl1_fx * pb_yz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xxy[j] + pa_yz[j] * fl1_fx * pb_xxz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxyz[j] + pa_yyz[j] * pb_xxyz[j]);

                t_yyz_xxyz[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.25 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - 0.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + fl3_fx * fl1_fz * pb_y[j] + 2.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_y[j] + 5.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_z[j] + 5.0 * pa_y[j] * fl2_fx * fl1_fz * pb_xx[j] - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_yz[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_yz[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xxy[j] - pa_yyz[j] * fl1_fz * fl1_fgb * pb_yz[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yz[j] + 2.5 * fl2_fx * fl1_fz * pb_xxy[j] + 6.0 * pa_yyz[j] * fl1_fz * fl1_fx * pb_yz[j] + 6.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xxz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxyz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxyz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_110_120(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (110,120)

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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyz = paDistances.data(19 * idx + 16);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyz_xxzz = primBuffer.data(150 * idx + 110);

            auto t_yyz_xyyy = primBuffer.data(150 * idx + 111);

            auto t_yyz_xyyz = primBuffer.data(150 * idx + 112);

            auto t_yyz_xyzz = primBuffer.data(150 * idx + 113);

            auto t_yyz_xzzz = primBuffer.data(150 * idx + 114);

            auto t_yyz_yyyy = primBuffer.data(150 * idx + 115);

            auto t_yyz_yyyz = primBuffer.data(150 * idx + 116);

            auto t_yyz_yyzz = primBuffer.data(150 * idx + 117);

            auto t_yyz_yzzz = primBuffer.data(150 * idx + 118);

            auto t_yyz_zzzz = primBuffer.data(150 * idx + 119);

            // Batch of Integrals (110,120)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, r_0_0, s_0_0, t_yyz_xxzz, t_yyz_xyyy, t_yyz_xyyz, t_yyz_xyzz, t_yyz_xzzz, \
                                     t_yyz_yyyy, t_yyz_yyyz, t_yyz_yyzz, t_yyz_yzzz, t_yyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yyz_xxzz[j] = fl_s_0_0 * (0.125 * fl3_fx * pa_z[j] + 0.25 * fl3_fx * pb_z[j] + 0.25 * pa_yyz[j] * fl2_fx + 0.5 * pa_yy[j] * fl2_fx * pb_z[j] + 0.25 * fl2_fx * pa_z[j] * pb_xx[j] + 0.25 * fl2_fx * pa_z[j] * pb_zz[j] + 0.5 * fl2_fx * pb_xxz[j] + 0.5 * pa_yyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyz[j] * fl1_fx * pb_zz[j] + pa_yy[j] * fl1_fx * pb_xxz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxzz[j] + pa_yyz[j] * pb_xxzz[j]);

                t_yyz_xxzz[j] += fl_r_0_0 * (-0.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb - 0.5 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.25 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 0.5 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + fl3_fx * fl1_fz * pa_z[j] + 2.0 * fl3_fx * fl1_fz * pb_z[j] + 2.5 * pa_yyz[j] * fl1_fz * fl2_fx + 5.0 * pa_yy[j] * fl1_fz * fl2_fx * pb_z[j] - 0.5 * fl1_fx * pa_z[j] * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xx[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_zz[j] - fl1_fz * fl1_fga * fl1_fx * pb_xxz[j] - pa_yyz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_yyz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xx[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_zz[j] + 5.0 * fl2_fx * fl1_fz * pb_xxz[j] + 6.0 * pa_yyz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_yyz[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xxz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxzz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xxzz[j]);

                t_yyz_xyyy[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl2_fx * pb_x[j] + 2.25 * fl2_fx * pa_z[j] * pb_xy[j] + 1.5 * pa_yyz[j] * pb_xy[j] * fl1_fx + 3.0 * pa_yz[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyyy[j] + pa_yyz[j] * pb_xyyy[j]);

                t_yyz_xyyy[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 15.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_x[j] + 22.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xy[j] - 1.5 * fl1_fx * pa_z[j] * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xy[j] * fl1_fx - 3.0 * pa_yyz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa_yyz[j] * fl1_fz * pb_xy[j] * fl1_fx + 36.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xyy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyyy[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyyy[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xyyy[j]);

                t_yyz_xyyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_x[j] + 0.25 * pa_yy[j] * fl2_fx * pb_x[j] + pa_y[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_z[j] * pb_xz[j] + 0.25 * fl2_fx * pb_xyy[j] + 0.5 * pa_yyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_yy[j] * fl1_fx * pb_xyy[j] + 2.0 * pa_yz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyyz[j] + pa_yyz[j] * pb_xyyz[j]);

                t_yyz_xyyz[j] += fl_r_0_0 * (3.0 * fl3_fx * fl1_fz * pb_x[j] - 0.25 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 0.5 * pa_yy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 2.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_x[j] + 10.0 * pa_y[j] * fl2_fx * fl1_fz * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xz[j] - 0.5 * fl1_fx * pa_z[j] * pb_xz[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xz[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xyy[j] - pa_yyz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_xyy[j] + 6.0 * pa_yyz[j] * fl1_fz * pb_xz[j] * fl1_fx + 6.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xyy[j] + 24.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xyz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyyz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyyz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xyyz[j]);

                t_yyz_xyzz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx * pb_x[j] + pa_y[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pa_z[j] * pb_xy[j] + 0.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_yyz[j] * pb_xy[j] * fl1_fx + pa_yy[j] * fl1_fx * pb_xyz[j] + pa_yz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyzz[j] + pa_yyz[j] * pb_xyzz[j]);

                t_yyz_xyzz[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 5.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_x[j] + 10.0 * pa_y[j] * fl2_fx * fl1_fz * pb_xz[j] - 0.5 * fl1_fx * pa_z[j] * pb_xy[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xy[j] * fl1_fx - fl1_fz * fl1_fga * fl1_fx * pb_xyz[j] - pa_yyz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xy[j] + 5.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_yyz[j] * fl1_fz * pb_xy[j] * fl1_fx + 12.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyzz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xyzz[j]);

                t_yyz_xzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_x[j] + 0.75 * pa_yy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_z[j] * pb_xz[j] + 0.75 * fl2_fx * pb_xzz[j] + 1.5 * pa_yyz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_yy[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xzzz[j] + pa_yyz[j] * pb_xzzz[j]);

                t_yyz_xzzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_x[j] - 1.5 * pa_yy[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_x[j] - 1.5 * fl1_fx * pa_z[j] * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_xz[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_xzz[j] - 3.0 * pa_yyz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pb_xzz[j] + 18.0 * pa_yyz[j] * fl1_fz * pb_xz[j] * fl1_fx + 18.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xzzz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_xzzz[j]);

                t_yyz_yyyy[j] = fl_s_0_0 * (1.875 * fl3_fx * pa_z[j] + 0.75 * pa_yyz[j] * fl2_fx + 6.0 * pa_yz[j] * fl2_fx * pb_y[j] + 4.5 * fl2_fx * pa_z[j] * pb_yy[j] + 3.0 * pa_yyz[j] * pb_yy[j] * fl1_fx + 4.0 * pa_yz[j] * fl1_fx * pb_yyy[j] + 0.5 * fl1_fx * pa_z[j] * pb_yyyy[j] + pa_yyz[j] * pb_yyyy[j]);

                t_yyz_yyyy[j] += fl_r_0_0 * (-4.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 15.0 * fl3_fx * fl1_fz * pa_z[j] - 0.75 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 3.0 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa_yz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 7.5 * pa_yyz[j] * fl1_fz * fl2_fx + 60.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_y[j] + 45.0 * fl2_fx * fl1_fz * pa_z[j] * pb_yy[j] - 3.0 * fl1_fx * pa_z[j] * pb_yy[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_z[j] * pb_yy[j] * fl1_fx - 6.0 * pa_yyz[j] * pb_yy[j] * fl1_fz * fl1_fgb + 36.0 * pa_yyz[j] * fl1_fz * pb_yy[j] * fl1_fx + 48.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_yyy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yyyy[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyyy[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_yyyy[j]);

                t_yyz_yyyz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 1.125 * fl3_fx * pb_y[j] + 0.75 * pa_yy[j] * fl2_fx * pb_y[j] + 1.5 * pa_yz[j] * fl2_fx * pb_z[j] + 1.5 * pa_y[j] * fl2_fx * pb_yy[j] + 2.25 * fl2_fx * pa_z[j] * pb_yz[j] + 0.25 * fl2_fx * pb_yyy[j] + 1.5 * pa_yyz[j] * pb_yz[j] * fl1_fx + 0.5 * pa_yy[j] * fl1_fx * pb_yyy[j] + 3.0 * pa_yz[j] * fl1_fx * pb_yyz[j] + 0.5 * fl1_fx * pa_z[j] * pb_yyyz[j] + pa_yyz[j] * pb_yyyz[j]);

                t_yyz_yyyz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl3_fx * fl1_fz + 9.0 * fl3_fx * fl1_fz * pb_y[j] - 0.75 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - 1.5 * pa_yy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 3.0 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 7.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_y[j] + 15.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_z[j] + 15.0 * pa_y[j] * fl2_fx * fl1_fz * pb_yy[j] + 22.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yz[j] - 1.5 * fl1_fx * pa_z[j] * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_yz[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_yyy[j] - 3.0 * pa_yyz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_yyy[j] + 18.0 * pa_yyz[j] * fl1_fz * pb_yz[j] * fl1_fx + 6.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yyy[j] + 36.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_yyz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yyyz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyyz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_yyyz[j]);

                t_yyz_yyzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 0.75 * fl3_fx * pb_z[j] + 0.25 * pa_yyz[j] * fl2_fx + 0.5 * pa_yy[j] * fl2_fx * pb_z[j] + pa_yz[j] * fl2_fx * pb_y[j] + 2.0 * pa_y[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pa_z[j] * pb_zz[j] + 0.25 * fl2_fx * pa_z[j] * pb_yy[j] + 0.5 * fl2_fx * pb_yyz[j] + 0.5 * pa_yyz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_yyz[j] * fl1_fx * pb_zz[j] + pa_yy[j] * fl1_fx * pb_yyz[j] + 2.0 * pa_yz[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_yyzz[j] + pa_yyz[j] * pb_yyzz[j]);

                t_yyz_yyzz[j] += fl_r_0_0 * (-fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] + 6.0 * fl3_fx * fl1_fz * pb_z[j] - 0.5 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.25 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 0.5 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 2.0 * pa_yz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 2.5 * pa_yyz[j] * fl1_fz * fl2_fx + 5.0 * pa_yy[j] * fl1_fz * fl2_fx * pb_z[j] + 10.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_y[j] + 20.0 * pa_y[j] * fl2_fx * fl1_fz * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_zz[j] - 0.5 * fl1_fx * pa_z[j] * pb_yy[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_yy[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_zz[j] - fl1_fz * fl1_fga * fl1_fx * pb_yyz[j] - pa_yyz[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_yyz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yy[j] + 5.0 * fl2_fx * fl1_fz * pb_yyz[j] + 6.0 * pa_yyz[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_yyz[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yyz[j] + 24.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_yzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yyzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyzz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_yyzz[j]);

                t_yyz_yzzz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 0.375 * fl3_fx * pb_y[j] + 0.75 * pa_yy[j] * fl2_fx * pb_y[j] + 1.5 * pa_yz[j] * fl2_fx * pb_z[j] + 1.5 * pa_y[j] * fl2_fx * pb_zz[j] + 0.75 * fl2_fx * pa_z[j] * pb_yz[j] + 0.75 * fl2_fx * pb_yzz[j] + 1.5 * pa_yyz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_yy[j] * fl1_fx * pb_yzz[j] + pa_yz[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_yzzz[j] + pa_yyz[j] * pb_yzzz[j]);

                t_yyz_yzzz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx * pb_y[j] - 1.5 * pa_yy[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 3.0 * pa_yz[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * pa_yy[j] * fl1_fz * fl2_fx * pb_y[j] + 15.0 * pa_yz[j] * fl2_fx * fl1_fz * pb_z[j] + 15.0 * pa_y[j] * fl2_fx * fl1_fz * pb_zz[j] - 1.5 * fl1_fx * pa_z[j] * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_yz[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_yzz[j] - 3.0 * pa_yyz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 7.5 * fl2_fx * fl1_fz * pa_z[j] * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pb_yzz[j] + 18.0 * pa_yyz[j] * fl1_fz * pb_yz[j] * fl1_fx + 18.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yzz[j] + 12.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_zzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yzzz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_yzzz[j]);

                t_yyz_zzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] + 1.5 * fl3_fx * pb_z[j] + 0.75 * pa_yyz[j] * fl2_fx + 3.0 * pa_yy[j] * fl2_fx * pb_z[j] + 1.5 * fl2_fx * pa_z[j] * pb_zz[j] + fl2_fx * pb_zzz[j] + 3.0 * pa_yyz[j] * pb_zz[j] * fl1_fx + 2.0 * pa_yy[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_z[j] * pb_zzzz[j] + pa_yyz[j] * pb_zzzz[j]);

                t_yyz_zzzz[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb - 3.0 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * pa_z[j] * fl2_fx - 3.0 * fl1_fz * fl1_fga * fl2_fx * pb_z[j] - 3.0 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa_yy[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pa_z[j] + 12.0 * fl3_fx * fl1_fz * pb_z[j] + 7.5 * pa_yyz[j] * fl1_fz * fl2_fx + 30.0 * pa_yy[j] * fl1_fz * fl2_fx * pb_z[j] - 3.0 * fl1_fx * pa_z[j] * pb_zz[j] * fl1_fz * fl1_fgb - 3.0 * fl1_fz * fl1_fga * pa_z[j] * pb_zz[j] * fl1_fx - 2.0 * fl1_fz * fl1_fga * fl1_fx * pb_zzz[j] - 6.0 * pa_yyz[j] * pb_zz[j] * fl1_fz * fl1_fgb + 15.0 * fl2_fx * fl1_fz * pa_z[j] * pb_zz[j] + 10.0 * fl2_fx * fl1_fz * pb_zzz[j] + 36.0 * pa_yyz[j] * fl1_fz * pb_zz[j] * fl1_fx + 24.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_zzz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_zzzz[j] + 6.0 * fl1_fx * fl1_fz * pa_z[j] * pb_zzzz[j] + 14.0 * pa_yyz[j] * fl1_fz * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_120_130(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (120,130)

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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yzz = paDistances.data(19 * idx + 17);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzz_xxxx = primBuffer.data(150 * idx + 120);

            auto t_yzz_xxxy = primBuffer.data(150 * idx + 121);

            auto t_yzz_xxxz = primBuffer.data(150 * idx + 122);

            auto t_yzz_xxyy = primBuffer.data(150 * idx + 123);

            auto t_yzz_xxyz = primBuffer.data(150 * idx + 124);

            auto t_yzz_xxzz = primBuffer.data(150 * idx + 125);

            auto t_yzz_xyyy = primBuffer.data(150 * idx + 126);

            auto t_yzz_xyyz = primBuffer.data(150 * idx + 127);

            auto t_yzz_xyzz = primBuffer.data(150 * idx + 128);

            auto t_yzz_xzzz = primBuffer.data(150 * idx + 129);

            // Batch of Integrals (120,130)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, \
                                     s_0_0, t_yzz_xxxx, t_yzz_xxxy, t_yzz_xxxz, t_yzz_xxyy, t_yzz_xxyz, t_yzz_xxzz, \
                                     t_yzz_xyyy, t_yzz_xyyz, t_yzz_xyzz, t_yzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yzz_xxxx[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa_yzz[j] * fl2_fx + 1.5 * pa_y[j] * fl2_fx * pb_xx[j] + 3.0 * pa_yzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pb_xxxx[j] + pa_yzz[j] * pb_xxxx[j]);

                t_yzz_xxxx[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl1_fz * fl3_fx + 7.5 * pa_yzz[j] * fl1_fz * fl2_fx - 3.0 * pa_y[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 6.0 * pa_yzz[j] * pb_xx[j] * fl1_fz * fl1_fgb + 15.0 * pa_y[j] * fl1_fz * fl2_fx * pb_xx[j] + 36.0 * pa_yzz[j] * fl1_fz * pb_xx[j] * fl1_fx - pa_y[j] * fl1_fz * fl1_fga * pb_xxxx[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxxx[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xxxx[j]);

                t_yzz_xxxy[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_x[j] + 0.75 * fl2_fx * pa_zz[j] * pb_x[j] + 0.75 * pa_y[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pb_xxx[j] + 1.5 * pa_yzz[j] * pb_xy[j] * fl1_fx + 0.5 * fl1_fx * pa_zz[j] * pb_xxx[j] + 0.5 * pa_y[j] * fl1_fx * pb_xxxy[j] + pa_yzz[j] * pb_xxxy[j]);

                t_yzz_xxxy[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 1.5 * fl1_fx * pa_zz[j] * pb_x[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_x[j] - 1.5 * pa_y[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxx[j] - 3.0 * pa_yzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xy[j] + 2.5 * fl2_fx * fl1_fz * pb_xxx[j] + 18.0 * pa_yzz[j] * fl1_fz * pb_xy[j] * fl1_fx + 6.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xxx[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xxxy[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxxy[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xxxy[j]);

                t_yzz_xxxz[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl2_fx * pb_x[j] + 0.75 * pa_y[j] * fl2_fx * pb_xz[j] + 1.5 * pa_yzz[j] * pb_xz[j] * fl1_fx + pa_yz[j] * fl1_fx * pb_xxx[j] + 0.5 * pa_y[j] * fl1_fx * pb_xxxz[j] + pa_yzz[j] * pb_xxxz[j]);

                t_yzz_xxxz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 15.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_x[j] - 1.5 * pa_y[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 3.0 * pa_yzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xz[j] + 18.0 * pa_yzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 12.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xxx[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xxxz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxxz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xxxz[j]);

                t_yzz_xxyy[j] = fl_s_0_0 * (0.125 * pa_y[j] * fl3_fx + 0.25 * fl3_fx * pb_y[j] + 0.25 * pa_yzz[j] * fl2_fx + 0.5 * fl2_fx * pa_zz[j] * pb_y[j] + 0.25 * pa_y[j] * fl2_fx * pb_xx[j] + 0.25 * pa_y[j] * fl2_fx * pb_yy[j] + 0.5 * fl2_fx * pb_xxy[j] + 0.5 * pa_yzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yzz[j] * fl1_fx * pb_yy[j] + fl1_fx * pa_zz[j] * pb_xxy[j] + 0.5 * pa_y[j] * fl1_fx * pb_xxyy[j] + pa_yzz[j] * pb_xxyy[j]);

                t_yzz_xxyy[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + pa_y[j] * fl1_fz * fl3_fx - fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb * pb_y[j] + 2.0 * fl3_fx * fl1_fz * pb_y[j] + 2.5 * pa_yzz[j] * fl1_fz * fl2_fx + 5.0 * fl2_fx * pa_zz[j] * fl1_fz * pb_y[j] - 0.5 * pa_y[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yy[j] - 0.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - fl1_fx * fl1_fz * fl1_fga * pb_xxy[j] - pa_yzz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_yzz[j] * fl1_fz * fl1_fgb * pb_yy[j] + 2.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xx[j] + 2.5 * pa_y[j] * fl1_fz * fl2_fx * pb_yy[j] + 5.0 * fl2_fx * fl1_fz * pb_xxy[j] + 6.0 * pa_yzz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_yzz[j] * fl1_fz * fl1_fx * pb_yy[j] + 12.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xxy[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xxyy[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxyy[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xxyy[j]);

                t_yzz_xxyz[j] = fl_s_0_0 * (0.25 * fl3_fx * pa_z[j] + 0.125 * fl3_fx * pb_z[j] + 0.5 * pa_yz[j] * fl2_fx * pb_y[j] + 0.25 * fl2_fx * pa_zz[j] * pb_z[j] + 0.5 * fl2_fx * pa_z[j] * pb_xx[j] + 0.25 * pa_y[j] * fl2_fx * pb_yz[j] + 0.25 * fl2_fx * pb_xxz[j] + 0.5 * pa_yzz[j] * fl1_fx * pb_yz[j] + pa_yz[j] * fl1_fx * pb_xxy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxz[j] + 0.5 * pa_y[j] * fl1_fx * pb_xxyz[j] + pa_yzz[j] * pb_xxyz[j]);

                t_yzz_xxyz[j] += fl_r_0_0 * (-0.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 2.0 * fl3_fx * pa_z[j] * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.25 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb * pb_z[j] + fl3_fx * fl1_fz * pb_z[j] + 5.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_y[j] + 2.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_z[j] + 5.0 * fl2_fx * pa_z[j] * fl1_fz * pb_xx[j] - 0.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yz[j] - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxz[j] - pa_yzz[j] * fl1_fz * fl1_fgb * pb_yz[j] + 2.5 * pa_y[j] * fl1_fz * fl2_fx * pb_yz[j] + 2.5 * fl2_fx * fl1_fz * pb_xxz[j] + 6.0 * pa_yzz[j] * fl1_fz * fl1_fx * pb_yz[j] + 12.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xxy[j] + 6.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xxz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xxyz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxyz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xxyz[j]);

                t_yzz_xxzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.25 * pa_yzz[j] * fl2_fx + pa_yz[j] * fl2_fx * pb_z[j] + 0.75 * pa_y[j] * fl2_fx * pb_xx[j] + 0.25 * pa_y[j] * fl2_fx * pb_zz[j] + 0.5 * pa_yzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_yz[j] * fl1_fx * pb_xxz[j] + 0.5 * pa_y[j] * fl1_fx * pb_xxzz[j] + pa_yzz[j] * pb_xxzz[j]);

                t_yzz_xxzz[j] += fl_r_0_0 * (-pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 2.5 * pa_yzz[j] * fl1_fz * fl2_fx + 10.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_z[j] + 7.5 * pa_y[j] * fl2_fx * fl1_fz * pb_xx[j] - 0.5 * pa_y[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - pa_yzz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_yzz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * pa_y[j] * fl1_fz * fl2_fx * pb_zz[j] + 6.0 * pa_yzz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_yzz[j] * fl1_fz * fl1_fx * pb_zz[j] + 24.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xxz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xxzz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxzz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xxzz[j]);

                t_yzz_xyyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_x[j] + 0.75 * fl2_fx * pa_zz[j] * pb_x[j] + 0.75 * pa_y[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pb_xyy[j] + 1.5 * pa_yzz[j] * pb_xy[j] * fl1_fx + 1.5 * fl1_fx * pa_zz[j] * pb_xyy[j] + 0.5 * pa_y[j] * fl1_fx * pb_xyyy[j] + pa_yzz[j] * pb_xyyy[j]);

                t_yzz_xyyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 1.5 * fl1_fx * pa_zz[j] * pb_x[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_x[j] - 1.5 * pa_y[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xyy[j] - 3.0 * pa_yzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xy[j] + 7.5 * fl2_fx * fl1_fz * pb_xyy[j] + 18.0 * pa_yzz[j] * fl1_fz * pb_xy[j] * fl1_fx + 18.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xyy[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xyyy[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyyy[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xyyy[j]);

                t_yzz_xyyz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx * pb_x[j] + fl2_fx * pa_z[j] * pb_xy[j] + 0.25 * pa_y[j] * fl2_fx * pb_xz[j] + 0.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_yzz[j] * pb_xz[j] * fl1_fx + pa_yz[j] * fl1_fx * pb_xyy[j] + fl1_fx * pa_zz[j] * pb_xyz[j] + 0.5 * pa_y[j] * fl1_fx * pb_xyyz[j] + pa_yzz[j] * pb_xyyz[j]);

                t_yzz_xyyz[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 5.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_x[j] + 10.0 * fl2_fx * pa_z[j] * fl1_fz * pb_xy[j] - 0.5 * pa_y[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - fl1_fx * fl1_fz * fl1_fga * pb_xyz[j] - pa_yzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 2.5 * pa_y[j] * fl1_fz * fl2_fx * pb_xz[j] + 5.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_yzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 12.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xyz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xyyz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyyz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xyyz[j]);

                t_yzz_xyzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_x[j] + 0.75 * pa_y[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pa_zz[j] * pb_x[j] + fl2_fx * pa_z[j] * pb_xz[j] + 0.25 * fl2_fx * pb_xzz[j] + 0.5 * pa_yzz[j] * pb_xy[j] * fl1_fx + 2.0 * pa_yz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xzz[j] + 0.5 * pa_y[j] * fl1_fx * pb_xyzz[j] + pa_yzz[j] * pb_xyzz[j]);

                t_yzz_xyzz[j] += fl_r_0_0 * (3.0 * fl3_fx * fl1_fz * pb_x[j] - 0.25 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 0.5 * fl1_fx * pa_zz[j] * pb_x[j] * fl1_fz * fl1_fgb + 7.5 * pa_y[j] * fl2_fx * fl1_fz * pb_xy[j] + 2.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_x[j] + 10.0 * fl2_fx * pa_z[j] * fl1_fz * pb_xz[j] - 0.5 * pa_y[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_xzz[j] - pa_yzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_xzz[j] + 6.0 * pa_yzz[j] * fl1_fz * pb_xy[j] * fl1_fx + 24.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xyz[j] + 6.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xzz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xyzz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyzz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xyzz[j]);

                t_yzz_xzzz[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl2_fx * pb_x[j] + 2.25 * pa_y[j] * fl2_fx * pb_xz[j] + 1.5 * pa_yzz[j] * pb_xz[j] * fl1_fx + 3.0 * pa_yz[j] * fl1_fx * pb_xzz[j] + 0.5 * pa_y[j] * fl1_fx * pb_xzzz[j] + pa_yzz[j] * pb_xzzz[j]);

                t_yzz_xzzz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 15.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_x[j] + 22.5 * pa_y[j] * fl2_fx * fl1_fz * pb_xz[j] - 1.5 * pa_y[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 3.0 * pa_yzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa_yzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 36.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xzz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xzzz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xzzz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_130_140(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (130,140)

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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

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

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzz_yyyy = primBuffer.data(150 * idx + 130);

            auto t_yzz_yyyz = primBuffer.data(150 * idx + 131);

            auto t_yzz_yyzz = primBuffer.data(150 * idx + 132);

            auto t_yzz_yzzz = primBuffer.data(150 * idx + 133);

            auto t_yzz_zzzz = primBuffer.data(150 * idx + 134);

            auto t_zzz_xxxx = primBuffer.data(150 * idx + 135);

            auto t_zzz_xxxy = primBuffer.data(150 * idx + 136);

            auto t_zzz_xxxz = primBuffer.data(150 * idx + 137);

            auto t_zzz_xxyy = primBuffer.data(150 * idx + 138);

            auto t_zzz_xxyz = primBuffer.data(150 * idx + 139);

            // Batch of Integrals (130,140)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     r_0_0, s_0_0, t_yzz_yyyy, t_yzz_yyyz, t_yzz_yyzz, t_yzz_yzzz, t_yzz_zzzz, \
                                     t_zzz_xxxx, t_zzz_xxxy, t_zzz_xxxz, t_zzz_xxyy, t_zzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_yzz_yyyy[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 1.5 * fl3_fx * pb_y[j] + 0.75 * pa_yzz[j] * fl2_fx + 3.0 * fl2_fx * pa_zz[j] * pb_y[j] + 1.5 * pa_y[j] * fl2_fx * pb_yy[j] + fl2_fx * pb_yyy[j] + 3.0 * pa_yzz[j] * pb_yy[j] * fl1_fx + 2.0 * fl1_fx * pa_zz[j] * pb_yyy[j] + 0.5 * pa_y[j] * fl1_fx * pb_yyyy[j] + pa_yzz[j] * pb_yyyy[j]);

                t_yzz_yyyy[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 3.0 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 3.0 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl1_fz * fl3_fx - 6.0 * fl1_fx * pa_zz[j] * pb_y[j] * fl1_fz * fl1_fgb + 12.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * pa_yzz[j] * fl1_fz * fl2_fx + 30.0 * fl2_fx * pa_zz[j] * fl1_fz * pb_y[j] - 3.0 * pa_y[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 2.0 * fl1_fx * fl1_fz * fl1_fga * pb_yyy[j] - 6.0 * pa_yzz[j] * pb_yy[j] * fl1_fz * fl1_fgb + 15.0 * pa_y[j] * fl1_fz * fl2_fx * pb_yy[j] + 10.0 * fl2_fx * fl1_fz * pb_yyy[j] + 36.0 * pa_yzz[j] * fl1_fz * pb_yy[j] * fl1_fx + 24.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yyy[j] - pa_y[j] * fl1_fz * fl1_fga * pb_yyyy[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyyy[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_yyyy[j]);

                t_yzz_yyyz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] + 0.375 * fl3_fx * pb_z[j] + 1.5 * pa_yz[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pa_zz[j] * pb_z[j] + 1.5 * fl2_fx * pa_z[j] * pb_yy[j] + 0.75 * pa_y[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pb_yyz[j] + 1.5 * pa_yzz[j] * pb_yz[j] * fl1_fx + pa_yz[j] * fl1_fx * pb_yyy[j] + 1.5 * fl1_fx * pa_zz[j] * pb_yyz[j] + 0.5 * pa_y[j] * fl1_fx * pb_yyyz[j] + pa_yzz[j] * pb_yyyz[j]);

                t_yzz_yyyz[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * pa_z[j] * fl1_fz - 0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 3.0 * pa_yz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb * pb_z[j] + 3.0 * fl3_fx * fl1_fz * pb_z[j] + 15.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_y[j] + 7.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_z[j] + 15.0 * fl2_fx * pa_z[j] * fl1_fz * pb_yy[j] - 1.5 * pa_y[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyz[j] - 3.0 * pa_yzz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa_y[j] * fl1_fz * fl2_fx * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pb_yyz[j] + 18.0 * pa_yzz[j] * fl1_fz * pb_yz[j] * fl1_fx + 12.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_yyy[j] + 18.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yyz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_yyyz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyyz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_yyyz[j]);

                t_yzz_yyzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * fl3_fx * pb_y[j] + 0.25 * pa_yzz[j] * fl2_fx + pa_yz[j] * fl2_fx * pb_z[j] + 0.75 * pa_y[j] * fl2_fx * pb_yy[j] + 0.5 * fl2_fx * pa_zz[j] * pb_y[j] + 2.0 * fl2_fx * pa_z[j] * pb_yz[j] + 0.25 * pa_y[j] * fl2_fx * pb_zz[j] + 0.5 * fl2_fx * pb_yzz[j] + 0.5 * pa_yzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_yzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_yz[j] * fl1_fx * pb_yyz[j] + fl1_fx * pa_zz[j] * pb_yzz[j] + 0.5 * pa_y[j] * fl1_fx * pb_yyzz[j] + pa_yzz[j] * pb_yyzz[j]);

                t_yzz_yyzz[j] += fl_r_0_0 * (-pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 6.0 * fl3_fx * fl1_fz * pb_y[j] - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.5 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - fl1_fx * pa_zz[j] * pb_y[j] * fl1_fz * fl1_fgb + 2.5 * pa_yzz[j] * fl1_fz * fl2_fx + 10.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_z[j] + 7.5 * pa_y[j] * fl2_fx * fl1_fz * pb_yy[j] + 5.0 * fl2_fx * pa_zz[j] * fl1_fz * pb_y[j] + 20.0 * fl2_fx * pa_z[j] * fl1_fz * pb_yz[j] - 0.5 * pa_y[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 0.5 * pa_y[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - fl1_fx * fl1_fz * fl1_fga * pb_yzz[j] - pa_yzz[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_yzz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 2.5 * pa_y[j] * fl1_fz * fl2_fx * pb_zz[j] + 5.0 * fl2_fx * fl1_fz * pb_yzz[j] + 6.0 * pa_yzz[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_yzz[j] * fl1_fz * fl1_fx * pb_zz[j] + 24.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_yyz[j] + 12.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yzz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_yyzz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyzz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_yyzz[j]);

                t_yzz_yzzz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] + 1.125 * fl3_fx * pb_z[j] + 1.5 * pa_yz[j] * fl2_fx * pb_y[j] + 2.25 * pa_y[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_z[j] + 1.5 * fl2_fx * pa_z[j] * pb_zz[j] + 0.25 * fl2_fx * pb_zzz[j] + 1.5 * pa_yzz[j] * pb_yz[j] * fl1_fx + 3.0 * pa_yz[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_zzz[j] + 0.5 * pa_y[j] * fl1_fx * pb_yzzz[j] + pa_yzz[j] * pb_yzzz[j]);

                t_yzz_yzzz[j] += fl_r_0_0 * (-1.5 * fl2_fx * pa_z[j] * fl1_fz * fl1_fgb + 6.0 * fl3_fx * pa_z[j] * fl1_fz + 9.0 * fl3_fx * fl1_fz * pb_z[j] - 0.75 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 3.0 * pa_yz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fx * pa_zz[j] * pb_z[j] * fl1_fz * fl1_fgb + 15.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_y[j] + 22.5 * pa_y[j] * fl2_fx * fl1_fz * pb_yz[j] + 7.5 * fl2_fx * pa_zz[j] * fl1_fz * pb_z[j] + 15.0 * fl2_fx * pa_z[j] * fl1_fz * pb_zz[j] - 1.5 * pa_y[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_zzz[j] - 3.0 * pa_yzz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 2.5 * fl2_fx * fl1_fz * pb_zzz[j] + 18.0 * pa_yzz[j] * fl1_fz * pb_yz[j] * fl1_fx + 36.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_yzz[j] + 6.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_zzz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_yzzz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yzzz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_yzzz[j]);

                t_yzz_zzzz[j] = fl_s_0_0 * (1.875 * pa_y[j] * fl3_fx + 0.75 * pa_yzz[j] * fl2_fx + 6.0 * pa_yz[j] * fl2_fx * pb_z[j] + 4.5 * pa_y[j] * fl2_fx * pb_zz[j] + 3.0 * pa_yzz[j] * pb_zz[j] * fl1_fx + 4.0 * pa_yz[j] * fl1_fx * pb_zzz[j] + 0.5 * pa_y[j] * fl1_fx * pb_zzzz[j] + pa_yzz[j] * pb_zzzz[j]);

                t_yzz_zzzz[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa_yz[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 7.5 * pa_yzz[j] * fl1_fz * fl2_fx + 60.0 * pa_yz[j] * fl1_fz * fl2_fx * pb_z[j] + 45.0 * pa_y[j] * fl2_fx * fl1_fz * pb_zz[j] - 3.0 * pa_y[j] * fl1_fx * pb_zz[j] * fl1_fz * fl1_fgb - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_zz[j] * fl1_fx - 6.0 * pa_yzz[j] * pb_zz[j] * fl1_fz * fl1_fgb + 36.0 * pa_yzz[j] * fl1_fz * pb_zz[j] * fl1_fx + 48.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_zzz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_zzzz[j] + 6.0 * pa_y[j] * fl1_fz * fl1_fx * pb_zzzz[j] + 14.0 * pa_yzz[j] * fl1_fz * pb_zzzz[j]);

                t_zzz_xxxx[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pa_zzz[j] * fl2_fx + 4.5 * pa_z[j] * fl2_fx * pb_xx[j] + 3.0 * pa_zzz[j] * pb_xx[j] * fl1_fx + 1.5 * pa_z[j] * fl1_fx * pb_xxxx[j] + pa_zzz[j] * pb_xxxx[j]);

                t_zzz_xxxx[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl1_fz * fl3_fx + 7.5 * pa_zzz[j] * fl1_fz * fl2_fx - 9.0 * pa_z[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 9.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 6.0 * pa_zzz[j] * pb_xx[j] * fl1_fz * fl1_fgb + 45.0 * pa_z[j] * fl1_fz * fl2_fx * pb_xx[j] + 36.0 * pa_zzz[j] * fl1_fz * pb_xx[j] * fl1_fx - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxxx[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxxx[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xxxx[j]);

                t_zzz_xxxy[j] = fl_s_0_0 * (2.25 * pa_z[j] * fl2_fx * pb_xy[j] + 1.5 * pa_zzz[j] * pb_xy[j] * fl1_fx + 1.5 * pa_z[j] * fl1_fx * pb_xxxy[j] + pa_zzz[j] * pb_xxxy[j]);

                t_zzz_xxxy[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 3.0 * pa_zzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 22.5 * pa_z[j] * fl1_fz * fl2_fx * pb_xy[j] + 18.0 * pa_zzz[j] * fl1_fz * pb_xy[j] * fl1_fx - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxxy[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxxy[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xxxy[j]);

                t_zzz_xxxz[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_x[j] + 2.25 * pa_zz[j] * fl2_fx * pb_x[j] + 2.25 * pa_z[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pb_xxx[j] + 1.5 * pa_zzz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_zz[j] * fl1_fx * pb_xxx[j] + 1.5 * pa_z[j] * fl1_fx * pb_xxxz[j] + pa_zzz[j] * pb_xxxz[j]);

                t_zzz_xxxz[j] += fl_r_0_0 * (-2.25 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 4.5 * pa_zz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 9.0 * fl3_fx * fl1_fz * pb_x[j] + 22.5 * pa_zz[j] * fl1_fz * fl2_fx * pb_x[j] - 4.5 * pa_z[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxx[j] - 3.0 * pa_zzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 22.5 * pa_z[j] * fl1_fz * fl2_fx * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pb_xxx[j] + 18.0 * pa_zzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 18.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xxx[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxxz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxxz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xxxz[j]);

                t_zzz_xxyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.25 * pa_zzz[j] * fl2_fx + 0.75 * pa_z[j] * fl2_fx * pb_xx[j] + 0.75 * pa_z[j] * fl2_fx * pb_yy[j] + 0.5 * pa_zzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_zzz[j] * fl1_fx * pb_yy[j] + 1.5 * pa_z[j] * fl1_fx * pb_xxyy[j] + pa_zzz[j] * pb_xxyy[j]);

                t_zzz_xxyy[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl1_fz * fl3_fx + 2.5 * pa_zzz[j] * fl1_fz * fl2_fx - 1.5 * pa_z[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yy[j] - 1.5 * pa_z[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - pa_zzz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_zzz[j] * fl1_fz * fl1_fgb * pb_yy[j] + 7.5 * pa_z[j] * fl1_fz * fl2_fx * pb_xx[j] + 7.5 * pa_z[j] * fl1_fz * fl2_fx * pb_yy[j] + 6.0 * pa_zzz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_zzz[j] * fl1_fz * fl1_fx * pb_yy[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxyy[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxyy[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xxyy[j]);

                t_zzz_xxyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_y[j] + 0.75 * pa_zz[j] * fl2_fx * pb_y[j] + 0.75 * pa_z[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pb_xxy[j] + 0.5 * pa_zzz[j] * fl1_fx * pb_yz[j] + 1.5 * pa_zz[j] * fl1_fx * pb_xxy[j] + 1.5 * pa_z[j] * fl1_fx * pb_xxyz[j] + pa_zzz[j] * pb_xxyz[j]);

                t_zzz_xxyz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 1.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] + 3.0 * fl3_fx * fl1_fz * pb_y[j] + 7.5 * pa_zz[j] * fl1_fz * fl2_fx * pb_y[j] - 1.5 * pa_z[j] * fl1_fx * fl1_fz * fl1_fgb * pb_yz[j] - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xxy[j] - pa_zzz[j] * fl1_fz * fl1_fgb * pb_yz[j] + 7.5 * pa_z[j] * fl1_fz * fl2_fx * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pb_xxy[j] + 6.0 * pa_zzz[j] * fl1_fz * fl1_fx * pb_yz[j] + 18.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xxy[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxyz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxyz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_140_150(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (140,150)

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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_zzz = paDistances.data(19 * idx + 18);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzz_xxzz = primBuffer.data(150 * idx + 140);

            auto t_zzz_xyyy = primBuffer.data(150 * idx + 141);

            auto t_zzz_xyyz = primBuffer.data(150 * idx + 142);

            auto t_zzz_xyzz = primBuffer.data(150 * idx + 143);

            auto t_zzz_xzzz = primBuffer.data(150 * idx + 144);

            auto t_zzz_yyyy = primBuffer.data(150 * idx + 145);

            auto t_zzz_yyyz = primBuffer.data(150 * idx + 146);

            auto t_zzz_yyzz = primBuffer.data(150 * idx + 147);

            auto t_zzz_yzzz = primBuffer.data(150 * idx + 148);

            auto t_zzz_zzzz = primBuffer.data(150 * idx + 149);

            // Batch of Integrals (140,150)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, \
                                     s_0_0, t_zzz_xxzz, t_zzz_xyyy, t_zzz_xyyz, t_zzz_xyzz, t_zzz_xzzz, t_zzz_yyyy, \
                                     t_zzz_yyyz, t_zzz_yyzz, t_zzz_yzzz, t_zzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fga = fga[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_zzz_xxzz[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * fl3_fx * pb_z[j] + 0.25 * pa_zzz[j] * fl2_fx + 1.5 * pa_zz[j] * fl2_fx * pb_z[j] + 2.25 * pa_z[j] * fl2_fx * pb_xx[j] + 0.75 * pa_z[j] * fl2_fx * pb_zz[j] + 1.5 * fl2_fx * pb_xxz[j] + 0.5 * pa_zzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_zzz[j] * fl1_fx * pb_zz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xxz[j] + 1.5 * pa_z[j] * fl1_fx * pb_xxzz[j] + pa_zzz[j] * pb_xxzz[j]);

                t_zzz_xxzz[j] += fl_r_0_0 * (-3.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 1.5 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 6.0 * fl3_fx * fl1_fz * pb_z[j] + 2.5 * pa_zzz[j] * fl1_fz * fl2_fx + 15.0 * pa_zz[j] * fl1_fz * fl2_fx * pb_z[j] + 22.5 * pa_z[j] * fl2_fx * fl1_fz * pb_xx[j] - 1.5 * pa_z[j] * fl1_fx * pb_xx[j] * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 1.5 * pa_z[j] * fl1_fz * fl1_fga * pb_xx[j] * fl1_fx - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xxz[j] - pa_zzz[j] * pb_xx[j] * fl1_fz * fl1_fgb - pa_zzz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 7.5 * pa_z[j] * fl1_fz * fl2_fx * pb_zz[j] + 15.0 * fl2_fx * fl1_fz * pb_xxz[j] + 6.0 * pa_zzz[j] * fl1_fz * pb_xx[j] * fl1_fx + 6.0 * pa_zzz[j] * fl1_fz * fl1_fx * pb_zz[j] + 36.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xxz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxzz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxzz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xxzz[j]);

                t_zzz_xyyy[j] = fl_s_0_0 * (2.25 * pa_z[j] * fl2_fx * pb_xy[j] + 1.5 * pa_zzz[j] * pb_xy[j] * fl1_fx + 1.5 * pa_z[j] * fl1_fx * pb_xyyy[j] + pa_zzz[j] * pb_xyyy[j]);

                t_zzz_xyyy[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 3.0 * pa_zzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 22.5 * pa_z[j] * fl1_fz * fl2_fx * pb_xy[j] + 18.0 * pa_zzz[j] * fl1_fz * pb_xy[j] * fl1_fx - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xyyy[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xyyy[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xyyy[j]);

                t_zzz_xyyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pb_x[j] + 0.75 * pa_zz[j] * fl2_fx * pb_x[j] + 0.75 * pa_z[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pb_xyy[j] + 0.5 * pa_zzz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_zz[j] * fl1_fx * pb_xyy[j] + 1.5 * pa_z[j] * fl1_fx * pb_xyyz[j] + pa_zzz[j] * pb_xyyz[j]);

                t_zzz_xyyz[j] += fl_r_0_0 * (-0.75 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 1.5 * pa_zz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 3.0 * fl3_fx * fl1_fz * pb_x[j] + 7.5 * pa_zz[j] * fl1_fz * fl2_fx * pb_x[j] - 1.5 * pa_z[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xyy[j] - pa_zzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa_z[j] * fl1_fz * fl2_fx * pb_xz[j] + 7.5 * fl2_fx * fl1_fz * pb_xyy[j] + 6.0 * pa_zzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 18.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xyy[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xyyz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xyyz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xyyz[j]);

                t_zzz_xyzz[j] = fl_s_0_0 * (2.25 * pa_z[j] * fl2_fx * pb_xy[j] + 1.5 * fl2_fx * pb_xyz[j] + 0.5 * pa_zzz[j] * pb_xy[j] * fl1_fx + 3.0 * pa_zz[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_z[j] * fl1_fx * pb_xyzz[j] + pa_zzz[j] * pb_xyzz[j]);

                t_zzz_xyzz[j] += fl_r_0_0 * (22.5 * pa_z[j] * fl2_fx * fl1_fz * pb_xy[j] - 1.5 * pa_z[j] * fl1_fx * pb_xy[j] * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fz * fl1_fga * pb_xy[j] * fl1_fx - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xyz[j] - pa_zzz[j] * pb_xy[j] * fl1_fz * fl1_fgb + 15.0 * fl2_fx * fl1_fz * pb_xyz[j] + 6.0 * pa_zzz[j] * fl1_fz * pb_xy[j] * fl1_fx + 36.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xyz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xyzz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xyzz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xyzz[j]);

                t_zzz_xzzz[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_x[j] + 2.25 * pa_zz[j] * fl2_fx * pb_x[j] + 6.75 * pa_z[j] * fl2_fx * pb_xz[j] + 2.25 * fl2_fx * pb_xzz[j] + 1.5 * pa_zzz[j] * pb_xz[j] * fl1_fx + 4.5 * pa_zz[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_z[j] * fl1_fx * pb_xzzz[j] + pa_zzz[j] * pb_xzzz[j]);

                t_zzz_xzzz[j] += fl_r_0_0 * (15.0 * fl3_fx * fl1_fz * pb_x[j] - 2.25 * fl2_fx * pb_x[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_x[j] - 4.5 * pa_zz[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb + 22.5 * pa_zz[j] * fl1_fz * fl2_fx * pb_x[j] + 67.5 * pa_z[j] * fl2_fx * fl1_fz * pb_xz[j] - 4.5 * pa_z[j] * fl1_fx * pb_xz[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_xz[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_xzz[j] - 3.0 * pa_zzz[j] * pb_xz[j] * fl1_fz * fl1_fgb + 22.5 * fl2_fx * fl1_fz * pb_xzz[j] + 18.0 * pa_zzz[j] * fl1_fz * pb_xz[j] * fl1_fx + 54.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xzz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xzzz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xzzz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_xzzz[j]);

                t_zzz_yyyy[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pa_zzz[j] * fl2_fx + 4.5 * pa_z[j] * fl2_fx * pb_yy[j] + 3.0 * pa_zzz[j] * pb_yy[j] * fl1_fx + 1.5 * pa_z[j] * fl1_fx * pb_yyyy[j] + pa_zzz[j] * pb_yyyy[j]);

                t_zzz_yyyy[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl1_fz * fl3_fx + 7.5 * pa_zzz[j] * fl1_fz * fl2_fx - 9.0 * pa_z[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 9.0 * pa_z[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 6.0 * pa_zzz[j] * pb_yy[j] * fl1_fz * fl1_fgb + 45.0 * pa_z[j] * fl1_fz * fl2_fx * pb_yy[j] + 36.0 * pa_zzz[j] * fl1_fz * pb_yy[j] * fl1_fx - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_yyyy[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_yyyy[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_yyyy[j]);

                t_zzz_yyyz[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_y[j] + 2.25 * pa_zz[j] * fl2_fx * pb_y[j] + 2.25 * pa_z[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pb_yyy[j] + 1.5 * pa_zzz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_zz[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_z[j] * fl1_fx * pb_yyyz[j] + pa_zzz[j] * pb_yyyz[j]);

                t_zzz_yyyz[j] += fl_r_0_0 * (-2.25 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 4.5 * pa_zz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 9.0 * fl3_fx * fl1_fz * pb_y[j] + 22.5 * pa_zz[j] * fl1_fz * fl2_fx * pb_y[j] - 4.5 * pa_z[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yyy[j] - 3.0 * pa_zzz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 22.5 * pa_z[j] * fl1_fz * fl2_fx * pb_yz[j] + 7.5 * fl2_fx * fl1_fz * pb_yyy[j] + 18.0 * pa_zzz[j] * fl1_fz * pb_yz[j] * fl1_fx + 18.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_yyy[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_yyyz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_yyyz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_yyyz[j]);

                t_zzz_yyzz[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * fl3_fx * pb_z[j] + 0.25 * pa_zzz[j] * fl2_fx + 1.5 * pa_zz[j] * fl2_fx * pb_z[j] + 2.25 * pa_z[j] * fl2_fx * pb_yy[j] + 0.75 * pa_z[j] * fl2_fx * pb_zz[j] + 1.5 * fl2_fx * pb_yyz[j] + 0.5 * pa_zzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_zzz[j] * fl1_fx * pb_zz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_z[j] * fl1_fx * pb_yyzz[j] + pa_zzz[j] * pb_yyzz[j]);

                t_zzz_yyzz[j] += fl_r_0_0 * (-3.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * fl2_fx * fl1_fz * fl1_fgb * pb_z[j] - 1.5 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] + 6.0 * fl3_fx * fl1_fz * pb_z[j] + 2.5 * pa_zzz[j] * fl1_fz * fl2_fx + 15.0 * pa_zz[j] * fl1_fz * fl2_fx * pb_z[j] + 22.5 * pa_z[j] * fl2_fx * fl1_fz * pb_yy[j] - 1.5 * pa_z[j] * fl1_fx * pb_yy[j] * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fx * fl1_fz * fl1_fgb * pb_zz[j] - 1.5 * pa_z[j] * fl1_fz * fl1_fga * pb_yy[j] * fl1_fx - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_yyz[j] - pa_zzz[j] * pb_yy[j] * fl1_fz * fl1_fgb - pa_zzz[j] * fl1_fz * fl1_fgb * pb_zz[j] + 7.5 * pa_z[j] * fl1_fz * fl2_fx * pb_zz[j] + 15.0 * fl2_fx * fl1_fz * pb_yyz[j] + 6.0 * pa_zzz[j] * fl1_fz * pb_yy[j] * fl1_fx + 6.0 * pa_zzz[j] * fl1_fz * fl1_fx * pb_zz[j] + 36.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_yyz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_yyzz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_yyzz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_yyzz[j]);

                t_zzz_yzzz[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_y[j] + 2.25 * pa_zz[j] * fl2_fx * pb_y[j] + 6.75 * pa_z[j] * fl2_fx * pb_yz[j] + 2.25 * fl2_fx * pb_yzz[j] + 1.5 * pa_zzz[j] * pb_yz[j] * fl1_fx + 4.5 * pa_zz[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_z[j] * fl1_fx * pb_yzzz[j] + pa_zzz[j] * pb_yzzz[j]);

                t_zzz_yzzz[j] += fl_r_0_0 * (15.0 * fl3_fx * fl1_fz * pb_y[j] - 2.25 * fl2_fx * pb_y[j] * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga * pb_y[j] - 4.5 * pa_zz[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb + 22.5 * pa_zz[j] * fl1_fz * fl2_fx * pb_y[j] + 67.5 * pa_z[j] * fl2_fx * fl1_fz * pb_yz[j] - 4.5 * pa_z[j] * fl1_fx * pb_yz[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_yz[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_yzz[j] - 3.0 * pa_zzz[j] * pb_yz[j] * fl1_fz * fl1_fgb + 22.5 * fl2_fx * fl1_fz * pb_yzz[j] + 18.0 * pa_zzz[j] * fl1_fz * pb_yz[j] * fl1_fx + 54.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_yzz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_yzzz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_yzzz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_yzzz[j]);

                t_zzz_zzzz[j] = fl_s_0_0 * (5.625 * pa_z[j] * fl3_fx + 7.5 * fl3_fx * pb_z[j] + 0.75 * pa_zzz[j] * fl2_fx + 9.0 * pa_zz[j] * fl2_fx * pb_z[j] + 13.5 * pa_z[j] * fl2_fx * pb_zz[j] + 3.0 * fl2_fx * pb_zzz[j] + 3.0 * pa_zzz[j] * pb_zz[j] * fl1_fx + 6.0 * pa_zz[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_z[j] * fl1_fx * pb_zzzz[j] + pa_zzz[j] * pb_zzzz[j]);

                t_zzz_zzzz[j] += fl_r_0_0 * (-13.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 45.0 * pa_z[j] * fl3_fx * fl1_fz + 60.0 * fl3_fx * fl1_fz * pb_z[j] - 2.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * fl2_fx * pb_z[j] * fl1_fz * fl1_fgb - 9.0 * fl2_fx * fl1_fz * fl1_fga * pb_z[j] - 3.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa_zz[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb + 7.5 * pa_zzz[j] * fl1_fz * fl2_fx + 90.0 * pa_zz[j] * fl1_fz * fl2_fx * pb_z[j] + 135.0 * pa_z[j] * fl2_fx * fl1_fz * pb_zz[j] - 9.0 * pa_z[j] * fl1_fx * pb_zz[j] * fl1_fz * fl1_fgb - 9.0 * pa_z[j] * fl1_fz * fl1_fga * pb_zz[j] * fl1_fx - 6.0 * fl1_fx * fl1_fz * fl1_fga * pb_zzz[j] - 6.0 * pa_zzz[j] * pb_zz[j] * fl1_fz * fl1_fgb + 30.0 * fl2_fx * fl1_fz * pb_zzz[j] + 36.0 * pa_zzz[j] * fl1_fz * pb_zz[j] * fl1_fx + 72.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_zzz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_zzzz[j] + 18.0 * pa_z[j] * fl1_fz * fl1_fx * pb_zzzz[j] + 14.0 * pa_zzz[j] * fl1_fz * pb_zzzz[j]);
            }

            idx++;
        }
    }


} // kinrecfunc namespace

