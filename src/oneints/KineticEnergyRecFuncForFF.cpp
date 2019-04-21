//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForFF.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForFF(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForFF_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_30_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_40_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_50_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_60_70(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_70_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_80_90(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_90_100(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                                  braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compKineticEnergyForFF_0_10(      CMemBlock2D<double>& primBuffer,
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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxx_xxx = primBuffer.data(100 * idx);

            auto t_xxx_xxy = primBuffer.data(100 * idx + 1);

            auto t_xxx_xxz = primBuffer.data(100 * idx + 2);

            auto t_xxx_xyy = primBuffer.data(100 * idx + 3);

            auto t_xxx_xyz = primBuffer.data(100 * idx + 4);

            auto t_xxx_xzz = primBuffer.data(100 * idx + 5);

            auto t_xxx_yyy = primBuffer.data(100 * idx + 6);

            auto t_xxx_yyz = primBuffer.data(100 * idx + 7);

            auto t_xxx_yzz = primBuffer.data(100 * idx + 8);

            auto t_xxx_zzz = primBuffer.data(100 * idx + 9);

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, r_0_0, s_0_0, t_xxx_xxx, t_xxx_xxy, t_xxx_xxz, t_xxx_xyy, t_xxx_xyz, \
                                     t_xxx_xzz, t_xxx_yyy, t_xxx_yyz, t_xxx_yzz, t_xxx_zzz: VLX_ALIGN)
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

                t_xxx_xxx[j] = fl_s_0_0 * (1.875 * fl3_fx + 2.25 * pa_xx[j] * fl2_fx + 6.75 * pa_x[j] * fl2_fx * pb_x[j] + 2.25 * fl2_fx * pb_xx[j] + 1.5 * pa_xxx[j] * pb_x[j] * fl1_fx + 4.5 * pa_xx[j] * fl1_fx * pb_xx[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxx[j] + pa_xxx[j] * pb_xxx[j]);

                t_xxx_xxx[j] += fl_r_0_0 * (11.25 * fl3_fx * fl1_fz - 2.25 * fl2_fx * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa_xx[j] * fl1_fz * fl2_fx + 54.0 * pa_x[j] * fl2_fx * fl1_fz * pb_x[j] - 4.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_xx[j] - 3.0 * pa_xxx[j] * pb_x[j] * fl1_fz * fl1_fgb + 18.0 * fl2_fx * fl1_fz * pb_xx[j] + 15.0 * pa_xxx[j] * fl1_fz * pb_x[j] * fl1_fx + 45.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xx[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxx[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxx[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_xxx[j]);

                t_xxx_xxy[j] = fl_s_0_0 * (2.25 * pa_x[j] * fl2_fx * pb_y[j] + 1.5 * fl2_fx * pb_xy[j] + 0.5 * pa_xxx[j] * fl1_fx * pb_y[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xy[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxy[j] + pa_xxx[j] * pb_xxy[j]);

                t_xxx_xxy[j] += fl_r_0_0 * (18.0 * pa_x[j] * fl2_fx * fl1_fz * pb_y[j] - 1.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 1.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_y[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xy[j] - pa_xxx[j] * fl1_fz * fl1_fgb * pb_y[j] + 12.0 * fl2_fx * fl1_fz * pb_xy[j] + 5.0 * pa_xxx[j] * fl1_fz * fl1_fx * pb_y[j] + 30.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xy[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxy[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_xxy[j]);

                t_xxx_xxz[j] = fl_s_0_0 * (2.25 * pa_x[j] * fl2_fx * pb_z[j] + 1.5 * fl2_fx * pb_xz[j] + 0.5 * pa_xxx[j] * fl1_fx * pb_z[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xxz[j] + pa_xxx[j] * pb_xxz[j]);

                t_xxx_xxz[j] += fl_r_0_0 * (18.0 * pa_x[j] * fl2_fx * fl1_fz * pb_z[j] - 1.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 1.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xz[j] - pa_xxx[j] * fl1_fz * fl1_fgb * pb_z[j] + 12.0 * fl2_fx * fl1_fz * pb_xz[j] + 5.0 * pa_xxx[j] * fl1_fz * fl1_fx * pb_z[j] + 30.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xxz[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_xxz[j]);

                t_xxx_xyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * pa_x[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pb_yy[j] + 0.5 * pa_xxx[j] * pb_x[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_yy[j] + 1.5 * pa_x[j] * fl1_fx * pb_xyy[j] + pa_xxx[j] * pb_xyy[j]);

                t_xxx_xyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_xx[j] * fl1_fz * fl2_fx - 1.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yy[j] - pa_xxx[j] * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl1_fz * fl2_fx * pb_x[j] + 6.0 * fl2_fx * fl1_fz * pb_yy[j] + 5.0 * pa_xxx[j] * fl1_fz * pb_x[j] * fl1_fx + 15.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yy[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xyy[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_xyy[j]);

                t_xxx_xyz[j] = fl_s_0_0 * (0.75 * fl2_fx * pb_yz[j] + 1.5 * pa_xx[j] * fl1_fx * pb_yz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xyz[j] + pa_xxx[j] * pb_xyz[j]);

                t_xxx_xyz[j] += fl_r_0_0 * (-1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yz[j] + 6.0 * fl2_fx * fl1_fz * pb_yz[j] + 15.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xyz[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_xyz[j]);

                t_xxx_xzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * pa_x[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pb_zz[j] + 0.5 * pa_xxx[j] * pb_x[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_zz[j] + 1.5 * pa_x[j] * fl1_fx * pb_xzz[j] + pa_xxx[j] * pb_xzz[j]);

                t_xxx_xzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_xx[j] * fl1_fz * fl2_fx - 1.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_zz[j] - pa_xxx[j] * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl1_fz * fl2_fx * pb_x[j] + 6.0 * fl2_fx * fl1_fz * pb_zz[j] + 5.0 * pa_xxx[j] * fl1_fz * pb_x[j] * fl1_fx + 15.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_zz[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_xzz[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xzz[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_xzz[j]);

                t_xxx_yyy[j] = fl_s_0_0 * (2.25 * pa_x[j] * fl2_fx * pb_y[j] + 1.5 * pa_xxx[j] * pb_y[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pb_yyy[j] + pa_xxx[j] * pb_yyy[j]);

                t_xxx_yyy[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 3.0 * pa_xxx[j] * pb_y[j] * fl1_fz * fl1_fgb + 18.0 * pa_x[j] * fl1_fz * fl2_fx * pb_y[j] + 15.0 * pa_xxx[j] * fl1_fz * pb_y[j] * fl1_fx - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yyy[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyy[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_yyy[j]);

                t_xxx_yyz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx * pb_z[j] + 0.5 * pa_xxx[j] * fl1_fx * pb_z[j] + 1.5 * pa_x[j] * fl1_fx * pb_yyz[j] + pa_xxx[j] * pb_yyz[j]);

                t_xxx_yyz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 1.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - pa_xxx[j] * fl1_fz * fl1_fgb * pb_z[j] + 6.0 * pa_x[j] * fl1_fz * fl2_fx * pb_z[j] + 5.0 * pa_xxx[j] * fl1_fz * fl1_fx * pb_z[j] - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yyz[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyz[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_yyz[j]);

                t_xxx_yzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx * pb_y[j] + 0.5 * pa_xxx[j] * pb_y[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pb_yzz[j] + pa_xxx[j] * pb_yzz[j]);

                t_xxx_yzz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - pa_xxx[j] * pb_y[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl1_fz * fl2_fx * pb_y[j] + 5.0 * pa_xxx[j] * fl1_fz * pb_y[j] * fl1_fx - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_yzz[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yzz[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_yzz[j]);

                t_xxx_zzz[j] = fl_s_0_0 * (2.25 * pa_x[j] * fl2_fx * pb_z[j] + 1.5 * pa_xxx[j] * pb_z[j] * fl1_fx + 1.5 * pa_x[j] * fl1_fx * pb_zzz[j] + pa_xxx[j] * pb_zzz[j]);

                t_xxx_zzz[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb - 4.5 * pa_x[j] * fl1_fz * fl1_fga * pb_z[j] * fl1_fx - 3.0 * pa_xxx[j] * pb_z[j] * fl1_fz * fl1_fgb + 18.0 * pa_x[j] * fl1_fz * fl2_fx * pb_z[j] + 15.0 * pa_xxx[j] * fl1_fz * pb_z[j] * fl1_fx - 3.0 * pa_x[j] * fl1_fz * fl1_fga * pb_zzz[j] + 15.0 * pa_x[j] * fl1_fz * fl1_fx * pb_zzz[j] + 12.0 * pa_xxx[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_10_20(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xxy = paDistances.data(19 * idx + 10);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxy_xxx = primBuffer.data(100 * idx + 10);

            auto t_xxy_xxy = primBuffer.data(100 * idx + 11);

            auto t_xxy_xxz = primBuffer.data(100 * idx + 12);

            auto t_xxy_xyy = primBuffer.data(100 * idx + 13);

            auto t_xxy_xyz = primBuffer.data(100 * idx + 14);

            auto t_xxy_xzz = primBuffer.data(100 * idx + 15);

            auto t_xxy_yyy = primBuffer.data(100 * idx + 16);

            auto t_xxy_yyz = primBuffer.data(100 * idx + 17);

            auto t_xxy_yzz = primBuffer.data(100 * idx + 18);

            auto t_xxy_zzz = primBuffer.data(100 * idx + 19);

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxy_xxx, t_xxy_xxy, t_xxy_xxz, t_xxy_xyy, \
                                     t_xxy_xyz, t_xxy_xzz, t_xxy_yyy, t_xxy_yyz, t_xxy_yzz, t_xxy_zzz: VLX_ALIGN)
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

                t_xxy_xxx[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl2_fx + 2.25 * fl2_fx * pa_y[j] * pb_x[j] + 1.5 * pa_xxy[j] * pb_x[j] * fl1_fx + 3.0 * pa_xy[j] * fl1_fx * pb_xx[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxx[j] + pa_xxy[j] * pb_xxx[j]);

                t_xxy_xxx[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_xy[j] * fl2_fx * fl1_fz + 18.0 * fl2_fx * fl1_fz * pa_y[j] * pb_x[j] - 1.5 * fl1_fx * pa_y[j] * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_x[j] * fl1_fx - 3.0 * pa_xxy[j] * pb_x[j] * fl1_fz * fl1_fgb + 15.0 * pa_xxy[j] * fl1_fz * pb_x[j] * fl1_fx + 30.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xx[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxx[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxx[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_xxx[j]);

                t_xxy_xxy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + pa_x[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_y[j] * pb_y[j] + 0.25 * fl2_fx * pb_xx[j] + 0.5 * pa_xxy[j] * fl1_fx * pb_y[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xx[j] + 2.0 * pa_xy[j] * fl1_fx * pb_xy[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxy[j] + pa_xxy[j] * pb_xxy[j]);

                t_xxy_xxy[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xx[j] * fl1_fz * fl2_fx + 8.0 * pa_x[j] * fl2_fx * fl1_fz * pb_x[j] + 6.0 * fl2_fx * fl1_fz * pa_y[j] * pb_y[j] - 0.5 * fl1_fx * pa_y[j] * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * fl1_fz * fl1_fga * pa_y[j] * fl1_fx * pb_y[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xx[j] - pa_xxy[j] * fl1_fz * fl1_fgb * pb_y[j] + 2.0 * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * pa_xxy[j] * fl1_fz * fl1_fx * pb_y[j] + 5.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xx[j] + 20.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xy[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxy[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxy[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_xxy[j]);

                t_xxy_xxz[j] = fl_s_0_0 * (0.75 * fl2_fx * pa_y[j] * pb_z[j] + 0.5 * pa_xxy[j] * fl1_fx * pb_z[j] + 2.0 * pa_xy[j] * fl1_fx * pb_xz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xxz[j] + pa_xxy[j] * pb_xxz[j]);

                t_xxy_xxz[j] += fl_r_0_0 * (6.0 * fl2_fx * fl1_fz * pa_y[j] * pb_z[j] - 0.5 * fl1_fx * pa_y[j] * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * fl1_fz * fl1_fga * pa_y[j] * fl1_fx * pb_z[j] - pa_xxy[j] * fl1_fz * fl1_fgb * pb_z[j] + 5.0 * pa_xxy[j] * fl1_fz * fl1_fx * pb_z[j] + 20.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_xz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xxz[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xxz[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_xxz[j]);

                t_xxy_xyy[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + pa_x[j] * fl2_fx * pb_y[j] + 0.25 * fl2_fx * pa_y[j] * pb_x[j] + 0.5 * fl2_fx * pb_xy[j] + 0.5 * pa_xxy[j] * pb_x[j] * fl1_fx + pa_xx[j] * fl1_fx * pb_xy[j] + pa_xy[j] * fl1_fx * pb_yy[j] + 0.5 * fl1_fx * pa_y[j] * pb_xyy[j] + pa_xxy[j] * pb_xyy[j]);

                t_xxy_xyy[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xy[j] * fl2_fx * fl1_fz + 8.0 * pa_x[j] * fl2_fx * fl1_fz * pb_y[j] - 0.5 * fl1_fx * pa_y[j] * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_y[j] * pb_x[j] * fl1_fx - fl1_fz * fl1_fga * fl1_fx * pb_xy[j] - pa_xxy[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pa_y[j] * pb_x[j] + 4.0 * fl2_fx * fl1_fz * pb_xy[j] + 5.0 * pa_xxy[j] * fl1_fz * pb_x[j] * fl1_fx + 10.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xy[j] + 10.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_yy[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xyy[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xyy[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_xyy[j]);

                t_xxy_xyz[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl2_fx * pb_z[j] + 0.25 * fl2_fx * pb_xz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xz[j] + pa_xy[j] * fl1_fx * pb_yz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xyz[j] + pa_xxy[j] * pb_xyz[j]);

                t_xxy_xyz[j] += fl_r_0_0 * (4.0 * pa_x[j] * fl2_fx * fl1_fz * pb_z[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xz[j] + 2.0 * fl2_fx * fl1_fz * pb_xz[j] + 5.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xz[j] + 10.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_yz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xyz[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xyz[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_xyz[j]);

                t_xxy_xzz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + 0.25 * fl2_fx * pa_y[j] * pb_x[j] + 0.5 * pa_xxy[j] * pb_x[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_zz[j] + 0.5 * fl1_fx * pa_y[j] * pb_xzz[j] + pa_xxy[j] * pb_xzz[j]);

                t_xxy_xzz[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xy[j] * fl2_fx * fl1_fz - 0.5 * fl1_fx * pa_y[j] * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_y[j] * pb_x[j] * fl1_fx - pa_xxy[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pa_y[j] * pb_x[j] + 5.0 * pa_xxy[j] * fl1_fz * pb_x[j] * fl1_fx + 10.0 * pa_xy[j] * fl1_fx * fl1_fz * pb_zz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_xzz[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_xzz[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_xzz[j]);

                t_xxy_yyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * fl2_fx * pa_y[j] * pb_y[j] + 0.75 * fl2_fx * pb_yy[j] + 1.5 * pa_xxy[j] * pb_y[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_yy[j] + 0.5 * fl1_fx * pa_y[j] * pb_yyy[j] + pa_xxy[j] * pb_yyy[j]);

                t_xxy_yyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_xx[j] * fl1_fz * fl2_fx - 1.5 * fl1_fx * pa_y[j] * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_y[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - 3.0 * pa_xxy[j] * pb_y[j] * fl1_fz * fl1_fgb + 6.0 * fl2_fx * fl1_fz * pa_y[j] * pb_y[j] + 6.0 * fl2_fx * fl1_fz * pb_yy[j] + 15.0 * pa_xxy[j] * fl1_fz * pb_y[j] * fl1_fx + 15.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yy[j] - fl1_fz * fl1_fga * pa_y[j] * pb_yyy[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_yyy[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_yyy[j]);

                t_xxy_yyz[j] = fl_s_0_0 * (0.25 * fl2_fx * pa_y[j] * pb_z[j] + 0.5 * fl2_fx * pb_yz[j] + 0.5 * pa_xxy[j] * fl1_fx * pb_z[j] + pa_xx[j] * fl1_fx * pb_yz[j] + 0.5 * fl1_fx * pa_y[j] * pb_yyz[j] + pa_xxy[j] * pb_yyz[j]);

                t_xxy_yyz[j] += fl_r_0_0 * (-0.5 * fl1_fx * pa_y[j] * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * fl1_fz * fl1_fga * pa_y[j] * fl1_fx * pb_z[j] - fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - pa_xxy[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * fl2_fx * fl1_fz * pa_y[j] * pb_z[j] + 4.0 * fl2_fx * fl1_fz * pb_yz[j] + 5.0 * pa_xxy[j] * fl1_fz * fl1_fx * pb_z[j] + 10.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_yyz[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_yyz[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_yyz[j]);

                t_xxy_yzz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + 0.25 * fl2_fx * pa_y[j] * pb_y[j] + 0.25 * fl2_fx * pb_zz[j] + 0.5 * pa_xxy[j] * pb_y[j] * fl1_fx + 0.5 * pa_xx[j] * fl1_fx * pb_zz[j] + 0.5 * fl1_fx * pa_y[j] * pb_yzz[j] + pa_xxy[j] * pb_yzz[j]);

                t_xxy_yzz[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_xx[j] * fl1_fz * fl2_fx - 0.5 * fl1_fx * pa_y[j] * pb_y[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_y[j] * pb_y[j] * fl1_fx - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - pa_xxy[j] * pb_y[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pa_y[j] * pb_y[j] + 2.0 * fl2_fx * fl1_fz * pb_zz[j] + 5.0 * pa_xxy[j] * fl1_fz * pb_y[j] * fl1_fx + 5.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_zz[j] - fl1_fz * fl1_fga * pa_y[j] * pb_yzz[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_yzz[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_yzz[j]);

                t_xxy_zzz[j] = fl_s_0_0 * (0.75 * fl2_fx * pa_y[j] * pb_z[j] + 1.5 * pa_xxy[j] * pb_z[j] * fl1_fx + 0.5 * fl1_fx * pa_y[j] * pb_zzz[j] + pa_xxy[j] * pb_zzz[j]);

                t_xxy_zzz[j] += fl_r_0_0 * (-1.5 * fl1_fx * pa_y[j] * pb_z[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_y[j] * pb_z[j] * fl1_fx - 3.0 * pa_xxy[j] * pb_z[j] * fl1_fz * fl1_fgb + 6.0 * fl2_fx * fl1_fz * pa_y[j] * pb_z[j] + 15.0 * pa_xxy[j] * fl1_fz * pb_z[j] * fl1_fx - fl1_fz * fl1_fga * pa_y[j] * pb_zzz[j] + 5.0 * fl1_fx * fl1_fz * pa_y[j] * pb_zzz[j] + 12.0 * pa_xxy[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_20_30(      CMemBlock2D<double>& primBuffer,
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

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xz = paDistances.data(19 * idx + 5);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxz = paDistances.data(19 * idx + 11);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxz_xxx = primBuffer.data(100 * idx + 20);

            auto t_xxz_xxy = primBuffer.data(100 * idx + 21);

            auto t_xxz_xxz = primBuffer.data(100 * idx + 22);

            auto t_xxz_xyy = primBuffer.data(100 * idx + 23);

            auto t_xxz_xyz = primBuffer.data(100 * idx + 24);

            auto t_xxz_xzz = primBuffer.data(100 * idx + 25);

            auto t_xxz_yyy = primBuffer.data(100 * idx + 26);

            auto t_xxz_yyz = primBuffer.data(100 * idx + 27);

            auto t_xxz_yzz = primBuffer.data(100 * idx + 28);

            auto t_xxz_zzz = primBuffer.data(100 * idx + 29);

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxz_xxx, t_xxz_xxy, t_xxz_xxz, t_xxz_xyy, \
                                     t_xxz_xyz, t_xxz_xzz, t_xxz_yyy, t_xxz_yyz, t_xxz_yzz, t_xxz_zzz: VLX_ALIGN)
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

                t_xxz_xxx[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl2_fx + 2.25 * fl2_fx * pa_z[j] * pb_x[j] + 1.5 * pa_xxz[j] * pb_x[j] * fl1_fx + 3.0 * pa_xz[j] * fl1_fx * pb_xx[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxx[j] + pa_xxz[j] * pb_xxx[j]);

                t_xxz_xxx[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_xz[j] * fl2_fx * fl1_fz + 18.0 * fl2_fx * fl1_fz * pa_z[j] * pb_x[j] - 1.5 * fl1_fx * pa_z[j] * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_x[j] * fl1_fx - 3.0 * pa_xxz[j] * pb_x[j] * fl1_fz * fl1_fgb + 15.0 * pa_xxz[j] * fl1_fz * pb_x[j] * fl1_fx + 30.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xx[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxx[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxx[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_xxx[j]);

                t_xxz_xxy[j] = fl_s_0_0 * (0.75 * fl2_fx * pa_z[j] * pb_y[j] + 0.5 * pa_xxz[j] * fl1_fx * pb_y[j] + 2.0 * pa_xz[j] * fl1_fx * pb_xy[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxy[j] + pa_xxz[j] * pb_xxy[j]);

                t_xxz_xxy[j] += fl_r_0_0 * (6.0 * fl2_fx * fl1_fz * pa_z[j] * pb_y[j] - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_y[j] - pa_xxz[j] * fl1_fz * fl1_fgb * pb_y[j] + 5.0 * pa_xxz[j] * fl1_fz * fl1_fx * pb_y[j] + 20.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxy[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxy[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_xxy[j]);

                t_xxz_xxz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + pa_x[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_z[j] * pb_z[j] + 0.25 * fl2_fx * pb_xx[j] + 0.5 * pa_xxz[j] * fl1_fx * pb_z[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xx[j] + 2.0 * pa_xz[j] * fl1_fx * pb_xz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxz[j] + pa_xxz[j] * pb_xxz[j]);

                t_xxz_xxz[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xx[j] * fl1_fz * fl2_fx + 8.0 * pa_x[j] * fl2_fx * fl1_fz * pb_x[j] + 6.0 * fl2_fx * fl1_fz * pa_z[j] * pb_z[j] - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xx[j] - pa_xxz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * pa_xxz[j] * fl1_fz * fl1_fx * pb_z[j] + 5.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xx[j] + 20.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxz[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_xxz[j]);

                t_xxz_xyy[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + 0.25 * fl2_fx * pa_z[j] * pb_x[j] + 0.5 * pa_xxz[j] * pb_x[j] * fl1_fx + pa_xz[j] * fl1_fx * pb_yy[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyy[j] + pa_xxz[j] * pb_xyy[j]);

                t_xxz_xyy[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xz[j] * fl2_fx * fl1_fz - 0.5 * fl1_fx * pa_z[j] * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_x[j] * fl1_fx - pa_xxz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_x[j] + 5.0 * pa_xxz[j] * fl1_fz * pb_x[j] * fl1_fx + 10.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyy[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyy[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_xyy[j]);

                t_xxz_xyz[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl2_fx * pb_y[j] + 0.25 * fl2_fx * pb_xy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xy[j] + pa_xz[j] * fl1_fx * pb_yz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyz[j] + pa_xxz[j] * pb_xyz[j]);

                t_xxz_xyz[j] += fl_r_0_0 * (4.0 * pa_x[j] * fl2_fx * fl1_fz * pb_y[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xy[j] + 2.0 * fl2_fx * fl1_fz * pb_xy[j] + 5.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xy[j] + 10.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyz[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_xyz[j]);

                t_xxz_xzz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + pa_x[j] * fl2_fx * pb_z[j] + 0.25 * fl2_fx * pa_z[j] * pb_x[j] + 0.5 * fl2_fx * pb_xz[j] + 0.5 * pa_xxz[j] * pb_x[j] * fl1_fx + pa_xx[j] * fl1_fx * pb_xz[j] + pa_xz[j] * fl1_fx * pb_zz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xzz[j] + pa_xxz[j] * pb_xzz[j]);

                t_xxz_xzz[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xz[j] * fl2_fx * fl1_fz + 8.0 * pa_x[j] * fl2_fx * fl1_fz * pb_z[j] - 0.5 * fl1_fx * pa_z[j] * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_x[j] * fl1_fx - fl1_fz * fl1_fga * fl1_fx * pb_xz[j] - pa_xxz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_x[j] + 4.0 * fl2_fx * fl1_fz * pb_xz[j] + 5.0 * pa_xxz[j] * fl1_fz * pb_x[j] * fl1_fx + 10.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_xz[j] + 10.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_zz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xzz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xzz[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_xzz[j]);

                t_xxz_yyy[j] = fl_s_0_0 * (0.75 * fl2_fx * pa_z[j] * pb_y[j] + 1.5 * pa_xxz[j] * pb_y[j] * fl1_fx + 0.5 * fl1_fx * pa_z[j] * pb_yyy[j] + pa_xxz[j] * pb_yyy[j]);

                t_xxz_yyy[j] += fl_r_0_0 * (-1.5 * fl1_fx * pa_z[j] * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_y[j] * fl1_fx - 3.0 * pa_xxz[j] * pb_y[j] * fl1_fz * fl1_fgb + 6.0 * fl2_fx * fl1_fz * pa_z[j] * pb_y[j] + 15.0 * pa_xxz[j] * fl1_fz * pb_y[j] * fl1_fx - fl1_fz * fl1_fga * pa_z[j] * pb_yyy[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyy[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_yyy[j]);

                t_xxz_yyz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + 0.25 * fl2_fx * pa_z[j] * pb_z[j] + 0.25 * fl2_fx * pb_yy[j] + 0.5 * pa_xxz[j] * fl1_fx * pb_z[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yy[j] + 0.5 * fl1_fx * pa_z[j] * pb_yyz[j] + pa_xxz[j] * pb_yyz[j]);

                t_xxz_yyz[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_xx[j] * fl1_fz * fl2_fx - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - pa_xxz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_z[j] + 2.0 * fl2_fx * fl1_fz * pb_yy[j] + 5.0 * pa_xxz[j] * fl1_fz * fl1_fx * pb_z[j] + 5.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yyz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyz[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_yyz[j]);

                t_xxz_yzz[j] = fl_s_0_0 * (0.25 * fl2_fx * pa_z[j] * pb_y[j] + 0.5 * fl2_fx * pb_yz[j] + 0.5 * pa_xxz[j] * pb_y[j] * fl1_fx + pa_xx[j] * fl1_fx * pb_yz[j] + 0.5 * fl1_fx * pa_z[j] * pb_yzz[j] + pa_xxz[j] * pb_yzz[j]);

                t_xxz_yzz[j] += fl_r_0_0 * (-0.5 * fl1_fx * pa_z[j] * pb_y[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_y[j] * fl1_fx - fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - pa_xxz[j] * pb_y[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_y[j] + 4.0 * fl2_fx * fl1_fz * pb_yz[j] + 5.0 * pa_xxz[j] * fl1_fz * pb_y[j] * fl1_fx + 10.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_yz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yzz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yzz[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_yzz[j]);

                t_xxz_zzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * fl2_fx * pa_z[j] * pb_z[j] + 0.75 * fl2_fx * pb_zz[j] + 1.5 * pa_xxz[j] * pb_z[j] * fl1_fx + 1.5 * pa_xx[j] * fl1_fx * pb_zz[j] + 0.5 * fl1_fx * pa_z[j] * pb_zzz[j] + pa_xxz[j] * pb_zzz[j]);

                t_xxz_zzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_xx[j] * fl1_fz * fl2_fx - 1.5 * fl1_fx * pa_z[j] * pb_z[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_z[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - 3.0 * pa_xxz[j] * pb_z[j] * fl1_fz * fl1_fgb + 6.0 * fl2_fx * fl1_fz * pa_z[j] * pb_z[j] + 6.0 * fl2_fx * fl1_fz * pb_zz[j] + 15.0 * pa_xxz[j] * fl1_fz * pb_z[j] * fl1_fx + 15.0 * pa_xx[j] * fl1_fz * fl1_fx * pb_zz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_zzz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_zzz[j] + 12.0 * pa_xxz[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_30_40(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(19 * idx + 1);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_yy = paDistances.data(19 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyy = paDistances.data(19 * idx + 12);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyy_xxx = primBuffer.data(100 * idx + 30);

            auto t_xyy_xxy = primBuffer.data(100 * idx + 31);

            auto t_xyy_xxz = primBuffer.data(100 * idx + 32);

            auto t_xyy_xyy = primBuffer.data(100 * idx + 33);

            auto t_xyy_xyz = primBuffer.data(100 * idx + 34);

            auto t_xyy_xzz = primBuffer.data(100 * idx + 35);

            auto t_xyy_yyy = primBuffer.data(100 * idx + 36);

            auto t_xyy_yyz = primBuffer.data(100 * idx + 37);

            auto t_xyy_yzz = primBuffer.data(100 * idx + 38);

            auto t_xyy_zzz = primBuffer.data(100 * idx + 39);

            // Batch of Integrals (30,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyy_xxx, t_xyy_xxy, t_xyy_xxz, t_xyy_xyy, \
                                     t_xyy_xyz, t_xyy_xzz, t_xyy_yyy, t_xyy_yyz, t_xyy_yzz, t_xyy_zzz: VLX_ALIGN)
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

                t_xyy_xxx[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * fl2_fx * pa_yy[j] + 0.75 * pa_x[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pb_xx[j] + 1.5 * pa_xyy[j] * pb_x[j] * fl1_fx + 1.5 * fl1_fx * pa_yy[j] * pb_xx[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxx[j] + pa_xyy[j] * pb_xxx[j]);

                t_xyy_xxx[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * fl1_fx * pa_yy[j] * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * fl2_fx * pa_yy[j] * fl1_fz - 1.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xx[j] - 3.0 * pa_xyy[j] * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl1_fz * fl2_fx * pb_x[j] + 6.0 * fl2_fx * fl1_fz * pb_xx[j] + 15.0 * pa_xyy[j] * fl1_fz * pb_x[j] * fl1_fx + 15.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xx[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxx[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxx[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_xxx[j]);

                t_xyy_xxy[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + fl2_fx * pa_y[j] * pb_x[j] + 0.25 * pa_x[j] * fl2_fx * pb_y[j] + 0.5 * fl2_fx * pb_xy[j] + 0.5 * pa_xyy[j] * fl1_fx * pb_y[j] + pa_xy[j] * fl1_fx * pb_xx[j] + fl1_fx * pa_yy[j] * pb_xy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxy[j] + pa_xyy[j] * pb_xxy[j]);

                t_xyy_xxy[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xy[j] * fl1_fz * fl2_fx + 8.0 * fl2_fx * pa_y[j] * fl1_fz * pb_x[j] - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_y[j] - fl1_fx * fl1_fz * fl1_fga * pb_xy[j] - pa_xyy[j] * fl1_fz * fl1_fgb * pb_y[j] + 2.0 * pa_x[j] * fl1_fz * fl2_fx * pb_y[j] + 4.0 * fl2_fx * fl1_fz * pb_xy[j] + 5.0 * pa_xyy[j] * fl1_fz * fl1_fx * pb_y[j] + 10.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xx[j] + 10.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxy[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_xxy[j]);

                t_xyy_xxz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx * pb_z[j] + 0.5 * fl2_fx * pb_xz[j] + 0.5 * pa_xyy[j] * fl1_fx * pb_z[j] + fl1_fx * pa_yy[j] * pb_xz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxz[j] + pa_xyy[j] * pb_xxz[j]);

                t_xyy_xxz[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - fl1_fx * fl1_fz * fl1_fga * pb_xz[j] - pa_xyy[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * pa_x[j] * fl1_fz * fl2_fx * pb_z[j] + 4.0 * fl2_fx * fl1_fz * pb_xz[j] + 5.0 * pa_xyy[j] * fl1_fz * fl1_fx * pb_z[j] + 10.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_xz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_xxz[j]);

                t_xyy_xyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_x[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pa_yy[j] + fl2_fx * pa_y[j] * pb_y[j] + 0.25 * fl2_fx * pb_yy[j] + 0.5 * pa_xyy[j] * pb_x[j] * fl1_fx + 2.0 * pa_xy[j] * fl1_fx * pb_xy[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyy[j] + pa_xyy[j] * pb_xyy[j]);

                t_xyy_xyy[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * fl1_fx * pa_yy[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl2_fx * fl1_fz * pb_x[j] + 2.0 * fl2_fx * pa_yy[j] * fl1_fz + 8.0 * fl2_fx * pa_y[j] * fl1_fz * pb_y[j] - 0.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yy[j] - pa_xyy[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pb_yy[j] + 5.0 * pa_xyy[j] * fl1_fz * pb_x[j] * fl1_fx + 20.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xy[j] + 5.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_yy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyy[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_xyy[j]);

                t_xyy_xyz[j] = fl_s_0_0 * (0.5 * fl2_fx * pa_y[j] * pb_z[j] + 0.25 * fl2_fx * pb_yz[j] + pa_xy[j] * fl1_fx * pb_xz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyz[j] + pa_xyy[j] * pb_xyz[j]);

                t_xyy_xyz[j] += fl_r_0_0 * (4.0 * fl2_fx * pa_y[j] * fl1_fz * pb_z[j] - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yz[j] + 2.0 * fl2_fx * fl1_fz * pb_yz[j] + 10.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xz[j] + 5.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_yz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_xyz[j]);

                t_xyy_xzz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * fl2_fx * pa_yy[j] + 0.25 * pa_x[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pb_zz[j] + 0.5 * pa_xyy[j] * pb_x[j] * fl1_fx + 0.5 * fl1_fx * pa_yy[j] * pb_zz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xzz[j] + pa_xyy[j] * pb_xzz[j]);

                t_xyy_xzz[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * fl1_fx * pa_yy[j] * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * fl2_fx * pa_yy[j] * fl1_fz - 0.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_zz[j] - pa_xyy[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl1_fz * fl2_fx * pb_x[j] + 2.0 * fl2_fx * fl1_fz * pb_zz[j] + 5.0 * pa_xyy[j] * fl1_fz * pb_x[j] * fl1_fx + 5.0 * fl1_fx * pa_yy[j] * fl1_fz * pb_zz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xzz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xzz[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_xzz[j]);

                t_xyy_yyy[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl2_fx + 2.25 * pa_x[j] * fl2_fx * pb_y[j] + 1.5 * pa_xyy[j] * pb_y[j] * fl1_fx + 3.0 * pa_xy[j] * fl1_fx * pb_yy[j] + 0.5 * pa_x[j] * fl1_fx * pb_yyy[j] + pa_xyy[j] * pb_yyy[j]);

                t_xyy_yyy[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_xy[j] * fl1_fz * fl2_fx + 18.0 * pa_x[j] * fl2_fx * fl1_fz * pb_y[j] - 1.5 * pa_x[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 3.0 * pa_xyy[j] * pb_y[j] * fl1_fz * fl1_fgb + 15.0 * pa_xyy[j] * fl1_fz * pb_y[j] * fl1_fx + 30.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yyy[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyy[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_yyy[j]);

                t_xyy_yyz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx * pb_z[j] + 0.5 * pa_xyy[j] * fl1_fx * pb_z[j] + 2.0 * pa_xy[j] * fl1_fx * pb_yz[j] + 0.5 * pa_x[j] * fl1_fx * pb_yyz[j] + pa_xyy[j] * pb_yyz[j]);

                t_xyy_yyz[j] += fl_r_0_0 * (6.0 * pa_x[j] * fl2_fx * fl1_fz * pb_z[j] - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - pa_xyy[j] * fl1_fz * fl1_fgb * pb_z[j] + 5.0 * pa_xyy[j] * fl1_fz * fl1_fx * pb_z[j] + 20.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yyz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyz[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_yyz[j]);

                t_xyy_yzz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + 0.25 * pa_x[j] * fl2_fx * pb_y[j] + 0.5 * pa_xyy[j] * pb_y[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_zz[j] + 0.5 * pa_x[j] * fl1_fx * pb_yzz[j] + pa_xyy[j] * pb_yzz[j]);

                t_xyy_yzz[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xy[j] * fl1_fz * fl2_fx - 0.5 * pa_x[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - pa_xyy[j] * pb_y[j] * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl1_fz * fl2_fx * pb_y[j] + 5.0 * pa_xyy[j] * fl1_fz * pb_y[j] * fl1_fx + 10.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_zz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yzz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yzz[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_yzz[j]);

                t_xyy_zzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx * pb_z[j] + 1.5 * pa_xyy[j] * pb_z[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_zzz[j] + pa_xyy[j] * pb_zzz[j]);

                t_xyy_zzz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_z[j] * fl1_fx - 3.0 * pa_xyy[j] * pb_z[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl1_fz * fl2_fx * pb_z[j] + 15.0 * pa_xyy[j] * fl1_fz * pb_z[j] * fl1_fx - pa_x[j] * fl1_fz * fl1_fga * pb_zzz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_zzz[j] + 12.0 * pa_xyy[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_40_50(      CMemBlock2D<double>& primBuffer,
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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyz_xxx = primBuffer.data(100 * idx + 40);

            auto t_xyz_xxy = primBuffer.data(100 * idx + 41);

            auto t_xyz_xxz = primBuffer.data(100 * idx + 42);

            auto t_xyz_xyy = primBuffer.data(100 * idx + 43);

            auto t_xyz_xyz = primBuffer.data(100 * idx + 44);

            auto t_xyz_xzz = primBuffer.data(100 * idx + 45);

            auto t_xyz_yyy = primBuffer.data(100 * idx + 46);

            auto t_xyz_yyz = primBuffer.data(100 * idx + 47);

            auto t_xyz_yzz = primBuffer.data(100 * idx + 48);

            auto t_xyz_zzz = primBuffer.data(100 * idx + 49);

            // Batch of Integrals (40,50)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyz_xxx, t_xyz_xxy, t_xyz_xxz, \
                                     t_xyz_xyy, t_xyz_xyz, t_xyz_xzz, t_xyz_yyy, t_xyz_yyz, t_xyz_yzz, t_xyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyz_xxx[j] = fl_s_0_0 * (0.75 * fl2_fx * pa_yz[j] + 1.5 * pa_xyz[j] * pb_x[j] * fl1_fx + 1.5 * fl1_fx * pa_yz[j] * pb_xx[j] + pa_xyz[j] * pb_xxx[j]);

                t_xyz_xxx[j] += fl_r_0_0 * (-1.5 * fl1_fx * pa_yz[j] * fl1_fz * fl1_fgb + 6.0 * fl2_fx * pa_yz[j] * fl1_fz - 3.0 * pa_xyz[j] * pb_x[j] * fl1_fz * fl1_fgb + 15.0 * pa_xyz[j] * fl1_fz * pb_x[j] * fl1_fx + 15.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xx[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_xxx[j]);

                t_xyz_xxy[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * fl2_fx * pa_z[j] * pb_x[j] + 0.5 * pa_xyz[j] * fl1_fx * pb_y[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xx[j] + fl1_fx * pa_yz[j] * pb_xy[j] + pa_xyz[j] * pb_xxy[j]);

                t_xyz_xxy[j] += fl_r_0_0 * (-0.5 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xz[j] * fl2_fx * fl1_fz + 4.0 * fl2_fx * fl1_fz * pa_z[j] * pb_x[j] - pa_xyz[j] * fl1_fz * fl1_fgb * pb_y[j] + 5.0 * pa_xyz[j] * fl1_fz * fl1_fx * pb_y[j] + 5.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xx[j] + 10.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xy[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_xxy[j]);

                t_xyz_xxz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * fl2_fx * pa_y[j] * pb_x[j] + 0.5 * pa_xyz[j] * fl1_fx * pb_z[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xx[j] + fl1_fx * pa_yz[j] * pb_xz[j] + pa_xyz[j] * pb_xxz[j]);

                t_xyz_xxz[j] += fl_r_0_0 * (-0.5 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xy[j] * fl1_fz * fl2_fx + 4.0 * fl2_fx * pa_y[j] * fl1_fz * pb_x[j] - pa_xyz[j] * fl1_fz * fl1_fgb * pb_z[j] + 5.0 * pa_xyz[j] * fl1_fz * fl1_fx * pb_z[j] + 5.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xx[j] + 10.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_xz[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_xxz[j]);

                t_xyz_xyy[j] = fl_s_0_0 * (0.25 * fl2_fx * pa_yz[j] + 0.5 * fl2_fx * pa_z[j] * pb_y[j] + 0.5 * pa_xyz[j] * pb_x[j] * fl1_fx + pa_xz[j] * fl1_fx * pb_xy[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yy[j] + pa_xyz[j] * pb_xyy[j]);

                t_xyz_xyy[j] += fl_r_0_0 * (-0.5 * fl1_fx * pa_yz[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * pa_yz[j] * fl1_fz + 4.0 * fl2_fx * fl1_fz * pa_z[j] * pb_y[j] - pa_xyz[j] * pb_x[j] * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl1_fz * pb_x[j] * fl1_fx + 10.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xy[j] + 5.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_yy[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_xyy[j]);

                t_xyz_xyz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_x[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pa_y[j] * pb_y[j] + 0.25 * fl2_fx * pa_z[j] * pb_z[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xy[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yz[j] + pa_xyz[j] * pb_xyz[j]);

                t_xyz_xyz[j] += fl_r_0_0 * (0.75 * fl3_fx * fl1_fz + 2.0 * pa_x[j] * fl2_fx * fl1_fz * pb_x[j] + 2.0 * fl2_fx * pa_y[j] * fl1_fz * pb_y[j] + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_z[j] + 5.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xy[j] + 5.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_xz[j] + 5.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_yz[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_xyz[j]);

                t_xyz_xzz[j] = fl_s_0_0 * (0.25 * fl2_fx * pa_yz[j] + 0.5 * fl2_fx * pa_y[j] * pb_z[j] + 0.5 * pa_xyz[j] * pb_x[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_xz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_zz[j] + pa_xyz[j] * pb_xzz[j]);

                t_xyz_xzz[j] += fl_r_0_0 * (-0.5 * fl1_fx * pa_yz[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * pa_yz[j] * fl1_fz + 4.0 * fl2_fx * pa_y[j] * fl1_fz * pb_z[j] - pa_xyz[j] * pb_x[j] * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl1_fz * pb_x[j] * fl1_fx + 10.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_xz[j] + 5.0 * fl1_fx * pa_yz[j] * fl1_fz * pb_zz[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_xzz[j]);

                t_xyz_yyy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 1.5 * pa_xyz[j] * pb_y[j] * fl1_fx + 1.5 * pa_xz[j] * fl1_fx * pb_yy[j] + pa_xyz[j] * pb_yyy[j]);

                t_xyz_yyy[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_xz[j] * fl2_fx * fl1_fz - 3.0 * pa_xyz[j] * pb_y[j] * fl1_fz * fl1_fgb + 15.0 * pa_xyz[j] * fl1_fz * pb_y[j] * fl1_fx + 15.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yy[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_yyy[j]);

                t_xyz_yyz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * pa_x[j] * fl2_fx * pb_y[j] + 0.5 * pa_xyz[j] * fl1_fx * pb_z[j] + 0.5 * pa_xy[j] * fl1_fx * pb_yy[j] + pa_xz[j] * fl1_fx * pb_yz[j] + pa_xyz[j] * pb_yyz[j]);

                t_xyz_yyz[j] += fl_r_0_0 * (-0.5 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xy[j] * fl1_fz * fl2_fx + 4.0 * pa_x[j] * fl2_fx * fl1_fz * pb_y[j] - pa_xyz[j] * fl1_fz * fl1_fgb * pb_z[j] + 5.0 * pa_xyz[j] * fl1_fz * fl1_fx * pb_z[j] + 5.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yy[j] + 10.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_yz[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_yyz[j]);

                t_xyz_yzz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * pa_x[j] * fl2_fx * pb_z[j] + 0.5 * pa_xyz[j] * pb_y[j] * fl1_fx + pa_xy[j] * fl1_fx * pb_yz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_zz[j] + pa_xyz[j] * pb_yzz[j]);

                t_xyz_yzz[j] += fl_r_0_0 * (-0.5 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xz[j] * fl2_fx * fl1_fz + 4.0 * pa_x[j] * fl2_fx * fl1_fz * pb_z[j] - pa_xyz[j] * pb_y[j] * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl1_fz * pb_y[j] * fl1_fx + 10.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_yz[j] + 5.0 * pa_xz[j] * fl1_fx * fl1_fz * pb_zz[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_yzz[j]);

                t_xyz_zzz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 1.5 * pa_xyz[j] * pb_z[j] * fl1_fx + 1.5 * pa_xy[j] * fl1_fx * pb_zz[j] + pa_xyz[j] * pb_zzz[j]);

                t_xyz_zzz[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_xy[j] * fl1_fz * fl2_fx - 3.0 * pa_xyz[j] * pb_z[j] * fl1_fz * fl1_fgb + 15.0 * pa_xyz[j] * fl1_fz * pb_z[j] * fl1_fx + 15.0 * pa_xy[j] * fl1_fz * fl1_fx * pb_zz[j] + 12.0 * pa_xyz[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_50_60(      CMemBlock2D<double>& primBuffer,
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

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xzz = paDistances.data(19 * idx + 14);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzz_xxx = primBuffer.data(100 * idx + 50);

            auto t_xzz_xxy = primBuffer.data(100 * idx + 51);

            auto t_xzz_xxz = primBuffer.data(100 * idx + 52);

            auto t_xzz_xyy = primBuffer.data(100 * idx + 53);

            auto t_xzz_xyz = primBuffer.data(100 * idx + 54);

            auto t_xzz_xzz = primBuffer.data(100 * idx + 55);

            auto t_xzz_yyy = primBuffer.data(100 * idx + 56);

            auto t_xzz_yyz = primBuffer.data(100 * idx + 57);

            auto t_xzz_yzz = primBuffer.data(100 * idx + 58);

            auto t_xzz_zzz = primBuffer.data(100 * idx + 59);

            // Batch of Integrals (50,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xzz_xxx, t_xzz_xxy, t_xzz_xxz, t_xzz_xyy, \
                                     t_xzz_xyz, t_xzz_xzz, t_xzz_yyy, t_xzz_yyz, t_xzz_yzz, t_xzz_zzz: VLX_ALIGN)
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

                t_xzz_xxx[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * fl2_fx * pa_zz[j] + 0.75 * pa_x[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pb_xx[j] + 1.5 * pa_xzz[j] * pb_x[j] * fl1_fx + 1.5 * fl1_fx * pa_zz[j] * pb_xx[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxx[j] + pa_xzz[j] * pb_xxx[j]);

                t_xzz_xxx[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * fl2_fx * pa_zz[j] * fl1_fz - 1.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xx[j] - 3.0 * pa_xzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl1_fz * fl2_fx * pb_x[j] + 6.0 * fl2_fx * fl1_fz * pb_xx[j] + 15.0 * pa_xzz[j] * fl1_fz * pb_x[j] * fl1_fx + 15.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xx[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxx[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxx[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_xxx[j]);

                t_xzz_xxy[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl2_fx * pb_y[j] + 0.5 * fl2_fx * pb_xy[j] + 0.5 * pa_xzz[j] * fl1_fx * pb_y[j] + fl1_fx * pa_zz[j] * pb_xy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxy[j] + pa_xzz[j] * pb_xxy[j]);

                t_xzz_xxy[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_y[j] - fl1_fx * fl1_fz * fl1_fga * pb_xy[j] - pa_xzz[j] * fl1_fz * fl1_fgb * pb_y[j] + 2.0 * pa_x[j] * fl1_fz * fl2_fx * pb_y[j] + 4.0 * fl2_fx * fl1_fz * pb_xy[j] + 5.0 * pa_xzz[j] * fl1_fz * fl1_fx * pb_y[j] + 10.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxy[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_xxy[j]);

                t_xzz_xxz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + fl2_fx * pa_z[j] * pb_x[j] + 0.25 * pa_x[j] * fl2_fx * pb_z[j] + 0.5 * fl2_fx * pb_xz[j] + 0.5 * pa_xzz[j] * fl1_fx * pb_z[j] + pa_xz[j] * fl1_fx * pb_xx[j] + fl1_fx * pa_zz[j] * pb_xz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xxz[j] + pa_xzz[j] * pb_xxz[j]);

                t_xzz_xxz[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xz[j] * fl1_fz * fl2_fx + 8.0 * fl2_fx * pa_z[j] * fl1_fz * pb_x[j] - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - fl1_fx * fl1_fz * fl1_fga * pb_xz[j] - pa_xzz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * pa_x[j] * fl1_fz * fl2_fx * pb_z[j] + 4.0 * fl2_fx * fl1_fz * pb_xz[j] + 5.0 * pa_xzz[j] * fl1_fz * fl1_fx * pb_z[j] + 10.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xx[j] + 10.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xxz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_xxz[j]);

                t_xzz_xyy[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * fl2_fx * pa_zz[j] + 0.25 * pa_x[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pb_yy[j] + 0.5 * pa_xzz[j] * pb_x[j] * fl1_fx + 0.5 * fl1_fx * pa_zz[j] * pb_yy[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyy[j] + pa_xzz[j] * pb_xyy[j]);

                t_xzz_xyy[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * fl2_fx * pa_zz[j] * fl1_fz - 0.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yy[j] - pa_xzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl1_fz * fl2_fx * pb_x[j] + 2.0 * fl2_fx * fl1_fz * pb_yy[j] + 5.0 * pa_xzz[j] * fl1_fz * pb_x[j] * fl1_fx + 5.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyy[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_xyy[j]);

                t_xzz_xyz[j] = fl_s_0_0 * (0.5 * fl2_fx * pa_z[j] * pb_y[j] + 0.25 * fl2_fx * pb_yz[j] + pa_xz[j] * fl1_fx * pb_xy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xyz[j] + pa_xzz[j] * pb_xyz[j]);

                t_xzz_xyz[j] += fl_r_0_0 * (4.0 * fl2_fx * pa_z[j] * fl1_fz * pb_y[j] - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_yz[j] + 2.0 * fl2_fx * fl1_fz * pb_yz[j] + 10.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xy[j] + 5.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xyz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_xyz[j]);

                t_xzz_xzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_x[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pa_zz[j] + fl2_fx * pa_z[j] * pb_z[j] + 0.25 * fl2_fx * pb_zz[j] + 0.5 * pa_xzz[j] * pb_x[j] * fl1_fx + 2.0 * pa_xz[j] * fl1_fx * pb_xz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_zz[j] + 0.5 * pa_x[j] * fl1_fx * pb_xzz[j] + pa_xzz[j] * pb_xzz[j]);

                t_xzz_xzz[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl2_fx * fl1_fz * pb_x[j] + 2.0 * fl2_fx * pa_zz[j] * fl1_fz + 8.0 * fl2_fx * pa_z[j] * fl1_fz * pb_z[j] - 0.5 * pa_x[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_zz[j] - pa_xzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pb_zz[j] + 5.0 * pa_xzz[j] * fl1_fz * pb_x[j] * fl1_fx + 20.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_xz[j] + 5.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_zz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_xzz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_xzz[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_xzz[j]);

                t_xzz_yyy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx * pb_y[j] + 1.5 * pa_xzz[j] * pb_y[j] * fl1_fx + 0.5 * pa_x[j] * fl1_fx * pb_yyy[j] + pa_xzz[j] * pb_yyy[j]);

                t_xzz_yyy[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 3.0 * pa_xzz[j] * pb_y[j] * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl1_fz * fl2_fx * pb_y[j] + 15.0 * pa_xzz[j] * fl1_fz * pb_y[j] * fl1_fx - pa_x[j] * fl1_fz * fl1_fga * pb_yyy[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyy[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_yyy[j]);

                t_xzz_yyz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + 0.25 * pa_x[j] * fl2_fx * pb_z[j] + 0.5 * pa_xzz[j] * fl1_fx * pb_z[j] + pa_xz[j] * fl1_fx * pb_yy[j] + 0.5 * pa_x[j] * fl1_fx * pb_yyz[j] + pa_xzz[j] * pb_yyz[j]);

                t_xzz_yyz[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xz[j] * fl1_fz * fl2_fx - 0.5 * pa_x[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * pa_x[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - pa_xzz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * pa_x[j] * fl1_fz * fl2_fx * pb_z[j] + 5.0 * pa_xzz[j] * fl1_fz * fl1_fx * pb_z[j] + 10.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_yy[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yyz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yyz[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_yyz[j]);

                t_xzz_yzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl2_fx * pb_y[j] + 0.5 * pa_xzz[j] * pb_y[j] * fl1_fx + 2.0 * pa_xz[j] * fl1_fx * pb_yz[j] + 0.5 * pa_x[j] * fl1_fx * pb_yzz[j] + pa_xzz[j] * pb_yzz[j]);

                t_xzz_yzz[j] += fl_r_0_0 * (6.0 * pa_x[j] * fl2_fx * fl1_fz * pb_y[j] - 0.5 * pa_x[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - pa_xzz[j] * pb_y[j] * fl1_fz * fl1_fgb + 5.0 * pa_xzz[j] * fl1_fz * pb_y[j] * fl1_fx + 20.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_yz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_yzz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_yzz[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_yzz[j]);

                t_xzz_zzz[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl2_fx + 2.25 * pa_x[j] * fl2_fx * pb_z[j] + 1.5 * pa_xzz[j] * pb_z[j] * fl1_fx + 3.0 * pa_xz[j] * fl1_fx * pb_zz[j] + 0.5 * pa_x[j] * fl1_fx * pb_zzz[j] + pa_xzz[j] * pb_zzz[j]);

                t_xzz_zzz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_xz[j] * fl1_fz * fl2_fx + 18.0 * pa_x[j] * fl2_fx * fl1_fz * pb_z[j] - 1.5 * pa_x[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl1_fz * fl1_fga * pb_z[j] * fl1_fx - 3.0 * pa_xzz[j] * pb_z[j] * fl1_fz * fl1_fgb + 15.0 * pa_xzz[j] * fl1_fz * pb_z[j] * fl1_fx + 30.0 * pa_xz[j] * fl1_fz * fl1_fx * pb_zz[j] - pa_x[j] * fl1_fz * fl1_fga * pb_zzz[j] + 5.0 * pa_x[j] * fl1_fz * fl1_fx * pb_zzz[j] + 12.0 * pa_xzz[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_60_70(      CMemBlock2D<double>& primBuffer,
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

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(19 * idx + 1);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(19 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(19 * idx + 15);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyy_xxx = primBuffer.data(100 * idx + 60);

            auto t_yyy_xxy = primBuffer.data(100 * idx + 61);

            auto t_yyy_xxz = primBuffer.data(100 * idx + 62);

            auto t_yyy_xyy = primBuffer.data(100 * idx + 63);

            auto t_yyy_xyz = primBuffer.data(100 * idx + 64);

            auto t_yyy_xzz = primBuffer.data(100 * idx + 65);

            auto t_yyy_yyy = primBuffer.data(100 * idx + 66);

            auto t_yyy_yyz = primBuffer.data(100 * idx + 67);

            auto t_yyy_yzz = primBuffer.data(100 * idx + 68);

            auto t_yyy_zzz = primBuffer.data(100 * idx + 69);

            // Batch of Integrals (60,70)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, r_0_0, s_0_0, t_yyy_xxx, t_yyy_xxy, t_yyy_xxz, t_yyy_xyy, t_yyy_xyz, \
                                     t_yyy_xzz, t_yyy_yyy, t_yyy_yyz, t_yyy_yzz, t_yyy_zzz: VLX_ALIGN)
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

                t_yyy_xxx[j] = fl_s_0_0 * (2.25 * pa_y[j] * fl2_fx * pb_x[j] + 1.5 * pa_yyy[j] * pb_x[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pb_xxx[j] + pa_yyy[j] * pb_xxx[j]);

                t_yyy_xxx[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 3.0 * pa_yyy[j] * pb_x[j] * fl1_fz * fl1_fgb + 18.0 * pa_y[j] * fl1_fz * fl2_fx * pb_x[j] + 15.0 * pa_yyy[j] * fl1_fz * pb_x[j] * fl1_fx - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxx[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxx[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_xxx[j]);

                t_yyy_xxy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * pa_y[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pb_xx[j] + 0.5 * pa_yyy[j] * fl1_fx * pb_y[j] + 1.5 * pa_yy[j] * fl1_fx * pb_xx[j] + 1.5 * pa_y[j] * fl1_fx * pb_xxy[j] + pa_yyy[j] * pb_xxy[j]);

                t_yyy_xxy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_yy[j] * fl1_fz * fl2_fx - 1.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 1.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_y[j] - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xx[j] - pa_yyy[j] * fl1_fz * fl1_fgb * pb_y[j] + 6.0 * pa_y[j] * fl1_fz * fl2_fx * pb_y[j] + 6.0 * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * pa_yyy[j] * fl1_fz * fl1_fx * pb_y[j] + 15.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xx[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxy[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_xxy[j]);

                t_yyy_xxz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx * pb_z[j] + 0.5 * pa_yyy[j] * fl1_fx * pb_z[j] + 1.5 * pa_y[j] * fl1_fx * pb_xxz[j] + pa_yyy[j] * pb_xxz[j]);

                t_yyy_xxz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 1.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - pa_yyy[j] * fl1_fz * fl1_fgb * pb_z[j] + 6.0 * pa_y[j] * fl1_fz * fl2_fx * pb_z[j] + 5.0 * pa_yyy[j] * fl1_fz * fl1_fx * pb_z[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xxz[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_xxz[j]);

                t_yyy_xyy[j] = fl_s_0_0 * (2.25 * pa_y[j] * fl2_fx * pb_x[j] + 1.5 * fl2_fx * pb_xy[j] + 0.5 * pa_yyy[j] * pb_x[j] * fl1_fx + 3.0 * pa_yy[j] * fl1_fx * pb_xy[j] + 1.5 * pa_y[j] * fl1_fx * pb_xyy[j] + pa_yyy[j] * pb_xyy[j]);

                t_yyy_xyy[j] += fl_r_0_0 * (18.0 * pa_y[j] * fl2_fx * fl1_fz * pb_x[j] - 1.5 * pa_y[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xy[j] - pa_yyy[j] * pb_x[j] * fl1_fz * fl1_fgb + 12.0 * fl2_fx * fl1_fz * pb_xy[j] + 5.0 * pa_yyy[j] * fl1_fz * pb_x[j] * fl1_fx + 30.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xy[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xyy[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_xyy[j]);

                t_yyy_xyz[j] = fl_s_0_0 * (0.75 * fl2_fx * pb_xz[j] + 1.5 * pa_yy[j] * fl1_fx * pb_xz[j] + 1.5 * pa_y[j] * fl1_fx * pb_xyz[j] + pa_yyy[j] * pb_xyz[j]);

                t_yyy_xyz[j] += fl_r_0_0 * (-1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xz[j] + 6.0 * fl2_fx * fl1_fz * pb_xz[j] + 15.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xyz[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_xyz[j]);

                t_yyy_xzz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx * pb_x[j] + 0.5 * pa_yyy[j] * pb_x[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pb_xzz[j] + pa_yyy[j] * pb_xzz[j]);

                t_yyy_xzz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - pa_yyy[j] * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl1_fz * fl2_fx * pb_x[j] + 5.0 * pa_yyy[j] * fl1_fz * pb_x[j] * fl1_fx - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_xzz[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xzz[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_xzz[j]);

                t_yyy_yyy[j] = fl_s_0_0 * (1.875 * fl3_fx + 2.25 * pa_yy[j] * fl2_fx + 6.75 * pa_y[j] * fl2_fx * pb_y[j] + 2.25 * fl2_fx * pb_yy[j] + 1.5 * pa_yyy[j] * pb_y[j] * fl1_fx + 4.5 * pa_yy[j] * fl1_fx * pb_yy[j] + 1.5 * pa_y[j] * fl1_fx * pb_yyy[j] + pa_yyy[j] * pb_yyy[j]);

                t_yyy_yyy[j] += fl_r_0_0 * (11.25 * fl3_fx * fl1_fz - 2.25 * fl2_fx * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa_yy[j] * fl1_fz * fl2_fx + 54.0 * pa_y[j] * fl2_fx * fl1_fz * pb_y[j] - 4.5 * pa_y[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_yy[j] - 3.0 * pa_yyy[j] * pb_y[j] * fl1_fz * fl1_fgb + 18.0 * fl2_fx * fl1_fz * pb_yy[j] + 15.0 * pa_yyy[j] * fl1_fz * pb_y[j] * fl1_fx + 45.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yy[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yyy[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyy[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_yyy[j]);

                t_yyy_yyz[j] = fl_s_0_0 * (2.25 * pa_y[j] * fl2_fx * pb_z[j] + 1.5 * fl2_fx * pb_yz[j] + 0.5 * pa_yyy[j] * fl1_fx * pb_z[j] + 3.0 * pa_yy[j] * fl1_fx * pb_yz[j] + 1.5 * pa_y[j] * fl1_fx * pb_yyz[j] + pa_yyy[j] * pb_yyz[j]);

                t_yyy_yyz[j] += fl_r_0_0 * (18.0 * pa_y[j] * fl2_fx * fl1_fz * pb_z[j] - 1.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 1.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_yz[j] - pa_yyy[j] * fl1_fz * fl1_fgb * pb_z[j] + 12.0 * fl2_fx * fl1_fz * pb_yz[j] + 5.0 * pa_yyy[j] * fl1_fz * fl1_fx * pb_z[j] + 30.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yyz[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyz[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_yyz[j]);

                t_yyy_yzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * pa_y[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pb_zz[j] + 0.5 * pa_yyy[j] * pb_y[j] * fl1_fx + 1.5 * pa_yy[j] * fl1_fx * pb_zz[j] + 1.5 * pa_y[j] * fl1_fx * pb_yzz[j] + pa_yyy[j] * pb_yzz[j]);

                t_yyy_yzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_yy[j] * fl1_fz * fl2_fx - 1.5 * pa_y[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_zz[j] - pa_yyy[j] * pb_y[j] * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl1_fz * fl2_fx * pb_y[j] + 6.0 * fl2_fx * fl1_fz * pb_zz[j] + 5.0 * pa_yyy[j] * fl1_fz * pb_y[j] * fl1_fx + 15.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_zz[j] - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_yzz[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yzz[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_yzz[j]);

                t_yyy_zzz[j] = fl_s_0_0 * (2.25 * pa_y[j] * fl2_fx * pb_z[j] + 1.5 * pa_yyy[j] * pb_z[j] * fl1_fx + 1.5 * pa_y[j] * fl1_fx * pb_zzz[j] + pa_yyy[j] * pb_zzz[j]);

                t_yyy_zzz[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb - 4.5 * pa_y[j] * fl1_fz * fl1_fga * pb_z[j] * fl1_fx - 3.0 * pa_yyy[j] * pb_z[j] * fl1_fz * fl1_fgb + 18.0 * pa_y[j] * fl1_fz * fl2_fx * pb_z[j] + 15.0 * pa_yyy[j] * fl1_fz * pb_z[j] * fl1_fx - 3.0 * pa_y[j] * fl1_fz * fl1_fga * pb_zzz[j] + 15.0 * pa_y[j] * fl1_fz * fl1_fx * pb_zzz[j] + 12.0 * pa_yyy[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_70_80(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyz = paDistances.data(19 * idx + 16);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyz_xxx = primBuffer.data(100 * idx + 70);

            auto t_yyz_xxy = primBuffer.data(100 * idx + 71);

            auto t_yyz_xxz = primBuffer.data(100 * idx + 72);

            auto t_yyz_xyy = primBuffer.data(100 * idx + 73);

            auto t_yyz_xyz = primBuffer.data(100 * idx + 74);

            auto t_yyz_xzz = primBuffer.data(100 * idx + 75);

            auto t_yyz_yyy = primBuffer.data(100 * idx + 76);

            auto t_yyz_yyz = primBuffer.data(100 * idx + 77);

            auto t_yyz_yzz = primBuffer.data(100 * idx + 78);

            auto t_yyz_zzz = primBuffer.data(100 * idx + 79);

            // Batch of Integrals (70,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_yyz_xxx, t_yyz_xxy, t_yyz_xxz, t_yyz_xyy, \
                                     t_yyz_xyz, t_yyz_xzz, t_yyz_yyy, t_yyz_yyz, t_yyz_yzz, t_yyz_zzz: VLX_ALIGN)
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

                t_yyz_xxx[j] = fl_s_0_0 * (0.75 * fl2_fx * pa_z[j] * pb_x[j] + 1.5 * pa_yyz[j] * pb_x[j] * fl1_fx + 0.5 * fl1_fx * pa_z[j] * pb_xxx[j] + pa_yyz[j] * pb_xxx[j]);

                t_yyz_xxx[j] += fl_r_0_0 * (-1.5 * fl1_fx * pa_z[j] * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_x[j] * fl1_fx - 3.0 * pa_yyz[j] * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * fl2_fx * fl1_fz * pa_z[j] * pb_x[j] + 15.0 * pa_yyz[j] * fl1_fz * pb_x[j] * fl1_fx - fl1_fz * fl1_fga * pa_z[j] * pb_xxx[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxx[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_xxx[j]);

                t_yyz_xxy[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + 0.25 * fl2_fx * pa_z[j] * pb_y[j] + 0.5 * pa_yyz[j] * fl1_fx * pb_y[j] + pa_yz[j] * fl1_fx * pb_xx[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxy[j] + pa_yyz[j] * pb_xxy[j]);

                t_yyz_xxy[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_yz[j] * fl2_fx * fl1_fz - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_y[j] - pa_yyz[j] * fl1_fz * fl1_fgb * pb_y[j] + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_y[j] + 5.0 * pa_yyz[j] * fl1_fz * fl1_fx * pb_y[j] + 10.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xx[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxy[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxy[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_xxy[j]);

                t_yyz_xxz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + 0.25 * fl2_fx * pa_z[j] * pb_z[j] + 0.25 * fl2_fx * pb_xx[j] + 0.5 * pa_yyz[j] * fl1_fx * pb_z[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xx[j] + 0.5 * fl1_fx * pa_z[j] * pb_xxz[j] + pa_yyz[j] * pb_xxz[j]);

                t_yyz_xxz[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_yy[j] * fl1_fz * fl2_fx - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xx[j] - pa_yyz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_z[j] + 2.0 * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * pa_yyz[j] * fl1_fz * fl1_fx * pb_z[j] + 5.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xx[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xxz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xxz[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_xxz[j]);

                t_yyz_xyy[j] = fl_s_0_0 * (0.75 * fl2_fx * pa_z[j] * pb_x[j] + 0.5 * pa_yyz[j] * pb_x[j] * fl1_fx + 2.0 * pa_yz[j] * fl1_fx * pb_xy[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyy[j] + pa_yyz[j] * pb_xyy[j]);

                t_yyz_xyy[j] += fl_r_0_0 * (6.0 * fl2_fx * fl1_fz * pa_z[j] * pb_x[j] - 0.5 * fl1_fx * pa_z[j] * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_x[j] * fl1_fx - pa_yyz[j] * pb_x[j] * fl1_fz * fl1_fgb + 5.0 * pa_yyz[j] * fl1_fz * pb_x[j] * fl1_fx + 20.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyy[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyy[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_xyy[j]);

                t_yyz_xyz[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pb_xy[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xy[j] + pa_yz[j] * fl1_fx * pb_xz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xyz[j] + pa_yyz[j] * pb_xyz[j]);

                t_yyz_xyz[j] += fl_r_0_0 * (4.0 * pa_y[j] * fl2_fx * fl1_fz * pb_x[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_xy[j] + 2.0 * fl2_fx * fl1_fz * pb_xy[j] + 5.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xy[j] + 10.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_xz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xyz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xyz[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_xyz[j]);

                t_yyz_xzz[j] = fl_s_0_0 * (0.25 * fl2_fx * pa_z[j] * pb_x[j] + 0.5 * fl2_fx * pb_xz[j] + 0.5 * pa_yyz[j] * pb_x[j] * fl1_fx + pa_yy[j] * fl1_fx * pb_xz[j] + 0.5 * fl1_fx * pa_z[j] * pb_xzz[j] + pa_yyz[j] * pb_xzz[j]);

                t_yyz_xzz[j] += fl_r_0_0 * (-0.5 * fl1_fx * pa_z[j] * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_x[j] * fl1_fx - fl1_fz * fl1_fga * fl1_fx * pb_xz[j] - pa_yyz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_x[j] + 4.0 * fl2_fx * fl1_fz * pb_xz[j] + 5.0 * pa_yyz[j] * fl1_fz * pb_x[j] * fl1_fx + 10.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_xz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_xzz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_xzz[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_xzz[j]);

                t_yyz_yyy[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl2_fx + 2.25 * fl2_fx * pa_z[j] * pb_y[j] + 1.5 * pa_yyz[j] * pb_y[j] * fl1_fx + 3.0 * pa_yz[j] * fl1_fx * pb_yy[j] + 0.5 * fl1_fx * pa_z[j] * pb_yyy[j] + pa_yyz[j] * pb_yyy[j]);

                t_yyz_yyy[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_yz[j] * fl2_fx * fl1_fz + 18.0 * fl2_fx * fl1_fz * pa_z[j] * pb_y[j] - 1.5 * fl1_fx * pa_z[j] * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_y[j] * fl1_fx - 3.0 * pa_yyz[j] * pb_y[j] * fl1_fz * fl1_fgb + 15.0 * pa_yyz[j] * fl1_fz * pb_y[j] * fl1_fx + 30.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_yy[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yyy[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyy[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_yyy[j]);

                t_yyz_yyz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + pa_y[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pa_z[j] * pb_z[j] + 0.25 * fl2_fx * pb_yy[j] + 0.5 * pa_yyz[j] * fl1_fx * pb_z[j] + 0.5 * pa_yy[j] * fl1_fx * pb_yy[j] + 2.0 * pa_yz[j] * fl1_fx * pb_yz[j] + 0.5 * fl1_fx * pa_z[j] * pb_yyz[j] + pa_yyz[j] * pb_yyz[j]);

                t_yyz_yyz[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_yy[j] * fl1_fz * fl2_fx + 8.0 * pa_y[j] * fl2_fx * fl1_fz * pb_y[j] + 6.0 * fl2_fx * fl1_fz * pa_z[j] * pb_z[j] - 0.5 * fl1_fx * pa_z[j] * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * fl1_fz * fl1_fga * pa_z[j] * fl1_fx * pb_z[j] - 0.5 * fl1_fz * fl1_fga * fl1_fx * pb_yy[j] - pa_yyz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * fl2_fx * fl1_fz * pb_yy[j] + 5.0 * pa_yyz[j] * fl1_fz * fl1_fx * pb_z[j] + 5.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yy[j] + 20.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_yz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yyz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yyz[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_yyz[j]);

                t_yyz_yzz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + pa_y[j] * fl2_fx * pb_z[j] + 0.25 * fl2_fx * pa_z[j] * pb_y[j] + 0.5 * fl2_fx * pb_yz[j] + 0.5 * pa_yyz[j] * pb_y[j] * fl1_fx + pa_yy[j] * fl1_fx * pb_yz[j] + pa_yz[j] * fl1_fx * pb_zz[j] + 0.5 * fl1_fx * pa_z[j] * pb_yzz[j] + pa_yyz[j] * pb_yzz[j]);

                t_yyz_yzz[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_yz[j] * fl2_fx * fl1_fz + 8.0 * pa_y[j] * fl2_fx * fl1_fz * pb_z[j] - 0.5 * fl1_fx * pa_z[j] * pb_y[j] * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * pa_z[j] * pb_y[j] * fl1_fx - fl1_fz * fl1_fga * fl1_fx * pb_yz[j] - pa_yyz[j] * pb_y[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pa_z[j] * pb_y[j] + 4.0 * fl2_fx * fl1_fz * pb_yz[j] + 5.0 * pa_yyz[j] * fl1_fz * pb_y[j] * fl1_fx + 10.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_yz[j] + 10.0 * pa_yz[j] * fl1_fx * fl1_fz * pb_zz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_yzz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_yzz[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_yzz[j]);

                t_yyz_zzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * fl2_fx * pa_z[j] * pb_z[j] + 0.75 * fl2_fx * pb_zz[j] + 1.5 * pa_yyz[j] * pb_z[j] * fl1_fx + 1.5 * pa_yy[j] * fl1_fx * pb_zz[j] + 0.5 * fl1_fx * pa_z[j] * pb_zzz[j] + pa_yyz[j] * pb_zzz[j]);

                t_yyz_zzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_yy[j] * fl1_fz * fl2_fx - 1.5 * fl1_fx * pa_z[j] * pb_z[j] * fl1_fz * fl1_fgb - 1.5 * fl1_fz * fl1_fga * pa_z[j] * pb_z[j] * fl1_fx - 1.5 * fl1_fz * fl1_fga * fl1_fx * pb_zz[j] - 3.0 * pa_yyz[j] * pb_z[j] * fl1_fz * fl1_fgb + 6.0 * fl2_fx * fl1_fz * pa_z[j] * pb_z[j] + 6.0 * fl2_fx * fl1_fz * pb_zz[j] + 15.0 * pa_yyz[j] * fl1_fz * pb_z[j] * fl1_fx + 15.0 * pa_yy[j] * fl1_fz * fl1_fx * pb_zz[j] - fl1_fz * fl1_fga * pa_z[j] * pb_zzz[j] + 5.0 * fl1_fx * fl1_fz * pa_z[j] * pb_zzz[j] + 12.0 * pa_yyz[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_80_90(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yzz = paDistances.data(19 * idx + 17);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzz_xxx = primBuffer.data(100 * idx + 80);

            auto t_yzz_xxy = primBuffer.data(100 * idx + 81);

            auto t_yzz_xxz = primBuffer.data(100 * idx + 82);

            auto t_yzz_xyy = primBuffer.data(100 * idx + 83);

            auto t_yzz_xyz = primBuffer.data(100 * idx + 84);

            auto t_yzz_xzz = primBuffer.data(100 * idx + 85);

            auto t_yzz_yyy = primBuffer.data(100 * idx + 86);

            auto t_yzz_yyz = primBuffer.data(100 * idx + 87);

            auto t_yzz_yzz = primBuffer.data(100 * idx + 88);

            auto t_yzz_zzz = primBuffer.data(100 * idx + 89);

            // Batch of Integrals (80,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_yzz_xxx, t_yzz_xxy, t_yzz_xxz, t_yzz_xyy, \
                                     t_yzz_xyz, t_yzz_xzz, t_yzz_yyy, t_yzz_yyz, t_yzz_yzz, t_yzz_zzz: VLX_ALIGN)
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

                t_yzz_xxx[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx * pb_x[j] + 1.5 * pa_yzz[j] * pb_x[j] * fl1_fx + 0.5 * pa_y[j] * fl1_fx * pb_xxx[j] + pa_yzz[j] * pb_xxx[j]);

                t_yzz_xxx[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 3.0 * pa_yzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl1_fz * fl2_fx * pb_x[j] + 15.0 * pa_yzz[j] * fl1_fz * pb_x[j] * fl1_fx - pa_y[j] * fl1_fz * fl1_fga * pb_xxx[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxx[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_xxx[j]);

                t_yzz_xxy[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * fl2_fx * pa_zz[j] + 0.25 * pa_y[j] * fl2_fx * pb_y[j] + 0.25 * fl2_fx * pb_xx[j] + 0.5 * pa_yzz[j] * fl1_fx * pb_y[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xx[j] + 0.5 * pa_y[j] * fl1_fx * pb_xxy[j] + pa_yzz[j] * pb_xxy[j]);

                t_yzz_xxy[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * fl2_fx * pa_zz[j] * fl1_fz - 0.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_y[j] - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_xx[j] - pa_yzz[j] * fl1_fz * fl1_fgb * pb_y[j] + 2.0 * pa_y[j] * fl1_fz * fl2_fx * pb_y[j] + 2.0 * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * pa_yzz[j] * fl1_fz * fl1_fx * pb_y[j] + 5.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xx[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xxy[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_xxy[j]);

                t_yzz_xxz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + 0.25 * pa_y[j] * fl2_fx * pb_z[j] + 0.5 * pa_yzz[j] * fl1_fx * pb_z[j] + pa_yz[j] * fl1_fx * pb_xx[j] + 0.5 * pa_y[j] * fl1_fx * pb_xxz[j] + pa_yzz[j] * pb_xxz[j]);

                t_yzz_xxz[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_yz[j] * fl1_fz * fl2_fx - 0.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - pa_yzz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * pa_y[j] * fl1_fz * fl2_fx * pb_z[j] + 5.0 * pa_yzz[j] * fl1_fz * fl1_fx * pb_z[j] + 10.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xx[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xxz[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_xxz[j]);

                t_yzz_xyy[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl2_fx * pb_x[j] + 0.5 * fl2_fx * pb_xy[j] + 0.5 * pa_yzz[j] * pb_x[j] * fl1_fx + fl1_fx * pa_zz[j] * pb_xy[j] + 0.5 * pa_y[j] * fl1_fx * pb_xyy[j] + pa_yzz[j] * pb_xyy[j]);

                t_yzz_xyy[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - fl1_fx * fl1_fz * fl1_fga * pb_xy[j] - pa_yzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 2.0 * pa_y[j] * fl1_fz * fl2_fx * pb_x[j] + 4.0 * fl2_fx * fl1_fz * pb_xy[j] + 5.0 * pa_yzz[j] * fl1_fz * pb_x[j] * fl1_fx + 10.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xy[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xyy[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_xyy[j]);

                t_yzz_xyz[j] = fl_s_0_0 * (0.5 * fl2_fx * pa_z[j] * pb_x[j] + 0.25 * fl2_fx * pb_xz[j] + pa_yz[j] * fl1_fx * pb_xy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xz[j] + 0.5 * pa_y[j] * fl1_fx * pb_xyz[j] + pa_yzz[j] * pb_xyz[j]);

                t_yzz_xyz[j] += fl_r_0_0 * (4.0 * fl2_fx * pa_z[j] * fl1_fz * pb_x[j] - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_xz[j] + 2.0 * fl2_fx * fl1_fz * pb_xz[j] + 10.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xy[j] + 5.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_xz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xyz[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_xyz[j]);

                t_yzz_xzz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl2_fx * pb_x[j] + 0.5 * pa_yzz[j] * pb_x[j] * fl1_fx + 2.0 * pa_yz[j] * fl1_fx * pb_xz[j] + 0.5 * pa_y[j] * fl1_fx * pb_xzz[j] + pa_yzz[j] * pb_xzz[j]);

                t_yzz_xzz[j] += fl_r_0_0 * (6.0 * pa_y[j] * fl2_fx * fl1_fz * pb_x[j] - 0.5 * pa_y[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - pa_yzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 5.0 * pa_yzz[j] * fl1_fz * pb_x[j] * fl1_fx + 20.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_xz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_xzz[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_xzz[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_xzz[j]);

                t_yzz_yyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * fl2_fx * pa_zz[j] + 0.75 * pa_y[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pb_yy[j] + 1.5 * pa_yzz[j] * pb_y[j] * fl1_fx + 1.5 * fl1_fx * pa_zz[j] * pb_yy[j] + 0.5 * pa_y[j] * fl1_fx * pb_yyy[j] + pa_yzz[j] * pb_yyy[j]);

                t_yzz_yyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * fl2_fx * pa_zz[j] * fl1_fz - 1.5 * pa_y[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yy[j] - 3.0 * pa_yzz[j] * pb_y[j] * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl1_fz * fl2_fx * pb_y[j] + 6.0 * fl2_fx * fl1_fz * pb_yy[j] + 15.0 * pa_yzz[j] * fl1_fz * pb_y[j] * fl1_fx + 15.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yy[j] - pa_y[j] * fl1_fz * fl1_fga * pb_yyy[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyy[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_yyy[j]);

                t_yzz_yyz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + fl2_fx * pa_z[j] * pb_y[j] + 0.25 * pa_y[j] * fl2_fx * pb_z[j] + 0.5 * fl2_fx * pb_yz[j] + 0.5 * pa_yzz[j] * fl1_fx * pb_z[j] + pa_yz[j] * fl1_fx * pb_yy[j] + fl1_fx * pa_zz[j] * pb_yz[j] + 0.5 * pa_y[j] * fl1_fx * pb_yyz[j] + pa_yzz[j] * pb_yyz[j]);

                t_yzz_yyz[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_yz[j] * fl1_fz * fl2_fx + 8.0 * fl2_fx * pa_z[j] * fl1_fz * pb_y[j] - 0.5 * pa_y[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - fl1_fx * fl1_fz * fl1_fga * pb_yz[j] - pa_yzz[j] * fl1_fz * fl1_fgb * pb_z[j] + 2.0 * pa_y[j] * fl1_fz * fl2_fx * pb_z[j] + 4.0 * fl2_fx * fl1_fz * pb_yz[j] + 5.0 * pa_yzz[j] * fl1_fz * fl1_fx * pb_z[j] + 10.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_yy[j] + 10.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_yz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_yyz[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yyz[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_yyz[j]);

                t_yzz_yzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_y[j] * fl2_fx * pb_y[j] + 0.25 * fl2_fx * pa_zz[j] + fl2_fx * pa_z[j] * pb_z[j] + 0.25 * fl2_fx * pb_zz[j] + 0.5 * pa_yzz[j] * pb_y[j] * fl1_fx + 2.0 * pa_yz[j] * fl1_fx * pb_yz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_zz[j] + 0.5 * pa_y[j] * fl1_fx * pb_yzz[j] + pa_yzz[j] * pb_yzz[j]);

                t_yzz_yzz[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * fl1_fx * pa_zz[j] * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl2_fx * fl1_fz * pb_y[j] + 2.0 * fl2_fx * pa_zz[j] * fl1_fz + 8.0 * fl2_fx * pa_z[j] * fl1_fz * pb_z[j] - 0.5 * pa_y[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 0.5 * fl1_fx * fl1_fz * fl1_fga * pb_zz[j] - pa_yzz[j] * pb_y[j] * fl1_fz * fl1_fgb + 2.0 * fl2_fx * fl1_fz * pb_zz[j] + 5.0 * pa_yzz[j] * fl1_fz * pb_y[j] * fl1_fx + 20.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_yz[j] + 5.0 * fl1_fx * pa_zz[j] * fl1_fz * pb_zz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_yzz[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_yzz[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_yzz[j]);

                t_yzz_zzz[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl2_fx + 2.25 * pa_y[j] * fl2_fx * pb_z[j] + 1.5 * pa_yzz[j] * pb_z[j] * fl1_fx + 3.0 * pa_yz[j] * fl1_fx * pb_zz[j] + 0.5 * pa_y[j] * fl1_fx * pb_zzz[j] + pa_yzz[j] * pb_zzz[j]);

                t_yzz_zzz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_yz[j] * fl1_fz * fl2_fx + 18.0 * pa_y[j] * fl2_fx * fl1_fz * pb_z[j] - 1.5 * pa_y[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * pb_z[j] * fl1_fx - 3.0 * pa_yzz[j] * pb_z[j] * fl1_fz * fl1_fgb + 15.0 * pa_yzz[j] * fl1_fz * pb_z[j] * fl1_fx + 30.0 * pa_yz[j] * fl1_fz * fl1_fx * pb_zz[j] - pa_y[j] * fl1_fz * fl1_fga * pb_zzz[j] + 5.0 * pa_y[j] * fl1_fz * fl1_fx * pb_zzz[j] + 12.0 * pa_yzz[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_90_100(      CMemBlock2D<double>& primBuffer,
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

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_zzz = paDistances.data(19 * idx + 18);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzz_xxx = primBuffer.data(100 * idx + 90);

            auto t_zzz_xxy = primBuffer.data(100 * idx + 91);

            auto t_zzz_xxz = primBuffer.data(100 * idx + 92);

            auto t_zzz_xyy = primBuffer.data(100 * idx + 93);

            auto t_zzz_xyz = primBuffer.data(100 * idx + 94);

            auto t_zzz_xzz = primBuffer.data(100 * idx + 95);

            auto t_zzz_yyy = primBuffer.data(100 * idx + 96);

            auto t_zzz_yyz = primBuffer.data(100 * idx + 97);

            auto t_zzz_yzz = primBuffer.data(100 * idx + 98);

            auto t_zzz_zzz = primBuffer.data(100 * idx + 99);

            // Batch of Integrals (90,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, r_0_0, s_0_0, t_zzz_xxx, t_zzz_xxy, t_zzz_xxz, t_zzz_xyy, t_zzz_xyz, \
                                     t_zzz_xzz, t_zzz_yyy, t_zzz_yyz, t_zzz_yzz, t_zzz_zzz: VLX_ALIGN)
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

                t_zzz_xxx[j] = fl_s_0_0 * (2.25 * pa_z[j] * fl2_fx * pb_x[j] + 1.5 * pa_zzz[j] * pb_x[j] * fl1_fx + 1.5 * pa_z[j] * fl1_fx * pb_xxx[j] + pa_zzz[j] * pb_xxx[j]);

                t_zzz_xxx[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 3.0 * pa_zzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 18.0 * pa_z[j] * fl1_fz * fl2_fx * pb_x[j] + 15.0 * pa_zzz[j] * fl1_fz * pb_x[j] * fl1_fx - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxx[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxx[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_xxx[j]);

                t_zzz_xxy[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx * pb_y[j] + 0.5 * pa_zzz[j] * fl1_fx * pb_y[j] + 1.5 * pa_z[j] * fl1_fx * pb_xxy[j] + pa_zzz[j] * pb_xxy[j]);

                t_zzz_xxy[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl1_fx * fl1_fz * fl1_fgb * pb_y[j] - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl1_fx * pb_y[j] - pa_zzz[j] * fl1_fz * fl1_fgb * pb_y[j] + 6.0 * pa_z[j] * fl1_fz * fl2_fx * pb_y[j] + 5.0 * pa_zzz[j] * fl1_fz * fl1_fx * pb_y[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxy[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxy[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_xxy[j]);

                t_zzz_xxz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.75 * pa_z[j] * fl2_fx * pb_z[j] + 0.75 * fl2_fx * pb_xx[j] + 0.5 * pa_zzz[j] * fl1_fx * pb_z[j] + 1.5 * pa_zz[j] * fl1_fx * pb_xx[j] + 1.5 * pa_z[j] * fl1_fx * pb_xxz[j] + pa_zzz[j] * pb_xxz[j]);

                t_zzz_xxz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_zz[j] * fl1_fz * fl2_fx - 1.5 * pa_z[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xx[j] - pa_zzz[j] * fl1_fz * fl1_fgb * pb_z[j] + 6.0 * pa_z[j] * fl1_fz * fl2_fx * pb_z[j] + 6.0 * fl2_fx * fl1_fz * pb_xx[j] + 5.0 * pa_zzz[j] * fl1_fz * fl1_fx * pb_z[j] + 15.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xx[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xxz[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xxz[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_xxz[j]);

                t_zzz_xyy[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl2_fx * pb_x[j] + 0.5 * pa_zzz[j] * pb_x[j] * fl1_fx + 1.5 * pa_z[j] * fl1_fx * pb_xyy[j] + pa_zzz[j] * pb_xyy[j]);

                t_zzz_xyy[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - pa_zzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 6.0 * pa_z[j] * fl1_fz * fl2_fx * pb_x[j] + 5.0 * pa_zzz[j] * fl1_fz * pb_x[j] * fl1_fx - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xyy[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xyy[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_xyy[j]);

                t_zzz_xyz[j] = fl_s_0_0 * (0.75 * fl2_fx * pb_xy[j] + 1.5 * pa_zz[j] * fl1_fx * pb_xy[j] + 1.5 * pa_z[j] * fl1_fx * pb_xyz[j] + pa_zzz[j] * pb_xyz[j]);

                t_zzz_xyz[j] += fl_r_0_0 * (-1.5 * fl1_fx * fl1_fz * fl1_fga * pb_xy[j] + 6.0 * fl2_fx * fl1_fz * pb_xy[j] + 15.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xy[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xyz[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xyz[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_xyz[j]);

                t_zzz_xzz[j] = fl_s_0_0 * (2.25 * pa_z[j] * fl2_fx * pb_x[j] + 1.5 * fl2_fx * pb_xz[j] + 0.5 * pa_zzz[j] * pb_x[j] * fl1_fx + 3.0 * pa_zz[j] * fl1_fx * pb_xz[j] + 1.5 * pa_z[j] * fl1_fx * pb_xzz[j] + pa_zzz[j] * pb_xzz[j]);

                t_zzz_xzz[j] += fl_r_0_0 * (18.0 * pa_z[j] * fl2_fx * fl1_fz * pb_x[j] - 1.5 * pa_z[j] * fl1_fx * pb_x[j] * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fz * fl1_fga * pb_x[j] * fl1_fx - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_xz[j] - pa_zzz[j] * pb_x[j] * fl1_fz * fl1_fgb + 12.0 * fl2_fx * fl1_fz * pb_xz[j] + 5.0 * pa_zzz[j] * fl1_fz * pb_x[j] * fl1_fx + 30.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_xz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_xzz[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_xzz[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_xzz[j]);

                t_zzz_yyy[j] = fl_s_0_0 * (2.25 * pa_z[j] * fl2_fx * pb_y[j] + 1.5 * pa_zzz[j] * pb_y[j] * fl1_fx + 1.5 * pa_z[j] * fl1_fx * pb_yyy[j] + pa_zzz[j] * pb_yyy[j]);

                t_zzz_yyy[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 3.0 * pa_zzz[j] * pb_y[j] * fl1_fz * fl1_fgb + 18.0 * pa_z[j] * fl1_fz * fl2_fx * pb_y[j] + 15.0 * pa_zzz[j] * fl1_fz * pb_y[j] * fl1_fx - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_yyy[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_yyy[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_yyy[j]);

                t_zzz_yyz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.75 * pa_z[j] * fl2_fx * pb_z[j] + 0.75 * fl2_fx * pb_yy[j] + 0.5 * pa_zzz[j] * fl1_fx * pb_z[j] + 1.5 * pa_zz[j] * fl1_fx * pb_yy[j] + 1.5 * pa_z[j] * fl1_fx * pb_yyz[j] + pa_zzz[j] * pb_yyz[j]);

                t_zzz_yyz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_zz[j] * fl1_fz * fl2_fx - 1.5 * pa_z[j] * fl1_fx * fl1_fz * fl1_fgb * pb_z[j] - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl1_fx * pb_z[j] - 1.5 * fl1_fx * fl1_fz * fl1_fga * pb_yy[j] - pa_zzz[j] * fl1_fz * fl1_fgb * pb_z[j] + 6.0 * pa_z[j] * fl1_fz * fl2_fx * pb_z[j] + 6.0 * fl2_fx * fl1_fz * pb_yy[j] + 5.0 * pa_zzz[j] * fl1_fz * fl1_fx * pb_z[j] + 15.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_yy[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_yyz[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_yyz[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_yyz[j]);

                t_zzz_yzz[j] = fl_s_0_0 * (2.25 * pa_z[j] * fl2_fx * pb_y[j] + 1.5 * fl2_fx * pb_yz[j] + 0.5 * pa_zzz[j] * pb_y[j] * fl1_fx + 3.0 * pa_zz[j] * fl1_fx * pb_yz[j] + 1.5 * pa_z[j] * fl1_fx * pb_yzz[j] + pa_zzz[j] * pb_yzz[j]);

                t_zzz_yzz[j] += fl_r_0_0 * (18.0 * pa_z[j] * fl2_fx * fl1_fz * pb_y[j] - 1.5 * pa_z[j] * fl1_fx * pb_y[j] * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fz * fl1_fga * pb_y[j] * fl1_fx - 3.0 * fl1_fx * fl1_fz * fl1_fga * pb_yz[j] - pa_zzz[j] * pb_y[j] * fl1_fz * fl1_fgb + 12.0 * fl2_fx * fl1_fz * pb_yz[j] + 5.0 * pa_zzz[j] * fl1_fz * pb_y[j] * fl1_fx + 30.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_yz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_yzz[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_yzz[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_yzz[j]);

                t_zzz_zzz[j] = fl_s_0_0 * (1.875 * fl3_fx + 2.25 * pa_zz[j] * fl2_fx + 6.75 * pa_z[j] * fl2_fx * pb_z[j] + 2.25 * fl2_fx * pb_zz[j] + 1.5 * pa_zzz[j] * pb_z[j] * fl1_fx + 4.5 * pa_zz[j] * fl1_fx * pb_zz[j] + 1.5 * pa_z[j] * fl1_fx * pb_zzz[j] + pa_zzz[j] * pb_zzz[j]);

                t_zzz_zzz[j] += fl_r_0_0 * (11.25 * fl3_fx * fl1_fz - 2.25 * fl2_fx * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa_zz[j] * fl1_fz * fl2_fx + 54.0 * pa_z[j] * fl2_fx * fl1_fz * pb_z[j] - 4.5 * pa_z[j] * fl1_fx * pb_z[j] * fl1_fz * fl1_fgb - 4.5 * pa_z[j] * fl1_fz * fl1_fga * pb_z[j] * fl1_fx - 4.5 * fl1_fx * fl1_fz * fl1_fga * pb_zz[j] - 3.0 * pa_zzz[j] * pb_z[j] * fl1_fz * fl1_fgb + 18.0 * fl2_fx * fl1_fz * pb_zz[j] + 15.0 * pa_zzz[j] * fl1_fz * pb_z[j] * fl1_fx + 45.0 * pa_zz[j] * fl1_fz * fl1_fx * pb_zz[j] - 3.0 * pa_z[j] * fl1_fz * fl1_fga * pb_zzz[j] + 15.0 * pa_z[j] * fl1_fz * fl1_fx * pb_zzz[j] + 12.0 * pa_zzz[j] * fl1_fz * pb_zzz[j]);
            }

            idx++;
        }
    }


} // kinrecfunc namespace

