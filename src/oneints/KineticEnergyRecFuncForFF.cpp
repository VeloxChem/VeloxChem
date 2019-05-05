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
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForFF_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_30_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_40_50(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_50_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_60_70(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_70_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_80_90(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFF_90_100(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                  braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compKineticEnergyForFF_0_10(      CMemBlock2D<double>& primBuffer,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xx = paDistances.data(19 * idx + 3);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(361 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(361 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(361 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(361 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(361 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(361 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(361 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(361 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(361 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(361 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(361 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(361 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(361 * idx + 18);

            auto pa2pb_xx_xx = pa2pbDistances.data(361 * idx + 60);

            auto pa2pb_xx_xy = pa2pbDistances.data(361 * idx + 61);

            auto pa2pb_xx_xz = pa2pbDistances.data(361 * idx + 62);

            auto pa2pb_xx_yy = pa2pbDistances.data(361 * idx + 63);

            auto pa2pb_xx_yz = pa2pbDistances.data(361 * idx + 64);

            auto pa2pb_xx_zz = pa2pbDistances.data(361 * idx + 65);

            auto pa2pb_xxx_x = pa2pbDistances.data(361 * idx + 171);

            auto pa2pb_xxx_y = pa2pbDistances.data(361 * idx + 172);

            auto pa2pb_xxx_z = pa2pbDistances.data(361 * idx + 173);

            auto pa2pb_xxx_xxx = pa2pbDistances.data(361 * idx + 180);

            auto pa2pb_xxx_xxy = pa2pbDistances.data(361 * idx + 181);

            auto pa2pb_xxx_xxz = pa2pbDistances.data(361 * idx + 182);

            auto pa2pb_xxx_xyy = pa2pbDistances.data(361 * idx + 183);

            auto pa2pb_xxx_xyz = pa2pbDistances.data(361 * idx + 184);

            auto pa2pb_xxx_xzz = pa2pbDistances.data(361 * idx + 185);

            auto pa2pb_xxx_yyy = pa2pbDistances.data(361 * idx + 186);

            auto pa2pb_xxx_yyz = pa2pbDistances.data(361 * idx + 187);

            auto pa2pb_xxx_yzz = pa2pbDistances.data(361 * idx + 188);

            auto pa2pb_xxx_zzz = pa2pbDistances.data(361 * idx + 189);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, \
                                     pa2pb_x_yzz, pa2pb_x_z, pa2pb_x_zzz, pa2pb_xx_xx, pa2pb_xx_xy, pa2pb_xx_xz, \
                                     pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa2pb_xxx_x, pa2pb_xxx_xxx, pa2pb_xxx_xxy, \
                                     pa2pb_xxx_xxz, pa2pb_xxx_xyy, pa2pb_xxx_xyz, pa2pb_xxx_xzz, pa2pb_xxx_y, \
                                     pa2pb_xxx_yyy, pa2pb_xxx_yyz, pa2pb_xxx_yzz, pa2pb_xxx_z, pa2pb_xxx_zzz, pa_xx, pb_xx, \
                                     pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, s_0_0, t_xxx_xxx, t_xxx_xxy, t_xxx_xxz, \
                                     t_xxx_xyy, t_xxx_xyz, t_xxx_xzz, t_xxx_yyy, t_xxx_yyz, t_xxx_yzz, t_xxx_zzz: VLX_ALIGN)
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

                t_xxx_xxx[j] = fl_s_0_0 * (1.875 * fl3_fx + 2.25 * pa_xx[j] * fl2_fx + 6.75 * pa2pb_x_x[j] * fl2_fx + 2.25 * pb_xx[j] * fl2_fx + 1.5 * pa2pb_xxx_x[j] * fl1_fx + 4.5 * pa2pb_xx_xx[j] * fl1_fx + 1.5 * pa2pb_x_xxx[j] * fl1_fx + pa2pb_xxx_xxx[j]);

                t_xxx_xxx[j] += fl_r_0_0 * (11.25 * fl3_fx * fl1_fz - 2.25 * fl2_fx * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa_xx[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxx_x[j] * fl1_fz * fl1_fgb + 18.0 * pb_xx[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xxx_x[j] * fl1_fz * fl1_fx + 45.0 * pa2pb_xx_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxx[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_xxx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_xxx[j] * fl1_fz);

                t_xxx_xxy[j] = fl_s_0_0 * (2.25 * pa2pb_x_y[j] * fl2_fx + 1.5 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xxx_y[j] * fl1_fx + 3.0 * pa2pb_xx_xy[j] * fl1_fx + 1.5 * pa2pb_x_xxy[j] * fl1_fx + pa2pb_xxx_xxy[j]);

                t_xxx_xxy[j] += fl_r_0_0 * (18.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_y[j] * fl1_fz * fl1_fgb + 12.0 * pb_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxx_y[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_xx_xy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_xxy[j] * fl1_fz);

                t_xxx_xxz[j] = fl_s_0_0 * (2.25 * pa2pb_x_z[j] * fl2_fx + 1.5 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xxx_z[j] * fl1_fx + 3.0 * pa2pb_xx_xz[j] * fl1_fx + 1.5 * pa2pb_x_xxz[j] * fl1_fx + pa2pb_xxx_xxz[j]);

                t_xxx_xxz[j] += fl_r_0_0 * (18.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_z[j] * fl1_fz * fl1_fgb + 12.0 * pb_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxx_z[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_xx_xz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_xxz[j] * fl1_fz);

                t_xxx_xyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_xxx_x[j] * fl1_fx + 1.5 * pa2pb_xx_yy[j] * fl1_fx + 1.5 * pa2pb_x_xyy[j] * fl1_fx + pa2pb_xxx_xyy[j]);

                t_xxx_xyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_xx[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_x[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_x[j] * fl1_fz * fl2_fx + 6.0 * pb_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxx_x[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_xx_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xyy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_xyy[j] * fl1_fz);

                t_xxx_xyz[j] = fl_s_0_0 * (0.75 * pb_yz[j] * fl2_fx + 1.5 * pa2pb_xx_yz[j] * fl1_fx + 1.5 * pa2pb_x_xyz[j] * fl1_fx + pa2pb_xxx_xyz[j]);

                t_xxx_xyz[j] += fl_r_0_0 * (-1.5 * pb_yz[j] * fl1_fx * fl1_fz * fl1_fga + 6.0 * pb_yz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xx_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xyz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_xyz[j] * fl1_fz);

                t_xxx_xzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_xxx_x[j] * fl1_fx + 1.5 * pa2pb_xx_zz[j] * fl1_fx + 1.5 * pa2pb_x_xzz[j] * fl1_fx + pa2pb_xxx_xzz[j]);

                t_xxx_xzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_xx[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_x[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_x[j] * fl1_fz * fl2_fx + 6.0 * pb_zz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxx_x[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_xx_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_xzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_xzz[j] * fl1_fz);

                t_xxx_yyy[j] = fl_s_0_0 * (2.25 * pa2pb_x_y[j] * fl2_fx + 1.5 * pa2pb_xxx_y[j] * fl1_fx + 1.5 * pa2pb_x_yyy[j] * fl1_fx + pa2pb_xxx_yyy[j]);

                t_xxx_yyy[j] += fl_r_0_0 * (-4.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxx_y[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_x_y[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xxx_y[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_yyy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_yyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_yyy[j] * fl1_fz);

                t_xxx_yyz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa2pb_xxx_z[j] * fl1_fx + 1.5 * pa2pb_x_yyz[j] * fl1_fx + pa2pb_xxx_yyz[j]);

                t_xxx_yyz[j] += fl_r_0_0 * (-1.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxx_z[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_z[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xxx_z[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_yyz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_yyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_yyz[j] * fl1_fz);

                t_xxx_yzz[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_xxx_y[j] * fl1_fx + 1.5 * pa2pb_x_yzz[j] * fl1_fx + pa2pb_xxx_yzz[j]);

                t_xxx_yzz[j] += fl_r_0_0 * (-1.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxx_y[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_y[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xxx_y[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_yzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_yzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_yzz[j] * fl1_fz);

                t_xxx_zzz[j] = fl_s_0_0 * (2.25 * pa2pb_x_z[j] * fl2_fx + 1.5 * pa2pb_xxx_z[j] * fl1_fx + 1.5 * pa2pb_x_zzz[j] * fl1_fx + pa2pb_xxx_zzz[j]);

                t_xxx_zzz[j] += fl_r_0_0 * (-4.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxx_z[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_x_z[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xxx_z[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_zzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_x_zzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_zzz[j] * fl1_fz);
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(361 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(361 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(361 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(361 * idx + 19);

            auto pa2pb_y_y = pa2pbDistances.data(361 * idx + 20);

            auto pa2pb_y_z = pa2pbDistances.data(361 * idx + 21);

            auto pa2pb_y_xxx = pa2pbDistances.data(361 * idx + 28);

            auto pa2pb_y_xxy = pa2pbDistances.data(361 * idx + 29);

            auto pa2pb_y_xxz = pa2pbDistances.data(361 * idx + 30);

            auto pa2pb_y_xyy = pa2pbDistances.data(361 * idx + 31);

            auto pa2pb_y_xyz = pa2pbDistances.data(361 * idx + 32);

            auto pa2pb_y_xzz = pa2pbDistances.data(361 * idx + 33);

            auto pa2pb_y_yyy = pa2pbDistances.data(361 * idx + 34);

            auto pa2pb_y_yyz = pa2pbDistances.data(361 * idx + 35);

            auto pa2pb_y_yzz = pa2pbDistances.data(361 * idx + 36);

            auto pa2pb_y_zzz = pa2pbDistances.data(361 * idx + 37);

            auto pa2pb_xx_xx = pa2pbDistances.data(361 * idx + 60);

            auto pa2pb_xx_xy = pa2pbDistances.data(361 * idx + 61);

            auto pa2pb_xx_xz = pa2pbDistances.data(361 * idx + 62);

            auto pa2pb_xx_yy = pa2pbDistances.data(361 * idx + 63);

            auto pa2pb_xx_yz = pa2pbDistances.data(361 * idx + 64);

            auto pa2pb_xx_zz = pa2pbDistances.data(361 * idx + 65);

            auto pa2pb_xy_xx = pa2pbDistances.data(361 * idx + 79);

            auto pa2pb_xy_xy = pa2pbDistances.data(361 * idx + 80);

            auto pa2pb_xy_xz = pa2pbDistances.data(361 * idx + 81);

            auto pa2pb_xy_yy = pa2pbDistances.data(361 * idx + 82);

            auto pa2pb_xy_yz = pa2pbDistances.data(361 * idx + 83);

            auto pa2pb_xy_zz = pa2pbDistances.data(361 * idx + 84);

            auto pa2pb_xxy_x = pa2pbDistances.data(361 * idx + 190);

            auto pa2pb_xxy_y = pa2pbDistances.data(361 * idx + 191);

            auto pa2pb_xxy_z = pa2pbDistances.data(361 * idx + 192);

            auto pa2pb_xxy_xxx = pa2pbDistances.data(361 * idx + 199);

            auto pa2pb_xxy_xxy = pa2pbDistances.data(361 * idx + 200);

            auto pa2pb_xxy_xxz = pa2pbDistances.data(361 * idx + 201);

            auto pa2pb_xxy_xyy = pa2pbDistances.data(361 * idx + 202);

            auto pa2pb_xxy_xyz = pa2pbDistances.data(361 * idx + 203);

            auto pa2pb_xxy_xzz = pa2pbDistances.data(361 * idx + 204);

            auto pa2pb_xxy_yyy = pa2pbDistances.data(361 * idx + 205);

            auto pa2pb_xxy_yyz = pa2pbDistances.data(361 * idx + 206);

            auto pa2pb_xxy_yzz = pa2pbDistances.data(361 * idx + 207);

            auto pa2pb_xxy_zzz = pa2pbDistances.data(361 * idx + 208);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xy, \
                                     pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa2pb_xxy_x, pa2pb_xxy_xxx, \
                                     pa2pb_xxy_xxy, pa2pb_xxy_xxz, pa2pb_xxy_xyy, pa2pb_xxy_xyz, pa2pb_xxy_xzz, \
                                     pa2pb_xxy_y, pa2pb_xxy_yyy, pa2pb_xxy_yyz, pa2pb_xxy_yzz, pa2pb_xxy_z, \
                                     pa2pb_xxy_zzz, pa2pb_xy_xx, pa2pb_xy_xy, pa2pb_xy_xz, pa2pb_xy_yy, pa2pb_xy_yz, \
                                     pa2pb_xy_zz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, pa2pb_y_xyy, \
                                     pa2pb_y_xyz, pa2pb_y_xzz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, \
                                     pa2pb_y_z, pa2pb_y_zzz, pa_xx, pa_xy, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, \
                                     s_0_0, t_xxy_xxx, t_xxy_xxy, t_xxy_xxz, t_xxy_xyy, t_xxy_xyz, t_xxy_xzz, \
                                     t_xxy_yyy, t_xxy_yyz, t_xxy_yzz, t_xxy_zzz: VLX_ALIGN)
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

                t_xxy_xxx[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl2_fx + 2.25 * pa2pb_y_x[j] * fl2_fx + 1.5 * pa2pb_xxy_x[j] * fl1_fx + 3.0 * pa2pb_xy_xx[j] * fl1_fx + 0.5 * pa2pb_y_xxx[j] * fl1_fx + pa2pb_xxy_xxx[j]);

                t_xxy_xxx[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_xy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xxy_x[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz - pa2pb_y_xxx[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xxx[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_xxx[j] * fl1_fz);

                t_xxy_xxy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_xxy_y[j] * fl1_fx + 0.5 * pa2pb_xx_xx[j] * fl1_fx + 2.0 * pa2pb_xy_xy[j] * fl1_fx + 0.5 * pa2pb_y_xxy[j] * fl1_fx + pa2pb_xxy_xxy[j]);

                t_xxy_xxy[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xx[j] * fl1_fz * fl2_fx + 8.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_y[j] * fl1_fz * fl1_fgb + 2.0 * pb_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxy_y[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xx_xx[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz - pa2pb_y_xxy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xxy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_xxy[j] * fl1_fz);

                t_xxy_xxz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa2pb_xxy_z[j] * fl1_fx + 2.0 * pa2pb_xy_xz[j] * fl1_fx + 0.5 * pa2pb_y_xxz[j] * fl1_fx + pa2pb_xxy_xxz[j]);

                t_xxy_xxz[j] += fl_r_0_0 * (6.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_z[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xxy_z[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz - pa2pb_y_xxz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xxz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_xxz[j] * fl1_fz);

                t_xxy_xyy[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + pa2pb_x_y[j] * fl2_fx + 0.25 * pa2pb_y_x[j] * fl2_fx + 0.5 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xxy_x[j] * fl1_fx + pa2pb_xx_xy[j] * fl1_fx + pa2pb_xy_yy[j] * fl1_fx + 0.5 * pa2pb_y_xyy[j] * fl1_fx + pa2pb_xxy_xyy[j]);

                t_xxy_xyy[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xy[j] * fl2_fx * fl1_fz + 8.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - pb_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_x[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz + 4.0 * pb_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxy_x[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xx_xy[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz - pa2pb_y_xyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xyy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_xyy[j] * fl1_fz);

                t_xxy_xyz[j] = fl_s_0_0 * (0.5 * pa2pb_x_z[j] * fl2_fx + 0.25 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xx_xz[j] * fl1_fx + pa2pb_xy_yz[j] * fl1_fx + 0.5 * pa2pb_y_xyz[j] * fl1_fx + pa2pb_xxy_xyz[j]);

                t_xxy_xyz[j] += fl_r_0_0 * (4.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz - 0.5 * pb_xz[j] * fl1_fz * fl1_fga * fl1_fx + 2.0 * pb_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xx_xz[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz - pa2pb_y_xyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xyz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_xyz[j] * fl1_fz);

                t_xxy_xzz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + 0.25 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa2pb_xxy_x[j] * fl1_fx + pa2pb_xy_zz[j] * fl1_fx + 0.5 * pa2pb_y_xzz[j] * fl1_fx + pa2pb_xxy_xzz[j]);

                t_xxy_xzz[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_x[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxy_x[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz - pa2pb_y_xzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_xzz[j] * fl1_fz);

                t_xxy_yyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pb_yy[j] * fl2_fx + 1.5 * pa2pb_xxy_y[j] * fl1_fx + 1.5 * pa2pb_xx_yy[j] * fl1_fx + 0.5 * pa2pb_y_yyy[j] * fl1_fx + pa2pb_xxy_yyy[j]);

                t_xxy_yyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_xx[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_y[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz + 6.0 * pb_yy[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xxy_y[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_xx_yy[j] * fl1_fz * fl1_fx - pa2pb_y_yyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_yyy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_yyy[j] * fl1_fz);

                t_xxy_yyz[j] = fl_s_0_0 * (0.25 * pa2pb_y_z[j] * fl2_fx + 0.5 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_xxy_z[j] * fl1_fx + pa2pb_xx_yz[j] * fl1_fx + 0.5 * pa2pb_y_yyz[j] * fl1_fx + pa2pb_xxy_yyz[j]);

                t_xxy_yyz[j] += fl_r_0_0 * (-0.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - pb_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_z[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz + 4.0 * pb_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxy_z[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xx_yz[j] * fl1_fz * fl1_fx - pa2pb_y_yyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_yyz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_yyz[j] * fl1_fz);

                t_xxy_yzz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + 0.25 * pa2pb_y_y[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_xxy_y[j] * fl1_fx + 0.5 * pa2pb_xx_zz[j] * fl1_fx + 0.5 * pa2pb_y_yzz[j] * fl1_fx + pa2pb_xxy_yzz[j]);

                t_xxy_yzz[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_xx[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_y[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz + 2.0 * pb_zz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxy_y[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xx_zz[j] * fl1_fz * fl1_fx - pa2pb_y_yzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_yzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_yzz[j] * fl1_fz);

                t_xxy_zzz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl2_fx + 1.5 * pa2pb_xxy_z[j] * fl1_fx + 0.5 * pa2pb_y_zzz[j] * fl1_fx + pa2pb_xxy_zzz[j]);

                t_xxy_zzz[j] += fl_r_0_0 * (-1.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_z[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xxy_z[j] * fl1_fz * fl1_fx - pa2pb_y_zzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_zzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxy_zzz[j] * fl1_fz);
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xz = paDistances.data(19 * idx + 5);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(361 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(361 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(361 * idx + 2);

            auto pa2pb_z_x = pa2pbDistances.data(361 * idx + 38);

            auto pa2pb_z_y = pa2pbDistances.data(361 * idx + 39);

            auto pa2pb_z_z = pa2pbDistances.data(361 * idx + 40);

            auto pa2pb_z_xxx = pa2pbDistances.data(361 * idx + 47);

            auto pa2pb_z_xxy = pa2pbDistances.data(361 * idx + 48);

            auto pa2pb_z_xxz = pa2pbDistances.data(361 * idx + 49);

            auto pa2pb_z_xyy = pa2pbDistances.data(361 * idx + 50);

            auto pa2pb_z_xyz = pa2pbDistances.data(361 * idx + 51);

            auto pa2pb_z_xzz = pa2pbDistances.data(361 * idx + 52);

            auto pa2pb_z_yyy = pa2pbDistances.data(361 * idx + 53);

            auto pa2pb_z_yyz = pa2pbDistances.data(361 * idx + 54);

            auto pa2pb_z_yzz = pa2pbDistances.data(361 * idx + 55);

            auto pa2pb_z_zzz = pa2pbDistances.data(361 * idx + 56);

            auto pa2pb_xx_xx = pa2pbDistances.data(361 * idx + 60);

            auto pa2pb_xx_xy = pa2pbDistances.data(361 * idx + 61);

            auto pa2pb_xx_xz = pa2pbDistances.data(361 * idx + 62);

            auto pa2pb_xx_yy = pa2pbDistances.data(361 * idx + 63);

            auto pa2pb_xx_yz = pa2pbDistances.data(361 * idx + 64);

            auto pa2pb_xx_zz = pa2pbDistances.data(361 * idx + 65);

            auto pa2pb_xz_xx = pa2pbDistances.data(361 * idx + 98);

            auto pa2pb_xz_xy = pa2pbDistances.data(361 * idx + 99);

            auto pa2pb_xz_xz = pa2pbDistances.data(361 * idx + 100);

            auto pa2pb_xz_yy = pa2pbDistances.data(361 * idx + 101);

            auto pa2pb_xz_yz = pa2pbDistances.data(361 * idx + 102);

            auto pa2pb_xz_zz = pa2pbDistances.data(361 * idx + 103);

            auto pa2pb_xxz_x = pa2pbDistances.data(361 * idx + 209);

            auto pa2pb_xxz_y = pa2pbDistances.data(361 * idx + 210);

            auto pa2pb_xxz_z = pa2pbDistances.data(361 * idx + 211);

            auto pa2pb_xxz_xxx = pa2pbDistances.data(361 * idx + 218);

            auto pa2pb_xxz_xxy = pa2pbDistances.data(361 * idx + 219);

            auto pa2pb_xxz_xxz = pa2pbDistances.data(361 * idx + 220);

            auto pa2pb_xxz_xyy = pa2pbDistances.data(361 * idx + 221);

            auto pa2pb_xxz_xyz = pa2pbDistances.data(361 * idx + 222);

            auto pa2pb_xxz_xzz = pa2pbDistances.data(361 * idx + 223);

            auto pa2pb_xxz_yyy = pa2pbDistances.data(361 * idx + 224);

            auto pa2pb_xxz_yyz = pa2pbDistances.data(361 * idx + 225);

            auto pa2pb_xxz_yzz = pa2pbDistances.data(361 * idx + 226);

            auto pa2pb_xxz_zzz = pa2pbDistances.data(361 * idx + 227);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xy, \
                                     pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa2pb_xxz_x, pa2pb_xxz_xxx, \
                                     pa2pb_xxz_xxy, pa2pb_xxz_xxz, pa2pb_xxz_xyy, pa2pb_xxz_xyz, pa2pb_xxz_xzz, \
                                     pa2pb_xxz_y, pa2pb_xxz_yyy, pa2pb_xxz_yyz, pa2pb_xxz_yzz, pa2pb_xxz_z, \
                                     pa2pb_xxz_zzz, pa2pb_xz_xx, pa2pb_xz_xy, pa2pb_xz_xz, pa2pb_xz_yy, pa2pb_xz_yz, \
                                     pa2pb_xz_zz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_xxz, pa2pb_z_xyy, \
                                     pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, \
                                     pa2pb_z_z, pa2pb_z_zzz, pa_xx, pa_xz, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, \
                                     s_0_0, t_xxz_xxx, t_xxz_xxy, t_xxz_xxz, t_xxz_xyy, t_xxz_xyz, t_xxz_xzz, \
                                     t_xxz_yyy, t_xxz_yyz, t_xxz_yzz, t_xxz_zzz: VLX_ALIGN)
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

                t_xxz_xxx[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl2_fx + 2.25 * pa2pb_z_x[j] * fl2_fx + 1.5 * pa2pb_xxz_x[j] * fl1_fx + 3.0 * pa2pb_xz_xx[j] * fl1_fx + 0.5 * pa2pb_z_xxx[j] * fl1_fx + pa2pb_xxz_xxx[j]);

                t_xxz_xxx[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_xz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xxz_x[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz - pa2pb_z_xxx[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xxx[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_xxx[j] * fl1_fz);

                t_xxz_xxy[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa2pb_xxz_y[j] * fl1_fx + 2.0 * pa2pb_xz_xy[j] * fl1_fx + 0.5 * pa2pb_z_xxy[j] * fl1_fx + pa2pb_xxz_xxy[j]);

                t_xxz_xxy[j] += fl_r_0_0 * (6.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_y[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xxz_y[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz - pa2pb_z_xxy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xxy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_xxy[j] * fl1_fz);

                t_xxz_xxz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + pa2pb_x_x[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_xxz_z[j] * fl1_fx + 0.5 * pa2pb_xx_xx[j] * fl1_fx + 2.0 * pa2pb_xz_xz[j] * fl1_fx + 0.5 * pa2pb_z_xxz[j] * fl1_fx + pa2pb_xxz_xxz[j]);

                t_xxz_xxz[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xx[j] * fl1_fz * fl2_fx + 8.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_z[j] * fl1_fz * fl1_fgb + 2.0 * pb_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxz_z[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xx_xx[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz - pa2pb_z_xxz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xxz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_xxz[j] * fl1_fz);

                t_xxz_xyy[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + 0.25 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa2pb_xxz_x[j] * fl1_fx + pa2pb_xz_yy[j] * fl1_fx + 0.5 * pa2pb_z_xyy[j] * fl1_fx + pa2pb_xxz_xyy[j]);

                t_xxz_xyy[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_x[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxz_x[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz - pa2pb_z_xyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xyy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_xyy[j] * fl1_fz);

                t_xxz_xyz[j] = fl_s_0_0 * (0.5 * pa2pb_x_y[j] * fl2_fx + 0.25 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xx_xy[j] * fl1_fx + pa2pb_xz_yz[j] * fl1_fx + 0.5 * pa2pb_z_xyz[j] * fl1_fx + pa2pb_xxz_xyz[j]);

                t_xxz_xyz[j] += fl_r_0_0 * (4.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz - 0.5 * pb_xy[j] * fl1_fz * fl1_fga * fl1_fx + 2.0 * pb_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xx_xy[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz - pa2pb_z_xyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xyz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_xyz[j] * fl1_fz);

                t_xxz_xzz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + pa2pb_x_z[j] * fl2_fx + 0.25 * pa2pb_z_x[j] * fl2_fx + 0.5 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xxz_x[j] * fl1_fx + pa2pb_xx_xz[j] * fl1_fx + pa2pb_xz_zz[j] * fl1_fx + 0.5 * pa2pb_z_xzz[j] * fl1_fx + pa2pb_xxz_xzz[j]);

                t_xxz_xzz[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xz[j] * fl2_fx * fl1_fz + 8.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - pb_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_x[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz + 4.0 * pb_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxz_x[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xx_xz[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz - pa2pb_z_xzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_xzz[j] * fl1_fz);

                t_xxz_yyy[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl2_fx + 1.5 * pa2pb_xxz_y[j] * fl1_fx + 0.5 * pa2pb_z_yyy[j] * fl1_fx + pa2pb_xxz_yyy[j]);

                t_xxz_yyy[j] += fl_r_0_0 * (-1.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_y[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xxz_y[j] * fl1_fz * fl1_fx - pa2pb_z_yyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_yyy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_yyy[j] * fl1_fz);

                t_xxz_yyz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_xx[j] * fl2_fx + 0.25 * pa2pb_z_z[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_xxz_z[j] * fl1_fx + 0.5 * pa2pb_xx_yy[j] * fl1_fx + 0.5 * pa2pb_z_yyz[j] * fl1_fx + pa2pb_xxz_yyz[j]);

                t_xxz_yyz[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_xx[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_z[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz + 2.0 * pb_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxz_z[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xx_yy[j] * fl1_fz * fl1_fx - pa2pb_z_yyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_yyz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_yyz[j] * fl1_fz);

                t_xxz_yzz[j] = fl_s_0_0 * (0.25 * pa2pb_z_y[j] * fl2_fx + 0.5 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_xxz_y[j] * fl1_fx + pa2pb_xx_yz[j] * fl1_fx + 0.5 * pa2pb_z_yzz[j] * fl1_fx + pa2pb_xxz_yzz[j]);

                t_xxz_yzz[j] += fl_r_0_0 * (-0.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - pb_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_y[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz + 4.0 * pb_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xxz_y[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xx_yz[j] * fl1_fz * fl1_fx - pa2pb_z_yzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_yzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_yzz[j] * fl1_fz);

                t_xxz_zzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_xx[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.75 * pb_zz[j] * fl2_fx + 1.5 * pa2pb_xxz_z[j] * fl1_fx + 1.5 * pa2pb_xx_zz[j] * fl1_fx + 0.5 * pa2pb_z_zzz[j] * fl1_fx + pa2pb_xxz_zzz[j]);

                t_xxz_zzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_xx[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_z[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz + 6.0 * pb_zz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xxz_z[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_xx_zz[j] * fl1_fz * fl1_fx - pa2pb_z_zzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_zzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xxz_zzz[j] * fl1_fz);
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_yy = paDistances.data(19 * idx + 6);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(361 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(361 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(361 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(361 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(361 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(361 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(361 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(361 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(361 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(361 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(361 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(361 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(361 * idx + 18);

            auto pa2pb_y_x = pa2pbDistances.data(361 * idx + 19);

            auto pa2pb_y_y = pa2pbDistances.data(361 * idx + 20);

            auto pa2pb_y_z = pa2pbDistances.data(361 * idx + 21);

            auto pa2pb_xy_xx = pa2pbDistances.data(361 * idx + 79);

            auto pa2pb_xy_xy = pa2pbDistances.data(361 * idx + 80);

            auto pa2pb_xy_xz = pa2pbDistances.data(361 * idx + 81);

            auto pa2pb_xy_yy = pa2pbDistances.data(361 * idx + 82);

            auto pa2pb_xy_yz = pa2pbDistances.data(361 * idx + 83);

            auto pa2pb_xy_zz = pa2pbDistances.data(361 * idx + 84);

            auto pa2pb_yy_xx = pa2pbDistances.data(361 * idx + 117);

            auto pa2pb_yy_xy = pa2pbDistances.data(361 * idx + 118);

            auto pa2pb_yy_xz = pa2pbDistances.data(361 * idx + 119);

            auto pa2pb_yy_yy = pa2pbDistances.data(361 * idx + 120);

            auto pa2pb_yy_yz = pa2pbDistances.data(361 * idx + 121);

            auto pa2pb_yy_zz = pa2pbDistances.data(361 * idx + 122);

            auto pa2pb_xyy_x = pa2pbDistances.data(361 * idx + 228);

            auto pa2pb_xyy_y = pa2pbDistances.data(361 * idx + 229);

            auto pa2pb_xyy_z = pa2pbDistances.data(361 * idx + 230);

            auto pa2pb_xyy_xxx = pa2pbDistances.data(361 * idx + 237);

            auto pa2pb_xyy_xxy = pa2pbDistances.data(361 * idx + 238);

            auto pa2pb_xyy_xxz = pa2pbDistances.data(361 * idx + 239);

            auto pa2pb_xyy_xyy = pa2pbDistances.data(361 * idx + 240);

            auto pa2pb_xyy_xyz = pa2pbDistances.data(361 * idx + 241);

            auto pa2pb_xyy_xzz = pa2pbDistances.data(361 * idx + 242);

            auto pa2pb_xyy_yyy = pa2pbDistances.data(361 * idx + 243);

            auto pa2pb_xyy_yyz = pa2pbDistances.data(361 * idx + 244);

            auto pa2pb_xyy_yzz = pa2pbDistances.data(361 * idx + 245);

            auto pa2pb_xyy_zzz = pa2pbDistances.data(361 * idx + 246);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, \
                                     pa2pb_x_yzz, pa2pb_x_z, pa2pb_x_zzz, pa2pb_xy_xx, pa2pb_xy_xy, pa2pb_xy_xz, \
                                     pa2pb_xy_yy, pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_xyy_x, pa2pb_xyy_xxx, pa2pb_xyy_xxy, \
                                     pa2pb_xyy_xxz, pa2pb_xyy_xyy, pa2pb_xyy_xyz, pa2pb_xyy_xzz, pa2pb_xyy_y, \
                                     pa2pb_xyy_yyy, pa2pb_xyy_yyz, pa2pb_xyy_yzz, pa2pb_xyy_z, pa2pb_xyy_zzz, pa2pb_y_x, \
                                     pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xy, pa2pb_yy_xz, pa2pb_yy_yy, \
                                     pa2pb_yy_yz, pa2pb_yy_zz, pa_xy, pa_yy, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, \
                                     s_0_0, t_xyy_xxx, t_xyy_xxy, t_xyy_xxz, t_xyy_xyy, t_xyy_xyz, t_xyy_xzz, \
                                     t_xyy_yyy, t_xyy_yyz, t_xyy_yzz, t_xyy_zzz: VLX_ALIGN)
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

                t_xyy_xxx[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pb_xx[j] * fl2_fx + 1.5 * pa2pb_xyy_x[j] * fl1_fx + 1.5 * pa2pb_yy_xx[j] * fl1_fx + 0.5 * pa2pb_x_xxx[j] * fl1_fx + pa2pb_xyy_xxx[j]);

                t_xyy_xxx[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyy_x[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_x[j] * fl1_fz * fl2_fx + 6.0 * pb_xx[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xyy_x[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz - pa2pb_x_xxx[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xxx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_xxx[j] * fl1_fz);

                t_xyy_xxy[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + pa2pb_y_x[j] * fl2_fx + 0.25 * pa2pb_x_y[j] * fl2_fx + 0.5 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xyy_y[j] * fl1_fx + pa2pb_xy_xx[j] * fl1_fx + pa2pb_yy_xy[j] * fl1_fx + 0.5 * pa2pb_x_xxy[j] * fl1_fx + pa2pb_xyy_xxy[j]);

                t_xyy_xxy[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xy[j] * fl1_fz * fl2_fx + 8.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - pb_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_y[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_x_y[j] * fl1_fz * fl2_fx + 4.0 * pb_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xyy_y[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xy_xx[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz - pa2pb_x_xxy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_xxy[j] * fl1_fz);

                t_xyy_xxz[j] = fl_s_0_0 * (0.25 * pa2pb_x_z[j] * fl2_fx + 0.5 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xyy_z[j] * fl1_fx + pa2pb_yy_xz[j] * fl1_fx + 0.5 * pa2pb_x_xxz[j] * fl1_fx + pa2pb_xyy_xxz[j]);

                t_xyy_xxz[j] += fl_r_0_0 * (-0.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - pb_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_z[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_x_z[j] * fl1_fz * fl2_fx + 4.0 * pb_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xyy_z[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz - pa2pb_x_xxz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_xxz[j] * fl1_fz);

                t_xyy_xyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.25 * pa_yy[j] * fl2_fx + pa2pb_y_y[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_xyy_x[j] * fl1_fx + 2.0 * pa2pb_xy_xy[j] * fl1_fx + 0.5 * pa2pb_yy_yy[j] * fl1_fx + 0.5 * pa2pb_x_xyy[j] * fl1_fx + pa2pb_xyy_xyy[j]);

                t_xyy_xyy[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz + 2.0 * pa_yy[j] * fl2_fx * fl1_fz + 8.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_x[j] * fl1_fz * fl1_fgb + 2.0 * pb_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xyy_x[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_xy_xy[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz - pa2pb_x_xyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_xyy[j] * fl1_fz);

                t_xyy_xyz[j] = fl_s_0_0 * (0.5 * pa2pb_y_z[j] * fl2_fx + 0.25 * pb_yz[j] * fl2_fx + pa2pb_xy_xz[j] * fl1_fx + 0.5 * pa2pb_yy_yz[j] * fl1_fx + 0.5 * pa2pb_x_xyz[j] * fl1_fx + pa2pb_xyy_xyz[j]);

                t_xyy_xyz[j] += fl_r_0_0 * (4.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz - 0.5 * pb_yz[j] * fl1_fx * fl1_fz * fl1_fga + 2.0 * pb_yz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_xy_xz[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz - pa2pb_x_xyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_xyz[j] * fl1_fz);

                t_xyy_xzz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + 0.25 * pa2pb_x_x[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_xyy_x[j] * fl1_fx + 0.5 * pa2pb_yy_zz[j] * fl1_fx + 0.5 * pa2pb_x_xzz[j] * fl1_fx + pa2pb_xyy_xzz[j]);

                t_xyy_xzz[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_yy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_x[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_x_x[j] * fl1_fz * fl2_fx + 2.0 * pb_zz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xyy_x[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz - pa2pb_x_xzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_xzz[j] * fl1_fz);

                t_xyy_yyy[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl2_fx + 2.25 * pa2pb_x_y[j] * fl2_fx + 1.5 * pa2pb_xyy_y[j] * fl1_fx + 3.0 * pa2pb_xy_yy[j] * fl1_fx + 0.5 * pa2pb_x_yyy[j] * fl1_fx + pa2pb_xyy_yyy[j]);

                t_xyy_yyy[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_xy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyy_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xyy_y[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_xy_yy[j] * fl1_fz * fl1_fx - pa2pb_x_yyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_yyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_yyy[j] * fl1_fz);

                t_xyy_yyz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa2pb_xyy_z[j] * fl1_fx + 2.0 * pa2pb_xy_yz[j] * fl1_fx + 0.5 * pa2pb_x_yyz[j] * fl1_fx + pa2pb_xyy_yyz[j]);

                t_xyy_yyz[j] += fl_r_0_0 * (6.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xyy_z[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xyy_z[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_xy_yz[j] * fl1_fz * fl1_fx - pa2pb_x_yyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_yyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_yyz[j] * fl1_fz);

                t_xyy_yzz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl2_fx + 0.25 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_xyy_y[j] * fl1_fx + pa2pb_xy_zz[j] * fl1_fx + 0.5 * pa2pb_x_yzz[j] * fl1_fx + pa2pb_xyy_yzz[j]);

                t_xyy_yzz[j] += fl_r_0_0 * (-pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xy[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xyy_y[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_x_y[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xyy_y[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xy_zz[j] * fl1_fz * fl1_fx - pa2pb_x_yzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_yzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_yzz[j] * fl1_fz);

                t_xyy_zzz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl2_fx + 1.5 * pa2pb_xyy_z[j] * fl1_fx + 0.5 * pa2pb_x_zzz[j] * fl1_fx + pa2pb_xyy_zzz[j]);

                t_xyy_zzz[j] += fl_r_0_0 * (-1.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyy_z[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_z[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xyy_z[j] * fl1_fz * fl1_fx - pa2pb_x_zzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_zzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFF_40_50(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yz = paDistances.data(19 * idx + 7);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(361 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(361 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(361 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(361 * idx + 19);

            auto pa2pb_y_y = pa2pbDistances.data(361 * idx + 20);

            auto pa2pb_y_z = pa2pbDistances.data(361 * idx + 21);

            auto pa2pb_z_x = pa2pbDistances.data(361 * idx + 38);

            auto pa2pb_z_y = pa2pbDistances.data(361 * idx + 39);

            auto pa2pb_z_z = pa2pbDistances.data(361 * idx + 40);

            auto pa2pb_xy_xx = pa2pbDistances.data(361 * idx + 79);

            auto pa2pb_xy_xy = pa2pbDistances.data(361 * idx + 80);

            auto pa2pb_xy_xz = pa2pbDistances.data(361 * idx + 81);

            auto pa2pb_xy_yy = pa2pbDistances.data(361 * idx + 82);

            auto pa2pb_xy_yz = pa2pbDistances.data(361 * idx + 83);

            auto pa2pb_xy_zz = pa2pbDistances.data(361 * idx + 84);

            auto pa2pb_xz_xx = pa2pbDistances.data(361 * idx + 98);

            auto pa2pb_xz_xy = pa2pbDistances.data(361 * idx + 99);

            auto pa2pb_xz_xz = pa2pbDistances.data(361 * idx + 100);

            auto pa2pb_xz_yy = pa2pbDistances.data(361 * idx + 101);

            auto pa2pb_xz_yz = pa2pbDistances.data(361 * idx + 102);

            auto pa2pb_xz_zz = pa2pbDistances.data(361 * idx + 103);

            auto pa2pb_yz_xx = pa2pbDistances.data(361 * idx + 136);

            auto pa2pb_yz_xy = pa2pbDistances.data(361 * idx + 137);

            auto pa2pb_yz_xz = pa2pbDistances.data(361 * idx + 138);

            auto pa2pb_yz_yy = pa2pbDistances.data(361 * idx + 139);

            auto pa2pb_yz_yz = pa2pbDistances.data(361 * idx + 140);

            auto pa2pb_yz_zz = pa2pbDistances.data(361 * idx + 141);

            auto pa2pb_xyz_x = pa2pbDistances.data(361 * idx + 247);

            auto pa2pb_xyz_y = pa2pbDistances.data(361 * idx + 248);

            auto pa2pb_xyz_z = pa2pbDistances.data(361 * idx + 249);

            auto pa2pb_xyz_xxx = pa2pbDistances.data(361 * idx + 256);

            auto pa2pb_xyz_xxy = pa2pbDistances.data(361 * idx + 257);

            auto pa2pb_xyz_xxz = pa2pbDistances.data(361 * idx + 258);

            auto pa2pb_xyz_xyy = pa2pbDistances.data(361 * idx + 259);

            auto pa2pb_xyz_xyz = pa2pbDistances.data(361 * idx + 260);

            auto pa2pb_xyz_xzz = pa2pbDistances.data(361 * idx + 261);

            auto pa2pb_xyz_yyy = pa2pbDistances.data(361 * idx + 262);

            auto pa2pb_xyz_yyz = pa2pbDistances.data(361 * idx + 263);

            auto pa2pb_xyz_yzz = pa2pbDistances.data(361 * idx + 264);

            auto pa2pb_xyz_zzz = pa2pbDistances.data(361 * idx + 265);

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

            #pragma omp simd aligned(fgb, fx, fz, pa2pb_x_x, pa2pb_x_y, pa2pb_x_z, pa2pb_xy_xx, pa2pb_xy_xy, \
                                     pa2pb_xy_xz, pa2pb_xy_yy, pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_xyz_x, pa2pb_xyz_xxx, \
                                     pa2pb_xyz_xxy, pa2pb_xyz_xxz, pa2pb_xyz_xyy, pa2pb_xyz_xyz, pa2pb_xyz_xzz, \
                                     pa2pb_xyz_y, pa2pb_xyz_yyy, pa2pb_xyz_yyz, pa2pb_xyz_yzz, pa2pb_xyz_z, \
                                     pa2pb_xyz_zzz, pa2pb_xz_xx, pa2pb_xz_xy, pa2pb_xz_xz, pa2pb_xz_yy, pa2pb_xz_yz, \
                                     pa2pb_xz_zz, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa2pb_yz_xx, pa2pb_yz_xy, pa2pb_yz_xz, \
                                     pa2pb_yz_yy, pa2pb_yz_yz, pa2pb_yz_zz, pa2pb_z_x, pa2pb_z_y, pa2pb_z_z, pa_xy, pa_xz, \
                                     pa_yz, r_0_0, s_0_0, t_xyz_xxx, t_xyz_xxy, t_xyz_xxz, t_xyz_xyy, t_xyz_xyz, \
                                     t_xyz_xzz, t_xyz_yyy, t_xyz_yyz, t_xyz_yzz, t_xyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyz_xxx[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl2_fx + 1.5 * pa2pb_xyz_x[j] * fl1_fx + 1.5 * pa2pb_yz_xx[j] * fl1_fx + pa2pb_xyz_xxx[j]);

                t_xyz_xxx[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_yz[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xyz_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xyz_x[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_xxx[j] * fl1_fz);

                t_xyz_xxy[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa2pb_xyz_y[j] * fl1_fx + 0.5 * pa2pb_xz_xx[j] * fl1_fx + pa2pb_yz_xy[j] * fl1_fx + pa2pb_xyz_xxy[j]);

                t_xyz_xxy[j] += fl_r_0_0 * (-0.5 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xz[j] * fl2_fx * fl1_fz + 4.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz - pa2pb_xyz_y[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xyz_y[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz + 10.0 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_xxy[j] * fl1_fz);

                t_xyz_xxz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa2pb_xyz_z[j] * fl1_fx + 0.5 * pa2pb_xy_xx[j] * fl1_fx + pa2pb_yz_xz[j] * fl1_fx + pa2pb_xyz_xxz[j]);

                t_xyz_xxz[j] += fl_r_0_0 * (-0.5 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xy[j] * fl1_fz * fl2_fx + 4.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz - pa2pb_xyz_z[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xyz_z[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xy_xx[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_xxz[j] * fl1_fz);

                t_xyz_xyy[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl2_fx + 0.5 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa2pb_xyz_x[j] * fl1_fx + pa2pb_xz_xy[j] * fl1_fx + 0.5 * pa2pb_yz_yy[j] * fl1_fx + pa2pb_xyz_xyy[j]);

                t_xyz_xyy[j] += fl_r_0_0 * (-0.5 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_yz[j] * fl2_fx * fl1_fz + 4.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz - pa2pb_xyz_x[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xyz_x[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz + 5.0 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_xyy[j] * fl1_fz);

                t_xyz_xyz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa2pb_x_x[j] * fl2_fx + 0.25 * pa2pb_y_y[j] * fl2_fx + 0.25 * pa2pb_z_z[j] * fl2_fx + 0.5 * pa2pb_xy_xy[j] * fl1_fx + 0.5 * pa2pb_xz_xz[j] * fl1_fx + 0.5 * pa2pb_yz_yz[j] * fl1_fx + pa2pb_xyz_xyz[j]);

                t_xyz_xyz[j] += fl_r_0_0 * (0.75 * fl3_fx * fl1_fz + 2.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz + 2.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz + 2.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xy_xy[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz + 5.0 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_xyz[j] * fl1_fz);

                t_xyz_xzz[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl2_fx + 0.5 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa2pb_xyz_x[j] * fl1_fx + pa2pb_xy_xz[j] * fl1_fx + 0.5 * pa2pb_yz_zz[j] * fl1_fx + pa2pb_xyz_xzz[j]);

                t_xyz_xzz[j] += fl_r_0_0 * (-0.5 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_yz[j] * fl2_fx * fl1_fz + 4.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz - pa2pb_xyz_x[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xyz_x[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xy_xz[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_xzz[j] * fl1_fz);

                t_xyz_yyy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl2_fx + 1.5 * pa2pb_xyz_y[j] * fl1_fx + 1.5 * pa2pb_xz_yy[j] * fl1_fx + pa2pb_xyz_yyy[j]);

                t_xyz_yyy[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_xz[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xyz_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xyz_y[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_yyy[j] * fl1_fz);

                t_xyz_yyz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl2_fx + 0.5 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_xyz_z[j] * fl1_fx + 0.5 * pa2pb_xy_yy[j] * fl1_fx + pa2pb_xz_yz[j] * fl1_fx + pa2pb_xyz_yyz[j]);

                t_xyz_yyz[j] += fl_r_0_0 * (-0.5 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xy[j] * fl1_fz * fl2_fx + 4.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz - pa2pb_xyz_z[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xyz_z[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xy_yy[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_yyz[j] * fl1_fz);

                t_xyz_yzz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl2_fx + 0.5 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa2pb_xyz_y[j] * fl1_fx + pa2pb_xy_yz[j] * fl1_fx + 0.5 * pa2pb_xz_zz[j] * fl1_fx + pa2pb_xyz_yzz[j]);

                t_xyz_yzz[j] += fl_r_0_0 * (-0.5 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_xz[j] * fl2_fx * fl1_fz + 4.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz - pa2pb_xyz_y[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xyz_y[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xy_yz[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_yzz[j] * fl1_fz);

                t_xyz_zzz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl2_fx + 1.5 * pa2pb_xyz_z[j] * fl1_fx + 1.5 * pa2pb_xy_zz[j] * fl1_fx + pa2pb_xyz_zzz[j]);

                t_xyz_zzz[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_xy[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_xyz_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xyz_z[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_xy_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyz_zzz[j] * fl1_fz);
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(361 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(361 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(361 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(361 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(361 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(361 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(361 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(361 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(361 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(361 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(361 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(361 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(361 * idx + 18);

            auto pa2pb_z_x = pa2pbDistances.data(361 * idx + 38);

            auto pa2pb_z_y = pa2pbDistances.data(361 * idx + 39);

            auto pa2pb_z_z = pa2pbDistances.data(361 * idx + 40);

            auto pa2pb_xz_xx = pa2pbDistances.data(361 * idx + 98);

            auto pa2pb_xz_xy = pa2pbDistances.data(361 * idx + 99);

            auto pa2pb_xz_xz = pa2pbDistances.data(361 * idx + 100);

            auto pa2pb_xz_yy = pa2pbDistances.data(361 * idx + 101);

            auto pa2pb_xz_yz = pa2pbDistances.data(361 * idx + 102);

            auto pa2pb_xz_zz = pa2pbDistances.data(361 * idx + 103);

            auto pa2pb_zz_xx = pa2pbDistances.data(361 * idx + 155);

            auto pa2pb_zz_xy = pa2pbDistances.data(361 * idx + 156);

            auto pa2pb_zz_xz = pa2pbDistances.data(361 * idx + 157);

            auto pa2pb_zz_yy = pa2pbDistances.data(361 * idx + 158);

            auto pa2pb_zz_yz = pa2pbDistances.data(361 * idx + 159);

            auto pa2pb_zz_zz = pa2pbDistances.data(361 * idx + 160);

            auto pa2pb_xzz_x = pa2pbDistances.data(361 * idx + 266);

            auto pa2pb_xzz_y = pa2pbDistances.data(361 * idx + 267);

            auto pa2pb_xzz_z = pa2pbDistances.data(361 * idx + 268);

            auto pa2pb_xzz_xxx = pa2pbDistances.data(361 * idx + 275);

            auto pa2pb_xzz_xxy = pa2pbDistances.data(361 * idx + 276);

            auto pa2pb_xzz_xxz = pa2pbDistances.data(361 * idx + 277);

            auto pa2pb_xzz_xyy = pa2pbDistances.data(361 * idx + 278);

            auto pa2pb_xzz_xyz = pa2pbDistances.data(361 * idx + 279);

            auto pa2pb_xzz_xzz = pa2pbDistances.data(361 * idx + 280);

            auto pa2pb_xzz_yyy = pa2pbDistances.data(361 * idx + 281);

            auto pa2pb_xzz_yyz = pa2pbDistances.data(361 * idx + 282);

            auto pa2pb_xzz_yzz = pa2pbDistances.data(361 * idx + 283);

            auto pa2pb_xzz_zzz = pa2pbDistances.data(361 * idx + 284);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, \
                                     pa2pb_x_yzz, pa2pb_x_z, pa2pb_x_zzz, pa2pb_xz_xx, pa2pb_xz_xy, pa2pb_xz_xz, \
                                     pa2pb_xz_yy, pa2pb_xz_yz, pa2pb_xz_zz, pa2pb_xzz_x, pa2pb_xzz_xxx, pa2pb_xzz_xxy, \
                                     pa2pb_xzz_xxz, pa2pb_xzz_xyy, pa2pb_xzz_xyz, pa2pb_xzz_xzz, pa2pb_xzz_y, \
                                     pa2pb_xzz_yyy, pa2pb_xzz_yyz, pa2pb_xzz_yzz, pa2pb_xzz_z, pa2pb_xzz_zzz, pa2pb_z_x, \
                                     pa2pb_z_y, pa2pb_z_z, pa2pb_zz_xx, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, \
                                     pa2pb_zz_yz, pa2pb_zz_zz, pa_xz, pa_zz, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, \
                                     s_0_0, t_xzz_xxx, t_xzz_xxy, t_xzz_xxz, t_xzz_xyy, t_xzz_xyz, t_xzz_xzz, \
                                     t_xzz_yyy, t_xzz_yyz, t_xzz_yzz, t_xzz_zzz: VLX_ALIGN)
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

                t_xzz_xxx[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.75 * pb_xx[j] * fl2_fx + 1.5 * pa2pb_xzz_x[j] * fl1_fx + 1.5 * pa2pb_zz_xx[j] * fl1_fx + 0.5 * pa2pb_x_xxx[j] * fl1_fx + pa2pb_xzz_xxx[j]);

                t_xzz_xxx[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzz_x[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_x[j] * fl1_fz * fl2_fx + 6.0 * pb_xx[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xzz_x[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz - pa2pb_x_xxx[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xxx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_xxx[j] * fl1_fz);

                t_xzz_xxy[j] = fl_s_0_0 * (0.25 * pa2pb_x_y[j] * fl2_fx + 0.5 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_xzz_y[j] * fl1_fx + pa2pb_zz_xy[j] * fl1_fx + 0.5 * pa2pb_x_xxy[j] * fl1_fx + pa2pb_xzz_xxy[j]);

                t_xzz_xxy[j] += fl_r_0_0 * (-0.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - pb_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_y[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_x_y[j] * fl1_fz * fl2_fx + 4.0 * pb_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xzz_y[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz - pa2pb_x_xxy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_xxy[j] * fl1_fz);

                t_xzz_xxz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + pa2pb_z_x[j] * fl2_fx + 0.25 * pa2pb_x_z[j] * fl2_fx + 0.5 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_xzz_z[j] * fl1_fx + pa2pb_xz_xx[j] * fl1_fx + pa2pb_zz_xz[j] * fl1_fx + 0.5 * pa2pb_x_xxz[j] * fl1_fx + pa2pb_xzz_xxz[j]);

                t_xzz_xxz[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xz[j] * fl1_fz * fl2_fx + 8.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - pb_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_z[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_x_z[j] * fl1_fz * fl2_fx + 4.0 * pb_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xzz_z[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xz_xx[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz - pa2pb_x_xxz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_xxz[j] * fl1_fz);

                t_xzz_xyy[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_zz[j] * fl2_fx + 0.25 * pa2pb_x_x[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_xzz_x[j] * fl1_fx + 0.5 * pa2pb_zz_yy[j] * fl1_fx + 0.5 * pa2pb_x_xyy[j] * fl1_fx + pa2pb_xzz_xyy[j]);

                t_xzz_xyy[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_x[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_x_x[j] * fl1_fz * fl2_fx + 2.0 * pb_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xzz_x[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz - pa2pb_x_xyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_xyy[j] * fl1_fz);

                t_xzz_xyz[j] = fl_s_0_0 * (0.5 * pa2pb_z_y[j] * fl2_fx + 0.25 * pb_yz[j] * fl2_fx + pa2pb_xz_xy[j] * fl1_fx + 0.5 * pa2pb_zz_yz[j] * fl1_fx + 0.5 * pa2pb_x_xyz[j] * fl1_fx + pa2pb_xzz_xyz[j]);

                t_xzz_xyz[j] += fl_r_0_0 * (4.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz - 0.5 * pb_yz[j] * fl1_fx * fl1_fz * fl1_fga + 2.0 * pb_yz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_xz_xy[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz - pa2pb_x_xyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_xyz[j] * fl1_fz);

                t_xzz_xzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_x_x[j] * fl2_fx + 0.25 * pa_zz[j] * fl2_fx + pa2pb_z_z[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_xzz_x[j] * fl1_fx + 2.0 * pa2pb_xz_xz[j] * fl1_fx + 0.5 * pa2pb_zz_zz[j] * fl1_fx + 0.5 * pa2pb_x_xzz[j] * fl1_fx + pa2pb_xzz_xzz[j]);

                t_xzz_xzz[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz + 2.0 * pa_zz[j] * fl2_fx * fl1_fz + 8.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_x[j] * fl1_fz * fl1_fgb + 2.0 * pb_zz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xzz_x[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_xz_xz[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz - pa2pb_x_xzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_xzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_xzz[j] * fl1_fz);

                t_xzz_yyy[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl2_fx + 1.5 * pa2pb_xzz_y[j] * fl1_fx + 0.5 * pa2pb_x_yyy[j] * fl1_fx + pa2pb_xzz_yyy[j]);

                t_xzz_yyy[j] += fl_r_0_0 * (-1.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xzz_y[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_y[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xzz_y[j] * fl1_fz * fl1_fx - pa2pb_x_yyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_yyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_yyy[j] * fl1_fz);

                t_xzz_yyz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl2_fx + 0.25 * pa2pb_x_z[j] * fl2_fx + 0.5 * pa2pb_xzz_z[j] * fl1_fx + pa2pb_xz_yy[j] * fl1_fx + 0.5 * pa2pb_x_yyz[j] * fl1_fx + pa2pb_xzz_yyz[j]);

                t_xzz_yyz[j] += fl_r_0_0 * (-pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_xz[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xzz_z[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_x_z[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xzz_z[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_xz_yy[j] * fl1_fz * fl1_fx - pa2pb_x_yyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_yyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_yyz[j] * fl1_fz);

                t_xzz_yzz[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl2_fx + 0.5 * pa2pb_xzz_y[j] * fl1_fx + 2.0 * pa2pb_xz_yz[j] * fl1_fx + 0.5 * pa2pb_x_yzz[j] * fl1_fx + pa2pb_xzz_yzz[j]);

                t_xzz_yzz[j] += fl_r_0_0 * (6.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xzz_y[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_xzz_y[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_xz_yz[j] * fl1_fz * fl1_fx - pa2pb_x_yzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_yzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_yzz[j] * fl1_fz);

                t_xzz_zzz[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl2_fx + 2.25 * pa2pb_x_z[j] * fl2_fx + 1.5 * pa2pb_xzz_z[j] * fl1_fx + 3.0 * pa2pb_xz_zz[j] * fl1_fx + 0.5 * pa2pb_x_zzz[j] * fl1_fx + pa2pb_xzz_zzz[j]);

                t_xzz_zzz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_xz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xzz_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xzz_z[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_xz_zz[j] * fl1_fz * fl1_fx - pa2pb_x_zzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_x_zzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_zzz[j] * fl1_fz);
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_yy = paDistances.data(19 * idx + 6);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(361 * idx + 19);

            auto pa2pb_y_y = pa2pbDistances.data(361 * idx + 20);

            auto pa2pb_y_z = pa2pbDistances.data(361 * idx + 21);

            auto pa2pb_y_xxx = pa2pbDistances.data(361 * idx + 28);

            auto pa2pb_y_xxy = pa2pbDistances.data(361 * idx + 29);

            auto pa2pb_y_xxz = pa2pbDistances.data(361 * idx + 30);

            auto pa2pb_y_xyy = pa2pbDistances.data(361 * idx + 31);

            auto pa2pb_y_xyz = pa2pbDistances.data(361 * idx + 32);

            auto pa2pb_y_xzz = pa2pbDistances.data(361 * idx + 33);

            auto pa2pb_y_yyy = pa2pbDistances.data(361 * idx + 34);

            auto pa2pb_y_yyz = pa2pbDistances.data(361 * idx + 35);

            auto pa2pb_y_yzz = pa2pbDistances.data(361 * idx + 36);

            auto pa2pb_y_zzz = pa2pbDistances.data(361 * idx + 37);

            auto pa2pb_yy_xx = pa2pbDistances.data(361 * idx + 117);

            auto pa2pb_yy_xy = pa2pbDistances.data(361 * idx + 118);

            auto pa2pb_yy_xz = pa2pbDistances.data(361 * idx + 119);

            auto pa2pb_yy_yy = pa2pbDistances.data(361 * idx + 120);

            auto pa2pb_yy_yz = pa2pbDistances.data(361 * idx + 121);

            auto pa2pb_yy_zz = pa2pbDistances.data(361 * idx + 122);

            auto pa2pb_yyy_x = pa2pbDistances.data(361 * idx + 285);

            auto pa2pb_yyy_y = pa2pbDistances.data(361 * idx + 286);

            auto pa2pb_yyy_z = pa2pbDistances.data(361 * idx + 287);

            auto pa2pb_yyy_xxx = pa2pbDistances.data(361 * idx + 294);

            auto pa2pb_yyy_xxy = pa2pbDistances.data(361 * idx + 295);

            auto pa2pb_yyy_xxz = pa2pbDistances.data(361 * idx + 296);

            auto pa2pb_yyy_xyy = pa2pbDistances.data(361 * idx + 297);

            auto pa2pb_yyy_xyz = pa2pbDistances.data(361 * idx + 298);

            auto pa2pb_yyy_xzz = pa2pbDistances.data(361 * idx + 299);

            auto pa2pb_yyy_yyy = pa2pbDistances.data(361 * idx + 300);

            auto pa2pb_yyy_yyz = pa2pbDistances.data(361 * idx + 301);

            auto pa2pb_yyy_yzz = pa2pbDistances.data(361 * idx + 302);

            auto pa2pb_yyy_zzz = pa2pbDistances.data(361 * idx + 303);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, \
                                     pa2pb_y_xyy, pa2pb_y_xyz, pa2pb_y_xzz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, \
                                     pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, pa2pb_yy_xx, pa2pb_yy_xy, pa2pb_yy_xz, \
                                     pa2pb_yy_yy, pa2pb_yy_yz, pa2pb_yy_zz, pa2pb_yyy_x, pa2pb_yyy_xxx, pa2pb_yyy_xxy, \
                                     pa2pb_yyy_xxz, pa2pb_yyy_xyy, pa2pb_yyy_xyz, pa2pb_yyy_xzz, pa2pb_yyy_y, \
                                     pa2pb_yyy_yyy, pa2pb_yyy_yyz, pa2pb_yyy_yzz, pa2pb_yyy_z, pa2pb_yyy_zzz, pa_yy, pb_xx, \
                                     pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, s_0_0, t_yyy_xxx, t_yyy_xxy, t_yyy_xxz, \
                                     t_yyy_xyy, t_yyy_xyz, t_yyy_xzz, t_yyy_yyy, t_yyy_yyz, t_yyy_yzz, t_yyy_zzz: VLX_ALIGN)
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

                t_yyy_xxx[j] = fl_s_0_0 * (2.25 * pa2pb_y_x[j] * fl2_fx + 1.5 * pa2pb_yyy_x[j] * fl1_fx + 1.5 * pa2pb_y_xxx[j] * fl1_fx + pa2pb_yyy_xxx[j]);

                t_yyy_xxx[j] += fl_r_0_0 * (-4.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyy_x[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_y_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yyy_x[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxx[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_xxx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xxx[j] * fl1_fz);

                t_yyy_xxy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_yyy_y[j] * fl1_fx + 1.5 * pa2pb_yy_xx[j] * fl1_fx + 1.5 * pa2pb_y_xxy[j] * fl1_fx + pa2pb_yyy_xxy[j]);

                t_yyy_xxy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_yy[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_y[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_y[j] * fl1_fz * fl2_fx + 6.0 * pb_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyy_y[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_yy_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xxy[j] * fl1_fz);

                t_yyy_xxz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa2pb_yyy_z[j] * fl1_fx + 1.5 * pa2pb_y_xxz[j] * fl1_fx + pa2pb_yyy_xxz[j]);

                t_yyy_xxz[j] += fl_r_0_0 * (-1.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyy_z[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_z[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yyy_z[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xxz[j] * fl1_fz);

                t_yyy_xyy[j] = fl_s_0_0 * (2.25 * pa2pb_y_x[j] * fl2_fx + 1.5 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_yyy_x[j] * fl1_fx + 3.0 * pa2pb_yy_xy[j] * fl1_fx + 1.5 * pa2pb_y_xyy[j] * fl1_fx + pa2pb_yyy_xyy[j]);

                t_yyy_xyy[j] += fl_r_0_0 * (18.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_x[j] * fl1_fz * fl1_fgb + 12.0 * pb_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyy_x[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_yy_xy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xyy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xyy[j] * fl1_fz);

                t_yyy_xyz[j] = fl_s_0_0 * (0.75 * pb_xz[j] * fl2_fx + 1.5 * pa2pb_yy_xz[j] * fl1_fx + 1.5 * pa2pb_y_xyz[j] * fl1_fx + pa2pb_yyy_xyz[j]);

                t_yyy_xyz[j] += fl_r_0_0 * (-1.5 * pb_xz[j] * fl1_fx * fl1_fz * fl1_fga + 6.0 * pb_xz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_yy_xz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xyz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xyz[j] * fl1_fz);

                t_yyy_xzz[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa2pb_yyy_x[j] * fl1_fx + 1.5 * pa2pb_y_xzz[j] * fl1_fx + pa2pb_yyy_xzz[j]);

                t_yyy_xzz[j] += fl_r_0_0 * (-1.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyy_x[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_x[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yyy_x[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_xzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xzz[j] * fl1_fz);

                t_yyy_yyy[j] = fl_s_0_0 * (1.875 * fl3_fx + 2.25 * pa_yy[j] * fl2_fx + 6.75 * pa2pb_y_y[j] * fl2_fx + 2.25 * pb_yy[j] * fl2_fx + 1.5 * pa2pb_yyy_y[j] * fl1_fx + 4.5 * pa2pb_yy_yy[j] * fl1_fx + 1.5 * pa2pb_y_yyy[j] * fl1_fx + pa2pb_yyy_yyy[j]);

                t_yyy_yyy[j] += fl_r_0_0 * (11.25 * fl3_fx * fl1_fz - 2.25 * fl2_fx * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa_yy[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_yy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyy_y[j] * fl1_fz * fl1_fgb + 18.0 * pb_yy[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_yyy_y[j] * fl1_fz * fl1_fx + 45.0 * pa2pb_yy_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_yyy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_yyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_yyy[j] * fl1_fz);

                t_yyy_yyz[j] = fl_s_0_0 * (2.25 * pa2pb_y_z[j] * fl2_fx + 1.5 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_yyy_z[j] * fl1_fx + 3.0 * pa2pb_yy_yz[j] * fl1_fx + 1.5 * pa2pb_y_yyz[j] * fl1_fx + pa2pb_yyy_yyz[j]);

                t_yyy_yyz[j] += fl_r_0_0 * (18.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_z[j] * fl1_fz * fl1_fgb + 12.0 * pb_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyy_z[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_yy_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_yyz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_yyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_yyz[j] * fl1_fz);

                t_yyy_yzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_yyy_y[j] * fl1_fx + 1.5 * pa2pb_yy_zz[j] * fl1_fx + 1.5 * pa2pb_y_yzz[j] * fl1_fx + pa2pb_yyy_yzz[j]);

                t_yyy_yzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_yy[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_y[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_y[j] * fl1_fz * fl2_fx + 6.0 * pb_zz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyy_y[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_yy_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_yzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_yzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_yzz[j] * fl1_fz);

                t_yyy_zzz[j] = fl_s_0_0 * (2.25 * pa2pb_y_z[j] * fl2_fx + 1.5 * pa2pb_yyy_z[j] * fl1_fx + 1.5 * pa2pb_y_zzz[j] * fl1_fx + pa2pb_yyy_zzz[j]);

                t_yyy_zzz[j] += fl_r_0_0 * (-4.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyy_z[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_y_z[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yyy_z[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_zzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_y_zzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_zzz[j] * fl1_fz);
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(361 * idx + 19);

            auto pa2pb_y_y = pa2pbDistances.data(361 * idx + 20);

            auto pa2pb_y_z = pa2pbDistances.data(361 * idx + 21);

            auto pa2pb_z_x = pa2pbDistances.data(361 * idx + 38);

            auto pa2pb_z_y = pa2pbDistances.data(361 * idx + 39);

            auto pa2pb_z_z = pa2pbDistances.data(361 * idx + 40);

            auto pa2pb_z_xxx = pa2pbDistances.data(361 * idx + 47);

            auto pa2pb_z_xxy = pa2pbDistances.data(361 * idx + 48);

            auto pa2pb_z_xxz = pa2pbDistances.data(361 * idx + 49);

            auto pa2pb_z_xyy = pa2pbDistances.data(361 * idx + 50);

            auto pa2pb_z_xyz = pa2pbDistances.data(361 * idx + 51);

            auto pa2pb_z_xzz = pa2pbDistances.data(361 * idx + 52);

            auto pa2pb_z_yyy = pa2pbDistances.data(361 * idx + 53);

            auto pa2pb_z_yyz = pa2pbDistances.data(361 * idx + 54);

            auto pa2pb_z_yzz = pa2pbDistances.data(361 * idx + 55);

            auto pa2pb_z_zzz = pa2pbDistances.data(361 * idx + 56);

            auto pa2pb_yy_xx = pa2pbDistances.data(361 * idx + 117);

            auto pa2pb_yy_xy = pa2pbDistances.data(361 * idx + 118);

            auto pa2pb_yy_xz = pa2pbDistances.data(361 * idx + 119);

            auto pa2pb_yy_yy = pa2pbDistances.data(361 * idx + 120);

            auto pa2pb_yy_yz = pa2pbDistances.data(361 * idx + 121);

            auto pa2pb_yy_zz = pa2pbDistances.data(361 * idx + 122);

            auto pa2pb_yz_xx = pa2pbDistances.data(361 * idx + 136);

            auto pa2pb_yz_xy = pa2pbDistances.data(361 * idx + 137);

            auto pa2pb_yz_xz = pa2pbDistances.data(361 * idx + 138);

            auto pa2pb_yz_yy = pa2pbDistances.data(361 * idx + 139);

            auto pa2pb_yz_yz = pa2pbDistances.data(361 * idx + 140);

            auto pa2pb_yz_zz = pa2pbDistances.data(361 * idx + 141);

            auto pa2pb_yyz_x = pa2pbDistances.data(361 * idx + 304);

            auto pa2pb_yyz_y = pa2pbDistances.data(361 * idx + 305);

            auto pa2pb_yyz_z = pa2pbDistances.data(361 * idx + 306);

            auto pa2pb_yyz_xxx = pa2pbDistances.data(361 * idx + 313);

            auto pa2pb_yyz_xxy = pa2pbDistances.data(361 * idx + 314);

            auto pa2pb_yyz_xxz = pa2pbDistances.data(361 * idx + 315);

            auto pa2pb_yyz_xyy = pa2pbDistances.data(361 * idx + 316);

            auto pa2pb_yyz_xyz = pa2pbDistances.data(361 * idx + 317);

            auto pa2pb_yyz_xzz = pa2pbDistances.data(361 * idx + 318);

            auto pa2pb_yyz_yyy = pa2pbDistances.data(361 * idx + 319);

            auto pa2pb_yyz_yyz = pa2pbDistances.data(361 * idx + 320);

            auto pa2pb_yyz_yzz = pa2pbDistances.data(361 * idx + 321);

            auto pa2pb_yyz_zzz = pa2pbDistances.data(361 * idx + 322);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xy, \
                                     pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, pa2pb_yy_zz, pa2pb_yyz_x, pa2pb_yyz_xxx, \
                                     pa2pb_yyz_xxy, pa2pb_yyz_xxz, pa2pb_yyz_xyy, pa2pb_yyz_xyz, pa2pb_yyz_xzz, \
                                     pa2pb_yyz_y, pa2pb_yyz_yyy, pa2pb_yyz_yyz, pa2pb_yyz_yzz, pa2pb_yyz_z, \
                                     pa2pb_yyz_zzz, pa2pb_yz_xx, pa2pb_yz_xy, pa2pb_yz_xz, pa2pb_yz_yy, pa2pb_yz_yz, \
                                     pa2pb_yz_zz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_xxz, pa2pb_z_xyy, \
                                     pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, \
                                     pa2pb_z_z, pa2pb_z_zzz, pa_yy, pa_yz, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, \
                                     s_0_0, t_yyz_xxx, t_yyz_xxy, t_yyz_xxz, t_yyz_xyy, t_yyz_xyz, t_yyz_xzz, \
                                     t_yyz_yyy, t_yyz_yyz, t_yyz_yzz, t_yyz_zzz: VLX_ALIGN)
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

                t_yyz_xxx[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl2_fx + 1.5 * pa2pb_yyz_x[j] * fl1_fx + 0.5 * pa2pb_z_xxx[j] * fl1_fx + pa2pb_yyz_xxx[j]);

                t_yyz_xxx[j] += fl_r_0_0 * (-1.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_x[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_yyz_x[j] * fl1_fz * fl1_fx - pa2pb_z_xxx[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xxx[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_xxx[j] * fl1_fz);

                t_yyz_xxy[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + 0.25 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa2pb_yyz_y[j] * fl1_fx + pa2pb_yz_xx[j] * fl1_fx + 0.5 * pa2pb_z_xxy[j] * fl1_fx + pa2pb_yyz_xxy[j]);

                t_yyz_xxy[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_y[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyz_y[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz - pa2pb_z_xxy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xxy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_xxy[j] * fl1_fz);

                t_yyz_xxz[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + 0.25 * pa2pb_z_z[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_yyz_z[j] * fl1_fx + 0.5 * pa2pb_yy_xx[j] * fl1_fx + 0.5 * pa2pb_z_xxz[j] * fl1_fx + pa2pb_yyz_xxz[j]);

                t_yyz_xxz[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_yy[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_z[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz + 2.0 * pb_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyz_z[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_yy_xx[j] * fl1_fz * fl1_fx - pa2pb_z_xxz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xxz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_xxz[j] * fl1_fz);

                t_yyz_xyy[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa2pb_yyz_x[j] * fl1_fx + 2.0 * pa2pb_yz_xy[j] * fl1_fx + 0.5 * pa2pb_z_xyy[j] * fl1_fx + pa2pb_yyz_xyy[j]);

                t_yyz_xyy[j] += fl_r_0_0 * (6.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_x[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_yyz_x[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz - pa2pb_z_xyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xyy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_xyy[j] * fl1_fz);

                t_yyz_xyz[j] = fl_s_0_0 * (0.5 * pa2pb_y_x[j] * fl2_fx + 0.25 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_yy_xy[j] * fl1_fx + pa2pb_yz_xz[j] * fl1_fx + 0.5 * pa2pb_z_xyz[j] * fl1_fx + pa2pb_yyz_xyz[j]);

                t_yyz_xyz[j] += fl_r_0_0 * (4.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz - 0.5 * pb_xy[j] * fl1_fz * fl1_fga * fl1_fx + 2.0 * pb_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yy_xy[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz - pa2pb_z_xyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xyz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_xyz[j] * fl1_fz);

                t_yyz_xzz[j] = fl_s_0_0 * (0.25 * pa2pb_z_x[j] * fl2_fx + 0.5 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_yyz_x[j] * fl1_fx + pa2pb_yy_xz[j] * fl1_fx + 0.5 * pa2pb_z_xzz[j] * fl1_fx + pa2pb_yyz_xzz[j]);

                t_yyz_xzz[j] += fl_r_0_0 * (-0.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - pb_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_x[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz + 4.0 * pb_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyz_x[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yy_xz[j] * fl1_fz * fl1_fx - pa2pb_z_xzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_xzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_xzz[j] * fl1_fz);

                t_yyz_yyy[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl2_fx + 2.25 * pa2pb_z_y[j] * fl2_fx + 1.5 * pa2pb_yyz_y[j] * fl1_fx + 3.0 * pa2pb_yz_yy[j] * fl1_fx + 0.5 * pa2pb_z_yyy[j] * fl1_fx + pa2pb_yyz_yyy[j]);

                t_yyz_yyy[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_yz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_yyz_y[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz - pa2pb_z_yyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_yyy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_yyy[j] * fl1_fz);

                t_yyz_yyz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.25 * pa_yy[j] * fl2_fx + pa2pb_y_y[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.25 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_yyz_z[j] * fl1_fx + 0.5 * pa2pb_yy_yy[j] * fl1_fx + 2.0 * pa2pb_yz_yz[j] * fl1_fx + 0.5 * pa2pb_z_yyz[j] * fl1_fx + pa2pb_yyz_yyz[j]);

                t_yyz_yyz[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_yy[j] * fl1_fz * fl2_fx + 8.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_z[j] * fl1_fz * fl1_fgb + 2.0 * pb_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyz_z[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_yy_yy[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz - pa2pb_z_yyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_yyz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_yyz[j] * fl1_fz);

                t_yyz_yzz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + pa2pb_y_z[j] * fl2_fx + 0.25 * pa2pb_z_y[j] * fl2_fx + 0.5 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_yyz_y[j] * fl1_fx + pa2pb_yy_yz[j] * fl1_fx + pa2pb_yz_zz[j] * fl1_fx + 0.5 * pa2pb_z_yzz[j] * fl1_fx + pa2pb_yyz_yzz[j]);

                t_yyz_yzz[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_yz[j] * fl2_fx * fl1_fz + 8.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - pb_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_y[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz + 4.0 * pb_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yyz_y[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yy_yz[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz - pa2pb_z_yzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_yzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_yzz[j] * fl1_fz);

                t_yyz_zzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_yy[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.75 * pb_zz[j] * fl2_fx + 1.5 * pa2pb_yyz_z[j] * fl1_fx + 1.5 * pa2pb_yy_zz[j] * fl1_fx + 0.5 * pa2pb_z_zzz[j] * fl1_fx + pa2pb_yyz_zzz[j]);

                t_yyz_zzz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_yy[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_z[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz + 6.0 * pb_zz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_yyz_z[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_yy_zz[j] * fl1_fz * fl1_fx - pa2pb_z_zzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_z_zzz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_zzz[j] * fl1_fz);
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(361 * idx + 19);

            auto pa2pb_y_y = pa2pbDistances.data(361 * idx + 20);

            auto pa2pb_y_z = pa2pbDistances.data(361 * idx + 21);

            auto pa2pb_y_xxx = pa2pbDistances.data(361 * idx + 28);

            auto pa2pb_y_xxy = pa2pbDistances.data(361 * idx + 29);

            auto pa2pb_y_xxz = pa2pbDistances.data(361 * idx + 30);

            auto pa2pb_y_xyy = pa2pbDistances.data(361 * idx + 31);

            auto pa2pb_y_xyz = pa2pbDistances.data(361 * idx + 32);

            auto pa2pb_y_xzz = pa2pbDistances.data(361 * idx + 33);

            auto pa2pb_y_yyy = pa2pbDistances.data(361 * idx + 34);

            auto pa2pb_y_yyz = pa2pbDistances.data(361 * idx + 35);

            auto pa2pb_y_yzz = pa2pbDistances.data(361 * idx + 36);

            auto pa2pb_y_zzz = pa2pbDistances.data(361 * idx + 37);

            auto pa2pb_z_x = pa2pbDistances.data(361 * idx + 38);

            auto pa2pb_z_y = pa2pbDistances.data(361 * idx + 39);

            auto pa2pb_z_z = pa2pbDistances.data(361 * idx + 40);

            auto pa2pb_yz_xx = pa2pbDistances.data(361 * idx + 136);

            auto pa2pb_yz_xy = pa2pbDistances.data(361 * idx + 137);

            auto pa2pb_yz_xz = pa2pbDistances.data(361 * idx + 138);

            auto pa2pb_yz_yy = pa2pbDistances.data(361 * idx + 139);

            auto pa2pb_yz_yz = pa2pbDistances.data(361 * idx + 140);

            auto pa2pb_yz_zz = pa2pbDistances.data(361 * idx + 141);

            auto pa2pb_zz_xx = pa2pbDistances.data(361 * idx + 155);

            auto pa2pb_zz_xy = pa2pbDistances.data(361 * idx + 156);

            auto pa2pb_zz_xz = pa2pbDistances.data(361 * idx + 157);

            auto pa2pb_zz_yy = pa2pbDistances.data(361 * idx + 158);

            auto pa2pb_zz_yz = pa2pbDistances.data(361 * idx + 159);

            auto pa2pb_zz_zz = pa2pbDistances.data(361 * idx + 160);

            auto pa2pb_yzz_x = pa2pbDistances.data(361 * idx + 323);

            auto pa2pb_yzz_y = pa2pbDistances.data(361 * idx + 324);

            auto pa2pb_yzz_z = pa2pbDistances.data(361 * idx + 325);

            auto pa2pb_yzz_xxx = pa2pbDistances.data(361 * idx + 332);

            auto pa2pb_yzz_xxy = pa2pbDistances.data(361 * idx + 333);

            auto pa2pb_yzz_xxz = pa2pbDistances.data(361 * idx + 334);

            auto pa2pb_yzz_xyy = pa2pbDistances.data(361 * idx + 335);

            auto pa2pb_yzz_xyz = pa2pbDistances.data(361 * idx + 336);

            auto pa2pb_yzz_xzz = pa2pbDistances.data(361 * idx + 337);

            auto pa2pb_yzz_yyy = pa2pbDistances.data(361 * idx + 338);

            auto pa2pb_yzz_yyz = pa2pbDistances.data(361 * idx + 339);

            auto pa2pb_yzz_yzz = pa2pbDistances.data(361 * idx + 340);

            auto pa2pb_yzz_zzz = pa2pbDistances.data(361 * idx + 341);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, \
                                     pa2pb_y_xyy, pa2pb_y_xyz, pa2pb_y_xzz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, \
                                     pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, pa2pb_yz_xx, pa2pb_yz_xy, pa2pb_yz_xz, \
                                     pa2pb_yz_yy, pa2pb_yz_yz, pa2pb_yz_zz, pa2pb_yzz_x, pa2pb_yzz_xxx, pa2pb_yzz_xxy, \
                                     pa2pb_yzz_xxz, pa2pb_yzz_xyy, pa2pb_yzz_xyz, pa2pb_yzz_xzz, pa2pb_yzz_y, \
                                     pa2pb_yzz_yyy, pa2pb_yzz_yyz, pa2pb_yzz_yzz, pa2pb_yzz_z, pa2pb_yzz_zzz, pa2pb_z_x, \
                                     pa2pb_z_y, pa2pb_z_z, pa2pb_zz_xx, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, \
                                     pa2pb_zz_yz, pa2pb_zz_zz, pa_yz, pa_zz, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, \
                                     s_0_0, t_yzz_xxx, t_yzz_xxy, t_yzz_xxz, t_yzz_xyy, t_yzz_xyz, t_yzz_xzz, \
                                     t_yzz_yyy, t_yzz_yyz, t_yzz_yzz, t_yzz_zzz: VLX_ALIGN)
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

                t_yzz_xxx[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl2_fx + 1.5 * pa2pb_yzz_x[j] * fl1_fx + 0.5 * pa2pb_y_xxx[j] * fl1_fx + pa2pb_yzz_xxx[j]);

                t_yzz_xxx[j] += fl_r_0_0 * (-1.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yzz_x[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yzz_x[j] * fl1_fz * fl1_fx - pa2pb_y_xxx[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xxx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xxx[j] * fl1_fz);

                t_yzz_xxy[j] = fl_s_0_0 * (0.125 * fl3_fx + 0.25 * pa_zz[j] * fl2_fx + 0.25 * pa2pb_y_y[j] * fl2_fx + 0.25 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_yzz_y[j] * fl1_fx + 0.5 * pa2pb_zz_xx[j] * fl1_fx + 0.5 * pa2pb_y_xxy[j] * fl1_fx + pa2pb_yzz_xxy[j]);

                t_yzz_xxy[j] += fl_r_0_0 * (-0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 0.75 * fl3_fx * fl1_fz + 2.0 * pa_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_y[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_y_y[j] * fl1_fz * fl2_fx + 2.0 * pb_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yzz_y[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz - pa2pb_y_xxy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xxy[j] * fl1_fz);

                t_yzz_xxz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + 0.25 * pa2pb_y_z[j] * fl2_fx + 0.5 * pa2pb_yzz_z[j] * fl1_fx + pa2pb_yz_xx[j] * fl1_fx + 0.5 * pa2pb_y_xxz[j] * fl1_fx + pa2pb_yzz_xxz[j]);

                t_yzz_xxz[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_yz[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yzz_z[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_y_z[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yzz_z[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yz_xx[j] * fl1_fz * fl1_fx - pa2pb_y_xxz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xxz[j] * fl1_fz);

                t_yzz_xyy[j] = fl_s_0_0 * (0.25 * pa2pb_y_x[j] * fl2_fx + 0.5 * pb_xy[j] * fl2_fx + 0.5 * pa2pb_yzz_x[j] * fl1_fx + pa2pb_zz_xy[j] * fl1_fx + 0.5 * pa2pb_y_xyy[j] * fl1_fx + pa2pb_yzz_xyy[j]);

                t_yzz_xyy[j] += fl_r_0_0 * (-0.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - pb_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_x[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_y_x[j] * fl1_fz * fl2_fx + 4.0 * pb_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yzz_x[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz - pa2pb_y_xyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xyy[j] * fl1_fz);

                t_yzz_xyz[j] = fl_s_0_0 * (0.5 * pa2pb_z_x[j] * fl2_fx + 0.25 * pb_xz[j] * fl2_fx + pa2pb_yz_xy[j] * fl1_fx + 0.5 * pa2pb_zz_xz[j] * fl1_fx + 0.5 * pa2pb_y_xyz[j] * fl1_fx + pa2pb_yzz_xyz[j]);

                t_yzz_xyz[j] += fl_r_0_0 * (4.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz - 0.5 * pb_xz[j] * fl1_fx * fl1_fz * fl1_fga + 2.0 * pb_xz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_yz_xy[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz - pa2pb_y_xyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xyz[j] * fl1_fz);

                t_yzz_xzz[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl2_fx + 0.5 * pa2pb_yzz_x[j] * fl1_fx + 2.0 * pa2pb_yz_xz[j] * fl1_fx + 0.5 * pa2pb_y_xzz[j] * fl1_fx + pa2pb_yzz_xzz[j]);

                t_yzz_xzz[j] += fl_r_0_0 * (6.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yzz_x[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_yzz_x[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_yz_xz[j] * fl1_fz * fl1_fx - pa2pb_y_xzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_xzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xzz[j] * fl1_fz);

                t_yzz_yyy[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.75 * pb_yy[j] * fl2_fx + 1.5 * pa2pb_yzz_y[j] * fl1_fx + 1.5 * pa2pb_zz_yy[j] * fl1_fx + 0.5 * pa2pb_y_yyy[j] * fl1_fx + pa2pb_yzz_yyy[j]);

                t_yzz_yyy[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzz_y[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_y[j] * fl1_fz * fl2_fx + 6.0 * pb_yy[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_yzz_y[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz - pa2pb_y_yyy[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_yyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_yyy[j] * fl1_fz);

                t_yzz_yyz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl2_fx + pa2pb_z_y[j] * fl2_fx + 0.25 * pa2pb_y_z[j] * fl2_fx + 0.5 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_yzz_z[j] * fl1_fx + pa2pb_yz_yy[j] * fl1_fx + pa2pb_zz_yz[j] * fl1_fx + 0.5 * pa2pb_y_yyz[j] * fl1_fx + pa2pb_yzz_yyz[j]);

                t_yzz_yyz[j] += fl_r_0_0 * (-pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 4.0 * pa_yz[j] * fl1_fz * fl2_fx + 8.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - pb_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_z[j] * fl1_fz * fl1_fgb + 2.0 * pa2pb_y_z[j] * fl1_fz * fl2_fx + 4.0 * pb_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yzz_z[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_yz_yy[j] * fl1_fz * fl1_fx + 10.0 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz - pa2pb_y_yyz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_yyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_yyz[j] * fl1_fz);

                t_yzz_yzz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa2pb_y_y[j] * fl2_fx + 0.25 * pa_zz[j] * fl2_fx + pa2pb_z_z[j] * fl2_fx + 0.25 * pb_zz[j] * fl2_fx + 0.5 * pa2pb_yzz_y[j] * fl1_fx + 2.0 * pa2pb_yz_yz[j] * fl1_fx + 0.5 * pa2pb_zz_zz[j] * fl1_fx + 0.5 * pa2pb_y_yzz[j] * fl1_fx + pa2pb_yzz_yzz[j]);

                t_yzz_yzz[j] += fl_r_0_0 * (2.25 * fl3_fx * fl1_fz - 0.25 * fl2_fx * fl1_fz * fl1_fgb - 0.25 * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz + 2.0 * pa_zz[j] * fl2_fx * fl1_fz + 8.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_y[j] * fl1_fz * fl1_fgb + 2.0 * pb_zz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yzz_y[j] * fl1_fz * fl1_fx + 20.0 * pa2pb_yz_yz[j] * fl1_fz * fl1_fx + 5.0 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz - pa2pb_y_yzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_yzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_yzz[j] * fl1_fz);

                t_yzz_zzz[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl2_fx + 2.25 * pa2pb_y_z[j] * fl2_fx + 1.5 * pa2pb_yzz_z[j] * fl1_fx + 3.0 * pa2pb_yz_zz[j] * fl1_fx + 0.5 * pa2pb_y_zzz[j] * fl1_fx + pa2pb_yzz_zzz[j]);

                t_yzz_zzz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_yz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yzz_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_yzz_z[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_yz_zz[j] * fl1_fz * fl1_fx - pa2pb_y_zzz[j] * fl1_fz * fl1_fga + 5.0 * pa2pb_y_zzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_zzz[j] * fl1_fz);
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
                                  const CMemBlock2D<double>& pa2pbDistances,
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

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_x = pa2pbDistances.data(361 * idx + 38);

            auto pa2pb_z_y = pa2pbDistances.data(361 * idx + 39);

            auto pa2pb_z_z = pa2pbDistances.data(361 * idx + 40);

            auto pa2pb_z_xxx = pa2pbDistances.data(361 * idx + 47);

            auto pa2pb_z_xxy = pa2pbDistances.data(361 * idx + 48);

            auto pa2pb_z_xxz = pa2pbDistances.data(361 * idx + 49);

            auto pa2pb_z_xyy = pa2pbDistances.data(361 * idx + 50);

            auto pa2pb_z_xyz = pa2pbDistances.data(361 * idx + 51);

            auto pa2pb_z_xzz = pa2pbDistances.data(361 * idx + 52);

            auto pa2pb_z_yyy = pa2pbDistances.data(361 * idx + 53);

            auto pa2pb_z_yyz = pa2pbDistances.data(361 * idx + 54);

            auto pa2pb_z_yzz = pa2pbDistances.data(361 * idx + 55);

            auto pa2pb_z_zzz = pa2pbDistances.data(361 * idx + 56);

            auto pa2pb_zz_xx = pa2pbDistances.data(361 * idx + 155);

            auto pa2pb_zz_xy = pa2pbDistances.data(361 * idx + 156);

            auto pa2pb_zz_xz = pa2pbDistances.data(361 * idx + 157);

            auto pa2pb_zz_yy = pa2pbDistances.data(361 * idx + 158);

            auto pa2pb_zz_yz = pa2pbDistances.data(361 * idx + 159);

            auto pa2pb_zz_zz = pa2pbDistances.data(361 * idx + 160);

            auto pa2pb_zzz_x = pa2pbDistances.data(361 * idx + 342);

            auto pa2pb_zzz_y = pa2pbDistances.data(361 * idx + 343);

            auto pa2pb_zzz_z = pa2pbDistances.data(361 * idx + 344);

            auto pa2pb_zzz_xxx = pa2pbDistances.data(361 * idx + 351);

            auto pa2pb_zzz_xxy = pa2pbDistances.data(361 * idx + 352);

            auto pa2pb_zzz_xxz = pa2pbDistances.data(361 * idx + 353);

            auto pa2pb_zzz_xyy = pa2pbDistances.data(361 * idx + 354);

            auto pa2pb_zzz_xyz = pa2pbDistances.data(361 * idx + 355);

            auto pa2pb_zzz_xzz = pa2pbDistances.data(361 * idx + 356);

            auto pa2pb_zzz_yyy = pa2pbDistances.data(361 * idx + 357);

            auto pa2pb_zzz_yyz = pa2pbDistances.data(361 * idx + 358);

            auto pa2pb_zzz_yzz = pa2pbDistances.data(361 * idx + 359);

            auto pa2pb_zzz_zzz = pa2pbDistances.data(361 * idx + 360);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_xxz, \
                                     pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, \
                                     pa2pb_z_yzz, pa2pb_z_z, pa2pb_z_zzz, pa2pb_zz_xx, pa2pb_zz_xy, pa2pb_zz_xz, \
                                     pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zz_zz, pa2pb_zzz_x, pa2pb_zzz_xxx, pa2pb_zzz_xxy, \
                                     pa2pb_zzz_xxz, pa2pb_zzz_xyy, pa2pb_zzz_xyz, pa2pb_zzz_xzz, pa2pb_zzz_y, \
                                     pa2pb_zzz_yyy, pa2pb_zzz_yyz, pa2pb_zzz_yzz, pa2pb_zzz_z, pa2pb_zzz_zzz, pa_zz, pb_xx, \
                                     pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, s_0_0, t_zzz_xxx, t_zzz_xxy, t_zzz_xxz, \
                                     t_zzz_xyy, t_zzz_xyz, t_zzz_xzz, t_zzz_yyy, t_zzz_yyz, t_zzz_yzz, t_zzz_zzz: VLX_ALIGN)
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

                t_zzz_xxx[j] = fl_s_0_0 * (2.25 * pa2pb_z_x[j] * fl2_fx + 1.5 * pa2pb_zzz_x[j] * fl1_fx + 1.5 * pa2pb_z_xxx[j] * fl1_fx + pa2pb_zzz_xxx[j]);

                t_zzz_xxx[j] += fl_r_0_0 * (-4.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zzz_x[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_z_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_zzz_x[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxx[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_xxx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xxx[j] * fl1_fz);

                t_zzz_xxy[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl2_fx + 0.5 * pa2pb_zzz_y[j] * fl1_fx + 1.5 * pa2pb_z_xxy[j] * fl1_fx + pa2pb_zzz_xxy[j]);

                t_zzz_xxy[j] += fl_r_0_0 * (-1.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_zzz_y[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_z_y[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_zzz_y[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xxy[j] * fl1_fz);

                t_zzz_xxz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.75 * pb_xx[j] * fl2_fx + 0.5 * pa2pb_zzz_z[j] * fl1_fx + 1.5 * pa2pb_zz_xx[j] * fl1_fx + 1.5 * pa2pb_z_xxz[j] * fl1_fx + pa2pb_zzz_xxz[j]);

                t_zzz_xxz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_zz[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_z[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_z_z[j] * fl1_fz * fl2_fx + 6.0 * pb_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_zzz_z[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_zz_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xxz[j] * fl1_fz);

                t_zzz_xyy[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl2_fx + 0.5 * pa2pb_zzz_x[j] * fl1_fx + 1.5 * pa2pb_z_xyy[j] * fl1_fx + pa2pb_zzz_xyy[j]);

                t_zzz_xyy[j] += fl_r_0_0 * (-1.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_zzz_x[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_z_x[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_zzz_x[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xyy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xyy[j] * fl1_fz);

                t_zzz_xyz[j] = fl_s_0_0 * (0.75 * pb_xy[j] * fl2_fx + 1.5 * pa2pb_zz_xy[j] * fl1_fx + 1.5 * pa2pb_z_xyz[j] * fl1_fx + pa2pb_zzz_xyz[j]);

                t_zzz_xyz[j] += fl_r_0_0 * (-1.5 * pb_xy[j] * fl1_fx * fl1_fz * fl1_fga + 6.0 * pb_xy[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_zz_xy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xyz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xyz[j] * fl1_fz);

                t_zzz_xzz[j] = fl_s_0_0 * (2.25 * pa2pb_z_x[j] * fl2_fx + 1.5 * pb_xz[j] * fl2_fx + 0.5 * pa2pb_zzz_x[j] * fl1_fx + 3.0 * pa2pb_zz_xz[j] * fl1_fx + 1.5 * pa2pb_z_xzz[j] * fl1_fx + pa2pb_zzz_xzz[j]);

                t_zzz_xzz[j] += fl_r_0_0 * (18.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_x[j] * fl1_fz * fl1_fgb + 12.0 * pb_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_zzz_x[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_zz_xz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_xzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xzz[j] * fl1_fz);

                t_zzz_yyy[j] = fl_s_0_0 * (2.25 * pa2pb_z_y[j] * fl2_fx + 1.5 * pa2pb_zzz_y[j] * fl1_fx + 1.5 * pa2pb_z_yyy[j] * fl1_fx + pa2pb_zzz_yyy[j]);

                t_zzz_yyy[j] += fl_r_0_0 * (-4.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zzz_y[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_z_y[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_zzz_y[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_yyy[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_yyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_yyy[j] * fl1_fz);

                t_zzz_yyz[j] = fl_s_0_0 * (0.375 * fl3_fx + 0.75 * pa_zz[j] * fl2_fx + 0.75 * pa2pb_z_z[j] * fl2_fx + 0.75 * pb_yy[j] * fl2_fx + 0.5 * pa2pb_zzz_z[j] * fl1_fx + 1.5 * pa2pb_zz_yy[j] * fl1_fx + 1.5 * pa2pb_z_yyz[j] * fl1_fx + pa2pb_zzz_yyz[j]);

                t_zzz_yyz[j] += fl_r_0_0 * (-0.75 * fl2_fx * fl1_fz * fl1_fgb - 0.75 * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.25 * fl3_fx * fl1_fz + 6.0 * pa_zz[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_z[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_z_z[j] * fl1_fz * fl2_fx + 6.0 * pb_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_zzz_z[j] * fl1_fz * fl1_fx + 15.0 * pa2pb_zz_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_yyz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_yyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_yyz[j] * fl1_fz);

                t_zzz_yzz[j] = fl_s_0_0 * (2.25 * pa2pb_z_y[j] * fl2_fx + 1.5 * pb_yz[j] * fl2_fx + 0.5 * pa2pb_zzz_y[j] * fl1_fx + 3.0 * pa2pb_zz_yz[j] * fl1_fx + 1.5 * pa2pb_z_yzz[j] * fl1_fx + pa2pb_zzz_yzz[j]);

                t_zzz_yzz[j] += fl_r_0_0 * (18.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_y[j] * fl1_fz * fl1_fgb + 12.0 * pb_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_zzz_y[j] * fl1_fz * fl1_fx + 30.0 * pa2pb_zz_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_yzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_yzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_yzz[j] * fl1_fz);

                t_zzz_zzz[j] = fl_s_0_0 * (1.875 * fl3_fx + 2.25 * pa_zz[j] * fl2_fx + 6.75 * pa2pb_z_z[j] * fl2_fx + 2.25 * pb_zz[j] * fl2_fx + 1.5 * pa2pb_zzz_z[j] * fl1_fx + 4.5 * pa2pb_zz_zz[j] * fl1_fx + 1.5 * pa2pb_z_zzz[j] * fl1_fx + pa2pb_zzz_zzz[j]);

                t_zzz_zzz[j] += fl_r_0_0 * (11.25 * fl3_fx * fl1_fz - 2.25 * fl2_fx * fl1_fz * fl1_fgb - 2.25 * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa_zz[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_z_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_zz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zzz_z[j] * fl1_fz * fl1_fgb + 18.0 * pb_zz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_zzz_z[j] * fl1_fz * fl1_fx + 45.0 * pa2pb_zz_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_zzz[j] * fl1_fz * fl1_fga + 15.0 * pa2pb_z_zzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }


} // kinrecfunc namespace

