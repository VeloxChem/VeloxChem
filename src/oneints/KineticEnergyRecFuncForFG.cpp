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
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForFG_0_5(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                               braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_5_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_10_15(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_15_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_20_25(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_25_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_30_35(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_35_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_40_45(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_45_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_50_55(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_55_60(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_60_65(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_65_70(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_70_75(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_75_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_80_85(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_85_90(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_90_95(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_95_100(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                  braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_100_105(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_105_110(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_110_115(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_115_120(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_120_125(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_125_130(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_130_135(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_135_140(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_140_145(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForFG_145_150(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compKineticEnergyForFG_0_5(      CMemBlock2D<double>& primBuffer,
                               const CMemBlock2D<double>& auxBuffer,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pbDistances,
                               const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(19 * idx);

            auto pa_xxx = paDistances.data(19 * idx + 9);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_xxxx = pa2pbDistances.data(646 * idx + 19);

            auto pa2pb_x_xxxy = pa2pbDistances.data(646 * idx + 20);

            auto pa2pb_x_xxxz = pa2pbDistances.data(646 * idx + 21);

            auto pa2pb_x_xxyy = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_x_xxyz = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 102);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 103);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xx_xxx = pa2pbDistances.data(646 * idx + 111);

            auto pa2pb_xx_xxy = pa2pbDistances.data(646 * idx + 112);

            auto pa2pb_xx_xxz = pa2pbDistances.data(646 * idx + 113);

            auto pa2pb_xx_xyy = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_xx_xyz = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_xxx_xx = pa2pbDistances.data(646 * idx + 309);

            auto pa2pb_xxx_xy = pa2pbDistances.data(646 * idx + 310);

            auto pa2pb_xxx_xz = pa2pbDistances.data(646 * idx + 311);

            auto pa2pb_xxx_yy = pa2pbDistances.data(646 * idx + 312);

            auto pa2pb_xxx_yz = pa2pbDistances.data(646 * idx + 313);

            auto pa2pb_xxx_xxxx = pa2pbDistances.data(646 * idx + 325);

            auto pa2pb_xxx_xxxy = pa2pbDistances.data(646 * idx + 326);

            auto pa2pb_xxx_xxxz = pa2pbDistances.data(646 * idx + 327);

            auto pa2pb_xxx_xxyy = pa2pbDistances.data(646 * idx + 328);

            auto pa2pb_xxx_xxyz = pa2pbDistances.data(646 * idx + 329);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxx_xxxx = primBuffer.data(150 * idx);

            auto t_xxx_xxxy = primBuffer.data(150 * idx + 1);

            auto t_xxx_xxxz = primBuffer.data(150 * idx + 2);

            auto t_xxx_xxyy = primBuffer.data(150 * idx + 3);

            auto t_xxx_xxyz = primBuffer.data(150 * idx + 4);

            // Batch of Integrals (0,5)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xxxx, pa2pb_x_xxxy, pa2pb_x_xxxz, \
                                     pa2pb_x_xxyy, pa2pb_x_xxyz, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_xx_x, pa2pb_xx_xxx, pa2pb_xx_xxy, pa2pb_xx_xxz, pa2pb_xx_xyy, pa2pb_xx_xyz, \
                                     pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxx_xx, pa2pb_xxx_xxxx, pa2pb_xxx_xxxy, \
                                     pa2pb_xxx_xxxz, pa2pb_xxx_xxyy, pa2pb_xxx_xxyz, pa2pb_xxx_xy, pa2pb_xxx_xz, \
                                     pa2pb_xxx_yy, pa2pb_xxx_yz, pa_x, pa_xxx, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_y, \
                                     pb_z, r_0_0, s_0_0, t_xxx_xxxx, t_xxx_xxxy, t_xxx_xxxz, t_xxx_xxyy, t_xxx_xxyz: VLX_ALIGN)
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

                t_xxx_xxxx[j] = fl_s_0_0 * (5.625 * pa_x[j] * fl3_fx + 7.5 * pb_x[j] * fl3_fx + 0.75 * pa_xxx[j] * fl2_fx + 9.0 * pa2pb_xx_x[j] * fl2_fx + 13.5 * pa2pb_x_xx[j] * fl2_fx + 3.0 * pb_xxx[j] * fl2_fx + 3.0 * pa2pb_xxx_xx[j] * fl1_fx + 6.0 * pa2pb_xx_xxx[j] * fl1_fx + 1.5 * pa2pb_x_xxxx[j] * fl1_fx + pa2pb_xxx_xxxx[j]);

                t_xxx_xxxx[j] += fl_r_0_0 * (-13.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 45.0 * pa_x[j] * fl3_fx * fl1_fz + 60.0 * pb_x[j] * fl3_fx * fl1_fz - 2.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xxx[j] * fl1_fz * fl2_fx + 90.0 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 135.0 * pa2pb_x_xx[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xxx_xx[j] * fl1_fz * fl1_fgb + 30.0 * pb_xxx[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xxx_xx[j] * fl1_fz * fl1_fx + 72.0 * pa2pb_xx_xxx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xxxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xxxx[j] * fl1_fz);

                t_xxx_xxxy[j] = fl_s_0_0 * (1.875 * pb_y[j] * fl3_fx + 2.25 * pa2pb_xx_y[j] * fl2_fx + 6.75 * pa2pb_x_xy[j] * fl2_fx + 2.25 * pb_xxy[j] * fl2_fx + 1.5 * pa2pb_xxx_xy[j] * fl1_fx + 4.5 * pa2pb_xx_xxy[j] * fl1_fx + 1.5 * pa2pb_x_xxxy[j] * fl1_fx + pa2pb_xxx_xxxy[j]);

                t_xxx_xxxy[j] += fl_r_0_0 * (15.0 * pb_y[j] * fl3_fx * fl1_fz - 2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_x_xy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxx_xy[j] * fl1_fz * fl1_fgb + 22.5 * pb_xxy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxx_xy[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_xx_xxy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xxxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xxxy[j] * fl1_fz);

                t_xxx_xxxz[j] = fl_s_0_0 * (1.875 * pb_z[j] * fl3_fx + 2.25 * pa2pb_xx_z[j] * fl2_fx + 6.75 * pa2pb_x_xz[j] * fl2_fx + 2.25 * pb_xxz[j] * fl2_fx + 1.5 * pa2pb_xxx_xz[j] * fl1_fx + 4.5 * pa2pb_xx_xxz[j] * fl1_fx + 1.5 * pa2pb_x_xxxz[j] * fl1_fx + pa2pb_xxx_xxxz[j]);

                t_xxx_xxxz[j] += fl_r_0_0 * (15.0 * pb_z[j] * fl3_fx * fl1_fz - 2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_x_xz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxx_xz[j] * fl1_fz * fl1_fgb + 22.5 * pb_xxz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxx_xz[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_xx_xxz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xxxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xxxz[j] * fl1_fz);

                t_xxx_xxyy[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pb_x[j] * fl3_fx + 0.25 * pa_xxx[j] * fl2_fx + 1.5 * pa2pb_xx_x[j] * fl2_fx + 2.25 * pa2pb_x_yy[j] * fl2_fx + 0.75 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_xxx_xx[j] * fl1_fx + 0.5 * pa2pb_xxx_yy[j] * fl1_fx + 3.0 * pa2pb_xx_xyy[j] * fl1_fx + 1.5 * pa2pb_x_xxyy[j] * fl1_fx + pa2pb_xxx_xxyy[j]);

                t_xxx_xxyy[j] += fl_r_0_0 * (-3.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxx_yy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xx[j] * fl1_fz * fl2_fx + 15.0 * pb_xyy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxx_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxx_yy[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xx_xyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xxyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xxyy[j] * fl1_fz);

                t_xxx_xxyz[j] = fl_s_0_0 * (2.25 * pa2pb_x_yz[j] * fl2_fx + 1.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_xxx_yz[j] * fl1_fx + 3.0 * pa2pb_xx_xyz[j] * fl1_fx + 1.5 * pa2pb_x_xxyz[j] * fl1_fx + pa2pb_xxx_xxyz[j]);

                t_xxx_xxyz[j] += fl_r_0_0 * (22.5 * pa2pb_x_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_yz[j] * fl1_fz * fl1_fgb + 15.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxx_yz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xx_xyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xxyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_5_10(      CMemBlock2D<double>& primBuffer,
                                const CMemBlock2D<double>& auxBuffer,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pbDistances,
                                const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(19 * idx);

            auto pa_xxx = paDistances.data(19 * idx + 9);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_x_xxzz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_x_xyyy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_x_xyyz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_x_xyzz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_x_xzzz = pa2pbDistances.data(646 * idx + 28);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 102);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 103);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xx_xzz = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_xx_yyy = pa2pbDistances.data(646 * idx + 117);

            auto pa2pb_xx_yyz = pa2pbDistances.data(646 * idx + 118);

            auto pa2pb_xx_yzz = pa2pbDistances.data(646 * idx + 119);

            auto pa2pb_xx_zzz = pa2pbDistances.data(646 * idx + 120);

            auto pa2pb_xxx_xx = pa2pbDistances.data(646 * idx + 309);

            auto pa2pb_xxx_xy = pa2pbDistances.data(646 * idx + 310);

            auto pa2pb_xxx_xz = pa2pbDistances.data(646 * idx + 311);

            auto pa2pb_xxx_zz = pa2pbDistances.data(646 * idx + 314);

            auto pa2pb_xxx_xxzz = pa2pbDistances.data(646 * idx + 330);

            auto pa2pb_xxx_xyyy = pa2pbDistances.data(646 * idx + 331);

            auto pa2pb_xxx_xyyz = pa2pbDistances.data(646 * idx + 332);

            auto pa2pb_xxx_xyzz = pa2pbDistances.data(646 * idx + 333);

            auto pa2pb_xxx_xzzz = pa2pbDistances.data(646 * idx + 334);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxx_xxzz = primBuffer.data(150 * idx + 5);

            auto t_xxx_xyyy = primBuffer.data(150 * idx + 6);

            auto t_xxx_xyyz = primBuffer.data(150 * idx + 7);

            auto t_xxx_xyzz = primBuffer.data(150 * idx + 8);

            auto t_xxx_xzzz = primBuffer.data(150 * idx + 9);

            // Batch of Integrals (5,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xxzz, pa2pb_x_xy, pa2pb_x_xyyy, \
                                     pa2pb_x_xyyz, pa2pb_x_xyzz, pa2pb_x_xz, pa2pb_x_xzzz, pa2pb_x_zz, pa2pb_xx_x, \
                                     pa2pb_xx_xzz, pa2pb_xx_y, pa2pb_xx_yyy, pa2pb_xx_yyz, pa2pb_xx_yzz, pa2pb_xx_z, \
                                     pa2pb_xx_zzz, pa2pb_xxx_xx, pa2pb_xxx_xxzz, pa2pb_xxx_xy, pa2pb_xxx_xyyy, \
                                     pa2pb_xxx_xyyz, pa2pb_xxx_xyzz, pa2pb_xxx_xz, pa2pb_xxx_xzzz, pa2pb_xxx_zz, pa_x, \
                                     pa_xxx, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, s_0_0, \
                                     t_xxx_xxzz, t_xxx_xyyy, t_xxx_xyyz, t_xxx_xyzz, t_xxx_xzzz: VLX_ALIGN)
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

                t_xxx_xxzz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pb_x[j] * fl3_fx + 0.25 * pa_xxx[j] * fl2_fx + 1.5 * pa2pb_xx_x[j] * fl2_fx + 2.25 * pa2pb_x_zz[j] * fl2_fx + 0.75 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_xxx_xx[j] * fl1_fx + 0.5 * pa2pb_xxx_zz[j] * fl1_fx + 3.0 * pa2pb_xx_xzz[j] * fl1_fx + 1.5 * pa2pb_x_xxzz[j] * fl1_fx + pa2pb_xxx_xxzz[j]);

                t_xxx_xxzz[j] += fl_r_0_0 * (-3.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxx_zz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xx[j] * fl1_fz * fl2_fx + 15.0 * pb_xzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxx_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxx_zz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xx_xzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xxzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xxzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xxzz[j] * fl1_fz);

                t_xxx_xyyy[j] = fl_s_0_0 * (1.125 * pb_y[j] * fl3_fx + 2.25 * pa2pb_xx_y[j] * fl2_fx + 2.25 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pb_yyy[j] * fl2_fx + 1.5 * pa2pb_xxx_xy[j] * fl1_fx + 1.5 * pa2pb_xx_yyy[j] * fl1_fx + 1.5 * pa2pb_x_xyyy[j] * fl1_fx + pa2pb_xxx_xyyy[j]);

                t_xxx_xyyy[j] += fl_r_0_0 * (-2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pb_y[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxx_xy[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_x_xy[j] * fl1_fz * fl2_fx + 7.5 * pb_yyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxx_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xx_yyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xyyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xyyy[j] * fl1_fz);

                t_xxx_xyyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_xx_z[j] * fl2_fx + 0.75 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_xxx_xz[j] * fl1_fx + 1.5 * pa2pb_xx_yyz[j] * fl1_fx + 1.5 * pa2pb_x_xyyz[j] * fl1_fx + pa2pb_xxx_xyyz[j]);

                t_xxx_xyyz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xz[j] * fl1_fz * fl2_fx + 7.5 * pb_yyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxx_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xx_yyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xyyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xyyz[j] * fl1_fz);

                t_xxx_xyzz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_xx_y[j] * fl2_fx + 0.75 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_xxx_xy[j] * fl1_fx + 1.5 * pa2pb_xx_yzz[j] * fl1_fx + 1.5 * pa2pb_x_xyzz[j] * fl1_fx + pa2pb_xxx_xyzz[j]);

                t_xxx_xyzz[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxx_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xy[j] * fl1_fz * fl2_fx + 7.5 * pb_yzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxx_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xx_yzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xyzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xyzz[j] * fl1_fz);

                t_xxx_xzzz[j] = fl_s_0_0 * (1.125 * pb_z[j] * fl3_fx + 2.25 * pa2pb_xx_z[j] * fl2_fx + 2.25 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pb_zzz[j] * fl2_fx + 1.5 * pa2pb_xxx_xz[j] * fl1_fx + 1.5 * pa2pb_xx_zzz[j] * fl1_fx + 1.5 * pa2pb_x_xzzz[j] * fl1_fx + pa2pb_xxx_xzzz[j]);

                t_xxx_xzzz[j] += fl_r_0_0 * (-2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pb_z[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxx_xz[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_x_xz[j] * fl1_fz * fl2_fx + 7.5 * pb_zzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxx_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xx_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_xzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_xzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_10_15(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(19 * idx);

            auto pa_xxx = paDistances.data(19 * idx + 9);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_x_yyyy = pa2pbDistances.data(646 * idx + 29);

            auto pa2pb_x_yyyz = pa2pbDistances.data(646 * idx + 30);

            auto pa2pb_x_yyzz = pa2pbDistances.data(646 * idx + 31);

            auto pa2pb_x_yzzz = pa2pbDistances.data(646 * idx + 32);

            auto pa2pb_x_zzzz = pa2pbDistances.data(646 * idx + 33);

            auto pa2pb_xxx_yy = pa2pbDistances.data(646 * idx + 312);

            auto pa2pb_xxx_yz = pa2pbDistances.data(646 * idx + 313);

            auto pa2pb_xxx_zz = pa2pbDistances.data(646 * idx + 314);

            auto pa2pb_xxx_yyyy = pa2pbDistances.data(646 * idx + 335);

            auto pa2pb_xxx_yyyz = pa2pbDistances.data(646 * idx + 336);

            auto pa2pb_xxx_yyzz = pa2pbDistances.data(646 * idx + 337);

            auto pa2pb_xxx_yzzz = pa2pbDistances.data(646 * idx + 338);

            auto pa2pb_xxx_zzzz = pa2pbDistances.data(646 * idx + 339);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxx_yyyy = primBuffer.data(150 * idx + 10);

            auto t_xxx_yyyz = primBuffer.data(150 * idx + 11);

            auto t_xxx_yyzz = primBuffer.data(150 * idx + 12);

            auto t_xxx_yzzz = primBuffer.data(150 * idx + 13);

            auto t_xxx_zzzz = primBuffer.data(150 * idx + 14);

            // Batch of Integrals (10,15)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_yy, pa2pb_x_yyyy, pa2pb_x_yyyz, pa2pb_x_yyzz, \
                                     pa2pb_x_yz, pa2pb_x_yzzz, pa2pb_x_zz, pa2pb_x_zzzz, pa2pb_xxx_yy, pa2pb_xxx_yyyy, \
                                     pa2pb_xxx_yyyz, pa2pb_xxx_yyzz, pa2pb_xxx_yz, pa2pb_xxx_yzzz, pa2pb_xxx_zz, \
                                     pa2pb_xxx_zzzz, pa_x, pa_xxx, r_0_0, s_0_0, t_xxx_yyyy, t_xxx_yyyz, t_xxx_yyzz, t_xxx_yzzz, \
                                     t_xxx_zzzz: VLX_ALIGN)
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

                t_xxx_yyyy[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pa_xxx[j] * fl2_fx + 4.5 * pa2pb_x_yy[j] * fl2_fx + 3.0 * pa2pb_xxx_yy[j] * fl1_fx + 1.5 * pa2pb_x_yyyy[j] * fl1_fx + pa2pb_xxx_yyyy[j]);

                t_xxx_yyyy[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl1_fz * fl3_fx + 7.5 * pa_xxx[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxx_yy[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_x_yy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xxx_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_yyyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_yyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_yyyy[j] * fl1_fz);

                t_xxx_yyyz[j] = fl_s_0_0 * (2.25 * pa2pb_x_yz[j] * fl2_fx + 1.5 * pa2pb_xxx_yz[j] * fl1_fx + 1.5 * pa2pb_x_yyyz[j] * fl1_fx + pa2pb_xxx_yyyz[j]);

                t_xxx_yyyz[j] += fl_r_0_0 * (-4.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxx_yz[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_x_yz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xxx_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_yyyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_yyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_yyyz[j] * fl1_fz);

                t_xxx_yyzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xxx[j] * fl2_fx + 0.75 * pa2pb_x_yy[j] * fl2_fx + 0.75 * pa2pb_x_zz[j] * fl2_fx + 0.5 * pa2pb_xxx_yy[j] * fl1_fx + 0.5 * pa2pb_xxx_zz[j] * fl1_fx + 1.5 * pa2pb_x_yyzz[j] * fl1_fx + pa2pb_xxx_yyzz[j]);

                t_xxx_yyzz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxx_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxx_zz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_yy[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_zz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxx_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxx_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_yyzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_yyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_yyzz[j] * fl1_fz);

                t_xxx_yzzz[j] = fl_s_0_0 * (2.25 * pa2pb_x_yz[j] * fl2_fx + 1.5 * pa2pb_xxx_yz[j] * fl1_fx + 1.5 * pa2pb_x_yzzz[j] * fl1_fx + pa2pb_xxx_yzzz[j]);

                t_xxx_yzzz[j] += fl_r_0_0 * (-4.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxx_yz[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_x_yz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xxx_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_yzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_yzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_yzzz[j] * fl1_fz);

                t_xxx_zzzz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pa_xxx[j] * fl2_fx + 4.5 * pa2pb_x_zz[j] * fl2_fx + 3.0 * pa2pb_xxx_zz[j] * fl1_fx + 1.5 * pa2pb_x_zzzz[j] * fl1_fx + pa2pb_xxx_zzzz[j]);

                t_xxx_zzzz[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl1_fz * fl3_fx + 7.5 * pa_xxx[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxx_zz[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_x_zz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xxx_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_x_zzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_x_zzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_15_20(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_xxxx = pa2pbDistances.data(646 * idx + 53);

            auto pa2pb_y_xxxy = pa2pbDistances.data(646 * idx + 54);

            auto pa2pb_y_xxxz = pa2pbDistances.data(646 * idx + 55);

            auto pa2pb_y_xxyy = pa2pbDistances.data(646 * idx + 56);

            auto pa2pb_y_xxyz = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 102);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 103);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xx_xxx = pa2pbDistances.data(646 * idx + 111);

            auto pa2pb_xx_xxy = pa2pbDistances.data(646 * idx + 112);

            auto pa2pb_xx_xxz = pa2pbDistances.data(646 * idx + 113);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 136);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 137);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 138);

            auto pa2pb_xy_xxx = pa2pbDistances.data(646 * idx + 145);

            auto pa2pb_xy_xxy = pa2pbDistances.data(646 * idx + 146);

            auto pa2pb_xy_xxz = pa2pbDistances.data(646 * idx + 147);

            auto pa2pb_xy_xyy = pa2pbDistances.data(646 * idx + 148);

            auto pa2pb_xy_xyz = pa2pbDistances.data(646 * idx + 149);

            auto pa2pb_xxy_xx = pa2pbDistances.data(646 * idx + 343);

            auto pa2pb_xxy_xy = pa2pbDistances.data(646 * idx + 344);

            auto pa2pb_xxy_xz = pa2pbDistances.data(646 * idx + 345);

            auto pa2pb_xxy_yy = pa2pbDistances.data(646 * idx + 346);

            auto pa2pb_xxy_yz = pa2pbDistances.data(646 * idx + 347);

            auto pa2pb_xxy_xxxx = pa2pbDistances.data(646 * idx + 359);

            auto pa2pb_xxy_xxxy = pa2pbDistances.data(646 * idx + 360);

            auto pa2pb_xxy_xxxz = pa2pbDistances.data(646 * idx + 361);

            auto pa2pb_xxy_xxyy = pa2pbDistances.data(646 * idx + 362);

            auto pa2pb_xxy_xxyz = pa2pbDistances.data(646 * idx + 363);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxy_xxxx = primBuffer.data(150 * idx + 15);

            auto t_xxy_xxxy = primBuffer.data(150 * idx + 16);

            auto t_xxy_xxxz = primBuffer.data(150 * idx + 17);

            auto t_xxy_xxyy = primBuffer.data(150 * idx + 18);

            auto t_xxy_xxyz = primBuffer.data(150 * idx + 19);

            // Batch of Integrals (15,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_xx_x, \
                                     pa2pb_xx_xxx, pa2pb_xx_xxy, pa2pb_xx_xxz, pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxy_xx, \
                                     pa2pb_xxy_xxxx, pa2pb_xxy_xxxy, pa2pb_xxy_xxxz, pa2pb_xxy_xxyy, pa2pb_xxy_xxyz, \
                                     pa2pb_xxy_xy, pa2pb_xxy_xz, pa2pb_xxy_yy, pa2pb_xxy_yz, pa2pb_xy_x, pa2pb_xy_xxx, \
                                     pa2pb_xy_xxy, pa2pb_xy_xxz, pa2pb_xy_xyy, pa2pb_xy_xyz, pa2pb_xy_y, pa2pb_xy_z, \
                                     pa2pb_y_xx, pa2pb_y_xxxx, pa2pb_y_xxxy, pa2pb_y_xxxz, pa2pb_y_xxyy, pa2pb_y_xxyz, \
                                     pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, pa_x, pa_xxy, pa_y, pb_x, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_y, pb_z, r_0_0, s_0_0, t_xxy_xxxx, t_xxy_xxxy, t_xxy_xxxz, t_xxy_xxyy, \
                                     t_xxy_xxyz: VLX_ALIGN)
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

                t_xxy_xxxx[j] = fl_s_0_0 * (1.875 * pa_y[j] * fl3_fx + 0.75 * pa_xxy[j] * fl2_fx + 6.0 * pa2pb_xy_x[j] * fl2_fx + 4.5 * pa2pb_y_xx[j] * fl2_fx + 3.0 * pa2pb_xxy_xx[j] * fl1_fx + 4.0 * pa2pb_xy_xxx[j] * fl1_fx + 0.5 * pa2pb_y_xxxx[j] * fl1_fx + pa2pb_xxy_xxxx[j]);

                t_xxy_xxxx[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xxy[j] * fl1_fz * fl2_fx + 60.0 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 45.0 * pa2pb_y_xx[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxy_xx[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_xxy_xx[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_xy_xxx[j] * fl1_fx * fl1_fz - pa2pb_y_xxxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xxxx[j] * fl1_fz);

                t_xxy_xxxy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 1.125 * pb_x[j] * fl3_fx + 0.75 * pa2pb_xx_x[j] * fl2_fx + 1.5 * pa2pb_xy_y[j] * fl2_fx + 1.5 * pa2pb_x_xx[j] * fl2_fx + 2.25 * pa2pb_y_xy[j] * fl2_fx + 0.25 * pb_xxx[j] * fl2_fx + 1.5 * pa2pb_xxy_xy[j] * fl1_fx + 0.5 * pa2pb_xx_xxx[j] * fl1_fx + 3.0 * pa2pb_xy_xxy[j] * fl1_fx + 0.5 * pa2pb_y_xxxy[j] * fl1_fx + pa2pb_xxy_xxxy[j]);

                t_xxy_xxxy[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz + 9.0 * pb_x[j] * fl3_fx * fl1_fz - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xy_y[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_y_xy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_xy[j] * fl1_fz * fl1_fgb + 2.5 * pb_xxx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxy_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xx_xxx[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xy_xxy[j] * fl1_fx * fl1_fz - pa2pb_y_xxxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xxxy[j] * fl1_fz);

                t_xxy_xxxz[j] = fl_s_0_0 * (1.5 * pa2pb_xy_z[j] * fl2_fx + 2.25 * pa2pb_y_xz[j] * fl2_fx + 1.5 * pa2pb_xxy_xz[j] * fl1_fx + 3.0 * pa2pb_xy_xxz[j] * fl1_fx + 0.5 * pa2pb_y_xxxz[j] * fl1_fx + pa2pb_xxy_xxxz[j]);

                t_xxy_xxxz[j] += fl_r_0_0 * (-3.0 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xy_z[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_y_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xxy_xz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xy_xxz[j] * fl1_fx * fl1_fz - pa2pb_y_xxxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xxxz[j] * fl1_fz);

                t_xxy_xxyy[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pb_y[j] * fl3_fx + 0.25 * pa_xxy[j] * fl2_fx + 0.5 * pa2pb_xx_y[j] * fl2_fx + pa2pb_xy_x[j] * fl2_fx + 2.0 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pa2pb_y_yy[j] * fl2_fx + 0.25 * pa2pb_y_xx[j] * fl2_fx + 0.5 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_xxy_xx[j] * fl1_fx + 0.5 * pa2pb_xxy_yy[j] * fl1_fx + pa2pb_xx_xxy[j] * fl1_fx + 2.0 * pa2pb_xy_xyy[j] * fl1_fx + 0.5 * pa2pb_y_xxyy[j] * fl1_fx + pa2pb_xxy_xxyy[j]);

                t_xxy_xxyy[j] += fl_r_0_0 * (-pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 6.0 * pb_y[j] * fl3_fx * fl1_fz - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xxy[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - pb_xxy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxy_yy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 5.0 * pb_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxy_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxy_yy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xx_xxy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xy_xyy[j] * fl1_fx * fl1_fz - pa2pb_y_xxyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xxyy[j] * fl1_fz);

                t_xxy_xxyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.25 * pa2pb_xx_z[j] * fl2_fx + pa2pb_x_xz[j] * fl2_fx + 0.75 * pa2pb_y_yz[j] * fl2_fx + 0.25 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_xxy_yz[j] * fl1_fx + 0.5 * pa2pb_xx_xxz[j] * fl1_fx + 2.0 * pa2pb_xy_xyz[j] * fl1_fx + 0.5 * pa2pb_y_xxyz[j] * fl1_fx + pa2pb_xxy_xxyz[j]);

                t_xxy_xxyz[j] += fl_r_0_0 * (3.0 * pb_z[j] * fl3_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_yz[j] * fl1_fz * fl1_fgb + 2.5 * pb_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxy_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xx_xxz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xy_xyz[j] * fl1_fx * fl1_fz - pa2pb_y_xxyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_20_25(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_y_xxzz = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_y_xyyy = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_y_xyyz = pa2pbDistances.data(646 * idx + 60);

            auto pa2pb_y_xyzz = pa2pbDistances.data(646 * idx + 61);

            auto pa2pb_y_xzzz = pa2pbDistances.data(646 * idx + 62);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 102);

            auto pa2pb_xx_xyy = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_xx_xyz = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_xx_xzz = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 136);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 137);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 138);

            auto pa2pb_xy_xzz = pa2pbDistances.data(646 * idx + 150);

            auto pa2pb_xy_yyy = pa2pbDistances.data(646 * idx + 151);

            auto pa2pb_xy_yyz = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_xy_yzz = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_xy_zzz = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_xxy_xx = pa2pbDistances.data(646 * idx + 343);

            auto pa2pb_xxy_xy = pa2pbDistances.data(646 * idx + 344);

            auto pa2pb_xxy_xz = pa2pbDistances.data(646 * idx + 345);

            auto pa2pb_xxy_zz = pa2pbDistances.data(646 * idx + 348);

            auto pa2pb_xxy_xxzz = pa2pbDistances.data(646 * idx + 364);

            auto pa2pb_xxy_xyyy = pa2pbDistances.data(646 * idx + 365);

            auto pa2pb_xxy_xyyz = pa2pbDistances.data(646 * idx + 366);

            auto pa2pb_xxy_xyzz = pa2pbDistances.data(646 * idx + 367);

            auto pa2pb_xxy_xzzz = pa2pbDistances.data(646 * idx + 368);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxy_xxzz = primBuffer.data(150 * idx + 20);

            auto t_xxy_xyyy = primBuffer.data(150 * idx + 21);

            auto t_xxy_xyyz = primBuffer.data(150 * idx + 22);

            auto t_xxy_xyzz = primBuffer.data(150 * idx + 23);

            auto t_xxy_xzzz = primBuffer.data(150 * idx + 24);

            // Batch of Integrals (20,25)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xx_x, \
                                     pa2pb_xx_xyy, pa2pb_xx_xyz, pa2pb_xx_xzz, pa2pb_xxy_xx, pa2pb_xxy_xxzz, \
                                     pa2pb_xxy_xy, pa2pb_xxy_xyyy, pa2pb_xxy_xyyz, pa2pb_xxy_xyzz, pa2pb_xxy_xz, \
                                     pa2pb_xxy_xzzz, pa2pb_xxy_zz, pa2pb_xy_x, pa2pb_xy_xzz, pa2pb_xy_y, pa2pb_xy_yyy, \
                                     pa2pb_xy_yyz, pa2pb_xy_yzz, pa2pb_xy_z, pa2pb_xy_zzz, pa2pb_y_xx, pa2pb_y_xxzz, \
                                     pa2pb_y_xy, pa2pb_y_xyyy, pa2pb_y_xyyz, pa2pb_y_xyzz, pa2pb_y_xz, pa2pb_y_xzzz, \
                                     pa2pb_y_zz, pa_x, pa_xxy, pa_y, pb_x, pb_xyy, pb_xyz, pb_xzz, r_0_0, s_0_0, t_xxy_xxzz, \
                                     t_xxy_xyyy, t_xxy_xyyz, t_xxy_xyzz, t_xxy_xzzz: VLX_ALIGN)
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

                t_xxy_xxzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.25 * pa_xxy[j] * fl2_fx + pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_y_zz[j] * fl2_fx + 0.25 * pa2pb_y_xx[j] * fl2_fx + 0.5 * pa2pb_xxy_xx[j] * fl1_fx + 0.5 * pa2pb_xxy_zz[j] * fl1_fx + 2.0 * pa2pb_xy_xzz[j] * fl1_fx + 0.5 * pa2pb_y_xxzz[j] * fl1_fx + pa2pb_xxy_xxzz[j]);

                t_xxy_xxzz[j] += fl_r_0_0 * (-pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xxy[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxy_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxy_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxy_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xy_xzz[j] * fl1_fx * fl1_fz - pa2pb_y_xxzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xxzz[j] * fl1_fz);

                t_xxy_xyyy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_xx_x[j] * fl2_fx + 1.5 * pa2pb_xy_y[j] * fl2_fx + 1.5 * pa2pb_x_yy[j] * fl2_fx + 0.75 * pa2pb_y_xy[j] * fl2_fx + 0.75 * pb_xyy[j] * fl2_fx + 1.5 * pa2pb_xxy_xy[j] * fl1_fx + 1.5 * pa2pb_xx_xyy[j] * fl1_fx + pa2pb_xy_yyy[j] * fl1_fx + 0.5 * pa2pb_y_xyyy[j] * fl1_fx + pa2pb_xxy_xyyy[j]);

                t_xxy_xyyy[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xy_y[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_x_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 7.5 * pb_xyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxy_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xx_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_yyy[j] * fl1_fx * fl1_fz - pa2pb_y_xyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xyyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xyyy[j] * fl1_fz);

                t_xxy_xyyz[j] = fl_s_0_0 * (0.5 * pa2pb_xy_z[j] * fl2_fx + pa2pb_x_yz[j] * fl2_fx + 0.25 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_xxy_xz[j] * fl1_fx + pa2pb_xx_xyz[j] * fl1_fx + pa2pb_xy_yyz[j] * fl1_fx + 0.5 * pa2pb_y_xyyz[j] * fl1_fx + pa2pb_xxy_xyyz[j]);

                t_xxy_xyyz[j] += fl_r_0_0 * (-pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_xy_z[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - pb_xyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_xz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 5.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxy_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xx_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_yyz[j] * fl1_fx * fl1_fz - pa2pb_y_xyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xyyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xyyz[j] * fl1_fz);

                t_xxy_xyzz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl3_fx + 0.125 * pb_x[j] * fl3_fx + 0.25 * pa2pb_xx_x[j] * fl2_fx + 0.5 * pa2pb_xy_y[j] * fl2_fx + 0.5 * pa2pb_x_zz[j] * fl2_fx + 0.25 * pa2pb_y_xy[j] * fl2_fx + 0.25 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_xxy_xy[j] * fl1_fx + 0.5 * pa2pb_xx_xzz[j] * fl1_fx + pa2pb_xy_yzz[j] * fl1_fx + 0.5 * pa2pb_y_xyzz[j] * fl1_fx + pa2pb_xxy_xyzz[j]);

                t_xxy_xyzz[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb + pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xy_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_x_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xzz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_xy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 2.5 * pb_xzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxy_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xx_xzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_yzz[j] * fl1_fx * fl1_fz - pa2pb_y_xyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xyzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xyzz[j] * fl1_fz);

                t_xxy_xzzz[j] = fl_s_0_0 * (1.5 * pa2pb_xy_z[j] * fl2_fx + 0.75 * pa2pb_y_xz[j] * fl2_fx + 1.5 * pa2pb_xxy_xz[j] * fl1_fx + pa2pb_xy_zzz[j] * fl1_fx + 0.5 * pa2pb_y_xzzz[j] * fl1_fx + pa2pb_xxy_xzzz[j]);

                t_xxy_xzzz[j] += fl_r_0_0 * (-3.0 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xy_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxy_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_zzz[j] * fl1_fx * fl1_fz - pa2pb_y_xzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_25_30(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_y_yyyy = pa2pbDistances.data(646 * idx + 63);

            auto pa2pb_y_yyyz = pa2pbDistances.data(646 * idx + 64);

            auto pa2pb_y_yyzz = pa2pbDistances.data(646 * idx + 65);

            auto pa2pb_y_yzzz = pa2pbDistances.data(646 * idx + 66);

            auto pa2pb_y_zzzz = pa2pbDistances.data(646 * idx + 67);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 103);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xx_yyy = pa2pbDistances.data(646 * idx + 117);

            auto pa2pb_xx_yyz = pa2pbDistances.data(646 * idx + 118);

            auto pa2pb_xx_yzz = pa2pbDistances.data(646 * idx + 119);

            auto pa2pb_xx_zzz = pa2pbDistances.data(646 * idx + 120);

            auto pa2pb_xxy_yy = pa2pbDistances.data(646 * idx + 346);

            auto pa2pb_xxy_yz = pa2pbDistances.data(646 * idx + 347);

            auto pa2pb_xxy_zz = pa2pbDistances.data(646 * idx + 348);

            auto pa2pb_xxy_yyyy = pa2pbDistances.data(646 * idx + 369);

            auto pa2pb_xxy_yyyz = pa2pbDistances.data(646 * idx + 370);

            auto pa2pb_xxy_yyzz = pa2pbDistances.data(646 * idx + 371);

            auto pa2pb_xxy_yzzz = pa2pbDistances.data(646 * idx + 372);

            auto pa2pb_xxy_zzzz = pa2pbDistances.data(646 * idx + 373);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxy_yyyy = primBuffer.data(150 * idx + 25);

            auto t_xxy_yyyz = primBuffer.data(150 * idx + 26);

            auto t_xxy_yyzz = primBuffer.data(150 * idx + 27);

            auto t_xxy_yzzz = primBuffer.data(150 * idx + 28);

            auto t_xxy_zzzz = primBuffer.data(150 * idx + 29);

            // Batch of Integrals (25,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_xx_y, pa2pb_xx_yyy, pa2pb_xx_yyz, pa2pb_xx_yzz, \
                                     pa2pb_xx_z, pa2pb_xx_zzz, pa2pb_xxy_yy, pa2pb_xxy_yyyy, pa2pb_xxy_yyyz, \
                                     pa2pb_xxy_yyzz, pa2pb_xxy_yz, pa2pb_xxy_yzzz, pa2pb_xxy_zz, pa2pb_xxy_zzzz, \
                                     pa2pb_y_yy, pa2pb_y_yyyy, pa2pb_y_yyyz, pa2pb_y_yyzz, pa2pb_y_yz, pa2pb_y_yzzz, \
                                     pa2pb_y_zz, pa2pb_y_zzzz, pa_xxy, pa_y, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, \
                                     s_0_0, t_xxy_yyyy, t_xxy_yyyz, t_xxy_yyzz, t_xxy_yzzz, t_xxy_zzzz: VLX_ALIGN)
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

                t_xxy_yyyy[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 1.5 * pb_y[j] * fl3_fx + 0.75 * pa_xxy[j] * fl2_fx + 3.0 * pa2pb_xx_y[j] * fl2_fx + 1.5 * pa2pb_y_yy[j] * fl2_fx + pb_yyy[j] * fl2_fx + 3.0 * pa2pb_xxy_yy[j] * fl1_fx + 2.0 * pa2pb_xx_yyy[j] * fl1_fx + 0.5 * pa2pb_y_yyyy[j] * fl1_fx + pa2pb_xxy_yyyy[j]);

                t_xxy_yyyy[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 12.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa_xxy[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_xx_y[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pb_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxy_yy[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 10.0 * pb_yyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xxy_yy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xx_yyy[j] * fl1_fz * fl1_fx - pa2pb_y_yyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_yyyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_yyyy[j] * fl1_fz);

                t_xxy_yyyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_xx_z[j] * fl2_fx + 0.75 * pa2pb_y_yz[j] * fl2_fx + 0.75 * pb_yyz[j] * fl2_fx + 1.5 * pa2pb_xxy_yz[j] * fl1_fx + 1.5 * pa2pb_xx_yyz[j] * fl1_fx + 0.5 * pa2pb_y_yyyz[j] * fl1_fx + pa2pb_xxy_yyyz[j]);

                t_xxy_yyyz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yyz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 7.5 * pb_yyz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxy_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xx_yyz[j] * fl1_fz * fl1_fx - pa2pb_y_yyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_yyyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_yyyz[j] * fl1_fz);

                t_xxy_yyzz[j] = fl_s_0_0 * (0.125 * pa_y[j] * fl3_fx + 0.25 * pb_y[j] * fl3_fx + 0.25 * pa_xxy[j] * fl2_fx + 0.5 * pa2pb_xx_y[j] * fl2_fx + 0.25 * pa2pb_y_yy[j] * fl2_fx + 0.25 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_xxy_yy[j] * fl1_fx + 0.5 * pa2pb_xxy_zz[j] * fl1_fx + pa2pb_xx_yzz[j] * fl1_fx + 0.5 * pa2pb_y_yyzz[j] * fl1_fx + pa2pb_xxy_yyzz[j]);

                t_xxy_yyzz[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb + pa_y[j] * fl3_fx * fl1_fz + 2.0 * pb_y[j] * fl3_fx * fl1_fz + 2.5 * pa_xxy[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xx_y[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - pb_yzz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxy_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxy_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 5.0 * pb_yzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxy_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxy_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xx_yzz[j] * fl1_fz * fl1_fx - pa2pb_y_yyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_yyzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_yyzz[j] * fl1_fz);

                t_xxy_yzzz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_xx_z[j] * fl2_fx + 0.75 * pa2pb_y_yz[j] * fl2_fx + 0.25 * pb_zzz[j] * fl2_fx + 1.5 * pa2pb_xxy_yz[j] * fl1_fx + 0.5 * pa2pb_xx_zzz[j] * fl1_fx + 0.5 * pa2pb_y_yzzz[j] * fl1_fx + pa2pb_xxy_yzzz[j]);

                t_xxy_yzzz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxy_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 2.5 * pb_zzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxy_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xx_zzz[j] * fl1_fz * fl1_fx - pa2pb_y_yzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_yzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_yzzz[j] * fl1_fz);

                t_xxy_zzzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa_xxy[j] * fl2_fx + 1.5 * pa2pb_y_zz[j] * fl2_fx + 3.0 * pa2pb_xxy_zz[j] * fl1_fx + 0.5 * pa2pb_y_zzzz[j] * fl1_fx + pa2pb_xxy_zzzz[j]);

                t_xxy_zzzz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 7.5 * pa_xxy[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxy_zz[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xxy_zz[j] * fl1_fz * fl1_fx - pa2pb_y_zzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_zzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxy_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_30_35(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(19 * idx);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_xxxx = pa2pbDistances.data(646 * idx + 87);

            auto pa2pb_z_xxxy = pa2pbDistances.data(646 * idx + 88);

            auto pa2pb_z_xxxz = pa2pbDistances.data(646 * idx + 89);

            auto pa2pb_z_xxyy = pa2pbDistances.data(646 * idx + 90);

            auto pa2pb_z_xxyz = pa2pbDistances.data(646 * idx + 91);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 102);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 103);

            auto pa2pb_xx_xxx = pa2pbDistances.data(646 * idx + 111);

            auto pa2pb_xx_xxy = pa2pbDistances.data(646 * idx + 112);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 171);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 172);

            auto pa2pb_xz_xxx = pa2pbDistances.data(646 * idx + 179);

            auto pa2pb_xz_xxy = pa2pbDistances.data(646 * idx + 180);

            auto pa2pb_xz_xxz = pa2pbDistances.data(646 * idx + 181);

            auto pa2pb_xz_xyy = pa2pbDistances.data(646 * idx + 182);

            auto pa2pb_xz_xyz = pa2pbDistances.data(646 * idx + 183);

            auto pa2pb_xxz_xx = pa2pbDistances.data(646 * idx + 377);

            auto pa2pb_xxz_xy = pa2pbDistances.data(646 * idx + 378);

            auto pa2pb_xxz_xz = pa2pbDistances.data(646 * idx + 379);

            auto pa2pb_xxz_yy = pa2pbDistances.data(646 * idx + 380);

            auto pa2pb_xxz_yz = pa2pbDistances.data(646 * idx + 381);

            auto pa2pb_xxz_xxxx = pa2pbDistances.data(646 * idx + 393);

            auto pa2pb_xxz_xxxy = pa2pbDistances.data(646 * idx + 394);

            auto pa2pb_xxz_xxxz = pa2pbDistances.data(646 * idx + 395);

            auto pa2pb_xxz_xxyy = pa2pbDistances.data(646 * idx + 396);

            auto pa2pb_xxz_xxyz = pa2pbDistances.data(646 * idx + 397);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxz_xxxx = primBuffer.data(150 * idx + 30);

            auto t_xxz_xxxy = primBuffer.data(150 * idx + 31);

            auto t_xxz_xxxz = primBuffer.data(150 * idx + 32);

            auto t_xxz_xxyy = primBuffer.data(150 * idx + 33);

            auto t_xxz_xxyz = primBuffer.data(150 * idx + 34);

            // Batch of Integrals (30,35)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_xx_x, pa2pb_xx_xxx, \
                                     pa2pb_xx_xxy, pa2pb_xx_y, pa2pb_xxz_xx, pa2pb_xxz_xxxx, pa2pb_xxz_xxxy, \
                                     pa2pb_xxz_xxxz, pa2pb_xxz_xxyy, pa2pb_xxz_xxyz, pa2pb_xxz_xy, pa2pb_xxz_xz, \
                                     pa2pb_xxz_yy, pa2pb_xxz_yz, pa2pb_xz_x, pa2pb_xz_xxx, pa2pb_xz_xxy, pa2pb_xz_xxz, \
                                     pa2pb_xz_xyy, pa2pb_xz_xyz, pa2pb_xz_y, pa2pb_xz_z, pa2pb_z_xx, pa2pb_z_xxxx, \
                                     pa2pb_z_xxxy, pa2pb_z_xxxz, pa2pb_z_xxyy, pa2pb_z_xxyz, pa2pb_z_xy, pa2pb_z_xz, \
                                     pa2pb_z_yy, pa2pb_z_yz, pa_x, pa_xxz, pa_z, pb_x, pb_xxx, pb_xxy, pb_y, r_0_0, s_0_0, \
                                     t_xxz_xxxx, t_xxz_xxxy, t_xxz_xxxz, t_xxz_xxyy, t_xxz_xxyz: VLX_ALIGN)
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

                t_xxz_xxxx[j] = fl_s_0_0 * (1.875 * pa_z[j] * fl3_fx + 0.75 * pa_xxz[j] * fl2_fx + 6.0 * pa2pb_xz_x[j] * fl2_fx + 4.5 * pa2pb_z_xx[j] * fl2_fx + 3.0 * pa2pb_xxz_xx[j] * fl1_fx + 4.0 * pa2pb_xz_xxx[j] * fl1_fx + 0.5 * pa2pb_z_xxxx[j] * fl1_fx + pa2pb_xxz_xxxx[j]);

                t_xxz_xxxx[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xxz[j] * fl1_fz * fl2_fx + 60.0 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 45.0 * pa2pb_z_xx[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxz_xx[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_xxz_xx[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_xz_xxx[j] * fl1_fx * fl1_fz - pa2pb_z_xxxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xxxx[j] * fl1_fz);

                t_xxz_xxxy[j] = fl_s_0_0 * (1.5 * pa2pb_xz_y[j] * fl2_fx + 2.25 * pa2pb_z_xy[j] * fl2_fx + 1.5 * pa2pb_xxz_xy[j] * fl1_fx + 3.0 * pa2pb_xz_xxy[j] * fl1_fx + 0.5 * pa2pb_z_xxxy[j] * fl1_fx + pa2pb_xxz_xxxy[j]);

                t_xxz_xxxy[j] += fl_r_0_0 * (-3.0 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_z_xy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xxz_xy[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xz_xxy[j] * fl1_fx * fl1_fz - pa2pb_z_xxxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xxxy[j] * fl1_fz);

                t_xxz_xxxz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 1.125 * pb_x[j] * fl3_fx + 0.75 * pa2pb_xx_x[j] * fl2_fx + 1.5 * pa2pb_xz_z[j] * fl2_fx + 1.5 * pa2pb_x_xx[j] * fl2_fx + 2.25 * pa2pb_z_xz[j] * fl2_fx + 0.25 * pb_xxx[j] * fl2_fx + 1.5 * pa2pb_xxz_xz[j] * fl1_fx + 0.5 * pa2pb_xx_xxx[j] * fl1_fx + 3.0 * pa2pb_xz_xxz[j] * fl1_fx + 0.5 * pa2pb_z_xxxz[j] * fl1_fx + pa2pb_xxz_xxxz[j]);

                t_xxz_xxxz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz + 9.0 * pb_x[j] * fl3_fx * fl1_fz - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_xz[j] * fl1_fz * fl1_fgb + 2.5 * pb_xxx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxz_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xx_xxx[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xz_xxz[j] * fl1_fx * fl1_fz - pa2pb_z_xxxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xxxz[j] * fl1_fz);

                t_xxz_xxyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.25 * pa_xxz[j] * fl2_fx + pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_z_yy[j] * fl2_fx + 0.25 * pa2pb_z_xx[j] * fl2_fx + 0.5 * pa2pb_xxz_xx[j] * fl1_fx + 0.5 * pa2pb_xxz_yy[j] * fl1_fx + 2.0 * pa2pb_xz_xyy[j] * fl1_fx + 0.5 * pa2pb_z_xxyy[j] * fl1_fx + pa2pb_xxz_xxyy[j]);

                t_xxz_xxyy[j] += fl_r_0_0 * (-pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz - 0.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xxz[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxz_yy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxz_yy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xz_xyy[j] * fl1_fx * fl1_fz - pa2pb_z_xxyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xxyy[j] * fl1_fz);

                t_xxz_xxyz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.25 * pa2pb_xx_y[j] * fl2_fx + pa2pb_x_xy[j] * fl2_fx + 0.75 * pa2pb_z_yz[j] * fl2_fx + 0.25 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_xxz_yz[j] * fl1_fx + 0.5 * pa2pb_xx_xxy[j] * fl1_fx + 2.0 * pa2pb_xz_xyz[j] * fl1_fx + 0.5 * pa2pb_z_xxyz[j] * fl1_fx + pa2pb_xxz_xxyz[j]);

                t_xxz_xxyz[j] += fl_r_0_0 * (3.0 * pb_y[j] * fl3_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xxy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_yz[j] * fl1_fz * fl1_fgb + 2.5 * pb_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxz_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xx_xxy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xz_xyz[j] * fl1_fx * fl1_fz - pa2pb_z_xxyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_35_40(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(19 * idx);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_z_xxzz = pa2pbDistances.data(646 * idx + 92);

            auto pa2pb_z_xyyy = pa2pbDistances.data(646 * idx + 93);

            auto pa2pb_z_xyyz = pa2pbDistances.data(646 * idx + 94);

            auto pa2pb_z_xyzz = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_z_xzzz = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 102);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xx_xxz = pa2pbDistances.data(646 * idx + 113);

            auto pa2pb_xx_xyy = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_xx_xyz = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_xx_xzz = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 171);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 172);

            auto pa2pb_xz_xzz = pa2pbDistances.data(646 * idx + 184);

            auto pa2pb_xz_yyy = pa2pbDistances.data(646 * idx + 185);

            auto pa2pb_xz_yyz = pa2pbDistances.data(646 * idx + 186);

            auto pa2pb_xz_yzz = pa2pbDistances.data(646 * idx + 187);

            auto pa2pb_xz_zzz = pa2pbDistances.data(646 * idx + 188);

            auto pa2pb_xxz_xx = pa2pbDistances.data(646 * idx + 377);

            auto pa2pb_xxz_xy = pa2pbDistances.data(646 * idx + 378);

            auto pa2pb_xxz_xz = pa2pbDistances.data(646 * idx + 379);

            auto pa2pb_xxz_zz = pa2pbDistances.data(646 * idx + 382);

            auto pa2pb_xxz_xxzz = pa2pbDistances.data(646 * idx + 398);

            auto pa2pb_xxz_xyyy = pa2pbDistances.data(646 * idx + 399);

            auto pa2pb_xxz_xyyz = pa2pbDistances.data(646 * idx + 400);

            auto pa2pb_xxz_xyzz = pa2pbDistances.data(646 * idx + 401);

            auto pa2pb_xxz_xzzz = pa2pbDistances.data(646 * idx + 402);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxz_xxzz = primBuffer.data(150 * idx + 35);

            auto t_xxz_xyyy = primBuffer.data(150 * idx + 36);

            auto t_xxz_xyyz = primBuffer.data(150 * idx + 37);

            auto t_xxz_xyzz = primBuffer.data(150 * idx + 38);

            auto t_xxz_xzzz = primBuffer.data(150 * idx + 39);

            // Batch of Integrals (35,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xx_x, \
                                     pa2pb_xx_xxz, pa2pb_xx_xyy, pa2pb_xx_xyz, pa2pb_xx_xzz, pa2pb_xx_z, pa2pb_xxz_xx, \
                                     pa2pb_xxz_xxzz, pa2pb_xxz_xy, pa2pb_xxz_xyyy, pa2pb_xxz_xyyz, pa2pb_xxz_xyzz, \
                                     pa2pb_xxz_xz, pa2pb_xxz_xzzz, pa2pb_xxz_zz, pa2pb_xz_x, pa2pb_xz_xzz, pa2pb_xz_y, \
                                     pa2pb_xz_yyy, pa2pb_xz_yyz, pa2pb_xz_yzz, pa2pb_xz_z, pa2pb_xz_zzz, pa2pb_z_xx, \
                                     pa2pb_z_xxzz, pa2pb_z_xy, pa2pb_z_xyyy, pa2pb_z_xyyz, pa2pb_z_xyzz, pa2pb_z_xz, \
                                     pa2pb_z_xzzz, pa2pb_z_zz, pa_x, pa_xxz, pa_z, pb_x, pb_xxz, pb_xyy, pb_xyz, pb_xzz, pb_z, \
                                     r_0_0, s_0_0, t_xxz_xxzz, t_xxz_xyyy, t_xxz_xyyz, t_xxz_xyzz, t_xxz_xzzz: VLX_ALIGN)
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

                t_xxz_xxzz[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pb_z[j] * fl3_fx + 0.25 * pa_xxz[j] * fl2_fx + 0.5 * pa2pb_xx_z[j] * fl2_fx + pa2pb_xz_x[j] * fl2_fx + 2.0 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pa2pb_z_zz[j] * fl2_fx + 0.25 * pa2pb_z_xx[j] * fl2_fx + 0.5 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_xxz_xx[j] * fl1_fx + 0.5 * pa2pb_xxz_zz[j] * fl1_fx + pa2pb_xx_xxz[j] * fl1_fx + 2.0 * pa2pb_xz_xzz[j] * fl1_fx + 0.5 * pa2pb_z_xxzz[j] * fl1_fx + pa2pb_xxz_xxzz[j]);

                t_xxz_xxzz[j] += fl_r_0_0 * (-pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 6.0 * pb_z[j] * fl3_fx * fl1_fz - 0.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xxz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - pb_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxz_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 5.0 * pb_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxz_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xx_xxz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xz_xzz[j] * fl1_fx * fl1_fz - pa2pb_z_xxzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xxzz[j] * fl1_fz);

                t_xxz_xyyy[j] = fl_s_0_0 * (1.5 * pa2pb_xz_y[j] * fl2_fx + 0.75 * pa2pb_z_xy[j] * fl2_fx + 1.5 * pa2pb_xxz_xy[j] * fl1_fx + pa2pb_xz_yyy[j] * fl1_fx + 0.5 * pa2pb_z_xyyy[j] * fl1_fx + pa2pb_xxz_xyyy[j]);

                t_xxz_xyyy[j] += fl_r_0_0 * (-3.0 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxz_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_yyy[j] * fl1_fx * fl1_fz - pa2pb_z_xyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xyyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xyyy[j] * fl1_fz);

                t_xxz_xyyz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl3_fx + 0.125 * pb_x[j] * fl3_fx + 0.25 * pa2pb_xx_x[j] * fl2_fx + 0.5 * pa2pb_xz_z[j] * fl2_fx + 0.5 * pa2pb_x_yy[j] * fl2_fx + 0.25 * pa2pb_z_xz[j] * fl2_fx + 0.25 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_xxz_xz[j] * fl1_fx + 0.5 * pa2pb_xx_xyy[j] * fl1_fx + pa2pb_xz_yyz[j] * fl1_fx + 0.5 * pa2pb_z_xyyz[j] * fl1_fx + pa2pb_xxz_xyyz[j]);

                t_xxz_xyyz[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb + pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_x_yy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xyy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_xz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 2.5 * pb_xyy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxz_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xx_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_yyz[j] * fl1_fx * fl1_fz - pa2pb_z_xyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xyyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xyyz[j] * fl1_fz);

                t_xxz_xyzz[j] = fl_s_0_0 * (0.5 * pa2pb_xz_y[j] * fl2_fx + pa2pb_x_yz[j] * fl2_fx + 0.25 * pa2pb_z_xy[j] * fl2_fx + 0.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_xxz_xy[j] * fl1_fx + pa2pb_xx_xyz[j] * fl1_fx + pa2pb_xz_yzz[j] * fl1_fx + 0.5 * pa2pb_z_xyzz[j] * fl1_fx + pa2pb_xxz_xyzz[j]);

                t_xxz_xyzz[j] += fl_r_0_0 * (-pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - pb_xyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_xy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 5.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxz_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xx_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_yzz[j] * fl1_fx * fl1_fz - pa2pb_z_xyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xyzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xyzz[j] * fl1_fz);

                t_xxz_xzzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_xx_x[j] * fl2_fx + 1.5 * pa2pb_xz_z[j] * fl2_fx + 1.5 * pa2pb_x_zz[j] * fl2_fx + 0.75 * pa2pb_z_xz[j] * fl2_fx + 0.75 * pb_xzz[j] * fl2_fx + 1.5 * pa2pb_xxz_xz[j] * fl1_fx + 1.5 * pa2pb_xx_xzz[j] * fl1_fx + pa2pb_xz_zzz[j] * fl1_fx + 0.5 * pa2pb_z_xzzz[j] * fl1_fx + pa2pb_xxz_xzzz[j]);

                t_xxz_xzzz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_x_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 7.5 * pb_xzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxz_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xx_xzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_zzz[j] * fl1_fx * fl1_fz - pa2pb_z_xzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_40_45(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_z_yyyy = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_z_yyyz = pa2pbDistances.data(646 * idx + 98);

            auto pa2pb_z_yyzz = pa2pbDistances.data(646 * idx + 99);

            auto pa2pb_z_yzzz = pa2pbDistances.data(646 * idx + 100);

            auto pa2pb_z_zzzz = pa2pbDistances.data(646 * idx + 101);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 103);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xx_yyy = pa2pbDistances.data(646 * idx + 117);

            auto pa2pb_xx_yyz = pa2pbDistances.data(646 * idx + 118);

            auto pa2pb_xx_yzz = pa2pbDistances.data(646 * idx + 119);

            auto pa2pb_xx_zzz = pa2pbDistances.data(646 * idx + 120);

            auto pa2pb_xxz_yy = pa2pbDistances.data(646 * idx + 380);

            auto pa2pb_xxz_yz = pa2pbDistances.data(646 * idx + 381);

            auto pa2pb_xxz_zz = pa2pbDistances.data(646 * idx + 382);

            auto pa2pb_xxz_yyyy = pa2pbDistances.data(646 * idx + 403);

            auto pa2pb_xxz_yyyz = pa2pbDistances.data(646 * idx + 404);

            auto pa2pb_xxz_yyzz = pa2pbDistances.data(646 * idx + 405);

            auto pa2pb_xxz_yzzz = pa2pbDistances.data(646 * idx + 406);

            auto pa2pb_xxz_zzzz = pa2pbDistances.data(646 * idx + 407);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxz_yyyy = primBuffer.data(150 * idx + 40);

            auto t_xxz_yyyz = primBuffer.data(150 * idx + 41);

            auto t_xxz_yyzz = primBuffer.data(150 * idx + 42);

            auto t_xxz_yzzz = primBuffer.data(150 * idx + 43);

            auto t_xxz_zzzz = primBuffer.data(150 * idx + 44);

            // Batch of Integrals (40,45)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_xx_y, pa2pb_xx_yyy, pa2pb_xx_yyz, pa2pb_xx_yzz, \
                                     pa2pb_xx_z, pa2pb_xx_zzz, pa2pb_xxz_yy, pa2pb_xxz_yyyy, pa2pb_xxz_yyyz, \
                                     pa2pb_xxz_yyzz, pa2pb_xxz_yz, pa2pb_xxz_yzzz, pa2pb_xxz_zz, pa2pb_xxz_zzzz, \
                                     pa2pb_z_yy, pa2pb_z_yyyy, pa2pb_z_yyyz, pa2pb_z_yyzz, pa2pb_z_yz, pa2pb_z_yzzz, \
                                     pa2pb_z_zz, pa2pb_z_zzzz, pa_xxz, pa_z, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, \
                                     s_0_0, t_xxz_yyyy, t_xxz_yyyz, t_xxz_yyzz, t_xxz_yzzz, t_xxz_zzzz: VLX_ALIGN)
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

                t_xxz_yyyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa_xxz[j] * fl2_fx + 1.5 * pa2pb_z_yy[j] * fl2_fx + 3.0 * pa2pb_xxz_yy[j] * fl1_fx + 0.5 * pa2pb_z_yyyy[j] * fl1_fx + pa2pb_xxz_yyyy[j]);

                t_xxz_yyyy[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 7.5 * pa_xxz[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxz_yy[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xxz_yy[j] * fl1_fz * fl1_fx - pa2pb_z_yyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_yyyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_yyyy[j] * fl1_fz);

                t_xxz_yyyz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_xx_y[j] * fl2_fx + 0.75 * pa2pb_z_yz[j] * fl2_fx + 0.25 * pb_yyy[j] * fl2_fx + 1.5 * pa2pb_xxz_yz[j] * fl1_fx + 0.5 * pa2pb_xx_yyy[j] * fl1_fx + 0.5 * pa2pb_z_yyyz[j] * fl1_fx + pa2pb_xxz_yyyz[j]);

                t_xxz_yyyz[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 2.5 * pb_yyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxz_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xx_yyy[j] * fl1_fz * fl1_fx - pa2pb_z_yyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_yyyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_yyyz[j] * fl1_fz);

                t_xxz_yyzz[j] = fl_s_0_0 * (0.125 * pa_z[j] * fl3_fx + 0.25 * pb_z[j] * fl3_fx + 0.25 * pa_xxz[j] * fl2_fx + 0.5 * pa2pb_xx_z[j] * fl2_fx + 0.25 * pa2pb_z_yy[j] * fl2_fx + 0.25 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_xxz_yy[j] * fl1_fx + 0.5 * pa2pb_xxz_zz[j] * fl1_fx + pa2pb_xx_yyz[j] * fl1_fx + 0.5 * pa2pb_z_yyzz[j] * fl1_fx + pa2pb_xxz_yyzz[j]);

                t_xxz_yyzz[j] += fl_r_0_0 * (-0.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb + pa_z[j] * fl3_fx * fl1_fz + 2.0 * pb_z[j] * fl3_fx * fl1_fz + 2.5 * pa_xxz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xx_z[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - pb_yyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxz_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 5.0 * pb_yyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxz_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxz_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xx_yyz[j] * fl1_fz * fl1_fx - pa2pb_z_yyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_yyzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_yyzz[j] * fl1_fz);

                t_xxz_yzzz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_xx_y[j] * fl2_fx + 0.75 * pa2pb_z_yz[j] * fl2_fx + 0.75 * pb_yzz[j] * fl2_fx + 1.5 * pa2pb_xxz_yz[j] * fl1_fx + 1.5 * pa2pb_xx_yzz[j] * fl1_fx + 0.5 * pa2pb_z_yzzz[j] * fl1_fx + pa2pb_xxz_yzzz[j]);

                t_xxz_yzzz[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxz_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 7.5 * pb_yzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxz_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xx_yzz[j] * fl1_fz * fl1_fx - pa2pb_z_yzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_yzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_yzzz[j] * fl1_fz);

                t_xxz_zzzz[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 1.5 * pb_z[j] * fl3_fx + 0.75 * pa_xxz[j] * fl2_fx + 3.0 * pa2pb_xx_z[j] * fl2_fx + 1.5 * pa2pb_z_zz[j] * fl2_fx + pb_zzz[j] * fl2_fx + 3.0 * pa2pb_xxz_zz[j] * fl1_fx + 2.0 * pa2pb_xx_zzz[j] * fl1_fx + 0.5 * pa2pb_z_zzzz[j] * fl1_fx + pa2pb_xxz_zzzz[j]);

                t_xxz_zzzz[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 12.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa_xxz[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_xx_z[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pb_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxz_zz[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 10.0 * pb_zzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xxz_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xx_zzz[j] * fl1_fz * fl1_fx - pa2pb_z_zzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_zzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_45_50(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (45,50)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_xxxx = pa2pbDistances.data(646 * idx + 19);

            auto pa2pb_x_xxxy = pa2pbDistances.data(646 * idx + 20);

            auto pa2pb_x_xxxz = pa2pbDistances.data(646 * idx + 21);

            auto pa2pb_x_xxyy = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_x_xxyz = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 136);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 137);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 138);

            auto pa2pb_xy_xxx = pa2pbDistances.data(646 * idx + 145);

            auto pa2pb_xy_xxy = pa2pbDistances.data(646 * idx + 146);

            auto pa2pb_xy_xxz = pa2pbDistances.data(646 * idx + 147);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 204);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 205);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 206);

            auto pa2pb_yy_xxx = pa2pbDistances.data(646 * idx + 213);

            auto pa2pb_yy_xxy = pa2pbDistances.data(646 * idx + 214);

            auto pa2pb_yy_xxz = pa2pbDistances.data(646 * idx + 215);

            auto pa2pb_yy_xyy = pa2pbDistances.data(646 * idx + 216);

            auto pa2pb_yy_xyz = pa2pbDistances.data(646 * idx + 217);

            auto pa2pb_xyy_xx = pa2pbDistances.data(646 * idx + 411);

            auto pa2pb_xyy_xy = pa2pbDistances.data(646 * idx + 412);

            auto pa2pb_xyy_xz = pa2pbDistances.data(646 * idx + 413);

            auto pa2pb_xyy_yy = pa2pbDistances.data(646 * idx + 414);

            auto pa2pb_xyy_yz = pa2pbDistances.data(646 * idx + 415);

            auto pa2pb_xyy_xxxx = pa2pbDistances.data(646 * idx + 427);

            auto pa2pb_xyy_xxxy = pa2pbDistances.data(646 * idx + 428);

            auto pa2pb_xyy_xxxz = pa2pbDistances.data(646 * idx + 429);

            auto pa2pb_xyy_xxyy = pa2pbDistances.data(646 * idx + 430);

            auto pa2pb_xyy_xxyz = pa2pbDistances.data(646 * idx + 431);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyy_xxxx = primBuffer.data(150 * idx + 45);

            auto t_xyy_xxxy = primBuffer.data(150 * idx + 46);

            auto t_xyy_xxxz = primBuffer.data(150 * idx + 47);

            auto t_xyy_xxyy = primBuffer.data(150 * idx + 48);

            auto t_xyy_xxyz = primBuffer.data(150 * idx + 49);

            // Batch of Integrals (45,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xxxx, pa2pb_x_xxxy, pa2pb_x_xxxz, \
                                     pa2pb_x_xxyy, pa2pb_x_xxyz, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_xy_x, pa2pb_xy_xxx, pa2pb_xy_xxy, pa2pb_xy_xxz, pa2pb_xy_y, pa2pb_xy_z, \
                                     pa2pb_xyy_xx, pa2pb_xyy_xxxx, pa2pb_xyy_xxxy, pa2pb_xyy_xxxz, pa2pb_xyy_xxyy, \
                                     pa2pb_xyy_xxyz, pa2pb_xyy_xy, pa2pb_xyy_xz, pa2pb_xyy_yy, pa2pb_xyy_yz, pa2pb_y_xx, \
                                     pa2pb_y_xy, pa2pb_y_xz, pa2pb_yy_x, pa2pb_yy_xxx, pa2pb_yy_xxy, pa2pb_yy_xxz, \
                                     pa2pb_yy_xyy, pa2pb_yy_xyz, pa2pb_yy_y, pa2pb_yy_z, pa_x, pa_xyy, pa_y, pb_x, pb_xxx, \
                                     pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_y, pb_z, r_0_0, s_0_0, t_xyy_xxxx, t_xyy_xxxy, \
                                     t_xyy_xxxz, t_xyy_xxyy, t_xyy_xxyz: VLX_ALIGN)
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

                t_xyy_xxxx[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 1.5 * pb_x[j] * fl3_fx + 0.75 * pa_xyy[j] * fl2_fx + 3.0 * pa2pb_yy_x[j] * fl2_fx + 1.5 * pa2pb_x_xx[j] * fl2_fx + pb_xxx[j] * fl2_fx + 3.0 * pa2pb_xyy_xx[j] * fl1_fx + 2.0 * pa2pb_yy_xxx[j] * fl1_fx + 0.5 * pa2pb_x_xxxx[j] * fl1_fx + pa2pb_xyy_xxxx[j]);

                t_xyy_xxxx[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx - 6.0 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_yy_x[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xyy_xx[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_xx[j] * fl1_fz * fl2_fx + 10.0 * pb_xxx[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xyy_xx[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yy_xxx[j] * fl1_fx * fl1_fz - pa2pb_x_xxxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xxxx[j] * fl1_fz);

                t_xyy_xxxy[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 0.375 * pb_y[j] * fl3_fx + 1.5 * pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_yy_y[j] * fl2_fx + 1.5 * pa2pb_y_xx[j] * fl2_fx + 0.75 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pb_xxy[j] * fl2_fx + 1.5 * pa2pb_xyy_xy[j] * fl1_fx + pa2pb_xy_xxx[j] * fl1_fx + 1.5 * pa2pb_yy_xxy[j] * fl1_fx + 0.5 * pa2pb_x_xxxy[j] * fl1_fx + pa2pb_xyy_xxxy[j]);

                t_xyy_xxxy[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_yy_y[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyy_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xy[j] * fl1_fz * fl2_fx + 7.5 * pb_xxy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyy_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_xxx[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yy_xxy[j] * fl1_fx * fl1_fz - pa2pb_x_xxxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xxxy[j] * fl1_fz);

                t_xyy_xxxz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_yy_z[j] * fl2_fx + 0.75 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pb_xxz[j] * fl2_fx + 1.5 * pa2pb_xyy_xz[j] * fl1_fx + 1.5 * pa2pb_yy_xxz[j] * fl1_fx + 0.5 * pa2pb_x_xxxz[j] * fl1_fx + pa2pb_xyy_xxxz[j]);

                t_xyy_xxxz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyy_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xz[j] * fl1_fz * fl2_fx + 7.5 * pb_xxz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyy_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yy_xxz[j] * fl1_fx * fl1_fz - pa2pb_x_xxxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xxxz[j] * fl1_fz);

                t_xyy_xxyy[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pb_x[j] * fl3_fx + 0.25 * pa_xyy[j] * fl2_fx + pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_x_xx[j] * fl2_fx + 0.5 * pa2pb_yy_x[j] * fl2_fx + 2.0 * pa2pb_y_xy[j] * fl2_fx + 0.25 * pa2pb_x_yy[j] * fl2_fx + 0.5 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_xyy_xx[j] * fl1_fx + 0.5 * pa2pb_xyy_yy[j] * fl1_fx + 2.0 * pa2pb_xy_xxy[j] * fl1_fx + pa2pb_yy_xyy[j] * fl1_fx + 0.5 * pa2pb_x_xxyy[j] * fl1_fx + pa2pb_xyy_xxyy[j]);

                t_xyy_xxyy[j] += fl_r_0_0 * (-pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 6.0 * pb_x[j] * fl3_fx * fl1_fz - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xyy[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yy_x[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyy_yy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_yy[j] * fl1_fz * fl2_fx + 5.0 * pb_xyy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyy_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xyy_yy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xy_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yy_xyy[j] * fl1_fx * fl1_fz - pa2pb_x_xxyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xxyy[j] * fl1_fz);

                t_xyy_xxyz[j] = fl_s_0_0 * (0.5 * pa2pb_xy_z[j] * fl2_fx + pa2pb_y_xz[j] * fl2_fx + 0.25 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_xyy_yz[j] * fl1_fx + pa2pb_xy_xxz[j] * fl1_fx + pa2pb_yy_xyz[j] * fl1_fx + 0.5 * pa2pb_x_xxyz[j] * fl1_fx + pa2pb_xyy_xxyz[j]);

                t_xyy_xxyz[j] += fl_r_0_0 * (-pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_xy_z[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_yz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_yz[j] * fl1_fz * fl2_fx + 5.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyy_yz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yy_xyz[j] * fl1_fx * fl1_fz - pa2pb_x_xxyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_50_55(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (50,55)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_x_xxzz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_x_xyyy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_x_xyyz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_x_xyzz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_x_xzzz = pa2pbDistances.data(646 * idx + 28);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 136);

            auto pa2pb_xy_xyy = pa2pbDistances.data(646 * idx + 148);

            auto pa2pb_xy_xyz = pa2pbDistances.data(646 * idx + 149);

            auto pa2pb_xy_xzz = pa2pbDistances.data(646 * idx + 150);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 204);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 205);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 206);

            auto pa2pb_yy_xzz = pa2pbDistances.data(646 * idx + 218);

            auto pa2pb_yy_yyy = pa2pbDistances.data(646 * idx + 219);

            auto pa2pb_yy_yyz = pa2pbDistances.data(646 * idx + 220);

            auto pa2pb_yy_yzz = pa2pbDistances.data(646 * idx + 221);

            auto pa2pb_yy_zzz = pa2pbDistances.data(646 * idx + 222);

            auto pa2pb_xyy_xx = pa2pbDistances.data(646 * idx + 411);

            auto pa2pb_xyy_xy = pa2pbDistances.data(646 * idx + 412);

            auto pa2pb_xyy_xz = pa2pbDistances.data(646 * idx + 413);

            auto pa2pb_xyy_zz = pa2pbDistances.data(646 * idx + 416);

            auto pa2pb_xyy_xxzz = pa2pbDistances.data(646 * idx + 432);

            auto pa2pb_xyy_xyyy = pa2pbDistances.data(646 * idx + 433);

            auto pa2pb_xyy_xyyz = pa2pbDistances.data(646 * idx + 434);

            auto pa2pb_xyy_xyzz = pa2pbDistances.data(646 * idx + 435);

            auto pa2pb_xyy_xzzz = pa2pbDistances.data(646 * idx + 436);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyy_xxzz = primBuffer.data(150 * idx + 50);

            auto t_xyy_xyyy = primBuffer.data(150 * idx + 51);

            auto t_xyy_xyyz = primBuffer.data(150 * idx + 52);

            auto t_xyy_xyzz = primBuffer.data(150 * idx + 53);

            auto t_xyy_xzzz = primBuffer.data(150 * idx + 54);

            // Batch of Integrals (50,55)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xxzz, pa2pb_x_xy, pa2pb_x_xyyy, \
                                     pa2pb_x_xyyz, pa2pb_x_xyzz, pa2pb_x_xz, pa2pb_x_xzzz, pa2pb_x_zz, pa2pb_xy_x, \
                                     pa2pb_xy_xyy, pa2pb_xy_xyz, pa2pb_xy_xzz, pa2pb_xyy_xx, pa2pb_xyy_xxzz, \
                                     pa2pb_xyy_xy, pa2pb_xyy_xyyy, pa2pb_xyy_xyyz, pa2pb_xyy_xyzz, pa2pb_xyy_xz, \
                                     pa2pb_xyy_xzzz, pa2pb_xyy_zz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yy_x, \
                                     pa2pb_yy_xzz, pa2pb_yy_y, pa2pb_yy_yyy, pa2pb_yy_yyz, pa2pb_yy_yzz, pa2pb_yy_z, \
                                     pa2pb_yy_zzz, pa_x, pa_xyy, pa_y, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, \
                                     r_0_0, s_0_0, t_xyy_xxzz, t_xyy_xyyy, t_xyy_xyyz, t_xyy_xyzz, t_xyy_xzzz: VLX_ALIGN)
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

                t_xyy_xxzz[j] = fl_s_0_0 * (0.125 * pa_x[j] * fl3_fx + 0.25 * pb_x[j] * fl3_fx + 0.25 * pa_xyy[j] * fl2_fx + 0.5 * pa2pb_yy_x[j] * fl2_fx + 0.25 * pa2pb_x_xx[j] * fl2_fx + 0.25 * pa2pb_x_zz[j] * fl2_fx + 0.5 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_xyy_xx[j] * fl1_fx + 0.5 * pa2pb_xyy_zz[j] * fl1_fx + pa2pb_yy_xzz[j] * fl1_fx + 0.5 * pa2pb_x_xxzz[j] * fl1_fx + pa2pb_xyy_xxzz[j]);

                t_xyy_xxzz[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + pa_x[j] * fl1_fz * fl3_fx - pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xyy[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yy_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyy_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_xx[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_x_zz[j] * fl1_fz * fl2_fx + 5.0 * pb_xzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyy_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xyy_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yy_xzz[j] * fl1_fx * fl1_fz - pa2pb_x_xxzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xxzz[j] * fl1_fz);

                t_xyy_xyyy[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 1.125 * pb_y[j] * fl3_fx + 1.5 * pa2pb_xy_x[j] * fl2_fx + 2.25 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pa2pb_yy_y[j] * fl2_fx + 1.5 * pa2pb_y_yy[j] * fl2_fx + 0.25 * pb_yyy[j] * fl2_fx + 1.5 * pa2pb_xyy_xy[j] * fl1_fx + 3.0 * pa2pb_xy_xyy[j] * fl1_fx + 0.5 * pa2pb_yy_yyy[j] * fl1_fx + 0.5 * pa2pb_x_xyyy[j] * fl1_fx + pa2pb_xyy_xyyy[j]);

                t_xyy_xyyy[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl3_fx * fl1_fz + 9.0 * pb_y[j] * fl3_fx * fl1_fz - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_yy_y[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyy_xy[j] * fl1_fz * fl1_fgb + 2.5 * pb_yyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyy_xy[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xy_xyy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_yyy[j] * fl1_fx * fl1_fz - pa2pb_x_xyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xyyy[j] * fl1_fz);

                t_xyy_xyyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_x_xz[j] * fl2_fx + 0.25 * pa2pb_yy_z[j] * fl2_fx + pa2pb_y_yz[j] * fl2_fx + 0.25 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_xyy_xz[j] * fl1_fx + 2.0 * pa2pb_xy_xyz[j] * fl1_fx + 0.5 * pa2pb_yy_yyz[j] * fl1_fx + 0.5 * pa2pb_x_xyyz[j] * fl1_fx + pa2pb_xyy_xyyz[j]);

                t_xyy_xyyz[j] += fl_r_0_0 * (3.0 * pb_z[j] * fl3_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_yy_z[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_xz[j] * fl1_fz * fl1_fgb + 2.5 * pb_yyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyy_xz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xy_xyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_yyz[j] * fl1_fx * fl1_fz - pa2pb_x_xyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xyyz[j] * fl1_fz);

                t_xyy_xyzz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl3_fx + 0.125 * pb_y[j] * fl3_fx + 0.5 * pa2pb_xy_x[j] * fl2_fx + 0.25 * pa2pb_yy_y[j] * fl2_fx + 0.5 * pa2pb_y_zz[j] * fl2_fx + 0.25 * pa2pb_x_xy[j] * fl2_fx + 0.25 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_xyy_xy[j] * fl1_fx + pa2pb_xy_xzz[j] * fl1_fx + 0.5 * pa2pb_yy_yzz[j] * fl1_fx + 0.5 * pa2pb_x_xyzz[j] * fl1_fx + pa2pb_xyy_xyzz[j]);

                t_xyy_xyzz[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb + pb_y[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_yy_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyy_xy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_xy[j] * fl1_fz * fl2_fx + 2.5 * pb_yzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyy_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_xzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_yzz[j] * fl1_fx * fl1_fz - pa2pb_x_xyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xyzz[j] * fl1_fz);

                t_xyy_xzzz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_yy_z[j] * fl2_fx + 0.75 * pa2pb_x_xz[j] * fl2_fx + 0.25 * pb_zzz[j] * fl2_fx + 1.5 * pa2pb_xyy_xz[j] * fl1_fx + 0.5 * pa2pb_yy_zzz[j] * fl1_fx + 0.5 * pa2pb_x_xzzz[j] * fl1_fx + pa2pb_xyy_xzzz[j]);

                t_xyy_xzzz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyy_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xz[j] * fl1_fz * fl2_fx + 2.5 * pb_zzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyy_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_zzz[j] * fl1_fx * fl1_fz - pa2pb_x_xzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_55_60(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (55,60)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_x_yyyy = pa2pbDistances.data(646 * idx + 29);

            auto pa2pb_x_yyyz = pa2pbDistances.data(646 * idx + 30);

            auto pa2pb_x_yyzz = pa2pbDistances.data(646 * idx + 31);

            auto pa2pb_x_yzzz = pa2pbDistances.data(646 * idx + 32);

            auto pa2pb_x_zzzz = pa2pbDistances.data(646 * idx + 33);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 137);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 138);

            auto pa2pb_xy_yyy = pa2pbDistances.data(646 * idx + 151);

            auto pa2pb_xy_yyz = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_xy_yzz = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_xy_zzz = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_xyy_yy = pa2pbDistances.data(646 * idx + 414);

            auto pa2pb_xyy_yz = pa2pbDistances.data(646 * idx + 415);

            auto pa2pb_xyy_zz = pa2pbDistances.data(646 * idx + 416);

            auto pa2pb_xyy_yyyy = pa2pbDistances.data(646 * idx + 437);

            auto pa2pb_xyy_yyyz = pa2pbDistances.data(646 * idx + 438);

            auto pa2pb_xyy_yyzz = pa2pbDistances.data(646 * idx + 439);

            auto pa2pb_xyy_yzzz = pa2pbDistances.data(646 * idx + 440);

            auto pa2pb_xyy_zzzz = pa2pbDistances.data(646 * idx + 441);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyy_yyyy = primBuffer.data(150 * idx + 55);

            auto t_xyy_yyyz = primBuffer.data(150 * idx + 56);

            auto t_xyy_yyzz = primBuffer.data(150 * idx + 57);

            auto t_xyy_yzzz = primBuffer.data(150 * idx + 58);

            auto t_xyy_zzzz = primBuffer.data(150 * idx + 59);

            // Batch of Integrals (55,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_yy, pa2pb_x_yyyy, pa2pb_x_yyyz, pa2pb_x_yyzz, \
                                     pa2pb_x_yz, pa2pb_x_yzzz, pa2pb_x_zz, pa2pb_x_zzzz, pa2pb_xy_y, pa2pb_xy_yyy, \
                                     pa2pb_xy_yyz, pa2pb_xy_yzz, pa2pb_xy_z, pa2pb_xy_zzz, pa2pb_xyy_yy, pa2pb_xyy_yyyy, \
                                     pa2pb_xyy_yyyz, pa2pb_xyy_yyzz, pa2pb_xyy_yz, pa2pb_xyy_yzzz, pa2pb_xyy_zz, \
                                     pa2pb_xyy_zzzz, pa_x, pa_xyy, r_0_0, s_0_0, t_xyy_yyyy, t_xyy_yyyz, t_xyy_yyzz, t_xyy_yzzz, \
                                     t_xyy_zzzz: VLX_ALIGN)
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

                t_xyy_yyyy[j] = fl_s_0_0 * (1.875 * pa_x[j] * fl3_fx + 0.75 * pa_xyy[j] * fl2_fx + 6.0 * pa2pb_xy_y[j] * fl2_fx + 4.5 * pa2pb_x_yy[j] * fl2_fx + 3.0 * pa2pb_xyy_yy[j] * fl1_fx + 4.0 * pa2pb_xy_yyy[j] * fl1_fx + 0.5 * pa2pb_x_yyyy[j] * fl1_fx + pa2pb_xyy_yyyy[j]);

                t_xyy_yyyy[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx + 60.0 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 45.0 * pa2pb_x_yy[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xyy_yy[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_xyy_yy[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_xy_yyy[j] * fl1_fz * fl1_fx - pa2pb_x_yyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_yyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_yyyy[j] * fl1_fz);

                t_xyy_yyyz[j] = fl_s_0_0 * (1.5 * pa2pb_xy_z[j] * fl2_fx + 2.25 * pa2pb_x_yz[j] * fl2_fx + 1.5 * pa2pb_xyy_yz[j] * fl1_fx + 3.0 * pa2pb_xy_yyz[j] * fl1_fx + 0.5 * pa2pb_x_yyyz[j] * fl1_fx + pa2pb_xyy_yyyz[j]);

                t_xyy_yyyz[j] += fl_r_0_0 * (-3.0 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xy_z[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_x_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyy_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyy_yz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xy_yyz[j] * fl1_fz * fl1_fx - pa2pb_x_yyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_yyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_yyyz[j] * fl1_fz);

                t_xyy_yyzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xyy[j] * fl2_fx + pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_x_zz[j] * fl2_fx + 0.25 * pa2pb_x_yy[j] * fl2_fx + 0.5 * pa2pb_xyy_yy[j] * fl1_fx + 0.5 * pa2pb_xyy_zz[j] * fl1_fx + 2.0 * pa2pb_xy_yzz[j] * fl1_fx + 0.5 * pa2pb_x_yyzz[j] * fl1_fx + pa2pb_xyy_yyzz[j]);

                t_xyy_yyzz[j] += fl_r_0_0 * (-pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xyy[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xyy_yy[j] * fl1_fz * fl1_fgb - pa2pb_xyy_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_yy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xyy_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xyy_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xy_yzz[j] * fl1_fz * fl1_fx - pa2pb_x_yyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_yyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_yyzz[j] * fl1_fz);

                t_xyy_yzzz[j] = fl_s_0_0 * (1.5 * pa2pb_xy_z[j] * fl2_fx + 0.75 * pa2pb_x_yz[j] * fl2_fx + 1.5 * pa2pb_xyy_yz[j] * fl1_fx + pa2pb_xy_zzz[j] * fl1_fx + 0.5 * pa2pb_x_yzzz[j] * fl1_fx + pa2pb_xyy_yzzz[j]);

                t_xyy_yzzz[j] += fl_r_0_0 * (-3.0 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xy_z[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyy_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_yz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyy_yz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_zzz[j] * fl1_fz * fl1_fx - pa2pb_x_yzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_yzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_yzzz[j] * fl1_fz);

                t_xyy_zzzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xyy[j] * fl2_fx + 1.5 * pa2pb_x_zz[j] * fl2_fx + 3.0 * pa2pb_xyy_zz[j] * fl1_fx + 0.5 * pa2pb_x_zzzz[j] * fl1_fx + pa2pb_xyy_zzzz[j]);

                t_xyy_zzzz[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xyy_zz[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_zz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xyy_zz[j] * fl1_fz * fl1_fx - pa2pb_x_zzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_zzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_60_65(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (60,65)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 136);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 137);

            auto pa2pb_xy_xxx = pa2pbDistances.data(646 * idx + 145);

            auto pa2pb_xy_xxy = pa2pbDistances.data(646 * idx + 146);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 171);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 172);

            auto pa2pb_xz_xxx = pa2pbDistances.data(646 * idx + 179);

            auto pa2pb_xz_xxy = pa2pbDistances.data(646 * idx + 180);

            auto pa2pb_xz_xxz = pa2pbDistances.data(646 * idx + 181);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 238);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 239);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 240);

            auto pa2pb_yz_xxx = pa2pbDistances.data(646 * idx + 247);

            auto pa2pb_yz_xxy = pa2pbDistances.data(646 * idx + 248);

            auto pa2pb_yz_xxz = pa2pbDistances.data(646 * idx + 249);

            auto pa2pb_yz_xyy = pa2pbDistances.data(646 * idx + 250);

            auto pa2pb_yz_xyz = pa2pbDistances.data(646 * idx + 251);

            auto pa2pb_xyz_xx = pa2pbDistances.data(646 * idx + 445);

            auto pa2pb_xyz_xy = pa2pbDistances.data(646 * idx + 446);

            auto pa2pb_xyz_xz = pa2pbDistances.data(646 * idx + 447);

            auto pa2pb_xyz_yy = pa2pbDistances.data(646 * idx + 448);

            auto pa2pb_xyz_yz = pa2pbDistances.data(646 * idx + 449);

            auto pa2pb_xyz_xxxx = pa2pbDistances.data(646 * idx + 461);

            auto pa2pb_xyz_xxxy = pa2pbDistances.data(646 * idx + 462);

            auto pa2pb_xyz_xxxz = pa2pbDistances.data(646 * idx + 463);

            auto pa2pb_xyz_xxyy = pa2pbDistances.data(646 * idx + 464);

            auto pa2pb_xyz_xxyz = pa2pbDistances.data(646 * idx + 465);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyz_xxxx = primBuffer.data(150 * idx + 60);

            auto t_xyz_xxxy = primBuffer.data(150 * idx + 61);

            auto t_xyz_xxxz = primBuffer.data(150 * idx + 62);

            auto t_xyz_xxyy = primBuffer.data(150 * idx + 63);

            auto t_xyz_xxyz = primBuffer.data(150 * idx + 64);

            // Batch of Integrals (60,65)

            #pragma omp simd aligned(fgb, fx, fz, pa2pb_x_xx, pa2pb_xy_x, pa2pb_xy_xxx, pa2pb_xy_xxy, \
                                     pa2pb_xy_y, pa2pb_xyz_xx, pa2pb_xyz_xxxx, pa2pb_xyz_xxxy, pa2pb_xyz_xxxz, \
                                     pa2pb_xyz_xxyy, pa2pb_xyz_xxyz, pa2pb_xyz_xy, pa2pb_xyz_xz, pa2pb_xyz_yy, \
                                     pa2pb_xyz_yz, pa2pb_xz_x, pa2pb_xz_xxx, pa2pb_xz_xxy, pa2pb_xz_xxz, pa2pb_xz_y, \
                                     pa2pb_xz_z, pa2pb_y_xx, pa2pb_y_xy, pa2pb_yz_x, pa2pb_yz_xxx, pa2pb_yz_xxy, \
                                     pa2pb_yz_xxz, pa2pb_yz_xyy, pa2pb_yz_xyz, pa2pb_yz_y, pa2pb_yz_z, pa2pb_z_xx, \
                                     pa2pb_z_xy, pa2pb_z_xz, pa_x, pa_xyz, pa_y, pa_z, pb_x, r_0_0, s_0_0, t_xyz_xxxx, \
                                     t_xyz_xxxy, t_xyz_xxxz, t_xyz_xxyy, t_xyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyz_xxxx[j] = fl_s_0_0 * (0.75 * pa_xyz[j] * fl2_fx + 3.0 * pa2pb_yz_x[j] * fl2_fx + 3.0 * pa2pb_xyz_xx[j] * fl1_fx + 2.0 * pa2pb_yz_xxx[j] * fl1_fx + pa2pb_xyz_xxxx[j]);

                t_xyz_xxxx[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xyz[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz - 6.0 * pa2pb_xyz_xx[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_xyz_xx[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yz_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xxxx[j] * fl1_fz);

                t_xyz_xxxy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_yz_y[j] * fl2_fx + 0.75 * pa2pb_z_xx[j] * fl2_fx + 1.5 * pa2pb_xyz_xy[j] * fl1_fx + 0.5 * pa2pb_xz_xxx[j] * fl1_fx + 1.5 * pa2pb_yz_xxy[j] * fl1_fx + pa2pb_xyz_xxxy[j]);

                t_xyz_xxxy[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xyz_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xz_xxx[j] * fl1_fx * fl1_fz + 18.0 * pa2pb_yz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xxxy[j] * fl1_fz);

                t_xyz_xxxz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_xx[j] * fl2_fx + 1.5 * pa2pb_xyz_xz[j] * fl1_fx + 0.5 * pa2pb_xy_xxx[j] * fl1_fx + 1.5 * pa2pb_yz_xxz[j] * fl1_fx + pa2pb_xyz_xxxz[j]);

                t_xyz_xxxz[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xyz_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xy_xxx[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xxxz[j] * fl1_fz);

                t_xyz_xxyy[j] = fl_s_0_0 * (0.25 * pa_xyz[j] * fl2_fx + 0.5 * pa2pb_xz_y[j] * fl2_fx + 0.5 * pa2pb_yz_x[j] * fl2_fx + pa2pb_z_xy[j] * fl2_fx + 0.5 * pa2pb_xyz_xx[j] * fl1_fx + 0.5 * pa2pb_xyz_yy[j] * fl1_fx + pa2pb_xz_xxy[j] * fl1_fx + pa2pb_yz_xyy[j] * fl1_fx + pa2pb_xyz_xxyy[j]);

                t_xyz_xxyy[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xyz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz - pa2pb_xyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyz_yy[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_xyz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xyz_yy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_xxy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yz_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xxyy[j] * fl1_fz);

                t_xyz_xxyz[j] = fl_s_0_0 * (0.125 * pa_x[j] * fl3_fx + 0.25 * pb_x[j] * fl3_fx + 0.25 * pa2pb_xy_y[j] * fl2_fx + 0.25 * pa2pb_xz_z[j] * fl2_fx + 0.25 * pa2pb_x_xx[j] * fl2_fx + 0.5 * pa2pb_y_xy[j] * fl2_fx + 0.5 * pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_xyz_yz[j] * fl1_fx + 0.5 * pa2pb_xy_xxy[j] * fl1_fx + 0.5 * pa2pb_xz_xxz[j] * fl1_fx + pa2pb_yz_xyz[j] * fl1_fx + pa2pb_xyz_xxyz[j]);

                t_xyz_xxyz[j] += fl_r_0_0 * (-0.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + pa_x[j] * fl3_fx * fl1_fz + 2.0 * pb_x[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz - pa2pb_xyz_yz[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_xyz_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xy_xxy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xz_xxz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yz_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_65_70(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (65,70)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 136);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 138);

            auto pa2pb_xy_xxz = pa2pbDistances.data(646 * idx + 147);

            auto pa2pb_xy_xyy = pa2pbDistances.data(646 * idx + 148);

            auto pa2pb_xy_xyz = pa2pbDistances.data(646 * idx + 149);

            auto pa2pb_xy_xzz = pa2pbDistances.data(646 * idx + 150);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_xz_xyy = pa2pbDistances.data(646 * idx + 182);

            auto pa2pb_xz_xyz = pa2pbDistances.data(646 * idx + 183);

            auto pa2pb_xz_xzz = pa2pbDistances.data(646 * idx + 184);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 238);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 239);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 240);

            auto pa2pb_yz_xzz = pa2pbDistances.data(646 * idx + 252);

            auto pa2pb_yz_yyy = pa2pbDistances.data(646 * idx + 253);

            auto pa2pb_yz_yyz = pa2pbDistances.data(646 * idx + 254);

            auto pa2pb_yz_yzz = pa2pbDistances.data(646 * idx + 255);

            auto pa2pb_yz_zzz = pa2pbDistances.data(646 * idx + 256);

            auto pa2pb_xyz_xx = pa2pbDistances.data(646 * idx + 445);

            auto pa2pb_xyz_xy = pa2pbDistances.data(646 * idx + 446);

            auto pa2pb_xyz_xz = pa2pbDistances.data(646 * idx + 447);

            auto pa2pb_xyz_zz = pa2pbDistances.data(646 * idx + 450);

            auto pa2pb_xyz_xxzz = pa2pbDistances.data(646 * idx + 466);

            auto pa2pb_xyz_xyyy = pa2pbDistances.data(646 * idx + 467);

            auto pa2pb_xyz_xyyz = pa2pbDistances.data(646 * idx + 468);

            auto pa2pb_xyz_xyzz = pa2pbDistances.data(646 * idx + 469);

            auto pa2pb_xyz_xzzz = pa2pbDistances.data(646 * idx + 470);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyz_xxzz = primBuffer.data(150 * idx + 65);

            auto t_xyz_xyyy = primBuffer.data(150 * idx + 66);

            auto t_xyz_xyyz = primBuffer.data(150 * idx + 67);

            auto t_xyz_xyzz = primBuffer.data(150 * idx + 68);

            auto t_xyz_xzzz = primBuffer.data(150 * idx + 69);

            // Batch of Integrals (65,70)

            #pragma omp simd aligned(fgb, fx, fz, pa2pb_x_xy, pa2pb_x_xz, pa2pb_xy_x, pa2pb_xy_xxz, \
                                     pa2pb_xy_xyy, pa2pb_xy_xyz, pa2pb_xy_xzz, pa2pb_xy_z, pa2pb_xyz_xx, pa2pb_xyz_xxzz, \
                                     pa2pb_xyz_xy, pa2pb_xyz_xyyy, pa2pb_xyz_xyyz, pa2pb_xyz_xyzz, pa2pb_xyz_xz, \
                                     pa2pb_xyz_xzzz, pa2pb_xyz_zz, pa2pb_xz_x, pa2pb_xz_xyy, pa2pb_xz_xyz, pa2pb_xz_xzz, \
                                     pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yz_x, pa2pb_yz_xzz, \
                                     pa2pb_yz_y, pa2pb_yz_yyy, pa2pb_yz_yyz, pa2pb_yz_yzz, pa2pb_yz_z, pa2pb_yz_zzz, \
                                     pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa_xyz, pa_y, pa_z, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_xyz_xxzz, t_xyz_xyyy, t_xyz_xyyz, t_xyz_xyzz, t_xyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyz_xxzz[j] = fl_s_0_0 * (0.25 * pa_xyz[j] * fl2_fx + 0.5 * pa2pb_xy_z[j] * fl2_fx + 0.5 * pa2pb_yz_x[j] * fl2_fx + pa2pb_y_xz[j] * fl2_fx + 0.5 * pa2pb_xyz_xx[j] * fl1_fx + 0.5 * pa2pb_xyz_zz[j] * fl1_fx + pa2pb_xy_xxz[j] * fl1_fx + pa2pb_yz_xzz[j] * fl1_fx + pa2pb_xyz_xxzz[j]);

                t_xyz_xxzz[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xyz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xy_z[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz - pa2pb_xyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyz_zz[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_xyz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xyz_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xxzz[j] * fl1_fz);

                t_xyz_xyyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_yz_y[j] * fl2_fx + 0.75 * pa2pb_z_yy[j] * fl2_fx + 1.5 * pa2pb_xyz_xy[j] * fl1_fx + 1.5 * pa2pb_xz_xyy[j] * fl1_fx + 0.5 * pa2pb_yz_yyy[j] * fl1_fx + pa2pb_xyz_xyyy[j]);

                t_xyz_xyyy[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xyz_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xz_xyy[j] * fl1_fx * fl1_fz + 6.0 * pa2pb_yz_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xyyy[j] * fl1_fz);

                t_xyz_xyyz[j] = fl_s_0_0 * (0.125 * pa_y[j] * fl3_fx + 0.25 * pb_y[j] * fl3_fx + 0.25 * pa2pb_xy_x[j] * fl2_fx + 0.5 * pa2pb_x_xy[j] * fl2_fx + 0.25 * pa2pb_yz_z[j] * fl2_fx + 0.25 * pa2pb_y_yy[j] * fl2_fx + 0.5 * pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_xyz_xz[j] * fl1_fx + 0.5 * pa2pb_xy_xyy[j] * fl1_fx + pa2pb_xz_xyz[j] * fl1_fx + 0.5 * pa2pb_yz_yyz[j] * fl1_fx + pa2pb_xyz_xyyz[j]);

                t_xyz_xyyz[j] += fl_r_0_0 * (-0.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + pa_y[j] * fl3_fx * fl1_fz + 2.0 * pb_y[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz - pa2pb_xyz_xz[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_xyz_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xy_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_xyz[j] * fl1_fx * fl1_fz + 6.0 * pa2pb_yz_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xyyz[j] * fl1_fz);

                t_xyz_xyzz[j] = fl_s_0_0 * (0.125 * pa_z[j] * fl3_fx + 0.25 * pb_z[j] * fl3_fx + 0.25 * pa2pb_xz_x[j] * fl2_fx + 0.5 * pa2pb_x_xz[j] * fl2_fx + 0.25 * pa2pb_yz_y[j] * fl2_fx + 0.5 * pa2pb_y_yz[j] * fl2_fx + 0.25 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pa2pb_xyz_xy[j] * fl1_fx + pa2pb_xy_xyz[j] * fl1_fx + 0.5 * pa2pb_xz_xzz[j] * fl1_fx + 0.5 * pa2pb_yz_yzz[j] * fl1_fx + pa2pb_xyz_xyzz[j]);

                t_xyz_xyzz[j] += fl_r_0_0 * (-0.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + pa_z[j] * fl3_fx * fl1_fz + 2.0 * pb_z[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz - pa2pb_xyz_xy[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_xyz_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_xyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xz_xzz[j] * fl1_fx * fl1_fz + 6.0 * pa2pb_yz_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xyzz[j] * fl1_fz);

                t_xyz_xzzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_zz[j] * fl2_fx + 1.5 * pa2pb_xyz_xz[j] * fl1_fx + 1.5 * pa2pb_xy_xzz[j] * fl1_fx + 0.5 * pa2pb_yz_zzz[j] * fl1_fx + pa2pb_xyz_xzzz[j]);

                t_xyz_xzzz[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xyz_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xy_xzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yz_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_70_75(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (70,75)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 137);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 138);

            auto pa2pb_xy_yyy = pa2pbDistances.data(646 * idx + 151);

            auto pa2pb_xy_yyz = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_xy_yzz = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_xy_zzz = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 171);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 172);

            auto pa2pb_xz_yyy = pa2pbDistances.data(646 * idx + 185);

            auto pa2pb_xz_yyz = pa2pbDistances.data(646 * idx + 186);

            auto pa2pb_xz_yzz = pa2pbDistances.data(646 * idx + 187);

            auto pa2pb_xz_zzz = pa2pbDistances.data(646 * idx + 188);

            auto pa2pb_xyz_yy = pa2pbDistances.data(646 * idx + 448);

            auto pa2pb_xyz_yz = pa2pbDistances.data(646 * idx + 449);

            auto pa2pb_xyz_zz = pa2pbDistances.data(646 * idx + 450);

            auto pa2pb_xyz_yyyy = pa2pbDistances.data(646 * idx + 471);

            auto pa2pb_xyz_yyyz = pa2pbDistances.data(646 * idx + 472);

            auto pa2pb_xyz_yyzz = pa2pbDistances.data(646 * idx + 473);

            auto pa2pb_xyz_yzzz = pa2pbDistances.data(646 * idx + 474);

            auto pa2pb_xyz_zzzz = pa2pbDistances.data(646 * idx + 475);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyz_yyyy = primBuffer.data(150 * idx + 70);

            auto t_xyz_yyyz = primBuffer.data(150 * idx + 71);

            auto t_xyz_yyzz = primBuffer.data(150 * idx + 72);

            auto t_xyz_yzzz = primBuffer.data(150 * idx + 73);

            auto t_xyz_zzzz = primBuffer.data(150 * idx + 74);

            // Batch of Integrals (70,75)

            #pragma omp simd aligned(fgb, fx, fz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xy_y, pa2pb_xy_yyy, \
                                     pa2pb_xy_yyz, pa2pb_xy_yzz, pa2pb_xy_z, pa2pb_xy_zzz, pa2pb_xyz_yy, pa2pb_xyz_yyyy, \
                                     pa2pb_xyz_yyyz, pa2pb_xyz_yyzz, pa2pb_xyz_yz, pa2pb_xyz_yzzz, pa2pb_xyz_zz, \
                                     pa2pb_xyz_zzzz, pa2pb_xz_y, pa2pb_xz_yyy, pa2pb_xz_yyz, pa2pb_xz_yzz, pa2pb_xz_z, \
                                     pa2pb_xz_zzz, pa_x, pa_xyz, r_0_0, s_0_0, t_xyz_yyyy, t_xyz_yyyz, t_xyz_yyzz, t_xyz_yzzz, \
                                     t_xyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_r_0_0 = r_0_0[j];

                double fl_s_0_0 = s_0_0[j];

                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                t_xyz_yyyy[j] = fl_s_0_0 * (0.75 * pa_xyz[j] * fl2_fx + 3.0 * pa2pb_xz_y[j] * fl2_fx + 3.0 * pa2pb_xyz_yy[j] * fl1_fx + 2.0 * pa2pb_xz_yyy[j] * fl1_fx + pa2pb_xyz_yyyy[j]);

                t_xyz_yyyy[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xyz[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz - 6.0 * pa2pb_xyz_yy[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_xyz_yy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xz_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_yyyy[j] * fl1_fz);

                t_xyz_yyyz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_yy[j] * fl2_fx + 1.5 * pa2pb_xyz_yz[j] * fl1_fx + 0.5 * pa2pb_xy_yyy[j] * fl1_fx + 1.5 * pa2pb_xz_yyz[j] * fl1_fx + pa2pb_xyz_yyyz[j]);

                t_xyz_yyyz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xyz_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xy_yyy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xz_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_yyyz[j] * fl1_fz);

                t_xyz_yyzz[j] = fl_s_0_0 * (0.25 * pa_xyz[j] * fl2_fx + 0.5 * pa2pb_xy_z[j] * fl2_fx + 0.5 * pa2pb_xz_y[j] * fl2_fx + pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_xyz_yy[j] * fl1_fx + 0.5 * pa2pb_xyz_zz[j] * fl1_fx + pa2pb_xy_yyz[j] * fl1_fx + pa2pb_xz_yzz[j] * fl1_fx + pa2pb_xyz_yyzz[j]);

                t_xyz_yyzz[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xyz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xy_z[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz - pa2pb_xyz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xyz_zz[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_xyz_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xyz_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xy_yyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_yyzz[j] * fl1_fz);

                t_xyz_yzzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_zz[j] * fl2_fx + 1.5 * pa2pb_xyz_yz[j] * fl1_fx + 1.5 * pa2pb_xy_yzz[j] * fl1_fx + 0.5 * pa2pb_xz_zzz[j] * fl1_fx + pa2pb_xyz_yzzz[j]);

                t_xyz_yzzz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xyz_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xy_yzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xz_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_yzzz[j] * fl1_fz);

                t_xyz_zzzz[j] = fl_s_0_0 * (0.75 * pa_xyz[j] * fl2_fx + 3.0 * pa2pb_xy_z[j] * fl2_fx + 3.0 * pa2pb_xyz_zz[j] * fl1_fx + 2.0 * pa2pb_xy_zzz[j] * fl1_fx + pa2pb_xyz_zzzz[j]);

                t_xyz_zzzz[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xyz[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_xy_z[j] * fl1_fz * fl2_fx - 6.0 * pa2pb_xyz_zz[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_xyz_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xy_zzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_75_80(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (75,80)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_xxxx = pa2pbDistances.data(646 * idx + 19);

            auto pa2pb_x_xxxy = pa2pbDistances.data(646 * idx + 20);

            auto pa2pb_x_xxxz = pa2pbDistances.data(646 * idx + 21);

            auto pa2pb_x_xxyy = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_x_xxyz = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 171);

            auto pa2pb_xz_xxx = pa2pbDistances.data(646 * idx + 179);

            auto pa2pb_xz_xxy = pa2pbDistances.data(646 * idx + 180);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_zz_xxx = pa2pbDistances.data(646 * idx + 281);

            auto pa2pb_zz_xxy = pa2pbDistances.data(646 * idx + 282);

            auto pa2pb_zz_xxz = pa2pbDistances.data(646 * idx + 283);

            auto pa2pb_zz_xyy = pa2pbDistances.data(646 * idx + 284);

            auto pa2pb_zz_xyz = pa2pbDistances.data(646 * idx + 285);

            auto pa2pb_xzz_xx = pa2pbDistances.data(646 * idx + 479);

            auto pa2pb_xzz_xy = pa2pbDistances.data(646 * idx + 480);

            auto pa2pb_xzz_xz = pa2pbDistances.data(646 * idx + 481);

            auto pa2pb_xzz_yy = pa2pbDistances.data(646 * idx + 482);

            auto pa2pb_xzz_yz = pa2pbDistances.data(646 * idx + 483);

            auto pa2pb_xzz_xxxx = pa2pbDistances.data(646 * idx + 495);

            auto pa2pb_xzz_xxxy = pa2pbDistances.data(646 * idx + 496);

            auto pa2pb_xzz_xxxz = pa2pbDistances.data(646 * idx + 497);

            auto pa2pb_xzz_xxyy = pa2pbDistances.data(646 * idx + 498);

            auto pa2pb_xzz_xxyz = pa2pbDistances.data(646 * idx + 499);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzz_xxxx = primBuffer.data(150 * idx + 75);

            auto t_xzz_xxxy = primBuffer.data(150 * idx + 76);

            auto t_xzz_xxxz = primBuffer.data(150 * idx + 77);

            auto t_xzz_xxyy = primBuffer.data(150 * idx + 78);

            auto t_xzz_xxyz = primBuffer.data(150 * idx + 79);

            // Batch of Integrals (75,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xxxx, pa2pb_x_xxxy, pa2pb_x_xxxz, \
                                     pa2pb_x_xxyy, pa2pb_x_xxyz, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_xz_x, pa2pb_xz_xxx, pa2pb_xz_xxy, pa2pb_xz_y, pa2pb_xzz_xx, pa2pb_xzz_xxxx, \
                                     pa2pb_xzz_xxxy, pa2pb_xzz_xxxz, pa2pb_xzz_xxyy, pa2pb_xzz_xxyz, pa2pb_xzz_xy, \
                                     pa2pb_xzz_xz, pa2pb_xzz_yy, pa2pb_xzz_yz, pa2pb_z_xx, pa2pb_z_xy, pa2pb_zz_x, \
                                     pa2pb_zz_xxx, pa2pb_zz_xxy, pa2pb_zz_xxz, pa2pb_zz_xyy, pa2pb_zz_xyz, pa2pb_zz_y, \
                                     pa2pb_zz_z, pa_x, pa_xzz, pa_z, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_y, pb_z, \
                                     r_0_0, s_0_0, t_xzz_xxxx, t_xzz_xxxy, t_xzz_xxxz, t_xzz_xxyy, t_xzz_xxyz: VLX_ALIGN)
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

                t_xzz_xxxx[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 1.5 * pb_x[j] * fl3_fx + 0.75 * pa_xzz[j] * fl2_fx + 3.0 * pa2pb_zz_x[j] * fl2_fx + 1.5 * pa2pb_x_xx[j] * fl2_fx + pb_xxx[j] * fl2_fx + 3.0 * pa2pb_xzz_xx[j] * fl1_fx + 2.0 * pa2pb_zz_xxx[j] * fl1_fx + 0.5 * pa2pb_x_xxxx[j] * fl1_fx + pa2pb_xzz_xxxx[j]);

                t_xzz_xxxx[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx - 6.0 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xzz[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xzz_xx[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_xx[j] * fl1_fz * fl2_fx + 10.0 * pb_xxx[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xzz_xx[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_zz_xxx[j] * fl1_fx * fl1_fz - pa2pb_x_xxxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xxxx[j] * fl1_fz);

                t_xzz_xxxy[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_zz_y[j] * fl2_fx + 0.75 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pb_xxy[j] * fl2_fx + 1.5 * pa2pb_xzz_xy[j] * fl1_fx + 1.5 * pa2pb_zz_xxy[j] * fl1_fx + 0.5 * pa2pb_x_xxxy[j] * fl1_fx + pa2pb_xzz_xxxy[j]);

                t_xzz_xxxy[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzz_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xy[j] * fl1_fz * fl2_fx + 7.5 * pb_xxy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xzz_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zz_xxy[j] * fl1_fx * fl1_fz - pa2pb_x_xxxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xxxy[j] * fl1_fz);

                t_xzz_xxxz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl3_fx + 0.375 * pb_z[j] * fl3_fx + 1.5 * pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_zz_z[j] * fl2_fx + 1.5 * pa2pb_z_xx[j] * fl2_fx + 0.75 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pb_xxz[j] * fl2_fx + 1.5 * pa2pb_xzz_xz[j] * fl1_fx + pa2pb_xz_xxx[j] * fl1_fx + 1.5 * pa2pb_zz_xxz[j] * fl1_fx + 0.5 * pa2pb_x_xxxz[j] * fl1_fx + pa2pb_xzz_xxxz[j]);

                t_xzz_xxxz[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_xz_x[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzz_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xz[j] * fl1_fz * fl2_fx + 7.5 * pb_xxz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xzz_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_xxx[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zz_xxz[j] * fl1_fx * fl1_fz - pa2pb_x_xxxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xxxz[j] * fl1_fz);

                t_xzz_xxyy[j] = fl_s_0_0 * (0.125 * pa_x[j] * fl3_fx + 0.25 * pb_x[j] * fl3_fx + 0.25 * pa_xzz[j] * fl2_fx + 0.5 * pa2pb_zz_x[j] * fl2_fx + 0.25 * pa2pb_x_xx[j] * fl2_fx + 0.25 * pa2pb_x_yy[j] * fl2_fx + 0.5 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_xzz_xx[j] * fl1_fx + 0.5 * pa2pb_xzz_yy[j] * fl1_fx + pa2pb_zz_xyy[j] * fl1_fx + 0.5 * pa2pb_x_xxyy[j] * fl1_fx + pa2pb_xzz_xxyy[j]);

                t_xzz_xxyy[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + pa_x[j] * fl1_fz * fl3_fx - pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xzz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xzz_yy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_xx[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_x_yy[j] * fl1_fz * fl2_fx + 5.0 * pb_xyy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xzz_yy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zz_xyy[j] * fl1_fx * fl1_fz - pa2pb_x_xxyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xxyy[j] * fl1_fz);

                t_xzz_xxyz[j] = fl_s_0_0 * (0.5 * pa2pb_xz_y[j] * fl2_fx + pa2pb_z_xy[j] * fl2_fx + 0.25 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_xzz_yz[j] * fl1_fx + pa2pb_xz_xxy[j] * fl1_fx + pa2pb_zz_xyz[j] * fl1_fx + 0.5 * pa2pb_x_xxyz[j] * fl1_fx + pa2pb_xzz_xxyz[j]);

                t_xzz_xxyz[j] += fl_r_0_0 * (-pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_xz_y[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_yz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_yz[j] * fl1_fz * fl2_fx + 5.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzz_yz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zz_xyz[j] * fl1_fx * fl1_fz - pa2pb_x_xxyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_80_85(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (80,85)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_x_xxzz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_x_xyyy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_x_xyyz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_x_xyzz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_x_xzzz = pa2pbDistances.data(646 * idx + 28);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 172);

            auto pa2pb_xz_xxz = pa2pbDistances.data(646 * idx + 181);

            auto pa2pb_xz_xyy = pa2pbDistances.data(646 * idx + 182);

            auto pa2pb_xz_xyz = pa2pbDistances.data(646 * idx + 183);

            auto pa2pb_xz_xzz = pa2pbDistances.data(646 * idx + 184);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_zz_xzz = pa2pbDistances.data(646 * idx + 286);

            auto pa2pb_zz_yyy = pa2pbDistances.data(646 * idx + 287);

            auto pa2pb_zz_yyz = pa2pbDistances.data(646 * idx + 288);

            auto pa2pb_zz_yzz = pa2pbDistances.data(646 * idx + 289);

            auto pa2pb_zz_zzz = pa2pbDistances.data(646 * idx + 290);

            auto pa2pb_xzz_xx = pa2pbDistances.data(646 * idx + 479);

            auto pa2pb_xzz_xy = pa2pbDistances.data(646 * idx + 480);

            auto pa2pb_xzz_xz = pa2pbDistances.data(646 * idx + 481);

            auto pa2pb_xzz_zz = pa2pbDistances.data(646 * idx + 484);

            auto pa2pb_xzz_xxzz = pa2pbDistances.data(646 * idx + 500);

            auto pa2pb_xzz_xyyy = pa2pbDistances.data(646 * idx + 501);

            auto pa2pb_xzz_xyyz = pa2pbDistances.data(646 * idx + 502);

            auto pa2pb_xzz_xyzz = pa2pbDistances.data(646 * idx + 503);

            auto pa2pb_xzz_xzzz = pa2pbDistances.data(646 * idx + 504);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzz_xxzz = primBuffer.data(150 * idx + 80);

            auto t_xzz_xyyy = primBuffer.data(150 * idx + 81);

            auto t_xzz_xyyz = primBuffer.data(150 * idx + 82);

            auto t_xzz_xyzz = primBuffer.data(150 * idx + 83);

            auto t_xzz_xzzz = primBuffer.data(150 * idx + 84);

            // Batch of Integrals (80,85)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xxzz, pa2pb_x_xy, pa2pb_x_xyyy, \
                                     pa2pb_x_xyyz, pa2pb_x_xyzz, pa2pb_x_xz, pa2pb_x_xzzz, pa2pb_x_zz, pa2pb_xz_x, \
                                     pa2pb_xz_xxz, pa2pb_xz_xyy, pa2pb_xz_xyz, pa2pb_xz_xzz, pa2pb_xz_z, pa2pb_xzz_xx, \
                                     pa2pb_xzz_xxzz, pa2pb_xzz_xy, pa2pb_xzz_xyyy, pa2pb_xzz_xyyz, pa2pb_xzz_xyzz, \
                                     pa2pb_xzz_xz, pa2pb_xzz_xzzz, pa2pb_xzz_zz, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, \
                                     pa2pb_z_zz, pa2pb_zz_x, pa2pb_zz_xzz, pa2pb_zz_y, pa2pb_zz_yyy, pa2pb_zz_yyz, \
                                     pa2pb_zz_yzz, pa2pb_zz_z, pa2pb_zz_zzz, pa_x, pa_xzz, pa_z, pb_x, pb_xzz, pb_y, pb_yyy, \
                                     pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, s_0_0, t_xzz_xxzz, t_xzz_xyyy, t_xzz_xyyz, \
                                     t_xzz_xyzz, t_xzz_xzzz: VLX_ALIGN)
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

                t_xzz_xxzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pb_x[j] * fl3_fx + 0.25 * pa_xzz[j] * fl2_fx + pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_xx[j] * fl2_fx + 0.5 * pa2pb_zz_x[j] * fl2_fx + 2.0 * pa2pb_z_xz[j] * fl2_fx + 0.25 * pa2pb_x_zz[j] * fl2_fx + 0.5 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_xzz_xx[j] * fl1_fx + 0.5 * pa2pb_xzz_zz[j] * fl1_fx + 2.0 * pa2pb_xz_xxz[j] * fl1_fx + pa2pb_zz_xzz[j] * fl1_fx + 0.5 * pa2pb_x_xxzz[j] * fl1_fx + pa2pb_xzz_xxzz[j]);

                t_xzz_xxzz[j] += fl_r_0_0 * (-pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 6.0 * pb_x[j] * fl3_fx * fl1_fz - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xzz[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xz_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_zz_x[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xzz_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_zz[j] * fl1_fz * fl2_fx + 5.0 * pb_xzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xzz_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xz_xxz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zz_xzz[j] * fl1_fx * fl1_fz - pa2pb_x_xxzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xxzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xxzz[j] * fl1_fz);

                t_xzz_xyyy[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_zz_y[j] * fl2_fx + 0.75 * pa2pb_x_xy[j] * fl2_fx + 0.25 * pb_yyy[j] * fl2_fx + 1.5 * pa2pb_xzz_xy[j] * fl1_fx + 0.5 * pa2pb_zz_yyy[j] * fl1_fx + 0.5 * pa2pb_x_xyyy[j] * fl1_fx + pa2pb_xzz_xyyy[j]);

                t_xzz_xyyy[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzz_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xy[j] * fl1_fz * fl2_fx + 2.5 * pb_yyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xzz_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yyy[j] * fl1_fx * fl1_fz - pa2pb_x_xyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xyyy[j] * fl1_fz);

                t_xzz_xyyz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl3_fx + 0.125 * pb_z[j] * fl3_fx + 0.5 * pa2pb_xz_x[j] * fl2_fx + 0.25 * pa2pb_zz_z[j] * fl2_fx + 0.5 * pa2pb_z_yy[j] * fl2_fx + 0.25 * pa2pb_x_xz[j] * fl2_fx + 0.25 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_xzz_xz[j] * fl1_fx + pa2pb_xz_xyy[j] * fl1_fx + 0.5 * pa2pb_zz_yyz[j] * fl1_fx + 0.5 * pa2pb_x_xyyz[j] * fl1_fx + pa2pb_xzz_xyyz[j]);

                t_xzz_xyyz[j] += fl_r_0_0 * (-0.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_z[j] * fl3_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + pb_z[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_xz_x[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_yy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_xz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_xz[j] * fl1_fz * fl2_fx + 2.5 * pb_yyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzz_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_xyy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yyz[j] * fl1_fx * fl1_fz - pa2pb_x_xyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xyyz[j] * fl1_fz);

                t_xzz_xyzz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_x_xy[j] * fl2_fx + 0.25 * pa2pb_zz_y[j] * fl2_fx + pa2pb_z_yz[j] * fl2_fx + 0.25 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_xzz_xy[j] * fl1_fx + 2.0 * pa2pb_xz_xyz[j] * fl1_fx + 0.5 * pa2pb_zz_yzz[j] * fl1_fx + 0.5 * pa2pb_x_xyzz[j] * fl1_fx + pa2pb_xzz_xyzz[j]);

                t_xzz_xyzz[j] += fl_r_0_0 * (3.0 * pb_y[j] * fl3_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzz_xy[j] * fl1_fz * fl1_fgb + 2.5 * pb_yzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzz_xy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xz_xyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yzz[j] * fl1_fx * fl1_fz - pa2pb_x_xyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xyzz[j] * fl1_fz);

                t_xzz_xzzz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl3_fx + 1.125 * pb_z[j] * fl3_fx + 1.5 * pa2pb_xz_x[j] * fl2_fx + 2.25 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pa2pb_zz_z[j] * fl2_fx + 1.5 * pa2pb_z_zz[j] * fl2_fx + 0.25 * pb_zzz[j] * fl2_fx + 1.5 * pa2pb_xzz_xz[j] * fl1_fx + 3.0 * pa2pb_xz_xzz[j] * fl1_fx + 0.5 * pa2pb_zz_zzz[j] * fl1_fx + 0.5 * pa2pb_x_xzzz[j] * fl1_fx + pa2pb_xzz_xzzz[j]);

                t_xzz_xzzz[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_z[j] * fl3_fx * fl1_fz + 9.0 * pb_z[j] * fl3_fx * fl1_fz - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xz_x[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzz_xz[j] * fl1_fz * fl1_fgb + 2.5 * pb_zzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xzz_xz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xz_xzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_zzz[j] * fl1_fx * fl1_fz - pa2pb_x_xzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_xzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_85_90(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (85,90)

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

            auto pa_x = paDistances.data(19 * idx);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_x_yyyy = pa2pbDistances.data(646 * idx + 29);

            auto pa2pb_x_yyyz = pa2pbDistances.data(646 * idx + 30);

            auto pa2pb_x_yyzz = pa2pbDistances.data(646 * idx + 31);

            auto pa2pb_x_yzzz = pa2pbDistances.data(646 * idx + 32);

            auto pa2pb_x_zzzz = pa2pbDistances.data(646 * idx + 33);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 171);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 172);

            auto pa2pb_xz_yyy = pa2pbDistances.data(646 * idx + 185);

            auto pa2pb_xz_yyz = pa2pbDistances.data(646 * idx + 186);

            auto pa2pb_xz_yzz = pa2pbDistances.data(646 * idx + 187);

            auto pa2pb_xz_zzz = pa2pbDistances.data(646 * idx + 188);

            auto pa2pb_xzz_yy = pa2pbDistances.data(646 * idx + 482);

            auto pa2pb_xzz_yz = pa2pbDistances.data(646 * idx + 483);

            auto pa2pb_xzz_zz = pa2pbDistances.data(646 * idx + 484);

            auto pa2pb_xzz_yyyy = pa2pbDistances.data(646 * idx + 505);

            auto pa2pb_xzz_yyyz = pa2pbDistances.data(646 * idx + 506);

            auto pa2pb_xzz_yyzz = pa2pbDistances.data(646 * idx + 507);

            auto pa2pb_xzz_yzzz = pa2pbDistances.data(646 * idx + 508);

            auto pa2pb_xzz_zzzz = pa2pbDistances.data(646 * idx + 509);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzz_yyyy = primBuffer.data(150 * idx + 85);

            auto t_xzz_yyyz = primBuffer.data(150 * idx + 86);

            auto t_xzz_yyzz = primBuffer.data(150 * idx + 87);

            auto t_xzz_yzzz = primBuffer.data(150 * idx + 88);

            auto t_xzz_zzzz = primBuffer.data(150 * idx + 89);

            // Batch of Integrals (85,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_yy, pa2pb_x_yyyy, pa2pb_x_yyyz, pa2pb_x_yyzz, \
                                     pa2pb_x_yz, pa2pb_x_yzzz, pa2pb_x_zz, pa2pb_x_zzzz, pa2pb_xz_y, pa2pb_xz_yyy, \
                                     pa2pb_xz_yyz, pa2pb_xz_yzz, pa2pb_xz_z, pa2pb_xz_zzz, pa2pb_xzz_yy, pa2pb_xzz_yyyy, \
                                     pa2pb_xzz_yyyz, pa2pb_xzz_yyzz, pa2pb_xzz_yz, pa2pb_xzz_yzzz, pa2pb_xzz_zz, \
                                     pa2pb_xzz_zzzz, pa_x, pa_xzz, r_0_0, s_0_0, t_xzz_yyyy, t_xzz_yyyz, t_xzz_yyzz, t_xzz_yzzz, \
                                     t_xzz_zzzz: VLX_ALIGN)
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

                t_xzz_yyyy[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xzz[j] * fl2_fx + 1.5 * pa2pb_x_yy[j] * fl2_fx + 3.0 * pa2pb_xzz_yy[j] * fl1_fx + 0.5 * pa2pb_x_yyyy[j] * fl1_fx + pa2pb_xzz_yyyy[j]);

                t_xzz_yyyy[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl1_fz * fl3_fx + 7.5 * pa_xzz[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xzz_yy[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_yy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xzz_yy[j] * fl1_fz * fl1_fx - pa2pb_x_yyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_yyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_yyyy[j] * fl1_fz);

                t_xzz_yyyz[j] = fl_s_0_0 * (1.5 * pa2pb_xz_y[j] * fl2_fx + 0.75 * pa2pb_x_yz[j] * fl2_fx + 1.5 * pa2pb_xzz_yz[j] * fl1_fx + pa2pb_xz_yyy[j] * fl1_fx + 0.5 * pa2pb_x_yyyz[j] * fl1_fx + pa2pb_xzz_yyyz[j]);

                t_xzz_yyyz[j] += fl_r_0_0 * (-3.0 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xz_y[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xzz_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_yz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xzz_yz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xz_yyy[j] * fl1_fz * fl1_fx - pa2pb_x_yyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_yyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_yyyz[j] * fl1_fz);

                t_xzz_yyzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xzz[j] * fl2_fx + pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_yy[j] * fl2_fx + 0.25 * pa2pb_x_zz[j] * fl2_fx + 0.5 * pa2pb_xzz_yy[j] * fl1_fx + 0.5 * pa2pb_xzz_zz[j] * fl1_fx + 2.0 * pa2pb_xz_yyz[j] * fl1_fx + 0.5 * pa2pb_x_yyzz[j] * fl1_fx + pa2pb_xzz_yyzz[j]);

                t_xzz_yyzz[j] += fl_r_0_0 * (-pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xzz[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xz_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xzz_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_zz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xzz_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xzz_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xz_yyz[j] * fl1_fz * fl1_fx - pa2pb_x_yyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_yyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_yyzz[j] * fl1_fz);

                t_xzz_yzzz[j] = fl_s_0_0 * (1.5 * pa2pb_xz_y[j] * fl2_fx + 2.25 * pa2pb_x_yz[j] * fl2_fx + 1.5 * pa2pb_xzz_yz[j] * fl1_fx + 3.0 * pa2pb_xz_yzz[j] * fl1_fx + 0.5 * pa2pb_x_yzzz[j] * fl1_fx + pa2pb_xzz_yzzz[j]);

                t_xzz_yzzz[j] += fl_r_0_0 * (-3.0 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_xz_y[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_x_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xzz_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xzz_yz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xz_yzz[j] * fl1_fz * fl1_fx - pa2pb_x_yzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_yzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_yzzz[j] * fl1_fz);

                t_xzz_zzzz[j] = fl_s_0_0 * (1.875 * pa_x[j] * fl3_fx + 0.75 * pa_xzz[j] * fl2_fx + 6.0 * pa2pb_xz_z[j] * fl2_fx + 4.5 * pa2pb_x_zz[j] * fl2_fx + 3.0 * pa2pb_xzz_zz[j] * fl1_fx + 4.0 * pa2pb_xz_zzz[j] * fl1_fx + 0.5 * pa2pb_x_zzzz[j] * fl1_fx + pa2pb_xzz_zzzz[j]);

                t_xzz_zzzz[j] += fl_r_0_0 * (-4.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xzz[j] * fl1_fz * fl2_fx + 60.0 * pa2pb_xz_z[j] * fl1_fz * fl2_fx + 45.0 * pa2pb_x_zz[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xzz_zz[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_xzz_zz[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_xz_zzz[j] * fl1_fz * fl1_fx - pa2pb_x_zzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_x_zzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_90_95(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
                                 const CMemBlock2D<double>& pa2pbDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (90,95)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_xxxx = pa2pbDistances.data(646 * idx + 53);

            auto pa2pb_y_xxxy = pa2pbDistances.data(646 * idx + 54);

            auto pa2pb_y_xxxz = pa2pbDistances.data(646 * idx + 55);

            auto pa2pb_y_xxyy = pa2pbDistances.data(646 * idx + 56);

            auto pa2pb_y_xxyz = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 204);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 205);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 206);

            auto pa2pb_yy_xxx = pa2pbDistances.data(646 * idx + 213);

            auto pa2pb_yy_xxy = pa2pbDistances.data(646 * idx + 214);

            auto pa2pb_yy_xxz = pa2pbDistances.data(646 * idx + 215);

            auto pa2pb_yyy_xx = pa2pbDistances.data(646 * idx + 513);

            auto pa2pb_yyy_xy = pa2pbDistances.data(646 * idx + 514);

            auto pa2pb_yyy_xz = pa2pbDistances.data(646 * idx + 515);

            auto pa2pb_yyy_yy = pa2pbDistances.data(646 * idx + 516);

            auto pa2pb_yyy_yz = pa2pbDistances.data(646 * idx + 517);

            auto pa2pb_yyy_xxxx = pa2pbDistances.data(646 * idx + 529);

            auto pa2pb_yyy_xxxy = pa2pbDistances.data(646 * idx + 530);

            auto pa2pb_yyy_xxxz = pa2pbDistances.data(646 * idx + 531);

            auto pa2pb_yyy_xxyy = pa2pbDistances.data(646 * idx + 532);

            auto pa2pb_yyy_xxyz = pa2pbDistances.data(646 * idx + 533);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyy_xxxx = primBuffer.data(150 * idx + 90);

            auto t_yyy_xxxy = primBuffer.data(150 * idx + 91);

            auto t_yyy_xxxz = primBuffer.data(150 * idx + 92);

            auto t_yyy_xxyy = primBuffer.data(150 * idx + 93);

            auto t_yyy_xxyz = primBuffer.data(150 * idx + 94);

            // Batch of Integrals (90,95)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_y_xxxx, pa2pb_y_xxxy, pa2pb_y_xxxz, \
                                     pa2pb_y_xxyy, pa2pb_y_xxyz, pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, \
                                     pa2pb_yy_x, pa2pb_yy_xxx, pa2pb_yy_xxy, pa2pb_yy_xxz, pa2pb_yy_y, pa2pb_yy_z, \
                                     pa2pb_yyy_xx, pa2pb_yyy_xxxx, pa2pb_yyy_xxxy, pa2pb_yyy_xxxz, pa2pb_yyy_xxyy, \
                                     pa2pb_yyy_xxyz, pa2pb_yyy_xy, pa2pb_yyy_xz, pa2pb_yyy_yy, pa2pb_yyy_yz, pa_y, pa_yyy, \
                                     pb_x, pb_xxx, pb_xxy, pb_xxz, pb_y, pb_z, r_0_0, s_0_0, t_yyy_xxxx, t_yyy_xxxy, \
                                     t_yyy_xxxz, t_yyy_xxyy, t_yyy_xxyz: VLX_ALIGN)
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

                t_yyy_xxxx[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pa_yyy[j] * fl2_fx + 4.5 * pa2pb_y_xx[j] * fl2_fx + 3.0 * pa2pb_yyy_xx[j] * fl1_fx + 1.5 * pa2pb_y_xxxx[j] * fl1_fx + pa2pb_yyy_xxxx[j]);

                t_yyy_xxxx[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl1_fz * fl3_fx + 7.5 * pa_yyy[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyy_xx[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_y_xx[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_yyy_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xxxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xxxx[j] * fl1_fz);

                t_yyy_xxxy[j] = fl_s_0_0 * (1.125 * pb_x[j] * fl3_fx + 2.25 * pa2pb_yy_x[j] * fl2_fx + 2.25 * pa2pb_y_xy[j] * fl2_fx + 0.75 * pb_xxx[j] * fl2_fx + 1.5 * pa2pb_yyy_xy[j] * fl1_fx + 1.5 * pa2pb_yy_xxx[j] * fl1_fx + 1.5 * pa2pb_y_xxxy[j] * fl1_fx + pa2pb_yyy_xxxy[j]);

                t_yyy_xxxy[j] += fl_r_0_0 * (-2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pb_x[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_yy_x[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyy_xy[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_y_xy[j] * fl1_fz * fl2_fx + 7.5 * pb_xxx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyy_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yy_xxx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xxxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xxxy[j] * fl1_fz);

                t_yyy_xxxz[j] = fl_s_0_0 * (2.25 * pa2pb_y_xz[j] * fl2_fx + 1.5 * pa2pb_yyy_xz[j] * fl1_fx + 1.5 * pa2pb_y_xxxz[j] * fl1_fx + pa2pb_yyy_xxxz[j]);

                t_yyy_xxxz[j] += fl_r_0_0 * (-4.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyy_xz[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_y_xz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yyy_xz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xxxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xxxz[j] * fl1_fz);

                t_yyy_xxyy[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pb_y[j] * fl3_fx + 0.25 * pa_yyy[j] * fl2_fx + 1.5 * pa2pb_yy_y[j] * fl2_fx + 2.25 * pa2pb_y_xx[j] * fl2_fx + 0.75 * pa2pb_y_yy[j] * fl2_fx + 1.5 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_yyy_xx[j] * fl1_fx + 0.5 * pa2pb_yyy_yy[j] * fl1_fx + 3.0 * pa2pb_yy_xxy[j] * fl1_fx + 1.5 * pa2pb_y_xxyy[j] * fl1_fx + pa2pb_yyy_xxyy[j]);

                t_yyy_xxyy[j] += fl_r_0_0 * (-3.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_y[j] * fl3_fx * fl1_fz + 2.5 * pa_yyy[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyy_yy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_yy[j] * fl1_fz * fl2_fx + 15.0 * pb_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyy_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyy_yy[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yy_xxy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xxyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xxyy[j] * fl1_fz);

                t_yyy_xxyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_yy_z[j] * fl2_fx + 0.75 * pa2pb_y_yz[j] * fl2_fx + 0.75 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_yyy_yz[j] * fl1_fx + 1.5 * pa2pb_yy_xxz[j] * fl1_fx + 1.5 * pa2pb_y_xxyz[j] * fl1_fx + pa2pb_yyy_xxyz[j]);

                t_yyy_xxyz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_z[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_yz[j] * fl1_fz * fl2_fx + 7.5 * pb_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyy_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yy_xxz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xxyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_95_100(      CMemBlock2D<double>& primBuffer,
                                  const CMemBlock2D<double>& auxBuffer,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pbDistances,
                                  const CMemBlock2D<double>& pa2pbDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (95,100)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_y_xxzz = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_y_xyyy = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_y_xyyz = pa2pbDistances.data(646 * idx + 60);

            auto pa2pb_y_xyzz = pa2pbDistances.data(646 * idx + 61);

            auto pa2pb_y_xzzz = pa2pbDistances.data(646 * idx + 62);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 204);

            auto pa2pb_yy_xyy = pa2pbDistances.data(646 * idx + 216);

            auto pa2pb_yy_xyz = pa2pbDistances.data(646 * idx + 217);

            auto pa2pb_yy_xzz = pa2pbDistances.data(646 * idx + 218);

            auto pa2pb_yyy_xx = pa2pbDistances.data(646 * idx + 513);

            auto pa2pb_yyy_xy = pa2pbDistances.data(646 * idx + 514);

            auto pa2pb_yyy_xz = pa2pbDistances.data(646 * idx + 515);

            auto pa2pb_yyy_zz = pa2pbDistances.data(646 * idx + 518);

            auto pa2pb_yyy_xxzz = pa2pbDistances.data(646 * idx + 534);

            auto pa2pb_yyy_xyyy = pa2pbDistances.data(646 * idx + 535);

            auto pa2pb_yyy_xyyz = pa2pbDistances.data(646 * idx + 536);

            auto pa2pb_yyy_xyzz = pa2pbDistances.data(646 * idx + 537);

            auto pa2pb_yyy_xzzz = pa2pbDistances.data(646 * idx + 538);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyy_xxzz = primBuffer.data(150 * idx + 95);

            auto t_yyy_xyyy = primBuffer.data(150 * idx + 96);

            auto t_yyy_xyyz = primBuffer.data(150 * idx + 97);

            auto t_yyy_xyzz = primBuffer.data(150 * idx + 98);

            auto t_yyy_xzzz = primBuffer.data(150 * idx + 99);

            // Batch of Integrals (95,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_y_xxzz, pa2pb_y_xy, pa2pb_y_xyyy, \
                                     pa2pb_y_xyyz, pa2pb_y_xyzz, pa2pb_y_xz, pa2pb_y_xzzz, pa2pb_y_zz, pa2pb_yy_x, \
                                     pa2pb_yy_xyy, pa2pb_yy_xyz, pa2pb_yy_xzz, pa2pb_yyy_xx, pa2pb_yyy_xxzz, \
                                     pa2pb_yyy_xy, pa2pb_yyy_xyyy, pa2pb_yyy_xyyz, pa2pb_yyy_xyzz, pa2pb_yyy_xz, \
                                     pa2pb_yyy_xzzz, pa2pb_yyy_zz, pa_y, pa_yyy, pb_x, pb_xyy, pb_xyz, pb_xzz, r_0_0, s_0_0, \
                                     t_yyy_xxzz, t_yyy_xyyy, t_yyy_xyyz, t_yyy_xyzz, t_yyy_xzzz: VLX_ALIGN)
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

                t_yyy_xxzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.25 * pa_yyy[j] * fl2_fx + 0.75 * pa2pb_y_xx[j] * fl2_fx + 0.75 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pa2pb_yyy_xx[j] * fl1_fx + 0.5 * pa2pb_yyy_zz[j] * fl1_fx + 1.5 * pa2pb_y_xxzz[j] * fl1_fx + pa2pb_yyy_xxzz[j]);

                t_yyy_xxzz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl1_fz * fl3_fx + 2.5 * pa_yyy[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyy_zz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_xx[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_y_zz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yyy_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyy_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xxzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xxzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xxzz[j] * fl1_fz);

                t_yyy_xyyy[j] = fl_s_0_0 * (1.875 * pb_x[j] * fl3_fx + 2.25 * pa2pb_yy_x[j] * fl2_fx + 6.75 * pa2pb_y_xy[j] * fl2_fx + 2.25 * pb_xyy[j] * fl2_fx + 1.5 * pa2pb_yyy_xy[j] * fl1_fx + 4.5 * pa2pb_yy_xyy[j] * fl1_fx + 1.5 * pa2pb_y_xyyy[j] * fl1_fx + pa2pb_yyy_xyyy[j]);

                t_yyy_xyyy[j] += fl_r_0_0 * (15.0 * pb_x[j] * fl3_fx * fl1_fz - 2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa2pb_yy_x[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_y_xy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyy_xy[j] * fl1_fz * fl1_fgb + 22.5 * pb_xyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyy_xy[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_yy_xyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xyyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xyyy[j] * fl1_fz);

                t_yyy_xyyz[j] = fl_s_0_0 * (2.25 * pa2pb_y_xz[j] * fl2_fx + 1.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_yyy_xz[j] * fl1_fx + 3.0 * pa2pb_yy_xyz[j] * fl1_fx + 1.5 * pa2pb_y_xyyz[j] * fl1_fx + pa2pb_yyy_xyyz[j]);

                t_yyy_xyyz[j] += fl_r_0_0 * (22.5 * pa2pb_y_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_xz[j] * fl1_fz * fl1_fgb + 15.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyy_xz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yy_xyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xyyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xyyz[j] * fl1_fz);

                t_yyy_xyzz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_yy_x[j] * fl2_fx + 0.75 * pa2pb_y_xy[j] * fl2_fx + 0.75 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_yyy_xy[j] * fl1_fx + 1.5 * pa2pb_yy_xzz[j] * fl1_fx + 1.5 * pa2pb_y_xyzz[j] * fl1_fx + pa2pb_yyy_xyzz[j]);

                t_yyy_xyzz[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_x[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_xy[j] * fl1_fz * fl2_fx + 7.5 * pb_xzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyy_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yy_xzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xyzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xyzz[j] * fl1_fz);

                t_yyy_xzzz[j] = fl_s_0_0 * (2.25 * pa2pb_y_xz[j] * fl2_fx + 1.5 * pa2pb_yyy_xz[j] * fl1_fx + 1.5 * pa2pb_y_xzzz[j] * fl1_fx + pa2pb_yyy_xzzz[j]);

                t_yyy_xzzz[j] += fl_r_0_0 * (-4.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyy_xz[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_y_xz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yyy_xz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_xzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_xzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_100_105(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (100,105)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_y_yyyy = pa2pbDistances.data(646 * idx + 63);

            auto pa2pb_y_yyyz = pa2pbDistances.data(646 * idx + 64);

            auto pa2pb_y_yyzz = pa2pbDistances.data(646 * idx + 65);

            auto pa2pb_y_yzzz = pa2pbDistances.data(646 * idx + 66);

            auto pa2pb_y_zzzz = pa2pbDistances.data(646 * idx + 67);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 205);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 206);

            auto pa2pb_yy_yyy = pa2pbDistances.data(646 * idx + 219);

            auto pa2pb_yy_yyz = pa2pbDistances.data(646 * idx + 220);

            auto pa2pb_yy_yzz = pa2pbDistances.data(646 * idx + 221);

            auto pa2pb_yy_zzz = pa2pbDistances.data(646 * idx + 222);

            auto pa2pb_yyy_yy = pa2pbDistances.data(646 * idx + 516);

            auto pa2pb_yyy_yz = pa2pbDistances.data(646 * idx + 517);

            auto pa2pb_yyy_zz = pa2pbDistances.data(646 * idx + 518);

            auto pa2pb_yyy_yyyy = pa2pbDistances.data(646 * idx + 539);

            auto pa2pb_yyy_yyyz = pa2pbDistances.data(646 * idx + 540);

            auto pa2pb_yyy_yyzz = pa2pbDistances.data(646 * idx + 541);

            auto pa2pb_yyy_yzzz = pa2pbDistances.data(646 * idx + 542);

            auto pa2pb_yyy_zzzz = pa2pbDistances.data(646 * idx + 543);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyy_yyyy = primBuffer.data(150 * idx + 100);

            auto t_yyy_yyyz = primBuffer.data(150 * idx + 101);

            auto t_yyy_yyzz = primBuffer.data(150 * idx + 102);

            auto t_yyy_yzzz = primBuffer.data(150 * idx + 103);

            auto t_yyy_zzzz = primBuffer.data(150 * idx + 104);

            // Batch of Integrals (100,105)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_yy, pa2pb_y_yyyy, pa2pb_y_yyyz, pa2pb_y_yyzz, \
                                     pa2pb_y_yz, pa2pb_y_yzzz, pa2pb_y_zz, pa2pb_y_zzzz, pa2pb_yy_y, pa2pb_yy_yyy, \
                                     pa2pb_yy_yyz, pa2pb_yy_yzz, pa2pb_yy_z, pa2pb_yy_zzz, pa2pb_yyy_yy, pa2pb_yyy_yyyy, \
                                     pa2pb_yyy_yyyz, pa2pb_yyy_yyzz, pa2pb_yyy_yz, pa2pb_yyy_yzzz, pa2pb_yyy_zz, \
                                     pa2pb_yyy_zzzz, pa_y, pa_yyy, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, s_0_0, \
                                     t_yyy_yyyy, t_yyy_yyyz, t_yyy_yyzz, t_yyy_yzzz, t_yyy_zzzz: VLX_ALIGN)
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

                t_yyy_yyyy[j] = fl_s_0_0 * (5.625 * pa_y[j] * fl3_fx + 7.5 * pb_y[j] * fl3_fx + 0.75 * pa_yyy[j] * fl2_fx + 9.0 * pa2pb_yy_y[j] * fl2_fx + 13.5 * pa2pb_y_yy[j] * fl2_fx + 3.0 * pb_yyy[j] * fl2_fx + 3.0 * pa2pb_yyy_yy[j] * fl1_fx + 6.0 * pa2pb_yy_yyy[j] * fl1_fx + 1.5 * pa2pb_y_yyyy[j] * fl1_fx + pa2pb_yyy_yyyy[j]);

                t_yyy_yyyy[j] += fl_r_0_0 * (-13.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 45.0 * pa_y[j] * fl3_fx * fl1_fz + 60.0 * pb_y[j] * fl3_fx * fl1_fz - 2.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_yyy[j] * fl1_fz * fl2_fx + 90.0 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 135.0 * pa2pb_y_yy[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yyy_yy[j] * fl1_fz * fl1_fgb + 30.0 * pb_yyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yyy_yy[j] * fl1_fz * fl1_fx + 72.0 * pa2pb_yy_yyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_yyyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_yyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_yyyy[j] * fl1_fz);

                t_yyy_yyyz[j] = fl_s_0_0 * (1.875 * pb_z[j] * fl3_fx + 2.25 * pa2pb_yy_z[j] * fl2_fx + 6.75 * pa2pb_y_yz[j] * fl2_fx + 2.25 * pb_yyz[j] * fl2_fx + 1.5 * pa2pb_yyy_yz[j] * fl1_fx + 4.5 * pa2pb_yy_yyz[j] * fl1_fx + 1.5 * pa2pb_y_yyyz[j] * fl1_fx + pa2pb_yyy_yyyz[j]);

                t_yyy_yyyz[j] += fl_r_0_0 * (15.0 * pb_z[j] * fl3_fx * fl1_fz - 2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa2pb_yy_z[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_y_yz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyy_yz[j] * fl1_fz * fl1_fgb + 22.5 * pb_yyz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyy_yz[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_yy_yyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_yyyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_yyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_yyyz[j] * fl1_fz);

                t_yyy_yyzz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pb_y[j] * fl3_fx + 0.25 * pa_yyy[j] * fl2_fx + 1.5 * pa2pb_yy_y[j] * fl2_fx + 2.25 * pa2pb_y_zz[j] * fl2_fx + 0.75 * pa2pb_y_yy[j] * fl2_fx + 1.5 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_yyy_yy[j] * fl1_fx + 0.5 * pa2pb_yyy_zz[j] * fl1_fx + 3.0 * pa2pb_yy_yzz[j] * fl1_fx + 1.5 * pa2pb_y_yyzz[j] * fl1_fx + pa2pb_yyy_yyzz[j]);

                t_yyy_yyzz[j] += fl_r_0_0 * (-3.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_y[j] * fl3_fx * fl1_fz + 2.5 * pa_yyy[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyy_yy[j] * fl1_fz * fl1_fgb - pa2pb_yyy_zz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_yy[j] * fl1_fz * fl2_fx + 15.0 * pb_yzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyy_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyy_zz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yy_yzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_yyzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_yyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_yyzz[j] * fl1_fz);

                t_yyy_yzzz[j] = fl_s_0_0 * (1.125 * pb_z[j] * fl3_fx + 2.25 * pa2pb_yy_z[j] * fl2_fx + 2.25 * pa2pb_y_yz[j] * fl2_fx + 0.75 * pb_zzz[j] * fl2_fx + 1.5 * pa2pb_yyy_yz[j] * fl1_fx + 1.5 * pa2pb_yy_zzz[j] * fl1_fx + 1.5 * pa2pb_y_yzzz[j] * fl1_fx + pa2pb_yyy_yzzz[j]);

                t_yyy_yzzz[j] += fl_r_0_0 * (-2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pb_z[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_yy_z[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyy_yz[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_y_yz[j] * fl1_fz * fl2_fx + 7.5 * pb_zzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyy_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yy_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_yzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_yzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_yzzz[j] * fl1_fz);

                t_yyy_zzzz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pa_yyy[j] * fl2_fx + 4.5 * pa2pb_y_zz[j] * fl2_fx + 3.0 * pa2pb_yyy_zz[j] * fl1_fx + 1.5 * pa2pb_y_zzzz[j] * fl1_fx + pa2pb_yyy_zzzz[j]);

                t_yyy_zzzz[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl1_fz * fl3_fx + 7.5 * pa_yyy[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyy_zz[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_y_zz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_yyy_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_y_zzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_y_zzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_105_110(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (105,110)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_xxxx = pa2pbDistances.data(646 * idx + 87);

            auto pa2pb_z_xxxy = pa2pbDistances.data(646 * idx + 88);

            auto pa2pb_z_xxxz = pa2pbDistances.data(646 * idx + 89);

            auto pa2pb_z_xxyy = pa2pbDistances.data(646 * idx + 90);

            auto pa2pb_z_xxyz = pa2pbDistances.data(646 * idx + 91);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 204);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 205);

            auto pa2pb_yy_xxx = pa2pbDistances.data(646 * idx + 213);

            auto pa2pb_yy_xxy = pa2pbDistances.data(646 * idx + 214);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 238);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 239);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 240);

            auto pa2pb_yz_xxx = pa2pbDistances.data(646 * idx + 247);

            auto pa2pb_yz_xxy = pa2pbDistances.data(646 * idx + 248);

            auto pa2pb_yz_xxz = pa2pbDistances.data(646 * idx + 249);

            auto pa2pb_yyz_xx = pa2pbDistances.data(646 * idx + 547);

            auto pa2pb_yyz_xy = pa2pbDistances.data(646 * idx + 548);

            auto pa2pb_yyz_xz = pa2pbDistances.data(646 * idx + 549);

            auto pa2pb_yyz_yy = pa2pbDistances.data(646 * idx + 550);

            auto pa2pb_yyz_yz = pa2pbDistances.data(646 * idx + 551);

            auto pa2pb_yyz_xxxx = pa2pbDistances.data(646 * idx + 563);

            auto pa2pb_yyz_xxxy = pa2pbDistances.data(646 * idx + 564);

            auto pa2pb_yyz_xxxz = pa2pbDistances.data(646 * idx + 565);

            auto pa2pb_yyz_xxyy = pa2pbDistances.data(646 * idx + 566);

            auto pa2pb_yyz_xxyz = pa2pbDistances.data(646 * idx + 567);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyz_xxxx = primBuffer.data(150 * idx + 105);

            auto t_yyz_xxxy = primBuffer.data(150 * idx + 106);

            auto t_yyz_xxxz = primBuffer.data(150 * idx + 107);

            auto t_yyz_xxyy = primBuffer.data(150 * idx + 108);

            auto t_yyz_xxyz = primBuffer.data(150 * idx + 109);

            // Batch of Integrals (105,110)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_yy_x, pa2pb_yy_xxx, pa2pb_yy_xxy, \
                                     pa2pb_yy_y, pa2pb_yyz_xx, pa2pb_yyz_xxxx, pa2pb_yyz_xxxy, pa2pb_yyz_xxxz, \
                                     pa2pb_yyz_xxyy, pa2pb_yyz_xxyz, pa2pb_yyz_xy, pa2pb_yyz_xz, pa2pb_yyz_yy, \
                                     pa2pb_yyz_yz, pa2pb_yz_x, pa2pb_yz_xxx, pa2pb_yz_xxy, pa2pb_yz_xxz, pa2pb_yz_y, \
                                     pa2pb_yz_z, pa2pb_z_xx, pa2pb_z_xxxx, pa2pb_z_xxxy, pa2pb_z_xxxz, pa2pb_z_xxyy, \
                                     pa2pb_z_xxyz, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, pa_y, pa_yyz, pa_z, pb_x, \
                                     pb_xxx, pb_xxy, pb_y, r_0_0, s_0_0, t_yyz_xxxx, t_yyz_xxxy, t_yyz_xxxz, t_yyz_xxyy, \
                                     t_yyz_xxyz: VLX_ALIGN)
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

                t_yyz_xxxx[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa_yyz[j] * fl2_fx + 1.5 * pa2pb_z_xx[j] * fl2_fx + 3.0 * pa2pb_yyz_xx[j] * fl1_fx + 0.5 * pa2pb_z_xxxx[j] * fl1_fx + pa2pb_yyz_xxxx[j]);

                t_yyz_xxxx[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 7.5 * pa_yyz[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyz_xx[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yyz_xx[j] * fl1_fz * fl1_fx - pa2pb_z_xxxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xxxx[j] * fl1_fz);

                t_yyz_xxxy[j] = fl_s_0_0 * (1.5 * pa2pb_yz_x[j] * fl2_fx + 0.75 * pa2pb_z_xy[j] * fl2_fx + 1.5 * pa2pb_yyz_xy[j] * fl1_fx + pa2pb_yz_xxx[j] * fl1_fx + 0.5 * pa2pb_z_xxxy[j] * fl1_fx + pa2pb_yyz_xxxy[j]);

                t_yyz_xxxy[j] += fl_r_0_0 * (-3.0 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyz_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_xxx[j] * fl1_fx * fl1_fz - pa2pb_z_xxxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xxxy[j] * fl1_fz);

                t_yyz_xxxz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_yy_x[j] * fl2_fx + 0.75 * pa2pb_z_xz[j] * fl2_fx + 0.25 * pb_xxx[j] * fl2_fx + 1.5 * pa2pb_yyz_xz[j] * fl1_fx + 0.5 * pa2pb_yy_xxx[j] * fl1_fx + 0.5 * pa2pb_z_xxxz[j] * fl1_fx + pa2pb_yyz_xxxz[j]);

                t_yyz_xxxz[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_x[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 2.5 * pb_xxx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyz_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xxx[j] * fl1_fz * fl1_fx - pa2pb_z_xxxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xxxz[j] * fl1_fz);

                t_yyz_xxyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.25 * pa_yyz[j] * fl2_fx + pa2pb_yz_y[j] * fl2_fx + 0.75 * pa2pb_z_xx[j] * fl2_fx + 0.25 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pa2pb_yyz_xx[j] * fl1_fx + 0.5 * pa2pb_yyz_yy[j] * fl1_fx + 2.0 * pa2pb_yz_xxy[j] * fl1_fx + 0.5 * pa2pb_z_xxyy[j] * fl1_fx + pa2pb_yyz_xxyy[j]);

                t_yyz_xxyy[j] += fl_r_0_0 * (-pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz - 0.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_yyz[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyz_yy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyz_yy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yz_xxy[j] * fl1_fx * fl1_fz - pa2pb_z_xxyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xxyy[j] * fl1_fz);

                t_yyz_xxyz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl3_fx + 0.125 * pb_y[j] * fl3_fx + 0.25 * pa2pb_yy_y[j] * fl2_fx + 0.5 * pa2pb_yz_z[j] * fl2_fx + 0.5 * pa2pb_y_xx[j] * fl2_fx + 0.25 * pa2pb_z_yz[j] * fl2_fx + 0.25 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_yyz_yz[j] * fl1_fx + 0.5 * pa2pb_yy_xxy[j] * fl1_fx + pa2pb_yz_xxz[j] * fl1_fx + 0.5 * pa2pb_z_xxyz[j] * fl1_fx + pa2pb_yyz_xxyz[j]);

                t_yyz_xxyz[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb + pb_y[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_xx[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xxy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_yz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 2.5 * pb_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyz_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xxy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_xxz[j] * fl1_fx * fl1_fz - pa2pb_z_xxyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_110_115(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (110,115)

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

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_z_xxzz = pa2pbDistances.data(646 * idx + 92);

            auto pa2pb_z_xyyy = pa2pbDistances.data(646 * idx + 93);

            auto pa2pb_z_xyyz = pa2pbDistances.data(646 * idx + 94);

            auto pa2pb_z_xyzz = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_z_xzzz = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 204);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 206);

            auto pa2pb_yy_xxz = pa2pbDistances.data(646 * idx + 215);

            auto pa2pb_yy_xyy = pa2pbDistances.data(646 * idx + 216);

            auto pa2pb_yy_xyz = pa2pbDistances.data(646 * idx + 217);

            auto pa2pb_yy_xzz = pa2pbDistances.data(646 * idx + 218);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 238);

            auto pa2pb_yz_xyy = pa2pbDistances.data(646 * idx + 250);

            auto pa2pb_yz_xyz = pa2pbDistances.data(646 * idx + 251);

            auto pa2pb_yz_xzz = pa2pbDistances.data(646 * idx + 252);

            auto pa2pb_yyz_xx = pa2pbDistances.data(646 * idx + 547);

            auto pa2pb_yyz_xy = pa2pbDistances.data(646 * idx + 548);

            auto pa2pb_yyz_xz = pa2pbDistances.data(646 * idx + 549);

            auto pa2pb_yyz_zz = pa2pbDistances.data(646 * idx + 552);

            auto pa2pb_yyz_xxzz = pa2pbDistances.data(646 * idx + 568);

            auto pa2pb_yyz_xyyy = pa2pbDistances.data(646 * idx + 569);

            auto pa2pb_yyz_xyyz = pa2pbDistances.data(646 * idx + 570);

            auto pa2pb_yyz_xyzz = pa2pbDistances.data(646 * idx + 571);

            auto pa2pb_yyz_xzzz = pa2pbDistances.data(646 * idx + 572);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyz_xxzz = primBuffer.data(150 * idx + 110);

            auto t_yyz_xyyy = primBuffer.data(150 * idx + 111);

            auto t_yyz_xyyz = primBuffer.data(150 * idx + 112);

            auto t_yyz_xyzz = primBuffer.data(150 * idx + 113);

            auto t_yyz_xzzz = primBuffer.data(150 * idx + 114);

            // Batch of Integrals (110,115)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xy, pa2pb_y_xz, pa2pb_yy_x, pa2pb_yy_xxz, \
                                     pa2pb_yy_xyy, pa2pb_yy_xyz, pa2pb_yy_xzz, pa2pb_yy_z, pa2pb_yyz_xx, pa2pb_yyz_xxzz, \
                                     pa2pb_yyz_xy, pa2pb_yyz_xyyy, pa2pb_yyz_xyyz, pa2pb_yyz_xyzz, pa2pb_yyz_xz, \
                                     pa2pb_yyz_xzzz, pa2pb_yyz_zz, pa2pb_yz_x, pa2pb_yz_xyy, pa2pb_yz_xyz, pa2pb_yz_xzz, \
                                     pa2pb_z_xx, pa2pb_z_xxzz, pa2pb_z_xy, pa2pb_z_xyyy, pa2pb_z_xyyz, pa2pb_z_xyzz, \
                                     pa2pb_z_xz, pa2pb_z_xzzz, pa2pb_z_zz, pa_yyz, pa_z, pb_x, pb_xxz, pb_xyy, pb_xyz, pb_xzz, \
                                     pb_z, r_0_0, s_0_0, t_yyz_xxzz, t_yyz_xyyy, t_yyz_xyyz, t_yyz_xyzz, t_yyz_xzzz: VLX_ALIGN)
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

                t_yyz_xxzz[j] = fl_s_0_0 * (0.125 * pa_z[j] * fl3_fx + 0.25 * pb_z[j] * fl3_fx + 0.25 * pa_yyz[j] * fl2_fx + 0.5 * pa2pb_yy_z[j] * fl2_fx + 0.25 * pa2pb_z_xx[j] * fl2_fx + 0.25 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_yyz_xx[j] * fl1_fx + 0.5 * pa2pb_yyz_zz[j] * fl1_fx + pa2pb_yy_xxz[j] * fl1_fx + 0.5 * pa2pb_z_xxzz[j] * fl1_fx + pa2pb_yyz_xxzz[j]);

                t_yyz_xxzz[j] += fl_r_0_0 * (-0.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + pa_z[j] * fl3_fx * fl1_fz + 2.0 * pb_z[j] * fl3_fx * fl1_fz + 2.5 * pa_yyz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yy_z[j] * fl1_fz * fl2_fx - 0.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - pb_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyz_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 5.0 * pb_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyz_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yy_xxz[j] * fl1_fz * fl1_fx - pa2pb_z_xxzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xxzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xxzz[j] * fl1_fz);

                t_yyz_xyyy[j] = fl_s_0_0 * (1.5 * pa2pb_yz_x[j] * fl2_fx + 2.25 * pa2pb_z_xy[j] * fl2_fx + 1.5 * pa2pb_yyz_xy[j] * fl1_fx + 3.0 * pa2pb_yz_xyy[j] * fl1_fx + 0.5 * pa2pb_z_xyyy[j] * fl1_fx + pa2pb_yyz_xyyy[j]);

                t_yyz_xyyy[j] += fl_r_0_0 * (-3.0 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_z_xy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yyz_xy[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yz_xyy[j] * fl1_fx * fl1_fz - pa2pb_z_xyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xyyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xyyy[j] * fl1_fz);

                t_yyz_xyyz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.25 * pa2pb_yy_x[j] * fl2_fx + pa2pb_y_xy[j] * fl2_fx + 0.75 * pa2pb_z_xz[j] * fl2_fx + 0.25 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_yyz_xz[j] * fl1_fx + 0.5 * pa2pb_yy_xyy[j] * fl1_fx + 2.0 * pa2pb_yz_xyz[j] * fl1_fx + 0.5 * pa2pb_z_xyyz[j] * fl1_fx + pa2pb_yyz_xyyz[j]);

                t_yyz_xyyz[j] += fl_r_0_0 * (3.0 * pb_x[j] * fl3_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_yy_x[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xyy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_xz[j] * fl1_fz * fl1_fgb + 2.5 * pb_xyy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyz_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xyy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yz_xyz[j] * fl1_fx * fl1_fz - pa2pb_z_xyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xyyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xyyz[j] * fl1_fz);

                t_yyz_xyzz[j] = fl_s_0_0 * (0.5 * pa2pb_yz_x[j] * fl2_fx + pa2pb_y_xz[j] * fl2_fx + 0.25 * pa2pb_z_xy[j] * fl2_fx + 0.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_yyz_xy[j] * fl1_fx + pa2pb_yy_xyz[j] * fl1_fx + pa2pb_yz_xzz[j] * fl1_fx + 0.5 * pa2pb_z_xyzz[j] * fl1_fx + pa2pb_yyz_xyzz[j]);

                t_yyz_xyzz[j] += fl_r_0_0 * (-pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - pb_xyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_xy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 5.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyz_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yy_xyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_xzz[j] * fl1_fx * fl1_fz - pa2pb_z_xyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xyzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xyzz[j] * fl1_fz);

                t_yyz_xzzz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_yy_x[j] * fl2_fx + 0.75 * pa2pb_z_xz[j] * fl2_fx + 0.75 * pb_xzz[j] * fl2_fx + 1.5 * pa2pb_yyz_xz[j] * fl1_fx + 1.5 * pa2pb_yy_xzz[j] * fl1_fx + 0.5 * pa2pb_z_xzzz[j] * fl1_fx + pa2pb_yyz_xzzz[j]);

                t_yyz_xzzz[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_x[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 7.5 * pb_xzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyz_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yy_xzz[j] * fl1_fz * fl1_fx - pa2pb_z_xzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_xzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_115_120(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (115,120)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_z_yyyy = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_z_yyyz = pa2pbDistances.data(646 * idx + 98);

            auto pa2pb_z_yyzz = pa2pbDistances.data(646 * idx + 99);

            auto pa2pb_z_yzzz = pa2pbDistances.data(646 * idx + 100);

            auto pa2pb_z_zzzz = pa2pbDistances.data(646 * idx + 101);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 205);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 206);

            auto pa2pb_yy_yyy = pa2pbDistances.data(646 * idx + 219);

            auto pa2pb_yy_yyz = pa2pbDistances.data(646 * idx + 220);

            auto pa2pb_yy_yzz = pa2pbDistances.data(646 * idx + 221);

            auto pa2pb_yy_zzz = pa2pbDistances.data(646 * idx + 222);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 239);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 240);

            auto pa2pb_yz_yyy = pa2pbDistances.data(646 * idx + 253);

            auto pa2pb_yz_yyz = pa2pbDistances.data(646 * idx + 254);

            auto pa2pb_yz_yzz = pa2pbDistances.data(646 * idx + 255);

            auto pa2pb_yz_zzz = pa2pbDistances.data(646 * idx + 256);

            auto pa2pb_yyz_yy = pa2pbDistances.data(646 * idx + 550);

            auto pa2pb_yyz_yz = pa2pbDistances.data(646 * idx + 551);

            auto pa2pb_yyz_zz = pa2pbDistances.data(646 * idx + 552);

            auto pa2pb_yyz_yyyy = pa2pbDistances.data(646 * idx + 573);

            auto pa2pb_yyz_yyyz = pa2pbDistances.data(646 * idx + 574);

            auto pa2pb_yyz_yyzz = pa2pbDistances.data(646 * idx + 575);

            auto pa2pb_yyz_yzzz = pa2pbDistances.data(646 * idx + 576);

            auto pa2pb_yyz_zzzz = pa2pbDistances.data(646 * idx + 577);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyz_yyyy = primBuffer.data(150 * idx + 115);

            auto t_yyz_yyyz = primBuffer.data(150 * idx + 116);

            auto t_yyz_yyzz = primBuffer.data(150 * idx + 117);

            auto t_yyz_yzzz = primBuffer.data(150 * idx + 118);

            auto t_yyz_zzzz = primBuffer.data(150 * idx + 119);

            // Batch of Integrals (115,120)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yy_y, \
                                     pa2pb_yy_yyy, pa2pb_yy_yyz, pa2pb_yy_yzz, pa2pb_yy_z, pa2pb_yy_zzz, pa2pb_yyz_yy, \
                                     pa2pb_yyz_yyyy, pa2pb_yyz_yyyz, pa2pb_yyz_yyzz, pa2pb_yyz_yz, pa2pb_yyz_yzzz, \
                                     pa2pb_yyz_zz, pa2pb_yyz_zzzz, pa2pb_yz_y, pa2pb_yz_yyy, pa2pb_yz_yyz, pa2pb_yz_yzz, \
                                     pa2pb_yz_z, pa2pb_yz_zzz, pa2pb_z_yy, pa2pb_z_yyyy, pa2pb_z_yyyz, pa2pb_z_yyzz, \
                                     pa2pb_z_yz, pa2pb_z_yzzz, pa2pb_z_zz, pa2pb_z_zzzz, pa_y, pa_yyz, pa_z, pb_y, pb_yyy, \
                                     pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, s_0_0, t_yyz_yyyy, t_yyz_yyyz, t_yyz_yyzz, \
                                     t_yyz_yzzz, t_yyz_zzzz: VLX_ALIGN)
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

                t_yyz_yyyy[j] = fl_s_0_0 * (1.875 * pa_z[j] * fl3_fx + 0.75 * pa_yyz[j] * fl2_fx + 6.0 * pa2pb_yz_y[j] * fl2_fx + 4.5 * pa2pb_z_yy[j] * fl2_fx + 3.0 * pa2pb_yyz_yy[j] * fl1_fx + 4.0 * pa2pb_yz_yyy[j] * fl1_fx + 0.5 * pa2pb_z_yyyy[j] * fl1_fx + pa2pb_yyz_yyyy[j]);

                t_yyz_yyyy[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_yyz[j] * fl1_fz * fl2_fx + 60.0 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 45.0 * pa2pb_z_yy[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyz_yy[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_yyz_yy[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_yz_yyy[j] * fl1_fx * fl1_fz - pa2pb_z_yyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_yyyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_yyyy[j] * fl1_fz);

                t_yyz_yyyz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 1.125 * pb_y[j] * fl3_fx + 0.75 * pa2pb_yy_y[j] * fl2_fx + 1.5 * pa2pb_yz_z[j] * fl2_fx + 1.5 * pa2pb_y_yy[j] * fl2_fx + 2.25 * pa2pb_z_yz[j] * fl2_fx + 0.25 * pb_yyy[j] * fl2_fx + 1.5 * pa2pb_yyz_yz[j] * fl1_fx + 0.5 * pa2pb_yy_yyy[j] * fl1_fx + 3.0 * pa2pb_yz_yyz[j] * fl1_fx + 0.5 * pa2pb_z_yyyz[j] * fl1_fx + pa2pb_yyz_yyyz[j]);

                t_yyz_yyyz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl3_fx * fl1_fz + 9.0 * pb_y[j] * fl3_fx * fl1_fz - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_yz[j] * fl1_fz * fl1_fgb + 2.5 * pb_yyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyz_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_yyy[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yz_yyz[j] * fl1_fx * fl1_fz - pa2pb_z_yyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_yyyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_yyyz[j] * fl1_fz);

                t_yyz_yyzz[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pb_z[j] * fl3_fx + 0.25 * pa_yyz[j] * fl2_fx + 0.5 * pa2pb_yy_z[j] * fl2_fx + pa2pb_yz_y[j] * fl2_fx + 2.0 * pa2pb_y_yz[j] * fl2_fx + 0.75 * pa2pb_z_zz[j] * fl2_fx + 0.25 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_yyz_yy[j] * fl1_fx + 0.5 * pa2pb_yyz_zz[j] * fl1_fx + pa2pb_yy_yyz[j] * fl1_fx + 2.0 * pa2pb_yz_yzz[j] * fl1_fx + 0.5 * pa2pb_z_yyzz[j] * fl1_fx + pa2pb_yyz_yyzz[j]);

                t_yyz_yyzz[j] += fl_r_0_0 * (-pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 6.0 * pb_z[j] * fl3_fx * fl1_fz - 0.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_yyz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yy_z[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - pb_yyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyz_yy[j] * fl1_fz * fl1_fgb - pa2pb_yyz_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 5.0 * pb_yyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyz_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyz_zz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yy_yyz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yz_yzz[j] * fl1_fx * fl1_fz - pa2pb_z_yyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_yyzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_yyzz[j] * fl1_fz);

                t_yyz_yzzz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_yy_y[j] * fl2_fx + 1.5 * pa2pb_yz_z[j] * fl2_fx + 1.5 * pa2pb_y_zz[j] * fl2_fx + 0.75 * pa2pb_z_yz[j] * fl2_fx + 0.75 * pb_yzz[j] * fl2_fx + 1.5 * pa2pb_yyz_yz[j] * fl1_fx + 1.5 * pa2pb_yy_yzz[j] * fl1_fx + pa2pb_yz_zzz[j] * fl1_fx + 0.5 * pa2pb_z_yzzz[j] * fl1_fx + pa2pb_yyz_yzzz[j]);

                t_yyz_yzzz[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyz_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 7.5 * pb_yzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyz_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yy_yzz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_zzz[j] * fl1_fx * fl1_fz - pa2pb_z_yzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_yzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_yzzz[j] * fl1_fz);

                t_yyz_zzzz[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 1.5 * pb_z[j] * fl3_fx + 0.75 * pa_yyz[j] * fl2_fx + 3.0 * pa2pb_yy_z[j] * fl2_fx + 1.5 * pa2pb_z_zz[j] * fl2_fx + pb_zzz[j] * fl2_fx + 3.0 * pa2pb_yyz_zz[j] * fl1_fx + 2.0 * pa2pb_yy_zzz[j] * fl1_fx + 0.5 * pa2pb_z_zzzz[j] * fl1_fx + pa2pb_yyz_zzzz[j]);

                t_yyz_zzzz[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 12.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa_yyz[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_yy_z[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pb_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyz_zz[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 10.0 * pb_zzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yyz_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yy_zzz[j] * fl1_fz * fl1_fx - pa2pb_z_zzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_z_zzzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_120_125(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (120,125)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_xxxx = pa2pbDistances.data(646 * idx + 53);

            auto pa2pb_y_xxxy = pa2pbDistances.data(646 * idx + 54);

            auto pa2pb_y_xxxz = pa2pbDistances.data(646 * idx + 55);

            auto pa2pb_y_xxyy = pa2pbDistances.data(646 * idx + 56);

            auto pa2pb_y_xxyz = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 238);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 239);

            auto pa2pb_yz_xxx = pa2pbDistances.data(646 * idx + 247);

            auto pa2pb_yz_xxy = pa2pbDistances.data(646 * idx + 248);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_zz_xxx = pa2pbDistances.data(646 * idx + 281);

            auto pa2pb_zz_xxy = pa2pbDistances.data(646 * idx + 282);

            auto pa2pb_zz_xxz = pa2pbDistances.data(646 * idx + 283);

            auto pa2pb_yzz_xx = pa2pbDistances.data(646 * idx + 581);

            auto pa2pb_yzz_xy = pa2pbDistances.data(646 * idx + 582);

            auto pa2pb_yzz_xz = pa2pbDistances.data(646 * idx + 583);

            auto pa2pb_yzz_yy = pa2pbDistances.data(646 * idx + 584);

            auto pa2pb_yzz_yz = pa2pbDistances.data(646 * idx + 585);

            auto pa2pb_yzz_xxxx = pa2pbDistances.data(646 * idx + 597);

            auto pa2pb_yzz_xxxy = pa2pbDistances.data(646 * idx + 598);

            auto pa2pb_yzz_xxxz = pa2pbDistances.data(646 * idx + 599);

            auto pa2pb_yzz_xxyy = pa2pbDistances.data(646 * idx + 600);

            auto pa2pb_yzz_xxyz = pa2pbDistances.data(646 * idx + 601);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzz_xxxx = primBuffer.data(150 * idx + 120);

            auto t_yzz_xxxy = primBuffer.data(150 * idx + 121);

            auto t_yzz_xxxz = primBuffer.data(150 * idx + 122);

            auto t_yzz_xxyy = primBuffer.data(150 * idx + 123);

            auto t_yzz_xxyz = primBuffer.data(150 * idx + 124);

            // Batch of Integrals (120,125)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_y_xxxx, pa2pb_y_xxxy, pa2pb_y_xxxz, \
                                     pa2pb_y_xxyy, pa2pb_y_xxyz, pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, \
                                     pa2pb_yz_x, pa2pb_yz_xxx, pa2pb_yz_xxy, pa2pb_yz_y, pa2pb_yzz_xx, pa2pb_yzz_xxxx, \
                                     pa2pb_yzz_xxxy, pa2pb_yzz_xxxz, pa2pb_yzz_xxyy, pa2pb_yzz_xxyz, pa2pb_yzz_xy, \
                                     pa2pb_yzz_xz, pa2pb_yzz_yy, pa2pb_yzz_yz, pa2pb_z_xx, pa2pb_zz_x, pa2pb_zz_xxx, \
                                     pa2pb_zz_xxy, pa2pb_zz_xxz, pa2pb_zz_y, pa2pb_zz_z, pa_y, pa_yzz, pa_z, pb_x, pb_xxx, \
                                     pb_xxy, pb_xxz, pb_y, pb_z, r_0_0, s_0_0, t_yzz_xxxx, t_yzz_xxxy, t_yzz_xxxz, \
                                     t_yzz_xxyy, t_yzz_xxyz: VLX_ALIGN)
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

                t_yzz_xxxx[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa_yzz[j] * fl2_fx + 1.5 * pa2pb_y_xx[j] * fl2_fx + 3.0 * pa2pb_yzz_xx[j] * fl1_fx + 0.5 * pa2pb_y_xxxx[j] * fl1_fx + pa2pb_yzz_xxxx[j]);

                t_yzz_xxxx[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl1_fz * fl3_fx + 7.5 * pa_yzz[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yzz_xx[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_xx[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_yzz_xx[j] * fl1_fz * fl1_fx - pa2pb_y_xxxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xxxx[j] * fl1_fz);

                t_yzz_xxxy[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_zz_x[j] * fl2_fx + 0.75 * pa2pb_y_xy[j] * fl2_fx + 0.25 * pb_xxx[j] * fl2_fx + 1.5 * pa2pb_yzz_xy[j] * fl1_fx + 0.5 * pa2pb_zz_xxx[j] * fl1_fx + 0.5 * pa2pb_y_xxxy[j] * fl1_fx + pa2pb_yzz_xxxy[j]);

                t_yzz_xxxy[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzz_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_xy[j] * fl1_fz * fl2_fx + 2.5 * pb_xxx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yzz_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xxx[j] * fl1_fx * fl1_fz - pa2pb_y_xxxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xxxy[j] * fl1_fz);

                t_yzz_xxxz[j] = fl_s_0_0 * (1.5 * pa2pb_yz_x[j] * fl2_fx + 0.75 * pa2pb_y_xz[j] * fl2_fx + 1.5 * pa2pb_yzz_xz[j] * fl1_fx + pa2pb_yz_xxx[j] * fl1_fx + 0.5 * pa2pb_y_xxxz[j] * fl1_fx + pa2pb_yzz_xxxz[j]);

                t_yzz_xxxz[j] += fl_r_0_0 * (-3.0 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_yz_x[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yzz_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_xz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yzz_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_xxx[j] * fl1_fz * fl1_fx - pa2pb_y_xxxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xxxz[j] * fl1_fz);

                t_yzz_xxyy[j] = fl_s_0_0 * (0.125 * pa_y[j] * fl3_fx + 0.25 * pb_y[j] * fl3_fx + 0.25 * pa_yzz[j] * fl2_fx + 0.5 * pa2pb_zz_y[j] * fl2_fx + 0.25 * pa2pb_y_xx[j] * fl2_fx + 0.25 * pa2pb_y_yy[j] * fl2_fx + 0.5 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_yzz_xx[j] * fl1_fx + 0.5 * pa2pb_yzz_yy[j] * fl1_fx + pa2pb_zz_xxy[j] * fl1_fx + 0.5 * pa2pb_y_xxyy[j] * fl1_fx + pa2pb_yzz_xxyy[j]);

                t_yzz_xxyy[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + pa_y[j] * fl1_fz * fl3_fx - pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pb_y[j] * fl3_fx * fl1_fz + 2.5 * pa_yzz[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yzz_yy[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_xx[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_y_yy[j] * fl1_fz * fl2_fx + 5.0 * pb_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yzz_yy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zz_xxy[j] * fl1_fx * fl1_fz - pa2pb_y_xxyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xxyy[j] * fl1_fz);

                t_yzz_xxyz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl3_fx + 0.125 * pb_z[j] * fl3_fx + 0.5 * pa2pb_yz_y[j] * fl2_fx + 0.25 * pa2pb_zz_z[j] * fl2_fx + 0.5 * pa2pb_z_xx[j] * fl2_fx + 0.25 * pa2pb_y_yz[j] * fl2_fx + 0.25 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_yzz_yz[j] * fl1_fx + pa2pb_yz_xxy[j] * fl1_fx + 0.5 * pa2pb_zz_xxz[j] * fl1_fx + 0.5 * pa2pb_y_xxyz[j] * fl1_fx + pa2pb_yzz_xxyz[j]);

                t_yzz_xxyz[j] += fl_r_0_0 * (-0.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 2.0 * pa_z[j] * fl3_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + pb_z[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_yz_y[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xx[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_yz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_yz[j] * fl1_fz * fl2_fx + 2.5 * pb_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzz_yz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_xxy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xxz[j] * fl1_fx * fl1_fz - pa2pb_y_xxyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_125_130(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (125,130)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 37);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 38);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 39);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_y_xxzz = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_y_xyyy = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_y_xyyz = pa2pbDistances.data(646 * idx + 60);

            auto pa2pb_y_xyzz = pa2pbDistances.data(646 * idx + 61);

            auto pa2pb_y_xzzz = pa2pbDistances.data(646 * idx + 62);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 238);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 240);

            auto pa2pb_yz_xxz = pa2pbDistances.data(646 * idx + 249);

            auto pa2pb_yz_xyy = pa2pbDistances.data(646 * idx + 250);

            auto pa2pb_yz_xyz = pa2pbDistances.data(646 * idx + 251);

            auto pa2pb_yz_xzz = pa2pbDistances.data(646 * idx + 252);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_zz_xyy = pa2pbDistances.data(646 * idx + 284);

            auto pa2pb_zz_xyz = pa2pbDistances.data(646 * idx + 285);

            auto pa2pb_zz_xzz = pa2pbDistances.data(646 * idx + 286);

            auto pa2pb_yzz_xx = pa2pbDistances.data(646 * idx + 581);

            auto pa2pb_yzz_xy = pa2pbDistances.data(646 * idx + 582);

            auto pa2pb_yzz_xz = pa2pbDistances.data(646 * idx + 583);

            auto pa2pb_yzz_zz = pa2pbDistances.data(646 * idx + 586);

            auto pa2pb_yzz_xxzz = pa2pbDistances.data(646 * idx + 602);

            auto pa2pb_yzz_xyyy = pa2pbDistances.data(646 * idx + 603);

            auto pa2pb_yzz_xyyz = pa2pbDistances.data(646 * idx + 604);

            auto pa2pb_yzz_xyzz = pa2pbDistances.data(646 * idx + 605);

            auto pa2pb_yzz_xzzz = pa2pbDistances.data(646 * idx + 606);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzz_xxzz = primBuffer.data(150 * idx + 125);

            auto t_yzz_xyyy = primBuffer.data(150 * idx + 126);

            auto t_yzz_xyyz = primBuffer.data(150 * idx + 127);

            auto t_yzz_xyzz = primBuffer.data(150 * idx + 128);

            auto t_yzz_xzzz = primBuffer.data(150 * idx + 129);

            // Batch of Integrals (125,130)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_y_xxzz, pa2pb_y_xy, pa2pb_y_xyyy, \
                                     pa2pb_y_xyyz, pa2pb_y_xyzz, pa2pb_y_xz, pa2pb_y_xzzz, pa2pb_y_zz, pa2pb_yz_x, \
                                     pa2pb_yz_xxz, pa2pb_yz_xyy, pa2pb_yz_xyz, pa2pb_yz_xzz, pa2pb_yz_z, pa2pb_yzz_xx, \
                                     pa2pb_yzz_xxzz, pa2pb_yzz_xy, pa2pb_yzz_xyyy, pa2pb_yzz_xyyz, pa2pb_yzz_xyzz, \
                                     pa2pb_yzz_xz, pa2pb_yzz_xzzz, pa2pb_yzz_zz, pa2pb_z_xy, pa2pb_z_xz, pa2pb_zz_x, \
                                     pa2pb_zz_xyy, pa2pb_zz_xyz, pa2pb_zz_xzz, pa_y, pa_yzz, pb_x, pb_xyy, pb_xyz, pb_xzz, \
                                     r_0_0, s_0_0, t_yzz_xxzz, t_yzz_xyyy, t_yzz_xyyz, t_yzz_xyzz, t_yzz_xzzz: VLX_ALIGN)
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

                t_yzz_xxzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.25 * pa_yzz[j] * fl2_fx + pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_xx[j] * fl2_fx + 0.25 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pa2pb_yzz_xx[j] * fl1_fx + 0.5 * pa2pb_yzz_zz[j] * fl1_fx + 2.0 * pa2pb_yz_xxz[j] * fl1_fx + 0.5 * pa2pb_y_xxzz[j] * fl1_fx + pa2pb_yzz_xxzz[j]);

                t_yzz_xxzz[j] += fl_r_0_0 * (-pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_yzz[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_yz_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yzz_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_zz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yzz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yzz_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yz_xxz[j] * fl1_fz * fl1_fx - pa2pb_y_xxzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xxzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xxzz[j] * fl1_fz);

                t_yzz_xyyy[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_zz_x[j] * fl2_fx + 0.75 * pa2pb_y_xy[j] * fl2_fx + 0.75 * pb_xyy[j] * fl2_fx + 1.5 * pa2pb_yzz_xy[j] * fl1_fx + 1.5 * pa2pb_zz_xyy[j] * fl1_fx + 0.5 * pa2pb_y_xyyy[j] * fl1_fx + pa2pb_yzz_xyyy[j]);

                t_yzz_xyyy[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzz_xy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_xy[j] * fl1_fz * fl2_fx + 7.5 * pb_xyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yzz_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zz_xyy[j] * fl1_fx * fl1_fz - pa2pb_y_xyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xyyy[j] * fl1_fz);

                t_yzz_xyyz[j] = fl_s_0_0 * (0.5 * pa2pb_yz_x[j] * fl2_fx + pa2pb_z_xy[j] * fl2_fx + 0.25 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_yzz_xz[j] * fl1_fx + pa2pb_yz_xyy[j] * fl1_fx + pa2pb_zz_xyz[j] * fl1_fx + 0.5 * pa2pb_y_xyyz[j] * fl1_fx + pa2pb_yzz_xyyz[j]);

                t_yzz_xyyz[j] += fl_r_0_0 * (-pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_yz_x[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_xz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_xz[j] * fl1_fz * fl2_fx + 5.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzz_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_xyy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zz_xyz[j] * fl1_fx * fl1_fz - pa2pb_y_xyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xyyz[j] * fl1_fz);

                t_yzz_xyzz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_y_xy[j] * fl2_fx + 0.25 * pa2pb_zz_x[j] * fl2_fx + pa2pb_z_xz[j] * fl2_fx + 0.25 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_yzz_xy[j] * fl1_fx + 2.0 * pa2pb_yz_xyz[j] * fl1_fx + 0.5 * pa2pb_zz_xzz[j] * fl1_fx + 0.5 * pa2pb_y_xyzz[j] * fl1_fx + pa2pb_yzz_xyzz[j]);

                t_yzz_xyzz[j] += fl_r_0_0 * (3.0 * pb_x[j] * fl3_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_xy[j] * fl1_fz * fl1_fgb + 2.5 * pb_xzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzz_xy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yz_xyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xzz[j] * fl1_fx * fl1_fz - pa2pb_y_xyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xyzz[j] * fl1_fz);

                t_yzz_xzzz[j] = fl_s_0_0 * (1.5 * pa2pb_yz_x[j] * fl2_fx + 2.25 * pa2pb_y_xz[j] * fl2_fx + 1.5 * pa2pb_yzz_xz[j] * fl1_fx + 3.0 * pa2pb_yz_xzz[j] * fl1_fx + 0.5 * pa2pb_y_xzzz[j] * fl1_fx + pa2pb_yzz_xzzz[j]);

                t_yzz_xzzz[j] += fl_r_0_0 * (-3.0 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_yz_x[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_y_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yzz_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yzz_xz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yz_xzz[j] * fl1_fz * fl1_fx - pa2pb_y_xzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_xzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_130_135(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (130,135)

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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 40);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_y_yyyy = pa2pbDistances.data(646 * idx + 63);

            auto pa2pb_y_yyyz = pa2pbDistances.data(646 * idx + 64);

            auto pa2pb_y_yyzz = pa2pbDistances.data(646 * idx + 65);

            auto pa2pb_y_yzzz = pa2pbDistances.data(646 * idx + 66);

            auto pa2pb_y_zzzz = pa2pbDistances.data(646 * idx + 67);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 239);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 240);

            auto pa2pb_yz_yyy = pa2pbDistances.data(646 * idx + 253);

            auto pa2pb_yz_yyz = pa2pbDistances.data(646 * idx + 254);

            auto pa2pb_yz_yzz = pa2pbDistances.data(646 * idx + 255);

            auto pa2pb_yz_zzz = pa2pbDistances.data(646 * idx + 256);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_zz_yyy = pa2pbDistances.data(646 * idx + 287);

            auto pa2pb_zz_yyz = pa2pbDistances.data(646 * idx + 288);

            auto pa2pb_zz_yzz = pa2pbDistances.data(646 * idx + 289);

            auto pa2pb_zz_zzz = pa2pbDistances.data(646 * idx + 290);

            auto pa2pb_yzz_yy = pa2pbDistances.data(646 * idx + 584);

            auto pa2pb_yzz_yz = pa2pbDistances.data(646 * idx + 585);

            auto pa2pb_yzz_zz = pa2pbDistances.data(646 * idx + 586);

            auto pa2pb_yzz_yyyy = pa2pbDistances.data(646 * idx + 607);

            auto pa2pb_yzz_yyyz = pa2pbDistances.data(646 * idx + 608);

            auto pa2pb_yzz_yyzz = pa2pbDistances.data(646 * idx + 609);

            auto pa2pb_yzz_yzzz = pa2pbDistances.data(646 * idx + 610);

            auto pa2pb_yzz_zzzz = pa2pbDistances.data(646 * idx + 611);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzz_yyyy = primBuffer.data(150 * idx + 130);

            auto t_yzz_yyyz = primBuffer.data(150 * idx + 131);

            auto t_yzz_yyzz = primBuffer.data(150 * idx + 132);

            auto t_yzz_yzzz = primBuffer.data(150 * idx + 133);

            auto t_yzz_zzzz = primBuffer.data(150 * idx + 134);

            // Batch of Integrals (130,135)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_yy, pa2pb_y_yyyy, pa2pb_y_yyyz, pa2pb_y_yyzz, \
                                     pa2pb_y_yz, pa2pb_y_yzzz, pa2pb_y_zz, pa2pb_y_zzzz, pa2pb_yz_y, pa2pb_yz_yyy, \
                                     pa2pb_yz_yyz, pa2pb_yz_yzz, pa2pb_yz_z, pa2pb_yz_zzz, pa2pb_yzz_yy, pa2pb_yzz_yyyy, \
                                     pa2pb_yzz_yyyz, pa2pb_yzz_yyzz, pa2pb_yzz_yz, pa2pb_yzz_yzzz, pa2pb_yzz_zz, \
                                     pa2pb_yzz_zzzz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa2pb_zz_y, pa2pb_zz_yyy, \
                                     pa2pb_zz_yyz, pa2pb_zz_yzz, pa2pb_zz_z, pa2pb_zz_zzz, pa_y, pa_yzz, pa_z, pb_y, pb_yyy, \
                                     pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, s_0_0, t_yzz_yyyy, t_yzz_yyyz, t_yzz_yyzz, \
                                     t_yzz_yzzz, t_yzz_zzzz: VLX_ALIGN)
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

                t_yzz_yyyy[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 1.5 * pb_y[j] * fl3_fx + 0.75 * pa_yzz[j] * fl2_fx + 3.0 * pa2pb_zz_y[j] * fl2_fx + 1.5 * pa2pb_y_yy[j] * fl2_fx + pb_yyy[j] * fl2_fx + 3.0 * pa2pb_yzz_yy[j] * fl1_fx + 2.0 * pa2pb_zz_yyy[j] * fl1_fx + 0.5 * pa2pb_y_yyyy[j] * fl1_fx + pa2pb_yzz_yyyy[j]);

                t_yzz_yyyy[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl1_fz * fl3_fx - 6.0 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa_yzz[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yzz_yy[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_yy[j] * fl1_fz * fl2_fx + 10.0 * pb_yyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yzz_yy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_zz_yyy[j] * fl1_fx * fl1_fz - pa2pb_y_yyyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_yyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_yyyy[j] * fl1_fz);

                t_yzz_yyyz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl3_fx + 0.375 * pb_z[j] * fl3_fx + 1.5 * pa2pb_yz_y[j] * fl2_fx + 0.75 * pa2pb_zz_z[j] * fl2_fx + 1.5 * pa2pb_z_yy[j] * fl2_fx + 0.75 * pa2pb_y_yz[j] * fl2_fx + 0.75 * pb_yyz[j] * fl2_fx + 1.5 * pa2pb_yzz_yz[j] * fl1_fx + pa2pb_yz_yyy[j] * fl1_fx + 1.5 * pa2pb_zz_yyz[j] * fl1_fx + 0.5 * pa2pb_y_yyyz[j] * fl1_fx + pa2pb_yzz_yyyz[j]);

                t_yzz_yyyz[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_yz_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzz_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_yz[j] * fl1_fz * fl2_fx + 7.5 * pb_yyz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yzz_yz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yz_yyy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zz_yyz[j] * fl1_fx * fl1_fz - pa2pb_y_yyyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_yyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_yyyz[j] * fl1_fz);

                t_yzz_yyzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pb_y[j] * fl3_fx + 0.25 * pa_yzz[j] * fl2_fx + pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_yy[j] * fl2_fx + 0.5 * pa2pb_zz_y[j] * fl2_fx + 2.0 * pa2pb_z_yz[j] * fl2_fx + 0.25 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_yzz_yy[j] * fl1_fx + 0.5 * pa2pb_yzz_zz[j] * fl1_fx + 2.0 * pa2pb_yz_yyz[j] * fl1_fx + pa2pb_zz_yzz[j] * fl1_fx + 0.5 * pa2pb_y_yyzz[j] * fl1_fx + pa2pb_yzz_yyzz[j]);

                t_yzz_yyzz[j] += fl_r_0_0 * (-pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 6.0 * pb_y[j] * fl3_fx * fl1_fz - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_yzz[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_yz_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_zz_y[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_yzz_zz[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_zz[j] * fl1_fz * fl2_fx + 5.0 * pb_yzz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzz_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yzz_zz[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yz_yyz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zz_yzz[j] * fl1_fx * fl1_fz - pa2pb_y_yyzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_yyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_yyzz[j] * fl1_fz);

                t_yzz_yzzz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl3_fx + 1.125 * pb_z[j] * fl3_fx + 1.5 * pa2pb_yz_y[j] * fl2_fx + 2.25 * pa2pb_y_yz[j] * fl2_fx + 0.75 * pa2pb_zz_z[j] * fl2_fx + 1.5 * pa2pb_z_zz[j] * fl2_fx + 0.25 * pb_zzz[j] * fl2_fx + 1.5 * pa2pb_yzz_yz[j] * fl1_fx + 3.0 * pa2pb_yz_yzz[j] * fl1_fx + 0.5 * pa2pb_zz_zzz[j] * fl1_fx + 0.5 * pa2pb_y_yzzz[j] * fl1_fx + pa2pb_yzz_yzzz[j]);

                t_yzz_yzzz[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 6.0 * pa_z[j] * fl3_fx * fl1_fz + 9.0 * pb_z[j] * fl3_fx * fl1_fz - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_yz_y[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzz_yz[j] * fl1_fz * fl1_fgb + 2.5 * pb_zzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yzz_yz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yz_yzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_zzz[j] * fl1_fx * fl1_fz - pa2pb_y_yzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_yzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_yzzz[j] * fl1_fz);

                t_yzz_zzzz[j] = fl_s_0_0 * (1.875 * pa_y[j] * fl3_fx + 0.75 * pa_yzz[j] * fl2_fx + 6.0 * pa2pb_yz_z[j] * fl2_fx + 4.5 * pa2pb_y_zz[j] * fl2_fx + 3.0 * pa2pb_yzz_zz[j] * fl1_fx + 4.0 * pa2pb_yz_zzz[j] * fl1_fx + 0.5 * pa2pb_y_zzzz[j] * fl1_fx + pa2pb_yzz_zzzz[j]);

                t_yzz_zzzz[j] += fl_r_0_0 * (-4.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_yzz[j] * fl1_fz * fl2_fx + 60.0 * pa2pb_yz_z[j] * fl1_fz * fl2_fx + 45.0 * pa2pb_y_zz[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yzz_zz[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_yzz_zz[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_yz_zzz[j] * fl1_fz * fl1_fx - pa2pb_y_zzzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_y_zzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_135_140(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (135,140)

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

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_xxxx = pa2pbDistances.data(646 * idx + 87);

            auto pa2pb_z_xxxy = pa2pbDistances.data(646 * idx + 88);

            auto pa2pb_z_xxxz = pa2pbDistances.data(646 * idx + 89);

            auto pa2pb_z_xxyy = pa2pbDistances.data(646 * idx + 90);

            auto pa2pb_z_xxyz = pa2pbDistances.data(646 * idx + 91);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_zz_xxx = pa2pbDistances.data(646 * idx + 281);

            auto pa2pb_zz_xxy = pa2pbDistances.data(646 * idx + 282);

            auto pa2pb_zzz_xx = pa2pbDistances.data(646 * idx + 615);

            auto pa2pb_zzz_xy = pa2pbDistances.data(646 * idx + 616);

            auto pa2pb_zzz_xz = pa2pbDistances.data(646 * idx + 617);

            auto pa2pb_zzz_yy = pa2pbDistances.data(646 * idx + 618);

            auto pa2pb_zzz_yz = pa2pbDistances.data(646 * idx + 619);

            auto pa2pb_zzz_xxxx = pa2pbDistances.data(646 * idx + 631);

            auto pa2pb_zzz_xxxy = pa2pbDistances.data(646 * idx + 632);

            auto pa2pb_zzz_xxxz = pa2pbDistances.data(646 * idx + 633);

            auto pa2pb_zzz_xxyy = pa2pbDistances.data(646 * idx + 634);

            auto pa2pb_zzz_xxyz = pa2pbDistances.data(646 * idx + 635);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzz_xxxx = primBuffer.data(150 * idx + 135);

            auto t_zzz_xxxy = primBuffer.data(150 * idx + 136);

            auto t_zzz_xxxz = primBuffer.data(150 * idx + 137);

            auto t_zzz_xxyy = primBuffer.data(150 * idx + 138);

            auto t_zzz_xxyz = primBuffer.data(150 * idx + 139);

            // Batch of Integrals (135,140)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_xx, pa2pb_z_xxxx, pa2pb_z_xxxy, pa2pb_z_xxxz, \
                                     pa2pb_z_xxyy, pa2pb_z_xxyz, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, \
                                     pa2pb_zz_x, pa2pb_zz_xxx, pa2pb_zz_xxy, pa2pb_zz_y, pa2pb_zzz_xx, pa2pb_zzz_xxxx, \
                                     pa2pb_zzz_xxxy, pa2pb_zzz_xxxz, pa2pb_zzz_xxyy, pa2pb_zzz_xxyz, pa2pb_zzz_xy, \
                                     pa2pb_zzz_xz, pa2pb_zzz_yy, pa2pb_zzz_yz, pa_z, pa_zzz, pb_x, pb_xxx, pb_xxy, pb_y, r_0_0, \
                                     s_0_0, t_zzz_xxxx, t_zzz_xxxy, t_zzz_xxxz, t_zzz_xxyy, t_zzz_xxyz: VLX_ALIGN)
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

                t_zzz_xxxx[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pa_zzz[j] * fl2_fx + 4.5 * pa2pb_z_xx[j] * fl2_fx + 3.0 * pa2pb_zzz_xx[j] * fl1_fx + 1.5 * pa2pb_z_xxxx[j] * fl1_fx + pa2pb_zzz_xxxx[j]);

                t_zzz_xxxx[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl1_fz * fl3_fx + 7.5 * pa_zzz[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_zzz_xx[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_z_xx[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_zzz_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xxxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xxxx[j] * fl1_fz);

                t_zzz_xxxy[j] = fl_s_0_0 * (2.25 * pa2pb_z_xy[j] * fl2_fx + 1.5 * pa2pb_zzz_xy[j] * fl1_fx + 1.5 * pa2pb_z_xxxy[j] * fl1_fx + pa2pb_zzz_xxxy[j]);

                t_zzz_xxxy[j] += fl_r_0_0 * (-4.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zzz_xy[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_z_xy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zzz_xy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xxxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xxxy[j] * fl1_fz);

                t_zzz_xxxz[j] = fl_s_0_0 * (1.125 * pb_x[j] * fl3_fx + 2.25 * pa2pb_zz_x[j] * fl2_fx + 2.25 * pa2pb_z_xz[j] * fl2_fx + 0.75 * pb_xxx[j] * fl2_fx + 1.5 * pa2pb_zzz_xz[j] * fl1_fx + 1.5 * pa2pb_zz_xxx[j] * fl1_fx + 1.5 * pa2pb_z_xxxz[j] * fl1_fx + pa2pb_zzz_xxxz[j]);

                t_zzz_xxxz[j] += fl_r_0_0 * (-2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pb_x[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_zz_x[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zzz_xz[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_z_xz[j] * fl1_fz * fl2_fx + 7.5 * pb_xxx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_zzz_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zz_xxx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xxxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xxxz[j] * fl1_fz);

                t_zzz_xxyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.25 * pa_zzz[j] * fl2_fx + 0.75 * pa2pb_z_xx[j] * fl2_fx + 0.75 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pa2pb_zzz_xx[j] * fl1_fx + 0.5 * pa2pb_zzz_yy[j] * fl1_fx + 1.5 * pa2pb_z_xxyy[j] * fl1_fx + pa2pb_zzz_xxyy[j]);

                t_zzz_xxyy[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl1_fz * fl3_fx + 2.5 * pa_zzz[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_zzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_zzz_yy[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_xx[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_z_yy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_zzz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xxyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xxyy[j] * fl1_fz);

                t_zzz_xxyz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_zz_y[j] * fl2_fx + 0.75 * pa2pb_z_yz[j] * fl2_fx + 0.75 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_zzz_yz[j] * fl1_fx + 1.5 * pa2pb_zz_xxy[j] * fl1_fx + 1.5 * pa2pb_z_xxyz[j] * fl1_fx + pa2pb_zzz_xxyz[j]);

                t_zzz_xxyz[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_y[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_yz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_yz[j] * fl1_fz * fl2_fx + 7.5 * pb_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzz_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zz_xxy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xxyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_140_145(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (140,145)

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

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_z_xxzz = pa2pbDistances.data(646 * idx + 92);

            auto pa2pb_z_xyyy = pa2pbDistances.data(646 * idx + 93);

            auto pa2pb_z_xyyz = pa2pbDistances.data(646 * idx + 94);

            auto pa2pb_z_xyzz = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_z_xzzz = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_zz_xxz = pa2pbDistances.data(646 * idx + 283);

            auto pa2pb_zz_xyy = pa2pbDistances.data(646 * idx + 284);

            auto pa2pb_zz_xyz = pa2pbDistances.data(646 * idx + 285);

            auto pa2pb_zz_xzz = pa2pbDistances.data(646 * idx + 286);

            auto pa2pb_zzz_xx = pa2pbDistances.data(646 * idx + 615);

            auto pa2pb_zzz_xy = pa2pbDistances.data(646 * idx + 616);

            auto pa2pb_zzz_xz = pa2pbDistances.data(646 * idx + 617);

            auto pa2pb_zzz_zz = pa2pbDistances.data(646 * idx + 620);

            auto pa2pb_zzz_xxzz = pa2pbDistances.data(646 * idx + 636);

            auto pa2pb_zzz_xyyy = pa2pbDistances.data(646 * idx + 637);

            auto pa2pb_zzz_xyyz = pa2pbDistances.data(646 * idx + 638);

            auto pa2pb_zzz_xyzz = pa2pbDistances.data(646 * idx + 639);

            auto pa2pb_zzz_xzzz = pa2pbDistances.data(646 * idx + 640);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzz_xxzz = primBuffer.data(150 * idx + 140);

            auto t_zzz_xyyy = primBuffer.data(150 * idx + 141);

            auto t_zzz_xyyz = primBuffer.data(150 * idx + 142);

            auto t_zzz_xyzz = primBuffer.data(150 * idx + 143);

            auto t_zzz_xzzz = primBuffer.data(150 * idx + 144);

            // Batch of Integrals (140,145)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_xx, pa2pb_z_xxzz, pa2pb_z_xy, pa2pb_z_xyyy, \
                                     pa2pb_z_xyyz, pa2pb_z_xyzz, pa2pb_z_xz, pa2pb_z_xzzz, pa2pb_z_zz, pa2pb_zz_x, \
                                     pa2pb_zz_xxz, pa2pb_zz_xyy, pa2pb_zz_xyz, pa2pb_zz_xzz, pa2pb_zz_z, pa2pb_zzz_xx, \
                                     pa2pb_zzz_xxzz, pa2pb_zzz_xy, pa2pb_zzz_xyyy, pa2pb_zzz_xyyz, pa2pb_zzz_xyzz, \
                                     pa2pb_zzz_xz, pa2pb_zzz_xzzz, pa2pb_zzz_zz, pa_z, pa_zzz, pb_x, pb_xxz, pb_xyy, pb_xyz, \
                                     pb_xzz, pb_z, r_0_0, s_0_0, t_zzz_xxzz, t_zzz_xyyy, t_zzz_xyyz, t_zzz_xyzz, \
                                     t_zzz_xzzz: VLX_ALIGN)
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

                t_zzz_xxzz[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pb_z[j] * fl3_fx + 0.25 * pa_zzz[j] * fl2_fx + 1.5 * pa2pb_zz_z[j] * fl2_fx + 2.25 * pa2pb_z_xx[j] * fl2_fx + 0.75 * pa2pb_z_zz[j] * fl2_fx + 1.5 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_zzz_xx[j] * fl1_fx + 0.5 * pa2pb_zzz_zz[j] * fl1_fx + 3.0 * pa2pb_zz_xxz[j] * fl1_fx + 1.5 * pa2pb_z_xxzz[j] * fl1_fx + pa2pb_zzz_xxzz[j]);

                t_zzz_xxzz[j] += fl_r_0_0 * (-3.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_z[j] * fl3_fx * fl1_fz + 2.5 * pa_zzz[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_zz_z[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_zzz_zz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_zz[j] * fl1_fz * fl2_fx + 15.0 * pb_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzz_xx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_zz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_zz_xxz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xxzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xxzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xxzz[j] * fl1_fz);

                t_zzz_xyyy[j] = fl_s_0_0 * (2.25 * pa2pb_z_xy[j] * fl2_fx + 1.5 * pa2pb_zzz_xy[j] * fl1_fx + 1.5 * pa2pb_z_xyyy[j] * fl1_fx + pa2pb_zzz_xyyy[j]);

                t_zzz_xyyy[j] += fl_r_0_0 * (-4.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zzz_xy[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_z_xy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zzz_xy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xyyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xyyy[j] * fl1_fz);

                t_zzz_xyyz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_zz_x[j] * fl2_fx + 0.75 * pa2pb_z_xz[j] * fl2_fx + 0.75 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_zzz_xz[j] * fl1_fx + 1.5 * pa2pb_zz_xyy[j] * fl1_fx + 1.5 * pa2pb_z_xyyz[j] * fl1_fx + pa2pb_zzz_xyyz[j]);

                t_zzz_xyyz[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_x[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_xz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_xz[j] * fl1_fz * fl2_fx + 7.5 * pb_xyy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzz_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zz_xyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xyyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xyyz[j] * fl1_fz);

                t_zzz_xyzz[j] = fl_s_0_0 * (2.25 * pa2pb_z_xy[j] * fl2_fx + 1.5 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_zzz_xy[j] * fl1_fx + 3.0 * pa2pb_zz_xyz[j] * fl1_fx + 1.5 * pa2pb_z_xyzz[j] * fl1_fx + pa2pb_zzz_xyzz[j]);

                t_zzz_xyzz[j] += fl_r_0_0 * (22.5 * pa2pb_z_xy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_xy[j] * fl1_fz * fl1_fgb + 15.0 * pb_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzz_xy[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_zz_xyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xyzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xyzz[j] * fl1_fz);

                t_zzz_xzzz[j] = fl_s_0_0 * (1.875 * pb_x[j] * fl3_fx + 2.25 * pa2pb_zz_x[j] * fl2_fx + 6.75 * pa2pb_z_xz[j] * fl2_fx + 2.25 * pb_xzz[j] * fl2_fx + 1.5 * pa2pb_zzz_xz[j] * fl1_fx + 4.5 * pa2pb_zz_xzz[j] * fl1_fx + 1.5 * pa2pb_z_xzzz[j] * fl1_fx + pa2pb_zzz_xzzz[j]);

                t_zzz_xzzz[j] += fl_r_0_0 * (15.0 * pb_x[j] * fl3_fx * fl1_fz - 2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa2pb_zz_x[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zzz_xz[j] * fl1_fz * fl1_fgb + 22.5 * pb_xzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_zzz_xz[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_zz_xzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_xzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_xzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFG_145_150(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (145,150)

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

            auto pa_z = paDistances.data(19 * idx + 2);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_z_yyyy = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_z_yyyz = pa2pbDistances.data(646 * idx + 98);

            auto pa2pb_z_yyzz = pa2pbDistances.data(646 * idx + 99);

            auto pa2pb_z_yzzz = pa2pbDistances.data(646 * idx + 100);

            auto pa2pb_z_zzzz = pa2pbDistances.data(646 * idx + 101);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_zz_yyy = pa2pbDistances.data(646 * idx + 287);

            auto pa2pb_zz_yyz = pa2pbDistances.data(646 * idx + 288);

            auto pa2pb_zz_yzz = pa2pbDistances.data(646 * idx + 289);

            auto pa2pb_zz_zzz = pa2pbDistances.data(646 * idx + 290);

            auto pa2pb_zzz_yy = pa2pbDistances.data(646 * idx + 618);

            auto pa2pb_zzz_yz = pa2pbDistances.data(646 * idx + 619);

            auto pa2pb_zzz_zz = pa2pbDistances.data(646 * idx + 620);

            auto pa2pb_zzz_yyyy = pa2pbDistances.data(646 * idx + 641);

            auto pa2pb_zzz_yyyz = pa2pbDistances.data(646 * idx + 642);

            auto pa2pb_zzz_yyzz = pa2pbDistances.data(646 * idx + 643);

            auto pa2pb_zzz_yzzz = pa2pbDistances.data(646 * idx + 644);

            auto pa2pb_zzz_zzzz = pa2pbDistances.data(646 * idx + 645);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzz_yyyy = primBuffer.data(150 * idx + 145);

            auto t_zzz_yyyz = primBuffer.data(150 * idx + 146);

            auto t_zzz_yyzz = primBuffer.data(150 * idx + 147);

            auto t_zzz_yzzz = primBuffer.data(150 * idx + 148);

            auto t_zzz_zzzz = primBuffer.data(150 * idx + 149);

            // Batch of Integrals (145,150)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_yy, pa2pb_z_yyyy, pa2pb_z_yyyz, pa2pb_z_yyzz, \
                                     pa2pb_z_yz, pa2pb_z_yzzz, pa2pb_z_zz, pa2pb_z_zzzz, pa2pb_zz_y, pa2pb_zz_yyy, \
                                     pa2pb_zz_yyz, pa2pb_zz_yzz, pa2pb_zz_z, pa2pb_zz_zzz, pa2pb_zzz_yy, pa2pb_zzz_yyyy, \
                                     pa2pb_zzz_yyyz, pa2pb_zzz_yyzz, pa2pb_zzz_yz, pa2pb_zzz_yzzz, pa2pb_zzz_zz, \
                                     pa2pb_zzz_zzzz, pa_z, pa_zzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, s_0_0, \
                                     t_zzz_yyyy, t_zzz_yyyz, t_zzz_yyzz, t_zzz_yzzz, t_zzz_zzzz: VLX_ALIGN)
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

                t_zzz_yyyy[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pa_zzz[j] * fl2_fx + 4.5 * pa2pb_z_yy[j] * fl2_fx + 3.0 * pa2pb_zzz_yy[j] * fl1_fx + 1.5 * pa2pb_z_yyyy[j] * fl1_fx + pa2pb_zzz_yyyy[j]);

                t_zzz_yyyy[j] += fl_r_0_0 * (-4.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl1_fz * fl3_fx + 7.5 * pa_zzz[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_zzz_yy[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_z_yy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_zzz_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_yyyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_yyyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_yyyy[j] * fl1_fz);

                t_zzz_yyyz[j] = fl_s_0_0 * (1.125 * pb_y[j] * fl3_fx + 2.25 * pa2pb_zz_y[j] * fl2_fx + 2.25 * pa2pb_z_yz[j] * fl2_fx + 0.75 * pb_yyy[j] * fl2_fx + 1.5 * pa2pb_zzz_yz[j] * fl1_fx + 1.5 * pa2pb_zz_yyy[j] * fl1_fx + 1.5 * pa2pb_z_yyyz[j] * fl1_fx + pa2pb_zzz_yyyz[j]);

                t_zzz_yyyz[j] += fl_r_0_0 * (-2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pb_y[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_zz_y[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zzz_yz[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_z_yz[j] * fl1_fz * fl2_fx + 7.5 * pb_yyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_zzz_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zz_yyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_yyyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_yyyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_yyyz[j] * fl1_fz);

                t_zzz_yyzz[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pb_z[j] * fl3_fx + 0.25 * pa_zzz[j] * fl2_fx + 1.5 * pa2pb_zz_z[j] * fl2_fx + 2.25 * pa2pb_z_yy[j] * fl2_fx + 0.75 * pa2pb_z_zz[j] * fl2_fx + 1.5 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_zzz_yy[j] * fl1_fx + 0.5 * pa2pb_zzz_zz[j] * fl1_fx + 3.0 * pa2pb_zz_yyz[j] * fl1_fx + 1.5 * pa2pb_z_yyzz[j] * fl1_fx + pa2pb_zzz_yyzz[j]);

                t_zzz_yyzz[j] += fl_r_0_0 * (-3.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_z[j] * fl3_fx * fl1_fz + 2.5 * pa_zzz[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_zz_z[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_zzz_zz[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_zz[j] * fl1_fz * fl2_fx + 15.0 * pb_yyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzz_yy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_zz[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_zz_yyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_yyzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_yyzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_yyzz[j] * fl1_fz);

                t_zzz_yzzz[j] = fl_s_0_0 * (1.875 * pb_y[j] * fl3_fx + 2.25 * pa2pb_zz_y[j] * fl2_fx + 6.75 * pa2pb_z_yz[j] * fl2_fx + 2.25 * pb_yzz[j] * fl2_fx + 1.5 * pa2pb_zzz_yz[j] * fl1_fx + 4.5 * pa2pb_zz_yzz[j] * fl1_fx + 1.5 * pa2pb_z_yzzz[j] * fl1_fx + pa2pb_zzz_yzzz[j]);

                t_zzz_yzzz[j] += fl_r_0_0 * (15.0 * pb_y[j] * fl3_fx * fl1_fz - 2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa2pb_zz_y[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zzz_yz[j] * fl1_fz * fl1_fgb + 22.5 * pb_yzz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_zzz_yz[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_zz_yzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_yzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_yzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_yzzz[j] * fl1_fz);

                t_zzz_zzzz[j] = fl_s_0_0 * (5.625 * pa_z[j] * fl3_fx + 7.5 * pb_z[j] * fl3_fx + 0.75 * pa_zzz[j] * fl2_fx + 9.0 * pa2pb_zz_z[j] * fl2_fx + 13.5 * pa2pb_z_zz[j] * fl2_fx + 3.0 * pb_zzz[j] * fl2_fx + 3.0 * pa2pb_zzz_zz[j] * fl1_fx + 6.0 * pa2pb_zz_zzz[j] * fl1_fx + 1.5 * pa2pb_z_zzzz[j] * fl1_fx + pa2pb_zzz_zzzz[j]);

                t_zzz_zzzz[j] += fl_r_0_0 * (-13.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb + 45.0 * pa_z[j] * fl3_fx * fl1_fz + 60.0 * pb_z[j] * fl3_fx * fl1_fz - 2.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_zzz[j] * fl1_fz * fl2_fx + 90.0 * pa2pb_zz_z[j] * fl1_fz * fl2_fx + 135.0 * pa2pb_z_zz[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zzz_zz[j] * fl1_fz * fl1_fgb + 30.0 * pb_zzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zzz_zz[j] * fl1_fz * fl1_fx + 72.0 * pa2pb_zz_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_z_zzzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_z_zzzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }


} // kinrecfunc namespace

