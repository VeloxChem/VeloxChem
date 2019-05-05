//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForGG.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForGG(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForGG_0_5(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                               braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_5_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_10_15(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_15_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_20_25(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_25_30(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_30_35(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_35_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_40_45(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_45_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_50_55(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_55_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_60_65(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_65_70(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_70_75(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_75_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_80_85(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_85_90(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_90_95(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_95_100(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                  braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_100_105(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_105_110(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_110_115(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_115_120(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_120_125(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_125_130(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_130_135(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_135_140(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_140_145(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_145_150(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_150_155(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_155_160(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_160_165(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_165_170(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_170_175(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_175_180(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_180_185(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_185_190(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_190_195(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_195_200(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_200_205(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_205_210(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_210_215(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_215_220(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGG_220_225(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compKineticEnergyForGG_0_5(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xxxx = paDistances.data(34 * idx + 19);

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

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_xxxx = pa2pbDistances.data(1156 * idx + 121);

            auto pa2pb_xx_xxxy = pa2pbDistances.data(1156 * idx + 122);

            auto pa2pb_xx_xxxz = pa2pbDistances.data(1156 * idx + 123);

            auto pa2pb_xx_xxyy = pa2pbDistances.data(1156 * idx + 124);

            auto pa2pb_xx_xxyz = pa2pbDistances.data(1156 * idx + 125);

            auto pa2pb_xxx_x = pa2pbDistances.data(1156 * idx + 306);

            auto pa2pb_xxx_y = pa2pbDistances.data(1156 * idx + 307);

            auto pa2pb_xxx_z = pa2pbDistances.data(1156 * idx + 308);

            auto pa2pb_xxx_xxx = pa2pbDistances.data(1156 * idx + 315);

            auto pa2pb_xxx_xxy = pa2pbDistances.data(1156 * idx + 316);

            auto pa2pb_xxx_xxz = pa2pbDistances.data(1156 * idx + 317);

            auto pa2pb_xxx_xyy = pa2pbDistances.data(1156 * idx + 318);

            auto pa2pb_xxx_xyz = pa2pbDistances.data(1156 * idx + 319);

            auto pa2pb_xxxx_xx = pa2pbDistances.data(1156 * idx + 649);

            auto pa2pb_xxxx_xy = pa2pbDistances.data(1156 * idx + 650);

            auto pa2pb_xxxx_xz = pa2pbDistances.data(1156 * idx + 651);

            auto pa2pb_xxxx_yy = pa2pbDistances.data(1156 * idx + 652);

            auto pa2pb_xxxx_yz = pa2pbDistances.data(1156 * idx + 653);

            auto pa2pb_xxxx_xxxx = pa2pbDistances.data(1156 * idx + 665);

            auto pa2pb_xxxx_xxxy = pa2pbDistances.data(1156 * idx + 666);

            auto pa2pb_xxxx_xxxz = pa2pbDistances.data(1156 * idx + 667);

            auto pa2pb_xxxx_xxyy = pa2pbDistances.data(1156 * idx + 668);

            auto pa2pb_xxxx_xxyz = pa2pbDistances.data(1156 * idx + 669);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxx_xxxx = primBuffer.data(225 * idx);

            auto t_xxxx_xxxy = primBuffer.data(225 * idx + 1);

            auto t_xxxx_xxxz = primBuffer.data(225 * idx + 2);

            auto t_xxxx_xxyy = primBuffer.data(225 * idx + 3);

            auto t_xxxx_xxyz = primBuffer.data(225 * idx + 4);

            // Batch of Integrals (0,5)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xxxx, \
                                     pa2pb_xx_xxxy, pa2pb_xx_xxxz, pa2pb_xx_xxyy, pa2pb_xx_xxyz, pa2pb_xx_xy, \
                                     pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xxx_x, pa2pb_xxx_xxx, pa2pb_xxx_xxy, \
                                     pa2pb_xxx_xxz, pa2pb_xxx_xyy, pa2pb_xxx_xyz, pa2pb_xxx_y, pa2pb_xxx_z, \
                                     pa2pb_xxxx_xx, pa2pb_xxxx_xxxx, pa2pb_xxxx_xxxy, pa2pb_xxxx_xxxz, pa2pb_xxxx_xxyy, \
                                     pa2pb_xxxx_xxyz, pa2pb_xxxx_xy, pa2pb_xxxx_xz, pa2pb_xxxx_yy, pa2pb_xxxx_yz, pa_xx, \
                                     pa_xxxx, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_yy, pb_yz, \
                                     r_0_0, s_0_0, t_xxxx_xxxx, t_xxxx_xxxy, t_xxxx_xxxz, t_xxxx_xxyy, t_xxxx_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxx_xxxx[j] = fl_s_0_0 * (6.5625 * fl4_fx + 11.25 * pa_xx[j] * fl3_fx + 30.0 * pa2pb_x_x[j] * fl3_fx + 11.25 * pb_xx[j] * fl3_fx + 0.75 * pa_xxxx[j] * fl2_fx + 12.0 * pa2pb_xxx_x[j] * fl2_fx + 27.0 * pa2pb_xx_xx[j] * fl2_fx + 12.0 * pa2pb_x_xxx[j] * fl2_fx + 3.0 * pa2pb_xxxx_xx[j] * fl1_fx + 8.0 * pa2pb_xxx_xxx[j] * fl1_fx + 0.75 * pb_xxxx[j] * fl2_fx + 3.0 * pa2pb_xx_xxxx[j] * fl1_fx + pa2pb_xxxx_xxxx[j]);

                t_xxxx_xxxx[j] += fl_r_0_0 * (52.5 * fl4_fx * fl1_fz - 11.25 * fl3_fx * fl1_fz * fl1_fgb - 11.25 * fl3_fx * fl1_fz * fl1_fga - 27.0 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 112.5 * pa_xx[j] * fl3_fx * fl1_fz + 300.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 4.5 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 36.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 36.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 27.0 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xxxx[j] * fl1_fx * fl1_fz * fl1_fgb - 24.0 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb + 112.5 * pb_xx[j] * fl3_fx * fl1_fz + 9.0 * pa_xxxx[j] * fl1_fz * fl2_fx + 144.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 324.0 * pa2pb_xx_xx[j] * fl2_fx * fl1_fz - 4.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 24.0 * pa2pb_x_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xxxx_xx[j] * fl1_fz * fl1_fgb + 144.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxxx_xx[j] * fl1_fz * fl1_fx + 112.0 * pa2pb_xxx_xxx[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxxx[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xxxx[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xxxx[j] * fl1_fz);

                t_xxxx_xxxy[j] = fl_s_0_0 * (7.5 * pa2pb_x_y[j] * fl3_fx + 5.625 * pb_xy[j] * fl3_fx + 3.0 * pa2pb_xxx_y[j] * fl2_fx + 13.5 * pa2pb_xx_xy[j] * fl2_fx + 9.0 * pa2pb_x_xxy[j] * fl2_fx + 1.5 * pa2pb_xxxx_xy[j] * fl1_fx + 6.0 * pa2pb_xxx_xxy[j] * fl1_fx + 0.75 * pb_xxxy[j] * fl2_fx + 3.0 * pa2pb_xx_xxxy[j] * fl1_fx + pa2pb_xxxx_xxxy[j]);

                t_xxxx_xxxy[j] += fl_r_0_0 * (75.0 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 9.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 56.25 * pb_xy[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_xx_xy[j] * fl2_fx * fl1_fz - 2.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_x_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxx_xy[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxx_xy[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_xxx_xxy[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxxy[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xxxy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xxxy[j] * fl1_fz);

                t_xxxx_xxxz[j] = fl_s_0_0 * (7.5 * pa2pb_x_z[j] * fl3_fx + 5.625 * pb_xz[j] * fl3_fx + 3.0 * pa2pb_xxx_z[j] * fl2_fx + 13.5 * pa2pb_xx_xz[j] * fl2_fx + 9.0 * pa2pb_x_xxz[j] * fl2_fx + 1.5 * pa2pb_xxxx_xz[j] * fl1_fx + 6.0 * pa2pb_xxx_xxz[j] * fl1_fx + 0.75 * pb_xxxz[j] * fl2_fx + 3.0 * pa2pb_xx_xxxz[j] * fl1_fx + pa2pb_xxxx_xxxz[j]);

                t_xxxx_xxxz[j] += fl_r_0_0 * (75.0 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 9.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 56.25 * pb_xz[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_xx_xz[j] * fl2_fx * fl1_fz - 2.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_x_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxx_xz[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxx_xz[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_xxx_xxz[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxxz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xxxz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xxxz[j] * fl1_fz);

                t_xxxx_xxyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_xx[j] * fl3_fx + 3.0 * pa2pb_x_x[j] * fl3_fx + 1.875 * pb_yy[j] * fl3_fx + 0.25 * pa_xxxx[j] * fl2_fx + 2.0 * pa2pb_xxx_x[j] * fl2_fx + 4.5 * pa2pb_xx_yy[j] * fl2_fx + 0.375 * pb_xx[j] * fl3_fx + 1.5 * pa2pb_xx_xx[j] * fl2_fx + 6.0 * pa2pb_x_xyy[j] * fl2_fx + 0.5 * pa2pb_xxxx_xx[j] * fl1_fx + 0.5 * pa2pb_xxxx_yy[j] * fl1_fx + 4.0 * pa2pb_xxx_xyy[j] * fl1_fx + 0.75 * pb_xxyy[j] * fl2_fx + 3.0 * pa2pb_xx_xxyy[j] * fl1_fx + pa2pb_xxxx_xxyy[j]);

                t_xxxx_xxyy[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 6.0 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 22.5 * pa_xx[j] * fl3_fx * fl1_fz - 1.5 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - pa_xxxx[j] * fl1_fx * fl1_fz * fl1_fgb - 4.0 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 18.75 * pb_yy[j] * fl3_fx * fl1_fz + 3.0 * pa_xxxx[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_xx_yy[j] * fl2_fx * fl1_fz - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_x_xyy[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_xx[j] * fl3_fx * fl1_fz - pa2pb_xxxx_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxxx_yy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_xx[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxx_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxx_yy[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xxx_xyy[j] * fl1_fz * fl1_fx - 3.0 * pb_xxyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxyy[j] * fl1_fz * fl1_fga + 9.0 * pb_xxyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xxyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xxyy[j] * fl1_fz);

                t_xxxx_xxyz[j] = fl_s_0_0 * (1.875 * pb_yz[j] * fl3_fx + 4.5 * pa2pb_xx_yz[j] * fl2_fx + 6.0 * pa2pb_x_xyz[j] * fl2_fx + 0.5 * pa2pb_xxxx_yz[j] * fl1_fx + 4.0 * pa2pb_xxx_xyz[j] * fl1_fx + 0.75 * pb_xxyz[j] * fl2_fx + 3.0 * pa2pb_xx_xxyz[j] * fl1_fx + pa2pb_xxxx_xxyz[j]);

                t_xxxx_xxyz[j] += fl_r_0_0 * (-4.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga + 18.75 * pb_yz[j] * fl3_fx * fl1_fz + 54.0 * pa2pb_xx_yz[j] * fl2_fx * fl1_fz - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_x_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxx_yz[j] * fl1_fz * fl1_fgb + 72.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxx_yz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xxx_xyz[j] * fl1_fz * fl1_fx - 3.0 * pb_xxyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxyz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xxyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_5_10(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_xx_xxzz = pa2pbDistances.data(1156 * idx + 126);

            auto pa2pb_xx_xyyy = pa2pbDistances.data(1156 * idx + 127);

            auto pa2pb_xx_xyyz = pa2pbDistances.data(1156 * idx + 128);

            auto pa2pb_xx_xyzz = pa2pbDistances.data(1156 * idx + 129);

            auto pa2pb_xx_xzzz = pa2pbDistances.data(1156 * idx + 130);

            auto pa2pb_xxx_x = pa2pbDistances.data(1156 * idx + 306);

            auto pa2pb_xxx_y = pa2pbDistances.data(1156 * idx + 307);

            auto pa2pb_xxx_z = pa2pbDistances.data(1156 * idx + 308);

            auto pa2pb_xxx_xzz = pa2pbDistances.data(1156 * idx + 320);

            auto pa2pb_xxx_yyy = pa2pbDistances.data(1156 * idx + 321);

            auto pa2pb_xxx_yyz = pa2pbDistances.data(1156 * idx + 322);

            auto pa2pb_xxx_yzz = pa2pbDistances.data(1156 * idx + 323);

            auto pa2pb_xxx_zzz = pa2pbDistances.data(1156 * idx + 324);

            auto pa2pb_xxxx_xx = pa2pbDistances.data(1156 * idx + 649);

            auto pa2pb_xxxx_xy = pa2pbDistances.data(1156 * idx + 650);

            auto pa2pb_xxxx_xz = pa2pbDistances.data(1156 * idx + 651);

            auto pa2pb_xxxx_zz = pa2pbDistances.data(1156 * idx + 654);

            auto pa2pb_xxxx_xxzz = pa2pbDistances.data(1156 * idx + 670);

            auto pa2pb_xxxx_xyyy = pa2pbDistances.data(1156 * idx + 671);

            auto pa2pb_xxxx_xyyz = pa2pbDistances.data(1156 * idx + 672);

            auto pa2pb_xxxx_xyzz = pa2pbDistances.data(1156 * idx + 673);

            auto pa2pb_xxxx_xzzz = pa2pbDistances.data(1156 * idx + 674);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxx_xxzz = primBuffer.data(225 * idx + 5);

            auto t_xxxx_xyyy = primBuffer.data(225 * idx + 6);

            auto t_xxxx_xyyz = primBuffer.data(225 * idx + 7);

            auto t_xxxx_xyzz = primBuffer.data(225 * idx + 8);

            auto t_xxxx_xzzz = primBuffer.data(225 * idx + 9);

            // Batch of Integrals (5,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_yyy, \
                                     pa2pb_x_yyz, pa2pb_x_yzz, pa2pb_x_z, pa2pb_x_zzz, pa2pb_xx_xx, pa2pb_xx_xxzz, \
                                     pa2pb_xx_xy, pa2pb_xx_xyyy, pa2pb_xx_xyyz, pa2pb_xx_xyzz, pa2pb_xx_xz, \
                                     pa2pb_xx_xzzz, pa2pb_xx_zz, pa2pb_xxx_x, pa2pb_xxx_xzz, pa2pb_xxx_y, pa2pb_xxx_yyy, \
                                     pa2pb_xxx_yyz, pa2pb_xxx_yzz, pa2pb_xxx_z, pa2pb_xxx_zzz, pa2pb_xxxx_xx, \
                                     pa2pb_xxxx_xxzz, pa2pb_xxxx_xy, pa2pb_xxxx_xyyy, pa2pb_xxxx_xyyz, pa2pb_xxxx_xyzz, \
                                     pa2pb_xxxx_xz, pa2pb_xxxx_xzzz, pa2pb_xxxx_zz, pa_xx, pa_xxxx, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, pb_xzzz, pb_zz, r_0_0, s_0_0, t_xxxx_xxzz, \
                                     t_xxxx_xyyy, t_xxxx_xyyz, t_xxxx_xyzz, t_xxxx_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxx_xxzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_xx[j] * fl3_fx + 3.0 * pa2pb_x_x[j] * fl3_fx + 1.875 * pb_zz[j] * fl3_fx + 0.25 * pa_xxxx[j] * fl2_fx + 2.0 * pa2pb_xxx_x[j] * fl2_fx + 4.5 * pa2pb_xx_zz[j] * fl2_fx + 0.375 * pb_xx[j] * fl3_fx + 1.5 * pa2pb_xx_xx[j] * fl2_fx + 6.0 * pa2pb_x_xzz[j] * fl2_fx + 0.5 * pa2pb_xxxx_xx[j] * fl1_fx + 0.5 * pa2pb_xxxx_zz[j] * fl1_fx + 4.0 * pa2pb_xxx_xzz[j] * fl1_fx + 0.75 * pb_xxzz[j] * fl2_fx + 3.0 * pa2pb_xx_xxzz[j] * fl1_fx + pa2pb_xxxx_xxzz[j]);

                t_xxxx_xxzz[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 6.0 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 22.5 * pa_xx[j] * fl3_fx * fl1_fz - 1.5 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - pa_xxxx[j] * fl1_fx * fl1_fz * fl1_fgb - 4.0 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 18.75 * pb_zz[j] * fl3_fx * fl1_fz + 3.0 * pa_xxxx[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_xx_zz[j] * fl2_fx * fl1_fz - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_x_xzz[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_xx[j] * fl3_fx * fl1_fz - pa2pb_xxxx_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxxx_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_xx[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxx_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxx_zz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xxx_xzz[j] * fl1_fz * fl1_fx - 3.0 * pb_xxzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xxzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xxzz[j] * fl1_fz);

                t_xxxx_xyyy[j] = fl_s_0_0 * (4.5 * pa2pb_x_y[j] * fl3_fx + 3.0 * pa2pb_xxx_y[j] * fl2_fx + 1.125 * pb_xy[j] * fl3_fx + 4.5 * pa2pb_xx_xy[j] * fl2_fx + 3.0 * pa2pb_x_yyy[j] * fl2_fx + 1.5 * pa2pb_xxxx_xy[j] * fl1_fx + 2.0 * pa2pb_xxx_yyy[j] * fl1_fx + 0.75 * pb_xyyy[j] * fl2_fx + 3.0 * pa2pb_xx_xyyy[j] * fl1_fx + pa2pb_xxxx_xyyy[j]);

                t_xxxx_xyyy[j] += fl_r_0_0 * (-9.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx - 2.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_yyy[j] * fl1_fx * fl1_fz * fl1_fga + 11.25 * pb_xy[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_xxxx_xy[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xx_xy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxx_xy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxx_yyy[j] * fl1_fz * fl1_fx - 3.0 * pb_xyyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xyyy[j] * fl1_fz * fl1_fga + 9.0 * pb_xyyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xyyy[j] * fl1_fz);

                t_xxxx_xyyz[j] = fl_s_0_0 * (1.5 * pa2pb_x_z[j] * fl3_fx + pa2pb_xxx_z[j] * fl2_fx + 0.375 * pb_xz[j] * fl3_fx + 1.5 * pa2pb_xx_xz[j] * fl2_fx + 3.0 * pa2pb_x_yyz[j] * fl2_fx + 0.5 * pa2pb_xxxx_xz[j] * fl1_fx + 2.0 * pa2pb_xxx_yyz[j] * fl1_fx + 0.75 * pb_xyyz[j] * fl2_fx + 3.0 * pa2pb_xx_xyyz[j] * fl1_fx + pa2pb_xxxx_xyyz[j]);

                t_xxxx_xyyz[j] += fl_r_0_0 * (-3.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 12.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_yyz[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_xz[j] * fl3_fx * fl1_fz - pa2pb_xxxx_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_xz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxx_xz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxx_yyz[j] * fl1_fz * fl1_fx - 3.0 * pb_xyyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xyyz[j] * fl1_fz * fl1_fga + 9.0 * pb_xyyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xyyz[j] * fl1_fz);

                t_xxxx_xyzz[j] = fl_s_0_0 * (1.5 * pa2pb_x_y[j] * fl3_fx + pa2pb_xxx_y[j] * fl2_fx + 0.375 * pb_xy[j] * fl3_fx + 1.5 * pa2pb_xx_xy[j] * fl2_fx + 3.0 * pa2pb_x_yzz[j] * fl2_fx + 0.5 * pa2pb_xxxx_xy[j] * fl1_fx + 2.0 * pa2pb_xxx_yzz[j] * fl1_fx + 0.75 * pb_xyzz[j] * fl2_fx + 3.0 * pa2pb_xx_xyzz[j] * fl1_fx + pa2pb_xxxx_xyzz[j]);

                t_xxxx_xyzz[j] += fl_r_0_0 * (-3.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 12.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_yzz[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_xy[j] * fl3_fx * fl1_fz - pa2pb_xxxx_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_xy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxx_xy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxx_yzz[j] * fl1_fz * fl1_fx - 3.0 * pb_xyzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xyzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xyzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xyzz[j] * fl1_fz);

                t_xxxx_xzzz[j] = fl_s_0_0 * (4.5 * pa2pb_x_z[j] * fl3_fx + 3.0 * pa2pb_xxx_z[j] * fl2_fx + 1.125 * pb_xz[j] * fl3_fx + 4.5 * pa2pb_xx_xz[j] * fl2_fx + 3.0 * pa2pb_x_zzz[j] * fl2_fx + 1.5 * pa2pb_xxxx_xz[j] * fl1_fx + 2.0 * pa2pb_xxx_zzz[j] * fl1_fx + 0.75 * pb_xzzz[j] * fl2_fx + 3.0 * pa2pb_xx_xzzz[j] * fl1_fx + pa2pb_xxxx_xzzz[j]);

                t_xxxx_xzzz[j] += fl_r_0_0 * (-9.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx - 2.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_zzz[j] * fl1_fx * fl1_fz * fl1_fga + 11.25 * pb_xz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_xxxx_xz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xx_xz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxx_xz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxx_zzz[j] * fl1_fz * fl1_fx - 3.0 * pb_xzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_xzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_10_15(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_xx_yyyy = pa2pbDistances.data(1156 * idx + 131);

            auto pa2pb_xx_yyyz = pa2pbDistances.data(1156 * idx + 132);

            auto pa2pb_xx_yyzz = pa2pbDistances.data(1156 * idx + 133);

            auto pa2pb_xx_yzzz = pa2pbDistances.data(1156 * idx + 134);

            auto pa2pb_xx_zzzz = pa2pbDistances.data(1156 * idx + 135);

            auto pa2pb_xxxx_yy = pa2pbDistances.data(1156 * idx + 652);

            auto pa2pb_xxxx_yz = pa2pbDistances.data(1156 * idx + 653);

            auto pa2pb_xxxx_zz = pa2pbDistances.data(1156 * idx + 654);

            auto pa2pb_xxxx_yyyy = pa2pbDistances.data(1156 * idx + 675);

            auto pa2pb_xxxx_yyyz = pa2pbDistances.data(1156 * idx + 676);

            auto pa2pb_xxxx_yyzz = pa2pbDistances.data(1156 * idx + 677);

            auto pa2pb_xxxx_yzzz = pa2pbDistances.data(1156 * idx + 678);

            auto pa2pb_xxxx_zzzz = pa2pbDistances.data(1156 * idx + 679);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxx_yyyy = primBuffer.data(225 * idx + 10);

            auto t_xxxx_yyyz = primBuffer.data(225 * idx + 11);

            auto t_xxxx_yyzz = primBuffer.data(225 * idx + 12);

            auto t_xxxx_yzzz = primBuffer.data(225 * idx + 13);

            auto t_xxxx_zzzz = primBuffer.data(225 * idx + 14);

            // Batch of Integrals (10,15)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_xx_yy, pa2pb_xx_yyyy, pa2pb_xx_yyyz, pa2pb_xx_yyzz, \
                                     pa2pb_xx_yz, pa2pb_xx_yzzz, pa2pb_xx_zz, pa2pb_xx_zzzz, pa2pb_xxxx_yy, \
                                     pa2pb_xxxx_yyyy, pa2pb_xxxx_yyyz, pa2pb_xxxx_yyzz, pa2pb_xxxx_yz, pa2pb_xxxx_yzzz, \
                                     pa2pb_xxxx_zz, pa2pb_xxxx_zzzz, pa_xx, pa_xxxx, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, \
                                     pb_yzzz, pb_zz, pb_zzzz, r_0_0, s_0_0, t_xxxx_yyyy, t_xxxx_yyyz, t_xxxx_yyzz, \
                                     t_xxxx_yzzz, t_xxxx_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxx_yyyy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_xx[j] * fl3_fx + 0.75 * pa_xxxx[j] * fl2_fx + 2.25 * pb_yy[j] * fl3_fx + 9.0 * pa2pb_xx_yy[j] * fl2_fx + 3.0 * pa2pb_xxxx_yy[j] * fl1_fx + 0.75 * pb_yyyy[j] * fl2_fx + 3.0 * pa2pb_xx_yyyy[j] * fl1_fx + pa2pb_xxxx_yyyy[j]);

                t_xxxx_yyyy[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 9.0 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx + 4.5 * fl4_fx * fl1_fz - 3.0 * pa_xxxx[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_xx[j] * fl1_fz * fl3_fx + 9.0 * pa_xxxx[j] * fl1_fz * fl2_fx - 4.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 18.0 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx + 22.5 * pb_yy[j] * fl3_fx * fl1_fz - 6.0 * pa2pb_xxxx_yy[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_xx_yy[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_xxxx_yy[j] * fl1_fz * fl1_fx - 3.0 * pb_yyyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_yyyy[j] * fl1_fz * fl1_fga + 9.0 * pb_yyyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_yyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_yyyy[j] * fl1_fz);

                t_xxxx_yyyz[j] = fl_s_0_0 * (1.125 * pb_yz[j] * fl3_fx + 4.5 * pa2pb_xx_yz[j] * fl2_fx + 1.5 * pa2pb_xxxx_yz[j] * fl1_fx + 0.75 * pb_yyyz[j] * fl2_fx + 3.0 * pa2pb_xx_yyyz[j] * fl1_fx + pa2pb_xxxx_yyyz[j]);

                t_xxxx_yyyz[j] += fl_r_0_0 * (-2.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx + 11.25 * pb_yz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_xxxx_yz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xx_yz[j] * fl1_fz * fl2_fx + 21.0 * pa2pb_xxxx_yz[j] * fl1_fz * fl1_fx - 3.0 * pb_yyyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_yyyz[j] * fl1_fz * fl1_fga + 9.0 * pb_yyyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_yyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_yyyz[j] * fl1_fz);

                t_xxxx_yyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.75 * pa_xx[j] * fl3_fx + 0.25 * pa_xxxx[j] * fl2_fx + 0.375 * pb_yy[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 1.5 * pa2pb_xx_yy[j] * fl2_fx + 1.5 * pa2pb_xx_zz[j] * fl2_fx + 0.5 * pa2pb_xxxx_yy[j] * fl1_fx + 0.5 * pa2pb_xxxx_zz[j] * fl1_fx + 0.75 * pb_yyzz[j] * fl2_fx + 3.0 * pa2pb_xx_yyzz[j] * fl1_fx + pa2pb_xxxx_yyzz[j]);

                t_xxxx_yyzz[j] += fl_r_0_0 * (-0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - 3.0 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx + 1.5 * fl4_fx * fl1_fz - pa_xxxx[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_xx[j] * fl1_fz * fl3_fx + 3.0 * pa_xxxx[j] * fl1_fz * fl2_fx - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz - pa2pb_xxxx_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxxx_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_yy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xx_zz[j] * fl1_fz * fl2_fx + 7.0 * pa2pb_xxxx_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxx_zz[j] * fl1_fz * fl1_fx - 3.0 * pb_yyzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_yyzz[j] * fl1_fz * fl1_fga + 9.0 * pb_yyzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_yyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_yyzz[j] * fl1_fz);

                t_xxxx_yzzz[j] = fl_s_0_0 * (1.125 * pb_yz[j] * fl3_fx + 4.5 * pa2pb_xx_yz[j] * fl2_fx + 1.5 * pa2pb_xxxx_yz[j] * fl1_fx + 0.75 * pb_yzzz[j] * fl2_fx + 3.0 * pa2pb_xx_yzzz[j] * fl1_fx + pa2pb_xxxx_yzzz[j]);

                t_xxxx_yzzz[j] += fl_r_0_0 * (-2.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx + 11.25 * pb_yz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_xxxx_yz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xx_yz[j] * fl1_fz * fl2_fx + 21.0 * pa2pb_xxxx_yz[j] * fl1_fz * fl1_fx - 3.0 * pb_yzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_yzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_yzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_yzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_yzzz[j] * fl1_fz);

                t_xxxx_zzzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_xx[j] * fl3_fx + 0.75 * pa_xxxx[j] * fl2_fx + 2.25 * pb_zz[j] * fl3_fx + 9.0 * pa2pb_xx_zz[j] * fl2_fx + 3.0 * pa2pb_xxxx_zz[j] * fl1_fx + 0.75 * pb_zzzz[j] * fl2_fx + 3.0 * pa2pb_xx_zzzz[j] * fl1_fx + pa2pb_xxxx_zzzz[j]);

                t_xxxx_zzzz[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 9.0 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx + 4.5 * fl4_fx * fl1_fz - 3.0 * pa_xxxx[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_xx[j] * fl1_fz * fl3_fx + 9.0 * pa_xxxx[j] * fl1_fz * fl2_fx - 4.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 18.0 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx + 22.5 * pb_zz[j] * fl3_fx * fl1_fz - 6.0 * pa2pb_xxxx_zz[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_xx_zz[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_xxxx_zz[j] * fl1_fz * fl1_fx - 3.0 * pb_zzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_zzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_zzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xx_zzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xxxx_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_15_20(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_xxxx = pa2pbDistances.data(1156 * idx + 155);

            auto pa2pb_xy_xxxy = pa2pbDistances.data(1156 * idx + 156);

            auto pa2pb_xy_xxxz = pa2pbDistances.data(1156 * idx + 157);

            auto pa2pb_xy_xxyy = pa2pbDistances.data(1156 * idx + 158);

            auto pa2pb_xy_xxyz = pa2pbDistances.data(1156 * idx + 159);

            auto pa2pb_xxx_x = pa2pbDistances.data(1156 * idx + 306);

            auto pa2pb_xxx_y = pa2pbDistances.data(1156 * idx + 307);

            auto pa2pb_xxx_z = pa2pbDistances.data(1156 * idx + 308);

            auto pa2pb_xxx_xxx = pa2pbDistances.data(1156 * idx + 315);

            auto pa2pb_xxx_xxy = pa2pbDistances.data(1156 * idx + 316);

            auto pa2pb_xxx_xxz = pa2pbDistances.data(1156 * idx + 317);

            auto pa2pb_xxy_x = pa2pbDistances.data(1156 * idx + 340);

            auto pa2pb_xxy_y = pa2pbDistances.data(1156 * idx + 341);

            auto pa2pb_xxy_z = pa2pbDistances.data(1156 * idx + 342);

            auto pa2pb_xxy_xxx = pa2pbDistances.data(1156 * idx + 349);

            auto pa2pb_xxy_xxy = pa2pbDistances.data(1156 * idx + 350);

            auto pa2pb_xxy_xxz = pa2pbDistances.data(1156 * idx + 351);

            auto pa2pb_xxy_xyy = pa2pbDistances.data(1156 * idx + 352);

            auto pa2pb_xxy_xyz = pa2pbDistances.data(1156 * idx + 353);

            auto pa2pb_xxxy_xx = pa2pbDistances.data(1156 * idx + 683);

            auto pa2pb_xxxy_xy = pa2pbDistances.data(1156 * idx + 684);

            auto pa2pb_xxxy_xz = pa2pbDistances.data(1156 * idx + 685);

            auto pa2pb_xxxy_yy = pa2pbDistances.data(1156 * idx + 686);

            auto pa2pb_xxxy_yz = pa2pbDistances.data(1156 * idx + 687);

            auto pa2pb_xxxy_xxxx = pa2pbDistances.data(1156 * idx + 699);

            auto pa2pb_xxxy_xxxy = pa2pbDistances.data(1156 * idx + 700);

            auto pa2pb_xxxy_xxxz = pa2pbDistances.data(1156 * idx + 701);

            auto pa2pb_xxxy_xxyy = pa2pbDistances.data(1156 * idx + 702);

            auto pa2pb_xxxy_xxyz = pa2pbDistances.data(1156 * idx + 703);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxy_xxxx = primBuffer.data(225 * idx + 15);

            auto t_xxxy_xxxy = primBuffer.data(225 * idx + 16);

            auto t_xxxy_xxxz = primBuffer.data(225 * idx + 17);

            auto t_xxxy_xxyy = primBuffer.data(225 * idx + 18);

            auto t_xxxy_xxyz = primBuffer.data(225 * idx + 19);

            // Batch of Integrals (15,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xy, pa2pb_xx_xz, pa2pb_xxx_x, \
                                     pa2pb_xxx_xxx, pa2pb_xxx_xxy, pa2pb_xxx_xxz, pa2pb_xxx_y, pa2pb_xxx_z, \
                                     pa2pb_xxxy_xx, pa2pb_xxxy_xxxx, pa2pb_xxxy_xxxy, pa2pb_xxxy_xxxz, pa2pb_xxxy_xxyy, \
                                     pa2pb_xxxy_xxyz, pa2pb_xxxy_xy, pa2pb_xxxy_xz, pa2pb_xxxy_yy, pa2pb_xxxy_yz, \
                                     pa2pb_xxy_x, pa2pb_xxy_xxx, pa2pb_xxy_xxy, pa2pb_xxy_xxz, pa2pb_xxy_xyy, \
                                     pa2pb_xxy_xyz, pa2pb_xxy_y, pa2pb_xxy_z, pa2pb_xy_xx, pa2pb_xy_xxxx, pa2pb_xy_xxxy, \
                                     pa2pb_xy_xxxz, pa2pb_xy_xxyy, pa2pb_xy_xxyz, pa2pb_xy_xy, pa2pb_xy_xz, pa2pb_xy_yy, \
                                     pa2pb_xy_yz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, pa2pb_y_xyy, \
                                     pa2pb_y_xyz, pa2pb_y_y, pa2pb_y_z, pa_xx, pa_xxxy, pa_xy, pb_xx, pb_xy, pb_xz, r_0_0, s_0_0, \
                                     t_xxxy_xxxx, t_xxxy_xxxy, t_xxxy_xxxz, t_xxxy_xxyy, t_xxxy_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxy_xxxx[j] = fl_s_0_0 * (5.625 * pa_xy[j] * fl3_fx + 7.5 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa_xxxy[j] * fl2_fx + 9.0 * pa2pb_xxy_x[j] * fl2_fx + 13.5 * pa2pb_xy_xx[j] * fl2_fx + 3.0 * pa2pb_y_xxx[j] * fl2_fx + 3.0 * pa2pb_xxxy_xx[j] * fl1_fx + 6.0 * pa2pb_xxy_xxx[j] * fl1_fx + 1.5 * pa2pb_xy_xxxx[j] * fl1_fx + pa2pb_xxxy_xxxx[j]);

                t_xxxy_xxxx[j] += fl_r_0_0 * (-13.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 56.25 * pa_xy[j] * fl3_fx * fl1_fz + 75.0 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 2.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xxxy[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_xxxy[j] * fl1_fz * fl2_fx + 108.0 * pa2pb_xxy_x[j] * fl2_fx * fl1_fz + 162.0 * pa2pb_xy_xx[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xxxy_xx[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxxy_xx[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_xxy_xxx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxxx[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxxx[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xxxx[j] * fl1_fz);

                t_xxxy_xxxy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.125 * pa_xx[j] * fl3_fx + 3.375 * pa2pb_x_x[j] * fl3_fx + 1.875 * pa2pb_y_y[j] * fl3_fx + 1.125 * pb_xx[j] * fl3_fx + 0.75 * pa2pb_xxx_x[j] * fl2_fx + 2.25 * pa2pb_xxy_y[j] * fl2_fx + 2.25 * pa2pb_xx_xx[j] * fl2_fx + 6.75 * pa2pb_xy_xy[j] * fl2_fx + 0.75 * pa2pb_x_xxx[j] * fl2_fx + 2.25 * pa2pb_y_xxy[j] * fl2_fx + 1.5 * pa2pb_xxxy_xy[j] * fl1_fx + 0.5 * pa2pb_xxx_xxx[j] * fl1_fx + 4.5 * pa2pb_xxy_xxy[j] * fl1_fx + 1.5 * pa2pb_xy_xxxy[j] * fl1_fx + pa2pb_xxxy_xxxy[j]);

                t_xxxy_xxxy[j] += fl_r_0_0 * (7.5 * fl4_fx * fl1_fz - 1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xx[j] * fl3_fx * fl1_fz + 33.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 18.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pb_xx[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xxy_y[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xx_xx[j] * fl2_fx * fl1_fz + 81.0 * pa2pb_xy_xy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxy_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxy_xy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxx_xxx[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_xxy_xxy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxxy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxxy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xxxy[j] * fl1_fz);

                t_xxxy_xxxz[j] = fl_s_0_0 * (1.875 * pa2pb_y_z[j] * fl3_fx + 2.25 * pa2pb_xxy_z[j] * fl2_fx + 6.75 * pa2pb_xy_xz[j] * fl2_fx + 2.25 * pa2pb_y_xxz[j] * fl2_fx + 1.5 * pa2pb_xxxy_xz[j] * fl1_fx + 4.5 * pa2pb_xxy_xxz[j] * fl1_fx + 1.5 * pa2pb_xy_xxxz[j] * fl1_fx + pa2pb_xxxy_xxxz[j]);

                t_xxxy_xxxz[j] += fl_r_0_0 * (18.75 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 27.0 * pa2pb_xxy_z[j] * fl2_fx * fl1_fz + 81.0 * pa2pb_xy_xz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxy_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxy_xz[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_xxy_xxz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxxz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxxz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xxxz[j] * fl1_fz);

                t_xxxy_xxyy[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 2.25 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa2pb_y_x[j] * fl3_fx + 1.5 * pb_xy[j] * fl3_fx + 0.25 * pa_xxxy[j] * fl2_fx + 0.5 * pa2pb_xxx_y[j] * fl2_fx + 1.5 * pa2pb_xxy_x[j] * fl2_fx + 3.0 * pa2pb_xx_xy[j] * fl2_fx + 2.25 * pa2pb_xy_yy[j] * fl2_fx + 0.75 * pa2pb_xy_xx[j] * fl2_fx + 1.5 * pa2pb_x_xxy[j] * fl2_fx + 1.5 * pa2pb_y_xyy[j] * fl2_fx + 0.5 * pa2pb_xxxy_xx[j] * fl1_fx + 0.5 * pa2pb_xxxy_yy[j] * fl1_fx + pa2pb_xxx_xxy[j] * fl1_fx + 3.0 * pa2pb_xxy_xyy[j] * fl1_fx + 1.5 * pa2pb_xy_xxyy[j] * fl1_fx + pa2pb_xxxy_xxyy[j]);

                t_xxxy_xxyy[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xy[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - pa_xxxy[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 15.0 * pb_xy[j] * fl3_fx * fl1_fz + 3.0 * pa_xxxy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xxy_x[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_xy[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xy_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xxy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxxy_yy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxy_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxy_yy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xxy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxy_xyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xxyy[j] * fl1_fz);

                t_xxxy_xxyz[j] = fl_s_0_0 * (1.125 * pa2pb_x_z[j] * fl3_fx + 0.75 * pb_xz[j] * fl3_fx + 0.25 * pa2pb_xxx_z[j] * fl2_fx + 1.5 * pa2pb_xx_xz[j] * fl2_fx + 2.25 * pa2pb_xy_yz[j] * fl2_fx + 0.75 * pa2pb_x_xxz[j] * fl2_fx + 1.5 * pa2pb_y_xyz[j] * fl2_fx + 0.5 * pa2pb_xxxy_yz[j] * fl1_fx + 0.5 * pa2pb_xxx_xxz[j] * fl1_fx + 3.0 * pa2pb_xxy_xyz[j] * fl1_fx + 1.5 * pa2pb_xy_xxyz[j] * fl1_fx + pa2pb_xxxy_xxyz[j]);

                t_xxxy_xxyz[j] += fl_r_0_0 * (11.25 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_xz[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xx_xz[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xy_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xxz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxy_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxx_xxz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxy_xyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_20_25(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_xy_xxzz = pa2pbDistances.data(1156 * idx + 160);

            auto pa2pb_xy_xyyy = pa2pbDistances.data(1156 * idx + 161);

            auto pa2pb_xy_xyyz = pa2pbDistances.data(1156 * idx + 162);

            auto pa2pb_xy_xyzz = pa2pbDistances.data(1156 * idx + 163);

            auto pa2pb_xy_xzzz = pa2pbDistances.data(1156 * idx + 164);

            auto pa2pb_xxx_x = pa2pbDistances.data(1156 * idx + 306);

            auto pa2pb_xxx_xyy = pa2pbDistances.data(1156 * idx + 318);

            auto pa2pb_xxx_xyz = pa2pbDistances.data(1156 * idx + 319);

            auto pa2pb_xxx_xzz = pa2pbDistances.data(1156 * idx + 320);

            auto pa2pb_xxy_x = pa2pbDistances.data(1156 * idx + 340);

            auto pa2pb_xxy_y = pa2pbDistances.data(1156 * idx + 341);

            auto pa2pb_xxy_z = pa2pbDistances.data(1156 * idx + 342);

            auto pa2pb_xxy_xzz = pa2pbDistances.data(1156 * idx + 354);

            auto pa2pb_xxy_yyy = pa2pbDistances.data(1156 * idx + 355);

            auto pa2pb_xxy_yyz = pa2pbDistances.data(1156 * idx + 356);

            auto pa2pb_xxy_yzz = pa2pbDistances.data(1156 * idx + 357);

            auto pa2pb_xxy_zzz = pa2pbDistances.data(1156 * idx + 358);

            auto pa2pb_xxxy_xx = pa2pbDistances.data(1156 * idx + 683);

            auto pa2pb_xxxy_xy = pa2pbDistances.data(1156 * idx + 684);

            auto pa2pb_xxxy_xz = pa2pbDistances.data(1156 * idx + 685);

            auto pa2pb_xxxy_zz = pa2pbDistances.data(1156 * idx + 688);

            auto pa2pb_xxxy_xxzz = pa2pbDistances.data(1156 * idx + 704);

            auto pa2pb_xxxy_xyyy = pa2pbDistances.data(1156 * idx + 705);

            auto pa2pb_xxxy_xyyz = pa2pbDistances.data(1156 * idx + 706);

            auto pa2pb_xxxy_xyzz = pa2pbDistances.data(1156 * idx + 707);

            auto pa2pb_xxxy_xzzz = pa2pbDistances.data(1156 * idx + 708);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxy_xxzz = primBuffer.data(225 * idx + 20);

            auto t_xxxy_xyyy = primBuffer.data(225 * idx + 21);

            auto t_xxxy_xyyz = primBuffer.data(225 * idx + 22);

            auto t_xxxy_xyzz = primBuffer.data(225 * idx + 23);

            auto t_xxxy_xzzz = primBuffer.data(225 * idx + 24);

            // Batch of Integrals (20,25)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_xzz, \
                                     pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa2pb_xxx_x, pa2pb_xxx_xyy, pa2pb_xxx_xyz, \
                                     pa2pb_xxx_xzz, pa2pb_xxxy_xx, pa2pb_xxxy_xxzz, pa2pb_xxxy_xy, pa2pb_xxxy_xyyy, \
                                     pa2pb_xxxy_xyyz, pa2pb_xxxy_xyzz, pa2pb_xxxy_xz, pa2pb_xxxy_xzzz, pa2pb_xxxy_zz, \
                                     pa2pb_xxy_x, pa2pb_xxy_xzz, pa2pb_xxy_y, pa2pb_xxy_yyy, pa2pb_xxy_yyz, \
                                     pa2pb_xxy_yzz, pa2pb_xxy_z, pa2pb_xxy_zzz, pa2pb_xy_xx, pa2pb_xy_xxzz, pa2pb_xy_xy, \
                                     pa2pb_xy_xyyy, pa2pb_xy_xyyz, pa2pb_xy_xyzz, pa2pb_xy_xz, pa2pb_xy_xzzz, \
                                     pa2pb_xy_zz, pa2pb_y_x, pa2pb_y_xzz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, \
                                     pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, pa_xx, pa_xxxy, pa_xy, pb_yy, pb_yz, pb_zz, r_0_0, \
                                     s_0_0, t_xxxy_xxzz, t_xxxy_xyyy, t_xxxy_xyyz, t_xxxy_xyzz, t_xxxy_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxy_xxzz[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 0.75 * pa2pb_y_x[j] * fl3_fx + 0.25 * pa_xxxy[j] * fl2_fx + 1.5 * pa2pb_xxy_x[j] * fl2_fx + 2.25 * pa2pb_xy_zz[j] * fl2_fx + 0.75 * pa2pb_xy_xx[j] * fl2_fx + 1.5 * pa2pb_y_xzz[j] * fl2_fx + 0.5 * pa2pb_xxxy_xx[j] * fl1_fx + 0.5 * pa2pb_xxxy_zz[j] * fl1_fx + 3.0 * pa2pb_xxy_xzz[j] * fl1_fx + 1.5 * pa2pb_xy_xxzz[j] * fl1_fx + pa2pb_xxxy_xxzz[j]);

                t_xxxy_xxzz[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xy[j] * fl3_fx * fl1_fz - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xxxy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 3.0 * pa_xxxy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xxy_x[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xy_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxxy_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxy_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxy_zz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxy_xzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xxzz[j] * fl1_fz);

                t_xxxy_xyyy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_xx[j] * fl3_fx + 1.125 * pa2pb_x_x[j] * fl3_fx + 1.125 * pa2pb_y_y[j] * fl3_fx + 1.125 * pb_yy[j] * fl3_fx + 0.75 * pa2pb_xxx_x[j] * fl2_fx + 2.25 * pa2pb_xxy_y[j] * fl2_fx + 2.25 * pa2pb_xx_yy[j] * fl2_fx + 2.25 * pa2pb_xy_xy[j] * fl2_fx + 2.25 * pa2pb_x_xyy[j] * fl2_fx + 0.75 * pa2pb_y_yyy[j] * fl2_fx + 1.5 * pa2pb_xxxy_xy[j] * fl1_fx + 1.5 * pa2pb_xxx_xyy[j] * fl1_fx + 1.5 * pa2pb_xxy_yyy[j] * fl1_fx + 1.5 * pa2pb_xy_xyyy[j] * fl1_fx + pa2pb_xxxy_xyyy[j]);

                t_xxxy_xyyy[j] += fl_r_0_0 * (-1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 4.5 * fl4_fx * fl1_fz + 11.25 * pa_xx[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 11.25 * pb_yy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xxy_y[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xx_yy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_xyy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxy_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xy_xy[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxy_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxx_xyy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxy_yyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xyyy[j] * fl1_fz);

                t_xxxy_xyyz[j] = fl_s_0_0 * (0.375 * pa2pb_y_z[j] * fl3_fx + 0.75 * pb_yz[j] * fl3_fx + 0.75 * pa2pb_xxy_z[j] * fl2_fx + 1.5 * pa2pb_xx_yz[j] * fl2_fx + 0.75 * pa2pb_xy_xz[j] * fl2_fx + 1.5 * pa2pb_x_xyz[j] * fl2_fx + 0.75 * pa2pb_y_yyz[j] * fl2_fx + 0.5 * pa2pb_xxxy_xz[j] * fl1_fx + pa2pb_xxx_xyz[j] * fl1_fx + 1.5 * pa2pb_xxy_yyz[j] * fl1_fx + 1.5 * pa2pb_xy_xyyz[j] * fl1_fx + pa2pb_xxxy_xyyz[j]);

                t_xxxy_xyyz[j] += fl_r_0_0 * (-0.75 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 7.5 * pb_yz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxy_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xx_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xyz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxy_xz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xyz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxy_yyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xyyz[j] * fl1_fz);

                t_xxxy_xyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 0.25 * pa2pb_xxx_x[j] * fl2_fx + 0.75 * pa2pb_xxy_y[j] * fl2_fx + 0.75 * pa2pb_xx_zz[j] * fl2_fx + 0.75 * pa2pb_xy_xy[j] * fl2_fx + 0.75 * pa2pb_x_xzz[j] * fl2_fx + 0.75 * pa2pb_y_yzz[j] * fl2_fx + 0.5 * pa2pb_xxxy_xy[j] * fl1_fx + 0.5 * pa2pb_xxx_xzz[j] * fl1_fx + 1.5 * pa2pb_xxy_yzz[j] * fl1_fx + 1.5 * pa2pb_xy_xyzz[j] * fl1_fx + pa2pb_xxxy_xyzz[j]);

                t_xxxy_xyzz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_xx[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xxy_y[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xx_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxy_xy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxx_xzz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxy_yzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xyzz[j] * fl1_fz);

                t_xxxy_xzzz[j] = fl_s_0_0 * (1.125 * pa2pb_y_z[j] * fl3_fx + 2.25 * pa2pb_xxy_z[j] * fl2_fx + 2.25 * pa2pb_xy_xz[j] * fl2_fx + 0.75 * pa2pb_y_zzz[j] * fl2_fx + 1.5 * pa2pb_xxxy_xz[j] * fl1_fx + 1.5 * pa2pb_xxy_zzz[j] * fl1_fx + 1.5 * pa2pb_xy_xzzz[j] * fl1_fx + pa2pb_xxxy_xzzz[j]);

                t_xxxy_xzzz[j] += fl_r_0_0 * (-2.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xxy_z[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxy_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xy_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxy_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxy_zzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_25_30(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_xy_yyyy = pa2pbDistances.data(1156 * idx + 165);

            auto pa2pb_xy_yyyz = pa2pbDistances.data(1156 * idx + 166);

            auto pa2pb_xy_yyzz = pa2pbDistances.data(1156 * idx + 167);

            auto pa2pb_xy_yzzz = pa2pbDistances.data(1156 * idx + 168);

            auto pa2pb_xy_zzzz = pa2pbDistances.data(1156 * idx + 169);

            auto pa2pb_xxx_y = pa2pbDistances.data(1156 * idx + 307);

            auto pa2pb_xxx_z = pa2pbDistances.data(1156 * idx + 308);

            auto pa2pb_xxx_yyy = pa2pbDistances.data(1156 * idx + 321);

            auto pa2pb_xxx_yyz = pa2pbDistances.data(1156 * idx + 322);

            auto pa2pb_xxx_yzz = pa2pbDistances.data(1156 * idx + 323);

            auto pa2pb_xxx_zzz = pa2pbDistances.data(1156 * idx + 324);

            auto pa2pb_xxxy_yy = pa2pbDistances.data(1156 * idx + 686);

            auto pa2pb_xxxy_yz = pa2pbDistances.data(1156 * idx + 687);

            auto pa2pb_xxxy_zz = pa2pbDistances.data(1156 * idx + 688);

            auto pa2pb_xxxy_yyyy = pa2pbDistances.data(1156 * idx + 709);

            auto pa2pb_xxxy_yyyz = pa2pbDistances.data(1156 * idx + 710);

            auto pa2pb_xxxy_yyzz = pa2pbDistances.data(1156 * idx + 711);

            auto pa2pb_xxxy_yzzz = pa2pbDistances.data(1156 * idx + 712);

            auto pa2pb_xxxy_zzzz = pa2pbDistances.data(1156 * idx + 713);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxy_yyyy = primBuffer.data(225 * idx + 25);

            auto t_xxxy_yyyz = primBuffer.data(225 * idx + 26);

            auto t_xxxy_yyzz = primBuffer.data(225 * idx + 27);

            auto t_xxxy_yzzz = primBuffer.data(225 * idx + 28);

            auto t_xxxy_zzzz = primBuffer.data(225 * idx + 29);

            // Batch of Integrals (25,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, \
                                     pa2pb_x_z, pa2pb_x_zzz, pa2pb_xxx_y, pa2pb_xxx_yyy, pa2pb_xxx_yyz, \
                                     pa2pb_xxx_yzz, pa2pb_xxx_z, pa2pb_xxx_zzz, pa2pb_xxxy_yy, pa2pb_xxxy_yyyy, \
                                     pa2pb_xxxy_yyyz, pa2pb_xxxy_yyzz, pa2pb_xxxy_yz, pa2pb_xxxy_yzzz, pa2pb_xxxy_zz, \
                                     pa2pb_xxxy_zzzz, pa2pb_xy_yy, pa2pb_xy_yyyy, pa2pb_xy_yyyz, pa2pb_xy_yyzz, \
                                     pa2pb_xy_yz, pa2pb_xy_yzzz, pa2pb_xy_zz, pa2pb_xy_zzzz, pa_xxxy, pa_xy, r_0_0, s_0_0, \
                                     t_xxxy_yyyy, t_xxxy_yyyz, t_xxxy_yyzz, t_xxxy_yzzz, t_xxxy_zzzz: VLX_ALIGN)
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

                t_xxxy_yyyy[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 4.5 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa_xxxy[j] * fl2_fx + 3.0 * pa2pb_xxx_y[j] * fl2_fx + 4.5 * pa2pb_xy_yy[j] * fl2_fx + 3.0 * pa2pb_x_yyy[j] * fl2_fx + 3.0 * pa2pb_xxxy_yy[j] * fl1_fx + 2.0 * pa2pb_xxx_yyy[j] * fl1_fx + 1.5 * pa2pb_xy_yyyy[j] * fl1_fx + pa2pb_xxxy_yyyy[j]);

                t_xxxy_yyyy[j] += fl_r_0_0 * (-4.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxxy[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_xy[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 9.0 * pa_xxxy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxxy_yy[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xy_yy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxxy_yy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxx_yyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_yyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_yyyy[j] * fl1_fz);

                t_xxxy_yyyz[j] = fl_s_0_0 * (1.125 * pa2pb_x_z[j] * fl3_fx + 0.75 * pa2pb_xxx_z[j] * fl2_fx + 2.25 * pa2pb_xy_yz[j] * fl2_fx + 2.25 * pa2pb_x_yyz[j] * fl2_fx + 1.5 * pa2pb_xxxy_yz[j] * fl1_fx + 1.5 * pa2pb_xxx_yyz[j] * fl1_fx + 1.5 * pa2pb_xy_yyyz[j] * fl1_fx + pa2pb_xxxy_yyyz[j]);

                t_xxxy_yyyz[j] += fl_r_0_0 * (-2.25 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_yyz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxxy_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xy_yz[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxy_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxx_yyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_yyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_yyyz[j] * fl1_fz);

                t_xxxy_yyzz[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 0.75 * pa2pb_x_y[j] * fl3_fx + 0.25 * pa_xxxy[j] * fl2_fx + 0.5 * pa2pb_xxx_y[j] * fl2_fx + 0.75 * pa2pb_xy_yy[j] * fl2_fx + 0.75 * pa2pb_xy_zz[j] * fl2_fx + 1.5 * pa2pb_x_yzz[j] * fl2_fx + 0.5 * pa2pb_xxxy_yy[j] * fl1_fx + 0.5 * pa2pb_xxxy_zz[j] * fl1_fx + pa2pb_xxx_yzz[j] * fl1_fx + 1.5 * pa2pb_xy_yyzz[j] * fl1_fx + pa2pb_xxxy_yyzz[j]);

                t_xxxy_yyzz[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxxy[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xy[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 3.0 * pa_xxxy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_yzz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxxy_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxxy_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_yy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xy_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxy_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxy_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_yzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_yyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_yyzz[j] * fl1_fz);

                t_xxxy_yzzz[j] = fl_s_0_0 * (1.125 * pa2pb_x_z[j] * fl3_fx + 0.75 * pa2pb_xxx_z[j] * fl2_fx + 2.25 * pa2pb_xy_yz[j] * fl2_fx + 0.75 * pa2pb_x_zzz[j] * fl2_fx + 1.5 * pa2pb_xxxy_yz[j] * fl1_fx + 0.5 * pa2pb_xxx_zzz[j] * fl1_fx + 1.5 * pa2pb_xy_yzzz[j] * fl1_fx + pa2pb_xxxy_yzzz[j]);

                t_xxxy_yzzz[j] += fl_r_0_0 * (-2.25 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxxy_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xy_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxy_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxx_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_yzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_yzzz[j] * fl1_fz);

                t_xxxy_zzzz[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 0.75 * pa_xxxy[j] * fl2_fx + 4.5 * pa2pb_xy_zz[j] * fl2_fx + 3.0 * pa2pb_xxxy_zz[j] * fl1_fx + 1.5 * pa2pb_xy_zzzz[j] * fl1_fx + pa2pb_xxxy_zzzz[j]);

                t_xxxy_zzzz[j] += fl_r_0_0 * (-4.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxxy[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_xy[j] * fl3_fx * fl1_fz + 9.0 * pa_xxxy[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxxy_zz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xy_zz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxxy_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_zzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_zzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxy_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_30_35(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_xxxx = pa2pbDistances.data(1156 * idx + 189);

            auto pa2pb_xz_xxxy = pa2pbDistances.data(1156 * idx + 190);

            auto pa2pb_xz_xxxz = pa2pbDistances.data(1156 * idx + 191);

            auto pa2pb_xz_xxyy = pa2pbDistances.data(1156 * idx + 192);

            auto pa2pb_xz_xxyz = pa2pbDistances.data(1156 * idx + 193);

            auto pa2pb_xxx_x = pa2pbDistances.data(1156 * idx + 306);

            auto pa2pb_xxx_y = pa2pbDistances.data(1156 * idx + 307);

            auto pa2pb_xxx_xxx = pa2pbDistances.data(1156 * idx + 315);

            auto pa2pb_xxx_xxy = pa2pbDistances.data(1156 * idx + 316);

            auto pa2pb_xxz_x = pa2pbDistances.data(1156 * idx + 374);

            auto pa2pb_xxz_y = pa2pbDistances.data(1156 * idx + 375);

            auto pa2pb_xxz_z = pa2pbDistances.data(1156 * idx + 376);

            auto pa2pb_xxz_xxx = pa2pbDistances.data(1156 * idx + 383);

            auto pa2pb_xxz_xxy = pa2pbDistances.data(1156 * idx + 384);

            auto pa2pb_xxz_xxz = pa2pbDistances.data(1156 * idx + 385);

            auto pa2pb_xxz_xyy = pa2pbDistances.data(1156 * idx + 386);

            auto pa2pb_xxz_xyz = pa2pbDistances.data(1156 * idx + 387);

            auto pa2pb_xxxz_xx = pa2pbDistances.data(1156 * idx + 717);

            auto pa2pb_xxxz_xy = pa2pbDistances.data(1156 * idx + 718);

            auto pa2pb_xxxz_xz = pa2pbDistances.data(1156 * idx + 719);

            auto pa2pb_xxxz_yy = pa2pbDistances.data(1156 * idx + 720);

            auto pa2pb_xxxz_yz = pa2pbDistances.data(1156 * idx + 721);

            auto pa2pb_xxxz_xxxx = pa2pbDistances.data(1156 * idx + 733);

            auto pa2pb_xxxz_xxxy = pa2pbDistances.data(1156 * idx + 734);

            auto pa2pb_xxxz_xxxz = pa2pbDistances.data(1156 * idx + 735);

            auto pa2pb_xxxz_xxyy = pa2pbDistances.data(1156 * idx + 736);

            auto pa2pb_xxxz_xxyz = pa2pbDistances.data(1156 * idx + 737);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxz_xxxx = primBuffer.data(225 * idx + 30);

            auto t_xxxz_xxxy = primBuffer.data(225 * idx + 31);

            auto t_xxxz_xxxz = primBuffer.data(225 * idx + 32);

            auto t_xxxz_xxyy = primBuffer.data(225 * idx + 33);

            auto t_xxxz_xxyz = primBuffer.data(225 * idx + 34);

            // Batch of Integrals (30,35)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_y, \
                                     pa2pb_xx_xx, pa2pb_xx_xy, pa2pb_xxx_x, pa2pb_xxx_xxx, pa2pb_xxx_xxy, pa2pb_xxx_y, \
                                     pa2pb_xxxz_xx, pa2pb_xxxz_xxxx, pa2pb_xxxz_xxxy, pa2pb_xxxz_xxxz, pa2pb_xxxz_xxyy, \
                                     pa2pb_xxxz_xxyz, pa2pb_xxxz_xy, pa2pb_xxxz_xz, pa2pb_xxxz_yy, pa2pb_xxxz_yz, \
                                     pa2pb_xxz_x, pa2pb_xxz_xxx, pa2pb_xxz_xxy, pa2pb_xxz_xxz, pa2pb_xxz_xyy, \
                                     pa2pb_xxz_xyz, pa2pb_xxz_y, pa2pb_xxz_z, pa2pb_xz_xx, pa2pb_xz_xxxx, pa2pb_xz_xxxy, \
                                     pa2pb_xz_xxxz, pa2pb_xz_xxyy, pa2pb_xz_xxyz, pa2pb_xz_xy, pa2pb_xz_xz, pa2pb_xz_yy, \
                                     pa2pb_xz_yz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_xxz, pa2pb_z_xyy, \
                                     pa2pb_z_xyz, pa2pb_z_y, pa2pb_z_z, pa_xx, pa_xxxz, pa_xz, pb_xx, pb_xy, r_0_0, s_0_0, \
                                     t_xxxz_xxxx, t_xxxz_xxxy, t_xxxz_xxxz, t_xxxz_xxyy, t_xxxz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxz_xxxx[j] = fl_s_0_0 * (5.625 * pa_xz[j] * fl3_fx + 7.5 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa_xxxz[j] * fl2_fx + 9.0 * pa2pb_xxz_x[j] * fl2_fx + 13.5 * pa2pb_xz_xx[j] * fl2_fx + 3.0 * pa2pb_z_xxx[j] * fl2_fx + 3.0 * pa2pb_xxxz_xx[j] * fl1_fx + 6.0 * pa2pb_xxz_xxx[j] * fl1_fx + 1.5 * pa2pb_xz_xxxx[j] * fl1_fx + pa2pb_xxxz_xxxx[j]);

                t_xxxz_xxxx[j] += fl_r_0_0 * (-13.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 56.25 * pa_xz[j] * fl3_fx * fl1_fz + 75.0 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 2.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xxxz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_xxxz[j] * fl1_fz * fl2_fx + 108.0 * pa2pb_xxz_x[j] * fl2_fx * fl1_fz + 162.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xxxz_xx[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxxz_xx[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_xxz_xxx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxxx[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxxx[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xxxx[j] * fl1_fz);

                t_xxxz_xxxy[j] = fl_s_0_0 * (1.875 * pa2pb_z_y[j] * fl3_fx + 2.25 * pa2pb_xxz_y[j] * fl2_fx + 6.75 * pa2pb_xz_xy[j] * fl2_fx + 2.25 * pa2pb_z_xxy[j] * fl2_fx + 1.5 * pa2pb_xxxz_xy[j] * fl1_fx + 4.5 * pa2pb_xxz_xxy[j] * fl1_fx + 1.5 * pa2pb_xz_xxxy[j] * fl1_fx + pa2pb_xxxz_xxxy[j]);

                t_xxxz_xxxy[j] += fl_r_0_0 * (18.75 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 27.0 * pa2pb_xxz_y[j] * fl2_fx * fl1_fz + 81.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxz_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxz_xy[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_xxz_xxy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxxy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxxy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xxxy[j] * fl1_fz);

                t_xxxz_xxxz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.125 * pa_xx[j] * fl3_fx + 3.375 * pa2pb_x_x[j] * fl3_fx + 1.875 * pa2pb_z_z[j] * fl3_fx + 1.125 * pb_xx[j] * fl3_fx + 0.75 * pa2pb_xxx_x[j] * fl2_fx + 2.25 * pa2pb_xxz_z[j] * fl2_fx + 2.25 * pa2pb_xx_xx[j] * fl2_fx + 6.75 * pa2pb_xz_xz[j] * fl2_fx + 0.75 * pa2pb_x_xxx[j] * fl2_fx + 2.25 * pa2pb_z_xxz[j] * fl2_fx + 1.5 * pa2pb_xxxz_xz[j] * fl1_fx + 0.5 * pa2pb_xxx_xxx[j] * fl1_fx + 4.5 * pa2pb_xxz_xxz[j] * fl1_fx + 1.5 * pa2pb_xz_xxxz[j] * fl1_fx + pa2pb_xxxz_xxxz[j]);

                t_xxxz_xxxz[j] += fl_r_0_0 * (7.5 * fl4_fx * fl1_fz - 1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xx[j] * fl3_fx * fl1_fz + 33.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 18.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pb_xx[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xxz_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xx_xx[j] * fl2_fx * fl1_fz + 81.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxz_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxx_xxx[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_xxz_xxz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxxz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxxz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xxxz[j] * fl1_fz);

                t_xxxz_xxyy[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 0.75 * pa2pb_z_x[j] * fl3_fx + 0.25 * pa_xxxz[j] * fl2_fx + 1.5 * pa2pb_xxz_x[j] * fl2_fx + 2.25 * pa2pb_xz_yy[j] * fl2_fx + 0.75 * pa2pb_xz_xx[j] * fl2_fx + 1.5 * pa2pb_z_xyy[j] * fl2_fx + 0.5 * pa2pb_xxxz_xx[j] * fl1_fx + 0.5 * pa2pb_xxxz_yy[j] * fl1_fx + 3.0 * pa2pb_xxz_xyy[j] * fl1_fx + 1.5 * pa2pb_xz_xxyy[j] * fl1_fx + pa2pb_xxxz_xxyy[j]);

                t_xxxz_xxyy[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xz[j] * fl3_fx * fl1_fz - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xxxz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 3.0 * pa_xxxz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xxz_x[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxxz_yy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxz_yy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxz_xyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xxyy[j] * fl1_fz);

                t_xxxz_xxyz[j] = fl_s_0_0 * (1.125 * pa2pb_x_y[j] * fl3_fx + 0.75 * pb_xy[j] * fl3_fx + 0.25 * pa2pb_xxx_y[j] * fl2_fx + 1.5 * pa2pb_xx_xy[j] * fl2_fx + 2.25 * pa2pb_xz_yz[j] * fl2_fx + 0.75 * pa2pb_x_xxy[j] * fl2_fx + 1.5 * pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_xxxz_yz[j] * fl1_fx + 0.5 * pa2pb_xxx_xxy[j] * fl1_fx + 3.0 * pa2pb_xxz_xyz[j] * fl1_fx + 1.5 * pa2pb_xz_xxyz[j] * fl1_fx + pa2pb_xxxz_xxyz[j]);

                t_xxxz_xxyz[j] += fl_r_0_0 * (11.25 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_xy[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xx_xy[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xxy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxz_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxx_xxy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxz_xyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_35_40(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_xz_xxzz = pa2pbDistances.data(1156 * idx + 194);

            auto pa2pb_xz_xyyy = pa2pbDistances.data(1156 * idx + 195);

            auto pa2pb_xz_xyyz = pa2pbDistances.data(1156 * idx + 196);

            auto pa2pb_xz_xyzz = pa2pbDistances.data(1156 * idx + 197);

            auto pa2pb_xz_xzzz = pa2pbDistances.data(1156 * idx + 198);

            auto pa2pb_xxx_x = pa2pbDistances.data(1156 * idx + 306);

            auto pa2pb_xxx_z = pa2pbDistances.data(1156 * idx + 308);

            auto pa2pb_xxx_xxz = pa2pbDistances.data(1156 * idx + 317);

            auto pa2pb_xxx_xyy = pa2pbDistances.data(1156 * idx + 318);

            auto pa2pb_xxx_xyz = pa2pbDistances.data(1156 * idx + 319);

            auto pa2pb_xxx_xzz = pa2pbDistances.data(1156 * idx + 320);

            auto pa2pb_xxz_x = pa2pbDistances.data(1156 * idx + 374);

            auto pa2pb_xxz_y = pa2pbDistances.data(1156 * idx + 375);

            auto pa2pb_xxz_z = pa2pbDistances.data(1156 * idx + 376);

            auto pa2pb_xxz_xzz = pa2pbDistances.data(1156 * idx + 388);

            auto pa2pb_xxz_yyy = pa2pbDistances.data(1156 * idx + 389);

            auto pa2pb_xxz_yyz = pa2pbDistances.data(1156 * idx + 390);

            auto pa2pb_xxz_yzz = pa2pbDistances.data(1156 * idx + 391);

            auto pa2pb_xxz_zzz = pa2pbDistances.data(1156 * idx + 392);

            auto pa2pb_xxxz_xx = pa2pbDistances.data(1156 * idx + 717);

            auto pa2pb_xxxz_xy = pa2pbDistances.data(1156 * idx + 718);

            auto pa2pb_xxxz_xz = pa2pbDistances.data(1156 * idx + 719);

            auto pa2pb_xxxz_zz = pa2pbDistances.data(1156 * idx + 722);

            auto pa2pb_xxxz_xxzz = pa2pbDistances.data(1156 * idx + 738);

            auto pa2pb_xxxz_xyyy = pa2pbDistances.data(1156 * idx + 739);

            auto pa2pb_xxxz_xyyz = pa2pbDistances.data(1156 * idx + 740);

            auto pa2pb_xxxz_xyzz = pa2pbDistances.data(1156 * idx + 741);

            auto pa2pb_xxxz_xzzz = pa2pbDistances.data(1156 * idx + 742);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxz_xxzz = primBuffer.data(225 * idx + 35);

            auto t_xxxz_xyyy = primBuffer.data(225 * idx + 36);

            auto t_xxxz_xyyz = primBuffer.data(225 * idx + 37);

            auto t_xxxz_xyzz = primBuffer.data(225 * idx + 38);

            auto t_xxxz_xzzz = primBuffer.data(225 * idx + 39);

            // Batch of Integrals (35,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxz, pa2pb_x_xyy, pa2pb_x_xyz, \
                                     pa2pb_x_xzz, pa2pb_x_z, pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, \
                                     pa2pb_xxx_x, pa2pb_xxx_xxz, pa2pb_xxx_xyy, pa2pb_xxx_xyz, pa2pb_xxx_xzz, \
                                     pa2pb_xxx_z, pa2pb_xxxz_xx, pa2pb_xxxz_xxzz, pa2pb_xxxz_xy, pa2pb_xxxz_xyyy, \
                                     pa2pb_xxxz_xyyz, pa2pb_xxxz_xyzz, pa2pb_xxxz_xz, pa2pb_xxxz_xzzz, pa2pb_xxxz_zz, \
                                     pa2pb_xxz_x, pa2pb_xxz_xzz, pa2pb_xxz_y, pa2pb_xxz_yyy, pa2pb_xxz_yyz, \
                                     pa2pb_xxz_yzz, pa2pb_xxz_z, pa2pb_xxz_zzz, pa2pb_xz_xx, pa2pb_xz_xxzz, pa2pb_xz_xy, \
                                     pa2pb_xz_xyyy, pa2pb_xz_xyyz, pa2pb_xz_xyzz, pa2pb_xz_xz, pa2pb_xz_xzzz, \
                                     pa2pb_xz_zz, pa2pb_z_x, pa2pb_z_xzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, \
                                     pa2pb_z_yzz, pa2pb_z_z, pa2pb_z_zzz, pa_xx, pa_xxxz, pa_xz, pb_xz, pb_yy, pb_yz, pb_zz, \
                                     r_0_0, s_0_0, t_xxxz_xxzz, t_xxxz_xyyy, t_xxxz_xyyz, t_xxxz_xyzz, t_xxxz_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxz_xxzz[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 2.25 * pa2pb_x_z[j] * fl3_fx + 0.75 * pa2pb_z_x[j] * fl3_fx + 1.5 * pb_xz[j] * fl3_fx + 0.25 * pa_xxxz[j] * fl2_fx + 0.5 * pa2pb_xxx_z[j] * fl2_fx + 1.5 * pa2pb_xxz_x[j] * fl2_fx + 3.0 * pa2pb_xx_xz[j] * fl2_fx + 2.25 * pa2pb_xz_zz[j] * fl2_fx + 0.75 * pa2pb_xz_xx[j] * fl2_fx + 1.5 * pa2pb_x_xxz[j] * fl2_fx + 1.5 * pa2pb_z_xzz[j] * fl2_fx + 0.5 * pa2pb_xxxz_xx[j] * fl1_fx + 0.5 * pa2pb_xxxz_zz[j] * fl1_fx + pa2pb_xxx_xxz[j] * fl1_fx + 3.0 * pa2pb_xxz_xzz[j] * fl1_fx + 1.5 * pa2pb_xz_xxzz[j] * fl1_fx + pa2pb_xxxz_xxzz[j]);

                t_xxxz_xxzz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xz[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - pa_xxxz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 15.0 * pb_xz[j] * fl3_fx * fl1_fz + 3.0 * pa_xxxz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xxz_x[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_xz[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xxz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxxz_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxz_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xxz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxz_xzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xxzz[j] * fl1_fz);

                t_xxxz_xyyy[j] = fl_s_0_0 * (1.125 * pa2pb_z_y[j] * fl3_fx + 2.25 * pa2pb_xxz_y[j] * fl2_fx + 2.25 * pa2pb_xz_xy[j] * fl2_fx + 0.75 * pa2pb_z_yyy[j] * fl2_fx + 1.5 * pa2pb_xxxz_xy[j] * fl1_fx + 1.5 * pa2pb_xxz_yyy[j] * fl1_fx + 1.5 * pa2pb_xz_xyyy[j] * fl1_fx + pa2pb_xxxz_xyyy[j]);

                t_xxxz_xyyy[j] += fl_r_0_0 * (-2.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xxz_y[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxz_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxz_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxz_yyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xyyy[j] * fl1_fz);

                t_xxxz_xyyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 0.25 * pa2pb_xxx_x[j] * fl2_fx + 0.75 * pa2pb_xxz_z[j] * fl2_fx + 0.75 * pa2pb_xx_yy[j] * fl2_fx + 0.75 * pa2pb_xz_xz[j] * fl2_fx + 0.75 * pa2pb_x_xyy[j] * fl2_fx + 0.75 * pa2pb_z_yyz[j] * fl2_fx + 0.5 * pa2pb_xxxz_xz[j] * fl1_fx + 0.5 * pa2pb_xxx_xyy[j] * fl1_fx + 1.5 * pa2pb_xxz_yyz[j] * fl1_fx + 1.5 * pa2pb_xz_xyyz[j] * fl1_fx + pa2pb_xxxz_xyyz[j]);

                t_xxxz_xyyz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_xx[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xxz_z[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xx_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xyy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxz_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxx_xyy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxz_yyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xyyz[j] * fl1_fz);

                t_xxxz_xyzz[j] = fl_s_0_0 * (0.375 * pa2pb_z_y[j] * fl3_fx + 0.75 * pb_yz[j] * fl3_fx + 0.75 * pa2pb_xxz_y[j] * fl2_fx + 1.5 * pa2pb_xx_yz[j] * fl2_fx + 0.75 * pa2pb_xz_xy[j] * fl2_fx + 1.5 * pa2pb_x_xyz[j] * fl2_fx + 0.75 * pa2pb_z_yzz[j] * fl2_fx + 0.5 * pa2pb_xxxz_xy[j] * fl1_fx + pa2pb_xxx_xyz[j] * fl1_fx + 1.5 * pa2pb_xxz_yzz[j] * fl1_fx + 1.5 * pa2pb_xz_xyzz[j] * fl1_fx + pa2pb_xxxz_xyzz[j]);

                t_xxxz_xyzz[j] += fl_r_0_0 * (-0.75 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 7.5 * pb_yz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxz_y[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xx_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xyz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxz_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_xyz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxz_yzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xyzz[j] * fl1_fz);

                t_xxxz_xzzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_xx[j] * fl3_fx + 1.125 * pa2pb_x_x[j] * fl3_fx + 1.125 * pa2pb_z_z[j] * fl3_fx + 1.125 * pb_zz[j] * fl3_fx + 0.75 * pa2pb_xxx_x[j] * fl2_fx + 2.25 * pa2pb_xxz_z[j] * fl2_fx + 2.25 * pa2pb_xx_zz[j] * fl2_fx + 2.25 * pa2pb_xz_xz[j] * fl2_fx + 2.25 * pa2pb_x_xzz[j] * fl2_fx + 0.75 * pa2pb_z_zzz[j] * fl2_fx + 1.5 * pa2pb_xxxz_xz[j] * fl1_fx + 1.5 * pa2pb_xxx_xzz[j] * fl1_fx + 1.5 * pa2pb_xxz_zzz[j] * fl1_fx + 1.5 * pa2pb_xz_xzzz[j] * fl1_fx + pa2pb_xxxz_xzzz[j]);

                t_xxxz_xzzz[j] += fl_r_0_0 * (-1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 4.5 * fl4_fx * fl1_fz + 11.25 * pa_xx[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xxx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 11.25 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxx_x[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xxz_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xx_zz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxz_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxz_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxx_xzz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxz_zzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_40_45(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
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

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_xz_yyyy = pa2pbDistances.data(1156 * idx + 199);

            auto pa2pb_xz_yyyz = pa2pbDistances.data(1156 * idx + 200);

            auto pa2pb_xz_yyzz = pa2pbDistances.data(1156 * idx + 201);

            auto pa2pb_xz_yzzz = pa2pbDistances.data(1156 * idx + 202);

            auto pa2pb_xz_zzzz = pa2pbDistances.data(1156 * idx + 203);

            auto pa2pb_xxx_y = pa2pbDistances.data(1156 * idx + 307);

            auto pa2pb_xxx_z = pa2pbDistances.data(1156 * idx + 308);

            auto pa2pb_xxx_yyy = pa2pbDistances.data(1156 * idx + 321);

            auto pa2pb_xxx_yyz = pa2pbDistances.data(1156 * idx + 322);

            auto pa2pb_xxx_yzz = pa2pbDistances.data(1156 * idx + 323);

            auto pa2pb_xxx_zzz = pa2pbDistances.data(1156 * idx + 324);

            auto pa2pb_xxxz_yy = pa2pbDistances.data(1156 * idx + 720);

            auto pa2pb_xxxz_yz = pa2pbDistances.data(1156 * idx + 721);

            auto pa2pb_xxxz_zz = pa2pbDistances.data(1156 * idx + 722);

            auto pa2pb_xxxz_yyyy = pa2pbDistances.data(1156 * idx + 743);

            auto pa2pb_xxxz_yyyz = pa2pbDistances.data(1156 * idx + 744);

            auto pa2pb_xxxz_yyzz = pa2pbDistances.data(1156 * idx + 745);

            auto pa2pb_xxxz_yzzz = pa2pbDistances.data(1156 * idx + 746);

            auto pa2pb_xxxz_zzzz = pa2pbDistances.data(1156 * idx + 747);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxz_yyyy = primBuffer.data(225 * idx + 40);

            auto t_xxxz_yyyz = primBuffer.data(225 * idx + 41);

            auto t_xxxz_yyzz = primBuffer.data(225 * idx + 42);

            auto t_xxxz_yzzz = primBuffer.data(225 * idx + 43);

            auto t_xxxz_zzzz = primBuffer.data(225 * idx + 44);

            // Batch of Integrals (40,45)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, \
                                     pa2pb_x_z, pa2pb_x_zzz, pa2pb_xxx_y, pa2pb_xxx_yyy, pa2pb_xxx_yyz, \
                                     pa2pb_xxx_yzz, pa2pb_xxx_z, pa2pb_xxx_zzz, pa2pb_xxxz_yy, pa2pb_xxxz_yyyy, \
                                     pa2pb_xxxz_yyyz, pa2pb_xxxz_yyzz, pa2pb_xxxz_yz, pa2pb_xxxz_yzzz, pa2pb_xxxz_zz, \
                                     pa2pb_xxxz_zzzz, pa2pb_xz_yy, pa2pb_xz_yyyy, pa2pb_xz_yyyz, pa2pb_xz_yyzz, \
                                     pa2pb_xz_yz, pa2pb_xz_yzzz, pa2pb_xz_zz, pa2pb_xz_zzzz, pa_xxxz, pa_xz, r_0_0, s_0_0, \
                                     t_xxxz_yyyy, t_xxxz_yyyz, t_xxxz_yyzz, t_xxxz_yzzz, t_xxxz_zzzz: VLX_ALIGN)
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

                t_xxxz_yyyy[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 0.75 * pa_xxxz[j] * fl2_fx + 4.5 * pa2pb_xz_yy[j] * fl2_fx + 3.0 * pa2pb_xxxz_yy[j] * fl1_fx + 1.5 * pa2pb_xz_yyyy[j] * fl1_fx + pa2pb_xxxz_yyyy[j]);

                t_xxxz_yyyy[j] += fl_r_0_0 * (-4.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxxz[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_xz[j] * fl3_fx * fl1_fz + 9.0 * pa_xxxz[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxxz_yy[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxxz_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_yyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_yyyy[j] * fl1_fz);

                t_xxxz_yyyz[j] = fl_s_0_0 * (1.125 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa2pb_xxx_y[j] * fl2_fx + 2.25 * pa2pb_xz_yz[j] * fl2_fx + 0.75 * pa2pb_x_yyy[j] * fl2_fx + 1.5 * pa2pb_xxxz_yz[j] * fl1_fx + 0.5 * pa2pb_xxx_yyy[j] * fl1_fx + 1.5 * pa2pb_xz_yyyz[j] * fl1_fx + pa2pb_xxxz_yyyz[j]);

                t_xxxz_yyyz[j] += fl_r_0_0 * (-2.25 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxxz_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxz_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxx_yyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_yyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_yyyz[j] * fl1_fz);

                t_xxxz_yyzz[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 0.75 * pa2pb_x_z[j] * fl3_fx + 0.25 * pa_xxxz[j] * fl2_fx + 0.5 * pa2pb_xxx_z[j] * fl2_fx + 0.75 * pa2pb_xz_yy[j] * fl2_fx + 0.75 * pa2pb_xz_zz[j] * fl2_fx + 1.5 * pa2pb_x_yyz[j] * fl2_fx + 0.5 * pa2pb_xxxz_yy[j] * fl1_fx + 0.5 * pa2pb_xxxz_zz[j] * fl1_fx + pa2pb_xxx_yyz[j] * fl1_fx + 1.5 * pa2pb_xz_yyzz[j] * fl1_fx + pa2pb_xxxz_yyzz[j]);

                t_xxxz_yyzz[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxxz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xz[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 3.0 * pa_xxxz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_yyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxxz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxxz_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxxz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxxz_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxx_yyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_yyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_yyzz[j] * fl1_fz);

                t_xxxz_yzzz[j] = fl_s_0_0 * (1.125 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa2pb_xxx_y[j] * fl2_fx + 2.25 * pa2pb_xz_yz[j] * fl2_fx + 2.25 * pa2pb_x_yzz[j] * fl2_fx + 1.5 * pa2pb_xxxz_yz[j] * fl1_fx + 1.5 * pa2pb_xxx_yzz[j] * fl1_fx + 1.5 * pa2pb_xz_yzzz[j] * fl1_fx + pa2pb_xxxz_yzzz[j]);

                t_xxxz_yzzz[j] += fl_r_0_0 * (-2.25 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxx_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxx_y[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_yzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxxz_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxxz_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxx_yzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_yzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_yzzz[j] * fl1_fz);

                t_xxxz_zzzz[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 4.5 * pa2pb_x_z[j] * fl3_fx + 0.75 * pa_xxxz[j] * fl2_fx + 3.0 * pa2pb_xxx_z[j] * fl2_fx + 4.5 * pa2pb_xz_zz[j] * fl2_fx + 3.0 * pa2pb_x_zzz[j] * fl2_fx + 3.0 * pa2pb_xxxz_zz[j] * fl1_fx + 2.0 * pa2pb_xxx_zzz[j] * fl1_fx + 1.5 * pa2pb_xz_zzzz[j] * fl1_fx + pa2pb_xxxz_zzzz[j]);

                t_xxxz_zzzz[j] += fl_r_0_0 * (-4.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxxz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xxx_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_xz[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 9.0 * pa_xxxz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xxx_z[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxxz_zz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxxz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxx_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_zzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_zzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxxz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_45_50(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

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

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_xxxx = pa2pbDistances.data(1156 * idx + 121);

            auto pa2pb_xx_xxxy = pa2pbDistances.data(1156 * idx + 122);

            auto pa2pb_xx_xxxz = pa2pbDistances.data(1156 * idx + 123);

            auto pa2pb_xx_xxyy = pa2pbDistances.data(1156 * idx + 124);

            auto pa2pb_xx_xxyz = pa2pbDistances.data(1156 * idx + 125);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_xxxx = pa2pbDistances.data(1156 * idx + 223);

            auto pa2pb_yy_xxxy = pa2pbDistances.data(1156 * idx + 224);

            auto pa2pb_yy_xxxz = pa2pbDistances.data(1156 * idx + 225);

            auto pa2pb_yy_xxyy = pa2pbDistances.data(1156 * idx + 226);

            auto pa2pb_yy_xxyz = pa2pbDistances.data(1156 * idx + 227);

            auto pa2pb_xxy_x = pa2pbDistances.data(1156 * idx + 340);

            auto pa2pb_xxy_y = pa2pbDistances.data(1156 * idx + 341);

            auto pa2pb_xxy_z = pa2pbDistances.data(1156 * idx + 342);

            auto pa2pb_xxy_xxx = pa2pbDistances.data(1156 * idx + 349);

            auto pa2pb_xxy_xxy = pa2pbDistances.data(1156 * idx + 350);

            auto pa2pb_xxy_xxz = pa2pbDistances.data(1156 * idx + 351);

            auto pa2pb_xyy_x = pa2pbDistances.data(1156 * idx + 408);

            auto pa2pb_xyy_y = pa2pbDistances.data(1156 * idx + 409);

            auto pa2pb_xyy_z = pa2pbDistances.data(1156 * idx + 410);

            auto pa2pb_xyy_xxx = pa2pbDistances.data(1156 * idx + 417);

            auto pa2pb_xyy_xxy = pa2pbDistances.data(1156 * idx + 418);

            auto pa2pb_xyy_xxz = pa2pbDistances.data(1156 * idx + 419);

            auto pa2pb_xyy_xyy = pa2pbDistances.data(1156 * idx + 420);

            auto pa2pb_xyy_xyz = pa2pbDistances.data(1156 * idx + 421);

            auto pa2pb_xxyy_xx = pa2pbDistances.data(1156 * idx + 751);

            auto pa2pb_xxyy_xy = pa2pbDistances.data(1156 * idx + 752);

            auto pa2pb_xxyy_xz = pa2pbDistances.data(1156 * idx + 753);

            auto pa2pb_xxyy_yy = pa2pbDistances.data(1156 * idx + 754);

            auto pa2pb_xxyy_yz = pa2pbDistances.data(1156 * idx + 755);

            auto pa2pb_xxyy_xxxx = pa2pbDistances.data(1156 * idx + 767);

            auto pa2pb_xxyy_xxxy = pa2pbDistances.data(1156 * idx + 768);

            auto pa2pb_xxyy_xxxz = pa2pbDistances.data(1156 * idx + 769);

            auto pa2pb_xxyy_xxyy = pa2pbDistances.data(1156 * idx + 770);

            auto pa2pb_xxyy_xxyz = pa2pbDistances.data(1156 * idx + 771);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyy_xxxx = primBuffer.data(225 * idx + 45);

            auto t_xxyy_xxxy = primBuffer.data(225 * idx + 46);

            auto t_xxyy_xxxz = primBuffer.data(225 * idx + 47);

            auto t_xxyy_xxyy = primBuffer.data(225 * idx + 48);

            auto t_xxyy_xxyz = primBuffer.data(225 * idx + 49);

            // Batch of Integrals (45,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xxxx, \
                                     pa2pb_xx_xxxy, pa2pb_xx_xxxz, pa2pb_xx_xxyy, pa2pb_xx_xxyz, pa2pb_xx_xy, \
                                     pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xxy_x, pa2pb_xxy_xxx, pa2pb_xxy_xxy, \
                                     pa2pb_xxy_xxz, pa2pb_xxy_y, pa2pb_xxy_z, pa2pb_xxyy_xx, pa2pb_xxyy_xxxx, \
                                     pa2pb_xxyy_xxxy, pa2pb_xxyy_xxxz, pa2pb_xxyy_xxyy, pa2pb_xxyy_xxyz, pa2pb_xxyy_xy, \
                                     pa2pb_xxyy_xz, pa2pb_xxyy_yy, pa2pb_xxyy_yz, pa2pb_xy_xx, pa2pb_xy_xy, pa2pb_xy_xz, \
                                     pa2pb_xyy_x, pa2pb_xyy_xxx, pa2pb_xyy_xxy, pa2pb_xyy_xxz, pa2pb_xyy_xyy, \
                                     pa2pb_xyy_xyz, pa2pb_xyy_y, pa2pb_xyy_z, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, \
                                     pa2pb_y_xxz, pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xxxx, pa2pb_yy_xxxy, \
                                     pa2pb_yy_xxxz, pa2pb_yy_xxyy, pa2pb_yy_xxyz, pa2pb_yy_xy, pa2pb_yy_xz, pa2pb_yy_yy, \
                                     pa2pb_yy_yz, pa_xx, pa_xxyy, pa_xy, pa_yy, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, \
                                     pb_xxyz, pb_xy, pb_xz, pb_yy, pb_yz, r_0_0, s_0_0, t_xxyy_xxxx, t_xxyy_xxxy, \
                                     t_xxyy_xxxz, t_xxyy_xxyy, t_xxyy_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxyy_xxxx[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_yy[j] * fl3_fx + 0.375 * pa_xx[j] * fl3_fx + 3.0 * pa2pb_x_x[j] * fl3_fx + 2.25 * pb_xx[j] * fl3_fx + 0.75 * pa_xxyy[j] * fl2_fx + 6.0 * pa2pb_xyy_x[j] * fl2_fx + 4.5 * pa2pb_yy_xx[j] * fl2_fx + 1.5 * pa2pb_xx_xx[j] * fl2_fx + 2.0 * pa2pb_x_xxx[j] * fl2_fx + 3.0 * pa2pb_xxyy_xx[j] * fl1_fx + 4.0 * pa2pb_xyy_xxx[j] * fl1_fx + 0.25 * pb_xxxx[j] * fl2_fx + 0.5 * pa2pb_xx_xxxx[j] * fl1_fx + 0.5 * pa2pb_yy_xxxx[j] * fl1_fx + pa2pb_xxyy_xxxx[j]);

                t_xxyy_xxxx[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 4.5 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 18.75 * pa_yy[j] * fl3_fx * fl1_fz - 1.5 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xx[j] * fl1_fz * fl3_fx - 12.0 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 22.5 * pb_xx[j] * fl3_fx * fl1_fz + 9.0 * pa_xxyy[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_xyy_x[j] * fl2_fx * fl1_fz + 54.0 * pa2pb_yy_xx[j] * fl2_fx * fl1_fz - 1.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 4.0 * pa2pb_x_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxyy_xx[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_xx[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxyy_xx[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xyy_xxx[j] * fl1_fx * fl1_fz - pb_xxxx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxxx[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxx[j] * fl2_fx * fl1_fz - pa2pb_yy_xxxx[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxxx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xxxx[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xxxx[j] * fl1_fz);

                t_xxyy_xxxy[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl3_fx + 2.25 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa2pb_x_y[j] * fl3_fx + 1.125 * pb_xy[j] * fl3_fx + 1.5 * pa2pb_xxy_x[j] * fl2_fx + 1.5 * pa2pb_xyy_y[j] * fl2_fx + 3.0 * pa2pb_xy_xx[j] * fl2_fx + 2.25 * pa2pb_yy_xy[j] * fl2_fx + 0.75 * pa2pb_xx_xy[j] * fl2_fx + 1.5 * pa2pb_x_xxy[j] * fl2_fx + 0.5 * pa2pb_y_xxx[j] * fl2_fx + 1.5 * pa2pb_xxyy_xy[j] * fl1_fx + pa2pb_xxy_xxx[j] * fl1_fx + 3.0 * pa2pb_xyy_xxy[j] * fl1_fx + 0.25 * pb_xxxy[j] * fl2_fx + 0.5 * pa2pb_xx_xxxy[j] * fl1_fx + 0.5 * pa2pb_yy_xxxy[j] * fl1_fx + pa2pb_xxyy_xxxy[j]);

                t_xxyy_xxxy[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_xy[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 11.25 * pb_xy[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xxy_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyy_y[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xy_xx[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yy_xy[j] * fl2_fx * fl1_fz - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xx_xy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyy_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxy_xxx[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyy_xxy[j] * fl1_fx * fl1_fz - pb_xxxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxxy[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxy[j] * fl2_fx * fl1_fz - pa2pb_yy_xxxy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxxy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xxxy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xxxy[j] * fl1_fz);

                t_xxyy_xxxz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl3_fx + 1.125 * pb_xz[j] * fl3_fx + 1.5 * pa2pb_xyy_z[j] * fl2_fx + 2.25 * pa2pb_yy_xz[j] * fl2_fx + 0.75 * pa2pb_xx_xz[j] * fl2_fx + 1.5 * pa2pb_x_xxz[j] * fl2_fx + 1.5 * pa2pb_xxyy_xz[j] * fl1_fx + 3.0 * pa2pb_xyy_xxz[j] * fl1_fx + 0.25 * pb_xxxz[j] * fl2_fx + 0.5 * pa2pb_xx_xxxz[j] * fl1_fx + 0.5 * pa2pb_yy_xxxz[j] * fl1_fx + pa2pb_xxyy_xxxz[j]);

                t_xxyy_xxxz[j] += fl_r_0_0 * (-1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 11.25 * pb_xz[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xyy_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yy_xz[j] * fl2_fx * fl1_fz - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xx_xz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyy_xz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyy_xxz[j] * fl1_fx * fl1_fz - pb_xxxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxxz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxz[j] * fl2_fx * fl1_fz - pa2pb_yy_xxxz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxxz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xxxz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xxxz[j] * fl1_fz);

                t_xxyy_xxyy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 1.5 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa_yy[j] * fl3_fx + 1.5 * pa2pb_y_y[j] * fl3_fx + 0.375 * pb_xx[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 0.25 * pa_xxyy[j] * fl2_fx + pa2pb_xxy_y[j] * fl2_fx + 0.75 * pa2pb_xx_xx[j] * fl2_fx + pa2pb_xyy_x[j] * fl2_fx + 4.0 * pa2pb_xy_xy[j] * fl2_fx + 0.75 * pa2pb_yy_yy[j] * fl2_fx + 0.25 * pa2pb_xx_yy[j] * fl2_fx + pa2pb_x_xyy[j] * fl2_fx + 0.25 * pa2pb_yy_xx[j] * fl2_fx + pa2pb_y_xxy[j] * fl2_fx + 0.5 * pa2pb_xxyy_xx[j] * fl1_fx + 0.5 * pa2pb_xxyy_yy[j] * fl1_fx + 2.0 * pa2pb_xxy_xxy[j] * fl1_fx + 2.0 * pa2pb_xyy_xyy[j] * fl1_fx + 0.25 * pb_xxyy[j] * fl2_fx + 0.5 * pa2pb_xx_xxyy[j] * fl1_fx + 0.5 * pa2pb_yy_xxyy[j] * fl1_fx + pa2pb_xxyy_xxyy[j]);

                t_xxyy_xxyy[j] += fl_r_0_0 * (4.5 * fl4_fx * fl1_fz - 0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_xx[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 3.75 * pa_yy[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz - 0.25 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - pb_xx[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxyy[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 3.0 * pa_xxyy[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xxy_y[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xx_xx[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xyy_x[j] * fl2_fx * fl1_fz + 48.0 * pa2pb_xy_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_yy[j] * fl2_fx * fl1_fz - 0.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_xxy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxyy_yy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_yy[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yy_xx[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyy_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxyy_yy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxy_xxy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyy_xyy[j] * fl1_fx * fl1_fz - pb_xxyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxyy[j] * fl1_fz * fl1_fga + 3.0 * pb_xxyy[j] * fl2_fx * fl1_fz - pa2pb_yy_xxyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xxyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xxyy[j] * fl1_fz);

                t_xxyy_xxyz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl3_fx + 0.375 * pb_yz[j] * fl3_fx + 0.5 * pa2pb_xxy_z[j] * fl2_fx + 2.0 * pa2pb_xy_xz[j] * fl2_fx + 0.75 * pa2pb_yy_yz[j] * fl2_fx + 0.25 * pa2pb_xx_yz[j] * fl2_fx + pa2pb_x_xyz[j] * fl2_fx + 0.5 * pa2pb_y_xxz[j] * fl2_fx + 0.5 * pa2pb_xxyy_yz[j] * fl1_fx + pa2pb_xxy_xxz[j] * fl1_fx + 2.0 * pa2pb_xyy_xyz[j] * fl1_fx + 0.25 * pb_xxyz[j] * fl2_fx + 0.5 * pa2pb_xx_xxyz[j] * fl1_fx + 0.5 * pa2pb_yy_xxyz[j] * fl1_fx + pa2pb_xxyy_xxyz[j]);

                t_xxyy_xxyz[j] += fl_r_0_0 * (7.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_yz[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_xxy_z[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_xy_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_yz[j] * fl2_fx * fl1_fz - 0.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_yz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_yz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyy_yz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxy_xxz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyy_xyz[j] * fl1_fx * fl1_fz - pb_xxyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxyz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxyz[j] * fl2_fx * fl1_fz - pa2pb_yy_xxyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xxyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_50_55(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_xx_xxzz = pa2pbDistances.data(1156 * idx + 126);

            auto pa2pb_xx_xyyy = pa2pbDistances.data(1156 * idx + 127);

            auto pa2pb_xx_xyyz = pa2pbDistances.data(1156 * idx + 128);

            auto pa2pb_xx_xyzz = pa2pbDistances.data(1156 * idx + 129);

            auto pa2pb_xx_xzzz = pa2pbDistances.data(1156 * idx + 130);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_yy_xxzz = pa2pbDistances.data(1156 * idx + 228);

            auto pa2pb_yy_xyyy = pa2pbDistances.data(1156 * idx + 229);

            auto pa2pb_yy_xyyz = pa2pbDistances.data(1156 * idx + 230);

            auto pa2pb_yy_xyzz = pa2pbDistances.data(1156 * idx + 231);

            auto pa2pb_yy_xzzz = pa2pbDistances.data(1156 * idx + 232);

            auto pa2pb_xxy_x = pa2pbDistances.data(1156 * idx + 340);

            auto pa2pb_xxy_xyy = pa2pbDistances.data(1156 * idx + 352);

            auto pa2pb_xxy_xyz = pa2pbDistances.data(1156 * idx + 353);

            auto pa2pb_xxy_xzz = pa2pbDistances.data(1156 * idx + 354);

            auto pa2pb_xyy_x = pa2pbDistances.data(1156 * idx + 408);

            auto pa2pb_xyy_y = pa2pbDistances.data(1156 * idx + 409);

            auto pa2pb_xyy_z = pa2pbDistances.data(1156 * idx + 410);

            auto pa2pb_xyy_xzz = pa2pbDistances.data(1156 * idx + 422);

            auto pa2pb_xyy_yyy = pa2pbDistances.data(1156 * idx + 423);

            auto pa2pb_xyy_yyz = pa2pbDistances.data(1156 * idx + 424);

            auto pa2pb_xyy_yzz = pa2pbDistances.data(1156 * idx + 425);

            auto pa2pb_xyy_zzz = pa2pbDistances.data(1156 * idx + 426);

            auto pa2pb_xxyy_xx = pa2pbDistances.data(1156 * idx + 751);

            auto pa2pb_xxyy_xy = pa2pbDistances.data(1156 * idx + 752);

            auto pa2pb_xxyy_xz = pa2pbDistances.data(1156 * idx + 753);

            auto pa2pb_xxyy_zz = pa2pbDistances.data(1156 * idx + 756);

            auto pa2pb_xxyy_xxzz = pa2pbDistances.data(1156 * idx + 772);

            auto pa2pb_xxyy_xyyy = pa2pbDistances.data(1156 * idx + 773);

            auto pa2pb_xxyy_xyyz = pa2pbDistances.data(1156 * idx + 774);

            auto pa2pb_xxyy_xyzz = pa2pbDistances.data(1156 * idx + 775);

            auto pa2pb_xxyy_xzzz = pa2pbDistances.data(1156 * idx + 776);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyy_xxzz = primBuffer.data(225 * idx + 50);

            auto t_xxyy_xyyy = primBuffer.data(225 * idx + 51);

            auto t_xxyy_xyyz = primBuffer.data(225 * idx + 52);

            auto t_xxyy_xyzz = primBuffer.data(225 * idx + 53);

            auto t_xxyy_xzzz = primBuffer.data(225 * idx + 54);

            // Batch of Integrals (50,55)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_yyy, \
                                     pa2pb_x_yyz, pa2pb_x_yzz, pa2pb_x_z, pa2pb_x_zzz, pa2pb_xx_xx, pa2pb_xx_xxzz, \
                                     pa2pb_xx_xy, pa2pb_xx_xyyy, pa2pb_xx_xyyz, pa2pb_xx_xyzz, pa2pb_xx_xz, \
                                     pa2pb_xx_xzzz, pa2pb_xx_zz, pa2pb_xxy_x, pa2pb_xxy_xyy, pa2pb_xxy_xyz, \
                                     pa2pb_xxy_xzz, pa2pb_xxyy_xx, pa2pb_xxyy_xxzz, pa2pb_xxyy_xy, pa2pb_xxyy_xyyy, \
                                     pa2pb_xxyy_xyyz, pa2pb_xxyy_xyzz, pa2pb_xxyy_xz, pa2pb_xxyy_xzzz, pa2pb_xxyy_zz, \
                                     pa2pb_xy_yy, pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_xyy_x, pa2pb_xyy_xzz, pa2pb_xyy_y, \
                                     pa2pb_xyy_yyy, pa2pb_xyy_yyz, pa2pb_xyy_yzz, pa2pb_xyy_z, pa2pb_xyy_zzz, pa2pb_y_x, \
                                     pa2pb_y_xyy, pa2pb_y_xyz, pa2pb_y_xzz, pa2pb_yy_xx, pa2pb_yy_xxzz, pa2pb_yy_xy, \
                                     pa2pb_yy_xyyy, pa2pb_yy_xyyz, pa2pb_yy_xyzz, pa2pb_yy_xz, pa2pb_yy_xzzz, \
                                     pa2pb_yy_zz, pa_xx, pa_xxyy, pa_xy, pa_yy, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, \
                                     pb_xz, pb_xzzz, pb_zz, r_0_0, s_0_0, t_xxyy_xxzz, t_xxyy_xyyy, t_xxyy_xyyz, \
                                     t_xxyy_xyzz, t_xxyy_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxyy_xxzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.125 * pa_xx[j] * fl3_fx + 0.5 * pa2pb_x_x[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 0.25 * pa_xxyy[j] * fl2_fx + pa2pb_xyy_x[j] * fl2_fx + 0.75 * pa2pb_yy_zz[j] * fl2_fx + 0.125 * pb_xx[j] * fl3_fx + 0.25 * pa2pb_xx_xx[j] * fl2_fx + 0.25 * pa2pb_xx_zz[j] * fl2_fx + pa2pb_x_xzz[j] * fl2_fx + 0.25 * pa2pb_yy_xx[j] * fl2_fx + 0.5 * pa2pb_xxyy_xx[j] * fl1_fx + 0.5 * pa2pb_xxyy_zz[j] * fl1_fx + 2.0 * pa2pb_xyy_xzz[j] * fl1_fx + 0.25 * pb_xxzz[j] * fl2_fx + 0.5 * pa2pb_xx_xxzz[j] * fl1_fx + 0.5 * pa2pb_yy_xxzz[j] * fl1_fx + pa2pb_xxyy_xxzz[j]);

                t_xxyy_xxzz[j] += fl_r_0_0 * (-0.5 * fl3_fx * fl1_fz * fl1_fgb - 0.5 * fl3_fx * fl1_fz * fl1_fga - pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_yy[j] * fl3_fx * fl1_fz - 0.5 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxyy[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_xx[j] * fl1_fz * fl3_fx - 2.0 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz + 3.0 * pa_xxyy[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xyy_x[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_zz[j] * fl2_fx * fl1_fz - 0.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_xx[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxyy_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_xx[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_xx_zz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yy_xx[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyy_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxyy_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyy_xzz[j] * fl1_fx * fl1_fz - pb_xxzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxzz[j] * fl2_fx * fl1_fz - pa2pb_yy_xxzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xxzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xxzz[j] * fl1_fz);

                t_xxyy_xyyy[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl3_fx + 2.25 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa2pb_y_x[j] * fl3_fx + 1.125 * pb_xy[j] * fl3_fx + 1.5 * pa2pb_xxy_x[j] * fl2_fx + 2.25 * pa2pb_xx_xy[j] * fl2_fx + 1.5 * pa2pb_xyy_y[j] * fl2_fx + 3.0 * pa2pb_xy_yy[j] * fl2_fx + 0.5 * pa2pb_x_yyy[j] * fl2_fx + 0.75 * pa2pb_yy_xy[j] * fl2_fx + 1.5 * pa2pb_y_xyy[j] * fl2_fx + 1.5 * pa2pb_xxyy_xy[j] * fl1_fx + 3.0 * pa2pb_xxy_xyy[j] * fl1_fx + pa2pb_xyy_yyy[j] * fl1_fx + 0.25 * pb_xyyy[j] * fl2_fx + 0.5 * pa2pb_xx_xyyy[j] * fl1_fx + 0.5 * pa2pb_yy_xyyy[j] * fl1_fx + pa2pb_xxyy_xyyy[j]);

                t_xxyy_xyyy[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_xy[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_xy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 11.25 * pb_xy[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xxy_x[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xx_xy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyy_y[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xy_yy[j] * fl2_fx * fl1_fz - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_xy[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_xy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyy_xy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxy_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_yyy[j] * fl1_fx * fl1_fz - pb_xyyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyyy[j] * fl1_fz * fl1_fga + 3.0 * pb_xyyy[j] * fl2_fx * fl1_fz - pa2pb_yy_xyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xyyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xyyy[j] * fl1_fz);

                t_xxyy_xyyz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl3_fx + 0.375 * pb_xz[j] * fl3_fx + 0.75 * pa2pb_xx_xz[j] * fl2_fx + 0.5 * pa2pb_xyy_z[j] * fl2_fx + 2.0 * pa2pb_xy_yz[j] * fl2_fx + 0.5 * pa2pb_x_yyz[j] * fl2_fx + 0.25 * pa2pb_yy_xz[j] * fl2_fx + pa2pb_y_xyz[j] * fl2_fx + 0.5 * pa2pb_xxyy_xz[j] * fl1_fx + 2.0 * pa2pb_xxy_xyz[j] * fl1_fx + pa2pb_xyy_yyz[j] * fl1_fx + 0.25 * pb_xyyz[j] * fl2_fx + 0.5 * pa2pb_xx_xyyz[j] * fl1_fx + 0.5 * pa2pb_yy_xyyz[j] * fl1_fx + pa2pb_xxyy_xyyz[j]);

                t_xxyy_xyyz[j] += fl_r_0_0 * (7.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - pb_xz[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xx_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyy_z[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_xy_yz[j] * fl2_fx * fl1_fz - 0.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_xyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_xz[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yy_xz[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyy_xz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxy_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_yyz[j] * fl1_fx * fl1_fz - pb_xyyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyyz[j] * fl1_fz * fl1_fga + 3.0 * pb_xyyz[j] * fl2_fx * fl1_fz - pa2pb_yy_xyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xyyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xyyz[j] * fl1_fz);

                t_xxyy_xyzz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl3_fx + 0.25 * pa2pb_x_y[j] * fl3_fx + 0.25 * pa2pb_y_x[j] * fl3_fx + 0.5 * pa2pb_xxy_x[j] * fl2_fx + 0.5 * pa2pb_xyy_y[j] * fl2_fx + pa2pb_xy_zz[j] * fl2_fx + 0.125 * pb_xy[j] * fl3_fx + 0.25 * pa2pb_xx_xy[j] * fl2_fx + 0.5 * pa2pb_x_yzz[j] * fl2_fx + 0.25 * pa2pb_yy_xy[j] * fl2_fx + 0.5 * pa2pb_y_xzz[j] * fl2_fx + 0.5 * pa2pb_xxyy_xy[j] * fl1_fx + pa2pb_xxy_xzz[j] * fl1_fx + pa2pb_xyy_yzz[j] * fl1_fx + 0.25 * pb_xyzz[j] * fl2_fx + 0.5 * pa2pb_xx_xyzz[j] * fl1_fx + 0.5 * pa2pb_yy_xyzz[j] * fl1_fx + pa2pb_xxyy_xyzz[j]);

                t_xxyy_xyzz[j] += fl_r_0_0 * (-pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 5.0 * pa_xy[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_xxy_x[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xyy_y[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xy_zz[j] * fl2_fx * fl1_fz - 0.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_xy[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xzz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_xy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_xy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yy_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyy_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxy_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_yzz[j] * fl1_fx * fl1_fz - pb_xyzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xyzz[j] * fl2_fx * fl1_fz - pa2pb_yy_xyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xyzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xyzz[j] * fl1_fz);

                t_xxyy_xzzz[j] = fl_s_0_0 * (0.75 * pa2pb_x_z[j] * fl3_fx + 1.5 * pa2pb_xyy_z[j] * fl2_fx + 0.375 * pb_xz[j] * fl3_fx + 0.75 * pa2pb_xx_xz[j] * fl2_fx + 0.5 * pa2pb_x_zzz[j] * fl2_fx + 0.75 * pa2pb_yy_xz[j] * fl2_fx + 1.5 * pa2pb_xxyy_xz[j] * fl1_fx + pa2pb_xyy_zzz[j] * fl1_fx + 0.25 * pb_xzzz[j] * fl2_fx + 0.5 * pa2pb_xx_xzzz[j] * fl1_fx + 0.5 * pa2pb_yy_xzzz[j] * fl1_fx + pa2pb_xxyy_xzzz[j]);

                t_xxyy_xzzz[j] += fl_r_0_0 * (-1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xyy_z[j] * fl2_fx * fl1_fz - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xz[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xx_xz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_xz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyy_xz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_zzz[j] * fl1_fx * fl1_fz - pb_xzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xzzz[j] * fl2_fx * fl1_fz - pa2pb_yy_xzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_xzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_55_60(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_xx_yyyy = pa2pbDistances.data(1156 * idx + 131);

            auto pa2pb_xx_yyyz = pa2pbDistances.data(1156 * idx + 132);

            auto pa2pb_xx_yyzz = pa2pbDistances.data(1156 * idx + 133);

            auto pa2pb_xx_yzzz = pa2pbDistances.data(1156 * idx + 134);

            auto pa2pb_xx_zzzz = pa2pbDistances.data(1156 * idx + 135);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_yy_yyyy = pa2pbDistances.data(1156 * idx + 233);

            auto pa2pb_yy_yyyz = pa2pbDistances.data(1156 * idx + 234);

            auto pa2pb_yy_yyzz = pa2pbDistances.data(1156 * idx + 235);

            auto pa2pb_yy_yzzz = pa2pbDistances.data(1156 * idx + 236);

            auto pa2pb_yy_zzzz = pa2pbDistances.data(1156 * idx + 237);

            auto pa2pb_xxy_y = pa2pbDistances.data(1156 * idx + 341);

            auto pa2pb_xxy_z = pa2pbDistances.data(1156 * idx + 342);

            auto pa2pb_xxy_yyy = pa2pbDistances.data(1156 * idx + 355);

            auto pa2pb_xxy_yyz = pa2pbDistances.data(1156 * idx + 356);

            auto pa2pb_xxy_yzz = pa2pbDistances.data(1156 * idx + 357);

            auto pa2pb_xxy_zzz = pa2pbDistances.data(1156 * idx + 358);

            auto pa2pb_xxyy_yy = pa2pbDistances.data(1156 * idx + 754);

            auto pa2pb_xxyy_yz = pa2pbDistances.data(1156 * idx + 755);

            auto pa2pb_xxyy_zz = pa2pbDistances.data(1156 * idx + 756);

            auto pa2pb_xxyy_yyyy = pa2pbDistances.data(1156 * idx + 777);

            auto pa2pb_xxyy_yyyz = pa2pbDistances.data(1156 * idx + 778);

            auto pa2pb_xxyy_yyzz = pa2pbDistances.data(1156 * idx + 779);

            auto pa2pb_xxyy_yzzz = pa2pbDistances.data(1156 * idx + 780);

            auto pa2pb_xxyy_zzzz = pa2pbDistances.data(1156 * idx + 781);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyy_yyyy = primBuffer.data(225 * idx + 55);

            auto t_xxyy_yyyz = primBuffer.data(225 * idx + 56);

            auto t_xxyy_yyzz = primBuffer.data(225 * idx + 57);

            auto t_xxyy_yzzz = primBuffer.data(225 * idx + 58);

            auto t_xxyy_zzzz = primBuffer.data(225 * idx + 59);

            // Batch of Integrals (55,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_xx_yy, pa2pb_xx_yyyy, pa2pb_xx_yyyz, pa2pb_xx_yyzz, \
                                     pa2pb_xx_yz, pa2pb_xx_yzzz, pa2pb_xx_zz, pa2pb_xx_zzzz, pa2pb_xxy_y, \
                                     pa2pb_xxy_yyy, pa2pb_xxy_yyz, pa2pb_xxy_yzz, pa2pb_xxy_z, pa2pb_xxy_zzz, \
                                     pa2pb_xxyy_yy, pa2pb_xxyy_yyyy, pa2pb_xxyy_yyyz, pa2pb_xxyy_yyzz, pa2pb_xxyy_yz, \
                                     pa2pb_xxyy_yzzz, pa2pb_xxyy_zz, pa2pb_xxyy_zzzz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, \
                                     pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, pa2pb_yy_yy, pa2pb_yy_yyyy, pa2pb_yy_yyyz, \
                                     pa2pb_yy_yyzz, pa2pb_yy_yz, pa2pb_yy_yzzz, pa2pb_yy_zz, pa2pb_yy_zzzz, pa_xx, pa_xxyy, \
                                     pa_yy, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, r_0_0, s_0_0, \
                                     t_xxyy_yyyy, t_xxyy_yyyz, t_xxyy_yyzz, t_xxyy_yzzz, t_xxyy_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxyy_yyyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_xx[j] * fl3_fx + 0.375 * pa_yy[j] * fl3_fx + 3.0 * pa2pb_y_y[j] * fl3_fx + 2.25 * pb_yy[j] * fl3_fx + 0.75 * pa_xxyy[j] * fl2_fx + 6.0 * pa2pb_xxy_y[j] * fl2_fx + 4.5 * pa2pb_xx_yy[j] * fl2_fx + 1.5 * pa2pb_yy_yy[j] * fl2_fx + 2.0 * pa2pb_y_yyy[j] * fl2_fx + 3.0 * pa2pb_xxyy_yy[j] * fl1_fx + 4.0 * pa2pb_xxy_yyy[j] * fl1_fx + 0.25 * pb_yyyy[j] * fl2_fx + 0.5 * pa2pb_xx_yyyy[j] * fl1_fx + 0.5 * pa2pb_yy_yyyy[j] * fl1_fx + pa2pb_xxyy_yyyy[j]);

                t_xxyy_yyyy[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl1_fz * fl1_fga * fl3_fx - 4.5 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 18.75 * pa_xx[j] * fl3_fx * fl1_fz - 0.75 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pb_yy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxyy[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_yy[j] * fl3_fx * fl1_fz + 30.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 22.5 * pb_yy[j] * fl3_fx * fl1_fz + 9.0 * pa_xxyy[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_xxy_y[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_xx_yy[j] * fl2_fx * fl1_fz - 1.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 4.0 * pa2pb_y_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxyy_yy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yy_yy[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxyy_yy[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xxy_yyy[j] * fl1_fz * fl1_fx - pb_yyyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyyy[j] * fl1_fz * fl1_fga + 3.0 * pb_yyyy[j] * fl2_fx * fl1_fz - pa2pb_yy_yyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_yyyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_yyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_yyyy[j] * fl1_fz);

                t_xxyy_yyyz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl3_fx + 1.125 * pb_yz[j] * fl3_fx + 1.5 * pa2pb_xxy_z[j] * fl2_fx + 2.25 * pa2pb_xx_yz[j] * fl2_fx + 0.75 * pa2pb_yy_yz[j] * fl2_fx + 1.5 * pa2pb_y_yyz[j] * fl2_fx + 1.5 * pa2pb_xxyy_yz[j] * fl1_fx + 3.0 * pa2pb_xxy_yyz[j] * fl1_fx + 0.25 * pb_yyyz[j] * fl2_fx + 0.5 * pa2pb_xx_yyyz[j] * fl1_fx + 0.5 * pa2pb_yy_yyyz[j] * fl1_fx + pa2pb_xxyy_yyyz[j]);

                t_xxyy_yyyz[j] += fl_r_0_0 * (-1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_yz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 11.25 * pb_yz[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xxy_z[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xx_yz[j] * fl2_fx * fl1_fz - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_yyz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yy_yz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyy_yz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxy_yyz[j] * fl1_fz * fl1_fx - pb_yyyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyyz[j] * fl1_fz * fl1_fga + 3.0 * pb_yyyz[j] * fl2_fx * fl1_fz - pa2pb_yy_yyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_yyyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_yyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_yyyz[j] * fl1_fz);

                t_xxyy_yyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.125 * pa_yy[j] * fl3_fx + 0.5 * pa2pb_y_y[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 0.25 * pa_xxyy[j] * fl2_fx + pa2pb_xxy_y[j] * fl2_fx + 0.75 * pa2pb_xx_zz[j] * fl2_fx + 0.125 * pb_yy[j] * fl3_fx + 0.25 * pa2pb_xx_yy[j] * fl2_fx + 0.25 * pa2pb_yy_yy[j] * fl2_fx + 0.25 * pa2pb_yy_zz[j] * fl2_fx + pa2pb_y_yzz[j] * fl2_fx + 0.5 * pa2pb_xxyy_yy[j] * fl1_fx + 0.5 * pa2pb_xxyy_zz[j] * fl1_fx + 2.0 * pa2pb_xxy_yzz[j] * fl1_fx + 0.25 * pb_yyzz[j] * fl2_fx + 0.5 * pa2pb_xx_yyzz[j] * fl1_fx + 0.5 * pa2pb_yy_yyzz[j] * fl1_fx + pa2pb_xxyy_yyzz[j]);

                t_xxyy_yyzz[j] += fl_r_0_0 * (-0.5 * fl3_fx * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * fl3_fx - pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_xx[j] * fl3_fx * fl1_fz - 0.25 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - pb_zz[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxyy[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_yy[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz + 3.0 * pa_xxyy[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xxy_y[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xx_zz[j] * fl2_fx * fl1_fz - 0.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_yy[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_yzz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxyy_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_yy[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_yy_yy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yy_zz[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyy_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxyy_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxy_yzz[j] * fl1_fz * fl1_fx - pb_yyzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyzz[j] * fl1_fz * fl1_fga + 3.0 * pb_yyzz[j] * fl2_fx * fl1_fz - pa2pb_yy_yyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_yyzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_yyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_yyzz[j] * fl1_fz);

                t_xxyy_yzzz[j] = fl_s_0_0 * (0.75 * pa2pb_y_z[j] * fl3_fx + 1.5 * pa2pb_xxy_z[j] * fl2_fx + 0.375 * pb_yz[j] * fl3_fx + 0.75 * pa2pb_xx_yz[j] * fl2_fx + 0.75 * pa2pb_yy_yz[j] * fl2_fx + 0.5 * pa2pb_y_zzz[j] * fl2_fx + 1.5 * pa2pb_xxyy_yz[j] * fl1_fx + pa2pb_xxy_zzz[j] * fl1_fx + 0.25 * pb_yzzz[j] * fl2_fx + 0.5 * pa2pb_xx_yzzz[j] * fl1_fx + 0.5 * pa2pb_yy_yzzz[j] * fl1_fx + pa2pb_xxyy_yzzz[j]);

                t_xxyy_yzzz[j] += fl_r_0_0 * (-1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xxy_z[j] * fl1_fz * fl2_fx - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_yz[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xx_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yy_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyy_yz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxy_zzz[j] * fl1_fz * fl1_fx - pb_yzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_yzzz[j] * fl2_fx * fl1_fz - pa2pb_yy_yzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_yzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_yzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_yzzz[j] * fl1_fz);

                t_xxyy_zzzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * pa_yy[j] * fl3_fx + 0.75 * pa_xxyy[j] * fl2_fx + 0.75 * pb_zz[j] * fl3_fx + 1.5 * pa2pb_xx_zz[j] * fl2_fx + 1.5 * pa2pb_yy_zz[j] * fl2_fx + 3.0 * pa2pb_xxyy_zz[j] * fl1_fx + 0.25 * pb_zzzz[j] * fl2_fx + 0.5 * pa2pb_xx_zzzz[j] * fl1_fx + 0.5 * pa2pb_yy_zzzz[j] * fl1_fx + pa2pb_xxyy_zzzz[j]);

                t_xxyy_zzzz[j] += fl_r_0_0 * (-0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - 1.5 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz - 0.75 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xx[j] * fl1_fz * fl3_fx + 3.75 * pa_yy[j] * fl3_fx * fl1_fz + 9.0 * pa_xxyy[j] * fl1_fz * fl2_fx - 1.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_zz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxyy_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_zz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yy_zz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxyy_zz[j] * fl1_fz * fl1_fx - pb_zzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_zzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_zzzz[j] * fl2_fx * fl1_fz - pa2pb_yy_zzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_zzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yy_zzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyy_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_60_65(      CMemBlock2D<double>& primBuffer,
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

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_xxxx = pa2pbDistances.data(1156 * idx + 257);

            auto pa2pb_yz_xxxy = pa2pbDistances.data(1156 * idx + 258);

            auto pa2pb_yz_xxxz = pa2pbDistances.data(1156 * idx + 259);

            auto pa2pb_yz_xxyy = pa2pbDistances.data(1156 * idx + 260);

            auto pa2pb_yz_xxyz = pa2pbDistances.data(1156 * idx + 261);

            auto pa2pb_xxy_x = pa2pbDistances.data(1156 * idx + 340);

            auto pa2pb_xxy_y = pa2pbDistances.data(1156 * idx + 341);

            auto pa2pb_xxy_xxx = pa2pbDistances.data(1156 * idx + 349);

            auto pa2pb_xxy_xxy = pa2pbDistances.data(1156 * idx + 350);

            auto pa2pb_xxz_x = pa2pbDistances.data(1156 * idx + 374);

            auto pa2pb_xxz_y = pa2pbDistances.data(1156 * idx + 375);

            auto pa2pb_xxz_z = pa2pbDistances.data(1156 * idx + 376);

            auto pa2pb_xxz_xxx = pa2pbDistances.data(1156 * idx + 383);

            auto pa2pb_xxz_xxy = pa2pbDistances.data(1156 * idx + 384);

            auto pa2pb_xxz_xxz = pa2pbDistances.data(1156 * idx + 385);

            auto pa2pb_xyz_x = pa2pbDistances.data(1156 * idx + 442);

            auto pa2pb_xyz_y = pa2pbDistances.data(1156 * idx + 443);

            auto pa2pb_xyz_z = pa2pbDistances.data(1156 * idx + 444);

            auto pa2pb_xyz_xxx = pa2pbDistances.data(1156 * idx + 451);

            auto pa2pb_xyz_xxy = pa2pbDistances.data(1156 * idx + 452);

            auto pa2pb_xyz_xxz = pa2pbDistances.data(1156 * idx + 453);

            auto pa2pb_xyz_xyy = pa2pbDistances.data(1156 * idx + 454);

            auto pa2pb_xyz_xyz = pa2pbDistances.data(1156 * idx + 455);

            auto pa2pb_xxyz_xx = pa2pbDistances.data(1156 * idx + 785);

            auto pa2pb_xxyz_xy = pa2pbDistances.data(1156 * idx + 786);

            auto pa2pb_xxyz_xz = pa2pbDistances.data(1156 * idx + 787);

            auto pa2pb_xxyz_yy = pa2pbDistances.data(1156 * idx + 788);

            auto pa2pb_xxyz_yz = pa2pbDistances.data(1156 * idx + 789);

            auto pa2pb_xxyz_xxxx = pa2pbDistances.data(1156 * idx + 801);

            auto pa2pb_xxyz_xxxy = pa2pbDistances.data(1156 * idx + 802);

            auto pa2pb_xxyz_xxxz = pa2pbDistances.data(1156 * idx + 803);

            auto pa2pb_xxyz_xxyy = pa2pbDistances.data(1156 * idx + 804);

            auto pa2pb_xxyz_xxyz = pa2pbDistances.data(1156 * idx + 805);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyz_xxxx = primBuffer.data(225 * idx + 60);

            auto t_xxyz_xxxy = primBuffer.data(225 * idx + 61);

            auto t_xxyz_xxxz = primBuffer.data(225 * idx + 62);

            auto t_xxyz_xxyy = primBuffer.data(225 * idx + 63);

            auto t_xxyz_xxyz = primBuffer.data(225 * idx + 64);

            // Batch of Integrals (60,65)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_xx_xx, pa2pb_xxy_x, pa2pb_xxy_xxx, \
                                     pa2pb_xxy_xxy, pa2pb_xxy_y, pa2pb_xxyz_xx, pa2pb_xxyz_xxxx, pa2pb_xxyz_xxxy, \
                                     pa2pb_xxyz_xxxz, pa2pb_xxyz_xxyy, pa2pb_xxyz_xxyz, pa2pb_xxyz_xy, pa2pb_xxyz_xz, \
                                     pa2pb_xxyz_yy, pa2pb_xxyz_yz, pa2pb_xxz_x, pa2pb_xxz_xxx, pa2pb_xxz_xxy, \
                                     pa2pb_xxz_xxz, pa2pb_xxz_y, pa2pb_xxz_z, pa2pb_xy_xx, pa2pb_xy_xy, pa2pb_xyz_x, \
                                     pa2pb_xyz_xxx, pa2pb_xyz_xxy, pa2pb_xyz_xxz, pa2pb_xyz_xyy, pa2pb_xyz_xyz, \
                                     pa2pb_xyz_y, pa2pb_xyz_z, pa2pb_xz_xx, pa2pb_xz_xy, pa2pb_xz_xz, pa2pb_y_x, \
                                     pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_y, pa2pb_yz_xx, pa2pb_yz_xxxx, pa2pb_yz_xxxy, \
                                     pa2pb_yz_xxxz, pa2pb_yz_xxyy, pa2pb_yz_xxyz, pa2pb_yz_xy, pa2pb_yz_xz, pa2pb_yz_yy, \
                                     pa2pb_yz_yz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_xxz, pa2pb_z_y, pa2pb_z_z, \
                                     pa_xx, pa_xxyz, pa_xy, pa_xz, pa_yz, pb_xx, r_0_0, s_0_0, t_xxyz_xxxx, t_xxyz_xxxy, \
                                     t_xxyz_xxxz, t_xxyz_xxyy, t_xxyz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxyz_xxxx[j] = fl_s_0_0 * (1.875 * pa_yz[j] * fl3_fx + 0.75 * pa_xxyz[j] * fl2_fx + 6.0 * pa2pb_xyz_x[j] * fl2_fx + 4.5 * pa2pb_yz_xx[j] * fl2_fx + 3.0 * pa2pb_xxyz_xx[j] * fl1_fx + 4.0 * pa2pb_xyz_xxx[j] * fl1_fx + 0.5 * pa2pb_yz_xxxx[j] * fl1_fx + pa2pb_xxyz_xxxx[j]);

                t_xxyz_xxxx[j] += fl_r_0_0 * (-4.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 18.75 * pa_yz[j] * fl3_fx * fl1_fz - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxyz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_xxyz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_xyz_x[j] * fl2_fx * fl1_fz + 54.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxyz_xx[j] * fl1_fz * fl1_fgb + 42.0 * pa2pb_xxyz_xx[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xyz_xxx[j] * fl1_fx * fl1_fz - pa2pb_yz_xxxx[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xxxx[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xxxx[j] * fl1_fz);

                t_xxyz_xxxy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl3_fx + 1.125 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa2pb_xxz_x[j] * fl2_fx + 1.5 * pa2pb_xyz_y[j] * fl2_fx + 1.5 * pa2pb_xz_xx[j] * fl2_fx + 2.25 * pa2pb_yz_xy[j] * fl2_fx + 0.25 * pa2pb_z_xxx[j] * fl2_fx + 1.5 * pa2pb_xxyz_xy[j] * fl1_fx + 0.5 * pa2pb_xxz_xxx[j] * fl1_fx + 3.0 * pa2pb_xyz_xxy[j] * fl1_fx + 0.5 * pa2pb_yz_xxxy[j] * fl1_fx + pa2pb_xxyz_xxxy[j]);

                t_xxyz_xxxy[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_xz[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa2pb_xxz_x[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyz_y[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_xy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyz_xy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxz_xxx[j] * fl1_fx * fl1_fz + 42.0 * pa2pb_xyz_xxy[j] * fl1_fx * fl1_fz - pa2pb_yz_xxxy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xxxy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xxxy[j] * fl1_fz);

                t_xxyz_xxxz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl3_fx + 1.125 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa2pb_xxy_x[j] * fl2_fx + 1.5 * pa2pb_xyz_z[j] * fl2_fx + 1.5 * pa2pb_xy_xx[j] * fl2_fx + 2.25 * pa2pb_yz_xz[j] * fl2_fx + 0.25 * pa2pb_y_xxx[j] * fl2_fx + 1.5 * pa2pb_xxyz_xz[j] * fl1_fx + 0.5 * pa2pb_xxy_xxx[j] * fl1_fx + 3.0 * pa2pb_xyz_xxz[j] * fl1_fx + 0.5 * pa2pb_yz_xxxz[j] * fl1_fx + pa2pb_xxyz_xxxz[j]);

                t_xxyz_xxxz[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_xy[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa2pb_xxy_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyz_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xy_xx[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_xz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyz_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxy_xxx[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyz_xxz[j] * fl1_fx * fl1_fz - pa2pb_yz_xxxz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xxxz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xxxz[j] * fl1_fz);

                t_xxyz_xxyy[j] = fl_s_0_0 * (0.375 * pa_yz[j] * fl3_fx + 0.75 * pa2pb_z_y[j] * fl3_fx + 0.25 * pa_xxyz[j] * fl2_fx + 0.5 * pa2pb_xxz_y[j] * fl2_fx + pa2pb_xyz_x[j] * fl2_fx + 2.0 * pa2pb_xz_xy[j] * fl2_fx + 0.75 * pa2pb_yz_yy[j] * fl2_fx + 0.25 * pa2pb_yz_xx[j] * fl2_fx + 0.5 * pa2pb_z_xxy[j] * fl2_fx + 0.5 * pa2pb_xxyz_xx[j] * fl1_fx + 0.5 * pa2pb_xxyz_yy[j] * fl1_fx + pa2pb_xxz_xxy[j] * fl1_fx + 2.0 * pa2pb_xyz_xyy[j] * fl1_fx + 0.5 * pa2pb_yz_xxyy[j] * fl1_fx + pa2pb_xxyz_xxyy[j]);

                t_xxyz_xxyy[j] += fl_r_0_0 * (-pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_yz[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_xxyz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxz_y[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xyz_x[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xxy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxyz_yy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxyz_yy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxz_xxy[j] * fl1_fx * fl1_fz + 28.0 * pa2pb_xyz_xyy[j] * fl1_fx * fl1_fz - pa2pb_yz_xxyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xxyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xxyy[j] * fl1_fz);

                t_xxyz_xxyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.125 * pa_xx[j] * fl3_fx + 0.5 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.125 * pb_xx[j] * fl3_fx + 0.25 * pa2pb_xxy_y[j] * fl2_fx + 0.25 * pa2pb_xxz_z[j] * fl2_fx + 0.25 * pa2pb_xx_xx[j] * fl2_fx + pa2pb_xy_xy[j] * fl2_fx + pa2pb_xz_xz[j] * fl2_fx + 0.75 * pa2pb_yz_yz[j] * fl2_fx + 0.25 * pa2pb_y_xxy[j] * fl2_fx + 0.25 * pa2pb_z_xxz[j] * fl2_fx + 0.5 * pa2pb_xxyz_yz[j] * fl1_fx + 0.5 * pa2pb_xxy_xxy[j] * fl1_fx + 0.5 * pa2pb_xxz_xxz[j] * fl1_fx + 2.0 * pa2pb_xyz_xyz[j] * fl1_fx + 0.5 * pa2pb_yz_xxyz[j] * fl1_fx + pa2pb_xxyz_xxyz[j]);

                t_xxyz_xxyz[j] += fl_r_0_0 * (1.5 * fl4_fx * fl1_fz - 0.125 * fl3_fx * fl1_fz * fl1_fgb - 0.125 * fl1_fz * fl1_fga * fl3_fx - 0.25 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 1.25 * pa_xx[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.25 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.25 * pb_xx[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_xx[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xxy_y[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_xxz_z[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_xx_xx[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xy_xy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_xxy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_yz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyz_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxy_xxy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxz_xxz[j] * fl1_fx * fl1_fz + 28.0 * pa2pb_xyz_xyz[j] * fl1_fx * fl1_fz - pa2pb_yz_xxyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xxyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_65_70(      CMemBlock2D<double>& primBuffer,
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

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_yz_xxzz = pa2pbDistances.data(1156 * idx + 262);

            auto pa2pb_yz_xyyy = pa2pbDistances.data(1156 * idx + 263);

            auto pa2pb_yz_xyyz = pa2pbDistances.data(1156 * idx + 264);

            auto pa2pb_yz_xyzz = pa2pbDistances.data(1156 * idx + 265);

            auto pa2pb_yz_xzzz = pa2pbDistances.data(1156 * idx + 266);

            auto pa2pb_xxy_x = pa2pbDistances.data(1156 * idx + 340);

            auto pa2pb_xxy_z = pa2pbDistances.data(1156 * idx + 342);

            auto pa2pb_xxy_xxz = pa2pbDistances.data(1156 * idx + 351);

            auto pa2pb_xxy_xyy = pa2pbDistances.data(1156 * idx + 352);

            auto pa2pb_xxy_xyz = pa2pbDistances.data(1156 * idx + 353);

            auto pa2pb_xxy_xzz = pa2pbDistances.data(1156 * idx + 354);

            auto pa2pb_xxz_x = pa2pbDistances.data(1156 * idx + 374);

            auto pa2pb_xxz_xyy = pa2pbDistances.data(1156 * idx + 386);

            auto pa2pb_xxz_xyz = pa2pbDistances.data(1156 * idx + 387);

            auto pa2pb_xxz_xzz = pa2pbDistances.data(1156 * idx + 388);

            auto pa2pb_xyz_x = pa2pbDistances.data(1156 * idx + 442);

            auto pa2pb_xyz_y = pa2pbDistances.data(1156 * idx + 443);

            auto pa2pb_xyz_z = pa2pbDistances.data(1156 * idx + 444);

            auto pa2pb_xyz_xzz = pa2pbDistances.data(1156 * idx + 456);

            auto pa2pb_xyz_yyy = pa2pbDistances.data(1156 * idx + 457);

            auto pa2pb_xyz_yyz = pa2pbDistances.data(1156 * idx + 458);

            auto pa2pb_xyz_yzz = pa2pbDistances.data(1156 * idx + 459);

            auto pa2pb_xyz_zzz = pa2pbDistances.data(1156 * idx + 460);

            auto pa2pb_xxyz_xx = pa2pbDistances.data(1156 * idx + 785);

            auto pa2pb_xxyz_xy = pa2pbDistances.data(1156 * idx + 786);

            auto pa2pb_xxyz_xz = pa2pbDistances.data(1156 * idx + 787);

            auto pa2pb_xxyz_zz = pa2pbDistances.data(1156 * idx + 790);

            auto pa2pb_xxyz_xxzz = pa2pbDistances.data(1156 * idx + 806);

            auto pa2pb_xxyz_xyyy = pa2pbDistances.data(1156 * idx + 807);

            auto pa2pb_xxyz_xyyz = pa2pbDistances.data(1156 * idx + 808);

            auto pa2pb_xxyz_xyzz = pa2pbDistances.data(1156 * idx + 809);

            auto pa2pb_xxyz_xzzz = pa2pbDistances.data(1156 * idx + 810);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyz_xxzz = primBuffer.data(225 * idx + 65);

            auto t_xxyz_xyyy = primBuffer.data(225 * idx + 66);

            auto t_xxyz_xyyz = primBuffer.data(225 * idx + 67);

            auto t_xxyz_xyzz = primBuffer.data(225 * idx + 68);

            auto t_xxyz_xzzz = primBuffer.data(225 * idx + 69);

            // Batch of Integrals (65,70)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xy, pa2pb_xx_xz, \
                                     pa2pb_xxy_x, pa2pb_xxy_xxz, pa2pb_xxy_xyy, pa2pb_xxy_xyz, pa2pb_xxy_xzz, \
                                     pa2pb_xxy_z, pa2pb_xxyz_xx, pa2pb_xxyz_xxzz, pa2pb_xxyz_xy, pa2pb_xxyz_xyyy, \
                                     pa2pb_xxyz_xyyz, pa2pb_xxyz_xyzz, pa2pb_xxyz_xz, pa2pb_xxyz_xzzz, pa2pb_xxyz_zz, \
                                     pa2pb_xxz_x, pa2pb_xxz_xyy, pa2pb_xxz_xyz, pa2pb_xxz_xzz, pa2pb_xy_xz, \
                                     pa2pb_xy_yy, pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_xyz_x, pa2pb_xyz_xzz, pa2pb_xyz_y, \
                                     pa2pb_xyz_yyy, pa2pb_xyz_yyz, pa2pb_xyz_yzz, pa2pb_xyz_z, pa2pb_xyz_zzz, \
                                     pa2pb_xz_yy, pa2pb_xz_yz, pa2pb_xz_zz, pa2pb_y_x, pa2pb_y_xxz, pa2pb_y_xyy, \
                                     pa2pb_y_xyz, pa2pb_y_xzz, pa2pb_y_z, pa2pb_yz_xx, pa2pb_yz_xxzz, pa2pb_yz_xy, \
                                     pa2pb_yz_xyyy, pa2pb_yz_xyyz, pa2pb_yz_xyzz, pa2pb_yz_xz, pa2pb_yz_xzzz, \
                                     pa2pb_yz_zz, pa2pb_z_x, pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, pa_xxyz, pa_xy, pa_xz, \
                                     pa_yz, pb_xy, pb_xz, r_0_0, s_0_0, t_xxyz_xxzz, t_xxyz_xyyy, t_xxyz_xyyz, \
                                     t_xxyz_xyzz, t_xxyz_xzzz: VLX_ALIGN)
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

                t_xxyz_xxzz[j] = fl_s_0_0 * (0.375 * pa_yz[j] * fl3_fx + 0.75 * pa2pb_y_z[j] * fl3_fx + 0.25 * pa_xxyz[j] * fl2_fx + 0.5 * pa2pb_xxy_z[j] * fl2_fx + pa2pb_xyz_x[j] * fl2_fx + 2.0 * pa2pb_xy_xz[j] * fl2_fx + 0.75 * pa2pb_yz_zz[j] * fl2_fx + 0.25 * pa2pb_yz_xx[j] * fl2_fx + 0.5 * pa2pb_y_xxz[j] * fl2_fx + 0.5 * pa2pb_xxyz_xx[j] * fl1_fx + 0.5 * pa2pb_xxyz_zz[j] * fl1_fx + pa2pb_xxy_xxz[j] * fl1_fx + 2.0 * pa2pb_xyz_xzz[j] * fl1_fx + 0.5 * pa2pb_yz_xxzz[j] * fl1_fx + pa2pb_xxyz_xxzz[j]);

                t_xxyz_xxzz[j] += fl_r_0_0 * (-pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_yz[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_xxyz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxy_z[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xyz_x[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_xy_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxyz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxyz_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxy_xxz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyz_xzz[j] * fl1_fx * fl1_fz - pa2pb_yz_xxzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xxzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xxzz[j] * fl1_fz);

                t_xxyz_xyyy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl3_fx + 0.375 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa2pb_xxz_x[j] * fl2_fx + 1.5 * pa2pb_xyz_y[j] * fl2_fx + 1.5 * pa2pb_xz_yy[j] * fl2_fx + 0.75 * pa2pb_yz_xy[j] * fl2_fx + 0.75 * pa2pb_z_xyy[j] * fl2_fx + 1.5 * pa2pb_xxyz_xy[j] * fl1_fx + 1.5 * pa2pb_xxz_xyy[j] * fl1_fx + pa2pb_xyz_yyy[j] * fl1_fx + 0.5 * pa2pb_yz_xyyy[j] * fl1_fx + pa2pb_xxyz_xyyy[j]);

                t_xxyz_xyyy[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_xz[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxz_x[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyz_y[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyz_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxz_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_yyy[j] * fl1_fx * fl1_fz - pa2pb_yz_xyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xyyy[j] * fl1_fz);

                t_xxyz_xyyz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl3_fx + 0.5 * pa2pb_x_y[j] * fl3_fx + 0.125 * pa2pb_y_x[j] * fl3_fx + 0.25 * pb_xy[j] * fl3_fx + 0.25 * pa2pb_xxy_x[j] * fl2_fx + 0.5 * pa2pb_xx_xy[j] * fl2_fx + 0.5 * pa2pb_xyz_z[j] * fl2_fx + 0.5 * pa2pb_xy_yy[j] * fl2_fx + pa2pb_xz_yz[j] * fl2_fx + 0.25 * pa2pb_yz_xz[j] * fl2_fx + 0.25 * pa2pb_y_xyy[j] * fl2_fx + 0.5 * pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_xxyz_xz[j] * fl1_fx + 0.5 * pa2pb_xxy_xyy[j] * fl1_fx + pa2pb_xxz_xyz[j] * fl1_fx + pa2pb_xyz_yyz[j] * fl1_fx + 0.5 * pa2pb_yz_xyyz[j] * fl1_fx + pa2pb_xxyz_xyyz[j]);

                t_xxyz_xyyz[j] += fl_r_0_0 * (-0.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 2.5 * pa_xy[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_xy[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 2.5 * pb_xy[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xxy_x[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xx_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyz_z[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xy_yy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_xyy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_xz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyz_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxy_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxz_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_yyz[j] * fl1_fx * fl1_fz - pa2pb_yz_xyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xyyz[j] * fl1_fz);

                t_xxyz_xyzz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl3_fx + 0.5 * pa2pb_x_z[j] * fl3_fx + 0.125 * pa2pb_z_x[j] * fl3_fx + 0.25 * pb_xz[j] * fl3_fx + 0.25 * pa2pb_xxz_x[j] * fl2_fx + 0.5 * pa2pb_xx_xz[j] * fl2_fx + 0.5 * pa2pb_xyz_y[j] * fl2_fx + pa2pb_xy_yz[j] * fl2_fx + 0.5 * pa2pb_xz_zz[j] * fl2_fx + 0.25 * pa2pb_yz_xy[j] * fl2_fx + 0.5 * pa2pb_y_xyz[j] * fl2_fx + 0.25 * pa2pb_z_xzz[j] * fl2_fx + 0.5 * pa2pb_xxyz_xy[j] * fl1_fx + pa2pb_xxy_xyz[j] * fl1_fx + 0.5 * pa2pb_xxz_xzz[j] * fl1_fx + pa2pb_xyz_yzz[j] * fl1_fx + 0.5 * pa2pb_yz_xyzz[j] * fl1_fx + pa2pb_xxyz_xyzz[j]);

                t_xxyz_xyzz[j] += fl_r_0_0 * (-0.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 2.5 * pa_xz[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_xz[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 2.5 * pb_xz[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xxz_x[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xx_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyz_y[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xy_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xyz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_xzz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_xy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyz_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxy_xyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxz_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyz_yzz[j] * fl1_fx * fl1_fz - pa2pb_yz_xyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xyzz[j] * fl1_fz);

                t_xxyz_xzzz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl3_fx + 0.375 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa2pb_xxy_x[j] * fl2_fx + 1.5 * pa2pb_xyz_z[j] * fl2_fx + 1.5 * pa2pb_xy_zz[j] * fl2_fx + 0.75 * pa2pb_yz_xz[j] * fl2_fx + 0.75 * pa2pb_y_xzz[j] * fl2_fx + 1.5 * pa2pb_xxyz_xz[j] * fl1_fx + 1.5 * pa2pb_xxy_xzz[j] * fl1_fx + pa2pb_xyz_zzz[j] * fl1_fx + 0.5 * pa2pb_yz_xzzz[j] * fl1_fx + pa2pb_xxyz_xzzz[j]);

                t_xxyz_xzzz[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_xy[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxy_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyz_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xy_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyz_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxy_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_zzz[j] * fl1_fx * fl1_fz - pa2pb_yz_xzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_xzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_70_75(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_yz_yyyy = pa2pbDistances.data(1156 * idx + 267);

            auto pa2pb_yz_yyyz = pa2pbDistances.data(1156 * idx + 268);

            auto pa2pb_yz_yyzz = pa2pbDistances.data(1156 * idx + 269);

            auto pa2pb_yz_yzzz = pa2pbDistances.data(1156 * idx + 270);

            auto pa2pb_yz_zzzz = pa2pbDistances.data(1156 * idx + 271);

            auto pa2pb_xxy_y = pa2pbDistances.data(1156 * idx + 341);

            auto pa2pb_xxy_z = pa2pbDistances.data(1156 * idx + 342);

            auto pa2pb_xxy_yyy = pa2pbDistances.data(1156 * idx + 355);

            auto pa2pb_xxy_yyz = pa2pbDistances.data(1156 * idx + 356);

            auto pa2pb_xxy_yzz = pa2pbDistances.data(1156 * idx + 357);

            auto pa2pb_xxy_zzz = pa2pbDistances.data(1156 * idx + 358);

            auto pa2pb_xxz_y = pa2pbDistances.data(1156 * idx + 375);

            auto pa2pb_xxz_z = pa2pbDistances.data(1156 * idx + 376);

            auto pa2pb_xxz_yyy = pa2pbDistances.data(1156 * idx + 389);

            auto pa2pb_xxz_yyz = pa2pbDistances.data(1156 * idx + 390);

            auto pa2pb_xxz_yzz = pa2pbDistances.data(1156 * idx + 391);

            auto pa2pb_xxz_zzz = pa2pbDistances.data(1156 * idx + 392);

            auto pa2pb_xxyz_yy = pa2pbDistances.data(1156 * idx + 788);

            auto pa2pb_xxyz_yz = pa2pbDistances.data(1156 * idx + 789);

            auto pa2pb_xxyz_zz = pa2pbDistances.data(1156 * idx + 790);

            auto pa2pb_xxyz_yyyy = pa2pbDistances.data(1156 * idx + 811);

            auto pa2pb_xxyz_yyyz = pa2pbDistances.data(1156 * idx + 812);

            auto pa2pb_xxyz_yyzz = pa2pbDistances.data(1156 * idx + 813);

            auto pa2pb_xxyz_yzzz = pa2pbDistances.data(1156 * idx + 814);

            auto pa2pb_xxyz_zzzz = pa2pbDistances.data(1156 * idx + 815);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyz_yyyy = primBuffer.data(225 * idx + 70);

            auto t_xxyz_yyyz = primBuffer.data(225 * idx + 71);

            auto t_xxyz_yyzz = primBuffer.data(225 * idx + 72);

            auto t_xxyz_yzzz = primBuffer.data(225 * idx + 73);

            auto t_xxyz_zzzz = primBuffer.data(225 * idx + 74);

            // Batch of Integrals (70,75)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xx_zz, pa2pb_xxy_y, \
                                     pa2pb_xxy_yyy, pa2pb_xxy_yyz, pa2pb_xxy_yzz, pa2pb_xxy_z, pa2pb_xxy_zzz, \
                                     pa2pb_xxyz_yy, pa2pb_xxyz_yyyy, pa2pb_xxyz_yyyz, pa2pb_xxyz_yyzz, pa2pb_xxyz_yz, \
                                     pa2pb_xxyz_yzzz, pa2pb_xxyz_zz, pa2pb_xxyz_zzzz, pa2pb_xxz_y, pa2pb_xxz_yyy, \
                                     pa2pb_xxz_yyz, pa2pb_xxz_yzz, pa2pb_xxz_z, pa2pb_xxz_zzz, pa2pb_y_y, pa2pb_y_yyy, \
                                     pa2pb_y_yyz, pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, pa2pb_yz_yy, pa2pb_yz_yyyy, \
                                     pa2pb_yz_yyyz, pa2pb_yz_yyzz, pa2pb_yz_yz, pa2pb_yz_yzzz, pa2pb_yz_zz, \
                                     pa2pb_yz_zzzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, pa2pb_z_z, \
                                     pa2pb_z_zzz, pa_xx, pa_xxyz, pa_yz, pb_yy, pb_yz, pb_zz, r_0_0, s_0_0, t_xxyz_yyyy, \
                                     t_xxyz_yyyz, t_xxyz_yyzz, t_xxyz_yzzz, t_xxyz_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxyz_yyyy[j] = fl_s_0_0 * (0.375 * pa_yz[j] * fl3_fx + 1.5 * pa2pb_z_y[j] * fl3_fx + 0.75 * pa_xxyz[j] * fl2_fx + 3.0 * pa2pb_xxz_y[j] * fl2_fx + 1.5 * pa2pb_yz_yy[j] * fl2_fx + pa2pb_z_yyy[j] * fl2_fx + 3.0 * pa2pb_xxyz_yy[j] * fl1_fx + 2.0 * pa2pb_xxz_yyy[j] * fl1_fx + 0.5 * pa2pb_yz_yyyy[j] * fl1_fx + pa2pb_xxyz_yyyy[j]);

                t_xxyz_yyyy[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_yz[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 9.0 * pa_xxyz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xxz_y[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxyz_yy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxyz_yy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxz_yyy[j] * fl1_fx * fl1_fz - pa2pb_yz_yyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_yyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_yyyy[j] * fl1_fz);

                t_xxyz_yyyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 0.75 * pa2pb_xxy_y[j] * fl2_fx + 0.75 * pa2pb_xxz_z[j] * fl2_fx + 0.75 * pa2pb_xx_yy[j] * fl2_fx + 0.75 * pa2pb_yz_yz[j] * fl2_fx + 0.25 * pa2pb_y_yyy[j] * fl2_fx + 0.75 * pa2pb_z_yyz[j] * fl2_fx + 1.5 * pa2pb_xxyz_yz[j] * fl1_fx + 0.5 * pa2pb_xxy_yyy[j] * fl1_fx + 1.5 * pa2pb_xxz_yyz[j] * fl1_fx + 0.5 * pa2pb_yz_yyyz[j] * fl1_fx + pa2pb_xxyz_yyyz[j]);

                t_xxyz_yyyz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl1_fz * fl1_fga * fl3_fx - 0.75 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_xx[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pb_yy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxy_y[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xxz_z[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xx_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yyz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyz_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxy_yyy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxz_yyz[j] * fl1_fx * fl1_fz - pa2pb_yz_yyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_yyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_yyyz[j] * fl1_fz);

                t_xxyz_yyzz[j] = fl_s_0_0 * (0.125 * pa_yz[j] * fl3_fx + 0.25 * pa2pb_y_z[j] * fl3_fx + 0.25 * pa2pb_z_y[j] * fl3_fx + 0.5 * pb_yz[j] * fl3_fx + 0.25 * pa_xxyz[j] * fl2_fx + 0.5 * pa2pb_xxy_z[j] * fl2_fx + 0.5 * pa2pb_xxz_y[j] * fl2_fx + pa2pb_xx_yz[j] * fl2_fx + 0.25 * pa2pb_yz_yy[j] * fl2_fx + 0.25 * pa2pb_yz_zz[j] * fl2_fx + 0.5 * pa2pb_y_yyz[j] * fl2_fx + 0.5 * pa2pb_z_yzz[j] * fl2_fx + 0.5 * pa2pb_xxyz_yy[j] * fl1_fx + 0.5 * pa2pb_xxyz_zz[j] * fl1_fx + pa2pb_xxy_yyz[j] * fl1_fx + pa2pb_xxz_yzz[j] * fl1_fx + 0.5 * pa2pb_yz_yyzz[j] * fl1_fx + pa2pb_xxyz_yyzz[j]);

                t_xxyz_yyzz[j] += fl_r_0_0 * (-0.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - pb_yz[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_yz[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 5.0 * pb_yz[j] * fl3_fx * fl1_fz + 3.0 * pa_xxyz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxy_z[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxz_y[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xx_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_yyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_yzz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxyz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxyz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxyz_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxy_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxz_yzz[j] * fl1_fx * fl1_fz - pa2pb_yz_yyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_yyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_yyzz[j] * fl1_fz);

                t_xxyz_yzzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 0.75 * pa2pb_xxy_y[j] * fl2_fx + 0.75 * pa2pb_xxz_z[j] * fl2_fx + 0.75 * pa2pb_xx_zz[j] * fl2_fx + 0.75 * pa2pb_yz_yz[j] * fl2_fx + 0.75 * pa2pb_y_yzz[j] * fl2_fx + 0.25 * pa2pb_z_zzz[j] * fl2_fx + 1.5 * pa2pb_xxyz_yz[j] * fl1_fx + 1.5 * pa2pb_xxy_yzz[j] * fl1_fx + 0.5 * pa2pb_xxz_zzz[j] * fl1_fx + 0.5 * pa2pb_yz_yzzz[j] * fl1_fx + pa2pb_xxyz_yzzz[j]);

                t_xxyz_yzzz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl1_fz * fl1_fga * fl3_fx - 0.75 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_xx[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pb_zz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xxy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xxy_y[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xxz_z[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xx_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yzz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxyz_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xxy_yzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxz_zzz[j] * fl1_fx * fl1_fz - pa2pb_yz_yzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_yzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_yzzz[j] * fl1_fz);

                t_xxyz_zzzz[j] = fl_s_0_0 * (0.375 * pa_yz[j] * fl3_fx + 1.5 * pa2pb_y_z[j] * fl3_fx + 0.75 * pa_xxyz[j] * fl2_fx + 3.0 * pa2pb_xxy_z[j] * fl2_fx + 1.5 * pa2pb_yz_zz[j] * fl2_fx + pa2pb_y_zzz[j] * fl2_fx + 3.0 * pa2pb_xxyz_zz[j] * fl1_fx + 2.0 * pa2pb_xxy_zzz[j] * fl1_fx + 0.5 * pa2pb_yz_zzzz[j] * fl1_fx + pa2pb_xxyz_zzzz[j]);

                t_xxyz_zzzz[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xxy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_yz[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 9.0 * pa_xxyz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xxy_z[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxyz_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxyz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxy_zzz[j] * fl1_fz * fl1_fx - pa2pb_yz_zzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yz_zzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxyz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_75_80(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

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

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_xxxx = pa2pbDistances.data(1156 * idx + 121);

            auto pa2pb_xx_xxxy = pa2pbDistances.data(1156 * idx + 122);

            auto pa2pb_xx_xxxz = pa2pbDistances.data(1156 * idx + 123);

            auto pa2pb_xx_xxyy = pa2pbDistances.data(1156 * idx + 124);

            auto pa2pb_xx_xxyz = pa2pbDistances.data(1156 * idx + 125);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_xxxx = pa2pbDistances.data(1156 * idx + 291);

            auto pa2pb_zz_xxxy = pa2pbDistances.data(1156 * idx + 292);

            auto pa2pb_zz_xxxz = pa2pbDistances.data(1156 * idx + 293);

            auto pa2pb_zz_xxyy = pa2pbDistances.data(1156 * idx + 294);

            auto pa2pb_zz_xxyz = pa2pbDistances.data(1156 * idx + 295);

            auto pa2pb_xxz_x = pa2pbDistances.data(1156 * idx + 374);

            auto pa2pb_xxz_y = pa2pbDistances.data(1156 * idx + 375);

            auto pa2pb_xxz_xxx = pa2pbDistances.data(1156 * idx + 383);

            auto pa2pb_xxz_xxy = pa2pbDistances.data(1156 * idx + 384);

            auto pa2pb_xzz_x = pa2pbDistances.data(1156 * idx + 476);

            auto pa2pb_xzz_y = pa2pbDistances.data(1156 * idx + 477);

            auto pa2pb_xzz_z = pa2pbDistances.data(1156 * idx + 478);

            auto pa2pb_xzz_xxx = pa2pbDistances.data(1156 * idx + 485);

            auto pa2pb_xzz_xxy = pa2pbDistances.data(1156 * idx + 486);

            auto pa2pb_xzz_xxz = pa2pbDistances.data(1156 * idx + 487);

            auto pa2pb_xzz_xyy = pa2pbDistances.data(1156 * idx + 488);

            auto pa2pb_xzz_xyz = pa2pbDistances.data(1156 * idx + 489);

            auto pa2pb_xxzz_xx = pa2pbDistances.data(1156 * idx + 819);

            auto pa2pb_xxzz_xy = pa2pbDistances.data(1156 * idx + 820);

            auto pa2pb_xxzz_xz = pa2pbDistances.data(1156 * idx + 821);

            auto pa2pb_xxzz_yy = pa2pbDistances.data(1156 * idx + 822);

            auto pa2pb_xxzz_yz = pa2pbDistances.data(1156 * idx + 823);

            auto pa2pb_xxzz_xxxx = pa2pbDistances.data(1156 * idx + 835);

            auto pa2pb_xxzz_xxxy = pa2pbDistances.data(1156 * idx + 836);

            auto pa2pb_xxzz_xxxz = pa2pbDistances.data(1156 * idx + 837);

            auto pa2pb_xxzz_xxyy = pa2pbDistances.data(1156 * idx + 838);

            auto pa2pb_xxzz_xxyz = pa2pbDistances.data(1156 * idx + 839);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxzz_xxxx = primBuffer.data(225 * idx + 75);

            auto t_xxzz_xxxy = primBuffer.data(225 * idx + 76);

            auto t_xxzz_xxxz = primBuffer.data(225 * idx + 77);

            auto t_xxzz_xxyy = primBuffer.data(225 * idx + 78);

            auto t_xxzz_xxyz = primBuffer.data(225 * idx + 79);

            // Batch of Integrals (75,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_y, pa2pb_x_z, pa2pb_xx_xx, pa2pb_xx_xxxx, \
                                     pa2pb_xx_xxxy, pa2pb_xx_xxxz, pa2pb_xx_xxyy, pa2pb_xx_xxyz, pa2pb_xx_xy, \
                                     pa2pb_xx_xz, pa2pb_xx_yy, pa2pb_xx_yz, pa2pb_xxz_x, pa2pb_xxz_xxx, pa2pb_xxz_xxy, \
                                     pa2pb_xxz_y, pa2pb_xxzz_xx, pa2pb_xxzz_xxxx, pa2pb_xxzz_xxxy, pa2pb_xxzz_xxxz, \
                                     pa2pb_xxzz_xxyy, pa2pb_xxzz_xxyz, pa2pb_xxzz_xy, pa2pb_xxzz_xz, pa2pb_xxzz_yy, \
                                     pa2pb_xxzz_yz, pa2pb_xz_xx, pa2pb_xz_xy, pa2pb_xzz_x, pa2pb_xzz_xxx, pa2pb_xzz_xxy, \
                                     pa2pb_xzz_xxz, pa2pb_xzz_xyy, pa2pb_xzz_xyz, pa2pb_xzz_y, pa2pb_xzz_z, pa2pb_z_x, \
                                     pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_y, pa2pb_zz_xx, pa2pb_zz_xxxx, pa2pb_zz_xxxy, \
                                     pa2pb_zz_xxxz, pa2pb_zz_xxyy, pa2pb_zz_xxyz, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, \
                                     pa2pb_zz_yz, pa_xx, pa_xxzz, pa_xz, pa_zz, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, \
                                     pb_xxyz, pb_xy, pb_xz, pb_yy, pb_yz, r_0_0, s_0_0, t_xxzz_xxxx, t_xxzz_xxxy, \
                                     t_xxzz_xxxz, t_xxzz_xxyy, t_xxzz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxzz_xxxx[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_zz[j] * fl3_fx + 0.375 * pa_xx[j] * fl3_fx + 3.0 * pa2pb_x_x[j] * fl3_fx + 2.25 * pb_xx[j] * fl3_fx + 0.75 * pa_xxzz[j] * fl2_fx + 6.0 * pa2pb_xzz_x[j] * fl2_fx + 4.5 * pa2pb_zz_xx[j] * fl2_fx + 1.5 * pa2pb_xx_xx[j] * fl2_fx + 2.0 * pa2pb_x_xxx[j] * fl2_fx + 3.0 * pa2pb_xxzz_xx[j] * fl1_fx + 4.0 * pa2pb_xzz_xxx[j] * fl1_fx + 0.25 * pb_xxxx[j] * fl2_fx + 0.5 * pa2pb_xx_xxxx[j] * fl1_fx + 0.5 * pa2pb_zz_xxxx[j] * fl1_fx + pa2pb_xxzz_xxxx[j]);

                t_xxzz_xxxx[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 4.5 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 18.75 * pa_zz[j] * fl3_fx * fl1_fz - 1.5 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xx[j] * fl1_fz * fl3_fx - 12.0 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 22.5 * pb_xx[j] * fl3_fx * fl1_fz + 9.0 * pa_xxzz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_xzz_x[j] * fl2_fx * fl1_fz + 54.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz - 1.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 4.0 * pa2pb_x_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxzz_xx[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_xx[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxzz_xx[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xzz_xxx[j] * fl1_fx * fl1_fz - pb_xxxx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxxx[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxx[j] * fl2_fx * fl1_fz - pa2pb_zz_xxxx[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxxx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxxx[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xxxx[j] * fl1_fz);

                t_xxzz_xxxy[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl3_fx + 1.125 * pb_xy[j] * fl3_fx + 1.5 * pa2pb_xzz_y[j] * fl2_fx + 2.25 * pa2pb_zz_xy[j] * fl2_fx + 0.75 * pa2pb_xx_xy[j] * fl2_fx + 1.5 * pa2pb_x_xxy[j] * fl2_fx + 1.5 * pa2pb_xxzz_xy[j] * fl1_fx + 3.0 * pa2pb_xzz_xxy[j] * fl1_fx + 0.25 * pb_xxxy[j] * fl2_fx + 0.5 * pa2pb_xx_xxxy[j] * fl1_fx + 0.5 * pa2pb_zz_xxxy[j] * fl1_fx + pa2pb_xxzz_xxxy[j]);

                t_xxzz_xxxy[j] += fl_r_0_0 * (-1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 11.25 * pb_xy[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xzz_y[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xx_xy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxzz_xy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xzz_xxy[j] * fl1_fx * fl1_fz - pb_xxxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxxy[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxy[j] * fl2_fx * fl1_fz - pa2pb_zz_xxxy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxxy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxxy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xxxy[j] * fl1_fz);

                t_xxzz_xxxz[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl3_fx + 2.25 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa2pb_x_z[j] * fl3_fx + 1.125 * pb_xz[j] * fl3_fx + 1.5 * pa2pb_xxz_x[j] * fl2_fx + 1.5 * pa2pb_xzz_z[j] * fl2_fx + 3.0 * pa2pb_xz_xx[j] * fl2_fx + 2.25 * pa2pb_zz_xz[j] * fl2_fx + 0.75 * pa2pb_xx_xz[j] * fl2_fx + 1.5 * pa2pb_x_xxz[j] * fl2_fx + 0.5 * pa2pb_z_xxx[j] * fl2_fx + 1.5 * pa2pb_xxzz_xz[j] * fl1_fx + pa2pb_xxz_xxx[j] * fl1_fx + 3.0 * pa2pb_xzz_xxz[j] * fl1_fx + 0.25 * pb_xxxz[j] * fl2_fx + 0.5 * pa2pb_xx_xxxz[j] * fl1_fx + 0.5 * pa2pb_zz_xxxz[j] * fl1_fx + pa2pb_xxzz_xxxz[j]);

                t_xxzz_xxxz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_xz[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 11.25 * pb_xz[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xxz_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xzz_z[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xx_xz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxzz_xz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxz_xxx[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xzz_xxz[j] * fl1_fx * fl1_fz - pb_xxxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxxz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxz[j] * fl2_fx * fl1_fz - pa2pb_zz_xxxz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxxz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxxz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xxxz[j] * fl1_fz);

                t_xxzz_xxyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_zz[j] * fl3_fx + 0.125 * pa_xx[j] * fl3_fx + 0.5 * pa2pb_x_x[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 0.25 * pa_xxzz[j] * fl2_fx + pa2pb_xzz_x[j] * fl2_fx + 0.75 * pa2pb_zz_yy[j] * fl2_fx + 0.125 * pb_xx[j] * fl3_fx + 0.25 * pa2pb_xx_xx[j] * fl2_fx + 0.25 * pa2pb_xx_yy[j] * fl2_fx + pa2pb_x_xyy[j] * fl2_fx + 0.25 * pa2pb_zz_xx[j] * fl2_fx + 0.5 * pa2pb_xxzz_xx[j] * fl1_fx + 0.5 * pa2pb_xxzz_yy[j] * fl1_fx + 2.0 * pa2pb_xzz_xyy[j] * fl1_fx + 0.25 * pb_xxyy[j] * fl2_fx + 0.5 * pa2pb_xx_xxyy[j] * fl1_fx + 0.5 * pa2pb_zz_xxyy[j] * fl1_fx + pa2pb_xxzz_xxyy[j]);

                t_xxzz_xxyy[j] += fl_r_0_0 * (-0.5 * fl3_fx * fl1_fz * fl1_fgb - 0.5 * fl3_fx * fl1_fz * fl1_fga - pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_zz[j] * fl3_fx * fl1_fz - 0.5 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxzz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_xx[j] * fl1_fz * fl3_fx - 2.0 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 3.0 * pa_xxzz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xzz_x[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz - 0.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_xx[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxzz_yy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_xx[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_xx_yy[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxzz_yy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xzz_xyy[j] * fl1_fx * fl1_fz - pb_xxyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxyy[j] * fl1_fz * fl1_fga + 3.0 * pb_xxyy[j] * fl2_fx * fl1_fz - pa2pb_zz_xxyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xxyy[j] * fl1_fz);

                t_xxzz_xxyz[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl3_fx + 0.375 * pb_yz[j] * fl3_fx + 0.5 * pa2pb_xxz_y[j] * fl2_fx + 2.0 * pa2pb_xz_xy[j] * fl2_fx + 0.75 * pa2pb_zz_yz[j] * fl2_fx + 0.25 * pa2pb_xx_yz[j] * fl2_fx + pa2pb_x_xyz[j] * fl2_fx + 0.5 * pa2pb_z_xxy[j] * fl2_fx + 0.5 * pa2pb_xxzz_yz[j] * fl1_fx + pa2pb_xxz_xxy[j] * fl1_fx + 2.0 * pa2pb_xzz_xyz[j] * fl1_fx + 0.25 * pb_xxyz[j] * fl2_fx + 0.5 * pa2pb_xx_xxyz[j] * fl1_fx + 0.5 * pa2pb_zz_xxyz[j] * fl1_fx + pa2pb_xxzz_xxyz[j]);

                t_xxzz_xxyz[j] += fl_r_0_0 * (7.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_yz[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_xxz_y[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz - 0.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xxy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_yz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_yz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxzz_yz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxz_xxy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xzz_xyz[j] * fl1_fx * fl1_fz - pb_xxyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxyz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxyz[j] * fl2_fx * fl1_fz - pa2pb_zz_xxyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_80_85(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_xx_xx = pa2pbDistances.data(1156 * idx + 105);

            auto pa2pb_xx_xy = pa2pbDistances.data(1156 * idx + 106);

            auto pa2pb_xx_xz = pa2pbDistances.data(1156 * idx + 107);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_xx_xxzz = pa2pbDistances.data(1156 * idx + 126);

            auto pa2pb_xx_xyyy = pa2pbDistances.data(1156 * idx + 127);

            auto pa2pb_xx_xyyz = pa2pbDistances.data(1156 * idx + 128);

            auto pa2pb_xx_xyzz = pa2pbDistances.data(1156 * idx + 129);

            auto pa2pb_xx_xzzz = pa2pbDistances.data(1156 * idx + 130);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_zz_xxzz = pa2pbDistances.data(1156 * idx + 296);

            auto pa2pb_zz_xyyy = pa2pbDistances.data(1156 * idx + 297);

            auto pa2pb_zz_xyyz = pa2pbDistances.data(1156 * idx + 298);

            auto pa2pb_zz_xyzz = pa2pbDistances.data(1156 * idx + 299);

            auto pa2pb_zz_xzzz = pa2pbDistances.data(1156 * idx + 300);

            auto pa2pb_xxz_x = pa2pbDistances.data(1156 * idx + 374);

            auto pa2pb_xxz_z = pa2pbDistances.data(1156 * idx + 376);

            auto pa2pb_xxz_xxz = pa2pbDistances.data(1156 * idx + 385);

            auto pa2pb_xxz_xyy = pa2pbDistances.data(1156 * idx + 386);

            auto pa2pb_xxz_xyz = pa2pbDistances.data(1156 * idx + 387);

            auto pa2pb_xxz_xzz = pa2pbDistances.data(1156 * idx + 388);

            auto pa2pb_xzz_x = pa2pbDistances.data(1156 * idx + 476);

            auto pa2pb_xzz_y = pa2pbDistances.data(1156 * idx + 477);

            auto pa2pb_xzz_z = pa2pbDistances.data(1156 * idx + 478);

            auto pa2pb_xzz_xzz = pa2pbDistances.data(1156 * idx + 490);

            auto pa2pb_xzz_yyy = pa2pbDistances.data(1156 * idx + 491);

            auto pa2pb_xzz_yyz = pa2pbDistances.data(1156 * idx + 492);

            auto pa2pb_xzz_yzz = pa2pbDistances.data(1156 * idx + 493);

            auto pa2pb_xzz_zzz = pa2pbDistances.data(1156 * idx + 494);

            auto pa2pb_xxzz_xx = pa2pbDistances.data(1156 * idx + 819);

            auto pa2pb_xxzz_xy = pa2pbDistances.data(1156 * idx + 820);

            auto pa2pb_xxzz_xz = pa2pbDistances.data(1156 * idx + 821);

            auto pa2pb_xxzz_zz = pa2pbDistances.data(1156 * idx + 824);

            auto pa2pb_xxzz_xxzz = pa2pbDistances.data(1156 * idx + 840);

            auto pa2pb_xxzz_xyyy = pa2pbDistances.data(1156 * idx + 841);

            auto pa2pb_xxzz_xyyz = pa2pbDistances.data(1156 * idx + 842);

            auto pa2pb_xxzz_xyzz = pa2pbDistances.data(1156 * idx + 843);

            auto pa2pb_xxzz_xzzz = pa2pbDistances.data(1156 * idx + 844);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxzz_xxzz = primBuffer.data(225 * idx + 80);

            auto t_xxzz_xyyy = primBuffer.data(225 * idx + 81);

            auto t_xxzz_xyyz = primBuffer.data(225 * idx + 82);

            auto t_xxzz_xyzz = primBuffer.data(225 * idx + 83);

            auto t_xxzz_xzzz = primBuffer.data(225 * idx + 84);

            // Batch of Integrals (80,85)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xzz, pa2pb_x_y, pa2pb_x_yyy, \
                                     pa2pb_x_yyz, pa2pb_x_yzz, pa2pb_x_z, pa2pb_x_zzz, pa2pb_xx_xx, pa2pb_xx_xxzz, \
                                     pa2pb_xx_xy, pa2pb_xx_xyyy, pa2pb_xx_xyyz, pa2pb_xx_xyzz, pa2pb_xx_xz, \
                                     pa2pb_xx_xzzz, pa2pb_xx_zz, pa2pb_xxz_x, pa2pb_xxz_xxz, pa2pb_xxz_xyy, \
                                     pa2pb_xxz_xyz, pa2pb_xxz_xzz, pa2pb_xxz_z, pa2pb_xxzz_xx, pa2pb_xxzz_xxzz, \
                                     pa2pb_xxzz_xy, pa2pb_xxzz_xyyy, pa2pb_xxzz_xyyz, pa2pb_xxzz_xyzz, pa2pb_xxzz_xz, \
                                     pa2pb_xxzz_xzzz, pa2pb_xxzz_zz, pa2pb_xz_xz, pa2pb_xz_yy, pa2pb_xz_yz, pa2pb_xz_zz, \
                                     pa2pb_xzz_x, pa2pb_xzz_xzz, pa2pb_xzz_y, pa2pb_xzz_yyy, pa2pb_xzz_yyz, \
                                     pa2pb_xzz_yzz, pa2pb_xzz_z, pa2pb_xzz_zzz, pa2pb_z_x, pa2pb_z_xxz, pa2pb_z_xyy, \
                                     pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_z_z, pa2pb_zz_xx, pa2pb_zz_xxzz, pa2pb_zz_xy, \
                                     pa2pb_zz_xyyy, pa2pb_zz_xyyz, pa2pb_zz_xyzz, pa2pb_zz_xz, pa2pb_zz_xzzz, \
                                     pa2pb_zz_zz, pa_xx, pa_xxzz, pa_xz, pa_zz, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, \
                                     pb_xz, pb_xzzz, pb_zz, r_0_0, s_0_0, t_xxzz_xxzz, t_xxzz_xyyy, t_xxzz_xyyz, \
                                     t_xxzz_xyzz, t_xxzz_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxzz_xxzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 1.5 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa_zz[j] * fl3_fx + 1.5 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_xx[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 0.25 * pa_xxzz[j] * fl2_fx + pa2pb_xxz_z[j] * fl2_fx + 0.75 * pa2pb_xx_xx[j] * fl2_fx + pa2pb_xzz_x[j] * fl2_fx + 4.0 * pa2pb_xz_xz[j] * fl2_fx + 0.75 * pa2pb_zz_zz[j] * fl2_fx + 0.25 * pa2pb_xx_zz[j] * fl2_fx + pa2pb_x_xzz[j] * fl2_fx + 0.25 * pa2pb_zz_xx[j] * fl2_fx + pa2pb_z_xxz[j] * fl2_fx + 0.5 * pa2pb_xxzz_xx[j] * fl1_fx + 0.5 * pa2pb_xxzz_zz[j] * fl1_fx + 2.0 * pa2pb_xxz_xxz[j] * fl1_fx + 2.0 * pa2pb_xzz_xzz[j] * fl1_fx + 0.25 * pb_xxzz[j] * fl2_fx + 0.5 * pa2pb_xx_xxzz[j] * fl1_fx + 0.5 * pa2pb_zz_xxzz[j] * fl1_fx + pa2pb_xxzz_xxzz[j]);

                t_xxzz_xxzz[j] += fl_r_0_0 * (4.5 * fl4_fx * fl1_fz - 0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_xx[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 3.75 * pa_zz[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 0.25 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - pb_xx[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz + 3.0 * pa_xxzz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xxz_z[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xx_xx[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xzz_x[j] * fl2_fx * fl1_fz + 48.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz - 0.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xxzz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_zz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxzz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxz_xxz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xzz_xzz[j] * fl1_fx * fl1_fz - pb_xxzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxzz[j] * fl2_fx * fl1_fz - pa2pb_zz_xxzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xxzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xxzz[j] * fl1_fz);

                t_xxzz_xyyy[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl3_fx + 1.5 * pa2pb_xzz_y[j] * fl2_fx + 0.375 * pb_xy[j] * fl3_fx + 0.75 * pa2pb_xx_xy[j] * fl2_fx + 0.5 * pa2pb_x_yyy[j] * fl2_fx + 0.75 * pa2pb_zz_xy[j] * fl2_fx + 1.5 * pa2pb_xxzz_xy[j] * fl1_fx + pa2pb_xzz_yyy[j] * fl1_fx + 0.25 * pb_xyyy[j] * fl2_fx + 0.5 * pa2pb_xx_xyyy[j] * fl1_fx + 0.5 * pa2pb_zz_xyyy[j] * fl1_fx + pa2pb_xxzz_xyyy[j]);

                t_xxzz_xyyy[j] += fl_r_0_0 * (-1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xzz_y[j] * fl2_fx * fl1_fz - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xy[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xx_xy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxzz_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_yyy[j] * fl1_fx * fl1_fz - pb_xyyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyyy[j] * fl1_fz * fl1_fga + 3.0 * pb_xyyy[j] * fl2_fx * fl1_fz - pa2pb_zz_xyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xyyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xyyy[j] * fl1_fz);

                t_xxzz_xyyz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl3_fx + 0.25 * pa2pb_x_z[j] * fl3_fx + 0.25 * pa2pb_z_x[j] * fl3_fx + 0.5 * pa2pb_xxz_x[j] * fl2_fx + 0.5 * pa2pb_xzz_z[j] * fl2_fx + pa2pb_xz_yy[j] * fl2_fx + 0.125 * pb_xz[j] * fl3_fx + 0.25 * pa2pb_xx_xz[j] * fl2_fx + 0.5 * pa2pb_x_yyz[j] * fl2_fx + 0.25 * pa2pb_zz_xz[j] * fl2_fx + 0.5 * pa2pb_z_xyy[j] * fl2_fx + 0.5 * pa2pb_xxzz_xz[j] * fl1_fx + pa2pb_xxz_xyy[j] * fl1_fx + pa2pb_xzz_yyz[j] * fl1_fx + 0.25 * pb_xyyz[j] * fl2_fx + 0.5 * pa2pb_xx_xyyz[j] * fl1_fx + 0.5 * pa2pb_zz_xyyz[j] * fl1_fx + pa2pb_xxzz_xyyz[j]);

                t_xxzz_xyyz[j] += fl_r_0_0 * (-pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 5.0 * pa_xz[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_xxz_x[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xzz_z[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz - 0.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_xz[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xyy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_xz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_xz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxzz_xz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxz_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_yyz[j] * fl1_fx * fl1_fz - pb_xyyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyyz[j] * fl1_fz * fl1_fga + 3.0 * pb_xyyz[j] * fl2_fx * fl1_fz - pa2pb_zz_xyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xyyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xyyz[j] * fl1_fz);

                t_xxzz_xyzz[j] = fl_s_0_0 * (0.75 * pa2pb_x_y[j] * fl3_fx + 0.375 * pb_xy[j] * fl3_fx + 0.75 * pa2pb_xx_xy[j] * fl2_fx + 0.5 * pa2pb_xzz_y[j] * fl2_fx + 2.0 * pa2pb_xz_yz[j] * fl2_fx + 0.5 * pa2pb_x_yzz[j] * fl2_fx + 0.25 * pa2pb_zz_xy[j] * fl2_fx + pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_xxzz_xy[j] * fl1_fx + 2.0 * pa2pb_xxz_xyz[j] * fl1_fx + pa2pb_xzz_yzz[j] * fl1_fx + 0.25 * pb_xyzz[j] * fl2_fx + 0.5 * pa2pb_xx_xyzz[j] * fl1_fx + 0.5 * pa2pb_zz_xyzz[j] * fl1_fx + pa2pb_xxzz_xyzz[j]);

                t_xxzz_xyzz[j] += fl_r_0_0 * (7.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - pb_xy[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xx_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzz_y[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz - 0.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_xyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_xy[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxzz_xy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxz_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_yzz[j] * fl1_fx * fl1_fz - pb_xyzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xyzz[j] * fl2_fx * fl1_fz - pa2pb_zz_xyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xyzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xyzz[j] * fl1_fz);

                t_xxzz_xzzz[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl3_fx + 2.25 * pa2pb_x_z[j] * fl3_fx + 0.75 * pa2pb_z_x[j] * fl3_fx + 1.125 * pb_xz[j] * fl3_fx + 1.5 * pa2pb_xxz_x[j] * fl2_fx + 2.25 * pa2pb_xx_xz[j] * fl2_fx + 1.5 * pa2pb_xzz_z[j] * fl2_fx + 3.0 * pa2pb_xz_zz[j] * fl2_fx + 0.5 * pa2pb_x_zzz[j] * fl2_fx + 0.75 * pa2pb_zz_xz[j] * fl2_fx + 1.5 * pa2pb_z_xzz[j] * fl2_fx + 1.5 * pa2pb_xxzz_xz[j] * fl1_fx + 3.0 * pa2pb_xxz_xzz[j] * fl1_fx + pa2pb_xzz_zzz[j] * fl1_fx + 0.25 * pb_xzzz[j] * fl2_fx + 0.5 * pa2pb_xx_xzzz[j] * fl1_fx + 0.5 * pa2pb_zz_xzzz[j] * fl1_fx + pa2pb_xxzz_xzzz[j]);

                t_xxzz_xzzz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_xz[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_xz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_xxz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 11.25 * pb_xz[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xxz_x[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xx_xz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xzz_z[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_xz[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxzz_xz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxz_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_zzz[j] * fl1_fx * fl1_fz - pb_xzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xzzz[j] * fl2_fx * fl1_fz - pa2pb_zz_xzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_xzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_85_90(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pbDistances,
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

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_xx_yy = pa2pbDistances.data(1156 * idx + 108);

            auto pa2pb_xx_yz = pa2pbDistances.data(1156 * idx + 109);

            auto pa2pb_xx_zz = pa2pbDistances.data(1156 * idx + 110);

            auto pa2pb_xx_yyyy = pa2pbDistances.data(1156 * idx + 131);

            auto pa2pb_xx_yyyz = pa2pbDistances.data(1156 * idx + 132);

            auto pa2pb_xx_yyzz = pa2pbDistances.data(1156 * idx + 133);

            auto pa2pb_xx_yzzz = pa2pbDistances.data(1156 * idx + 134);

            auto pa2pb_xx_zzzz = pa2pbDistances.data(1156 * idx + 135);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_zz_yyyy = pa2pbDistances.data(1156 * idx + 301);

            auto pa2pb_zz_yyyz = pa2pbDistances.data(1156 * idx + 302);

            auto pa2pb_zz_yyzz = pa2pbDistances.data(1156 * idx + 303);

            auto pa2pb_zz_yzzz = pa2pbDistances.data(1156 * idx + 304);

            auto pa2pb_zz_zzzz = pa2pbDistances.data(1156 * idx + 305);

            auto pa2pb_xxz_y = pa2pbDistances.data(1156 * idx + 375);

            auto pa2pb_xxz_z = pa2pbDistances.data(1156 * idx + 376);

            auto pa2pb_xxz_yyy = pa2pbDistances.data(1156 * idx + 389);

            auto pa2pb_xxz_yyz = pa2pbDistances.data(1156 * idx + 390);

            auto pa2pb_xxz_yzz = pa2pbDistances.data(1156 * idx + 391);

            auto pa2pb_xxz_zzz = pa2pbDistances.data(1156 * idx + 392);

            auto pa2pb_xxzz_yy = pa2pbDistances.data(1156 * idx + 822);

            auto pa2pb_xxzz_yz = pa2pbDistances.data(1156 * idx + 823);

            auto pa2pb_xxzz_zz = pa2pbDistances.data(1156 * idx + 824);

            auto pa2pb_xxzz_yyyy = pa2pbDistances.data(1156 * idx + 845);

            auto pa2pb_xxzz_yyyz = pa2pbDistances.data(1156 * idx + 846);

            auto pa2pb_xxzz_yyzz = pa2pbDistances.data(1156 * idx + 847);

            auto pa2pb_xxzz_yzzz = pa2pbDistances.data(1156 * idx + 848);

            auto pa2pb_xxzz_zzzz = pa2pbDistances.data(1156 * idx + 849);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxzz_yyyy = primBuffer.data(225 * idx + 85);

            auto t_xxzz_yyyz = primBuffer.data(225 * idx + 86);

            auto t_xxzz_yyzz = primBuffer.data(225 * idx + 87);

            auto t_xxzz_yzzz = primBuffer.data(225 * idx + 88);

            auto t_xxzz_zzzz = primBuffer.data(225 * idx + 89);

            // Batch of Integrals (85,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_xx_yy, pa2pb_xx_yyyy, pa2pb_xx_yyyz, pa2pb_xx_yyzz, \
                                     pa2pb_xx_yz, pa2pb_xx_yzzz, pa2pb_xx_zz, pa2pb_xx_zzzz, pa2pb_xxz_y, \
                                     pa2pb_xxz_yyy, pa2pb_xxz_yyz, pa2pb_xxz_yzz, pa2pb_xxz_z, pa2pb_xxz_zzz, \
                                     pa2pb_xxzz_yy, pa2pb_xxzz_yyyy, pa2pb_xxzz_yyyz, pa2pb_xxzz_yyzz, pa2pb_xxzz_yz, \
                                     pa2pb_xxzz_yzzz, pa2pb_xxzz_zz, pa2pb_xxzz_zzzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, \
                                     pa2pb_z_yzz, pa2pb_z_z, pa2pb_z_zzz, pa2pb_zz_yy, pa2pb_zz_yyyy, pa2pb_zz_yyyz, \
                                     pa2pb_zz_yyzz, pa2pb_zz_yz, pa2pb_zz_yzzz, pa2pb_zz_zz, pa2pb_zz_zzzz, pa_xx, pa_xxzz, \
                                     pa_zz, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, r_0_0, s_0_0, \
                                     t_xxzz_yyyy, t_xxzz_yyyz, t_xxzz_yyzz, t_xxzz_yzzz, t_xxzz_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxzz_yyyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * pa_zz[j] * fl3_fx + 0.75 * pa_xxzz[j] * fl2_fx + 0.75 * pb_yy[j] * fl3_fx + 1.5 * pa2pb_xx_yy[j] * fl2_fx + 1.5 * pa2pb_zz_yy[j] * fl2_fx + 3.0 * pa2pb_xxzz_yy[j] * fl1_fx + 0.25 * pb_yyyy[j] * fl2_fx + 0.5 * pa2pb_xx_yyyy[j] * fl1_fx + 0.5 * pa2pb_zz_yyyy[j] * fl1_fx + pa2pb_xxzz_yyyy[j]);

                t_xxzz_yyyy[j] += fl_r_0_0 * (-0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - 1.5 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz - 0.75 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xx[j] * fl1_fz * fl3_fx + 3.75 * pa_zz[j] * fl3_fx * fl1_fz + 9.0 * pa_xxzz[j] * fl1_fz * fl2_fx - 1.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_yy[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxzz_yy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xx_yy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxzz_yy[j] * fl1_fz * fl1_fx - pb_yyyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyyy[j] * fl1_fz * fl1_fga + 3.0 * pb_yyyy[j] * fl2_fx * fl1_fz - pa2pb_zz_yyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_yyyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_yyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_yyyy[j] * fl1_fz);

                t_xxzz_yyyz[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl3_fx + 1.5 * pa2pb_xxz_y[j] * fl2_fx + 0.375 * pb_yz[j] * fl3_fx + 0.75 * pa2pb_xx_yz[j] * fl2_fx + 0.75 * pa2pb_zz_yz[j] * fl2_fx + 0.5 * pa2pb_z_yyy[j] * fl2_fx + 1.5 * pa2pb_xxzz_yz[j] * fl1_fx + pa2pb_xxz_yyy[j] * fl1_fx + 0.25 * pb_yyyz[j] * fl2_fx + 0.5 * pa2pb_xx_yyyz[j] * fl1_fx + 0.5 * pa2pb_zz_yyyz[j] * fl1_fx + pa2pb_xxzz_yyyz[j]);

                t_xxzz_yyyz[j] += fl_r_0_0 * (-1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xxz_y[j] * fl1_fz * fl2_fx - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_yz[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xx_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxzz_yz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxz_yyy[j] * fl1_fz * fl1_fx - pb_yyyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyyz[j] * fl1_fz * fl1_fga + 3.0 * pb_yyyz[j] * fl2_fx * fl1_fz - pa2pb_zz_yyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_yyyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_yyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_yyyz[j] * fl1_fz);

                t_xxzz_yyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.125 * pa_zz[j] * fl3_fx + 0.5 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 0.25 * pa_xxzz[j] * fl2_fx + pa2pb_xxz_z[j] * fl2_fx + 0.75 * pa2pb_xx_yy[j] * fl2_fx + 0.125 * pb_zz[j] * fl3_fx + 0.25 * pa2pb_xx_zz[j] * fl2_fx + 0.25 * pa2pb_zz_yy[j] * fl2_fx + 0.25 * pa2pb_zz_zz[j] * fl2_fx + pa2pb_z_yyz[j] * fl2_fx + 0.5 * pa2pb_xxzz_yy[j] * fl1_fx + 0.5 * pa2pb_xxzz_zz[j] * fl1_fx + 2.0 * pa2pb_xxz_yyz[j] * fl1_fx + 0.25 * pb_yyzz[j] * fl2_fx + 0.5 * pa2pb_xx_yyzz[j] * fl1_fx + 0.5 * pa2pb_zz_yyzz[j] * fl1_fx + pa2pb_xxzz_yyzz[j]);

                t_xxzz_yyzz[j] += fl_r_0_0 * (-0.5 * fl3_fx * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * fl3_fx - pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_xx[j] * fl3_fx * fl1_fz - 0.25 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - pb_yy[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_zz[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 3.0 * pa_xxzz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xxz_z[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xx_yy[j] * fl2_fx * fl1_fz - 0.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_zz[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_yyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xxzz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xx_zz[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xxzz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xxzz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xxz_yyz[j] * fl1_fz * fl1_fx - pb_yyzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyzz[j] * fl1_fz * fl1_fga + 3.0 * pb_yyzz[j] * fl2_fx * fl1_fz - pa2pb_zz_yyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_yyzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_yyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_yyzz[j] * fl1_fz);

                t_xxzz_yzzz[j] = fl_s_0_0 * (0.75 * pa2pb_z_y[j] * fl3_fx + 1.125 * pb_yz[j] * fl3_fx + 1.5 * pa2pb_xxz_y[j] * fl2_fx + 2.25 * pa2pb_xx_yz[j] * fl2_fx + 0.75 * pa2pb_zz_yz[j] * fl2_fx + 1.5 * pa2pb_z_yzz[j] * fl2_fx + 1.5 * pa2pb_xxzz_yz[j] * fl1_fx + 3.0 * pa2pb_xxz_yzz[j] * fl1_fx + 0.25 * pb_yzzz[j] * fl2_fx + 0.5 * pa2pb_xx_yzzz[j] * fl1_fx + 0.5 * pa2pb_zz_yzzz[j] * fl1_fx + pa2pb_xxzz_yzzz[j]);

                t_xxzz_yzzz[j] += fl_r_0_0 * (-1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_yz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_xxz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 11.25 * pb_yz[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xxz_y[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xx_yz[j] * fl2_fx * fl1_fz - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_yzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xxzz_yz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xxz_yzz[j] * fl1_fz * fl1_fx - pb_yzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_yzzz[j] * fl2_fx * fl1_fz - pa2pb_zz_yzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_yzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_yzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_yzzz[j] * fl1_fz);

                t_xxzz_zzzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_xx[j] * fl3_fx + 0.375 * pa_zz[j] * fl3_fx + 3.0 * pa2pb_z_z[j] * fl3_fx + 2.25 * pb_zz[j] * fl3_fx + 0.75 * pa_xxzz[j] * fl2_fx + 6.0 * pa2pb_xxz_z[j] * fl2_fx + 4.5 * pa2pb_xx_zz[j] * fl2_fx + 1.5 * pa2pb_zz_zz[j] * fl2_fx + 2.0 * pa2pb_z_zzz[j] * fl2_fx + 3.0 * pa2pb_xxzz_zz[j] * fl1_fx + 4.0 * pa2pb_xxz_zzz[j] * fl1_fx + 0.25 * pb_zzzz[j] * fl2_fx + 0.5 * pa2pb_xx_zzzz[j] * fl1_fx + 0.5 * pa2pb_zz_zzzz[j] * fl1_fx + pa2pb_xxzz_zzzz[j]);

                t_xxzz_zzzz[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl1_fz * fl1_fga * fl3_fx - 4.5 * pa_xx[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 18.75 * pa_xx[j] * fl3_fx * fl1_fz - 0.75 * pa_xx[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pb_zz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxzz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xxz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_zz[j] * fl3_fx * fl1_fz + 30.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 22.5 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa_xxzz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_xxz_z[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_xx_zz[j] * fl2_fx * fl1_fz - 1.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 4.0 * pa2pb_z_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xxzz_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xxzz_zz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xxz_zzz[j] * fl1_fz * fl1_fx - pb_zzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_zzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_zzzz[j] * fl2_fx * fl1_fz - pa2pb_zz_zzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xx_zzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_zzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xxzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_90_95(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_xxxx = pa2pbDistances.data(1156 * idx + 155);

            auto pa2pb_xy_xxxy = pa2pbDistances.data(1156 * idx + 156);

            auto pa2pb_xy_xxxz = pa2pbDistances.data(1156 * idx + 157);

            auto pa2pb_xy_xxyy = pa2pbDistances.data(1156 * idx + 158);

            auto pa2pb_xy_xxyz = pa2pbDistances.data(1156 * idx + 159);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_xyy_x = pa2pbDistances.data(1156 * idx + 408);

            auto pa2pb_xyy_y = pa2pbDistances.data(1156 * idx + 409);

            auto pa2pb_xyy_z = pa2pbDistances.data(1156 * idx + 410);

            auto pa2pb_xyy_xxx = pa2pbDistances.data(1156 * idx + 417);

            auto pa2pb_xyy_xxy = pa2pbDistances.data(1156 * idx + 418);

            auto pa2pb_xyy_xxz = pa2pbDistances.data(1156 * idx + 419);

            auto pa2pb_yyy_x = pa2pbDistances.data(1156 * idx + 510);

            auto pa2pb_yyy_y = pa2pbDistances.data(1156 * idx + 511);

            auto pa2pb_yyy_z = pa2pbDistances.data(1156 * idx + 512);

            auto pa2pb_yyy_xxx = pa2pbDistances.data(1156 * idx + 519);

            auto pa2pb_yyy_xxy = pa2pbDistances.data(1156 * idx + 520);

            auto pa2pb_yyy_xxz = pa2pbDistances.data(1156 * idx + 521);

            auto pa2pb_yyy_xyy = pa2pbDistances.data(1156 * idx + 522);

            auto pa2pb_yyy_xyz = pa2pbDistances.data(1156 * idx + 523);

            auto pa2pb_xyyy_xx = pa2pbDistances.data(1156 * idx + 853);

            auto pa2pb_xyyy_xy = pa2pbDistances.data(1156 * idx + 854);

            auto pa2pb_xyyy_xz = pa2pbDistances.data(1156 * idx + 855);

            auto pa2pb_xyyy_yy = pa2pbDistances.data(1156 * idx + 856);

            auto pa2pb_xyyy_yz = pa2pbDistances.data(1156 * idx + 857);

            auto pa2pb_xyyy_xxxx = pa2pbDistances.data(1156 * idx + 869);

            auto pa2pb_xyyy_xxxy = pa2pbDistances.data(1156 * idx + 870);

            auto pa2pb_xyyy_xxxz = pa2pbDistances.data(1156 * idx + 871);

            auto pa2pb_xyyy_xxyy = pa2pbDistances.data(1156 * idx + 872);

            auto pa2pb_xyyy_xxyz = pa2pbDistances.data(1156 * idx + 873);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyy_xxxx = primBuffer.data(225 * idx + 90);

            auto t_xyyy_xxxy = primBuffer.data(225 * idx + 91);

            auto t_xyyy_xxxz = primBuffer.data(225 * idx + 92);

            auto t_xyyy_xxyy = primBuffer.data(225 * idx + 93);

            auto t_xyyy_xxyz = primBuffer.data(225 * idx + 94);

            // Batch of Integrals (90,95)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_y, pa2pb_x_z, pa2pb_xy_xx, pa2pb_xy_xxxx, pa2pb_xy_xxxy, pa2pb_xy_xxxz, \
                                     pa2pb_xy_xxyy, pa2pb_xy_xxyz, pa2pb_xy_xy, pa2pb_xy_xz, pa2pb_xy_yy, pa2pb_xy_yz, \
                                     pa2pb_xyy_x, pa2pb_xyy_xxx, pa2pb_xyy_xxy, pa2pb_xyy_xxz, pa2pb_xyy_y, \
                                     pa2pb_xyy_z, pa2pb_xyyy_xx, pa2pb_xyyy_xxxx, pa2pb_xyyy_xxxy, pa2pb_xyyy_xxxz, \
                                     pa2pb_xyyy_xxyy, pa2pb_xyyy_xxyz, pa2pb_xyyy_xy, pa2pb_xyyy_xz, pa2pb_xyyy_yy, \
                                     pa2pb_xyyy_yz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, pa2pb_y_xyy, \
                                     pa2pb_y_xyz, pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xy, pa2pb_yy_xz, \
                                     pa2pb_yyy_x, pa2pb_yyy_xxx, pa2pb_yyy_xxy, pa2pb_yyy_xxz, pa2pb_yyy_xyy, \
                                     pa2pb_yyy_xyz, pa2pb_yyy_y, pa2pb_yyy_z, pa_xy, pa_xyyy, pa_yy, pb_xx, pb_xy, pb_xz, r_0_0, \
                                     s_0_0, t_xyyy_xxxx, t_xyyy_xxxy, t_xyyy_xxxz, t_xyyy_xxyy, t_xyyy_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyyy_xxxx[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 4.5 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa_xyyy[j] * fl2_fx + 3.0 * pa2pb_yyy_x[j] * fl2_fx + 4.5 * pa2pb_xy_xx[j] * fl2_fx + 3.0 * pa2pb_y_xxx[j] * fl2_fx + 3.0 * pa2pb_xyyy_xx[j] * fl1_fx + 2.0 * pa2pb_yyy_xxx[j] * fl1_fx + 1.5 * pa2pb_xy_xxxx[j] * fl1_fx + pa2pb_xyyy_xxxx[j]);

                t_xyyy_xxxx[j] += fl_r_0_0 * (-4.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xyyy[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_xy[j] * fl1_fz * fl3_fx - 6.0 * pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 9.0 * pa_xyyy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_yyy_x[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xyyy_xx[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xy_xx[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xyyy_xx[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyy_xxx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxxx[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxxx[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xxxx[j] * fl1_fz);

                t_xyyy_xxxy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_yy[j] * fl3_fx + 1.125 * pa2pb_x_x[j] * fl3_fx + 1.125 * pa2pb_y_y[j] * fl3_fx + 1.125 * pb_xx[j] * fl3_fx + 2.25 * pa2pb_xyy_x[j] * fl2_fx + 0.75 * pa2pb_yyy_y[j] * fl2_fx + 2.25 * pa2pb_yy_xx[j] * fl2_fx + 2.25 * pa2pb_xy_xy[j] * fl2_fx + 0.75 * pa2pb_x_xxx[j] * fl2_fx + 2.25 * pa2pb_y_xxy[j] * fl2_fx + 1.5 * pa2pb_xyyy_xy[j] * fl1_fx + 1.5 * pa2pb_xyy_xxx[j] * fl1_fx + 1.5 * pa2pb_yyy_xxy[j] * fl1_fx + 1.5 * pa2pb_xy_xxxy[j] * fl1_fx + pa2pb_xyyy_xxxy[j]);

                t_xyyy_xxxy[j] += fl_r_0_0 * (-1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 4.5 * fl4_fx * fl1_fz + 11.25 * pa_yy[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 11.25 * pb_xx[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xyy_x[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yyy_y[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yy_xx[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_y_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyy_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xy_xy[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyy_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xyy_xxx[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyy_xxy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxxy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxxy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xxxy[j] * fl1_fz);

                t_xyyy_xxxz[j] = fl_s_0_0 * (1.125 * pa2pb_y_z[j] * fl3_fx + 0.75 * pa2pb_yyy_z[j] * fl2_fx + 2.25 * pa2pb_xy_xz[j] * fl2_fx + 2.25 * pa2pb_y_xxz[j] * fl2_fx + 1.5 * pa2pb_xyyy_xz[j] * fl1_fx + 1.5 * pa2pb_yyy_xxz[j] * fl1_fx + 1.5 * pa2pb_xy_xxxz[j] * fl1_fx + pa2pb_xyyy_xxxz[j]);

                t_xyyy_xxxz[j] += fl_r_0_0 * (-2.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yyy_z[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyy_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xy_xz[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyy_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyy_xxz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxxz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxxz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xxxz[j] * fl1_fz);

                t_xyyy_xxyy[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 2.25 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa2pb_x_y[j] * fl3_fx + 1.5 * pb_xy[j] * fl3_fx + 0.25 * pa_xyyy[j] * fl2_fx + 1.5 * pa2pb_xyy_y[j] * fl2_fx + 2.25 * pa2pb_xy_xx[j] * fl2_fx + 0.5 * pa2pb_yyy_x[j] * fl2_fx + 3.0 * pa2pb_yy_xy[j] * fl2_fx + 0.75 * pa2pb_xy_yy[j] * fl2_fx + 1.5 * pa2pb_x_xxy[j] * fl2_fx + 1.5 * pa2pb_y_xyy[j] * fl2_fx + 0.5 * pa2pb_xyyy_xx[j] * fl1_fx + 0.5 * pa2pb_xyyy_yy[j] * fl1_fx + 3.0 * pa2pb_xyy_xxy[j] * fl1_fx + pa2pb_yyy_xyy[j] * fl1_fx + 1.5 * pa2pb_xy_xxyy[j] * fl1_fx + pa2pb_xyyy_xxyy[j]);

                t_xyyy_xxyy[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xy[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyyy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz - pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pb_xy[j] * fl3_fx * fl1_fz + 3.0 * pa_xyyy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyy_y[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xy_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyy_x[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_xy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_y_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyyy_yy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_yy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyy_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyyy_yy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyy_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xxyy[j] * fl1_fz);

                t_xyyy_xxyz[j] = fl_s_0_0 * (0.375 * pa2pb_x_z[j] * fl3_fx + 0.75 * pb_xz[j] * fl3_fx + 0.75 * pa2pb_xyy_z[j] * fl2_fx + 1.5 * pa2pb_yy_xz[j] * fl2_fx + 0.75 * pa2pb_xy_yz[j] * fl2_fx + 0.75 * pa2pb_x_xxz[j] * fl2_fx + 1.5 * pa2pb_y_xyz[j] * fl2_fx + 0.5 * pa2pb_xyyy_yz[j] * fl1_fx + 1.5 * pa2pb_xyy_xxz[j] * fl1_fx + pa2pb_yyy_xyz[j] * fl1_fx + 1.5 * pa2pb_xy_xxyz[j] * fl1_fx + pa2pb_xyyy_xxyz[j]);

                t_xyyy_xxyz[j] += fl_r_0_0 * (-0.75 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 7.5 * pb_xz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xyy_z[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yy_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_y_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyy_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xyy_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_95_100(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_xy_xxzz = pa2pbDistances.data(1156 * idx + 160);

            auto pa2pb_xy_xyyy = pa2pbDistances.data(1156 * idx + 161);

            auto pa2pb_xy_xyyz = pa2pbDistances.data(1156 * idx + 162);

            auto pa2pb_xy_xyzz = pa2pbDistances.data(1156 * idx + 163);

            auto pa2pb_xy_xzzz = pa2pbDistances.data(1156 * idx + 164);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_xyy_x = pa2pbDistances.data(1156 * idx + 408);

            auto pa2pb_xyy_xyy = pa2pbDistances.data(1156 * idx + 420);

            auto pa2pb_xyy_xyz = pa2pbDistances.data(1156 * idx + 421);

            auto pa2pb_xyy_xzz = pa2pbDistances.data(1156 * idx + 422);

            auto pa2pb_yyy_x = pa2pbDistances.data(1156 * idx + 510);

            auto pa2pb_yyy_y = pa2pbDistances.data(1156 * idx + 511);

            auto pa2pb_yyy_z = pa2pbDistances.data(1156 * idx + 512);

            auto pa2pb_yyy_xzz = pa2pbDistances.data(1156 * idx + 524);

            auto pa2pb_yyy_yyy = pa2pbDistances.data(1156 * idx + 525);

            auto pa2pb_yyy_yyz = pa2pbDistances.data(1156 * idx + 526);

            auto pa2pb_yyy_yzz = pa2pbDistances.data(1156 * idx + 527);

            auto pa2pb_yyy_zzz = pa2pbDistances.data(1156 * idx + 528);

            auto pa2pb_xyyy_xx = pa2pbDistances.data(1156 * idx + 853);

            auto pa2pb_xyyy_xy = pa2pbDistances.data(1156 * idx + 854);

            auto pa2pb_xyyy_xz = pa2pbDistances.data(1156 * idx + 855);

            auto pa2pb_xyyy_zz = pa2pbDistances.data(1156 * idx + 858);

            auto pa2pb_xyyy_xxzz = pa2pbDistances.data(1156 * idx + 874);

            auto pa2pb_xyyy_xyyy = pa2pbDistances.data(1156 * idx + 875);

            auto pa2pb_xyyy_xyyz = pa2pbDistances.data(1156 * idx + 876);

            auto pa2pb_xyyy_xyzz = pa2pbDistances.data(1156 * idx + 877);

            auto pa2pb_xyyy_xzzz = pa2pbDistances.data(1156 * idx + 878);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyy_xxzz = primBuffer.data(225 * idx + 95);

            auto t_xyyy_xyyy = primBuffer.data(225 * idx + 96);

            auto t_xyyy_xyyz = primBuffer.data(225 * idx + 97);

            auto t_xyyy_xyzz = primBuffer.data(225 * idx + 98);

            auto t_xyyy_xzzz = primBuffer.data(225 * idx + 99);

            // Batch of Integrals (95,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_xzz, \
                                     pa2pb_xy_xx, pa2pb_xy_xxzz, pa2pb_xy_xy, pa2pb_xy_xyyy, pa2pb_xy_xyyz, \
                                     pa2pb_xy_xyzz, pa2pb_xy_xz, pa2pb_xy_xzzz, pa2pb_xy_zz, pa2pb_xyy_x, pa2pb_xyy_xyy, \
                                     pa2pb_xyy_xyz, pa2pb_xyy_xzz, pa2pb_xyyy_xx, pa2pb_xyyy_xxzz, pa2pb_xyyy_xy, \
                                     pa2pb_xyyy_xyyy, pa2pb_xyyy_xyyz, pa2pb_xyyy_xyzz, pa2pb_xyyy_xz, pa2pb_xyyy_xzzz, \
                                     pa2pb_xyyy_zz, pa2pb_y_x, pa2pb_y_xzz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, \
                                     pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, pa2pb_yy_yy, pa2pb_yy_yz, pa2pb_yy_zz, \
                                     pa2pb_yyy_x, pa2pb_yyy_xzz, pa2pb_yyy_y, pa2pb_yyy_yyy, pa2pb_yyy_yyz, \
                                     pa2pb_yyy_yzz, pa2pb_yyy_z, pa2pb_yyy_zzz, pa_xy, pa_xyyy, pa_yy, pb_yy, pb_yz, pb_zz, \
                                     r_0_0, s_0_0, t_xyyy_xxzz, t_xyyy_xyyy, t_xyyy_xyyz, t_xyyy_xyzz, t_xyyy_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyyy_xxzz[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 0.75 * pa2pb_y_x[j] * fl3_fx + 0.25 * pa_xyyy[j] * fl2_fx + 0.5 * pa2pb_yyy_x[j] * fl2_fx + 0.75 * pa2pb_xy_xx[j] * fl2_fx + 0.75 * pa2pb_xy_zz[j] * fl2_fx + 1.5 * pa2pb_y_xzz[j] * fl2_fx + 0.5 * pa2pb_xyyy_xx[j] * fl1_fx + 0.5 * pa2pb_xyyy_zz[j] * fl1_fx + pa2pb_yyy_xzz[j] * fl1_fx + 1.5 * pa2pb_xy_xxzz[j] * fl1_fx + pa2pb_xyyy_xxzz[j]);

                t_xyyy_xxzz[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xy[j] * fl1_fz * fl3_fx - pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 3.0 * pa_xyyy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yyy_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyyy_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xx[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xy_zz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyy_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyyy_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xxzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xxzz[j] * fl1_fz);

                t_xyyy_xyyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa2pb_x_x[j] * fl3_fx + 1.125 * pa_yy[j] * fl3_fx + 3.375 * pa2pb_y_y[j] * fl3_fx + 1.125 * pb_yy[j] * fl3_fx + 2.25 * pa2pb_xyy_x[j] * fl2_fx + 6.75 * pa2pb_xy_xy[j] * fl2_fx + 0.75 * pa2pb_yyy_y[j] * fl2_fx + 2.25 * pa2pb_yy_yy[j] * fl2_fx + 2.25 * pa2pb_x_xyy[j] * fl2_fx + 0.75 * pa2pb_y_yyy[j] * fl2_fx + 1.5 * pa2pb_xyyy_xy[j] * fl1_fx + 4.5 * pa2pb_xyy_xyy[j] * fl1_fx + 0.5 * pa2pb_yyy_yyy[j] * fl1_fx + 1.5 * pa2pb_xy_xyyy[j] * fl1_fx + pa2pb_xyyy_xyyy[j]);

                t_xyyy_xyyy[j] += fl_r_0_0 * (7.5 * fl4_fx * fl1_fz - 1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 18.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 11.25 * pa_yy[j] * fl3_fx * fl1_fz + 33.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pb_yy[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xyy_x[j] * fl1_fz * fl2_fx + 81.0 * pa2pb_xy_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yyy_y[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yy_yy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyy_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyy_xy[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_xyy_xyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyy_yyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xyyy[j] * fl1_fz);

                t_xyyy_xyyz[j] = fl_s_0_0 * (1.125 * pa2pb_y_z[j] * fl3_fx + 0.75 * pb_yz[j] * fl3_fx + 2.25 * pa2pb_xy_xz[j] * fl2_fx + 0.25 * pa2pb_yyy_z[j] * fl2_fx + 1.5 * pa2pb_yy_yz[j] * fl2_fx + 1.5 * pa2pb_x_xyz[j] * fl2_fx + 0.75 * pa2pb_y_yyz[j] * fl2_fx + 0.5 * pa2pb_xyyy_xz[j] * fl1_fx + 3.0 * pa2pb_xyy_xyz[j] * fl1_fx + 0.5 * pa2pb_yyy_yyz[j] * fl1_fx + 1.5 * pa2pb_xy_xyyz[j] * fl1_fx + pa2pb_xyyy_xyyz[j]);

                t_xyyy_xyyz[j] += fl_r_0_0 * (11.25 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_yz[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xy_xz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yyy_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yy_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyy_xz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyy_xyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyy_yyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xyyz[j] * fl1_fz);

                t_xyyy_xyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 0.75 * pa2pb_xyy_x[j] * fl2_fx + 0.25 * pa2pb_yyy_y[j] * fl2_fx + 0.75 * pa2pb_yy_zz[j] * fl2_fx + 0.75 * pa2pb_xy_xy[j] * fl2_fx + 0.75 * pa2pb_x_xzz[j] * fl2_fx + 0.75 * pa2pb_y_yzz[j] * fl2_fx + 0.5 * pa2pb_xyyy_xy[j] * fl1_fx + 1.5 * pa2pb_xyy_xzz[j] * fl1_fx + 0.5 * pa2pb_yyy_yzz[j] * fl1_fx + 1.5 * pa2pb_xy_xyzz[j] * fl1_fx + pa2pb_xyyy_xyzz[j]);

                t_xyyy_xyzz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_yy[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xyy_x[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_yyy_y[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xy[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyy_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xyy_xzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyy_yzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xyzz[j] * fl1_fz);

                t_xyyy_xzzz[j] = fl_s_0_0 * (1.125 * pa2pb_y_z[j] * fl3_fx + 0.75 * pa2pb_yyy_z[j] * fl2_fx + 2.25 * pa2pb_xy_xz[j] * fl2_fx + 0.75 * pa2pb_y_zzz[j] * fl2_fx + 1.5 * pa2pb_xyyy_xz[j] * fl1_fx + 0.5 * pa2pb_yyy_zzz[j] * fl1_fx + 1.5 * pa2pb_xy_xzzz[j] * fl1_fx + pa2pb_xyyy_xzzz[j]);

                t_xyyy_xzzz[j] += fl_r_0_0 * (-2.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yyy_z[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyy_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xy_xz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyy_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyy_zzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_xzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_100_105(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_xy_yyyy = pa2pbDistances.data(1156 * idx + 165);

            auto pa2pb_xy_yyyz = pa2pbDistances.data(1156 * idx + 166);

            auto pa2pb_xy_yyzz = pa2pbDistances.data(1156 * idx + 167);

            auto pa2pb_xy_yzzz = pa2pbDistances.data(1156 * idx + 168);

            auto pa2pb_xy_zzzz = pa2pbDistances.data(1156 * idx + 169);

            auto pa2pb_xyy_y = pa2pbDistances.data(1156 * idx + 409);

            auto pa2pb_xyy_z = pa2pbDistances.data(1156 * idx + 410);

            auto pa2pb_xyy_yyy = pa2pbDistances.data(1156 * idx + 423);

            auto pa2pb_xyy_yyz = pa2pbDistances.data(1156 * idx + 424);

            auto pa2pb_xyy_yzz = pa2pbDistances.data(1156 * idx + 425);

            auto pa2pb_xyy_zzz = pa2pbDistances.data(1156 * idx + 426);

            auto pa2pb_xyyy_yy = pa2pbDistances.data(1156 * idx + 856);

            auto pa2pb_xyyy_yz = pa2pbDistances.data(1156 * idx + 857);

            auto pa2pb_xyyy_zz = pa2pbDistances.data(1156 * idx + 858);

            auto pa2pb_xyyy_yyyy = pa2pbDistances.data(1156 * idx + 879);

            auto pa2pb_xyyy_yyyz = pa2pbDistances.data(1156 * idx + 880);

            auto pa2pb_xyyy_yyzz = pa2pbDistances.data(1156 * idx + 881);

            auto pa2pb_xyyy_yzzz = pa2pbDistances.data(1156 * idx + 882);

            auto pa2pb_xyyy_zzzz = pa2pbDistances.data(1156 * idx + 883);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyy_yyyy = primBuffer.data(225 * idx + 100);

            auto t_xyyy_yyyz = primBuffer.data(225 * idx + 101);

            auto t_xyyy_yyzz = primBuffer.data(225 * idx + 102);

            auto t_xyyy_yzzz = primBuffer.data(225 * idx + 103);

            auto t_xyyy_zzzz = primBuffer.data(225 * idx + 104);

            // Batch of Integrals (100,105)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, \
                                     pa2pb_x_z, pa2pb_x_zzz, pa2pb_xy_yy, pa2pb_xy_yyyy, pa2pb_xy_yyyz, \
                                     pa2pb_xy_yyzz, pa2pb_xy_yz, pa2pb_xy_yzzz, pa2pb_xy_zz, pa2pb_xy_zzzz, pa2pb_xyy_y, \
                                     pa2pb_xyy_yyy, pa2pb_xyy_yyz, pa2pb_xyy_yzz, pa2pb_xyy_z, pa2pb_xyy_zzz, \
                                     pa2pb_xyyy_yy, pa2pb_xyyy_yyyy, pa2pb_xyyy_yyyz, pa2pb_xyyy_yyzz, pa2pb_xyyy_yz, \
                                     pa2pb_xyyy_yzzz, pa2pb_xyyy_zz, pa2pb_xyyy_zzzz, pa_xy, pa_xyyy, r_0_0, s_0_0, \
                                     t_xyyy_yyyy, t_xyyy_yyyz, t_xyyy_yyzz, t_xyyy_yzzz, t_xyyy_zzzz: VLX_ALIGN)
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

                t_xyyy_yyyy[j] = fl_s_0_0 * (5.625 * pa_xy[j] * fl3_fx + 7.5 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa_xyyy[j] * fl2_fx + 9.0 * pa2pb_xyy_y[j] * fl2_fx + 13.5 * pa2pb_xy_yy[j] * fl2_fx + 3.0 * pa2pb_x_yyy[j] * fl2_fx + 3.0 * pa2pb_xyyy_yy[j] * fl1_fx + 6.0 * pa2pb_xyy_yyy[j] * fl1_fx + 1.5 * pa2pb_xy_yyyy[j] * fl1_fx + pa2pb_xyyy_yyyy[j]);

                t_xyyy_yyyy[j] += fl_r_0_0 * (-13.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 56.25 * pa_xy[j] * fl3_fx * fl1_fz + 75.0 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 2.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xyyy[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_xyyy[j] * fl1_fz * fl2_fx + 108.0 * pa2pb_xyy_y[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_xy_yy[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xyyy_yy[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xyyy_yy[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_xyy_yyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_yyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_yyyy[j] * fl1_fz);

                t_xyyy_yyyz[j] = fl_s_0_0 * (1.875 * pa2pb_x_z[j] * fl3_fx + 2.25 * pa2pb_xyy_z[j] * fl2_fx + 6.75 * pa2pb_xy_yz[j] * fl2_fx + 2.25 * pa2pb_x_yyz[j] * fl2_fx + 1.5 * pa2pb_xyyy_yz[j] * fl1_fx + 4.5 * pa2pb_xyy_yyz[j] * fl1_fx + 1.5 * pa2pb_xy_yyyz[j] * fl1_fx + pa2pb_xyyy_yyyz[j]);

                t_xyyy_yyyz[j] += fl_r_0_0 * (18.75 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 27.0 * pa2pb_xyy_z[j] * fl1_fz * fl2_fx + 81.0 * pa2pb_xy_yz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyy_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyy_yz[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_xyy_yyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_yyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_yyyz[j] * fl1_fz);

                t_xyyy_yyzz[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 0.75 * pa2pb_x_y[j] * fl3_fx + 0.25 * pa_xyyy[j] * fl2_fx + 1.5 * pa2pb_xyy_y[j] * fl2_fx + 2.25 * pa2pb_xy_zz[j] * fl2_fx + 0.75 * pa2pb_xy_yy[j] * fl2_fx + 1.5 * pa2pb_x_yzz[j] * fl2_fx + 0.5 * pa2pb_xyyy_yy[j] * fl1_fx + 0.5 * pa2pb_xyyy_zz[j] * fl1_fx + 3.0 * pa2pb_xyy_yzz[j] * fl1_fx + 1.5 * pa2pb_xy_yyzz[j] * fl1_fx + pa2pb_xyyy_yyzz[j]);

                t_xyyy_yyzz[j] += fl_r_0_0 * (-3.0 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xy[j] * fl3_fx * fl1_fz - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyyy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 3.0 * pa_xyyy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyy_y[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xy_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_yy[j] * fl1_fz * fl1_fgb - pa2pb_xyyy_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_yy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyy_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyyy_zz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyy_yzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_yyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_yyzz[j] * fl1_fz);

                t_xyyy_yzzz[j] = fl_s_0_0 * (1.125 * pa2pb_x_z[j] * fl3_fx + 2.25 * pa2pb_xyy_z[j] * fl2_fx + 2.25 * pa2pb_xy_yz[j] * fl2_fx + 0.75 * pa2pb_x_zzz[j] * fl2_fx + 1.5 * pa2pb_xyyy_yz[j] * fl1_fx + 1.5 * pa2pb_xyy_zzz[j] * fl1_fx + 1.5 * pa2pb_xy_yzzz[j] * fl1_fx + pa2pb_xyyy_yzzz[j]);

                t_xyyy_yzzz[j] += fl_r_0_0 * (-2.25 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xyy_z[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyy_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xy_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyy_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xyy_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_yzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_yzzz[j] * fl1_fz);

                t_xyyy_zzzz[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 0.75 * pa_xyyy[j] * fl2_fx + 4.5 * pa2pb_xy_zz[j] * fl2_fx + 3.0 * pa2pb_xyyy_zz[j] * fl1_fx + 1.5 * pa2pb_xy_zzzz[j] * fl1_fx + pa2pb_xyyy_zzzz[j]);

                t_xyyy_zzzz[j] += fl_r_0_0 * (-4.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xyyy[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_xy[j] * fl1_fz * fl3_fx + 9.0 * pa_xyyy[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xyyy_zz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xy_zz[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_xyyy_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_zzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xy_zzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyyy_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_105_110(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_xxxx = pa2pbDistances.data(1156 * idx + 189);

            auto pa2pb_xz_xxxy = pa2pbDistances.data(1156 * idx + 190);

            auto pa2pb_xz_xxxz = pa2pbDistances.data(1156 * idx + 191);

            auto pa2pb_xz_xxyy = pa2pbDistances.data(1156 * idx + 192);

            auto pa2pb_xz_xxyz = pa2pbDistances.data(1156 * idx + 193);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_xyy_x = pa2pbDistances.data(1156 * idx + 408);

            auto pa2pb_xyy_y = pa2pbDistances.data(1156 * idx + 409);

            auto pa2pb_xyy_xxx = pa2pbDistances.data(1156 * idx + 417);

            auto pa2pb_xyy_xxy = pa2pbDistances.data(1156 * idx + 418);

            auto pa2pb_xyz_x = pa2pbDistances.data(1156 * idx + 442);

            auto pa2pb_xyz_y = pa2pbDistances.data(1156 * idx + 443);

            auto pa2pb_xyz_z = pa2pbDistances.data(1156 * idx + 444);

            auto pa2pb_xyz_xxx = pa2pbDistances.data(1156 * idx + 451);

            auto pa2pb_xyz_xxy = pa2pbDistances.data(1156 * idx + 452);

            auto pa2pb_xyz_xxz = pa2pbDistances.data(1156 * idx + 453);

            auto pa2pb_yyz_x = pa2pbDistances.data(1156 * idx + 544);

            auto pa2pb_yyz_y = pa2pbDistances.data(1156 * idx + 545);

            auto pa2pb_yyz_z = pa2pbDistances.data(1156 * idx + 546);

            auto pa2pb_yyz_xxx = pa2pbDistances.data(1156 * idx + 553);

            auto pa2pb_yyz_xxy = pa2pbDistances.data(1156 * idx + 554);

            auto pa2pb_yyz_xxz = pa2pbDistances.data(1156 * idx + 555);

            auto pa2pb_yyz_xyy = pa2pbDistances.data(1156 * idx + 556);

            auto pa2pb_yyz_xyz = pa2pbDistances.data(1156 * idx + 557);

            auto pa2pb_xyyz_xx = pa2pbDistances.data(1156 * idx + 887);

            auto pa2pb_xyyz_xy = pa2pbDistances.data(1156 * idx + 888);

            auto pa2pb_xyyz_xz = pa2pbDistances.data(1156 * idx + 889);

            auto pa2pb_xyyz_yy = pa2pbDistances.data(1156 * idx + 890);

            auto pa2pb_xyyz_yz = pa2pbDistances.data(1156 * idx + 891);

            auto pa2pb_xyyz_xxxx = pa2pbDistances.data(1156 * idx + 903);

            auto pa2pb_xyyz_xxxy = pa2pbDistances.data(1156 * idx + 904);

            auto pa2pb_xyyz_xxxz = pa2pbDistances.data(1156 * idx + 905);

            auto pa2pb_xyyz_xxyy = pa2pbDistances.data(1156 * idx + 906);

            auto pa2pb_xyyz_xxyz = pa2pbDistances.data(1156 * idx + 907);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyz_xxxx = primBuffer.data(225 * idx + 105);

            auto t_xyyz_xxxy = primBuffer.data(225 * idx + 106);

            auto t_xyyz_xxxz = primBuffer.data(225 * idx + 107);

            auto t_xyyz_xxyy = primBuffer.data(225 * idx + 108);

            auto t_xyyz_xxyz = primBuffer.data(225 * idx + 109);

            // Batch of Integrals (105,110)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_y, \
                                     pa2pb_xy_xx, pa2pb_xyy_x, pa2pb_xyy_xxx, pa2pb_xyy_xxy, pa2pb_xyy_y, \
                                     pa2pb_xyyz_xx, pa2pb_xyyz_xxxx, pa2pb_xyyz_xxxy, pa2pb_xyyz_xxxz, pa2pb_xyyz_xxyy, \
                                     pa2pb_xyyz_xxyz, pa2pb_xyyz_xy, pa2pb_xyyz_xz, pa2pb_xyyz_yy, pa2pb_xyyz_yz, \
                                     pa2pb_xyz_x, pa2pb_xyz_xxx, pa2pb_xyz_xxy, pa2pb_xyz_xxz, pa2pb_xyz_y, \
                                     pa2pb_xyz_z, pa2pb_xz_xx, pa2pb_xz_xxxx, pa2pb_xz_xxxy, pa2pb_xz_xxxz, \
                                     pa2pb_xz_xxyy, pa2pb_xz_xxyz, pa2pb_xz_xy, pa2pb_xz_xz, pa2pb_xz_yy, pa2pb_xz_yz, \
                                     pa2pb_y_x, pa2pb_yy_xx, pa2pb_yy_xy, pa2pb_yyz_x, pa2pb_yyz_xxx, pa2pb_yyz_xxy, \
                                     pa2pb_yyz_xxz, pa2pb_yyz_xyy, pa2pb_yyz_xyz, pa2pb_yyz_y, pa2pb_yyz_z, pa2pb_yz_xx, \
                                     pa2pb_yz_xy, pa2pb_yz_xz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_xxz, \
                                     pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_y, pa2pb_z_z, pa_xy, pa_xyyz, pa_xz, pa_yy, pa_yz, \
                                     pb_xx, pb_xy, r_0_0, s_0_0, t_xyyz_xxxx, t_xyyz_xxxy, t_xyyz_xxxz, t_xyyz_xxyy, \
                                     t_xyyz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyyz_xxxx[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 1.5 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa_xyyz[j] * fl2_fx + 3.0 * pa2pb_yyz_x[j] * fl2_fx + 1.5 * pa2pb_xz_xx[j] * fl2_fx + pa2pb_z_xxx[j] * fl2_fx + 3.0 * pa2pb_xyyz_xx[j] * fl1_fx + 2.0 * pa2pb_yyz_xxx[j] * fl1_fx + 0.5 * pa2pb_xz_xxxx[j] * fl1_fx + pa2pb_xyyz_xxxx[j]);

                t_xyyz_xxxx[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xyyz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xz[j] * fl3_fx * fl1_fz - 6.0 * pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 9.0 * pa_xyyz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_yyz_x[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xyyz_xx[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xyyz_xx[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyz_xxx[j] * fl1_fx * fl1_fz - pa2pb_xz_xxxx[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xxxx[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xxxx[j] * fl1_fz);

                t_xyyz_xxxy[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl3_fx + 0.375 * pa2pb_z_y[j] * fl3_fx + 1.5 * pa2pb_xyz_x[j] * fl2_fx + 0.75 * pa2pb_yyz_y[j] * fl2_fx + 1.5 * pa2pb_yz_xx[j] * fl2_fx + 0.75 * pa2pb_xz_xy[j] * fl2_fx + 0.75 * pa2pb_z_xxy[j] * fl2_fx + 1.5 * pa2pb_xyyz_xy[j] * fl1_fx + pa2pb_xyz_xxx[j] * fl1_fx + 1.5 * pa2pb_yyz_xxy[j] * fl1_fx + 0.5 * pa2pb_xz_xxxy[j] * fl1_fx + pa2pb_xyyz_xxxy[j]);

                t_xyyz_xxxy[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_yz[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xyz_x[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yyz_y[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyz_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_xxx[j] * fl1_fx * fl1_fz + 21.0 * pa2pb_yyz_xxy[j] * fl1_fx * fl1_fz - pa2pb_xz_xxxy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xxxy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xxxy[j] * fl1_fz);

                t_xyyz_xxxz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_xx[j] * fl3_fx + 0.75 * pa2pb_xyy_x[j] * fl2_fx + 0.75 * pa2pb_yyz_z[j] * fl2_fx + 0.75 * pa2pb_yy_xx[j] * fl2_fx + 0.75 * pa2pb_xz_xz[j] * fl2_fx + 0.25 * pa2pb_x_xxx[j] * fl2_fx + 0.75 * pa2pb_z_xxz[j] * fl2_fx + 1.5 * pa2pb_xyyz_xz[j] * fl1_fx + 0.5 * pa2pb_xyy_xxx[j] * fl1_fx + 1.5 * pa2pb_yyz_xxz[j] * fl1_fx + 0.5 * pa2pb_xz_xxxz[j] * fl1_fx + pa2pb_xyyz_xxxz[j]);

                t_xyyz_xxxz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_yy[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xyy_x[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yyz_z[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyz_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyy_xxx[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyz_xxz[j] * fl1_fx * fl1_fz - pa2pb_xz_xxxz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xxxz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xxxz[j] * fl1_fz);

                t_xyyz_xxyy[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 0.75 * pa2pb_z_x[j] * fl3_fx + 0.25 * pa_xyyz[j] * fl2_fx + pa2pb_xyz_y[j] * fl2_fx + 0.75 * pa2pb_xz_xx[j] * fl2_fx + 0.5 * pa2pb_yyz_x[j] * fl2_fx + 2.0 * pa2pb_yz_xy[j] * fl2_fx + 0.25 * pa2pb_xz_yy[j] * fl2_fx + 0.5 * pa2pb_z_xyy[j] * fl2_fx + 0.5 * pa2pb_xyyz_xx[j] * fl1_fx + 0.5 * pa2pb_xyyz_yy[j] * fl1_fx + 2.0 * pa2pb_xyz_xxy[j] * fl1_fx + pa2pb_yyz_xyy[j] * fl1_fx + 0.5 * pa2pb_xz_xxyy[j] * fl1_fx + pa2pb_xyyz_xxyy[j]);

                t_xyyz_xxyy[j] += fl_r_0_0 * (-pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_xz[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 0.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyyz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_xyyz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xyz_y[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyz_x[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyyz_yy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyyz_yy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xyy[j] * fl1_fx * fl1_fz - pa2pb_xz_xxyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xxyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xxyy[j] * fl1_fz);

                t_xyyz_xxyz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl3_fx + 0.5 * pa2pb_y_x[j] * fl3_fx + 0.125 * pa2pb_x_y[j] * fl3_fx + 0.25 * pb_xy[j] * fl3_fx + 0.25 * pa2pb_xyy_y[j] * fl2_fx + 0.5 * pa2pb_xyz_z[j] * fl2_fx + 0.5 * pa2pb_xy_xx[j] * fl2_fx + 0.5 * pa2pb_yy_xy[j] * fl2_fx + pa2pb_yz_xz[j] * fl2_fx + 0.25 * pa2pb_xz_yz[j] * fl2_fx + 0.25 * pa2pb_x_xxy[j] * fl2_fx + 0.5 * pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_xyyz_yz[j] * fl1_fx + 0.5 * pa2pb_xyy_xxy[j] * fl1_fx + pa2pb_xyz_xxz[j] * fl1_fx + pa2pb_yyz_xyz[j] * fl1_fx + 0.5 * pa2pb_xz_xxyz[j] * fl1_fx + pa2pb_xyyz_xxyz[j]);

                t_xyyz_xxyz[j] += fl_r_0_0 * (-0.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 2.5 * pa_xy[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 2.5 * pb_xy[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xyy_y[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xyz_z[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xy_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yy_xy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_xxy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_yz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyz_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyy_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyz_xyz[j] * fl1_fx * fl1_fz - pa2pb_xz_xxyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xxyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_110_115(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_xz_xxzz = pa2pbDistances.data(1156 * idx + 194);

            auto pa2pb_xz_xyyy = pa2pbDistances.data(1156 * idx + 195);

            auto pa2pb_xz_xyyz = pa2pbDistances.data(1156 * idx + 196);

            auto pa2pb_xz_xyzz = pa2pbDistances.data(1156 * idx + 197);

            auto pa2pb_xz_xzzz = pa2pbDistances.data(1156 * idx + 198);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_xyy_x = pa2pbDistances.data(1156 * idx + 408);

            auto pa2pb_xyy_z = pa2pbDistances.data(1156 * idx + 410);

            auto pa2pb_xyy_xxz = pa2pbDistances.data(1156 * idx + 419);

            auto pa2pb_xyy_xyy = pa2pbDistances.data(1156 * idx + 420);

            auto pa2pb_xyy_xyz = pa2pbDistances.data(1156 * idx + 421);

            auto pa2pb_xyy_xzz = pa2pbDistances.data(1156 * idx + 422);

            auto pa2pb_xyz_x = pa2pbDistances.data(1156 * idx + 442);

            auto pa2pb_xyz_xyy = pa2pbDistances.data(1156 * idx + 454);

            auto pa2pb_xyz_xyz = pa2pbDistances.data(1156 * idx + 455);

            auto pa2pb_xyz_xzz = pa2pbDistances.data(1156 * idx + 456);

            auto pa2pb_yyz_x = pa2pbDistances.data(1156 * idx + 544);

            auto pa2pb_yyz_y = pa2pbDistances.data(1156 * idx + 545);

            auto pa2pb_yyz_z = pa2pbDistances.data(1156 * idx + 546);

            auto pa2pb_yyz_xzz = pa2pbDistances.data(1156 * idx + 558);

            auto pa2pb_yyz_yyy = pa2pbDistances.data(1156 * idx + 559);

            auto pa2pb_yyz_yyz = pa2pbDistances.data(1156 * idx + 560);

            auto pa2pb_yyz_yzz = pa2pbDistances.data(1156 * idx + 561);

            auto pa2pb_yyz_zzz = pa2pbDistances.data(1156 * idx + 562);

            auto pa2pb_xyyz_xx = pa2pbDistances.data(1156 * idx + 887);

            auto pa2pb_xyyz_xy = pa2pbDistances.data(1156 * idx + 888);

            auto pa2pb_xyyz_xz = pa2pbDistances.data(1156 * idx + 889);

            auto pa2pb_xyyz_zz = pa2pbDistances.data(1156 * idx + 892);

            auto pa2pb_xyyz_xxzz = pa2pbDistances.data(1156 * idx + 908);

            auto pa2pb_xyyz_xyyy = pa2pbDistances.data(1156 * idx + 909);

            auto pa2pb_xyyz_xyyz = pa2pbDistances.data(1156 * idx + 910);

            auto pa2pb_xyyz_xyzz = pa2pbDistances.data(1156 * idx + 911);

            auto pa2pb_xyyz_xzzz = pa2pbDistances.data(1156 * idx + 912);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyz_xxzz = primBuffer.data(225 * idx + 110);

            auto t_xyyz_xyyy = primBuffer.data(225 * idx + 111);

            auto t_xyyz_xyyz = primBuffer.data(225 * idx + 112);

            auto t_xyyz_xyzz = primBuffer.data(225 * idx + 113);

            auto t_xyyz_xzzz = primBuffer.data(225 * idx + 114);

            // Batch of Integrals (110,115)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxz, pa2pb_x_xyy, pa2pb_x_xyz, \
                                     pa2pb_x_xzz, pa2pb_x_z, pa2pb_xy_xy, pa2pb_xy_xz, pa2pb_xyy_x, pa2pb_xyy_xxz, \
                                     pa2pb_xyy_xyy, pa2pb_xyy_xyz, pa2pb_xyy_xzz, pa2pb_xyy_z, pa2pb_xyyz_xx, \
                                     pa2pb_xyyz_xxzz, pa2pb_xyyz_xy, pa2pb_xyyz_xyyy, pa2pb_xyyz_xyyz, pa2pb_xyyz_xyzz, \
                                     pa2pb_xyyz_xz, pa2pb_xyyz_xzzz, pa2pb_xyyz_zz, pa2pb_xyz_x, pa2pb_xyz_xyy, \
                                     pa2pb_xyz_xyz, pa2pb_xyz_xzz, pa2pb_xz_xx, pa2pb_xz_xxzz, pa2pb_xz_xy, \
                                     pa2pb_xz_xyyy, pa2pb_xz_xyyz, pa2pb_xz_xyzz, pa2pb_xz_xz, pa2pb_xz_xzzz, \
                                     pa2pb_xz_zz, pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, \
                                     pa2pb_yy_zz, pa2pb_yyz_x, pa2pb_yyz_xzz, pa2pb_yyz_y, pa2pb_yyz_yyy, \
                                     pa2pb_yyz_yyz, pa2pb_yyz_yzz, pa2pb_yyz_z, pa2pb_yyz_zzz, pa2pb_yz_yy, pa2pb_yz_yz, \
                                     pa2pb_yz_zz, pa2pb_z_x, pa2pb_z_xzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, \
                                     pa2pb_z_yzz, pa2pb_z_z, pa2pb_z_zzz, pa_xyyz, pa_xz, pa_yy, pa_yz, pb_xz, pb_yy, pb_yz, \
                                     pb_zz, r_0_0, s_0_0, t_xyyz_xxzz, t_xyyz_xyyy, t_xyyz_xyyz, t_xyyz_xyzz, \
                                     t_xyyz_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyyz_xxzz[j] = fl_s_0_0 * (0.125 * pa_xz[j] * fl3_fx + 0.25 * pa2pb_x_z[j] * fl3_fx + 0.25 * pa2pb_z_x[j] * fl3_fx + 0.5 * pb_xz[j] * fl3_fx + 0.25 * pa_xyyz[j] * fl2_fx + 0.5 * pa2pb_xyy_z[j] * fl2_fx + 0.5 * pa2pb_yyz_x[j] * fl2_fx + pa2pb_yy_xz[j] * fl2_fx + 0.25 * pa2pb_xz_xx[j] * fl2_fx + 0.25 * pa2pb_xz_zz[j] * fl2_fx + 0.5 * pa2pb_x_xxz[j] * fl2_fx + 0.5 * pa2pb_z_xzz[j] * fl2_fx + 0.5 * pa2pb_xyyz_xx[j] * fl1_fx + 0.5 * pa2pb_xyyz_zz[j] * fl1_fx + pa2pb_xyy_xxz[j] * fl1_fx + pa2pb_yyz_xzz[j] * fl1_fx + 0.5 * pa2pb_xz_xxzz[j] * fl1_fx + pa2pb_xyyz_xxzz[j]);

                t_xyyz_xxzz[j] += fl_r_0_0 * (-0.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_xz[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz - pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 5.0 * pb_xz[j] * fl3_fx * fl1_fz + 3.0 * pa_xyyz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xyy_z[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yyz_x[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_yy_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyyz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyyz_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyz_xzz[j] * fl1_fx * fl1_fz - pa2pb_xz_xxzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xxzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xxzz[j] * fl1_fz);

                t_xyyz_xyyy[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl3_fx + 1.125 * pa2pb_z_y[j] * fl3_fx + 1.5 * pa2pb_xyz_x[j] * fl2_fx + 2.25 * pa2pb_xz_xy[j] * fl2_fx + 0.75 * pa2pb_yyz_y[j] * fl2_fx + 1.5 * pa2pb_yz_yy[j] * fl2_fx + 0.25 * pa2pb_z_yyy[j] * fl2_fx + 1.5 * pa2pb_xyyz_xy[j] * fl1_fx + 3.0 * pa2pb_xyz_xyy[j] * fl1_fx + 0.5 * pa2pb_yyz_yyy[j] * fl1_fx + 0.5 * pa2pb_xz_xyyy[j] * fl1_fx + pa2pb_xyyz_xyyy[j]);

                t_xyyz_xyyy[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_yz[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_x[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yyz_y[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyz_xy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyz_xy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyz_xyy[j] * fl1_fx * fl1_fz + 7.0 * pa2pb_yyz_yyy[j] * fl1_fx * fl1_fz - pa2pb_xz_xyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xyyy[j] * fl1_fz);

                t_xyyz_xyyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.125 * pa_yy[j] * fl3_fx + 0.5 * pa2pb_y_y[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.125 * pb_yy[j] * fl3_fx + 0.25 * pa2pb_xyy_x[j] * fl2_fx + pa2pb_xy_xy[j] * fl2_fx + 0.75 * pa2pb_xz_xz[j] * fl2_fx + 0.25 * pa2pb_yyz_z[j] * fl2_fx + 0.25 * pa2pb_yy_yy[j] * fl2_fx + pa2pb_yz_yz[j] * fl2_fx + 0.25 * pa2pb_x_xyy[j] * fl2_fx + 0.25 * pa2pb_z_yyz[j] * fl2_fx + 0.5 * pa2pb_xyyz_xz[j] * fl1_fx + 0.5 * pa2pb_xyy_xyy[j] * fl1_fx + 2.0 * pa2pb_xyz_xyz[j] * fl1_fx + 0.5 * pa2pb_yyz_yyz[j] * fl1_fx + 0.5 * pa2pb_xz_xyyz[j] * fl1_fx + pa2pb_xyyz_xyyz[j]);

                t_xyyz_xyyz[j] += fl_r_0_0 * (1.5 * fl4_fx * fl1_fz - 0.125 * fl3_fx * fl1_fz * fl1_fgb - 0.125 * fl3_fx * fl1_fz * fl1_fga - 0.25 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 1.25 * pa_yy[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_yy[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_xyy_x[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xy_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yyz_z[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yy_yy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_xyy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_xz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyz_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyy_xyy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyz_xyz[j] * fl1_fx * fl1_fz + 7.0 * pa2pb_yyz_yyz[j] * fl1_fx * fl1_fz - pa2pb_xz_xyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xyyz[j] * fl1_fz);

                t_xyyz_xyzz[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl3_fx + 0.5 * pa2pb_y_z[j] * fl3_fx + 0.125 * pa2pb_z_y[j] * fl3_fx + 0.25 * pb_yz[j] * fl3_fx + 0.5 * pa2pb_xyz_x[j] * fl2_fx + pa2pb_xy_xz[j] * fl2_fx + 0.25 * pa2pb_yyz_y[j] * fl2_fx + 0.5 * pa2pb_yy_yz[j] * fl2_fx + 0.5 * pa2pb_yz_zz[j] * fl2_fx + 0.25 * pa2pb_xz_xy[j] * fl2_fx + 0.5 * pa2pb_x_xyz[j] * fl2_fx + 0.25 * pa2pb_z_yzz[j] * fl2_fx + 0.5 * pa2pb_xyyz_xy[j] * fl1_fx + pa2pb_xyy_xyz[j] * fl1_fx + pa2pb_xyz_xzz[j] * fl1_fx + 0.5 * pa2pb_yyz_yzz[j] * fl1_fx + 0.5 * pa2pb_xz_xyzz[j] * fl1_fx + pa2pb_xyyz_xyzz[j]);

                t_xyyz_xyzz[j] += fl_r_0_0 * (-0.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 2.5 * pa_yz[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 2.5 * pb_yz[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_xyz_x[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xy_xz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yyz_y[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yy_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_xyz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_xy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyz_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_xzz[j] * fl1_fx * fl1_fz + 7.0 * pa2pb_yyz_yzz[j] * fl1_fx * fl1_fz - pa2pb_xz_xyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xyzz[j] * fl1_fz);

                t_xyyz_xzzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 0.75 * pa2pb_xyy_x[j] * fl2_fx + 0.75 * pa2pb_yyz_z[j] * fl2_fx + 0.75 * pa2pb_yy_zz[j] * fl2_fx + 0.75 * pa2pb_xz_xz[j] * fl2_fx + 0.75 * pa2pb_x_xzz[j] * fl2_fx + 0.25 * pa2pb_z_zzz[j] * fl2_fx + 1.5 * pa2pb_xyyz_xz[j] * fl1_fx + 1.5 * pa2pb_xyy_xzz[j] * fl1_fx + 0.5 * pa2pb_yyz_zzz[j] * fl1_fx + 0.5 * pa2pb_xz_xzzz[j] * fl1_fx + pa2pb_xyyz_xzzz[j]);

                t_xyyz_xzzz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_yy[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xyy_x[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yyz_z[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyz_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xyy_xzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyz_zzz[j] * fl1_fx * fl1_fz - pa2pb_xz_xzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_xzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_115_120(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_xz_yyyy = pa2pbDistances.data(1156 * idx + 199);

            auto pa2pb_xz_yyyz = pa2pbDistances.data(1156 * idx + 200);

            auto pa2pb_xz_yyzz = pa2pbDistances.data(1156 * idx + 201);

            auto pa2pb_xz_yzzz = pa2pbDistances.data(1156 * idx + 202);

            auto pa2pb_xz_zzzz = pa2pbDistances.data(1156 * idx + 203);

            auto pa2pb_xyy_y = pa2pbDistances.data(1156 * idx + 409);

            auto pa2pb_xyy_z = pa2pbDistances.data(1156 * idx + 410);

            auto pa2pb_xyy_yyy = pa2pbDistances.data(1156 * idx + 423);

            auto pa2pb_xyy_yyz = pa2pbDistances.data(1156 * idx + 424);

            auto pa2pb_xyy_yzz = pa2pbDistances.data(1156 * idx + 425);

            auto pa2pb_xyy_zzz = pa2pbDistances.data(1156 * idx + 426);

            auto pa2pb_xyz_y = pa2pbDistances.data(1156 * idx + 443);

            auto pa2pb_xyz_z = pa2pbDistances.data(1156 * idx + 444);

            auto pa2pb_xyz_yyy = pa2pbDistances.data(1156 * idx + 457);

            auto pa2pb_xyz_yyz = pa2pbDistances.data(1156 * idx + 458);

            auto pa2pb_xyz_yzz = pa2pbDistances.data(1156 * idx + 459);

            auto pa2pb_xyz_zzz = pa2pbDistances.data(1156 * idx + 460);

            auto pa2pb_xyyz_yy = pa2pbDistances.data(1156 * idx + 890);

            auto pa2pb_xyyz_yz = pa2pbDistances.data(1156 * idx + 891);

            auto pa2pb_xyyz_zz = pa2pbDistances.data(1156 * idx + 892);

            auto pa2pb_xyyz_yyyy = pa2pbDistances.data(1156 * idx + 913);

            auto pa2pb_xyyz_yyyz = pa2pbDistances.data(1156 * idx + 914);

            auto pa2pb_xyyz_yyzz = pa2pbDistances.data(1156 * idx + 915);

            auto pa2pb_xyyz_yzzz = pa2pbDistances.data(1156 * idx + 916);

            auto pa2pb_xyyz_zzzz = pa2pbDistances.data(1156 * idx + 917);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyz_yyyy = primBuffer.data(225 * idx + 115);

            auto t_xyyz_yyyz = primBuffer.data(225 * idx + 116);

            auto t_xyyz_yyzz = primBuffer.data(225 * idx + 117);

            auto t_xyyz_yzzz = primBuffer.data(225 * idx + 118);

            auto t_xyyz_zzzz = primBuffer.data(225 * idx + 119);

            // Batch of Integrals (115,120)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, \
                                     pa2pb_x_z, pa2pb_x_zzz, pa2pb_xy_yy, pa2pb_xy_yz, pa2pb_xy_zz, pa2pb_xyy_y, \
                                     pa2pb_xyy_yyy, pa2pb_xyy_yyz, pa2pb_xyy_yzz, pa2pb_xyy_z, pa2pb_xyy_zzz, \
                                     pa2pb_xyyz_yy, pa2pb_xyyz_yyyy, pa2pb_xyyz_yyyz, pa2pb_xyyz_yyzz, pa2pb_xyyz_yz, \
                                     pa2pb_xyyz_yzzz, pa2pb_xyyz_zz, pa2pb_xyyz_zzzz, pa2pb_xyz_y, pa2pb_xyz_yyy, \
                                     pa2pb_xyz_yyz, pa2pb_xyz_yzz, pa2pb_xyz_z, pa2pb_xyz_zzz, pa2pb_xz_yy, \
                                     pa2pb_xz_yyyy, pa2pb_xz_yyyz, pa2pb_xz_yyzz, pa2pb_xz_yz, pa2pb_xz_yzzz, \
                                     pa2pb_xz_zz, pa2pb_xz_zzzz, pa_xy, pa_xyyz, pa_xz, r_0_0, s_0_0, t_xyyz_yyyy, \
                                     t_xyyz_yyyz, t_xyyz_yyzz, t_xyyz_yzzz, t_xyyz_zzzz: VLX_ALIGN)
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

                t_xyyz_yyyy[j] = fl_s_0_0 * (1.875 * pa_xz[j] * fl3_fx + 0.75 * pa_xyyz[j] * fl2_fx + 6.0 * pa2pb_xyz_y[j] * fl2_fx + 4.5 * pa2pb_xz_yy[j] * fl2_fx + 3.0 * pa2pb_xyyz_yy[j] * fl1_fx + 4.0 * pa2pb_xyz_yyy[j] * fl1_fx + 0.5 * pa2pb_xz_yyyy[j] * fl1_fx + pa2pb_xyyz_yyyy[j]);

                t_xyyz_yyyy[j] += fl_r_0_0 * (-4.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 18.75 * pa_xz[j] * fl3_fx * fl1_fz - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xyyz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_xyyz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_xyz_y[j] * fl2_fx * fl1_fz + 54.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xyyz_yy[j] * fl1_fz * fl1_fgb + 42.0 * pa2pb_xyyz_yy[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xyz_yyy[j] * fl1_fx * fl1_fz - pa2pb_xz_yyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_yyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_yyyy[j] * fl1_fz);

                t_xyyz_yyyz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl3_fx + 1.125 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa2pb_xyy_y[j] * fl2_fx + 1.5 * pa2pb_xyz_z[j] * fl2_fx + 1.5 * pa2pb_xy_yy[j] * fl2_fx + 2.25 * pa2pb_xz_yz[j] * fl2_fx + 0.25 * pa2pb_x_yyy[j] * fl2_fx + 1.5 * pa2pb_xyyz_yz[j] * fl1_fx + 0.5 * pa2pb_xyy_yyy[j] * fl1_fx + 3.0 * pa2pb_xyz_yyz[j] * fl1_fx + 0.5 * pa2pb_xz_yyyz[j] * fl1_fx + pa2pb_xyyz_yyyz[j]);

                t_xyyz_yyyz[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_xy[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa2pb_xyy_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyz_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xy_yy[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyyz_yz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyz_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyy_yyy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyz_yyz[j] * fl1_fx * fl1_fz - pa2pb_xz_yyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_yyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_yyyz[j] * fl1_fz);

                t_xyyz_yyzz[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 0.75 * pa2pb_x_z[j] * fl3_fx + 0.25 * pa_xyyz[j] * fl2_fx + 0.5 * pa2pb_xyy_z[j] * fl2_fx + pa2pb_xyz_y[j] * fl2_fx + 2.0 * pa2pb_xy_yz[j] * fl2_fx + 0.75 * pa2pb_xz_zz[j] * fl2_fx + 0.25 * pa2pb_xz_yy[j] * fl2_fx + 0.5 * pa2pb_x_yyz[j] * fl2_fx + 0.5 * pa2pb_xyyz_yy[j] * fl1_fx + 0.5 * pa2pb_xyyz_zz[j] * fl1_fx + pa2pb_xyy_yyz[j] * fl1_fx + 2.0 * pa2pb_xyz_yzz[j] * fl1_fx + 0.5 * pa2pb_xz_yyzz[j] * fl1_fx + pa2pb_xyyz_yyzz[j]);

                t_xyyz_yyzz[j] += fl_r_0_0 * (-pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_xz[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_xyyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_xyyz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xyy_z[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xyz_y[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_xy_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xyyz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xyyz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyyz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyyz_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyy_yyz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyz_yzz[j] * fl1_fx * fl1_fz - pa2pb_xz_yyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_yyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_yyzz[j] * fl1_fz);

                t_xyyz_yzzz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl3_fx + 0.375 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa2pb_xyy_y[j] * fl2_fx + 1.5 * pa2pb_xyz_z[j] * fl2_fx + 1.5 * pa2pb_xy_zz[j] * fl2_fx + 0.75 * pa2pb_xz_yz[j] * fl2_fx + 0.75 * pa2pb_x_yzz[j] * fl2_fx + 1.5 * pa2pb_xyyz_yz[j] * fl1_fx + 1.5 * pa2pb_xyy_yzz[j] * fl1_fx + pa2pb_xyz_zzz[j] * fl1_fx + 0.5 * pa2pb_xz_yzzz[j] * fl1_fx + pa2pb_xyyz_yzzz[j]);

                t_xyyz_yzzz[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_xy[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_xyy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xyy_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyz_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xy_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_yzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyyz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyyz_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xyy_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_zzz[j] * fl1_fx * fl1_fz - pa2pb_xz_yzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_yzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_yzzz[j] * fl1_fz);

                t_xyyz_zzzz[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 1.5 * pa2pb_x_z[j] * fl3_fx + 0.75 * pa_xyyz[j] * fl2_fx + 3.0 * pa2pb_xyy_z[j] * fl2_fx + 1.5 * pa2pb_xz_zz[j] * fl2_fx + pa2pb_x_zzz[j] * fl2_fx + 3.0 * pa2pb_xyyz_zz[j] * fl1_fx + 2.0 * pa2pb_xyy_zzz[j] * fl1_fx + 0.5 * pa2pb_xz_zzzz[j] * fl1_fx + pa2pb_xyyz_zzzz[j]);

                t_xyyz_zzzz[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_x_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xyyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_xyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xz[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 9.0 * pa_xyyz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xyy_z[j] * fl1_fz * fl2_fx - 3.0 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xyyz_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xyyz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyy_zzz[j] * fl1_fz * fl1_fx - pa2pb_xz_zzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xz_zzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_xyyz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_120_125(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_xxxx = pa2pbDistances.data(1156 * idx + 155);

            auto pa2pb_xy_xxxy = pa2pbDistances.data(1156 * idx + 156);

            auto pa2pb_xy_xxxz = pa2pbDistances.data(1156 * idx + 157);

            auto pa2pb_xy_xxyy = pa2pbDistances.data(1156 * idx + 158);

            auto pa2pb_xy_xxyz = pa2pbDistances.data(1156 * idx + 159);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_xyz_x = pa2pbDistances.data(1156 * idx + 442);

            auto pa2pb_xyz_y = pa2pbDistances.data(1156 * idx + 443);

            auto pa2pb_xyz_xxx = pa2pbDistances.data(1156 * idx + 451);

            auto pa2pb_xyz_xxy = pa2pbDistances.data(1156 * idx + 452);

            auto pa2pb_xzz_x = pa2pbDistances.data(1156 * idx + 476);

            auto pa2pb_xzz_y = pa2pbDistances.data(1156 * idx + 477);

            auto pa2pb_xzz_z = pa2pbDistances.data(1156 * idx + 478);

            auto pa2pb_xzz_xxx = pa2pbDistances.data(1156 * idx + 485);

            auto pa2pb_xzz_xxy = pa2pbDistances.data(1156 * idx + 486);

            auto pa2pb_xzz_xxz = pa2pbDistances.data(1156 * idx + 487);

            auto pa2pb_yzz_x = pa2pbDistances.data(1156 * idx + 578);

            auto pa2pb_yzz_y = pa2pbDistances.data(1156 * idx + 579);

            auto pa2pb_yzz_z = pa2pbDistances.data(1156 * idx + 580);

            auto pa2pb_yzz_xxx = pa2pbDistances.data(1156 * idx + 587);

            auto pa2pb_yzz_xxy = pa2pbDistances.data(1156 * idx + 588);

            auto pa2pb_yzz_xxz = pa2pbDistances.data(1156 * idx + 589);

            auto pa2pb_yzz_xyy = pa2pbDistances.data(1156 * idx + 590);

            auto pa2pb_yzz_xyz = pa2pbDistances.data(1156 * idx + 591);

            auto pa2pb_xyzz_xx = pa2pbDistances.data(1156 * idx + 921);

            auto pa2pb_xyzz_xy = pa2pbDistances.data(1156 * idx + 922);

            auto pa2pb_xyzz_xz = pa2pbDistances.data(1156 * idx + 923);

            auto pa2pb_xyzz_yy = pa2pbDistances.data(1156 * idx + 924);

            auto pa2pb_xyzz_yz = pa2pbDistances.data(1156 * idx + 925);

            auto pa2pb_xyzz_xxxx = pa2pbDistances.data(1156 * idx + 937);

            auto pa2pb_xyzz_xxxy = pa2pbDistances.data(1156 * idx + 938);

            auto pa2pb_xyzz_xxxz = pa2pbDistances.data(1156 * idx + 939);

            auto pa2pb_xyzz_xxyy = pa2pbDistances.data(1156 * idx + 940);

            auto pa2pb_xyzz_xxyz = pa2pbDistances.data(1156 * idx + 941);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyzz_xxxx = primBuffer.data(225 * idx + 120);

            auto t_xyzz_xxxy = primBuffer.data(225 * idx + 121);

            auto t_xyzz_xxxz = primBuffer.data(225 * idx + 122);

            auto t_xyzz_xxyy = primBuffer.data(225 * idx + 123);

            auto t_xyzz_xxyz = primBuffer.data(225 * idx + 124);

            // Batch of Integrals (120,125)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_xxz, \
                                     pa2pb_x_y, pa2pb_x_z, pa2pb_xy_xx, pa2pb_xy_xxxx, pa2pb_xy_xxxy, pa2pb_xy_xxxz, \
                                     pa2pb_xy_xxyy, pa2pb_xy_xxyz, pa2pb_xy_xy, pa2pb_xy_xz, pa2pb_xy_yy, pa2pb_xy_yz, \
                                     pa2pb_xyz_x, pa2pb_xyz_xxx, pa2pb_xyz_xxy, pa2pb_xyz_y, pa2pb_xyzz_xx, \
                                     pa2pb_xyzz_xxxx, pa2pb_xyzz_xxxy, pa2pb_xyzz_xxxz, pa2pb_xyzz_xxyy, pa2pb_xyzz_xxyz, \
                                     pa2pb_xyzz_xy, pa2pb_xyzz_xz, pa2pb_xyzz_yy, pa2pb_xyzz_yz, pa2pb_xz_xx, \
                                     pa2pb_xzz_x, pa2pb_xzz_xxx, pa2pb_xzz_xxy, pa2pb_xzz_xxz, pa2pb_xzz_y, \
                                     pa2pb_xzz_z, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, pa2pb_y_xyy, \
                                     pa2pb_y_xyz, pa2pb_y_y, pa2pb_y_z, pa2pb_yz_xx, pa2pb_yz_xy, pa2pb_yzz_x, \
                                     pa2pb_yzz_xxx, pa2pb_yzz_xxy, pa2pb_yzz_xxz, pa2pb_yzz_xyy, pa2pb_yzz_xyz, \
                                     pa2pb_yzz_y, pa2pb_yzz_z, pa2pb_z_x, pa2pb_zz_xx, pa2pb_zz_xy, pa2pb_zz_xz, pa_xy, \
                                     pa_xyzz, pa_xz, pa_yz, pa_zz, pb_xx, pb_xy, pb_xz, r_0_0, s_0_0, t_xyzz_xxxx, \
                                     t_xyzz_xxxy, t_xyzz_xxxz, t_xyzz_xxyy, t_xyzz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyzz_xxxx[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 1.5 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa_xyzz[j] * fl2_fx + 3.0 * pa2pb_yzz_x[j] * fl2_fx + 1.5 * pa2pb_xy_xx[j] * fl2_fx + pa2pb_y_xxx[j] * fl2_fx + 3.0 * pa2pb_xyzz_xx[j] * fl1_fx + 2.0 * pa2pb_yzz_xxx[j] * fl1_fx + 0.5 * pa2pb_xy_xxxx[j] * fl1_fx + pa2pb_xyzz_xxxx[j]);

                t_xyzz_xxxx[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xyzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xy[j] * fl1_fz * fl3_fx - 6.0 * pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 9.0 * pa_xyzz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_yzz_x[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xyzz_xx[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xy_xx[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xyzz_xx[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yzz_xxx[j] * fl1_fx * fl1_fz - pa2pb_xy_xxxx[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xxxx[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xxxx[j] * fl1_fz);

                t_xyzz_xxxy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_zz[j] * fl3_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pb_xx[j] * fl3_fx + 0.75 * pa2pb_xzz_x[j] * fl2_fx + 0.75 * pa2pb_yzz_y[j] * fl2_fx + 0.75 * pa2pb_zz_xx[j] * fl2_fx + 0.75 * pa2pb_xy_xy[j] * fl2_fx + 0.25 * pa2pb_x_xxx[j] * fl2_fx + 0.75 * pa2pb_y_xxy[j] * fl2_fx + 1.5 * pa2pb_xyzz_xy[j] * fl1_fx + 0.5 * pa2pb_xzz_xxx[j] * fl1_fx + 1.5 * pa2pb_yzz_xxy[j] * fl1_fx + 0.5 * pa2pb_xy_xxxy[j] * fl1_fx + pa2pb_xyzz_xxxy[j]);

                t_xyzz_xxxy[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_zz[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xzz_x[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yzz_y[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyzz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xy[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyzz_xy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xzz_xxx[j] * fl1_fx * fl1_fz + 21.0 * pa2pb_yzz_xxy[j] * fl1_fx * fl1_fz - pa2pb_xy_xxxy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xxxy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xxxy[j] * fl1_fz);

                t_xyzz_xxxz[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl3_fx + 0.375 * pa2pb_y_z[j] * fl3_fx + 1.5 * pa2pb_xyz_x[j] * fl2_fx + 0.75 * pa2pb_yzz_z[j] * fl2_fx + 1.5 * pa2pb_yz_xx[j] * fl2_fx + 0.75 * pa2pb_xy_xz[j] * fl2_fx + 0.75 * pa2pb_y_xxz[j] * fl2_fx + 1.5 * pa2pb_xyzz_xz[j] * fl1_fx + pa2pb_xyz_xxx[j] * fl1_fx + 1.5 * pa2pb_yzz_xxz[j] * fl1_fx + 0.5 * pa2pb_xy_xxxz[j] * fl1_fx + pa2pb_xyzz_xxxz[j]);

                t_xyzz_xxxz[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_yz[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xyz_x[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yzz_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyzz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyzz_xz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_xxx[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yzz_xxz[j] * fl1_fx * fl1_fz - pa2pb_xy_xxxz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xxxz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xxxz[j] * fl1_fz);

                t_xyzz_xxyy[j] = fl_s_0_0 * (0.125 * pa_xy[j] * fl3_fx + 0.25 * pa2pb_x_y[j] * fl3_fx + 0.25 * pa2pb_y_x[j] * fl3_fx + 0.5 * pb_xy[j] * fl3_fx + 0.25 * pa_xyzz[j] * fl2_fx + 0.5 * pa2pb_xzz_y[j] * fl2_fx + 0.5 * pa2pb_yzz_x[j] * fl2_fx + pa2pb_zz_xy[j] * fl2_fx + 0.25 * pa2pb_xy_xx[j] * fl2_fx + 0.25 * pa2pb_xy_yy[j] * fl2_fx + 0.5 * pa2pb_x_xxy[j] * fl2_fx + 0.5 * pa2pb_y_xyy[j] * fl2_fx + 0.5 * pa2pb_xyzz_xx[j] * fl1_fx + 0.5 * pa2pb_xyzz_yy[j] * fl1_fx + pa2pb_xzz_xxy[j] * fl1_fx + pa2pb_yzz_xyy[j] * fl1_fx + 0.5 * pa2pb_xy_xxyy[j] * fl1_fx + pa2pb_xyzz_xxyy[j]);

                t_xyzz_xxyy[j] += fl_r_0_0 * (-0.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyzz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_xy[j] * fl1_fz * fl3_fx - pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz - pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 5.0 * pb_xy[j] * fl3_fx * fl1_fz + 3.0 * pa_xyzz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xzz_y[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzz_x[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_y_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyzz_yy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xy_xx[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_xy_yy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyzz_yy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yzz_xyy[j] * fl1_fx * fl1_fz - pa2pb_xy_xxyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xxyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xxyy[j] * fl1_fz);

                t_xyzz_xxyz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl3_fx + 0.5 * pa2pb_z_x[j] * fl3_fx + 0.125 * pa2pb_x_z[j] * fl3_fx + 0.25 * pb_xz[j] * fl3_fx + 0.5 * pa2pb_xyz_y[j] * fl2_fx + 0.25 * pa2pb_xzz_z[j] * fl2_fx + 0.5 * pa2pb_xz_xx[j] * fl2_fx + pa2pb_yz_xy[j] * fl2_fx + 0.5 * pa2pb_zz_xz[j] * fl2_fx + 0.25 * pa2pb_xy_yz[j] * fl2_fx + 0.25 * pa2pb_x_xxz[j] * fl2_fx + 0.5 * pa2pb_y_xyz[j] * fl2_fx + 0.5 * pa2pb_xyzz_yz[j] * fl1_fx + pa2pb_xyz_xxy[j] * fl1_fx + 0.5 * pa2pb_xzz_xxz[j] * fl1_fx + pa2pb_yzz_xyz[j] * fl1_fx + 0.5 * pa2pb_xy_xxyz[j] * fl1_fx + pa2pb_xyzz_xxyz[j]);

                t_xyzz_xxyz[j] += fl_r_0_0 * (-0.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 2.5 * pa_xz[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 2.5 * pb_xz[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_xyz_y[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_xzz_z[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_y_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_yz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xy_yz[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyzz_yz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_xxy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xzz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yzz_xyz[j] * fl1_fx * fl1_fz - pa2pb_xy_xxyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xxyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_125_130(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_xy_xx = pa2pbDistances.data(1156 * idx + 139);

            auto pa2pb_xy_xy = pa2pbDistances.data(1156 * idx + 140);

            auto pa2pb_xy_xz = pa2pbDistances.data(1156 * idx + 141);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_xy_xxzz = pa2pbDistances.data(1156 * idx + 160);

            auto pa2pb_xy_xyyy = pa2pbDistances.data(1156 * idx + 161);

            auto pa2pb_xy_xyyz = pa2pbDistances.data(1156 * idx + 162);

            auto pa2pb_xy_xyzz = pa2pbDistances.data(1156 * idx + 163);

            auto pa2pb_xy_xzzz = pa2pbDistances.data(1156 * idx + 164);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_xyz_x = pa2pbDistances.data(1156 * idx + 442);

            auto pa2pb_xyz_z = pa2pbDistances.data(1156 * idx + 444);

            auto pa2pb_xyz_xxz = pa2pbDistances.data(1156 * idx + 453);

            auto pa2pb_xyz_xyy = pa2pbDistances.data(1156 * idx + 454);

            auto pa2pb_xyz_xyz = pa2pbDistances.data(1156 * idx + 455);

            auto pa2pb_xyz_xzz = pa2pbDistances.data(1156 * idx + 456);

            auto pa2pb_xzz_x = pa2pbDistances.data(1156 * idx + 476);

            auto pa2pb_xzz_xyy = pa2pbDistances.data(1156 * idx + 488);

            auto pa2pb_xzz_xyz = pa2pbDistances.data(1156 * idx + 489);

            auto pa2pb_xzz_xzz = pa2pbDistances.data(1156 * idx + 490);

            auto pa2pb_yzz_x = pa2pbDistances.data(1156 * idx + 578);

            auto pa2pb_yzz_y = pa2pbDistances.data(1156 * idx + 579);

            auto pa2pb_yzz_z = pa2pbDistances.data(1156 * idx + 580);

            auto pa2pb_yzz_xzz = pa2pbDistances.data(1156 * idx + 592);

            auto pa2pb_yzz_yyy = pa2pbDistances.data(1156 * idx + 593);

            auto pa2pb_yzz_yyz = pa2pbDistances.data(1156 * idx + 594);

            auto pa2pb_yzz_yzz = pa2pbDistances.data(1156 * idx + 595);

            auto pa2pb_yzz_zzz = pa2pbDistances.data(1156 * idx + 596);

            auto pa2pb_xyzz_xx = pa2pbDistances.data(1156 * idx + 921);

            auto pa2pb_xyzz_xy = pa2pbDistances.data(1156 * idx + 922);

            auto pa2pb_xyzz_xz = pa2pbDistances.data(1156 * idx + 923);

            auto pa2pb_xyzz_zz = pa2pbDistances.data(1156 * idx + 926);

            auto pa2pb_xyzz_xxzz = pa2pbDistances.data(1156 * idx + 942);

            auto pa2pb_xyzz_xyyy = pa2pbDistances.data(1156 * idx + 943);

            auto pa2pb_xyzz_xyyz = pa2pbDistances.data(1156 * idx + 944);

            auto pa2pb_xyzz_xyzz = pa2pbDistances.data(1156 * idx + 945);

            auto pa2pb_xyzz_xzzz = pa2pbDistances.data(1156 * idx + 946);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyzz_xxzz = primBuffer.data(225 * idx + 125);

            auto t_xyzz_xyyy = primBuffer.data(225 * idx + 126);

            auto t_xyzz_xyyz = primBuffer.data(225 * idx + 127);

            auto t_xyzz_xyzz = primBuffer.data(225 * idx + 128);

            auto t_xyzz_xzzz = primBuffer.data(225 * idx + 129);

            // Batch of Integrals (125,130)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xyy, pa2pb_x_xyz, pa2pb_x_xzz, \
                                     pa2pb_xy_xx, pa2pb_xy_xxzz, pa2pb_xy_xy, pa2pb_xy_xyyy, pa2pb_xy_xyyz, \
                                     pa2pb_xy_xyzz, pa2pb_xy_xz, pa2pb_xy_xzzz, pa2pb_xy_zz, pa2pb_xyz_x, pa2pb_xyz_xxz, \
                                     pa2pb_xyz_xyy, pa2pb_xyz_xyz, pa2pb_xyz_xzz, pa2pb_xyz_z, pa2pb_xyzz_xx, \
                                     pa2pb_xyzz_xxzz, pa2pb_xyzz_xy, pa2pb_xyzz_xyyy, pa2pb_xyzz_xyyz, pa2pb_xyzz_xyzz, \
                                     pa2pb_xyzz_xz, pa2pb_xyzz_xzzz, pa2pb_xyzz_zz, pa2pb_xz_xy, pa2pb_xz_xz, \
                                     pa2pb_xzz_x, pa2pb_xzz_xyy, pa2pb_xzz_xyz, pa2pb_xzz_xzz, pa2pb_y_x, pa2pb_y_xzz, \
                                     pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, pa2pb_y_z, pa2pb_y_zzz, \
                                     pa2pb_yz_xz, pa2pb_yz_yy, pa2pb_yz_yz, pa2pb_yz_zz, pa2pb_yzz_x, pa2pb_yzz_xzz, \
                                     pa2pb_yzz_y, pa2pb_yzz_yyy, pa2pb_yzz_yyz, pa2pb_yzz_yzz, pa2pb_yzz_z, \
                                     pa2pb_yzz_zzz, pa2pb_z_y, pa2pb_z_z, pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zz_zz, pa_xy, \
                                     pa_xyzz, pa_yz, pa_zz, pb_yy, pb_yz, pb_zz, r_0_0, s_0_0, t_xyzz_xxzz, t_xyzz_xyyy, \
                                     t_xyzz_xyyz, t_xyzz_xyzz, t_xyzz_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyzz_xxzz[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 0.75 * pa2pb_y_x[j] * fl3_fx + 0.25 * pa_xyzz[j] * fl2_fx + pa2pb_xyz_z[j] * fl2_fx + 0.75 * pa2pb_xy_xx[j] * fl2_fx + 0.5 * pa2pb_yzz_x[j] * fl2_fx + 2.0 * pa2pb_yz_xz[j] * fl2_fx + 0.25 * pa2pb_xy_zz[j] * fl2_fx + 0.5 * pa2pb_y_xzz[j] * fl2_fx + 0.5 * pa2pb_xyzz_xx[j] * fl1_fx + 0.5 * pa2pb_xyzz_zz[j] * fl1_fx + 2.0 * pa2pb_xyz_xxz[j] * fl1_fx + pa2pb_yzz_xzz[j] * fl1_fx + 0.5 * pa2pb_xy_xxzz[j] * fl1_fx + pa2pb_xyzz_xxzz[j]);

                t_xyzz_xxzz[j] += fl_r_0_0 * (-pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_xy[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 0.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_xyzz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xyz_z[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xy_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzz_x[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xyzz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xy_zz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyzz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyz_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xzz[j] * fl1_fx * fl1_fz - pa2pb_xy_xxzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xxzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xxzz[j] * fl1_fz);

                t_xyzz_xyyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_zz[j] * fl3_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 0.75 * pa2pb_xzz_x[j] * fl2_fx + 0.75 * pa2pb_yzz_y[j] * fl2_fx + 0.75 * pa2pb_zz_yy[j] * fl2_fx + 0.75 * pa2pb_xy_xy[j] * fl2_fx + 0.75 * pa2pb_x_xyy[j] * fl2_fx + 0.25 * pa2pb_y_yyy[j] * fl2_fx + 1.5 * pa2pb_xyzz_xy[j] * fl1_fx + 1.5 * pa2pb_xzz_xyy[j] * fl1_fx + 0.5 * pa2pb_yzz_yyy[j] * fl1_fx + 0.5 * pa2pb_xy_xyyy[j] * fl1_fx + pa2pb_xyzz_xyyy[j]);

                t_xyzz_xyyy[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_zz[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xzz_x[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yzz_y[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_y_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyzz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_xy[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyzz_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xzz_xyy[j] * fl1_fx * fl1_fz + 7.0 * pa2pb_yzz_yyy[j] * fl1_fx * fl1_fz - pa2pb_xy_xyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xyyy[j] * fl1_fz);

                t_xyzz_xyyz[j] = fl_s_0_0 * (0.25 * pa_yz[j] * fl3_fx + 0.5 * pa2pb_z_y[j] * fl3_fx + 0.125 * pa2pb_y_z[j] * fl3_fx + 0.25 * pb_yz[j] * fl3_fx + 0.5 * pa2pb_xyz_x[j] * fl2_fx + pa2pb_xz_xy[j] * fl2_fx + 0.25 * pa2pb_yzz_z[j] * fl2_fx + 0.5 * pa2pb_yz_yy[j] * fl2_fx + 0.5 * pa2pb_zz_yz[j] * fl2_fx + 0.25 * pa2pb_xy_xz[j] * fl2_fx + 0.5 * pa2pb_x_xyz[j] * fl2_fx + 0.25 * pa2pb_y_yyz[j] * fl2_fx + 0.5 * pa2pb_xyzz_xz[j] * fl1_fx + pa2pb_xyz_xyy[j] * fl1_fx + pa2pb_xzz_xyz[j] * fl1_fx + 0.5 * pa2pb_yzz_yyz[j] * fl1_fx + 0.5 * pa2pb_xy_xyyz[j] * fl1_fx + pa2pb_xyzz_xyyz[j]);

                t_xyzz_xyyz[j] += fl_r_0_0 * (-0.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 2.5 * pa_yz[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 2.5 * pb_yz[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_xyz_x[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yzz_z[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_y_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_xz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xy_xz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyzz_xz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_xyz[j] * fl1_fx * fl1_fz + 7.0 * pa2pb_yzz_yyz[j] * fl1_fx * fl1_fz - pa2pb_xy_xyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xyyz[j] * fl1_fz);

                t_xyzz_xyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.125 * pa_zz[j] * fl3_fx + 0.5 * pa2pb_z_z[j] * fl3_fx + 0.125 * pb_zz[j] * fl3_fx + 0.75 * pa2pb_xy_xy[j] * fl2_fx + 0.25 * pa2pb_xzz_x[j] * fl2_fx + pa2pb_xz_xz[j] * fl2_fx + 0.25 * pa2pb_yzz_y[j] * fl2_fx + pa2pb_yz_yz[j] * fl2_fx + 0.25 * pa2pb_zz_zz[j] * fl2_fx + 0.25 * pa2pb_x_xzz[j] * fl2_fx + 0.25 * pa2pb_y_yzz[j] * fl2_fx + 0.5 * pa2pb_xyzz_xy[j] * fl1_fx + 2.0 * pa2pb_xyz_xyz[j] * fl1_fx + 0.5 * pa2pb_xzz_xzz[j] * fl1_fx + 0.5 * pa2pb_yzz_yzz[j] * fl1_fx + 0.5 * pa2pb_xy_xyzz[j] * fl1_fx + pa2pb_xyzz_xyzz[j]);

                t_xyzz_xyzz[j] += fl_r_0_0 * (1.5 * fl4_fx * fl1_fz - 0.125 * fl3_fx * fl1_fz * fl1_fgb - 0.125 * fl3_fx * fl1_fz * fl1_fga - 0.25 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 1.25 * pa_zz[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 0.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xy_xy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_xzz_x[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_yzz_y[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_y_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_xy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyzz_xy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyz_xyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xzz_xzz[j] * fl1_fx * fl1_fz + 7.0 * pa2pb_yzz_yzz[j] * fl1_fx * fl1_fz - pa2pb_xy_xyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xyzz[j] * fl1_fz);

                t_xyzz_xzzz[j] = fl_s_0_0 * (0.75 * pa_yz[j] * fl3_fx + 1.125 * pa2pb_y_z[j] * fl3_fx + 1.5 * pa2pb_xyz_x[j] * fl2_fx + 2.25 * pa2pb_xy_xz[j] * fl2_fx + 0.75 * pa2pb_yzz_z[j] * fl2_fx + 1.5 * pa2pb_yz_zz[j] * fl2_fx + 0.25 * pa2pb_y_zzz[j] * fl2_fx + 1.5 * pa2pb_xyzz_xz[j] * fl1_fx + 3.0 * pa2pb_xyz_xzz[j] * fl1_fx + 0.5 * pa2pb_yzz_zzz[j] * fl1_fx + 0.5 * pa2pb_xy_xzzz[j] * fl1_fx + pa2pb_xyzz_xzzz[j]);

                t_xyzz_xzzz[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_yz[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_x[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xy_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yzz_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyzz_xz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyzz_xz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyz_xzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yzz_zzz[j] * fl1_fx * fl1_fz - pa2pb_xy_xzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_xzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_130_135(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_xy_yy = pa2pbDistances.data(1156 * idx + 142);

            auto pa2pb_xy_yz = pa2pbDistances.data(1156 * idx + 143);

            auto pa2pb_xy_zz = pa2pbDistances.data(1156 * idx + 144);

            auto pa2pb_xy_yyyy = pa2pbDistances.data(1156 * idx + 165);

            auto pa2pb_xy_yyyz = pa2pbDistances.data(1156 * idx + 166);

            auto pa2pb_xy_yyzz = pa2pbDistances.data(1156 * idx + 167);

            auto pa2pb_xy_yzzz = pa2pbDistances.data(1156 * idx + 168);

            auto pa2pb_xy_zzzz = pa2pbDistances.data(1156 * idx + 169);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_xyz_y = pa2pbDistances.data(1156 * idx + 443);

            auto pa2pb_xyz_z = pa2pbDistances.data(1156 * idx + 444);

            auto pa2pb_xyz_yyy = pa2pbDistances.data(1156 * idx + 457);

            auto pa2pb_xyz_yyz = pa2pbDistances.data(1156 * idx + 458);

            auto pa2pb_xyz_yzz = pa2pbDistances.data(1156 * idx + 459);

            auto pa2pb_xyz_zzz = pa2pbDistances.data(1156 * idx + 460);

            auto pa2pb_xzz_y = pa2pbDistances.data(1156 * idx + 477);

            auto pa2pb_xzz_z = pa2pbDistances.data(1156 * idx + 478);

            auto pa2pb_xzz_yyy = pa2pbDistances.data(1156 * idx + 491);

            auto pa2pb_xzz_yyz = pa2pbDistances.data(1156 * idx + 492);

            auto pa2pb_xzz_yzz = pa2pbDistances.data(1156 * idx + 493);

            auto pa2pb_xzz_zzz = pa2pbDistances.data(1156 * idx + 494);

            auto pa2pb_xyzz_yy = pa2pbDistances.data(1156 * idx + 924);

            auto pa2pb_xyzz_yz = pa2pbDistances.data(1156 * idx + 925);

            auto pa2pb_xyzz_zz = pa2pbDistances.data(1156 * idx + 926);

            auto pa2pb_xyzz_yyyy = pa2pbDistances.data(1156 * idx + 947);

            auto pa2pb_xyzz_yyyz = pa2pbDistances.data(1156 * idx + 948);

            auto pa2pb_xyzz_yyzz = pa2pbDistances.data(1156 * idx + 949);

            auto pa2pb_xyzz_yzzz = pa2pbDistances.data(1156 * idx + 950);

            auto pa2pb_xyzz_zzzz = pa2pbDistances.data(1156 * idx + 951);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyzz_yyyy = primBuffer.data(225 * idx + 130);

            auto t_xyzz_yyyz = primBuffer.data(225 * idx + 131);

            auto t_xyzz_yyzz = primBuffer.data(225 * idx + 132);

            auto t_xyzz_yzzz = primBuffer.data(225 * idx + 133);

            auto t_xyzz_zzzz = primBuffer.data(225 * idx + 134);

            // Batch of Integrals (130,135)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, \
                                     pa2pb_x_z, pa2pb_x_zzz, pa2pb_xy_yy, pa2pb_xy_yyyy, pa2pb_xy_yyyz, \
                                     pa2pb_xy_yyzz, pa2pb_xy_yz, pa2pb_xy_yzzz, pa2pb_xy_zz, pa2pb_xy_zzzz, pa2pb_xyz_y, \
                                     pa2pb_xyz_yyy, pa2pb_xyz_yyz, pa2pb_xyz_yzz, pa2pb_xyz_z, pa2pb_xyz_zzz, \
                                     pa2pb_xyzz_yy, pa2pb_xyzz_yyyy, pa2pb_xyzz_yyyz, pa2pb_xyzz_yyzz, pa2pb_xyzz_yz, \
                                     pa2pb_xyzz_yzzz, pa2pb_xyzz_zz, pa2pb_xyzz_zzzz, pa2pb_xz_yy, pa2pb_xz_yz, \
                                     pa2pb_xz_zz, pa2pb_xzz_y, pa2pb_xzz_yyy, pa2pb_xzz_yyz, pa2pb_xzz_yzz, \
                                     pa2pb_xzz_z, pa2pb_xzz_zzz, pa_xy, pa_xyzz, pa_xz, r_0_0, s_0_0, t_xyzz_yyyy, \
                                     t_xyzz_yyyz, t_xyzz_yyzz, t_xyzz_yzzz, t_xyzz_zzzz: VLX_ALIGN)
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

                t_xyzz_yyyy[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 1.5 * pa2pb_x_y[j] * fl3_fx + 0.75 * pa_xyzz[j] * fl2_fx + 3.0 * pa2pb_xzz_y[j] * fl2_fx + 1.5 * pa2pb_xy_yy[j] * fl2_fx + pa2pb_x_yyy[j] * fl2_fx + 3.0 * pa2pb_xyzz_yy[j] * fl1_fx + 2.0 * pa2pb_xzz_yyy[j] * fl1_fx + 0.5 * pa2pb_xy_yyyy[j] * fl1_fx + pa2pb_xyzz_yyyy[j]);

                t_xyzz_yyyy[j] += fl_r_0_0 * (-1.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xyzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xy[j] * fl1_fz * fl3_fx - 6.0 * pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 9.0 * pa_xyzz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_xzz_y[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xyzz_yy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xy_yy[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xyzz_yy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xzz_yyy[j] * fl1_fx * fl1_fz - pa2pb_xy_yyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_yyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_yyyy[j] * fl1_fz);

                t_xyzz_yyyz[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl3_fx + 0.375 * pa2pb_x_z[j] * fl3_fx + 1.5 * pa2pb_xyz_y[j] * fl2_fx + 0.75 * pa2pb_xzz_z[j] * fl2_fx + 1.5 * pa2pb_xz_yy[j] * fl2_fx + 0.75 * pa2pb_xy_yz[j] * fl2_fx + 0.75 * pa2pb_x_yyz[j] * fl2_fx + 1.5 * pa2pb_xyzz_yz[j] * fl1_fx + pa2pb_xyz_yyy[j] * fl1_fx + 1.5 * pa2pb_xzz_yyz[j] * fl1_fx + 0.5 * pa2pb_xy_yyyz[j] * fl1_fx + pa2pb_xyzz_yyyz[j]);

                t_xyzz_yyyz[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_xz[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_xyz_y[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xzz_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyzz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xy_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyzz_yz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyz_yyy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xzz_yyz[j] * fl1_fx * fl1_fz - pa2pb_xy_yyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_yyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_yyyz[j] * fl1_fz);

                t_xyzz_yyzz[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 0.75 * pa2pb_x_y[j] * fl3_fx + 0.25 * pa_xyzz[j] * fl2_fx + pa2pb_xyz_z[j] * fl2_fx + 0.75 * pa2pb_xy_yy[j] * fl2_fx + 0.5 * pa2pb_xzz_y[j] * fl2_fx + 2.0 * pa2pb_xz_yz[j] * fl2_fx + 0.25 * pa2pb_xy_zz[j] * fl2_fx + 0.5 * pa2pb_x_yzz[j] * fl2_fx + 0.5 * pa2pb_xyzz_yy[j] * fl1_fx + 0.5 * pa2pb_xyzz_zz[j] * fl1_fx + 2.0 * pa2pb_xyz_yyz[j] * fl1_fx + pa2pb_xzz_yzz[j] * fl1_fx + 0.5 * pa2pb_xy_yyzz[j] * fl1_fx + pa2pb_xyzz_yyzz[j]);

                t_xyzz_yyzz[j] += fl_r_0_0 * (-pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_xy[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 0.25 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_xyzz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_xyz_z[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xy_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzz_y[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xyzz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_xy_zz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xyzz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xyzz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_xyz_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzz_yzz[j] * fl1_fx * fl1_fz - pa2pb_xy_yyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_yyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_yyzz[j] * fl1_fz);

                t_xyzz_yzzz[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl3_fx + 1.125 * pa2pb_x_z[j] * fl3_fx + 1.5 * pa2pb_xyz_y[j] * fl2_fx + 2.25 * pa2pb_xy_yz[j] * fl2_fx + 0.75 * pa2pb_xzz_z[j] * fl2_fx + 1.5 * pa2pb_xz_zz[j] * fl2_fx + 0.25 * pa2pb_x_zzz[j] * fl2_fx + 1.5 * pa2pb_xyzz_yz[j] * fl1_fx + 3.0 * pa2pb_xyz_yzz[j] * fl1_fx + 0.5 * pa2pb_xzz_zzz[j] * fl1_fx + 0.5 * pa2pb_xy_yzzz[j] * fl1_fx + pa2pb_xyzz_yzzz[j]);

                t_xyzz_yzzz[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * pa_xz[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyz_y[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xy_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_xzz_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyzz_yz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xyzz_yz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xyz_yzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xzz_zzz[j] * fl1_fx * fl1_fz - pa2pb_xy_yzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_yzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_yzzz[j] * fl1_fz);

                t_xyzz_zzzz[j] = fl_s_0_0 * (1.875 * pa_xy[j] * fl3_fx + 0.75 * pa_xyzz[j] * fl2_fx + 6.0 * pa2pb_xyz_z[j] * fl2_fx + 4.5 * pa2pb_xy_zz[j] * fl2_fx + 3.0 * pa2pb_xyzz_zz[j] * fl1_fx + 4.0 * pa2pb_xyz_zzz[j] * fl1_fx + 0.5 * pa2pb_xy_zzzz[j] * fl1_fx + pa2pb_xyzz_zzzz[j]);

                t_xyzz_zzzz[j] += fl_r_0_0 * (-4.5 * pa_xy[j] * fl2_fx * fl1_fz * fl1_fgb + 18.75 * pa_xy[j] * fl3_fx * fl1_fz - 0.75 * pa_xy[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xyzz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_xyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_xyzz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_xyz_z[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_xy_zz[j] * fl2_fx * fl1_fz - 3.0 * pa2pb_xy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xyzz_zz[j] * fl1_fz * fl1_fgb + 42.0 * pa2pb_xyzz_zz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_xyz_zzz[j] * fl1_fz * fl1_fx - pa2pb_xy_zzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_xy_zzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xyzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_135_140(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_xxx = pa2pbDistances.data(1156 * idx + 9);

            auto pa2pb_x_xxy = pa2pbDistances.data(1156 * idx + 10);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_xxxx = pa2pbDistances.data(1156 * idx + 189);

            auto pa2pb_xz_xxxy = pa2pbDistances.data(1156 * idx + 190);

            auto pa2pb_xz_xxxz = pa2pbDistances.data(1156 * idx + 191);

            auto pa2pb_xz_xxyy = pa2pbDistances.data(1156 * idx + 192);

            auto pa2pb_xz_xxyz = pa2pbDistances.data(1156 * idx + 193);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_xzz_x = pa2pbDistances.data(1156 * idx + 476);

            auto pa2pb_xzz_y = pa2pbDistances.data(1156 * idx + 477);

            auto pa2pb_xzz_xxx = pa2pbDistances.data(1156 * idx + 485);

            auto pa2pb_xzz_xxy = pa2pbDistances.data(1156 * idx + 486);

            auto pa2pb_zzz_x = pa2pbDistances.data(1156 * idx + 612);

            auto pa2pb_zzz_y = pa2pbDistances.data(1156 * idx + 613);

            auto pa2pb_zzz_z = pa2pbDistances.data(1156 * idx + 614);

            auto pa2pb_zzz_xxx = pa2pbDistances.data(1156 * idx + 621);

            auto pa2pb_zzz_xxy = pa2pbDistances.data(1156 * idx + 622);

            auto pa2pb_zzz_xxz = pa2pbDistances.data(1156 * idx + 623);

            auto pa2pb_zzz_xyy = pa2pbDistances.data(1156 * idx + 624);

            auto pa2pb_zzz_xyz = pa2pbDistances.data(1156 * idx + 625);

            auto pa2pb_xzzz_xx = pa2pbDistances.data(1156 * idx + 955);

            auto pa2pb_xzzz_xy = pa2pbDistances.data(1156 * idx + 956);

            auto pa2pb_xzzz_xz = pa2pbDistances.data(1156 * idx + 957);

            auto pa2pb_xzzz_yy = pa2pbDistances.data(1156 * idx + 958);

            auto pa2pb_xzzz_yz = pa2pbDistances.data(1156 * idx + 959);

            auto pa2pb_xzzz_xxxx = pa2pbDistances.data(1156 * idx + 971);

            auto pa2pb_xzzz_xxxy = pa2pbDistances.data(1156 * idx + 972);

            auto pa2pb_xzzz_xxxz = pa2pbDistances.data(1156 * idx + 973);

            auto pa2pb_xzzz_xxyy = pa2pbDistances.data(1156 * idx + 974);

            auto pa2pb_xzzz_xxyz = pa2pbDistances.data(1156 * idx + 975);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzzz_xxxx = primBuffer.data(225 * idx + 135);

            auto t_xzzz_xxxy = primBuffer.data(225 * idx + 136);

            auto t_xzzz_xxxz = primBuffer.data(225 * idx + 137);

            auto t_xzzz_xxyy = primBuffer.data(225 * idx + 138);

            auto t_xzzz_xxyz = primBuffer.data(225 * idx + 139);

            // Batch of Integrals (135,140)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxx, pa2pb_x_xxy, pa2pb_x_y, \
                                     pa2pb_xz_xx, pa2pb_xz_xxxx, pa2pb_xz_xxxy, pa2pb_xz_xxxz, pa2pb_xz_xxyy, \
                                     pa2pb_xz_xxyz, pa2pb_xz_xy, pa2pb_xz_xz, pa2pb_xz_yy, pa2pb_xz_yz, pa2pb_xzz_x, \
                                     pa2pb_xzz_xxx, pa2pb_xzz_xxy, pa2pb_xzz_y, pa2pb_xzzz_xx, pa2pb_xzzz_xxxx, \
                                     pa2pb_xzzz_xxxy, pa2pb_xzzz_xxxz, pa2pb_xzzz_xxyy, pa2pb_xzzz_xxyz, pa2pb_xzzz_xy, \
                                     pa2pb_xzzz_xz, pa2pb_xzzz_yy, pa2pb_xzzz_yz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, \
                                     pa2pb_z_xxz, pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_y, pa2pb_z_z, pa2pb_zz_xx, \
                                     pa2pb_zz_xy, pa2pb_zzz_x, pa2pb_zzz_xxx, pa2pb_zzz_xxy, pa2pb_zzz_xxz, \
                                     pa2pb_zzz_xyy, pa2pb_zzz_xyz, pa2pb_zzz_y, pa2pb_zzz_z, pa_xz, pa_xzzz, pa_zz, pb_xx, \
                                     pb_xy, r_0_0, s_0_0, t_xzzz_xxxx, t_xzzz_xxxy, t_xzzz_xxxz, t_xzzz_xxyy, \
                                     t_xzzz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xzzz_xxxx[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 4.5 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa_xzzz[j] * fl2_fx + 3.0 * pa2pb_zzz_x[j] * fl2_fx + 4.5 * pa2pb_xz_xx[j] * fl2_fx + 3.0 * pa2pb_z_xxx[j] * fl2_fx + 3.0 * pa2pb_xzzz_xx[j] * fl1_fx + 2.0 * pa2pb_zzz_xxx[j] * fl1_fx + 1.5 * pa2pb_xz_xxxx[j] * fl1_fx + pa2pb_xzzz_xxxx[j]);

                t_xzzz_xxxx[j] += fl_r_0_0 * (-4.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_xz[j] * fl1_fz * fl3_fx - 6.0 * pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 9.0 * pa_xzzz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_zzz_x[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xzzz_xx[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xz_xx[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xzzz_xx[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_zzz_xxx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxxx[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxxx[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xxxx[j] * fl1_fz);

                t_xzzz_xxxy[j] = fl_s_0_0 * (1.125 * pa2pb_z_y[j] * fl3_fx + 0.75 * pa2pb_zzz_y[j] * fl2_fx + 2.25 * pa2pb_xz_xy[j] * fl2_fx + 2.25 * pa2pb_z_xxy[j] * fl2_fx + 1.5 * pa2pb_xzzz_xy[j] * fl1_fx + 1.5 * pa2pb_zzz_xxy[j] * fl1_fx + 1.5 * pa2pb_xz_xxxy[j] * fl1_fx + pa2pb_xzzz_xxxy[j]);

                t_xzzz_xxxy[j] += fl_r_0_0 * (-2.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_zzz_y[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzzz_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xz_xy[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xzzz_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_zzz_xxy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxxy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxxy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xxxy[j] * fl1_fz);

                t_xzzz_xxxz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_zz[j] * fl3_fx + 1.125 * pa2pb_x_x[j] * fl3_fx + 1.125 * pa2pb_z_z[j] * fl3_fx + 1.125 * pb_xx[j] * fl3_fx + 2.25 * pa2pb_xzz_x[j] * fl2_fx + 0.75 * pa2pb_zzz_z[j] * fl2_fx + 2.25 * pa2pb_zz_xx[j] * fl2_fx + 2.25 * pa2pb_xz_xz[j] * fl2_fx + 0.75 * pa2pb_x_xxx[j] * fl2_fx + 2.25 * pa2pb_z_xxz[j] * fl2_fx + 1.5 * pa2pb_xzzz_xz[j] * fl1_fx + 1.5 * pa2pb_xzz_xxx[j] * fl1_fx + 1.5 * pa2pb_zzz_xxz[j] * fl1_fx + 1.5 * pa2pb_xz_xxxz[j] * fl1_fx + pa2pb_xzzz_xxxz[j]);

                t_xzzz_xxxz[j] += fl_r_0_0 * (-1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 4.5 * fl4_fx * fl1_fz + 11.25 * pa_zz[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 11.25 * pb_xx[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xzz_x[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_zzz_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_z_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzzz_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xz_xz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_xxx[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xzzz_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xzz_xxx[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_zzz_xxz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxxz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxxz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xxxz[j] * fl1_fz);

                t_xzzz_xxyy[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 0.75 * pa2pb_z_x[j] * fl3_fx + 0.25 * pa_xzzz[j] * fl2_fx + 0.5 * pa2pb_zzz_x[j] * fl2_fx + 0.75 * pa2pb_xz_xx[j] * fl2_fx + 0.75 * pa2pb_xz_yy[j] * fl2_fx + 1.5 * pa2pb_z_xyy[j] * fl2_fx + 0.5 * pa2pb_xzzz_xx[j] * fl1_fx + 0.5 * pa2pb_xzzz_yy[j] * fl1_fx + pa2pb_zzz_xyy[j] * fl1_fx + 1.5 * pa2pb_xz_xxyy[j] * fl1_fx + pa2pb_xzzz_xxyy[j]);

                t_xzzz_xxyy[j] += fl_r_0_0 * (-1.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_xz[j] * fl1_fz * fl3_fx - pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 3.0 * pa_xzzz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_zzz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xzzz_yy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xx[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_xz_yy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xzzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xzzz_yy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xxyy[j] * fl1_fz);

                t_xzzz_xxyz[j] = fl_s_0_0 * (0.375 * pa2pb_x_y[j] * fl3_fx + 0.75 * pb_xy[j] * fl3_fx + 0.75 * pa2pb_xzz_y[j] * fl2_fx + 1.5 * pa2pb_zz_xy[j] * fl2_fx + 0.75 * pa2pb_xz_yz[j] * fl2_fx + 0.75 * pa2pb_x_xxy[j] * fl2_fx + 1.5 * pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_xzzz_yz[j] * fl1_fx + 1.5 * pa2pb_xzz_xxy[j] * fl1_fx + pa2pb_zzz_xyz[j] * fl1_fx + 1.5 * pa2pb_xz_xxyz[j] * fl1_fx + pa2pb_xzzz_xxyz[j]);

                t_xzzz_xxyz[j] += fl_r_0_0 * (-0.75 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 7.5 * pb_xy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xzz_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_z_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_xxy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xzzz_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xzz_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_140_145(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_x = pa2pbDistances.data(1156 * idx);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_xxz = pa2pbDistances.data(1156 * idx + 11);

            auto pa2pb_x_xyy = pa2pbDistances.data(1156 * idx + 12);

            auto pa2pb_x_xyz = pa2pbDistances.data(1156 * idx + 13);

            auto pa2pb_x_xzz = pa2pbDistances.data(1156 * idx + 14);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_xz_xx = pa2pbDistances.data(1156 * idx + 173);

            auto pa2pb_xz_xy = pa2pbDistances.data(1156 * idx + 174);

            auto pa2pb_xz_xz = pa2pbDistances.data(1156 * idx + 175);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_xz_xxzz = pa2pbDistances.data(1156 * idx + 194);

            auto pa2pb_xz_xyyy = pa2pbDistances.data(1156 * idx + 195);

            auto pa2pb_xz_xyyz = pa2pbDistances.data(1156 * idx + 196);

            auto pa2pb_xz_xyzz = pa2pbDistances.data(1156 * idx + 197);

            auto pa2pb_xz_xzzz = pa2pbDistances.data(1156 * idx + 198);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_xzz_x = pa2pbDistances.data(1156 * idx + 476);

            auto pa2pb_xzz_z = pa2pbDistances.data(1156 * idx + 478);

            auto pa2pb_xzz_xxz = pa2pbDistances.data(1156 * idx + 487);

            auto pa2pb_xzz_xyy = pa2pbDistances.data(1156 * idx + 488);

            auto pa2pb_xzz_xyz = pa2pbDistances.data(1156 * idx + 489);

            auto pa2pb_xzz_xzz = pa2pbDistances.data(1156 * idx + 490);

            auto pa2pb_zzz_x = pa2pbDistances.data(1156 * idx + 612);

            auto pa2pb_zzz_y = pa2pbDistances.data(1156 * idx + 613);

            auto pa2pb_zzz_z = pa2pbDistances.data(1156 * idx + 614);

            auto pa2pb_zzz_xzz = pa2pbDistances.data(1156 * idx + 626);

            auto pa2pb_zzz_yyy = pa2pbDistances.data(1156 * idx + 627);

            auto pa2pb_zzz_yyz = pa2pbDistances.data(1156 * idx + 628);

            auto pa2pb_zzz_yzz = pa2pbDistances.data(1156 * idx + 629);

            auto pa2pb_zzz_zzz = pa2pbDistances.data(1156 * idx + 630);

            auto pa2pb_xzzz_xx = pa2pbDistances.data(1156 * idx + 955);

            auto pa2pb_xzzz_xy = pa2pbDistances.data(1156 * idx + 956);

            auto pa2pb_xzzz_xz = pa2pbDistances.data(1156 * idx + 957);

            auto pa2pb_xzzz_zz = pa2pbDistances.data(1156 * idx + 960);

            auto pa2pb_xzzz_xxzz = pa2pbDistances.data(1156 * idx + 976);

            auto pa2pb_xzzz_xyyy = pa2pbDistances.data(1156 * idx + 977);

            auto pa2pb_xzzz_xyyz = pa2pbDistances.data(1156 * idx + 978);

            auto pa2pb_xzzz_xyzz = pa2pbDistances.data(1156 * idx + 979);

            auto pa2pb_xzzz_xzzz = pa2pbDistances.data(1156 * idx + 980);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzzz_xxzz = primBuffer.data(225 * idx + 140);

            auto t_xzzz_xyyy = primBuffer.data(225 * idx + 141);

            auto t_xzzz_xyyz = primBuffer.data(225 * idx + 142);

            auto t_xzzz_xyzz = primBuffer.data(225 * idx + 143);

            auto t_xzzz_xzzz = primBuffer.data(225 * idx + 144);

            // Batch of Integrals (140,145)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_x, pa2pb_x_xxz, pa2pb_x_xyy, pa2pb_x_xyz, \
                                     pa2pb_x_xzz, pa2pb_x_z, pa2pb_xz_xx, pa2pb_xz_xxzz, pa2pb_xz_xy, pa2pb_xz_xyyy, \
                                     pa2pb_xz_xyyz, pa2pb_xz_xyzz, pa2pb_xz_xz, pa2pb_xz_xzzz, pa2pb_xz_zz, pa2pb_xzz_x, \
                                     pa2pb_xzz_xxz, pa2pb_xzz_xyy, pa2pb_xzz_xyz, pa2pb_xzz_xzz, pa2pb_xzz_z, \
                                     pa2pb_xzzz_xx, pa2pb_xzzz_xxzz, pa2pb_xzzz_xy, pa2pb_xzzz_xyyy, pa2pb_xzzz_xyyz, \
                                     pa2pb_xzzz_xyzz, pa2pb_xzzz_xz, pa2pb_xzzz_xzzz, pa2pb_xzzz_zz, pa2pb_z_x, \
                                     pa2pb_z_xzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, pa2pb_z_z, \
                                     pa2pb_z_zzz, pa2pb_zz_xz, pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zz_zz, pa2pb_zzz_x, \
                                     pa2pb_zzz_xzz, pa2pb_zzz_y, pa2pb_zzz_yyy, pa2pb_zzz_yyz, pa2pb_zzz_yzz, \
                                     pa2pb_zzz_z, pa2pb_zzz_zzz, pa_xz, pa_xzzz, pa_zz, pb_xz, pb_yy, pb_yz, pb_zz, r_0_0, s_0_0, \
                                     t_xzzz_xxzz, t_xzzz_xyyy, t_xzzz_xyyz, t_xzzz_xyzz, t_xzzz_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xzzz_xxzz[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 2.25 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa2pb_x_z[j] * fl3_fx + 1.5 * pb_xz[j] * fl3_fx + 0.25 * pa_xzzz[j] * fl2_fx + 1.5 * pa2pb_xzz_z[j] * fl2_fx + 2.25 * pa2pb_xz_xx[j] * fl2_fx + 0.5 * pa2pb_zzz_x[j] * fl2_fx + 3.0 * pa2pb_zz_xz[j] * fl2_fx + 0.75 * pa2pb_xz_zz[j] * fl2_fx + 1.5 * pa2pb_x_xxz[j] * fl2_fx + 1.5 * pa2pb_z_xzz[j] * fl2_fx + 0.5 * pa2pb_xzzz_xx[j] * fl1_fx + 0.5 * pa2pb_xzzz_zz[j] * fl1_fx + 3.0 * pa2pb_xzz_xxz[j] * fl1_fx + pa2pb_zzz_xzz[j] * fl1_fx + 1.5 * pa2pb_xz_xxzz[j] * fl1_fx + pa2pb_xzzz_xxzz[j]);

                t_xzzz_xxzz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xz[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - pa_xzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz - pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pb_xz[j] * fl3_fx * fl1_fz + 3.0 * pa_xzzz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xzz_z[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xz_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzz_x[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_z_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_xzzz_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_zz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_xxz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xzzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xzzz_zz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xzz_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xxzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xxzz[j] * fl1_fz);

                t_xzzz_xyyy[j] = fl_s_0_0 * (1.125 * pa2pb_z_y[j] * fl3_fx + 0.75 * pa2pb_zzz_y[j] * fl2_fx + 2.25 * pa2pb_xz_xy[j] * fl2_fx + 0.75 * pa2pb_z_yyy[j] * fl2_fx + 1.5 * pa2pb_xzzz_xy[j] * fl1_fx + 0.5 * pa2pb_zzz_yyy[j] * fl1_fx + 1.5 * pa2pb_xz_xyyy[j] * fl1_fx + pa2pb_xzzz_xyyy[j]);

                t_xzzz_xyyy[j] += fl_r_0_0 * (-2.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_zzz_y[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzzz_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xz_xy[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xzzz_xy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzz_yyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xyyy[j] * fl1_fz);

                t_xzzz_xyyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_zz[j] * fl3_fx + 0.375 * pa2pb_x_x[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 0.75 * pa2pb_xzz_x[j] * fl2_fx + 0.25 * pa2pb_zzz_z[j] * fl2_fx + 0.75 * pa2pb_zz_yy[j] * fl2_fx + 0.75 * pa2pb_xz_xz[j] * fl2_fx + 0.75 * pa2pb_x_xyy[j] * fl2_fx + 0.75 * pa2pb_z_yyz[j] * fl2_fx + 0.5 * pa2pb_xzzz_xz[j] * fl1_fx + 1.5 * pa2pb_xzz_xyy[j] * fl1_fx + 0.5 * pa2pb_zzz_yyz[j] * fl1_fx + 1.5 * pa2pb_xz_xyyz[j] * fl1_fx + pa2pb_xzzz_xyyz[j]);

                t_xzzz_xyyz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_zz[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_xzz_x[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_zzz_z[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_xz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_xyy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xzzz_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xzz_xyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzz_yyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xyyz[j] * fl1_fz);

                t_xzzz_xyzz[j] = fl_s_0_0 * (1.125 * pa2pb_z_y[j] * fl3_fx + 0.75 * pb_yz[j] * fl3_fx + 2.25 * pa2pb_xz_xy[j] * fl2_fx + 0.25 * pa2pb_zzz_y[j] * fl2_fx + 1.5 * pa2pb_zz_yz[j] * fl2_fx + 1.5 * pa2pb_x_xyz[j] * fl2_fx + 0.75 * pa2pb_z_yzz[j] * fl2_fx + 0.5 * pa2pb_xzzz_xy[j] * fl1_fx + 3.0 * pa2pb_xzz_xyz[j] * fl1_fx + 0.5 * pa2pb_zzz_yzz[j] * fl1_fx + 1.5 * pa2pb_xz_xyzz[j] * fl1_fx + pa2pb_xzzz_xyzz[j]);

                t_xzzz_xyzz[j] += fl_r_0_0 * (11.25 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_yz[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xz_xy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zzz_y[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_x_xyz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xzzz_xy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xzz_xyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzz_yzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xyzz[j] * fl1_fz);

                t_xzzz_xzzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa2pb_x_x[j] * fl3_fx + 1.125 * pa_zz[j] * fl3_fx + 3.375 * pa2pb_z_z[j] * fl3_fx + 1.125 * pb_zz[j] * fl3_fx + 2.25 * pa2pb_xzz_x[j] * fl2_fx + 6.75 * pa2pb_xz_xz[j] * fl2_fx + 0.75 * pa2pb_zzz_z[j] * fl2_fx + 2.25 * pa2pb_zz_zz[j] * fl2_fx + 2.25 * pa2pb_x_xzz[j] * fl2_fx + 0.75 * pa2pb_z_zzz[j] * fl2_fx + 1.5 * pa2pb_xzzz_xz[j] * fl1_fx + 4.5 * pa2pb_xzz_xzz[j] * fl1_fx + 0.5 * pa2pb_zzz_zzz[j] * fl1_fx + 1.5 * pa2pb_xz_xzzz[j] * fl1_fx + pa2pb_xzzz_xzzz[j]);

                t_xzzz_xzzz[j] += fl_r_0_0 * (7.5 * fl4_fx * fl1_fz - 1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 18.75 * pa2pb_x_x[j] * fl3_fx * fl1_fz + 11.25 * pa_zz[j] * fl3_fx * fl1_fz + 33.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_x[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xzz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pb_zz[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xzz_x[j] * fl1_fz * fl2_fx + 81.0 * pa2pb_xz_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zzz_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzzz_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_x_xzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xzzz_xz[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_xzz_xzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzz_zzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_xzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_145_150(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_y = pa2pbDistances.data(1156 * idx + 1);

            auto pa2pb_x_z = pa2pbDistances.data(1156 * idx + 2);

            auto pa2pb_x_yyy = pa2pbDistances.data(1156 * idx + 15);

            auto pa2pb_x_yyz = pa2pbDistances.data(1156 * idx + 16);

            auto pa2pb_x_yzz = pa2pbDistances.data(1156 * idx + 17);

            auto pa2pb_x_zzz = pa2pbDistances.data(1156 * idx + 18);

            auto pa2pb_xz_yy = pa2pbDistances.data(1156 * idx + 176);

            auto pa2pb_xz_yz = pa2pbDistances.data(1156 * idx + 177);

            auto pa2pb_xz_zz = pa2pbDistances.data(1156 * idx + 178);

            auto pa2pb_xz_yyyy = pa2pbDistances.data(1156 * idx + 199);

            auto pa2pb_xz_yyyz = pa2pbDistances.data(1156 * idx + 200);

            auto pa2pb_xz_yyzz = pa2pbDistances.data(1156 * idx + 201);

            auto pa2pb_xz_yzzz = pa2pbDistances.data(1156 * idx + 202);

            auto pa2pb_xz_zzzz = pa2pbDistances.data(1156 * idx + 203);

            auto pa2pb_xzz_y = pa2pbDistances.data(1156 * idx + 477);

            auto pa2pb_xzz_z = pa2pbDistances.data(1156 * idx + 478);

            auto pa2pb_xzz_yyy = pa2pbDistances.data(1156 * idx + 491);

            auto pa2pb_xzz_yyz = pa2pbDistances.data(1156 * idx + 492);

            auto pa2pb_xzz_yzz = pa2pbDistances.data(1156 * idx + 493);

            auto pa2pb_xzz_zzz = pa2pbDistances.data(1156 * idx + 494);

            auto pa2pb_xzzz_yy = pa2pbDistances.data(1156 * idx + 958);

            auto pa2pb_xzzz_yz = pa2pbDistances.data(1156 * idx + 959);

            auto pa2pb_xzzz_zz = pa2pbDistances.data(1156 * idx + 960);

            auto pa2pb_xzzz_yyyy = pa2pbDistances.data(1156 * idx + 981);

            auto pa2pb_xzzz_yyyz = pa2pbDistances.data(1156 * idx + 982);

            auto pa2pb_xzzz_yyzz = pa2pbDistances.data(1156 * idx + 983);

            auto pa2pb_xzzz_yzzz = pa2pbDistances.data(1156 * idx + 984);

            auto pa2pb_xzzz_zzzz = pa2pbDistances.data(1156 * idx + 985);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzzz_yyyy = primBuffer.data(225 * idx + 145);

            auto t_xzzz_yyyz = primBuffer.data(225 * idx + 146);

            auto t_xzzz_yyzz = primBuffer.data(225 * idx + 147);

            auto t_xzzz_yzzz = primBuffer.data(225 * idx + 148);

            auto t_xzzz_zzzz = primBuffer.data(225 * idx + 149);

            // Batch of Integrals (145,150)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_y, pa2pb_x_yyy, pa2pb_x_yyz, pa2pb_x_yzz, \
                                     pa2pb_x_z, pa2pb_x_zzz, pa2pb_xz_yy, pa2pb_xz_yyyy, pa2pb_xz_yyyz, \
                                     pa2pb_xz_yyzz, pa2pb_xz_yz, pa2pb_xz_yzzz, pa2pb_xz_zz, pa2pb_xz_zzzz, pa2pb_xzz_y, \
                                     pa2pb_xzz_yyy, pa2pb_xzz_yyz, pa2pb_xzz_yzz, pa2pb_xzz_z, pa2pb_xzz_zzz, \
                                     pa2pb_xzzz_yy, pa2pb_xzzz_yyyy, pa2pb_xzzz_yyyz, pa2pb_xzzz_yyzz, pa2pb_xzzz_yz, \
                                     pa2pb_xzzz_yzzz, pa2pb_xzzz_zz, pa2pb_xzzz_zzzz, pa_xz, pa_xzzz, r_0_0, s_0_0, \
                                     t_xzzz_yyyy, t_xzzz_yyyz, t_xzzz_yyzz, t_xzzz_yzzz, t_xzzz_zzzz: VLX_ALIGN)
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

                t_xzzz_yyyy[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 0.75 * pa_xzzz[j] * fl2_fx + 4.5 * pa2pb_xz_yy[j] * fl2_fx + 3.0 * pa2pb_xzzz_yy[j] * fl1_fx + 1.5 * pa2pb_xz_yyyy[j] * fl1_fx + pa2pb_xzzz_yyyy[j]);

                t_xzzz_yyyy[j] += fl_r_0_0 * (-4.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_xz[j] * fl1_fz * fl3_fx + 9.0 * pa_xzzz[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_xzzz_yy[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_xz_yy[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_xzzz_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_yyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_yyyy[j] * fl1_fz);

                t_xzzz_yyyz[j] = fl_s_0_0 * (1.125 * pa2pb_x_y[j] * fl3_fx + 2.25 * pa2pb_xzz_y[j] * fl2_fx + 2.25 * pa2pb_xz_yz[j] * fl2_fx + 0.75 * pa2pb_x_yyy[j] * fl2_fx + 1.5 * pa2pb_xzzz_yz[j] * fl1_fx + 1.5 * pa2pb_xzz_yyy[j] * fl1_fx + 1.5 * pa2pb_xz_yyyz[j] * fl1_fx + pa2pb_xzzz_yyyz[j]);

                t_xzzz_yyyz[j] += fl_r_0_0 * (-2.25 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_x_y[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_xzz_y[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzzz_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_xz_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_x_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xzzz_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_xzz_yyy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_yyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_yyyz[j] * fl1_fz);

                t_xzzz_yyzz[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 0.75 * pa2pb_x_z[j] * fl3_fx + 0.25 * pa_xzzz[j] * fl2_fx + 1.5 * pa2pb_xzz_z[j] * fl2_fx + 2.25 * pa2pb_xz_yy[j] * fl2_fx + 0.75 * pa2pb_xz_zz[j] * fl2_fx + 1.5 * pa2pb_x_yyz[j] * fl2_fx + 0.5 * pa2pb_xzzz_yy[j] * fl1_fx + 0.5 * pa2pb_xzzz_zz[j] * fl1_fx + 3.0 * pa2pb_xzz_yyz[j] * fl1_fx + 1.5 * pa2pb_xz_yyzz[j] * fl1_fx + pa2pb_xzzz_yyzz[j]);

                t_xzzz_yyzz[j] += fl_r_0_0 * (-3.0 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_xz[j] * fl3_fx * fl1_fz - 0.75 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - pa_xzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_z[j] * fl3_fx * fl1_fz + 3.0 * pa_xzzz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xzz_z[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_xz_yy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_xzzz_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_xz_zz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_x_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_xzzz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_xzzz_zz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_xzz_yyz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_yyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_yyzz[j] * fl1_fz);

                t_xzzz_yzzz[j] = fl_s_0_0 * (1.875 * pa2pb_x_y[j] * fl3_fx + 2.25 * pa2pb_xzz_y[j] * fl2_fx + 6.75 * pa2pb_xz_yz[j] * fl2_fx + 2.25 * pa2pb_x_yzz[j] * fl2_fx + 1.5 * pa2pb_xzzz_yz[j] * fl1_fx + 4.5 * pa2pb_xzz_yzz[j] * fl1_fx + 1.5 * pa2pb_xz_yzzz[j] * fl1_fx + pa2pb_xzzz_yzzz[j]);

                t_xzzz_yzzz[j] += fl_r_0_0 * (18.75 * pa2pb_x_y[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_x_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_xzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 27.0 * pa2pb_xzz_y[j] * fl1_fz * fl2_fx + 81.0 * pa2pb_xz_yz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzzz_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_x_yzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_xzzz_yz[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_xzz_yzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_yzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_yzzz[j] * fl1_fz);

                t_xzzz_zzzz[j] = fl_s_0_0 * (5.625 * pa_xz[j] * fl3_fx + 7.5 * pa2pb_x_z[j] * fl3_fx + 0.75 * pa_xzzz[j] * fl2_fx + 9.0 * pa2pb_xzz_z[j] * fl2_fx + 13.5 * pa2pb_xz_zz[j] * fl2_fx + 3.0 * pa2pb_x_zzz[j] * fl2_fx + 3.0 * pa2pb_xzzz_zz[j] * fl1_fx + 6.0 * pa2pb_xzz_zzz[j] * fl1_fx + 1.5 * pa2pb_xz_zzzz[j] * fl1_fx + pa2pb_xzzz_zzzz[j]);

                t_xzzz_zzzz[j] += fl_r_0_0 * (-13.5 * pa_xz[j] * fl2_fx * fl1_fz * fl1_fgb + 56.25 * pa_xz[j] * fl3_fx * fl1_fz + 75.0 * pa2pb_x_z[j] * fl3_fx * fl1_fz - 2.25 * pa_xz[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_x_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_xzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_xzzz[j] * fl1_fz * fl2_fx + 108.0 * pa2pb_xzz_z[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_xz_zz[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_xz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xzzz_zz[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_x_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_xzzz_zz[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_xzz_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_zzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_xz_zzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_xzzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_150_155(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (150,155)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

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

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_xxxx = pa2pbDistances.data(1156 * idx + 223);

            auto pa2pb_yy_xxxy = pa2pbDistances.data(1156 * idx + 224);

            auto pa2pb_yy_xxxz = pa2pbDistances.data(1156 * idx + 225);

            auto pa2pb_yy_xxyy = pa2pbDistances.data(1156 * idx + 226);

            auto pa2pb_yy_xxyz = pa2pbDistances.data(1156 * idx + 227);

            auto pa2pb_yyy_x = pa2pbDistances.data(1156 * idx + 510);

            auto pa2pb_yyy_y = pa2pbDistances.data(1156 * idx + 511);

            auto pa2pb_yyy_z = pa2pbDistances.data(1156 * idx + 512);

            auto pa2pb_yyy_xxx = pa2pbDistances.data(1156 * idx + 519);

            auto pa2pb_yyy_xxy = pa2pbDistances.data(1156 * idx + 520);

            auto pa2pb_yyy_xxz = pa2pbDistances.data(1156 * idx + 521);

            auto pa2pb_yyyy_xx = pa2pbDistances.data(1156 * idx + 989);

            auto pa2pb_yyyy_xy = pa2pbDistances.data(1156 * idx + 990);

            auto pa2pb_yyyy_xz = pa2pbDistances.data(1156 * idx + 991);

            auto pa2pb_yyyy_yy = pa2pbDistances.data(1156 * idx + 992);

            auto pa2pb_yyyy_yz = pa2pbDistances.data(1156 * idx + 993);

            auto pa2pb_yyyy_xxxx = pa2pbDistances.data(1156 * idx + 1005);

            auto pa2pb_yyyy_xxxy = pa2pbDistances.data(1156 * idx + 1006);

            auto pa2pb_yyyy_xxxz = pa2pbDistances.data(1156 * idx + 1007);

            auto pa2pb_yyyy_xxyy = pa2pbDistances.data(1156 * idx + 1008);

            auto pa2pb_yyyy_xxyz = pa2pbDistances.data(1156 * idx + 1009);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyy_xxxx = primBuffer.data(225 * idx + 150);

            auto t_yyyy_xxxy = primBuffer.data(225 * idx + 151);

            auto t_yyyy_xxxz = primBuffer.data(225 * idx + 152);

            auto t_yyyy_xxyy = primBuffer.data(225 * idx + 153);

            auto t_yyyy_xxyz = primBuffer.data(225 * idx + 154);

            // Batch of Integrals (150,155)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, \
                                     pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xxxx, pa2pb_yy_xxxy, pa2pb_yy_xxxz, \
                                     pa2pb_yy_xxyy, pa2pb_yy_xxyz, pa2pb_yy_xy, pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, \
                                     pa2pb_yyy_x, pa2pb_yyy_xxx, pa2pb_yyy_xxy, pa2pb_yyy_xxz, pa2pb_yyy_y, \
                                     pa2pb_yyy_z, pa2pb_yyyy_xx, pa2pb_yyyy_xxxx, pa2pb_yyyy_xxxy, pa2pb_yyyy_xxxz, \
                                     pa2pb_yyyy_xxyy, pa2pb_yyyy_xxyz, pa2pb_yyyy_xy, pa2pb_yyyy_xz, pa2pb_yyyy_yy, \
                                     pa2pb_yyyy_yz, pa_yy, pa_yyyy, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, pb_xxyz, pb_xy, \
                                     pb_xz, pb_yy, pb_yz, r_0_0, s_0_0, t_yyyy_xxxx, t_yyyy_xxxy, t_yyyy_xxxz, \
                                     t_yyyy_xxyy, t_yyyy_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyyy_xxxx[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_yy[j] * fl3_fx + 0.75 * pa_yyyy[j] * fl2_fx + 2.25 * pb_xx[j] * fl3_fx + 9.0 * pa2pb_yy_xx[j] * fl2_fx + 3.0 * pa2pb_yyyy_xx[j] * fl1_fx + 0.75 * pb_xxxx[j] * fl2_fx + 3.0 * pa2pb_yy_xxxx[j] * fl1_fx + pa2pb_yyyy_xxxx[j]);

                t_yyyy_xxxx[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 9.0 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx + 4.5 * fl4_fx * fl1_fz - 3.0 * pa_yyyy[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_yy[j] * fl1_fz * fl3_fx + 9.0 * pa_yyyy[j] * fl1_fz * fl2_fx - 4.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 18.0 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx + 22.5 * pb_xx[j] * fl3_fx * fl1_fz - 6.0 * pa2pb_yyyy_xx[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_yy_xx[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_yyyy_xx[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxxx[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xxxx[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xxxx[j] * fl1_fz);

                t_yyyy_xxxy[j] = fl_s_0_0 * (4.5 * pa2pb_y_x[j] * fl3_fx + 3.0 * pa2pb_yyy_x[j] * fl2_fx + 1.125 * pb_xy[j] * fl3_fx + 4.5 * pa2pb_yy_xy[j] * fl2_fx + 3.0 * pa2pb_y_xxx[j] * fl2_fx + 1.5 * pa2pb_yyyy_xy[j] * fl1_fx + 2.0 * pa2pb_yyy_xxx[j] * fl1_fx + 0.75 * pb_xxxy[j] * fl2_fx + 3.0 * pa2pb_yy_xxxy[j] * fl1_fx + pa2pb_yyyy_xxxy[j]);

                t_yyyy_xxxy[j] += fl_r_0_0 * (-9.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_yyy_x[j] * fl1_fz * fl2_fx - 2.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_xxx[j] * fl1_fx * fl1_fz * fl1_fga + 11.25 * pb_xy[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_yyyy_xy[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_yy_xy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyy_xy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyy_xxx[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxxy[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xxxy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xxxy[j] * fl1_fz);

                t_yyyy_xxxz[j] = fl_s_0_0 * (1.125 * pb_xz[j] * fl3_fx + 4.5 * pa2pb_yy_xz[j] * fl2_fx + 1.5 * pa2pb_yyyy_xz[j] * fl1_fx + 0.75 * pb_xxxz[j] * fl2_fx + 3.0 * pa2pb_yy_xxxz[j] * fl1_fx + pa2pb_yyyy_xxxz[j]);

                t_yyyy_xxxz[j] += fl_r_0_0 * (-2.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx + 11.25 * pb_xz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_yyyy_xz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_yy_xz[j] * fl1_fz * fl2_fx + 21.0 * pa2pb_yyyy_xz[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxxz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xxxz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xxxz[j] * fl1_fz);

                t_yyyy_xxyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_yy[j] * fl3_fx + 3.0 * pa2pb_y_y[j] * fl3_fx + 1.875 * pb_xx[j] * fl3_fx + 0.25 * pa_yyyy[j] * fl2_fx + 2.0 * pa2pb_yyy_y[j] * fl2_fx + 4.5 * pa2pb_yy_xx[j] * fl2_fx + 0.375 * pb_yy[j] * fl3_fx + 1.5 * pa2pb_yy_yy[j] * fl2_fx + 6.0 * pa2pb_y_xxy[j] * fl2_fx + 0.5 * pa2pb_yyyy_xx[j] * fl1_fx + 0.5 * pa2pb_yyyy_yy[j] * fl1_fx + 4.0 * pa2pb_yyy_xxy[j] * fl1_fx + 0.75 * pb_xxyy[j] * fl2_fx + 3.0 * pa2pb_yy_xxyy[j] * fl1_fx + pa2pb_yyyy_xxyy[j]);

                t_yyyy_xxyy[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 6.0 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 22.5 * pa_yy[j] * fl3_fx * fl1_fz - 1.5 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - pa_yyyy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.0 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 18.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.0 * pa_yyyy[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_yyy_y[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_yy_xx[j] * fl2_fx * fl1_fz - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_y_xxy[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_yy[j] * fl3_fx * fl1_fz - pa2pb_yyyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyyy_yy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yy_yy[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyy_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyyy_yy[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_yyy_xxy[j] * fl1_fz * fl1_fx - 3.0 * pb_xxyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxyy[j] * fl1_fz * fl1_fga + 9.0 * pb_xxyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xxyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xxyy[j] * fl1_fz);

                t_yyyy_xxyz[j] = fl_s_0_0 * (1.5 * pa2pb_y_z[j] * fl3_fx + pa2pb_yyy_z[j] * fl2_fx + 0.375 * pb_yz[j] * fl3_fx + 1.5 * pa2pb_yy_yz[j] * fl2_fx + 3.0 * pa2pb_y_xxz[j] * fl2_fx + 0.5 * pa2pb_yyyy_yz[j] * fl1_fx + 2.0 * pa2pb_yyy_xxz[j] * fl1_fx + 0.75 * pb_xxyz[j] * fl2_fx + 3.0 * pa2pb_yy_xxyz[j] * fl1_fx + pa2pb_yyyy_xxyz[j]);

                t_yyyy_xxyz[j] += fl_r_0_0 * (-3.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 12.0 * pa2pb_yyy_z[j] * fl1_fz * fl2_fx - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_xxz[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_yz[j] * fl3_fx * fl1_fz - pa2pb_yyyy_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yy_yz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyy_yz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyy_xxz[j] * fl1_fz * fl1_fx - 3.0 * pb_xxyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxyz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xxyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_155_160(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (155,160)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_yy_xxzz = pa2pbDistances.data(1156 * idx + 228);

            auto pa2pb_yy_xyyy = pa2pbDistances.data(1156 * idx + 229);

            auto pa2pb_yy_xyyz = pa2pbDistances.data(1156 * idx + 230);

            auto pa2pb_yy_xyzz = pa2pbDistances.data(1156 * idx + 231);

            auto pa2pb_yy_xzzz = pa2pbDistances.data(1156 * idx + 232);

            auto pa2pb_yyy_x = pa2pbDistances.data(1156 * idx + 510);

            auto pa2pb_yyy_xyy = pa2pbDistances.data(1156 * idx + 522);

            auto pa2pb_yyy_xyz = pa2pbDistances.data(1156 * idx + 523);

            auto pa2pb_yyy_xzz = pa2pbDistances.data(1156 * idx + 524);

            auto pa2pb_yyyy_xx = pa2pbDistances.data(1156 * idx + 989);

            auto pa2pb_yyyy_xy = pa2pbDistances.data(1156 * idx + 990);

            auto pa2pb_yyyy_xz = pa2pbDistances.data(1156 * idx + 991);

            auto pa2pb_yyyy_zz = pa2pbDistances.data(1156 * idx + 994);

            auto pa2pb_yyyy_xxzz = pa2pbDistances.data(1156 * idx + 1010);

            auto pa2pb_yyyy_xyyy = pa2pbDistances.data(1156 * idx + 1011);

            auto pa2pb_yyyy_xyyz = pa2pbDistances.data(1156 * idx + 1012);

            auto pa2pb_yyyy_xyzz = pa2pbDistances.data(1156 * idx + 1013);

            auto pa2pb_yyyy_xzzz = pa2pbDistances.data(1156 * idx + 1014);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyy_xxzz = primBuffer.data(225 * idx + 155);

            auto t_yyyy_xyyy = primBuffer.data(225 * idx + 156);

            auto t_yyyy_xyyz = primBuffer.data(225 * idx + 157);

            auto t_yyyy_xyzz = primBuffer.data(225 * idx + 158);

            auto t_yyyy_xzzz = primBuffer.data(225 * idx + 159);

            // Batch of Integrals (155,160)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xyy, pa2pb_y_xyz, pa2pb_y_xzz, \
                                     pa2pb_yy_xx, pa2pb_yy_xxzz, pa2pb_yy_xy, pa2pb_yy_xyyy, pa2pb_yy_xyyz, \
                                     pa2pb_yy_xyzz, pa2pb_yy_xz, pa2pb_yy_xzzz, pa2pb_yy_zz, pa2pb_yyy_x, pa2pb_yyy_xyy, \
                                     pa2pb_yyy_xyz, pa2pb_yyy_xzz, pa2pb_yyyy_xx, pa2pb_yyyy_xxzz, pa2pb_yyyy_xy, \
                                     pa2pb_yyyy_xyyy, pa2pb_yyyy_xyyz, pa2pb_yyyy_xyzz, pa2pb_yyyy_xz, pa2pb_yyyy_xzzz, \
                                     pa2pb_yyyy_zz, pa_yy, pa_yyyy, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, \
                                     pb_xzzz, pb_zz, r_0_0, s_0_0, t_yyyy_xxzz, t_yyyy_xyyy, t_yyyy_xyyz, t_yyyy_xyzz, \
                                     t_yyyy_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyyy_xxzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.75 * pa_yy[j] * fl3_fx + 0.25 * pa_yyyy[j] * fl2_fx + 0.375 * pb_xx[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 1.5 * pa2pb_yy_xx[j] * fl2_fx + 1.5 * pa2pb_yy_zz[j] * fl2_fx + 0.5 * pa2pb_yyyy_xx[j] * fl1_fx + 0.5 * pa2pb_yyyy_zz[j] * fl1_fx + 0.75 * pb_xxzz[j] * fl2_fx + 3.0 * pa2pb_yy_xxzz[j] * fl1_fx + pa2pb_yyyy_xxzz[j]);

                t_yyyy_xxzz[j] += fl_r_0_0 * (-0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - 3.0 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx + 1.5 * fl4_fx * fl1_fz - pa_yyyy[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_yy[j] * fl1_fz * fl3_fx + 3.0 * pa_yyyy[j] * fl1_fz * fl2_fx - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz - pa2pb_yyyy_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyyy_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yy_xx[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yy_zz[j] * fl1_fz * fl2_fx + 7.0 * pa2pb_yyyy_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyyy_zz[j] * fl1_fz * fl1_fx - 3.0 * pb_xxzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xxzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xxzz[j] * fl1_fz);

                t_yyyy_xyyy[j] = fl_s_0_0 * (7.5 * pa2pb_y_x[j] * fl3_fx + 5.625 * pb_xy[j] * fl3_fx + 3.0 * pa2pb_yyy_x[j] * fl2_fx + 13.5 * pa2pb_yy_xy[j] * fl2_fx + 9.0 * pa2pb_y_xyy[j] * fl2_fx + 1.5 * pa2pb_yyyy_xy[j] * fl1_fx + 6.0 * pa2pb_yyy_xyy[j] * fl1_fx + 0.75 * pb_xyyy[j] * fl2_fx + 3.0 * pa2pb_yy_xyyy[j] * fl1_fx + pa2pb_yyyy_xyyy[j]);

                t_yyyy_xyyy[j] += fl_r_0_0 * (75.0 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 9.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 56.25 * pb_xy[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_yyy_x[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_yy_xy[j] * fl2_fx * fl1_fz - 2.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_y_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyyy_xy[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyy_xy[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_yyy_xyy[j] * fl1_fz * fl1_fx - 3.0 * pb_xyyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xyyy[j] * fl1_fz * fl1_fga + 9.0 * pb_xyyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xyyy[j] * fl1_fz);

                t_yyyy_xyyz[j] = fl_s_0_0 * (1.875 * pb_xz[j] * fl3_fx + 4.5 * pa2pb_yy_xz[j] * fl2_fx + 6.0 * pa2pb_y_xyz[j] * fl2_fx + 0.5 * pa2pb_yyyy_xz[j] * fl1_fx + 4.0 * pa2pb_yyy_xyz[j] * fl1_fx + 0.75 * pb_xyyz[j] * fl2_fx + 3.0 * pa2pb_yy_xyyz[j] * fl1_fx + pa2pb_yyyy_xyyz[j]);

                t_yyyy_xyyz[j] += fl_r_0_0 * (-4.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga + 18.75 * pb_xz[j] * fl3_fx * fl1_fz + 54.0 * pa2pb_yy_xz[j] * fl2_fx * fl1_fz - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_y_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyy_xz[j] * fl1_fz * fl1_fgb + 72.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyy_xz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_yyy_xyz[j] * fl1_fz * fl1_fx - 3.0 * pb_xyyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xyyz[j] * fl1_fz * fl1_fga + 9.0 * pb_xyyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xyyz[j] * fl1_fz);

                t_yyyy_xyzz[j] = fl_s_0_0 * (1.5 * pa2pb_y_x[j] * fl3_fx + pa2pb_yyy_x[j] * fl2_fx + 0.375 * pb_xy[j] * fl3_fx + 1.5 * pa2pb_yy_xy[j] * fl2_fx + 3.0 * pa2pb_y_xzz[j] * fl2_fx + 0.5 * pa2pb_yyyy_xy[j] * fl1_fx + 2.0 * pa2pb_yyy_xzz[j] * fl1_fx + 0.75 * pb_xyzz[j] * fl2_fx + 3.0 * pa2pb_yy_xyzz[j] * fl1_fx + pa2pb_yyyy_xyzz[j]);

                t_yyyy_xyzz[j] += fl_r_0_0 * (-3.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 12.0 * pa2pb_yyy_x[j] * fl1_fz * fl2_fx - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_xzz[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_xy[j] * fl3_fx * fl1_fz - pa2pb_yyyy_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yy_xy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyy_xy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyy_xzz[j] * fl1_fz * fl1_fx - 3.0 * pb_xyzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xyzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xyzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xyzz[j] * fl1_fz);

                t_yyyy_xzzz[j] = fl_s_0_0 * (1.125 * pb_xz[j] * fl3_fx + 4.5 * pa2pb_yy_xz[j] * fl2_fx + 1.5 * pa2pb_yyyy_xz[j] * fl1_fx + 0.75 * pb_xzzz[j] * fl2_fx + 3.0 * pa2pb_yy_xzzz[j] * fl1_fx + pa2pb_yyyy_xzzz[j]);

                t_yyyy_xzzz[j] += fl_r_0_0 * (-2.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx + 11.25 * pb_xz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_yyyy_xz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_yy_xz[j] * fl1_fz * fl2_fx + 21.0 * pa2pb_yyyy_xz[j] * fl1_fz * fl1_fx - 3.0 * pb_xzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_xzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_160_165(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (160,165)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_yy_yyyy = pa2pbDistances.data(1156 * idx + 233);

            auto pa2pb_yy_yyyz = pa2pbDistances.data(1156 * idx + 234);

            auto pa2pb_yy_yyzz = pa2pbDistances.data(1156 * idx + 235);

            auto pa2pb_yy_yzzz = pa2pbDistances.data(1156 * idx + 236);

            auto pa2pb_yy_zzzz = pa2pbDistances.data(1156 * idx + 237);

            auto pa2pb_yyy_y = pa2pbDistances.data(1156 * idx + 511);

            auto pa2pb_yyy_z = pa2pbDistances.data(1156 * idx + 512);

            auto pa2pb_yyy_yyy = pa2pbDistances.data(1156 * idx + 525);

            auto pa2pb_yyy_yyz = pa2pbDistances.data(1156 * idx + 526);

            auto pa2pb_yyy_yzz = pa2pbDistances.data(1156 * idx + 527);

            auto pa2pb_yyy_zzz = pa2pbDistances.data(1156 * idx + 528);

            auto pa2pb_yyyy_yy = pa2pbDistances.data(1156 * idx + 992);

            auto pa2pb_yyyy_yz = pa2pbDistances.data(1156 * idx + 993);

            auto pa2pb_yyyy_zz = pa2pbDistances.data(1156 * idx + 994);

            auto pa2pb_yyyy_yyyy = pa2pbDistances.data(1156 * idx + 1015);

            auto pa2pb_yyyy_yyyz = pa2pbDistances.data(1156 * idx + 1016);

            auto pa2pb_yyyy_yyzz = pa2pbDistances.data(1156 * idx + 1017);

            auto pa2pb_yyyy_yzzz = pa2pbDistances.data(1156 * idx + 1018);

            auto pa2pb_yyyy_zzzz = pa2pbDistances.data(1156 * idx + 1019);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyy_yyyy = primBuffer.data(225 * idx + 160);

            auto t_yyyy_yyyz = primBuffer.data(225 * idx + 161);

            auto t_yyyy_yyzz = primBuffer.data(225 * idx + 162);

            auto t_yyyy_yzzz = primBuffer.data(225 * idx + 163);

            auto t_yyyy_zzzz = primBuffer.data(225 * idx + 164);

            // Batch of Integrals (160,165)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, \
                                     pa2pb_y_z, pa2pb_y_zzz, pa2pb_yy_yy, pa2pb_yy_yyyy, pa2pb_yy_yyyz, \
                                     pa2pb_yy_yyzz, pa2pb_yy_yz, pa2pb_yy_yzzz, pa2pb_yy_zz, pa2pb_yy_zzzz, pa2pb_yyy_y, \
                                     pa2pb_yyy_yyy, pa2pb_yyy_yyz, pa2pb_yyy_yzz, pa2pb_yyy_z, pa2pb_yyy_zzz, \
                                     pa2pb_yyyy_yy, pa2pb_yyyy_yyyy, pa2pb_yyyy_yyyz, pa2pb_yyyy_yyzz, pa2pb_yyyy_yz, \
                                     pa2pb_yyyy_yzzz, pa2pb_yyyy_zz, pa2pb_yyyy_zzzz, pa_yy, pa_yyyy, pb_yy, pb_yyyy, pb_yyyz, \
                                     pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, r_0_0, s_0_0, t_yyyy_yyyy, t_yyyy_yyyz, \
                                     t_yyyy_yyzz, t_yyyy_yzzz, t_yyyy_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyyy_yyyy[j] = fl_s_0_0 * (6.5625 * fl4_fx + 11.25 * pa_yy[j] * fl3_fx + 30.0 * pa2pb_y_y[j] * fl3_fx + 11.25 * pb_yy[j] * fl3_fx + 0.75 * pa_yyyy[j] * fl2_fx + 12.0 * pa2pb_yyy_y[j] * fl2_fx + 27.0 * pa2pb_yy_yy[j] * fl2_fx + 12.0 * pa2pb_y_yyy[j] * fl2_fx + 3.0 * pa2pb_yyyy_yy[j] * fl1_fx + 8.0 * pa2pb_yyy_yyy[j] * fl1_fx + 0.75 * pb_yyyy[j] * fl2_fx + 3.0 * pa2pb_yy_yyyy[j] * fl1_fx + pa2pb_yyyy_yyyy[j]);

                t_yyyy_yyyy[j] += fl_r_0_0 * (52.5 * fl4_fx * fl1_fz - 11.25 * fl3_fx * fl1_fz * fl1_fgb - 11.25 * fl3_fx * fl1_fz * fl1_fga - 27.0 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 112.5 * pa_yy[j] * fl3_fx * fl1_fz + 300.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz - 4.5 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 36.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 36.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 27.0 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_yyyy[j] * fl1_fx * fl1_fz * fl1_fgb - 24.0 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 112.5 * pb_yy[j] * fl3_fx * fl1_fz + 9.0 * pa_yyyy[j] * fl1_fz * fl2_fx + 144.0 * pa2pb_yyy_y[j] * fl1_fz * fl2_fx + 324.0 * pa2pb_yy_yy[j] * fl2_fx * fl1_fz - 4.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 24.0 * pa2pb_y_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yyyy_yy[j] * fl1_fz * fl1_fgb + 144.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yyyy_yy[j] * fl1_fz * fl1_fx + 112.0 * pa2pb_yyy_yyy[j] * fl1_fz * fl1_fx - 3.0 * pb_yyyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_yyyy[j] * fl1_fz * fl1_fga + 9.0 * pb_yyyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_yyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_yyyy[j] * fl1_fz);

                t_yyyy_yyyz[j] = fl_s_0_0 * (7.5 * pa2pb_y_z[j] * fl3_fx + 5.625 * pb_yz[j] * fl3_fx + 3.0 * pa2pb_yyy_z[j] * fl2_fx + 13.5 * pa2pb_yy_yz[j] * fl2_fx + 9.0 * pa2pb_y_yyz[j] * fl2_fx + 1.5 * pa2pb_yyyy_yz[j] * fl1_fx + 6.0 * pa2pb_yyy_yyz[j] * fl1_fx + 0.75 * pb_yyyz[j] * fl2_fx + 3.0 * pa2pb_yy_yyyz[j] * fl1_fx + pa2pb_yyyy_yyyz[j]);

                t_yyyy_yyyz[j] += fl_r_0_0 * (75.0 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 9.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 56.25 * pb_yz[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_yyy_z[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_yy_yz[j] * fl2_fx * fl1_fz - 2.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_y_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyyy_yz[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyy_yz[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_yyy_yyz[j] * fl1_fz * fl1_fx - 3.0 * pb_yyyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_yyyz[j] * fl1_fz * fl1_fga + 9.0 * pb_yyyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_yyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_yyyz[j] * fl1_fz);

                t_yyyy_yyzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_yy[j] * fl3_fx + 3.0 * pa2pb_y_y[j] * fl3_fx + 1.875 * pb_zz[j] * fl3_fx + 0.25 * pa_yyyy[j] * fl2_fx + 2.0 * pa2pb_yyy_y[j] * fl2_fx + 4.5 * pa2pb_yy_zz[j] * fl2_fx + 0.375 * pb_yy[j] * fl3_fx + 1.5 * pa2pb_yy_yy[j] * fl2_fx + 6.0 * pa2pb_y_yzz[j] * fl2_fx + 0.5 * pa2pb_yyyy_yy[j] * fl1_fx + 0.5 * pa2pb_yyyy_zz[j] * fl1_fx + 4.0 * pa2pb_yyy_yzz[j] * fl1_fx + 0.75 * pb_yyzz[j] * fl2_fx + 3.0 * pa2pb_yy_yyzz[j] * fl1_fx + pa2pb_yyyy_yyzz[j]);

                t_yyyy_yyzz[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 6.0 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 22.5 * pa_yy[j] * fl3_fx * fl1_fz - 1.5 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - pa_yyyy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.0 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 18.75 * pb_zz[j] * fl3_fx * fl1_fz + 3.0 * pa_yyyy[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_yyy_y[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_yy_zz[j] * fl2_fx * fl1_fz - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_y_yzz[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_yy[j] * fl3_fx * fl1_fz - pa2pb_yyyy_yy[j] * fl1_fz * fl1_fgb - pa2pb_yyyy_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yy_yy[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyy_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyyy_zz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_yyy_yzz[j] * fl1_fz * fl1_fx - 3.0 * pb_yyzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_yyzz[j] * fl1_fz * fl1_fga + 9.0 * pb_yyzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_yyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_yyzz[j] * fl1_fz);

                t_yyyy_yzzz[j] = fl_s_0_0 * (4.5 * pa2pb_y_z[j] * fl3_fx + 3.0 * pa2pb_yyy_z[j] * fl2_fx + 1.125 * pb_yz[j] * fl3_fx + 4.5 * pa2pb_yy_yz[j] * fl2_fx + 3.0 * pa2pb_y_zzz[j] * fl2_fx + 1.5 * pa2pb_yyyy_yz[j] * fl1_fx + 2.0 * pa2pb_yyy_zzz[j] * fl1_fx + 0.75 * pb_yzzz[j] * fl2_fx + 3.0 * pa2pb_yy_yzzz[j] * fl1_fx + pa2pb_yyyy_yzzz[j]);

                t_yyyy_yzzz[j] += fl_r_0_0 * (-9.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_yyy_z[j] * fl1_fz * fl2_fx - 2.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_zzz[j] * fl1_fx * fl1_fz * fl1_fga + 11.25 * pb_yz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_yyyy_yz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_yy_yz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyy_yz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyy_zzz[j] * fl1_fz * fl1_fx - 3.0 * pb_yzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_yzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_yzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_yzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_yzzz[j] * fl1_fz);

                t_yyyy_zzzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_yy[j] * fl3_fx + 0.75 * pa_yyyy[j] * fl2_fx + 2.25 * pb_zz[j] * fl3_fx + 9.0 * pa2pb_yy_zz[j] * fl2_fx + 3.0 * pa2pb_yyyy_zz[j] * fl1_fx + 0.75 * pb_zzzz[j] * fl2_fx + 3.0 * pa2pb_yy_zzzz[j] * fl1_fx + pa2pb_yyyy_zzzz[j]);

                t_yyyy_zzzz[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 9.0 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx + 4.5 * fl4_fx * fl1_fz - 3.0 * pa_yyyy[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_yy[j] * fl1_fz * fl3_fx + 9.0 * pa_yyyy[j] * fl1_fz * fl2_fx - 4.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 18.0 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx + 22.5 * pb_zz[j] * fl3_fx * fl1_fz - 6.0 * pa2pb_yyyy_zz[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_yy_zz[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_yyyy_zz[j] * fl1_fz * fl1_fx - 3.0 * pb_zzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_zzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_zzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yy_zzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yyyy_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_165_170(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (165,170)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_xxxx = pa2pbDistances.data(1156 * idx + 257);

            auto pa2pb_yz_xxxy = pa2pbDistances.data(1156 * idx + 258);

            auto pa2pb_yz_xxxz = pa2pbDistances.data(1156 * idx + 259);

            auto pa2pb_yz_xxyy = pa2pbDistances.data(1156 * idx + 260);

            auto pa2pb_yz_xxyz = pa2pbDistances.data(1156 * idx + 261);

            auto pa2pb_yyy_x = pa2pbDistances.data(1156 * idx + 510);

            auto pa2pb_yyy_y = pa2pbDistances.data(1156 * idx + 511);

            auto pa2pb_yyy_xxx = pa2pbDistances.data(1156 * idx + 519);

            auto pa2pb_yyy_xxy = pa2pbDistances.data(1156 * idx + 520);

            auto pa2pb_yyz_x = pa2pbDistances.data(1156 * idx + 544);

            auto pa2pb_yyz_y = pa2pbDistances.data(1156 * idx + 545);

            auto pa2pb_yyz_z = pa2pbDistances.data(1156 * idx + 546);

            auto pa2pb_yyz_xxx = pa2pbDistances.data(1156 * idx + 553);

            auto pa2pb_yyz_xxy = pa2pbDistances.data(1156 * idx + 554);

            auto pa2pb_yyz_xxz = pa2pbDistances.data(1156 * idx + 555);

            auto pa2pb_yyyz_xx = pa2pbDistances.data(1156 * idx + 1023);

            auto pa2pb_yyyz_xy = pa2pbDistances.data(1156 * idx + 1024);

            auto pa2pb_yyyz_xz = pa2pbDistances.data(1156 * idx + 1025);

            auto pa2pb_yyyz_yy = pa2pbDistances.data(1156 * idx + 1026);

            auto pa2pb_yyyz_yz = pa2pbDistances.data(1156 * idx + 1027);

            auto pa2pb_yyyz_xxxx = pa2pbDistances.data(1156 * idx + 1039);

            auto pa2pb_yyyz_xxxy = pa2pbDistances.data(1156 * idx + 1040);

            auto pa2pb_yyyz_xxxz = pa2pbDistances.data(1156 * idx + 1041);

            auto pa2pb_yyyz_xxyy = pa2pbDistances.data(1156 * idx + 1042);

            auto pa2pb_yyyz_xxyz = pa2pbDistances.data(1156 * idx + 1043);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyz_xxxx = primBuffer.data(225 * idx + 165);

            auto t_yyyz_xxxy = primBuffer.data(225 * idx + 166);

            auto t_yyyz_xxxz = primBuffer.data(225 * idx + 167);

            auto t_yyyz_xxyy = primBuffer.data(225 * idx + 168);

            auto t_yyyz_xxyz = primBuffer.data(225 * idx + 169);

            // Batch of Integrals (165,170)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_y, \
                                     pa2pb_yy_xx, pa2pb_yyy_x, pa2pb_yyy_xxx, pa2pb_yyy_xxy, pa2pb_yyy_y, \
                                     pa2pb_yyyz_xx, pa2pb_yyyz_xxxx, pa2pb_yyyz_xxxy, pa2pb_yyyz_xxxz, pa2pb_yyyz_xxyy, \
                                     pa2pb_yyyz_xxyz, pa2pb_yyyz_xy, pa2pb_yyyz_xz, pa2pb_yyyz_yy, pa2pb_yyyz_yz, \
                                     pa2pb_yyz_x, pa2pb_yyz_xxx, pa2pb_yyz_xxy, pa2pb_yyz_xxz, pa2pb_yyz_y, \
                                     pa2pb_yyz_z, pa2pb_yz_xx, pa2pb_yz_xxxx, pa2pb_yz_xxxy, pa2pb_yz_xxxz, \
                                     pa2pb_yz_xxyy, pa2pb_yz_xxyz, pa2pb_yz_xy, pa2pb_yz_xz, pa2pb_yz_yy, pa2pb_yz_yz, \
                                     pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_xxz, pa2pb_z_y, pa2pb_z_z, pa_yy, \
                                     pa_yyyz, pa_yz, pb_xx, r_0_0, s_0_0, t_yyyz_xxxx, t_yyyz_xxxy, t_yyyz_xxxz, \
                                     t_yyyz_xxyy, t_yyyz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyyz_xxxx[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 0.75 * pa_yyyz[j] * fl2_fx + 4.5 * pa2pb_yz_xx[j] * fl2_fx + 3.0 * pa2pb_yyyz_xx[j] * fl1_fx + 1.5 * pa2pb_yz_xxxx[j] * fl1_fx + pa2pb_yyyz_xxxx[j]);

                t_yyyz_xxxx[j] += fl_r_0_0 * (-4.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyyz[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_yz[j] * fl3_fx * fl1_fz + 9.0 * pa_yyyz[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyyz_xx[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yyyz_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxxx[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxxx[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xxxx[j] * fl1_fz);

                t_yyyz_xxxy[j] = fl_s_0_0 * (1.125 * pa2pb_z_x[j] * fl3_fx + 2.25 * pa2pb_yyz_x[j] * fl2_fx + 2.25 * pa2pb_yz_xy[j] * fl2_fx + 0.75 * pa2pb_z_xxx[j] * fl2_fx + 1.5 * pa2pb_yyyz_xy[j] * fl1_fx + 1.5 * pa2pb_yyz_xxx[j] * fl1_fx + 1.5 * pa2pb_yz_xxxy[j] * fl1_fx + pa2pb_yyyz_xxxy[j]);

                t_yyyz_xxxy[j] += fl_r_0_0 * (-2.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_yyz_x[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyyz_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyz_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyz_xxx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xxxy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxxy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xxxy[j] * fl1_fz);

                t_yyyz_xxxz[j] = fl_s_0_0 * (1.125 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa2pb_yyy_x[j] * fl2_fx + 2.25 * pa2pb_yz_xz[j] * fl2_fx + 0.75 * pa2pb_y_xxx[j] * fl2_fx + 1.5 * pa2pb_yyyz_xz[j] * fl1_fx + 0.5 * pa2pb_yyy_xxx[j] * fl1_fx + 1.5 * pa2pb_yz_xxxz[j] * fl1_fx + pa2pb_yyyz_xxxz[j]);

                t_yyyz_xxxz[j] += fl_r_0_0 * (-2.25 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yyy_x[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyyz_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyz_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyy_xxx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxxz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxxz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xxxz[j] * fl1_fz);

                t_yyyz_xxyy[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 0.75 * pa2pb_z_y[j] * fl3_fx + 0.25 * pa_yyyz[j] * fl2_fx + 1.5 * pa2pb_yyz_y[j] * fl2_fx + 2.25 * pa2pb_yz_xx[j] * fl2_fx + 0.75 * pa2pb_yz_yy[j] * fl2_fx + 1.5 * pa2pb_z_xxy[j] * fl2_fx + 0.5 * pa2pb_yyyz_xx[j] * fl1_fx + 0.5 * pa2pb_yyyz_yy[j] * fl1_fx + 3.0 * pa2pb_yyz_xxy[j] * fl1_fx + 1.5 * pa2pb_yz_xxyy[j] * fl1_fx + pa2pb_yyyz_xxyy[j]);

                t_yyyz_xxyy[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_yz[j] * fl3_fx * fl1_fz - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_yyyz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 3.0 * pa_yyyz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yyz_y[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyyz_yy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyyz_yy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yyz_xxy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xxyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xxyy[j] * fl1_fz);

                t_yyyz_xxyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_xx[j] * fl3_fx + 0.25 * pa2pb_yyy_y[j] * fl2_fx + 0.75 * pa2pb_yyz_z[j] * fl2_fx + 0.75 * pa2pb_yy_xx[j] * fl2_fx + 0.75 * pa2pb_yz_yz[j] * fl2_fx + 0.75 * pa2pb_y_xxy[j] * fl2_fx + 0.75 * pa2pb_z_xxz[j] * fl2_fx + 0.5 * pa2pb_yyyz_yz[j] * fl1_fx + 0.5 * pa2pb_yyy_xxy[j] * fl1_fx + 1.5 * pa2pb_yyz_xxz[j] * fl1_fx + 1.5 * pa2pb_yz_xxyz[j] * fl1_fx + pa2pb_yyyz_xxyz[j]);

                t_yyyz_xxyz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_yy[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_yyy_y[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yyz_z[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yy_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xxy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyz_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyy_xxy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyz_xxz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xxyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_170_175(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (170,175)

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

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_yz_xxzz = pa2pbDistances.data(1156 * idx + 262);

            auto pa2pb_yz_xyyy = pa2pbDistances.data(1156 * idx + 263);

            auto pa2pb_yz_xyyz = pa2pbDistances.data(1156 * idx + 264);

            auto pa2pb_yz_xyzz = pa2pbDistances.data(1156 * idx + 265);

            auto pa2pb_yz_xzzz = pa2pbDistances.data(1156 * idx + 266);

            auto pa2pb_yyy_x = pa2pbDistances.data(1156 * idx + 510);

            auto pa2pb_yyy_z = pa2pbDistances.data(1156 * idx + 512);

            auto pa2pb_yyy_xxz = pa2pbDistances.data(1156 * idx + 521);

            auto pa2pb_yyy_xyy = pa2pbDistances.data(1156 * idx + 522);

            auto pa2pb_yyy_xyz = pa2pbDistances.data(1156 * idx + 523);

            auto pa2pb_yyy_xzz = pa2pbDistances.data(1156 * idx + 524);

            auto pa2pb_yyz_x = pa2pbDistances.data(1156 * idx + 544);

            auto pa2pb_yyz_xyy = pa2pbDistances.data(1156 * idx + 556);

            auto pa2pb_yyz_xyz = pa2pbDistances.data(1156 * idx + 557);

            auto pa2pb_yyz_xzz = pa2pbDistances.data(1156 * idx + 558);

            auto pa2pb_yyyz_xx = pa2pbDistances.data(1156 * idx + 1023);

            auto pa2pb_yyyz_xy = pa2pbDistances.data(1156 * idx + 1024);

            auto pa2pb_yyyz_xz = pa2pbDistances.data(1156 * idx + 1025);

            auto pa2pb_yyyz_zz = pa2pbDistances.data(1156 * idx + 1028);

            auto pa2pb_yyyz_xxzz = pa2pbDistances.data(1156 * idx + 1044);

            auto pa2pb_yyyz_xyyy = pa2pbDistances.data(1156 * idx + 1045);

            auto pa2pb_yyyz_xyyz = pa2pbDistances.data(1156 * idx + 1046);

            auto pa2pb_yyyz_xyzz = pa2pbDistances.data(1156 * idx + 1047);

            auto pa2pb_yyyz_xzzz = pa2pbDistances.data(1156 * idx + 1048);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyz_xxzz = primBuffer.data(225 * idx + 170);

            auto t_yyyz_xyyy = primBuffer.data(225 * idx + 171);

            auto t_yyyz_xyyz = primBuffer.data(225 * idx + 172);

            auto t_yyyz_xyzz = primBuffer.data(225 * idx + 173);

            auto t_yyyz_xzzz = primBuffer.data(225 * idx + 174);

            // Batch of Integrals (170,175)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xxz, pa2pb_y_xyy, pa2pb_y_xyz, \
                                     pa2pb_y_xzz, pa2pb_y_z, pa2pb_yy_xy, pa2pb_yy_xz, pa2pb_yyy_x, pa2pb_yyy_xxz, \
                                     pa2pb_yyy_xyy, pa2pb_yyy_xyz, pa2pb_yyy_xzz, pa2pb_yyy_z, pa2pb_yyyz_xx, \
                                     pa2pb_yyyz_xxzz, pa2pb_yyyz_xy, pa2pb_yyyz_xyyy, pa2pb_yyyz_xyyz, pa2pb_yyyz_xyzz, \
                                     pa2pb_yyyz_xz, pa2pb_yyyz_xzzz, pa2pb_yyyz_zz, pa2pb_yyz_x, pa2pb_yyz_xyy, \
                                     pa2pb_yyz_xyz, pa2pb_yyz_xzz, pa2pb_yz_xx, pa2pb_yz_xxzz, pa2pb_yz_xy, \
                                     pa2pb_yz_xyyy, pa2pb_yz_xyyz, pa2pb_yz_xyzz, pa2pb_yz_xz, pa2pb_yz_xzzz, \
                                     pa2pb_yz_zz, pa2pb_z_x, pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, pa_yyyz, pa_yz, pb_xy, \
                                     pb_xz, r_0_0, s_0_0, t_yyyz_xxzz, t_yyyz_xyyy, t_yyyz_xyyz, t_yyyz_xyzz, \
                                     t_yyyz_xzzz: VLX_ALIGN)
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

                t_yyyz_xxzz[j] = fl_s_0_0 * (0.375 * pa_yz[j] * fl3_fx + 0.75 * pa2pb_y_z[j] * fl3_fx + 0.25 * pa_yyyz[j] * fl2_fx + 0.5 * pa2pb_yyy_z[j] * fl2_fx + 0.75 * pa2pb_yz_xx[j] * fl2_fx + 0.75 * pa2pb_yz_zz[j] * fl2_fx + 1.5 * pa2pb_y_xxz[j] * fl2_fx + 0.5 * pa2pb_yyyz_xx[j] * fl1_fx + 0.5 * pa2pb_yyyz_zz[j] * fl1_fx + pa2pb_yyy_xxz[j] * fl1_fx + 1.5 * pa2pb_yz_xxzz[j] * fl1_fx + pa2pb_yyyz_xxzz[j]);

                t_yyyz_xxzz[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_yz[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 3.0 * pa_yyyz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yyy_z[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyyz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyyz_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyyz_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xxz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xxzz[j] * fl1_fz);

                t_yyyz_xyyy[j] = fl_s_0_0 * (1.875 * pa2pb_z_x[j] * fl3_fx + 2.25 * pa2pb_yyz_x[j] * fl2_fx + 6.75 * pa2pb_yz_xy[j] * fl2_fx + 2.25 * pa2pb_z_xyy[j] * fl2_fx + 1.5 * pa2pb_yyyz_xy[j] * fl1_fx + 4.5 * pa2pb_yyz_xyy[j] * fl1_fx + 1.5 * pa2pb_yz_xyyy[j] * fl1_fx + pa2pb_yyyz_xyyy[j]);

                t_yyyz_xyyy[j] += fl_r_0_0 * (18.75 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 27.0 * pa2pb_yyz_x[j] * fl2_fx * fl1_fz + 81.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyyz_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyz_xy[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_yyz_xyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xyyy[j] * fl1_fz);

                t_yyyz_xyyz[j] = fl_s_0_0 * (1.125 * pa2pb_y_x[j] * fl3_fx + 0.75 * pb_xy[j] * fl3_fx + 0.25 * pa2pb_yyy_x[j] * fl2_fx + 1.5 * pa2pb_yy_xy[j] * fl2_fx + 2.25 * pa2pb_yz_xz[j] * fl2_fx + 0.75 * pa2pb_y_xyy[j] * fl2_fx + 1.5 * pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_yyyz_xz[j] * fl1_fx + 0.5 * pa2pb_yyy_xyy[j] * fl1_fx + 3.0 * pa2pb_yyz_xyz[j] * fl1_fx + 1.5 * pa2pb_yz_xyyz[j] * fl1_fx + pa2pb_yyyz_xyyz[j]);

                t_yyyz_xyyz[j] += fl_r_0_0 * (11.25 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_xy[j] * fl3_fx * fl1_fz + 3.0 * pa2pb_yyy_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yy_xy[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyz_xz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyy_xyy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yyz_xyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xyyz[j] * fl1_fz);

                t_yyyz_xyzz[j] = fl_s_0_0 * (0.375 * pa2pb_z_x[j] * fl3_fx + 0.75 * pb_xz[j] * fl3_fx + 0.75 * pa2pb_yyz_x[j] * fl2_fx + 1.5 * pa2pb_yy_xz[j] * fl2_fx + 0.75 * pa2pb_yz_xy[j] * fl2_fx + 1.5 * pa2pb_y_xyz[j] * fl2_fx + 0.75 * pa2pb_z_xzz[j] * fl2_fx + 0.5 * pa2pb_yyyz_xy[j] * fl1_fx + pa2pb_yyy_xyz[j] * fl1_fx + 1.5 * pa2pb_yyz_xzz[j] * fl1_fx + 1.5 * pa2pb_yz_xyzz[j] * fl1_fx + pa2pb_yyyz_xyzz[j]);

                t_yyyz_xyzz[j] += fl_r_0_0 * (-0.75 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 7.5 * pb_xz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yyz_x[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yy_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xyz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyz_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_xyz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyz_xzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xyzz[j] * fl1_fz);

                t_yyyz_xzzz[j] = fl_s_0_0 * (1.125 * pa2pb_y_x[j] * fl3_fx + 0.75 * pa2pb_yyy_x[j] * fl2_fx + 2.25 * pa2pb_yz_xz[j] * fl2_fx + 2.25 * pa2pb_y_xzz[j] * fl2_fx + 1.5 * pa2pb_yyyz_xz[j] * fl1_fx + 1.5 * pa2pb_yyy_xzz[j] * fl1_fx + 1.5 * pa2pb_yz_xzzz[j] * fl1_fx + pa2pb_yyyz_xzzz[j]);

                t_yyyz_xzzz[j] += fl_r_0_0 * (-2.25 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_yyy_x[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yyy_x[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyyz_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyz_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyy_xzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_175_180(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (175,180)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_yz_yyyy = pa2pbDistances.data(1156 * idx + 267);

            auto pa2pb_yz_yyyz = pa2pbDistances.data(1156 * idx + 268);

            auto pa2pb_yz_yyzz = pa2pbDistances.data(1156 * idx + 269);

            auto pa2pb_yz_yzzz = pa2pbDistances.data(1156 * idx + 270);

            auto pa2pb_yz_zzzz = pa2pbDistances.data(1156 * idx + 271);

            auto pa2pb_yyy_y = pa2pbDistances.data(1156 * idx + 511);

            auto pa2pb_yyy_z = pa2pbDistances.data(1156 * idx + 512);

            auto pa2pb_yyy_yyy = pa2pbDistances.data(1156 * idx + 525);

            auto pa2pb_yyy_yyz = pa2pbDistances.data(1156 * idx + 526);

            auto pa2pb_yyy_yzz = pa2pbDistances.data(1156 * idx + 527);

            auto pa2pb_yyy_zzz = pa2pbDistances.data(1156 * idx + 528);

            auto pa2pb_yyz_y = pa2pbDistances.data(1156 * idx + 545);

            auto pa2pb_yyz_z = pa2pbDistances.data(1156 * idx + 546);

            auto pa2pb_yyz_yyy = pa2pbDistances.data(1156 * idx + 559);

            auto pa2pb_yyz_yyz = pa2pbDistances.data(1156 * idx + 560);

            auto pa2pb_yyz_yzz = pa2pbDistances.data(1156 * idx + 561);

            auto pa2pb_yyz_zzz = pa2pbDistances.data(1156 * idx + 562);

            auto pa2pb_yyyz_yy = pa2pbDistances.data(1156 * idx + 1026);

            auto pa2pb_yyyz_yz = pa2pbDistances.data(1156 * idx + 1027);

            auto pa2pb_yyyz_zz = pa2pbDistances.data(1156 * idx + 1028);

            auto pa2pb_yyyz_yyyy = pa2pbDistances.data(1156 * idx + 1049);

            auto pa2pb_yyyz_yyyz = pa2pbDistances.data(1156 * idx + 1050);

            auto pa2pb_yyyz_yyzz = pa2pbDistances.data(1156 * idx + 1051);

            auto pa2pb_yyyz_yzzz = pa2pbDistances.data(1156 * idx + 1052);

            auto pa2pb_yyyz_zzzz = pa2pbDistances.data(1156 * idx + 1053);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyz_yyyy = primBuffer.data(225 * idx + 175);

            auto t_yyyz_yyyz = primBuffer.data(225 * idx + 176);

            auto t_yyyz_yyzz = primBuffer.data(225 * idx + 177);

            auto t_yyyz_yzzz = primBuffer.data(225 * idx + 178);

            auto t_yyyz_zzzz = primBuffer.data(225 * idx + 179);

            // Batch of Integrals (175,180)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, \
                                     pa2pb_y_z, pa2pb_y_zzz, pa2pb_yy_yy, pa2pb_yy_yz, pa2pb_yy_zz, pa2pb_yyy_y, \
                                     pa2pb_yyy_yyy, pa2pb_yyy_yyz, pa2pb_yyy_yzz, pa2pb_yyy_z, pa2pb_yyy_zzz, \
                                     pa2pb_yyyz_yy, pa2pb_yyyz_yyyy, pa2pb_yyyz_yyyz, pa2pb_yyyz_yyzz, pa2pb_yyyz_yz, \
                                     pa2pb_yyyz_yzzz, pa2pb_yyyz_zz, pa2pb_yyyz_zzzz, pa2pb_yyz_y, pa2pb_yyz_yyy, \
                                     pa2pb_yyz_yyz, pa2pb_yyz_yzz, pa2pb_yyz_z, pa2pb_yyz_zzz, pa2pb_yz_yy, \
                                     pa2pb_yz_yyyy, pa2pb_yz_yyyz, pa2pb_yz_yyzz, pa2pb_yz_yz, pa2pb_yz_yzzz, \
                                     pa2pb_yz_zz, pa2pb_yz_zzzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, \
                                     pa2pb_z_z, pa2pb_z_zzz, pa_yy, pa_yyyz, pa_yz, pb_yy, pb_yz, pb_zz, r_0_0, s_0_0, \
                                     t_yyyz_yyyy, t_yyyz_yyyz, t_yyyz_yyzz, t_yyyz_yzzz, t_yyyz_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyyz_yyyy[j] = fl_s_0_0 * (5.625 * pa_yz[j] * fl3_fx + 7.5 * pa2pb_z_y[j] * fl3_fx + 0.75 * pa_yyyz[j] * fl2_fx + 9.0 * pa2pb_yyz_y[j] * fl2_fx + 13.5 * pa2pb_yz_yy[j] * fl2_fx + 3.0 * pa2pb_z_yyy[j] * fl2_fx + 3.0 * pa2pb_yyyz_yy[j] * fl1_fx + 6.0 * pa2pb_yyz_yyy[j] * fl1_fx + 1.5 * pa2pb_yz_yyyy[j] * fl1_fx + pa2pb_yyyz_yyyy[j]);

                t_yyyz_yyyy[j] += fl_r_0_0 * (-13.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 56.25 * pa_yz[j] * fl3_fx * fl1_fz + 75.0 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 2.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_yyyz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_yyyz[j] * fl1_fz * fl2_fx + 108.0 * pa2pb_yyz_y[j] * fl2_fx * fl1_fz + 162.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yyyz_yy[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yyyz_yy[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_yyz_yyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_yyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_yyyy[j] * fl1_fz);

                t_yyyz_yyyz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.125 * pa_yy[j] * fl3_fx + 3.375 * pa2pb_y_y[j] * fl3_fx + 1.875 * pa2pb_z_z[j] * fl3_fx + 1.125 * pb_yy[j] * fl3_fx + 0.75 * pa2pb_yyy_y[j] * fl2_fx + 2.25 * pa2pb_yyz_z[j] * fl2_fx + 2.25 * pa2pb_yy_yy[j] * fl2_fx + 6.75 * pa2pb_yz_yz[j] * fl2_fx + 0.75 * pa2pb_y_yyy[j] * fl2_fx + 2.25 * pa2pb_z_yyz[j] * fl2_fx + 1.5 * pa2pb_yyyz_yz[j] * fl1_fx + 0.5 * pa2pb_yyy_yyy[j] * fl1_fx + 4.5 * pa2pb_yyz_yyz[j] * fl1_fx + 1.5 * pa2pb_yz_yyyz[j] * fl1_fx + pa2pb_yyyz_yyyz[j]);

                t_yyyz_yyyz[j] += fl_r_0_0 * (7.5 * fl4_fx * fl1_fz - 1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_yy[j] * fl3_fx * fl1_fz + 33.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 18.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pb_yy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yyy_y[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_yyz_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yy_yy[j] * fl2_fx * fl1_fz + 81.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyyz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyz_yz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyy_yyy[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_yyz_yyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_yyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_yyyz[j] * fl1_fz);

                t_yyyz_yyzz[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 2.25 * pa2pb_y_z[j] * fl3_fx + 0.75 * pa2pb_z_y[j] * fl3_fx + 1.5 * pb_yz[j] * fl3_fx + 0.25 * pa_yyyz[j] * fl2_fx + 0.5 * pa2pb_yyy_z[j] * fl2_fx + 1.5 * pa2pb_yyz_y[j] * fl2_fx + 3.0 * pa2pb_yy_yz[j] * fl2_fx + 2.25 * pa2pb_yz_zz[j] * fl2_fx + 0.75 * pa2pb_yz_yy[j] * fl2_fx + 1.5 * pa2pb_y_yyz[j] * fl2_fx + 1.5 * pa2pb_z_yzz[j] * fl2_fx + 0.5 * pa2pb_yyyz_yy[j] * fl1_fx + 0.5 * pa2pb_yyyz_zz[j] * fl1_fx + pa2pb_yyy_yyz[j] * fl1_fx + 3.0 * pa2pb_yyz_yzz[j] * fl1_fx + 1.5 * pa2pb_yz_yyzz[j] * fl1_fx + pa2pb_yyyz_yyzz[j]);

                t_yyyz_yyzz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_yz[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - pa_yyyz[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 15.0 * pb_yz[j] * fl3_fx * fl1_fz + 3.0 * pa_yyyz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yyy_z[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yyz_y[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_yz[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_yyz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_yy[j] * fl1_fz * fl1_fgb - pa2pb_yyyz_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyyz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyyz_zz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyy_yyz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yyz_yzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_yyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_yyzz[j] * fl1_fz);

                t_yyyz_yzzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_yy[j] * fl3_fx + 1.125 * pa2pb_y_y[j] * fl3_fx + 1.125 * pa2pb_z_z[j] * fl3_fx + 1.125 * pb_zz[j] * fl3_fx + 0.75 * pa2pb_yyy_y[j] * fl2_fx + 2.25 * pa2pb_yyz_z[j] * fl2_fx + 2.25 * pa2pb_yy_zz[j] * fl2_fx + 2.25 * pa2pb_yz_yz[j] * fl2_fx + 2.25 * pa2pb_y_yzz[j] * fl2_fx + 0.75 * pa2pb_z_zzz[j] * fl2_fx + 1.5 * pa2pb_yyyz_yz[j] * fl1_fx + 1.5 * pa2pb_yyy_yzz[j] * fl1_fx + 1.5 * pa2pb_yyz_zzz[j] * fl1_fx + 1.5 * pa2pb_yz_yzzz[j] * fl1_fx + pa2pb_yyyz_yzzz[j]);

                t_yyyz_yzzz[j] += fl_r_0_0 * (-1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 4.5 * fl4_fx * fl1_fz + 11.25 * pa_yy[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_y[j] * fl1_fz * fl1_fga * fl2_fx - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yyy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 11.25 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 11.25 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yyy_y[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_yyz_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_yy_zz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_yzz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyyz_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyyz_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyy_yzz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yyz_zzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_yzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_yzzz[j] * fl1_fz);

                t_yyyz_zzzz[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 4.5 * pa2pb_y_z[j] * fl3_fx + 0.75 * pa_yyyz[j] * fl2_fx + 3.0 * pa2pb_yyy_z[j] * fl2_fx + 4.5 * pa2pb_yz_zz[j] * fl2_fx + 3.0 * pa2pb_y_zzz[j] * fl2_fx + 3.0 * pa2pb_yyyz_zz[j] * fl1_fx + 2.0 * pa2pb_yyy_zzz[j] * fl1_fx + 1.5 * pa2pb_yz_zzzz[j] * fl1_fx + pa2pb_yyyz_zzzz[j]);

                t_yyyz_zzzz[j] += fl_r_0_0 * (-4.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_y_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyyz[j] * fl1_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_yyy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_yz[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 9.0 * pa_yyyz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_yyy_z[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyyz_zz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yyyz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyy_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_zzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_zzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyyz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_180_185(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (180,185)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

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

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_xxxx = pa2pbDistances.data(1156 * idx + 223);

            auto pa2pb_yy_xxxy = pa2pbDistances.data(1156 * idx + 224);

            auto pa2pb_yy_xxxz = pa2pbDistances.data(1156 * idx + 225);

            auto pa2pb_yy_xxyy = pa2pbDistances.data(1156 * idx + 226);

            auto pa2pb_yy_xxyz = pa2pbDistances.data(1156 * idx + 227);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_xxxx = pa2pbDistances.data(1156 * idx + 291);

            auto pa2pb_zz_xxxy = pa2pbDistances.data(1156 * idx + 292);

            auto pa2pb_zz_xxxz = pa2pbDistances.data(1156 * idx + 293);

            auto pa2pb_zz_xxyy = pa2pbDistances.data(1156 * idx + 294);

            auto pa2pb_zz_xxyz = pa2pbDistances.data(1156 * idx + 295);

            auto pa2pb_yyz_x = pa2pbDistances.data(1156 * idx + 544);

            auto pa2pb_yyz_y = pa2pbDistances.data(1156 * idx + 545);

            auto pa2pb_yyz_xxx = pa2pbDistances.data(1156 * idx + 553);

            auto pa2pb_yyz_xxy = pa2pbDistances.data(1156 * idx + 554);

            auto pa2pb_yzz_x = pa2pbDistances.data(1156 * idx + 578);

            auto pa2pb_yzz_y = pa2pbDistances.data(1156 * idx + 579);

            auto pa2pb_yzz_z = pa2pbDistances.data(1156 * idx + 580);

            auto pa2pb_yzz_xxx = pa2pbDistances.data(1156 * idx + 587);

            auto pa2pb_yzz_xxy = pa2pbDistances.data(1156 * idx + 588);

            auto pa2pb_yzz_xxz = pa2pbDistances.data(1156 * idx + 589);

            auto pa2pb_yyzz_xx = pa2pbDistances.data(1156 * idx + 1057);

            auto pa2pb_yyzz_xy = pa2pbDistances.data(1156 * idx + 1058);

            auto pa2pb_yyzz_xz = pa2pbDistances.data(1156 * idx + 1059);

            auto pa2pb_yyzz_yy = pa2pbDistances.data(1156 * idx + 1060);

            auto pa2pb_yyzz_yz = pa2pbDistances.data(1156 * idx + 1061);

            auto pa2pb_yyzz_xxxx = pa2pbDistances.data(1156 * idx + 1073);

            auto pa2pb_yyzz_xxxy = pa2pbDistances.data(1156 * idx + 1074);

            auto pa2pb_yyzz_xxxz = pa2pbDistances.data(1156 * idx + 1075);

            auto pa2pb_yyzz_xxyy = pa2pbDistances.data(1156 * idx + 1076);

            auto pa2pb_yyzz_xxyz = pa2pbDistances.data(1156 * idx + 1077);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyzz_xxxx = primBuffer.data(225 * idx + 180);

            auto t_yyzz_xxxy = primBuffer.data(225 * idx + 181);

            auto t_yyzz_xxxz = primBuffer.data(225 * idx + 182);

            auto t_yyzz_xxyy = primBuffer.data(225 * idx + 183);

            auto t_yyzz_xxyz = primBuffer.data(225 * idx + 184);

            // Batch of Integrals (180,185)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_xxz, \
                                     pa2pb_y_y, pa2pb_y_z, pa2pb_yy_xx, pa2pb_yy_xxxx, pa2pb_yy_xxxy, pa2pb_yy_xxxz, \
                                     pa2pb_yy_xxyy, pa2pb_yy_xxyz, pa2pb_yy_xy, pa2pb_yy_xz, pa2pb_yy_yy, pa2pb_yy_yz, \
                                     pa2pb_yyz_x, pa2pb_yyz_xxx, pa2pb_yyz_xxy, pa2pb_yyz_y, pa2pb_yyzz_xx, \
                                     pa2pb_yyzz_xxxx, pa2pb_yyzz_xxxy, pa2pb_yyzz_xxxz, pa2pb_yyzz_xxyy, pa2pb_yyzz_xxyz, \
                                     pa2pb_yyzz_xy, pa2pb_yyzz_xz, pa2pb_yyzz_yy, pa2pb_yyzz_yz, pa2pb_yz_xx, \
                                     pa2pb_yzz_x, pa2pb_yzz_xxx, pa2pb_yzz_xxy, pa2pb_yzz_xxz, pa2pb_yzz_y, \
                                     pa2pb_yzz_z, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_y, pa2pb_zz_xx, \
                                     pa2pb_zz_xxxx, pa2pb_zz_xxxy, pa2pb_zz_xxxz, pa2pb_zz_xxyy, pa2pb_zz_xxyz, \
                                     pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, pa2pb_zz_yz, pa_yy, pa_yyzz, pa_yz, pa_zz, pb_xx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_yy, pb_yz, r_0_0, s_0_0, \
                                     t_yyzz_xxxx, t_yyzz_xxxy, t_yyzz_xxxz, t_yyzz_xxyy, t_yyzz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyzz_xxxx[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.375 * pa_zz[j] * fl3_fx + 0.75 * pa_yyzz[j] * fl2_fx + 0.75 * pb_xx[j] * fl3_fx + 1.5 * pa2pb_yy_xx[j] * fl2_fx + 1.5 * pa2pb_zz_xx[j] * fl2_fx + 3.0 * pa2pb_yyzz_xx[j] * fl1_fx + 0.25 * pb_xxxx[j] * fl2_fx + 0.5 * pa2pb_yy_xxxx[j] * fl1_fx + 0.5 * pa2pb_zz_xxxx[j] * fl1_fx + pa2pb_yyzz_xxxx[j]);

                t_yyzz_xxxx[j] += fl_r_0_0 * (-0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - 1.5 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz - 0.75 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_yy[j] * fl1_fz * fl3_fx + 3.75 * pa_zz[j] * fl3_fx * fl1_fz + 9.0 * pa_yyzz[j] * fl1_fz * fl2_fx - 1.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_xx[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyzz_xx[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yy_xx[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yyzz_xx[j] * fl1_fz * fl1_fx - pb_xxxx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxxx[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxx[j] * fl2_fx * fl1_fz - pa2pb_zz_xxxx[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xxxx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxxx[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xxxx[j] * fl1_fz);

                t_yyzz_xxxy[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl3_fx + 1.5 * pa2pb_yzz_x[j] * fl2_fx + 0.375 * pb_xy[j] * fl3_fx + 0.75 * pa2pb_yy_xy[j] * fl2_fx + 0.5 * pa2pb_y_xxx[j] * fl2_fx + 0.75 * pa2pb_zz_xy[j] * fl2_fx + 1.5 * pa2pb_yyzz_xy[j] * fl1_fx + pa2pb_yzz_xxx[j] * fl1_fx + 0.25 * pb_xxxy[j] * fl2_fx + 0.5 * pa2pb_yy_xxxy[j] * fl1_fx + 0.5 * pa2pb_zz_xxxy[j] * fl1_fx + pa2pb_yyzz_xxxy[j]);

                t_yyzz_xxxy[j] += fl_r_0_0 * (-1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_yzz_x[j] * fl2_fx * fl1_fz - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xy[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yy_xy[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyzz_xy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xxx[j] * fl1_fx * fl1_fz - pb_xxxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxxy[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxy[j] * fl2_fx * fl1_fz - pa2pb_zz_xxxy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xxxy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxxy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xxxy[j] * fl1_fz);

                t_yyzz_xxxz[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl3_fx + 1.5 * pa2pb_yyz_x[j] * fl2_fx + 0.375 * pb_xz[j] * fl3_fx + 0.75 * pa2pb_yy_xz[j] * fl2_fx + 0.75 * pa2pb_zz_xz[j] * fl2_fx + 0.5 * pa2pb_z_xxx[j] * fl2_fx + 1.5 * pa2pb_yyzz_xz[j] * fl1_fx + pa2pb_yyz_xxx[j] * fl1_fx + 0.25 * pb_xxxz[j] * fl2_fx + 0.5 * pa2pb_yy_xxxz[j] * fl1_fx + 0.5 * pa2pb_zz_xxxz[j] * fl1_fx + pa2pb_yyzz_xxxz[j]);

                t_yyzz_xxxz[j] += fl_r_0_0 * (-1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_yyz_x[j] * fl1_fz * fl2_fx - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xz[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xxx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yy_xz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyzz_xz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyz_xxx[j] * fl1_fz * fl1_fx - pb_xxxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxxz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxxz[j] * fl2_fx * fl1_fz - pa2pb_zz_xxxz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xxxz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxxz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xxxz[j] * fl1_fz);

                t_yyzz_xxyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_zz[j] * fl3_fx + 0.125 * pa_yy[j] * fl3_fx + 0.5 * pa2pb_y_y[j] * fl3_fx + 0.375 * pb_xx[j] * fl3_fx + 0.25 * pa_yyzz[j] * fl2_fx + pa2pb_yzz_y[j] * fl2_fx + 0.75 * pa2pb_zz_xx[j] * fl2_fx + 0.125 * pb_yy[j] * fl3_fx + 0.25 * pa2pb_yy_xx[j] * fl2_fx + 0.25 * pa2pb_yy_yy[j] * fl2_fx + pa2pb_y_xxy[j] * fl2_fx + 0.25 * pa2pb_zz_yy[j] * fl2_fx + 0.5 * pa2pb_yyzz_xx[j] * fl1_fx + 0.5 * pa2pb_yyzz_yy[j] * fl1_fx + 2.0 * pa2pb_yzz_xxy[j] * fl1_fx + 0.25 * pb_xxyy[j] * fl2_fx + 0.5 * pa2pb_yy_xxyy[j] * fl1_fx + 0.5 * pa2pb_zz_xxyy[j] * fl1_fx + pa2pb_yyzz_xxyy[j]);

                t_yyzz_xxyy[j] += fl_r_0_0 * (-0.5 * fl3_fx * fl1_fz * fl1_fgb - 0.5 * fl3_fx * fl1_fz * fl1_fga - pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_zz[j] * fl3_fx * fl1_fz - 0.5 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyzz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_yy[j] * fl1_fz * fl3_fx - 2.0 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.0 * pa_yyzz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_yzz_y[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz - 0.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_yy[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyzz_yy[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yy_xx[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_yy_yy[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyzz_yy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yzz_xxy[j] * fl1_fx * fl1_fz - pb_xxyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxyy[j] * fl1_fz * fl1_fga + 3.0 * pb_xxyy[j] * fl2_fx * fl1_fz - pa2pb_zz_xxyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xxyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xxyy[j] * fl1_fz);

                t_yyzz_xxyz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl3_fx + 0.25 * pa2pb_y_z[j] * fl3_fx + 0.25 * pa2pb_z_y[j] * fl3_fx + 0.5 * pa2pb_yyz_y[j] * fl2_fx + 0.5 * pa2pb_yzz_z[j] * fl2_fx + pa2pb_yz_xx[j] * fl2_fx + 0.125 * pb_yz[j] * fl3_fx + 0.25 * pa2pb_yy_yz[j] * fl2_fx + 0.5 * pa2pb_y_xxz[j] * fl2_fx + 0.25 * pa2pb_zz_yz[j] * fl2_fx + 0.5 * pa2pb_z_xxy[j] * fl2_fx + 0.5 * pa2pb_yyzz_yz[j] * fl1_fx + pa2pb_yyz_xxy[j] * fl1_fx + pa2pb_yzz_xxz[j] * fl1_fx + 0.25 * pb_xxyz[j] * fl2_fx + 0.5 * pa2pb_yy_xxyz[j] * fl1_fx + 0.5 * pa2pb_zz_xxyz[j] * fl1_fx + pa2pb_yyzz_xxyz[j]);

                t_yyzz_xxyz[j] += fl_r_0_0 * (-pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 5.0 * pa_yz[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb - pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_yyz_y[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yzz_z[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz - 0.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_yz[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xxy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_yz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yy_yz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyzz_yz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyz_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xxz[j] * fl1_fx * fl1_fz - pb_xxyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxyz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxyz[j] * fl2_fx * fl1_fz - pa2pb_zz_xxyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xxyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_185_190(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (185,190)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_yy_xx = pa2pbDistances.data(1156 * idx + 207);

            auto pa2pb_yy_xy = pa2pbDistances.data(1156 * idx + 208);

            auto pa2pb_yy_xz = pa2pbDistances.data(1156 * idx + 209);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_yy_xxzz = pa2pbDistances.data(1156 * idx + 228);

            auto pa2pb_yy_xyyy = pa2pbDistances.data(1156 * idx + 229);

            auto pa2pb_yy_xyyz = pa2pbDistances.data(1156 * idx + 230);

            auto pa2pb_yy_xyzz = pa2pbDistances.data(1156 * idx + 231);

            auto pa2pb_yy_xzzz = pa2pbDistances.data(1156 * idx + 232);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_zz_xxzz = pa2pbDistances.data(1156 * idx + 296);

            auto pa2pb_zz_xyyy = pa2pbDistances.data(1156 * idx + 297);

            auto pa2pb_zz_xyyz = pa2pbDistances.data(1156 * idx + 298);

            auto pa2pb_zz_xyzz = pa2pbDistances.data(1156 * idx + 299);

            auto pa2pb_zz_xzzz = pa2pbDistances.data(1156 * idx + 300);

            auto pa2pb_yyz_x = pa2pbDistances.data(1156 * idx + 544);

            auto pa2pb_yyz_z = pa2pbDistances.data(1156 * idx + 546);

            auto pa2pb_yyz_xxz = pa2pbDistances.data(1156 * idx + 555);

            auto pa2pb_yyz_xyy = pa2pbDistances.data(1156 * idx + 556);

            auto pa2pb_yyz_xyz = pa2pbDistances.data(1156 * idx + 557);

            auto pa2pb_yyz_xzz = pa2pbDistances.data(1156 * idx + 558);

            auto pa2pb_yzz_x = pa2pbDistances.data(1156 * idx + 578);

            auto pa2pb_yzz_xyy = pa2pbDistances.data(1156 * idx + 590);

            auto pa2pb_yzz_xyz = pa2pbDistances.data(1156 * idx + 591);

            auto pa2pb_yzz_xzz = pa2pbDistances.data(1156 * idx + 592);

            auto pa2pb_yyzz_xx = pa2pbDistances.data(1156 * idx + 1057);

            auto pa2pb_yyzz_xy = pa2pbDistances.data(1156 * idx + 1058);

            auto pa2pb_yyzz_xz = pa2pbDistances.data(1156 * idx + 1059);

            auto pa2pb_yyzz_zz = pa2pbDistances.data(1156 * idx + 1062);

            auto pa2pb_yyzz_xxzz = pa2pbDistances.data(1156 * idx + 1078);

            auto pa2pb_yyzz_xyyy = pa2pbDistances.data(1156 * idx + 1079);

            auto pa2pb_yyzz_xyyz = pa2pbDistances.data(1156 * idx + 1080);

            auto pa2pb_yyzz_xyzz = pa2pbDistances.data(1156 * idx + 1081);

            auto pa2pb_yyzz_xzzz = pa2pbDistances.data(1156 * idx + 1082);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyzz_xxzz = primBuffer.data(225 * idx + 185);

            auto t_yyzz_xyyy = primBuffer.data(225 * idx + 186);

            auto t_yyzz_xyyz = primBuffer.data(225 * idx + 187);

            auto t_yyzz_xyzz = primBuffer.data(225 * idx + 188);

            auto t_yyzz_xzzz = primBuffer.data(225 * idx + 189);

            // Batch of Integrals (185,190)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xyy, pa2pb_y_xyz, pa2pb_y_xzz, \
                                     pa2pb_yy_xx, pa2pb_yy_xxzz, pa2pb_yy_xy, pa2pb_yy_xyyy, pa2pb_yy_xyyz, \
                                     pa2pb_yy_xyzz, pa2pb_yy_xz, pa2pb_yy_xzzz, pa2pb_yy_zz, pa2pb_yyz_x, pa2pb_yyz_xxz, \
                                     pa2pb_yyz_xyy, pa2pb_yyz_xyz, pa2pb_yyz_xzz, pa2pb_yyz_z, pa2pb_yyzz_xx, \
                                     pa2pb_yyzz_xxzz, pa2pb_yyzz_xy, pa2pb_yyzz_xyyy, pa2pb_yyzz_xyyz, pa2pb_yyzz_xyzz, \
                                     pa2pb_yyzz_xz, pa2pb_yyzz_xzzz, pa2pb_yyzz_zz, pa2pb_yz_xy, pa2pb_yz_xz, \
                                     pa2pb_yzz_x, pa2pb_yzz_xyy, pa2pb_yzz_xyz, pa2pb_yzz_xzz, pa2pb_z_x, pa2pb_z_xxz, \
                                     pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_z_z, pa2pb_zz_xx, pa2pb_zz_xxzz, \
                                     pa2pb_zz_xy, pa2pb_zz_xyyy, pa2pb_zz_xyyz, pa2pb_zz_xyzz, pa2pb_zz_xz, \
                                     pa2pb_zz_xzzz, pa2pb_zz_zz, pa_yy, pa_yyzz, pa_zz, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, \
                                     pb_xyzz, pb_xz, pb_xzzz, pb_zz, r_0_0, s_0_0, t_yyzz_xxzz, t_yyzz_xyyy, t_yyzz_xyyz, \
                                     t_yyzz_xyzz, t_yyzz_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyzz_xxzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.125 * pa_zz[j] * fl3_fx + 0.5 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_xx[j] * fl3_fx + 0.25 * pa_yyzz[j] * fl2_fx + pa2pb_yyz_z[j] * fl2_fx + 0.75 * pa2pb_yy_xx[j] * fl2_fx + 0.125 * pb_zz[j] * fl3_fx + 0.25 * pa2pb_yy_zz[j] * fl2_fx + 0.25 * pa2pb_zz_xx[j] * fl2_fx + 0.25 * pa2pb_zz_zz[j] * fl2_fx + pa2pb_z_xxz[j] * fl2_fx + 0.5 * pa2pb_yyzz_xx[j] * fl1_fx + 0.5 * pa2pb_yyzz_zz[j] * fl1_fx + 2.0 * pa2pb_yyz_xxz[j] * fl1_fx + 0.25 * pb_xxzz[j] * fl2_fx + 0.5 * pa2pb_yy_xxzz[j] * fl1_fx + 0.5 * pa2pb_zz_xxzz[j] * fl1_fx + pa2pb_yyzz_xxzz[j]);

                t_yyzz_xxzz[j] += fl_r_0_0 * (-0.5 * fl3_fx * fl1_fz * fl1_fgb - 0.5 * fl1_fz * fl1_fga * fl3_fx - pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_yy[j] * fl3_fx * fl1_fz - 0.25 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - pb_xx[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pa_zz[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.0 * pa_yyzz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_yyz_z[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yy_xx[j] * fl2_fx * fl1_fz - 0.25 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb + 1.25 * pb_zz[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_xxz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yyzz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yy_zz[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyzz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyz_xxz[j] * fl1_fz * fl1_fx - pb_xxzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xxzz[j] * fl2_fx * fl1_fz - pa2pb_zz_xxzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xxzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xxzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xxzz[j] * fl1_fz);

                t_yyzz_xyyy[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl3_fx + 1.125 * pb_xy[j] * fl3_fx + 1.5 * pa2pb_yzz_x[j] * fl2_fx + 2.25 * pa2pb_zz_xy[j] * fl2_fx + 0.75 * pa2pb_yy_xy[j] * fl2_fx + 1.5 * pa2pb_y_xyy[j] * fl2_fx + 1.5 * pa2pb_yyzz_xy[j] * fl1_fx + 3.0 * pa2pb_yzz_xyy[j] * fl1_fx + 0.25 * pb_xyyy[j] * fl2_fx + 0.5 * pa2pb_yy_xyyy[j] * fl1_fx + 0.5 * pa2pb_zz_xyyy[j] * fl1_fx + pa2pb_yyzz_xyyy[j]);

                t_yyzz_xyyy[j] += fl_r_0_0 * (-1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 11.25 * pb_xy[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_yzz_x[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_xy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yy_xy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyzz_xy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yzz_xyy[j] * fl1_fx * fl1_fz - pb_xyyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xyyy[j] * fl1_fz * fl1_fga + 3.0 * pb_xyyy[j] * fl2_fx * fl1_fz - pa2pb_zz_xyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xyyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xyyy[j] * fl1_fz);

                t_yyzz_xyyz[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl3_fx + 0.375 * pb_xz[j] * fl3_fx + 0.5 * pa2pb_yyz_x[j] * fl2_fx + 2.0 * pa2pb_yz_xy[j] * fl2_fx + 0.75 * pa2pb_zz_xz[j] * fl2_fx + 0.25 * pa2pb_yy_xz[j] * fl2_fx + pa2pb_y_xyz[j] * fl2_fx + 0.5 * pa2pb_z_xyy[j] * fl2_fx + 0.5 * pa2pb_yyzz_xz[j] * fl1_fx + pa2pb_yyz_xyy[j] * fl1_fx + 2.0 * pa2pb_yzz_xyz[j] * fl1_fx + 0.25 * pb_xyyz[j] * fl2_fx + 0.5 * pa2pb_yy_xyyz[j] * fl1_fx + 0.5 * pa2pb_zz_xyyz[j] * fl1_fx + pa2pb_yyzz_xyyz[j]);

                t_yyzz_xyyz[j] += fl_r_0_0 * (7.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xz[j] * fl3_fx * fl1_fz + 6.0 * pa2pb_yyz_x[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz - 0.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xyy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_xz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yy_xz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyzz_xz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyz_xyy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yzz_xyz[j] * fl1_fx * fl1_fz - pb_xyyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xyyz[j] * fl1_fz * fl1_fga + 3.0 * pb_xyyz[j] * fl2_fx * fl1_fz - pa2pb_zz_xyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xyyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xyyz[j] * fl1_fz);

                t_yyzz_xyzz[j] = fl_s_0_0 * (0.75 * pa2pb_y_x[j] * fl3_fx + 0.375 * pb_xy[j] * fl3_fx + 0.75 * pa2pb_yy_xy[j] * fl2_fx + 0.5 * pa2pb_yzz_x[j] * fl2_fx + 2.0 * pa2pb_yz_xz[j] * fl2_fx + 0.5 * pa2pb_y_xzz[j] * fl2_fx + 0.25 * pa2pb_zz_xy[j] * fl2_fx + pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_yyzz_xy[j] * fl1_fx + 2.0 * pa2pb_yyz_xyz[j] * fl1_fx + pa2pb_yzz_xzz[j] * fl1_fx + 0.25 * pb_xyzz[j] * fl2_fx + 0.5 * pa2pb_yy_xyzz[j] * fl1_fx + 0.5 * pa2pb_zz_xyzz[j] * fl1_fx + pa2pb_yyzz_xyzz[j]);

                t_yyzz_xyzz[j] += fl_r_0_0 * (7.5 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - pb_xy[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_xy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yy_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzz_x[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz - 0.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_xyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_xy[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyzz_xy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyz_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_xzz[j] * fl1_fx * fl1_fz - pb_xyzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xyzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xyzz[j] * fl2_fx * fl1_fz - pa2pb_zz_xyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xyzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xyzz[j] * fl1_fz);

                t_yyzz_xzzz[j] = fl_s_0_0 * (0.75 * pa2pb_z_x[j] * fl3_fx + 1.125 * pb_xz[j] * fl3_fx + 1.5 * pa2pb_yyz_x[j] * fl2_fx + 2.25 * pa2pb_yy_xz[j] * fl2_fx + 0.75 * pa2pb_zz_xz[j] * fl2_fx + 1.5 * pa2pb_z_xzz[j] * fl2_fx + 1.5 * pa2pb_yyzz_xz[j] * fl1_fx + 3.0 * pa2pb_yyz_xzz[j] * fl1_fx + 0.25 * pb_xzzz[j] * fl2_fx + 0.5 * pa2pb_yy_xzzz[j] * fl1_fx + 0.5 * pa2pb_zz_xzzz[j] * fl1_fx + pa2pb_yyzz_xzzz[j]);

                t_yyzz_xzzz[j] += fl_r_0_0 * (-1.5 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_x[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_xz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_yyz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 11.25 * pb_xz[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_yyz_x[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_yy_xz[j] * fl2_fx * fl1_fz - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyzz_xz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yyz_xzz[j] * fl1_fz * fl1_fx - pb_xzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_xzzz[j] * fl2_fx * fl1_fz - pa2pb_zz_xzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_xzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_xzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_190_195(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (190,195)

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

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_yy_yy = pa2pbDistances.data(1156 * idx + 210);

            auto pa2pb_yy_yz = pa2pbDistances.data(1156 * idx + 211);

            auto pa2pb_yy_zz = pa2pbDistances.data(1156 * idx + 212);

            auto pa2pb_yy_yyyy = pa2pbDistances.data(1156 * idx + 233);

            auto pa2pb_yy_yyyz = pa2pbDistances.data(1156 * idx + 234);

            auto pa2pb_yy_yyzz = pa2pbDistances.data(1156 * idx + 235);

            auto pa2pb_yy_yzzz = pa2pbDistances.data(1156 * idx + 236);

            auto pa2pb_yy_zzzz = pa2pbDistances.data(1156 * idx + 237);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_zz_yyyy = pa2pbDistances.data(1156 * idx + 301);

            auto pa2pb_zz_yyyz = pa2pbDistances.data(1156 * idx + 302);

            auto pa2pb_zz_yyzz = pa2pbDistances.data(1156 * idx + 303);

            auto pa2pb_zz_yzzz = pa2pbDistances.data(1156 * idx + 304);

            auto pa2pb_zz_zzzz = pa2pbDistances.data(1156 * idx + 305);

            auto pa2pb_yyz_y = pa2pbDistances.data(1156 * idx + 545);

            auto pa2pb_yyz_z = pa2pbDistances.data(1156 * idx + 546);

            auto pa2pb_yyz_yyy = pa2pbDistances.data(1156 * idx + 559);

            auto pa2pb_yyz_yyz = pa2pbDistances.data(1156 * idx + 560);

            auto pa2pb_yyz_yzz = pa2pbDistances.data(1156 * idx + 561);

            auto pa2pb_yyz_zzz = pa2pbDistances.data(1156 * idx + 562);

            auto pa2pb_yzz_y = pa2pbDistances.data(1156 * idx + 579);

            auto pa2pb_yzz_z = pa2pbDistances.data(1156 * idx + 580);

            auto pa2pb_yzz_yyy = pa2pbDistances.data(1156 * idx + 593);

            auto pa2pb_yzz_yyz = pa2pbDistances.data(1156 * idx + 594);

            auto pa2pb_yzz_yzz = pa2pbDistances.data(1156 * idx + 595);

            auto pa2pb_yzz_zzz = pa2pbDistances.data(1156 * idx + 596);

            auto pa2pb_yyzz_yy = pa2pbDistances.data(1156 * idx + 1060);

            auto pa2pb_yyzz_yz = pa2pbDistances.data(1156 * idx + 1061);

            auto pa2pb_yyzz_zz = pa2pbDistances.data(1156 * idx + 1062);

            auto pa2pb_yyzz_yyyy = pa2pbDistances.data(1156 * idx + 1083);

            auto pa2pb_yyzz_yyyz = pa2pbDistances.data(1156 * idx + 1084);

            auto pa2pb_yyzz_yyzz = pa2pbDistances.data(1156 * idx + 1085);

            auto pa2pb_yyzz_yzzz = pa2pbDistances.data(1156 * idx + 1086);

            auto pa2pb_yyzz_zzzz = pa2pbDistances.data(1156 * idx + 1087);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyzz_yyyy = primBuffer.data(225 * idx + 190);

            auto t_yyzz_yyyz = primBuffer.data(225 * idx + 191);

            auto t_yyzz_yyzz = primBuffer.data(225 * idx + 192);

            auto t_yyzz_yzzz = primBuffer.data(225 * idx + 193);

            auto t_yyzz_zzzz = primBuffer.data(225 * idx + 194);

            // Batch of Integrals (190,195)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, \
                                     pa2pb_y_z, pa2pb_y_zzz, pa2pb_yy_yy, pa2pb_yy_yyyy, pa2pb_yy_yyyz, \
                                     pa2pb_yy_yyzz, pa2pb_yy_yz, pa2pb_yy_yzzz, pa2pb_yy_zz, pa2pb_yy_zzzz, pa2pb_yyz_y, \
                                     pa2pb_yyz_yyy, pa2pb_yyz_yyz, pa2pb_yyz_yzz, pa2pb_yyz_z, pa2pb_yyz_zzz, \
                                     pa2pb_yyzz_yy, pa2pb_yyzz_yyyy, pa2pb_yyzz_yyyz, pa2pb_yyzz_yyzz, pa2pb_yyzz_yz, \
                                     pa2pb_yyzz_yzzz, pa2pb_yyzz_zz, pa2pb_yyzz_zzzz, pa2pb_yz_yy, pa2pb_yz_yz, \
                                     pa2pb_yz_zz, pa2pb_yzz_y, pa2pb_yzz_yyy, pa2pb_yzz_yyz, pa2pb_yzz_yzz, \
                                     pa2pb_yzz_z, pa2pb_yzz_zzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, \
                                     pa2pb_z_z, pa2pb_z_zzz, pa2pb_zz_yy, pa2pb_zz_yyyy, pa2pb_zz_yyyz, \
                                     pa2pb_zz_yyzz, pa2pb_zz_yz, pa2pb_zz_yzzz, pa2pb_zz_zz, pa2pb_zz_zzzz, pa_yy, pa_yyzz, \
                                     pa_yz, pa_zz, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, r_0_0, \
                                     s_0_0, t_yyzz_yyyy, t_yyzz_yyyz, t_yyzz_yyzz, t_yyzz_yzzz, t_yyzz_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyzz_yyyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_zz[j] * fl3_fx + 0.375 * pa_yy[j] * fl3_fx + 3.0 * pa2pb_y_y[j] * fl3_fx + 2.25 * pb_yy[j] * fl3_fx + 0.75 * pa_yyzz[j] * fl2_fx + 6.0 * pa2pb_yzz_y[j] * fl2_fx + 4.5 * pa2pb_zz_yy[j] * fl2_fx + 1.5 * pa2pb_yy_yy[j] * fl2_fx + 2.0 * pa2pb_y_yyy[j] * fl2_fx + 3.0 * pa2pb_yyzz_yy[j] * fl1_fx + 4.0 * pa2pb_yzz_yyy[j] * fl1_fx + 0.25 * pb_yyyy[j] * fl2_fx + 0.5 * pa2pb_yy_yyyy[j] * fl1_fx + 0.5 * pa2pb_zz_yyyy[j] * fl1_fx + pa2pb_yyzz_yyyy[j]);

                t_yyzz_yyyy[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 4.5 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 18.75 * pa_zz[j] * fl3_fx * fl1_fz - 1.5 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_yy[j] * fl1_fz * fl3_fx - 12.0 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 22.5 * pb_yy[j] * fl3_fx * fl1_fz + 9.0 * pa_yyzz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_yzz_y[j] * fl2_fx * fl1_fz + 54.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz - 1.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 4.0 * pa2pb_y_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyzz_yy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_yy_yy[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yyzz_yy[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_yzz_yyy[j] * fl1_fx * fl1_fz - pb_yyyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_yyyy[j] * fl1_fz * fl1_fga + 3.0 * pb_yyyy[j] * fl2_fx * fl1_fz - pa2pb_zz_yyyy[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_yyyy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_yyyy[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_yyyy[j] * fl1_fz);

                t_yyzz_yyyz[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl3_fx + 2.25 * pa2pb_z_y[j] * fl3_fx + 0.75 * pa2pb_y_z[j] * fl3_fx + 1.125 * pb_yz[j] * fl3_fx + 1.5 * pa2pb_yyz_y[j] * fl2_fx + 1.5 * pa2pb_yzz_z[j] * fl2_fx + 3.0 * pa2pb_yz_yy[j] * fl2_fx + 2.25 * pa2pb_zz_yz[j] * fl2_fx + 0.75 * pa2pb_yy_yz[j] * fl2_fx + 1.5 * pa2pb_y_yyz[j] * fl2_fx + 0.5 * pa2pb_z_yyy[j] * fl2_fx + 1.5 * pa2pb_yyzz_yz[j] * fl1_fx + pa2pb_yyz_yyy[j] * fl1_fx + 3.0 * pa2pb_yzz_yyz[j] * fl1_fx + 0.25 * pb_yyyz[j] * fl2_fx + 0.5 * pa2pb_yy_yyyz[j] * fl1_fx + 0.5 * pa2pb_zz_yyyz[j] * fl1_fx + pa2pb_yyzz_yyyz[j]);

                t_yyzz_yyyz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_yz[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 11.25 * pb_yz[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_yyz_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yzz_z[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_yyy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yy_yz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyzz_yz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyz_yyy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yzz_yyz[j] * fl1_fx * fl1_fz - pb_yyyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_yyyz[j] * fl1_fz * fl1_fga + 3.0 * pb_yyyz[j] * fl2_fx * fl1_fz - pa2pb_zz_yyyz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_yyyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_yyyz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_yyyz[j] * fl1_fz);

                t_yyzz_yyzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 1.5 * pa2pb_y_y[j] * fl3_fx + 0.375 * pa_zz[j] * fl3_fx + 1.5 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 0.375 * pb_zz[j] * fl3_fx + 0.25 * pa_yyzz[j] * fl2_fx + pa2pb_yyz_z[j] * fl2_fx + 0.75 * pa2pb_yy_yy[j] * fl2_fx + pa2pb_yzz_y[j] * fl2_fx + 4.0 * pa2pb_yz_yz[j] * fl2_fx + 0.75 * pa2pb_zz_zz[j] * fl2_fx + 0.25 * pa2pb_yy_zz[j] * fl2_fx + pa2pb_y_yzz[j] * fl2_fx + 0.25 * pa2pb_zz_yy[j] * fl2_fx + pa2pb_z_yyz[j] * fl2_fx + 0.5 * pa2pb_yyzz_yy[j] * fl1_fx + 0.5 * pa2pb_yyzz_zz[j] * fl1_fx + 2.0 * pa2pb_yyz_yyz[j] * fl1_fx + 2.0 * pa2pb_yzz_yzz[j] * fl1_fx + 0.25 * pb_yyzz[j] * fl2_fx + 0.5 * pa2pb_yy_yyzz[j] * fl1_fx + 0.5 * pa2pb_zz_yyzz[j] * fl1_fx + pa2pb_yyzz_yyzz[j]);

                t_yyzz_yyzz[j] += fl_r_0_0 * (4.5 * fl4_fx * fl1_fz - 0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb - pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 3.75 * pa_yy[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 3.75 * pa_zz[j] * fl3_fx * fl1_fz + 15.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 0.25 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 0.25 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - pb_yy[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyzz[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 2.0 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pb_yy[j] * fl3_fx * fl1_fz + 3.75 * pb_zz[j] * fl3_fx * fl1_fz + 3.0 * pa_yyzz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_yyz_z[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yy_yy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_yzz_y[j] * fl2_fx * fl1_fz + 48.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz - 0.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_yyz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_yyzz_zz[j] * fl1_fz * fl1_fgb + 3.0 * pa2pb_yy_zz[j] * fl1_fz * fl2_fx + 12.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yyzz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yyzz_zz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yyz_yyz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_yzz_yzz[j] * fl1_fx * fl1_fz - pb_yyzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_yyzz[j] * fl1_fz * fl1_fga + 3.0 * pb_yyzz[j] * fl2_fx * fl1_fz - pa2pb_zz_yyzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_yyzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_yyzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_yyzz[j] * fl1_fz);

                t_yyzz_yzzz[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl3_fx + 2.25 * pa2pb_y_z[j] * fl3_fx + 0.75 * pa2pb_z_y[j] * fl3_fx + 1.125 * pb_yz[j] * fl3_fx + 1.5 * pa2pb_yyz_y[j] * fl2_fx + 2.25 * pa2pb_yy_yz[j] * fl2_fx + 1.5 * pa2pb_yzz_z[j] * fl2_fx + 3.0 * pa2pb_yz_zz[j] * fl2_fx + 0.5 * pa2pb_y_zzz[j] * fl2_fx + 0.75 * pa2pb_zz_yz[j] * fl2_fx + 1.5 * pa2pb_z_yzz[j] * fl2_fx + 1.5 * pa2pb_yyzz_yz[j] * fl1_fx + 3.0 * pa2pb_yyz_yzz[j] * fl1_fx + pa2pb_yzz_zzz[j] * fl1_fx + 0.25 * pb_yzzz[j] * fl2_fx + 0.5 * pa2pb_yy_yzzz[j] * fl1_fx + 0.5 * pa2pb_zz_yzzz[j] * fl1_fx + pa2pb_yyzz_yzzz[j]);

                t_yyzz_yzzz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 15.0 * pa_yz[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_yz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa2pb_yyz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 11.25 * pb_yz[j] * fl3_fx * fl1_fz + 18.0 * pa2pb_yyz_y[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_yy_yz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yzz_z[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_yzz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_yz[j] * fl1_fz * fl1_fgb + 6.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yyzz_yz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yyz_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzz_zzz[j] * fl1_fx * fl1_fz - pb_yzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_yzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_yzzz[j] * fl2_fx * fl1_fz - pa2pb_zz_yzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_yzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_yzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_yzzz[j] * fl1_fz);

                t_yyzz_zzzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_yy[j] * fl3_fx + 0.375 * pa_zz[j] * fl3_fx + 3.0 * pa2pb_z_z[j] * fl3_fx + 2.25 * pb_zz[j] * fl3_fx + 0.75 * pa_yyzz[j] * fl2_fx + 6.0 * pa2pb_yyz_z[j] * fl2_fx + 4.5 * pa2pb_yy_zz[j] * fl2_fx + 1.5 * pa2pb_zz_zz[j] * fl2_fx + 2.0 * pa2pb_z_zzz[j] * fl2_fx + 3.0 * pa2pb_yyzz_zz[j] * fl1_fx + 4.0 * pa2pb_yyz_zzz[j] * fl1_fx + 0.25 * pb_zzzz[j] * fl2_fx + 0.5 * pa2pb_yy_zzzz[j] * fl1_fx + 0.5 * pa2pb_zz_zzzz[j] * fl1_fx + pa2pb_yyzz_zzzz[j]);

                t_yyzz_zzzz[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl1_fz * fl1_fga * fl3_fx - 4.5 * pa_yy[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 18.75 * pa_yy[j] * fl3_fx * fl1_fz - 0.75 * pa_yy[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_z_z[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pb_zz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyzz[j] * fl1_fx * fl1_fz * fl1_fgb - 12.0 * pa2pb_yyz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_zz[j] * fl3_fx * fl1_fz + 30.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 22.5 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa_yyzz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_yyz_z[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_yy_zz[j] * fl2_fx * fl1_fz - 1.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 4.0 * pa2pb_z_zzz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yyzz_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yyzz_zz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_yyz_zzz[j] * fl1_fz * fl1_fx - pb_zzzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_zzzz[j] * fl1_fz * fl1_fga + 3.0 * pb_zzzz[j] * fl2_fx * fl1_fz - pa2pb_zz_zzzz[j] * fl1_fz * fl1_fga + 7.0 * pa2pb_yy_zzzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zz_zzzz[j] * fl1_fx * fl1_fz + 16.0 * pa2pb_yyzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_195_200(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (195,200)

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

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_xxx = pa2pbDistances.data(1156 * idx + 43);

            auto pa2pb_y_xxy = pa2pbDistances.data(1156 * idx + 44);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_xxxx = pa2pbDistances.data(1156 * idx + 257);

            auto pa2pb_yz_xxxy = pa2pbDistances.data(1156 * idx + 258);

            auto pa2pb_yz_xxxz = pa2pbDistances.data(1156 * idx + 259);

            auto pa2pb_yz_xxyy = pa2pbDistances.data(1156 * idx + 260);

            auto pa2pb_yz_xxyz = pa2pbDistances.data(1156 * idx + 261);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_yzz_x = pa2pbDistances.data(1156 * idx + 578);

            auto pa2pb_yzz_y = pa2pbDistances.data(1156 * idx + 579);

            auto pa2pb_yzz_xxx = pa2pbDistances.data(1156 * idx + 587);

            auto pa2pb_yzz_xxy = pa2pbDistances.data(1156 * idx + 588);

            auto pa2pb_zzz_x = pa2pbDistances.data(1156 * idx + 612);

            auto pa2pb_zzz_y = pa2pbDistances.data(1156 * idx + 613);

            auto pa2pb_zzz_z = pa2pbDistances.data(1156 * idx + 614);

            auto pa2pb_zzz_xxx = pa2pbDistances.data(1156 * idx + 621);

            auto pa2pb_zzz_xxy = pa2pbDistances.data(1156 * idx + 622);

            auto pa2pb_zzz_xxz = pa2pbDistances.data(1156 * idx + 623);

            auto pa2pb_yzzz_xx = pa2pbDistances.data(1156 * idx + 1091);

            auto pa2pb_yzzz_xy = pa2pbDistances.data(1156 * idx + 1092);

            auto pa2pb_yzzz_xz = pa2pbDistances.data(1156 * idx + 1093);

            auto pa2pb_yzzz_yy = pa2pbDistances.data(1156 * idx + 1094);

            auto pa2pb_yzzz_yz = pa2pbDistances.data(1156 * idx + 1095);

            auto pa2pb_yzzz_xxxx = pa2pbDistances.data(1156 * idx + 1107);

            auto pa2pb_yzzz_xxxy = pa2pbDistances.data(1156 * idx + 1108);

            auto pa2pb_yzzz_xxxz = pa2pbDistances.data(1156 * idx + 1109);

            auto pa2pb_yzzz_xxyy = pa2pbDistances.data(1156 * idx + 1110);

            auto pa2pb_yzzz_xxyz = pa2pbDistances.data(1156 * idx + 1111);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzzz_xxxx = primBuffer.data(225 * idx + 195);

            auto t_yzzz_xxxy = primBuffer.data(225 * idx + 196);

            auto t_yzzz_xxxz = primBuffer.data(225 * idx + 197);

            auto t_yzzz_xxyy = primBuffer.data(225 * idx + 198);

            auto t_yzzz_xxyz = primBuffer.data(225 * idx + 199);

            // Batch of Integrals (195,200)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xxx, pa2pb_y_xxy, pa2pb_y_y, \
                                     pa2pb_yz_xx, pa2pb_yz_xxxx, pa2pb_yz_xxxy, pa2pb_yz_xxxz, pa2pb_yz_xxyy, \
                                     pa2pb_yz_xxyz, pa2pb_yz_xy, pa2pb_yz_xz, pa2pb_yz_yy, pa2pb_yz_yz, pa2pb_yzz_x, \
                                     pa2pb_yzz_xxx, pa2pb_yzz_xxy, pa2pb_yzz_y, pa2pb_yzzz_xx, pa2pb_yzzz_xxxx, \
                                     pa2pb_yzzz_xxxy, pa2pb_yzzz_xxxz, pa2pb_yzzz_xxyy, pa2pb_yzzz_xxyz, pa2pb_yzzz_xy, \
                                     pa2pb_yzzz_xz, pa2pb_yzzz_yy, pa2pb_yzzz_yz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, \
                                     pa2pb_z_xxz, pa2pb_z_y, pa2pb_z_z, pa2pb_zz_xx, pa2pb_zzz_x, pa2pb_zzz_xxx, \
                                     pa2pb_zzz_xxy, pa2pb_zzz_xxz, pa2pb_zzz_y, pa2pb_zzz_z, pa_yz, pa_yzzz, pa_zz, pb_xx, \
                                     r_0_0, s_0_0, t_yzzz_xxxx, t_yzzz_xxxy, t_yzzz_xxxz, t_yzzz_xxyy, t_yzzz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yzzz_xxxx[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 0.75 * pa_yzzz[j] * fl2_fx + 4.5 * pa2pb_yz_xx[j] * fl2_fx + 3.0 * pa2pb_yzzz_xx[j] * fl1_fx + 1.5 * pa2pb_yz_xxxx[j] * fl1_fx + pa2pb_yzzz_xxxx[j]);

                t_yzzz_xxxx[j] += fl_r_0_0 * (-4.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_yz[j] * fl1_fz * fl3_fx + 9.0 * pa_yzzz[j] * fl1_fz * fl2_fx - 9.0 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_yzzz_xx[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_yz_xx[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_yzzz_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxxx[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxxx[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xxxx[j] * fl1_fz);

                t_yzzz_xxxy[j] = fl_s_0_0 * (1.125 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa2pb_zzz_x[j] * fl2_fx + 2.25 * pa2pb_yz_xy[j] * fl2_fx + 0.75 * pa2pb_z_xxx[j] * fl2_fx + 1.5 * pa2pb_yzzz_xy[j] * fl1_fx + 0.5 * pa2pb_zzz_xxx[j] * fl1_fx + 1.5 * pa2pb_yz_xxxy[j] * fl1_fx + pa2pb_yzzz_xxxy[j]);

                t_yzzz_xxxy[j] += fl_r_0_0 * (-2.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_zzz_x[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzzz_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_yz_xy[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yzzz_xy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzz_xxx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xxxy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxxy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xxxy[j] * fl1_fz);

                t_yzzz_xxxz[j] = fl_s_0_0 * (1.125 * pa2pb_y_x[j] * fl3_fx + 2.25 * pa2pb_yzz_x[j] * fl2_fx + 2.25 * pa2pb_yz_xz[j] * fl2_fx + 0.75 * pa2pb_y_xxx[j] * fl2_fx + 1.5 * pa2pb_yzzz_xz[j] * fl1_fx + 1.5 * pa2pb_yzz_xxx[j] * fl1_fx + 1.5 * pa2pb_yz_xxxz[j] * fl1_fx + pa2pb_yzzz_xxxz[j]);

                t_yzzz_xxxz[j] += fl_r_0_0 * (-2.25 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_yzz_x[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzzz_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_yz_xz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_y_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yzzz_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yzz_xxx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxxz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxxz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xxxz[j] * fl1_fz);

                t_yzzz_xxyy[j] = fl_s_0_0 * (0.375 * pa_yz[j] * fl3_fx + 0.75 * pa2pb_z_y[j] * fl3_fx + 0.25 * pa_yzzz[j] * fl2_fx + 0.5 * pa2pb_zzz_y[j] * fl2_fx + 0.75 * pa2pb_yz_xx[j] * fl2_fx + 0.75 * pa2pb_yz_yy[j] * fl2_fx + 1.5 * pa2pb_z_xxy[j] * fl2_fx + 0.5 * pa2pb_yzzz_xx[j] * fl1_fx + 0.5 * pa2pb_yzzz_yy[j] * fl1_fx + pa2pb_zzz_xxy[j] * fl1_fx + 1.5 * pa2pb_yz_xxyy[j] * fl1_fx + pa2pb_yzzz_xxyy[j]);

                t_yzzz_xxyy[j] += fl_r_0_0 * (-1.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_yzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa_yz[j] * fl1_fz * fl3_fx - pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 3.0 * pa_yzzz[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_zzz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yzzz_yy[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_xx[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_yz_yy[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yzzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yzzz_yy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xxy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xxyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xxyy[j] * fl1_fz);

                t_yzzz_xxyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_zz[j] * fl3_fx + 0.375 * pa2pb_y_y[j] * fl3_fx + 0.375 * pa2pb_z_z[j] * fl3_fx + 0.375 * pb_xx[j] * fl3_fx + 0.75 * pa2pb_yzz_y[j] * fl2_fx + 0.25 * pa2pb_zzz_z[j] * fl2_fx + 0.75 * pa2pb_zz_xx[j] * fl2_fx + 0.75 * pa2pb_yz_yz[j] * fl2_fx + 0.75 * pa2pb_y_xxy[j] * fl2_fx + 0.75 * pa2pb_z_xxz[j] * fl2_fx + 0.5 * pa2pb_yzzz_yz[j] * fl1_fx + 1.5 * pa2pb_yzz_xxy[j] * fl1_fx + 0.5 * pa2pb_zzz_xxz[j] * fl1_fx + 1.5 * pa2pb_yz_xxyz[j] * fl1_fx + pa2pb_yzzz_xxyz[j]);

                t_yzzz_xxyz[j] += fl_r_0_0 * (-0.375 * fl3_fx * fl1_fz * fl1_fgb - 0.375 * fl3_fx * fl1_fz * fl1_fga - 0.75 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 1.5 * fl4_fx * fl1_fz + 3.75 * pa_zz[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yzz_y[j] * fl1_fz * fl2_fx + 3.0 * pa2pb_zzz_z[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_yz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_y_xxy[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yzzz_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yzz_xxy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzz_xxz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xxyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_200_205(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (200,205)

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

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_x = pa2pbDistances.data(1156 * idx + 34);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_xxz = pa2pbDistances.data(1156 * idx + 45);

            auto pa2pb_y_xyy = pa2pbDistances.data(1156 * idx + 46);

            auto pa2pb_y_xyz = pa2pbDistances.data(1156 * idx + 47);

            auto pa2pb_y_xzz = pa2pbDistances.data(1156 * idx + 48);

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_yz_xx = pa2pbDistances.data(1156 * idx + 241);

            auto pa2pb_yz_xy = pa2pbDistances.data(1156 * idx + 242);

            auto pa2pb_yz_xz = pa2pbDistances.data(1156 * idx + 243);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_yz_xxzz = pa2pbDistances.data(1156 * idx + 262);

            auto pa2pb_yz_xyyy = pa2pbDistances.data(1156 * idx + 263);

            auto pa2pb_yz_xyyz = pa2pbDistances.data(1156 * idx + 264);

            auto pa2pb_yz_xyzz = pa2pbDistances.data(1156 * idx + 265);

            auto pa2pb_yz_xzzz = pa2pbDistances.data(1156 * idx + 266);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_yzz_x = pa2pbDistances.data(1156 * idx + 578);

            auto pa2pb_yzz_z = pa2pbDistances.data(1156 * idx + 580);

            auto pa2pb_yzz_xxz = pa2pbDistances.data(1156 * idx + 589);

            auto pa2pb_yzz_xyy = pa2pbDistances.data(1156 * idx + 590);

            auto pa2pb_yzz_xyz = pa2pbDistances.data(1156 * idx + 591);

            auto pa2pb_yzz_xzz = pa2pbDistances.data(1156 * idx + 592);

            auto pa2pb_zzz_x = pa2pbDistances.data(1156 * idx + 612);

            auto pa2pb_zzz_xyy = pa2pbDistances.data(1156 * idx + 624);

            auto pa2pb_zzz_xyz = pa2pbDistances.data(1156 * idx + 625);

            auto pa2pb_zzz_xzz = pa2pbDistances.data(1156 * idx + 626);

            auto pa2pb_yzzz_xx = pa2pbDistances.data(1156 * idx + 1091);

            auto pa2pb_yzzz_xy = pa2pbDistances.data(1156 * idx + 1092);

            auto pa2pb_yzzz_xz = pa2pbDistances.data(1156 * idx + 1093);

            auto pa2pb_yzzz_zz = pa2pbDistances.data(1156 * idx + 1096);

            auto pa2pb_yzzz_xxzz = pa2pbDistances.data(1156 * idx + 1112);

            auto pa2pb_yzzz_xyyy = pa2pbDistances.data(1156 * idx + 1113);

            auto pa2pb_yzzz_xyyz = pa2pbDistances.data(1156 * idx + 1114);

            auto pa2pb_yzzz_xyzz = pa2pbDistances.data(1156 * idx + 1115);

            auto pa2pb_yzzz_xzzz = pa2pbDistances.data(1156 * idx + 1116);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzzz_xxzz = primBuffer.data(225 * idx + 200);

            auto t_yzzz_xyyy = primBuffer.data(225 * idx + 201);

            auto t_yzzz_xyyz = primBuffer.data(225 * idx + 202);

            auto t_yzzz_xyzz = primBuffer.data(225 * idx + 203);

            auto t_yzzz_xzzz = primBuffer.data(225 * idx + 204);

            // Batch of Integrals (200,205)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_x, pa2pb_y_xxz, pa2pb_y_xyy, pa2pb_y_xyz, \
                                     pa2pb_y_xzz, pa2pb_y_z, pa2pb_yz_xx, pa2pb_yz_xxzz, pa2pb_yz_xy, pa2pb_yz_xyyy, \
                                     pa2pb_yz_xyyz, pa2pb_yz_xyzz, pa2pb_yz_xz, pa2pb_yz_xzzz, pa2pb_yz_zz, pa2pb_yzz_x, \
                                     pa2pb_yzz_xxz, pa2pb_yzz_xyy, pa2pb_yzz_xyz, pa2pb_yzz_xzz, pa2pb_yzz_z, \
                                     pa2pb_yzzz_xx, pa2pb_yzzz_xxzz, pa2pb_yzzz_xy, pa2pb_yzzz_xyyy, pa2pb_yzzz_xyyz, \
                                     pa2pb_yzzz_xyzz, pa2pb_yzzz_xz, pa2pb_yzzz_xzzz, pa2pb_yzzz_zz, pa2pb_z_x, \
                                     pa2pb_z_xyy, pa2pb_z_xyz, pa2pb_z_xzz, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zzz_x, \
                                     pa2pb_zzz_xyy, pa2pb_zzz_xyz, pa2pb_zzz_xzz, pa_yz, pa_yzzz, pb_xy, pb_xz, r_0_0, s_0_0, \
                                     t_yzzz_xxzz, t_yzzz_xyyy, t_yzzz_xyyz, t_yzzz_xyzz, t_yzzz_xzzz: VLX_ALIGN)
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

                t_yzzz_xxzz[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 0.75 * pa2pb_y_z[j] * fl3_fx + 0.25 * pa_yzzz[j] * fl2_fx + 1.5 * pa2pb_yzz_z[j] * fl2_fx + 2.25 * pa2pb_yz_xx[j] * fl2_fx + 0.75 * pa2pb_yz_zz[j] * fl2_fx + 1.5 * pa2pb_y_xxz[j] * fl2_fx + 0.5 * pa2pb_yzzz_xx[j] * fl1_fx + 0.5 * pa2pb_yzzz_zz[j] * fl1_fx + 3.0 * pa2pb_yzz_xxz[j] * fl1_fx + 1.5 * pa2pb_yz_xxzz[j] * fl1_fx + pa2pb_yzzz_xxzz[j]);

                t_yzzz_xxzz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_yz[j] * fl3_fx * fl1_fz - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - pa_yzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz + 3.0 * pa_yzzz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yzz_z[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_yz_xx[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_yzzz_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_zz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_y_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yzzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yzzz_zz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yzz_xxz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xxzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xxzz[j] * fl1_fz);

                t_yzzz_xyyy[j] = fl_s_0_0 * (1.125 * pa2pb_z_x[j] * fl3_fx + 0.75 * pa2pb_zzz_x[j] * fl2_fx + 2.25 * pa2pb_yz_xy[j] * fl2_fx + 2.25 * pa2pb_z_xyy[j] * fl2_fx + 1.5 * pa2pb_yzzz_xy[j] * fl1_fx + 1.5 * pa2pb_zzz_xyy[j] * fl1_fx + 1.5 * pa2pb_yz_xyyy[j] * fl1_fx + pa2pb_yzzz_xyyy[j]);

                t_yzzz_xyyy[j] += fl_r_0_0 * (-2.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_zzz_x[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzzz_xy[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_yz_xy[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yzzz_xy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_zzz_xyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xyyy[j] * fl1_fz);

                t_yzzz_xyyz[j] = fl_s_0_0 * (0.375 * pa2pb_y_x[j] * fl3_fx + 0.75 * pb_xy[j] * fl3_fx + 0.75 * pa2pb_yzz_x[j] * fl2_fx + 1.5 * pa2pb_zz_xy[j] * fl2_fx + 0.75 * pa2pb_yz_xz[j] * fl2_fx + 0.75 * pa2pb_y_xyy[j] * fl2_fx + 1.5 * pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_yzzz_xz[j] * fl1_fx + 1.5 * pa2pb_yzz_xyy[j] * fl1_fx + pa2pb_zzz_xyz[j] * fl1_fx + 1.5 * pa2pb_yz_xyyz[j] * fl1_fx + pa2pb_yzzz_xyyz[j]);

                t_yzzz_xyyz[j] += fl_r_0_0 * (-0.75 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.75 * pa2pb_y_x[j] * fl3_fx * fl1_fz + 7.5 * pb_xy[j] * fl3_fx * fl1_fz + 9.0 * pa2pb_yzz_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_z_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_xz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_xz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_y_xyy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yzzz_xz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yzz_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_xyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xyyz[j] * fl1_fz);

                t_yzzz_xyzz[j] = fl_s_0_0 * (1.125 * pa2pb_z_x[j] * fl3_fx + 0.75 * pb_xz[j] * fl3_fx + 2.25 * pa2pb_yz_xy[j] * fl2_fx + 0.25 * pa2pb_zzz_x[j] * fl2_fx + 1.5 * pa2pb_zz_xz[j] * fl2_fx + 1.5 * pa2pb_y_xyz[j] * fl2_fx + 0.75 * pa2pb_z_xzz[j] * fl2_fx + 0.5 * pa2pb_yzzz_xy[j] * fl1_fx + 3.0 * pa2pb_yzz_xyz[j] * fl1_fx + 0.5 * pa2pb_zzz_xzz[j] * fl1_fx + 1.5 * pa2pb_yz_xyzz[j] * fl1_fx + pa2pb_yzzz_xyzz[j]);

                t_yzzz_xyzz[j] += fl_r_0_0 * (11.25 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 0.75 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pb_xz[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_yz_xy[j] * fl2_fx * fl1_fz + 3.0 * pa2pb_zzz_x[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_xy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_y_xyz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yzzz_xy[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yzz_xyz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzz_xzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xyzz[j] * fl1_fz);

                t_yzzz_xzzz[j] = fl_s_0_0 * (1.875 * pa2pb_y_x[j] * fl3_fx + 2.25 * pa2pb_yzz_x[j] * fl2_fx + 6.75 * pa2pb_yz_xz[j] * fl2_fx + 2.25 * pa2pb_y_xzz[j] * fl2_fx + 1.5 * pa2pb_yzzz_xz[j] * fl1_fx + 4.5 * pa2pb_yzz_xzz[j] * fl1_fx + 1.5 * pa2pb_yz_xzzz[j] * fl1_fx + pa2pb_yzzz_xzzz[j]);

                t_yzzz_xzzz[j] += fl_r_0_0 * (18.75 * pa2pb_y_x[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 27.0 * pa2pb_yzz_x[j] * fl1_fz * fl2_fx + 81.0 * pa2pb_yz_xz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzzz_xz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_y_xzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yzzz_xz[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_yzz_xzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_xzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_205_210(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (205,210)

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

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_y = pa2pbDistances.data(1156 * idx + 35);

            auto pa2pb_y_z = pa2pbDistances.data(1156 * idx + 36);

            auto pa2pb_y_yyy = pa2pbDistances.data(1156 * idx + 49);

            auto pa2pb_y_yyz = pa2pbDistances.data(1156 * idx + 50);

            auto pa2pb_y_yzz = pa2pbDistances.data(1156 * idx + 51);

            auto pa2pb_y_zzz = pa2pbDistances.data(1156 * idx + 52);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_yz_yy = pa2pbDistances.data(1156 * idx + 244);

            auto pa2pb_yz_yz = pa2pbDistances.data(1156 * idx + 245);

            auto pa2pb_yz_zz = pa2pbDistances.data(1156 * idx + 246);

            auto pa2pb_yz_yyyy = pa2pbDistances.data(1156 * idx + 267);

            auto pa2pb_yz_yyyz = pa2pbDistances.data(1156 * idx + 268);

            auto pa2pb_yz_yyzz = pa2pbDistances.data(1156 * idx + 269);

            auto pa2pb_yz_yzzz = pa2pbDistances.data(1156 * idx + 270);

            auto pa2pb_yz_zzzz = pa2pbDistances.data(1156 * idx + 271);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_yzz_y = pa2pbDistances.data(1156 * idx + 579);

            auto pa2pb_yzz_z = pa2pbDistances.data(1156 * idx + 580);

            auto pa2pb_yzz_yyy = pa2pbDistances.data(1156 * idx + 593);

            auto pa2pb_yzz_yyz = pa2pbDistances.data(1156 * idx + 594);

            auto pa2pb_yzz_yzz = pa2pbDistances.data(1156 * idx + 595);

            auto pa2pb_yzz_zzz = pa2pbDistances.data(1156 * idx + 596);

            auto pa2pb_zzz_y = pa2pbDistances.data(1156 * idx + 613);

            auto pa2pb_zzz_z = pa2pbDistances.data(1156 * idx + 614);

            auto pa2pb_zzz_yyy = pa2pbDistances.data(1156 * idx + 627);

            auto pa2pb_zzz_yyz = pa2pbDistances.data(1156 * idx + 628);

            auto pa2pb_zzz_yzz = pa2pbDistances.data(1156 * idx + 629);

            auto pa2pb_zzz_zzz = pa2pbDistances.data(1156 * idx + 630);

            auto pa2pb_yzzz_yy = pa2pbDistances.data(1156 * idx + 1094);

            auto pa2pb_yzzz_yz = pa2pbDistances.data(1156 * idx + 1095);

            auto pa2pb_yzzz_zz = pa2pbDistances.data(1156 * idx + 1096);

            auto pa2pb_yzzz_yyyy = pa2pbDistances.data(1156 * idx + 1117);

            auto pa2pb_yzzz_yyyz = pa2pbDistances.data(1156 * idx + 1118);

            auto pa2pb_yzzz_yyzz = pa2pbDistances.data(1156 * idx + 1119);

            auto pa2pb_yzzz_yzzz = pa2pbDistances.data(1156 * idx + 1120);

            auto pa2pb_yzzz_zzzz = pa2pbDistances.data(1156 * idx + 1121);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzzz_yyyy = primBuffer.data(225 * idx + 205);

            auto t_yzzz_yyyz = primBuffer.data(225 * idx + 206);

            auto t_yzzz_yyzz = primBuffer.data(225 * idx + 207);

            auto t_yzzz_yzzz = primBuffer.data(225 * idx + 208);

            auto t_yzzz_zzzz = primBuffer.data(225 * idx + 209);

            // Batch of Integrals (205,210)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_y, pa2pb_y_yyy, pa2pb_y_yyz, pa2pb_y_yzz, \
                                     pa2pb_y_z, pa2pb_y_zzz, pa2pb_yz_yy, pa2pb_yz_yyyy, pa2pb_yz_yyyz, \
                                     pa2pb_yz_yyzz, pa2pb_yz_yz, pa2pb_yz_yzzz, pa2pb_yz_zz, pa2pb_yz_zzzz, pa2pb_yzz_y, \
                                     pa2pb_yzz_yyy, pa2pb_yzz_yyz, pa2pb_yzz_yzz, pa2pb_yzz_z, pa2pb_yzz_zzz, \
                                     pa2pb_yzzz_yy, pa2pb_yzzz_yyyy, pa2pb_yzzz_yyyz, pa2pb_yzzz_yyzz, pa2pb_yzzz_yz, \
                                     pa2pb_yzzz_yzzz, pa2pb_yzzz_zz, pa2pb_yzzz_zzzz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, \
                                     pa2pb_z_yzz, pa2pb_z_z, pa2pb_z_zzz, pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zz_zz, \
                                     pa2pb_zzz_y, pa2pb_zzz_yyy, pa2pb_zzz_yyz, pa2pb_zzz_yzz, pa2pb_zzz_z, \
                                     pa2pb_zzz_zzz, pa_yz, pa_yzzz, pa_zz, pb_yy, pb_yz, pb_zz, r_0_0, s_0_0, t_yzzz_yyyy, \
                                     t_yzzz_yyyz, t_yzzz_yyzz, t_yzzz_yzzz, t_yzzz_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yzzz_yyyy[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 4.5 * pa2pb_z_y[j] * fl3_fx + 0.75 * pa_yzzz[j] * fl2_fx + 3.0 * pa2pb_zzz_y[j] * fl2_fx + 4.5 * pa2pb_yz_yy[j] * fl2_fx + 3.0 * pa2pb_z_yyy[j] * fl2_fx + 3.0 * pa2pb_yzzz_yy[j] * fl1_fx + 2.0 * pa2pb_zzz_yyy[j] * fl1_fx + 1.5 * pa2pb_yz_yyyy[j] * fl1_fx + pa2pb_yzzz_yyyy[j]);

                t_yzzz_yyyy[j] += fl_r_0_0 * (-4.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_yzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa_yz[j] * fl1_fz * fl3_fx - 6.0 * pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 9.0 * pa_yzzz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_zzz_y[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yzzz_yy[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_yz_yy[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yzzz_yy[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_zzz_yyy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyyy[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_yyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_yyyy[j] * fl1_fz);

                t_yzzz_yyyz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_zz[j] * fl3_fx + 1.125 * pa2pb_y_y[j] * fl3_fx + 1.125 * pa2pb_z_z[j] * fl3_fx + 1.125 * pb_yy[j] * fl3_fx + 2.25 * pa2pb_yzz_y[j] * fl2_fx + 0.75 * pa2pb_zzz_z[j] * fl2_fx + 2.25 * pa2pb_zz_yy[j] * fl2_fx + 2.25 * pa2pb_yz_yz[j] * fl2_fx + 0.75 * pa2pb_y_yyy[j] * fl2_fx + 2.25 * pa2pb_z_yyz[j] * fl2_fx + 1.5 * pa2pb_yzzz_yz[j] * fl1_fx + 1.5 * pa2pb_yzz_yyy[j] * fl1_fx + 1.5 * pa2pb_zzz_yyz[j] * fl1_fx + 1.5 * pa2pb_yz_yyyz[j] * fl1_fx + pa2pb_yzzz_yyyz[j]);

                t_yzzz_yyyz[j] += fl_r_0_0 * (-1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 4.5 * fl4_fx * fl1_fz + 11.25 * pa_zz[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_y_y[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 11.25 * pb_yy[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_yzz_y[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_zzz_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_z_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzzz_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_yz_yz[j] * fl1_fz * fl2_fx + 9.0 * pa2pb_y_yyy[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yzzz_yz[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_yzz_yyy[j] * fl1_fz * fl1_fx + 21.0 * pa2pb_zzz_yyz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyyz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_yyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_yyyz[j] * fl1_fz);

                t_yzzz_yyzz[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 2.25 * pa2pb_z_y[j] * fl3_fx + 0.75 * pa2pb_y_z[j] * fl3_fx + 1.5 * pb_yz[j] * fl3_fx + 0.25 * pa_yzzz[j] * fl2_fx + 1.5 * pa2pb_yzz_z[j] * fl2_fx + 2.25 * pa2pb_yz_yy[j] * fl2_fx + 0.5 * pa2pb_zzz_y[j] * fl2_fx + 3.0 * pa2pb_zz_yz[j] * fl2_fx + 0.75 * pa2pb_yz_zz[j] * fl2_fx + 1.5 * pa2pb_y_yyz[j] * fl2_fx + 1.5 * pa2pb_z_yzz[j] * fl2_fx + 0.5 * pa2pb_yzzz_yy[j] * fl1_fx + 0.5 * pa2pb_yzzz_zz[j] * fl1_fx + 3.0 * pa2pb_yzz_yyz[j] * fl1_fx + pa2pb_zzz_yzz[j] * fl1_fx + 1.5 * pa2pb_yz_yyzz[j] * fl1_fx + pa2pb_yzzz_yyzz[j]);

                t_yzzz_yyzz[j] += fl_r_0_0 * (-3.0 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 11.25 * pa_yz[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 0.75 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - pa_yzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_z[j] * fl3_fx * fl1_fz - pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pb_yz[j] * fl3_fx * fl1_fz + 3.0 * pa_yzzz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yzz_z[j] * fl1_fz * fl2_fx + 27.0 * pa2pb_yz_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzz_y[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_z_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_yzzz_zz[j] * fl1_fz * fl1_fgb + 9.0 * pa2pb_yz_zz[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_y_yyz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_yzzz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_yzzz_zz[j] * fl1_fz * fl1_fx + 42.0 * pa2pb_yzz_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzz_yzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_yyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_yyzz[j] * fl1_fz);

                t_yzzz_yzzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa2pb_y_y[j] * fl3_fx + 1.125 * pa_zz[j] * fl3_fx + 3.375 * pa2pb_z_z[j] * fl3_fx + 1.125 * pb_zz[j] * fl3_fx + 2.25 * pa2pb_yzz_y[j] * fl2_fx + 6.75 * pa2pb_yz_yz[j] * fl2_fx + 0.75 * pa2pb_zzz_z[j] * fl2_fx + 2.25 * pa2pb_zz_zz[j] * fl2_fx + 2.25 * pa2pb_y_yzz[j] * fl2_fx + 0.75 * pa2pb_z_zzz[j] * fl2_fx + 1.5 * pa2pb_yzzz_yz[j] * fl1_fx + 4.5 * pa2pb_yzz_yzz[j] * fl1_fx + 0.5 * pa2pb_zzz_zzz[j] * fl1_fx + 1.5 * pa2pb_yz_yzzz[j] * fl1_fx + pa2pb_yzzz_yzzz[j]);

                t_yzzz_yzzz[j] += fl_r_0_0 * (7.5 * fl4_fx * fl1_fz - 1.125 * fl3_fx * fl1_fz * fl1_fgb - 1.125 * fl3_fx * fl1_fz * fl1_fga - 2.25 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 18.75 * pa2pb_y_y[j] * fl3_fx * fl1_fz + 11.25 * pa_zz[j] * fl3_fx * fl1_fz + 33.75 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_y_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.25 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa2pb_yzz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 11.25 * pb_zz[j] * fl3_fx * fl1_fz + 27.0 * pa2pb_yzz_y[j] * fl1_fz * fl2_fx + 81.0 * pa2pb_yz_yz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_zzz_z[j] * fl2_fx * fl1_fz + 27.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzzz_yz[j] * fl1_fz * fl1_fgb + 27.0 * pa2pb_y_yzz[j] * fl2_fx * fl1_fz + 9.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_yzzz_yz[j] * fl1_fz * fl1_fx + 63.0 * pa2pb_yzz_yzz[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzz_zzz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_yzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_yzzz[j] * fl1_fz);

                t_yzzz_zzzz[j] = fl_s_0_0 * (5.625 * pa_yz[j] * fl3_fx + 7.5 * pa2pb_y_z[j] * fl3_fx + 0.75 * pa_yzzz[j] * fl2_fx + 9.0 * pa2pb_yzz_z[j] * fl2_fx + 13.5 * pa2pb_yz_zz[j] * fl2_fx + 3.0 * pa2pb_y_zzz[j] * fl2_fx + 3.0 * pa2pb_yzzz_zz[j] * fl1_fx + 6.0 * pa2pb_yzz_zzz[j] * fl1_fx + 1.5 * pa2pb_yz_zzzz[j] * fl1_fx + pa2pb_yzzz_zzzz[j]);

                t_yzzz_zzzz[j] += fl_r_0_0 * (-13.5 * pa_yz[j] * fl2_fx * fl1_fz * fl1_fgb + 56.25 * pa_yz[j] * fl3_fx * fl1_fz + 75.0 * pa2pb_y_z[j] * fl3_fx * fl1_fz - 2.25 * pa_yz[j] * fl1_fz * fl1_fga * fl2_fx - 9.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_y_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_yzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_yzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_yzzz[j] * fl1_fz * fl2_fx + 108.0 * pa2pb_yzz_z[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_yz_zz[j] * fl2_fx * fl1_fz - 9.0 * pa2pb_yz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yzzz_zz[j] * fl1_fz * fl1_fgb + 36.0 * pa2pb_y_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_yzzz_zz[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_yzz_zzz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_zzzz[j] * fl1_fz * fl1_fga + 21.0 * pa2pb_yz_zzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_yzzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_210_215(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (210,215)

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

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

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

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_xxx = pa2pbDistances.data(1156 * idx + 77);

            auto pa2pb_z_xxy = pa2pbDistances.data(1156 * idx + 78);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_xxxx = pa2pbDistances.data(1156 * idx + 291);

            auto pa2pb_zz_xxxy = pa2pbDistances.data(1156 * idx + 292);

            auto pa2pb_zz_xxxz = pa2pbDistances.data(1156 * idx + 293);

            auto pa2pb_zz_xxyy = pa2pbDistances.data(1156 * idx + 294);

            auto pa2pb_zz_xxyz = pa2pbDistances.data(1156 * idx + 295);

            auto pa2pb_zzz_x = pa2pbDistances.data(1156 * idx + 612);

            auto pa2pb_zzz_y = pa2pbDistances.data(1156 * idx + 613);

            auto pa2pb_zzz_xxx = pa2pbDistances.data(1156 * idx + 621);

            auto pa2pb_zzz_xxy = pa2pbDistances.data(1156 * idx + 622);

            auto pa2pb_zzzz_xx = pa2pbDistances.data(1156 * idx + 1125);

            auto pa2pb_zzzz_xy = pa2pbDistances.data(1156 * idx + 1126);

            auto pa2pb_zzzz_xz = pa2pbDistances.data(1156 * idx + 1127);

            auto pa2pb_zzzz_yy = pa2pbDistances.data(1156 * idx + 1128);

            auto pa2pb_zzzz_yz = pa2pbDistances.data(1156 * idx + 1129);

            auto pa2pb_zzzz_xxxx = pa2pbDistances.data(1156 * idx + 1141);

            auto pa2pb_zzzz_xxxy = pa2pbDistances.data(1156 * idx + 1142);

            auto pa2pb_zzzz_xxxz = pa2pbDistances.data(1156 * idx + 1143);

            auto pa2pb_zzzz_xxyy = pa2pbDistances.data(1156 * idx + 1144);

            auto pa2pb_zzzz_xxyz = pa2pbDistances.data(1156 * idx + 1145);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzzz_xxxx = primBuffer.data(225 * idx + 210);

            auto t_zzzz_xxxy = primBuffer.data(225 * idx + 211);

            auto t_zzzz_xxxz = primBuffer.data(225 * idx + 212);

            auto t_zzzz_xxyy = primBuffer.data(225 * idx + 213);

            auto t_zzzz_xxyz = primBuffer.data(225 * idx + 214);

            // Batch of Integrals (210,215)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_x, pa2pb_z_xxx, pa2pb_z_xxy, pa2pb_z_y, \
                                     pa2pb_zz_xx, pa2pb_zz_xxxx, pa2pb_zz_xxxy, pa2pb_zz_xxxz, pa2pb_zz_xxyy, \
                                     pa2pb_zz_xxyz, pa2pb_zz_xy, pa2pb_zz_xz, pa2pb_zz_yy, pa2pb_zz_yz, pa2pb_zzz_x, \
                                     pa2pb_zzz_xxx, pa2pb_zzz_xxy, pa2pb_zzz_y, pa2pb_zzzz_xx, pa2pb_zzzz_xxxx, \
                                     pa2pb_zzzz_xxxy, pa2pb_zzzz_xxxz, pa2pb_zzzz_xxyy, pa2pb_zzzz_xxyz, pa2pb_zzzz_xy, \
                                     pa2pb_zzzz_xz, pa2pb_zzzz_yy, pa2pb_zzzz_yz, pa_zz, pa_zzzz, pb_xx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_yy, pb_yz, r_0_0, s_0_0, t_zzzz_xxxx, \
                                     t_zzzz_xxxy, t_zzzz_xxxz, t_zzzz_xxyy, t_zzzz_xxyz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_zzzz_xxxx[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_zz[j] * fl3_fx + 0.75 * pa_zzzz[j] * fl2_fx + 2.25 * pb_xx[j] * fl3_fx + 9.0 * pa2pb_zz_xx[j] * fl2_fx + 3.0 * pa2pb_zzzz_xx[j] * fl1_fx + 0.75 * pb_xxxx[j] * fl2_fx + 3.0 * pa2pb_zz_xxxx[j] * fl1_fx + pa2pb_zzzz_xxxx[j]);

                t_zzzz_xxxx[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 9.0 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx + 4.5 * fl4_fx * fl1_fz - 3.0 * pa_zzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_zz[j] * fl1_fz * fl3_fx + 9.0 * pa_zzzz[j] * fl1_fz * fl2_fx - 4.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 18.0 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx + 22.5 * pb_xx[j] * fl3_fx * fl1_fz - 6.0 * pa2pb_zzzz_xx[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_zz_xx[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_zzzz_xx[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxxx[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxx[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xxxx[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xxxx[j] * fl1_fz);

                t_zzzz_xxxy[j] = fl_s_0_0 * (1.125 * pb_xy[j] * fl3_fx + 4.5 * pa2pb_zz_xy[j] * fl2_fx + 1.5 * pa2pb_zzzz_xy[j] * fl1_fx + 0.75 * pb_xxxy[j] * fl2_fx + 3.0 * pa2pb_zz_xxxy[j] * fl1_fx + pa2pb_zzzz_xxxy[j]);

                t_zzzz_xxxy[j] += fl_r_0_0 * (-2.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx + 11.25 * pb_xy[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_zzzz_xy[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_zz_xy[j] * fl1_fz * fl2_fx + 21.0 * pa2pb_zzzz_xy[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxxy[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xxxy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xxxy[j] * fl1_fz);

                t_zzzz_xxxz[j] = fl_s_0_0 * (4.5 * pa2pb_z_x[j] * fl3_fx + 3.0 * pa2pb_zzz_x[j] * fl2_fx + 1.125 * pb_xz[j] * fl3_fx + 4.5 * pa2pb_zz_xz[j] * fl2_fx + 3.0 * pa2pb_z_xxx[j] * fl2_fx + 1.5 * pa2pb_zzzz_xz[j] * fl1_fx + 2.0 * pa2pb_zzz_xxx[j] * fl1_fx + 0.75 * pb_xxxz[j] * fl2_fx + 3.0 * pa2pb_zz_xxxz[j] * fl1_fx + pa2pb_zzzz_xxxz[j]);

                t_zzzz_xxxz[j] += fl_r_0_0 * (-9.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_zzz_x[j] * fl1_fz * fl2_fx - 2.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_xxx[j] * fl1_fx * fl1_fz * fl1_fga + 11.25 * pb_xz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_zzzz_xz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_zz_xz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_z_xxx[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_zzzz_xz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_zzz_xxx[j] * fl1_fz * fl1_fx - 3.0 * pb_xxxz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxxz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxxz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xxxz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xxxz[j] * fl1_fz);

                t_zzzz_xxyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.75 * pa_zz[j] * fl3_fx + 0.25 * pa_zzzz[j] * fl2_fx + 0.375 * pb_xx[j] * fl3_fx + 0.375 * pb_yy[j] * fl3_fx + 1.5 * pa2pb_zz_xx[j] * fl2_fx + 1.5 * pa2pb_zz_yy[j] * fl2_fx + 0.5 * pa2pb_zzzz_xx[j] * fl1_fx + 0.5 * pa2pb_zzzz_yy[j] * fl1_fx + 0.75 * pb_xxyy[j] * fl2_fx + 3.0 * pa2pb_zz_xxyy[j] * fl1_fx + pa2pb_zzzz_xxyy[j]);

                t_zzzz_xxyy[j] += fl_r_0_0 * (-0.75 * fl3_fx * fl1_fz * fl1_fgb - 0.75 * fl3_fx * fl1_fz * fl1_fga - 3.0 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx + 1.5 * fl4_fx * fl1_fz - pa_zzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa_zz[j] * fl1_fz * fl3_fx + 3.0 * pa_zzzz[j] * fl1_fz * fl2_fx - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx + 3.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.75 * pb_yy[j] * fl3_fx * fl1_fz - pa2pb_zzzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_zzzz_yy[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_zz_xx[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zz_yy[j] * fl1_fz * fl2_fx + 7.0 * pa2pb_zzzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzzz_yy[j] * fl1_fz * fl1_fx - 3.0 * pb_xxyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxyy[j] * fl1_fz * fl1_fga + 9.0 * pb_xxyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xxyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xxyy[j] * fl1_fz);

                t_zzzz_xxyz[j] = fl_s_0_0 * (1.5 * pa2pb_z_y[j] * fl3_fx + pa2pb_zzz_y[j] * fl2_fx + 0.375 * pb_yz[j] * fl3_fx + 1.5 * pa2pb_zz_yz[j] * fl2_fx + 3.0 * pa2pb_z_xxy[j] * fl2_fx + 0.5 * pa2pb_zzzz_yz[j] * fl1_fx + 2.0 * pa2pb_zzz_xxy[j] * fl1_fx + 0.75 * pb_xxyz[j] * fl2_fx + 3.0 * pa2pb_zz_xxyz[j] * fl1_fx + pa2pb_zzzz_xxyz[j]);

                t_zzzz_xxyz[j] += fl_r_0_0 * (-3.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 12.0 * pa2pb_zzz_y[j] * fl1_fz * fl2_fx - 0.75 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_xxy[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_yz[j] * fl3_fx * fl1_fz - pa2pb_zzzz_yz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_zz_yz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_z_xxy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_zzzz_yz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_zzz_xxy[j] * fl1_fz * fl1_fx - 3.0 * pb_xxyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxyz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xxyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xxyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_215_220(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (215,220)

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

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_x = pa2pbDistances.data(1156 * idx + 68);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_xxz = pa2pbDistances.data(1156 * idx + 79);

            auto pa2pb_z_xyy = pa2pbDistances.data(1156 * idx + 80);

            auto pa2pb_z_xyz = pa2pbDistances.data(1156 * idx + 81);

            auto pa2pb_z_xzz = pa2pbDistances.data(1156 * idx + 82);

            auto pa2pb_zz_xx = pa2pbDistances.data(1156 * idx + 275);

            auto pa2pb_zz_xy = pa2pbDistances.data(1156 * idx + 276);

            auto pa2pb_zz_xz = pa2pbDistances.data(1156 * idx + 277);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_zz_xxzz = pa2pbDistances.data(1156 * idx + 296);

            auto pa2pb_zz_xyyy = pa2pbDistances.data(1156 * idx + 297);

            auto pa2pb_zz_xyyz = pa2pbDistances.data(1156 * idx + 298);

            auto pa2pb_zz_xyzz = pa2pbDistances.data(1156 * idx + 299);

            auto pa2pb_zz_xzzz = pa2pbDistances.data(1156 * idx + 300);

            auto pa2pb_zzz_x = pa2pbDistances.data(1156 * idx + 612);

            auto pa2pb_zzz_z = pa2pbDistances.data(1156 * idx + 614);

            auto pa2pb_zzz_xxz = pa2pbDistances.data(1156 * idx + 623);

            auto pa2pb_zzz_xyy = pa2pbDistances.data(1156 * idx + 624);

            auto pa2pb_zzz_xyz = pa2pbDistances.data(1156 * idx + 625);

            auto pa2pb_zzz_xzz = pa2pbDistances.data(1156 * idx + 626);

            auto pa2pb_zzzz_xx = pa2pbDistances.data(1156 * idx + 1125);

            auto pa2pb_zzzz_xy = pa2pbDistances.data(1156 * idx + 1126);

            auto pa2pb_zzzz_xz = pa2pbDistances.data(1156 * idx + 1127);

            auto pa2pb_zzzz_zz = pa2pbDistances.data(1156 * idx + 1130);

            auto pa2pb_zzzz_xxzz = pa2pbDistances.data(1156 * idx + 1146);

            auto pa2pb_zzzz_xyyy = pa2pbDistances.data(1156 * idx + 1147);

            auto pa2pb_zzzz_xyyz = pa2pbDistances.data(1156 * idx + 1148);

            auto pa2pb_zzzz_xyzz = pa2pbDistances.data(1156 * idx + 1149);

            auto pa2pb_zzzz_xzzz = pa2pbDistances.data(1156 * idx + 1150);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzzz_xxzz = primBuffer.data(225 * idx + 215);

            auto t_zzzz_xyyy = primBuffer.data(225 * idx + 216);

            auto t_zzzz_xyyz = primBuffer.data(225 * idx + 217);

            auto t_zzzz_xyzz = primBuffer.data(225 * idx + 218);

            auto t_zzzz_xzzz = primBuffer.data(225 * idx + 219);

            // Batch of Integrals (215,220)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_x, pa2pb_z_xxz, pa2pb_z_xyy, pa2pb_z_xyz, \
                                     pa2pb_z_xzz, pa2pb_z_z, pa2pb_zz_xx, pa2pb_zz_xxzz, pa2pb_zz_xy, pa2pb_zz_xyyy, \
                                     pa2pb_zz_xyyz, pa2pb_zz_xyzz, pa2pb_zz_xz, pa2pb_zz_xzzz, pa2pb_zz_zz, pa2pb_zzz_x, \
                                     pa2pb_zzz_xxz, pa2pb_zzz_xyy, pa2pb_zzz_xyz, pa2pb_zzz_xzz, pa2pb_zzz_z, \
                                     pa2pb_zzzz_xx, pa2pb_zzzz_xxzz, pa2pb_zzzz_xy, pa2pb_zzzz_xyyy, pa2pb_zzzz_xyyz, \
                                     pa2pb_zzzz_xyzz, pa2pb_zzzz_xz, pa2pb_zzzz_xzzz, pa2pb_zzzz_zz, pa_zz, pa_zzzz, pb_xx, \
                                     pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, pb_xzzz, pb_zz, r_0_0, s_0_0, \
                                     t_zzzz_xxzz, t_zzzz_xyyy, t_zzzz_xyyz, t_zzzz_xyzz, t_zzzz_xzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_zzzz_xxzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_zz[j] * fl3_fx + 3.0 * pa2pb_z_z[j] * fl3_fx + 1.875 * pb_xx[j] * fl3_fx + 0.25 * pa_zzzz[j] * fl2_fx + 2.0 * pa2pb_zzz_z[j] * fl2_fx + 4.5 * pa2pb_zz_xx[j] * fl2_fx + 0.375 * pb_zz[j] * fl3_fx + 1.5 * pa2pb_zz_zz[j] * fl2_fx + 6.0 * pa2pb_z_xxz[j] * fl2_fx + 0.5 * pa2pb_zzzz_xx[j] * fl1_fx + 0.5 * pa2pb_zzzz_zz[j] * fl1_fx + 4.0 * pa2pb_zzz_xxz[j] * fl1_fx + 0.75 * pb_xxzz[j] * fl2_fx + 3.0 * pa2pb_zz_xxzz[j] * fl1_fx + pa2pb_zzzz_xxzz[j]);

                t_zzzz_xxzz[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 6.0 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 22.5 * pa_zz[j] * fl3_fx * fl1_fz - 1.5 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fga - pa_zzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.0 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 18.75 * pb_xx[j] * fl3_fx * fl1_fz + 3.0 * pa_zzzz[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_zzz_z[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_zz_xx[j] * fl2_fx * fl1_fz - 0.75 * pb_xx[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_xx[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_z_xxz[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_zz[j] * fl3_fx * fl1_fz - pa2pb_zzzz_xx[j] * fl1_fz * fl1_fgb - pa2pb_zzzz_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_zz_zz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_z_xxz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_zzzz_xx[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzzz_zz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_zzz_xxz[j] * fl1_fz * fl1_fx - 3.0 * pb_xxzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xxzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xxzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xxzz[j] * fl1_fz);

                t_zzzz_xyyy[j] = fl_s_0_0 * (1.125 * pb_xy[j] * fl3_fx + 4.5 * pa2pb_zz_xy[j] * fl2_fx + 1.5 * pa2pb_zzzz_xy[j] * fl1_fx + 0.75 * pb_xyyy[j] * fl2_fx + 3.0 * pa2pb_zz_xyyy[j] * fl1_fx + pa2pb_zzzz_xyyy[j]);

                t_zzzz_xyyy[j] += fl_r_0_0 * (-2.25 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx + 11.25 * pb_xy[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_zzzz_xy[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_zz_xy[j] * fl1_fz * fl2_fx + 21.0 * pa2pb_zzzz_xy[j] * fl1_fz * fl1_fx - 3.0 * pb_xyyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xyyy[j] * fl1_fz * fl1_fga + 9.0 * pb_xyyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xyyy[j] * fl1_fz);

                t_zzzz_xyyz[j] = fl_s_0_0 * (1.5 * pa2pb_z_x[j] * fl3_fx + pa2pb_zzz_x[j] * fl2_fx + 0.375 * pb_xz[j] * fl3_fx + 1.5 * pa2pb_zz_xz[j] * fl2_fx + 3.0 * pa2pb_z_xyy[j] * fl2_fx + 0.5 * pa2pb_zzzz_xz[j] * fl1_fx + 2.0 * pa2pb_zzz_xyy[j] * fl1_fx + 0.75 * pb_xyyz[j] * fl2_fx + 3.0 * pa2pb_zz_xyyz[j] * fl1_fx + pa2pb_zzzz_xyyz[j]);

                t_zzzz_xyyz[j] += fl_r_0_0 * (-3.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_x[j] * fl3_fx * fl1_fz + 12.0 * pa2pb_zzz_x[j] * fl1_fz * fl2_fx - 0.75 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_xyy[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_xz[j] * fl3_fx * fl1_fz - pa2pb_zzzz_xz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_zz_xz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_z_xyy[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_zzzz_xz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_zzz_xyy[j] * fl1_fz * fl1_fx - 3.0 * pb_xyyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xyyz[j] * fl1_fz * fl1_fga + 9.0 * pb_xyyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xyyz[j] * fl1_fz);

                t_zzzz_xyzz[j] = fl_s_0_0 * (1.875 * pb_xy[j] * fl3_fx + 4.5 * pa2pb_zz_xy[j] * fl2_fx + 6.0 * pa2pb_z_xyz[j] * fl2_fx + 0.5 * pa2pb_zzzz_xy[j] * fl1_fx + 4.0 * pa2pb_zzz_xyz[j] * fl1_fx + 0.75 * pb_xyzz[j] * fl2_fx + 3.0 * pa2pb_zz_xyzz[j] * fl1_fx + pa2pb_zzzz_xyzz[j]);

                t_zzzz_xyzz[j] += fl_r_0_0 * (-4.5 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fga + 18.75 * pb_xy[j] * fl3_fx * fl1_fz + 54.0 * pa2pb_zz_xy[j] * fl2_fx * fl1_fz - 0.75 * pb_xy[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_xy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_xy[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_z_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzzz_xy[j] * fl1_fz * fl1_fgb + 72.0 * pa2pb_z_xyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_zzzz_xy[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_zzz_xyz[j] * fl1_fz * fl1_fx - 3.0 * pb_xyzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xyzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xyzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xyzz[j] * fl1_fz);

                t_zzzz_xzzz[j] = fl_s_0_0 * (7.5 * pa2pb_z_x[j] * fl3_fx + 5.625 * pb_xz[j] * fl3_fx + 3.0 * pa2pb_zzz_x[j] * fl2_fx + 13.5 * pa2pb_zz_xz[j] * fl2_fx + 9.0 * pa2pb_z_xzz[j] * fl2_fx + 1.5 * pa2pb_zzzz_xz[j] * fl1_fx + 6.0 * pa2pb_zzz_xzz[j] * fl1_fx + 0.75 * pb_xzzz[j] * fl2_fx + 3.0 * pa2pb_zz_xzzz[j] * fl1_fx + pa2pb_zzzz_xzzz[j]);

                t_zzzz_xzzz[j] += fl_r_0_0 * (75.0 * pa2pb_z_x[j] * fl3_fx * fl1_fz - 9.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_x[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zzz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 56.25 * pb_xz[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_zzz_x[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_zz_xz[j] * fl2_fx * fl1_fz - 2.25 * pb_xz[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_xz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_xz[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_z_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zzzz_xz[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_z_xzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_zzzz_xz[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_zzz_xzz[j] * fl1_fz * fl1_fx - 3.0 * pb_xzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_xzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_xzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_xzzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGG_220_225(      CMemBlock2D<double>& primBuffer,
                                   const CMemBlock2D<double>& auxBuffer,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pbDistances,
                                   const CMemBlock2D<double>& pa2pbDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (220,225)

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

            auto pa_zz = paDistances.data(34 * idx + 8);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_y = pa2pbDistances.data(1156 * idx + 69);

            auto pa2pb_z_z = pa2pbDistances.data(1156 * idx + 70);

            auto pa2pb_z_yyy = pa2pbDistances.data(1156 * idx + 83);

            auto pa2pb_z_yyz = pa2pbDistances.data(1156 * idx + 84);

            auto pa2pb_z_yzz = pa2pbDistances.data(1156 * idx + 85);

            auto pa2pb_z_zzz = pa2pbDistances.data(1156 * idx + 86);

            auto pa2pb_zz_yy = pa2pbDistances.data(1156 * idx + 278);

            auto pa2pb_zz_yz = pa2pbDistances.data(1156 * idx + 279);

            auto pa2pb_zz_zz = pa2pbDistances.data(1156 * idx + 280);

            auto pa2pb_zz_yyyy = pa2pbDistances.data(1156 * idx + 301);

            auto pa2pb_zz_yyyz = pa2pbDistances.data(1156 * idx + 302);

            auto pa2pb_zz_yyzz = pa2pbDistances.data(1156 * idx + 303);

            auto pa2pb_zz_yzzz = pa2pbDistances.data(1156 * idx + 304);

            auto pa2pb_zz_zzzz = pa2pbDistances.data(1156 * idx + 305);

            auto pa2pb_zzz_y = pa2pbDistances.data(1156 * idx + 613);

            auto pa2pb_zzz_z = pa2pbDistances.data(1156 * idx + 614);

            auto pa2pb_zzz_yyy = pa2pbDistances.data(1156 * idx + 627);

            auto pa2pb_zzz_yyz = pa2pbDistances.data(1156 * idx + 628);

            auto pa2pb_zzz_yzz = pa2pbDistances.data(1156 * idx + 629);

            auto pa2pb_zzz_zzz = pa2pbDistances.data(1156 * idx + 630);

            auto pa2pb_zzzz_yy = pa2pbDistances.data(1156 * idx + 1128);

            auto pa2pb_zzzz_yz = pa2pbDistances.data(1156 * idx + 1129);

            auto pa2pb_zzzz_zz = pa2pbDistances.data(1156 * idx + 1130);

            auto pa2pb_zzzz_yyyy = pa2pbDistances.data(1156 * idx + 1151);

            auto pa2pb_zzzz_yyyz = pa2pbDistances.data(1156 * idx + 1152);

            auto pa2pb_zzzz_yyzz = pa2pbDistances.data(1156 * idx + 1153);

            auto pa2pb_zzzz_yzzz = pa2pbDistances.data(1156 * idx + 1154);

            auto pa2pb_zzzz_zzzz = pa2pbDistances.data(1156 * idx + 1155);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzzz_yyyy = primBuffer.data(225 * idx + 220);

            auto t_zzzz_yyyz = primBuffer.data(225 * idx + 221);

            auto t_zzzz_yyzz = primBuffer.data(225 * idx + 222);

            auto t_zzzz_yzzz = primBuffer.data(225 * idx + 223);

            auto t_zzzz_zzzz = primBuffer.data(225 * idx + 224);

            // Batch of Integrals (220,225)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_y, pa2pb_z_yyy, pa2pb_z_yyz, pa2pb_z_yzz, \
                                     pa2pb_z_z, pa2pb_z_zzz, pa2pb_zz_yy, pa2pb_zz_yyyy, pa2pb_zz_yyyz, \
                                     pa2pb_zz_yyzz, pa2pb_zz_yz, pa2pb_zz_yzzz, pa2pb_zz_zz, pa2pb_zz_zzzz, pa2pb_zzz_y, \
                                     pa2pb_zzz_yyy, pa2pb_zzz_yyz, pa2pb_zzz_yzz, pa2pb_zzz_z, pa2pb_zzz_zzz, \
                                     pa2pb_zzzz_yy, pa2pb_zzzz_yyyy, pa2pb_zzzz_yyyz, pa2pb_zzzz_yyzz, pa2pb_zzzz_yz, \
                                     pa2pb_zzzz_yzzz, pa2pb_zzzz_zz, pa2pb_zzzz_zzzz, pa_zz, pa_zzzz, pb_yy, pb_yyyy, pb_yyyz, \
                                     pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, r_0_0, s_0_0, t_zzzz_yyyy, t_zzzz_yyyz, \
                                     t_zzzz_yyzz, t_zzzz_yzzz, t_zzzz_zzzz: VLX_ALIGN)
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

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_zzzz_yyyy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_zz[j] * fl3_fx + 0.75 * pa_zzzz[j] * fl2_fx + 2.25 * pb_yy[j] * fl3_fx + 9.0 * pa2pb_zz_yy[j] * fl2_fx + 3.0 * pa2pb_zzzz_yy[j] * fl1_fx + 0.75 * pb_yyyy[j] * fl2_fx + 3.0 * pa2pb_zz_yyyy[j] * fl1_fx + pa2pb_zzzz_yyyy[j]);

                t_zzzz_yyyy[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 9.0 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx + 4.5 * fl4_fx * fl1_fz - 3.0 * pa_zzzz[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_zz[j] * fl1_fz * fl3_fx + 9.0 * pa_zzzz[j] * fl1_fz * fl2_fx - 4.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - 18.0 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx + 22.5 * pb_yy[j] * fl3_fx * fl1_fz - 6.0 * pa2pb_zzzz_yy[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_zz_yy[j] * fl1_fz * fl2_fx + 42.0 * pa2pb_zzzz_yy[j] * fl1_fz * fl1_fx - 3.0 * pb_yyyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_yyyy[j] * fl1_fz * fl1_fga + 9.0 * pb_yyyy[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_yyyy[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_yyyy[j] * fl1_fz);

                t_zzzz_yyyz[j] = fl_s_0_0 * (4.5 * pa2pb_z_y[j] * fl3_fx + 3.0 * pa2pb_zzz_y[j] * fl2_fx + 1.125 * pb_yz[j] * fl3_fx + 4.5 * pa2pb_zz_yz[j] * fl2_fx + 3.0 * pa2pb_z_yyy[j] * fl2_fx + 1.5 * pa2pb_zzzz_yz[j] * fl1_fx + 2.0 * pa2pb_zzz_yyy[j] * fl1_fx + 0.75 * pb_yyyz[j] * fl2_fx + 3.0 * pa2pb_zz_yyyz[j] * fl1_fx + pa2pb_zzzz_yyyz[j]);

                t_zzzz_yyyz[j] += fl_r_0_0 * (-9.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pa2pb_z_y[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_zzz_y[j] * fl1_fz * fl2_fx - 2.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_yyy[j] * fl1_fx * fl1_fz * fl1_fga + 11.25 * pb_yz[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_zzzz_yz[j] * fl1_fz * fl1_fgb + 54.0 * pa2pb_zz_yz[j] * fl1_fz * fl2_fx + 36.0 * pa2pb_z_yyy[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_zzzz_yz[j] * fl1_fz * fl1_fx + 28.0 * pa2pb_zzz_yyy[j] * fl1_fz * fl1_fx - 3.0 * pb_yyyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_yyyz[j] * fl1_fz * fl1_fga + 9.0 * pb_yyyz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_yyyz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_yyyz[j] * fl1_fz);

                t_zzzz_yyzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_zz[j] * fl3_fx + 3.0 * pa2pb_z_z[j] * fl3_fx + 1.875 * pb_yy[j] * fl3_fx + 0.25 * pa_zzzz[j] * fl2_fx + 2.0 * pa2pb_zzz_z[j] * fl2_fx + 4.5 * pa2pb_zz_yy[j] * fl2_fx + 0.375 * pb_zz[j] * fl3_fx + 1.5 * pa2pb_zz_zz[j] * fl2_fx + 6.0 * pa2pb_z_yyz[j] * fl2_fx + 0.5 * pa2pb_zzzz_yy[j] * fl1_fx + 0.5 * pa2pb_zzzz_zz[j] * fl1_fx + 4.0 * pa2pb_zzz_yyz[j] * fl1_fx + 0.75 * pb_yyzz[j] * fl2_fx + 3.0 * pa2pb_zz_yyzz[j] * fl1_fx + pa2pb_zzzz_yyzz[j]);

                t_zzzz_yyzz[j] += fl_r_0_0 * (-2.25 * fl3_fx * fl1_fz * fl1_fgb - 2.25 * fl3_fx * fl1_fz * fl1_fga - 6.0 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 7.5 * fl4_fx * fl1_fz + 22.5 * pa_zz[j] * fl3_fx * fl1_fz - 1.5 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 6.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fga - pa_zzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 4.0 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 30.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz + 18.75 * pb_yy[j] * fl3_fx * fl1_fz + 3.0 * pa_zzzz[j] * fl1_fz * fl2_fx + 24.0 * pa2pb_zzz_z[j] * fl1_fz * fl2_fx + 54.0 * pa2pb_zz_yy[j] * fl2_fx * fl1_fz - 0.75 * pb_yy[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_yy[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_z_yyz[j] * fl1_fx * fl1_fz * fl1_fga + 3.75 * pb_zz[j] * fl3_fx * fl1_fz - pa2pb_zzzz_yy[j] * fl1_fz * fl1_fgb - pa2pb_zzzz_zz[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_zz_zz[j] * fl1_fz * fl2_fx + 72.0 * pa2pb_z_yyz[j] * fl2_fx * fl1_fz + 7.0 * pa2pb_zzzz_yy[j] * fl1_fz * fl1_fx + 7.0 * pa2pb_zzzz_zz[j] * fl1_fz * fl1_fx + 56.0 * pa2pb_zzz_yyz[j] * fl1_fz * fl1_fx - 3.0 * pb_yyzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_yyzz[j] * fl1_fz * fl1_fga + 9.0 * pb_yyzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_yyzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_yyzz[j] * fl1_fz);

                t_zzzz_yzzz[j] = fl_s_0_0 * (7.5 * pa2pb_z_y[j] * fl3_fx + 5.625 * pb_yz[j] * fl3_fx + 3.0 * pa2pb_zzz_y[j] * fl2_fx + 13.5 * pa2pb_zz_yz[j] * fl2_fx + 9.0 * pa2pb_z_yzz[j] * fl2_fx + 1.5 * pa2pb_zzzz_yz[j] * fl1_fx + 6.0 * pa2pb_zzz_yzz[j] * fl1_fx + 0.75 * pb_yzzz[j] * fl2_fx + 3.0 * pa2pb_zz_yzzz[j] * fl1_fx + pa2pb_zzzz_yzzz[j]);

                t_zzzz_yzzz[j] += fl_r_0_0 * (75.0 * pa2pb_z_y[j] * fl3_fx * fl1_fz - 9.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_z_y[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zzz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 56.25 * pb_yz[j] * fl3_fx * fl1_fz + 36.0 * pa2pb_zzz_y[j] * fl1_fz * fl2_fx + 162.0 * pa2pb_zz_yz[j] * fl2_fx * fl1_fz - 2.25 * pb_yz[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_yz[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_yz[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_z_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zzzz_yz[j] * fl1_fz * fl1_fgb + 108.0 * pa2pb_z_yzz[j] * fl2_fx * fl1_fz + 21.0 * pa2pb_zzzz_yz[j] * fl1_fz * fl1_fx + 84.0 * pa2pb_zzz_yzz[j] * fl1_fz * fl1_fx - 3.0 * pb_yzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_yzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_yzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_yzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_yzzz[j] * fl1_fz);

                t_zzzz_zzzz[j] = fl_s_0_0 * (6.5625 * fl4_fx + 11.25 * pa_zz[j] * fl3_fx + 30.0 * pa2pb_z_z[j] * fl3_fx + 11.25 * pb_zz[j] * fl3_fx + 0.75 * pa_zzzz[j] * fl2_fx + 12.0 * pa2pb_zzz_z[j] * fl2_fx + 27.0 * pa2pb_zz_zz[j] * fl2_fx + 12.0 * pa2pb_z_zzz[j] * fl2_fx + 3.0 * pa2pb_zzzz_zz[j] * fl1_fx + 8.0 * pa2pb_zzz_zzz[j] * fl1_fx + 0.75 * pb_zzzz[j] * fl2_fx + 3.0 * pa2pb_zz_zzzz[j] * fl1_fx + pa2pb_zzzz_zzzz[j]);

                t_zzzz_zzzz[j] += fl_r_0_0 * (52.5 * fl4_fx * fl1_fz - 11.25 * fl3_fx * fl1_fz * fl1_fgb - 11.25 * fl3_fx * fl1_fz * fl1_fga - 27.0 * pa_zz[j] * fl2_fx * fl1_fz * fl1_fgb + 112.5 * pa_zz[j] * fl3_fx * fl1_fz + 300.0 * pa2pb_z_z[j] * fl3_fx * fl1_fz - 4.5 * pa_zz[j] * fl1_fz * fl1_fga * fl2_fx - 36.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fgb - 36.0 * pa2pb_z_z[j] * fl2_fx * fl1_fz * fl1_fga - 27.0 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_zzzz[j] * fl1_fx * fl1_fz * fl1_fgb - 24.0 * pa2pb_zzz_z[j] * fl1_fx * fl1_fz * fl1_fgb + 112.5 * pb_zz[j] * fl3_fx * fl1_fz + 9.0 * pa_zzzz[j] * fl1_fz * fl2_fx + 144.0 * pa2pb_zzz_z[j] * fl1_fz * fl2_fx + 324.0 * pa2pb_zz_zz[j] * fl2_fx * fl1_fz - 4.5 * pb_zz[j] * fl2_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_zz_zz[j] * fl1_fx * fl1_fz * fl1_fgb - 18.0 * pa2pb_zz_zz[j] * fl1_fz * fl1_fga * fl1_fx - 24.0 * pa2pb_z_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zzzz_zz[j] * fl1_fz * fl1_fgb + 144.0 * pa2pb_z_zzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zzzz_zz[j] * fl1_fz * fl1_fx + 112.0 * pa2pb_zzz_zzz[j] * fl1_fz * fl1_fx - 3.0 * pb_zzzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_zzzz[j] * fl1_fz * fl1_fga + 9.0 * pb_zzzz[j] * fl2_fx * fl1_fz + 42.0 * pa2pb_zz_zzzz[j] * fl1_fz * fl1_fx + 16.0 * pa2pb_zzzz_zzzz[j] * fl1_fz);
            }

            idx++;
        }
    }


} // kinrecfunc namespace

