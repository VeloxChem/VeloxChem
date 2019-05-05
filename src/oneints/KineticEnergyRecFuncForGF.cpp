//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForGF.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForGF(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CMemBlock2D<double>& pa2pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForGF_0_5(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                               braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_5_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_10_15(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_15_20(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_20_25(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_25_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_30_35(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_35_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_40_45(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_45_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_50_55(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_55_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_60_65(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_65_70(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_70_75(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_75_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_80_85(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_85_90(primBuffer, auxBuffer, osFactors, paDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_90_95(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                 braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_95_100(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                  braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_100_105(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_105_110(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_110_115(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_115_120(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_120_125(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_125_130(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_130_135(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_135_140(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_140_145(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 

        kinrecfunc::compKineticEnergyForGF_145_150(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, pa2pbDistances, 
                                                   braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compKineticEnergyForGF_0_5(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_xxx = paDistances.data(34 * idx + 9);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xx_xxx = pa2pbDistances.data(646 * idx + 66);

            auto pa2pb_xx_xxy = pa2pbDistances.data(646 * idx + 67);

            auto pa2pb_xx_xxz = pa2pbDistances.data(646 * idx + 68);

            auto pa2pb_xx_xyy = pa2pbDistances.data(646 * idx + 69);

            auto pa2pb_xx_xyz = pa2pbDistances.data(646 * idx + 70);

            auto pa2pb_xxx_xx = pa2pbDistances.data(646 * idx + 174);

            auto pa2pb_xxx_xy = pa2pbDistances.data(646 * idx + 175);

            auto pa2pb_xxx_xz = pa2pbDistances.data(646 * idx + 176);

            auto pa2pb_xxx_yy = pa2pbDistances.data(646 * idx + 177);

            auto pa2pb_xxx_yz = pa2pbDistances.data(646 * idx + 178);

            auto pa2pb_xxxx_x = pa2pbDistances.data(646 * idx + 361);

            auto pa2pb_xxxx_y = pa2pbDistances.data(646 * idx + 362);

            auto pa2pb_xxxx_z = pa2pbDistances.data(646 * idx + 363);

            auto pa2pb_xxxx_xxx = pa2pbDistances.data(646 * idx + 370);

            auto pa2pb_xxxx_xxy = pa2pbDistances.data(646 * idx + 371);

            auto pa2pb_xxxx_xxz = pa2pbDistances.data(646 * idx + 372);

            auto pa2pb_xxxx_xyy = pa2pbDistances.data(646 * idx + 373);

            auto pa2pb_xxxx_xyz = pa2pbDistances.data(646 * idx + 374);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxx_xxx = primBuffer.data(150 * idx);

            auto t_xxxx_xxy = primBuffer.data(150 * idx + 1);

            auto t_xxxx_xxz = primBuffer.data(150 * idx + 2);

            auto t_xxxx_xyy = primBuffer.data(150 * idx + 3);

            auto t_xxxx_xyz = primBuffer.data(150 * idx + 4);

            // Batch of Integrals (0,5)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_xx_x, pa2pb_xx_xxx, pa2pb_xx_xxy, pa2pb_xx_xxz, pa2pb_xx_xyy, pa2pb_xx_xyz, \
                                     pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxx_xx, pa2pb_xxx_xy, pa2pb_xxx_xz, pa2pb_xxx_yy, \
                                     pa2pb_xxx_yz, pa2pb_xxxx_x, pa2pb_xxxx_xxx, pa2pb_xxxx_xxy, pa2pb_xxxx_xxz, \
                                     pa2pb_xxxx_xyy, pa2pb_xxxx_xyz, pa2pb_xxxx_y, pa2pb_xxxx_z, pa_x, pa_xxx, pb_x, pb_xxx, \
                                     pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_y, pb_z, r_0_0, s_0_0, t_xxxx_xxx, t_xxxx_xxy, \
                                     t_xxxx_xxz, t_xxxx_xyy, t_xxxx_xyz: VLX_ALIGN)
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

                t_xxxx_xxx[j] = fl_s_0_0 * (7.5 * pa_x[j] * fl3_fx + 5.625 * pb_x[j] * fl3_fx + 3.0 * pa_xxx[j] * fl2_fx + 13.5 * pa2pb_xx_x[j] * fl2_fx + 9.0 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pa2pb_xxxx_x[j] * fl1_fx + 6.0 * pa2pb_xxx_xx[j] * fl1_fx + 0.75 * pb_xxx[j] * fl2_fx + 3.0 * pa2pb_xx_xxx[j] * fl1_fx + pa2pb_xxxx_xxx[j]);

                t_xxxx_xxx[j] += fl_r_0_0 * (60.0 * pa_x[j] * fl3_fx * fl1_fz - 9.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pb_x[j] * fl3_fx * fl1_fz + 30.0 * pa_xxx[j] * fl1_fz * fl2_fx + 135.0 * pa2pb_xx_x[j] * fl2_fx * fl1_fz - 2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxx_x[j] * fl1_fz * fl1_fgb + 90.0 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxxx_x[j] * fl1_fz * fl1_fx + 72.0 * pa2pb_xxx_xx[j] * fl1_fz * fl1_fx - 3.0 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxx[j] * fl1_fz * fl1_fga + 7.5 * pb_xxx[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_xxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_xxx[j] * fl1_fz);

                t_xxxx_xxy[j] = fl_s_0_0 * (1.875 * pb_y[j] * fl3_fx + 4.5 * pa2pb_xx_y[j] * fl2_fx + 6.0 * pa2pb_x_xy[j] * fl2_fx + 0.5 * pa2pb_xxxx_y[j] * fl1_fx + 4.0 * pa2pb_xxx_xy[j] * fl1_fx + 0.75 * pb_xxy[j] * fl2_fx + 3.0 * pa2pb_xx_xxy[j] * fl1_fx + pa2pb_xxxx_xxy[j]);

                t_xxxx_xxy[j] += fl_r_0_0 * (-4.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga + 15.0 * pb_y[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_xx_y[j] * fl2_fx * fl1_fz - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxx_y[j] * fl1_fz * fl1_fgb + 60.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxx_y[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_xxx_xy[j] * fl1_fz * fl1_fx - 3.0 * pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxy[j] * fl1_fz * fl1_fga + 7.5 * pb_xxy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_xxy[j] * fl1_fz);

                t_xxxx_xxz[j] = fl_s_0_0 * (1.875 * pb_z[j] * fl3_fx + 4.5 * pa2pb_xx_z[j] * fl2_fx + 6.0 * pa2pb_x_xz[j] * fl2_fx + 0.5 * pa2pb_xxxx_z[j] * fl1_fx + 4.0 * pa2pb_xxx_xz[j] * fl1_fx + 0.75 * pb_xxz[j] * fl2_fx + 3.0 * pa2pb_xx_xxz[j] * fl1_fx + pa2pb_xxxx_xxz[j]);

                t_xxxx_xxz[j] += fl_r_0_0 * (-4.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga + 15.0 * pb_z[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_xx_z[j] * fl2_fx * fl1_fz - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxx_z[j] * fl1_fz * fl1_fgb + 60.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxx_z[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_xxx_xz[j] * fl1_fz * fl1_fx - 3.0 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xxz[j] * fl1_fz * fl1_fga + 7.5 * pb_xxz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_xxz[j] * fl1_fz);

                t_xxxx_xyy[j] = fl_s_0_0 * (1.5 * pa_x[j] * fl3_fx + pa_xxx[j] * fl2_fx + 0.375 * pb_x[j] * fl3_fx + 1.5 * pa2pb_xx_x[j] * fl2_fx + 3.0 * pa2pb_x_yy[j] * fl2_fx + 0.5 * pa2pb_xxxx_x[j] * fl1_fx + 2.0 * pa2pb_xxx_yy[j] * fl1_fx + 0.75 * pb_xyy[j] * fl2_fx + 3.0 * pa2pb_xx_xyy[j] * fl1_fx + pa2pb_xxxx_xyy[j]);

                t_xxxx_xyy[j] += fl_r_0_0 * (-3.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_x[j] * fl3_fx * fl1_fz + 10.0 * pa_xxx[j] * fl1_fz * fl2_fx - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fga + 3.0 * pb_x[j] * fl3_fx * fl1_fz - pa2pb_xxxx_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxx_x[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xxx_yy[j] * fl1_fz * fl1_fx - 3.0 * pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xyy[j] * fl1_fz * fl1_fga + 7.5 * pb_xyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_xyy[j] * fl1_fz);

                t_xxxx_xyz[j] = fl_s_0_0 * (3.0 * pa2pb_x_yz[j] * fl2_fx + 2.0 * pa2pb_xxx_yz[j] * fl1_fx + 0.75 * pb_xyz[j] * fl2_fx + 3.0 * pa2pb_xx_xyz[j] * fl1_fx + pa2pb_xxxx_xyz[j]);

                t_xxxx_xyz[j] += fl_r_0_0 * (-6.0 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fga + 30.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_xxx_yz[j] * fl1_fz * fl1_fx - 3.0 * pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xyz[j] * fl1_fz * fl1_fga + 7.5 * pb_xyz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_5_10(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_xxx = paDistances.data(34 * idx + 9);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xx_xzz = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_xx_yyy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_xx_yyz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_xx_yzz = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_xx_zzz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_xxx_zz = pa2pbDistances.data(646 * idx + 179);

            auto pa2pb_xxxx_x = pa2pbDistances.data(646 * idx + 361);

            auto pa2pb_xxxx_y = pa2pbDistances.data(646 * idx + 362);

            auto pa2pb_xxxx_z = pa2pbDistances.data(646 * idx + 363);

            auto pa2pb_xxxx_xzz = pa2pbDistances.data(646 * idx + 375);

            auto pa2pb_xxxx_yyy = pa2pbDistances.data(646 * idx + 376);

            auto pa2pb_xxxx_yyz = pa2pbDistances.data(646 * idx + 377);

            auto pa2pb_xxxx_yzz = pa2pbDistances.data(646 * idx + 378);

            auto pa2pb_xxxx_zzz = pa2pbDistances.data(646 * idx + 379);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxx_xzz = primBuffer.data(150 * idx + 5);

            auto t_xxxx_yyy = primBuffer.data(150 * idx + 6);

            auto t_xxxx_yyz = primBuffer.data(150 * idx + 7);

            auto t_xxxx_yzz = primBuffer.data(150 * idx + 8);

            auto t_xxxx_zzz = primBuffer.data(150 * idx + 9);

            // Batch of Integrals (5,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_zz, pa2pb_xx_x, pa2pb_xx_xzz, pa2pb_xx_y, \
                                     pa2pb_xx_yyy, pa2pb_xx_yyz, pa2pb_xx_yzz, pa2pb_xx_z, pa2pb_xx_zzz, pa2pb_xxx_zz, \
                                     pa2pb_xxxx_x, pa2pb_xxxx_xzz, pa2pb_xxxx_y, pa2pb_xxxx_yyy, pa2pb_xxxx_yyz, \
                                     pa2pb_xxxx_yzz, pa2pb_xxxx_z, pa2pb_xxxx_zzz, pa_x, pa_xxx, pb_x, pb_xzz, pb_y, pb_yyy, \
                                     pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, s_0_0, t_xxxx_xzz, t_xxxx_yyy, t_xxxx_yyz, \
                                     t_xxxx_yzz, t_xxxx_zzz: VLX_ALIGN)
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

                t_xxxx_xzz[j] = fl_s_0_0 * (1.5 * pa_x[j] * fl3_fx + pa_xxx[j] * fl2_fx + 0.375 * pb_x[j] * fl3_fx + 1.5 * pa2pb_xx_x[j] * fl2_fx + 3.0 * pa2pb_x_zz[j] * fl2_fx + 0.5 * pa2pb_xxxx_x[j] * fl1_fx + 2.0 * pa2pb_xxx_zz[j] * fl1_fx + 0.75 * pb_xzz[j] * fl2_fx + 3.0 * pa2pb_xx_xzz[j] * fl1_fx + pa2pb_xxxx_xzz[j]);

                t_xxxx_xzz[j] += fl_r_0_0 * (-3.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_x[j] * fl3_fx * fl1_fz + 10.0 * pa_xxx[j] * fl1_fz * fl2_fx - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fga + 3.0 * pb_x[j] * fl3_fx * fl1_fz - pa2pb_xxxx_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxx_x[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xxx_zz[j] * fl1_fz * fl1_fx - 3.0 * pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_xzz[j] * fl1_fz * fl1_fga + 7.5 * pb_xzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_xzz[j] * fl1_fz);

                t_xxxx_yyy[j] = fl_s_0_0 * (1.125 * pb_y[j] * fl3_fx + 4.5 * pa2pb_xx_y[j] * fl2_fx + 1.5 * pa2pb_xxxx_y[j] * fl1_fx + 0.75 * pb_yyy[j] * fl2_fx + 3.0 * pa2pb_xx_yyy[j] * fl1_fx + pa2pb_xxxx_yyy[j]);

                t_xxxx_yyy[j] += fl_r_0_0 * (-2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx + 9.0 * pb_y[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_xxxx_y[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xxxx_y[j] * fl1_fz * fl1_fx - 3.0 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_yyy[j] * fl1_fz * fl1_fga + 7.5 * pb_yyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_yyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_yyy[j] * fl1_fz);

                t_xxxx_yyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 1.5 * pa2pb_xx_z[j] * fl2_fx + 0.5 * pa2pb_xxxx_z[j] * fl1_fx + 0.75 * pb_yyz[j] * fl2_fx + 3.0 * pa2pb_xx_yyz[j] * fl1_fx + pa2pb_xxxx_yyz[j]);

                t_xxxx_yyz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx + 3.0 * pb_z[j] * fl3_fx * fl1_fz - pa2pb_xxxx_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxxx_z[j] * fl1_fz * fl1_fx - 3.0 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_yyz[j] * fl1_fz * fl1_fga + 7.5 * pb_yyz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_yyz[j] * fl1_fz);

                t_xxxx_yzz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 1.5 * pa2pb_xx_y[j] * fl2_fx + 0.5 * pa2pb_xxxx_y[j] * fl1_fx + 0.75 * pb_yzz[j] * fl2_fx + 3.0 * pa2pb_xx_yzz[j] * fl1_fx + pa2pb_xxxx_yzz[j]);

                t_xxxx_yzz[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx + 3.0 * pb_y[j] * fl3_fx * fl1_fz - pa2pb_xxxx_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_xxxx_y[j] * fl1_fz * fl1_fx - 3.0 * pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_yzz[j] * fl1_fz * fl1_fga + 7.5 * pb_yzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_yzz[j] * fl1_fz);

                t_xxxx_zzz[j] = fl_s_0_0 * (1.125 * pb_z[j] * fl3_fx + 4.5 * pa2pb_xx_z[j] * fl2_fx + 1.5 * pa2pb_xxxx_z[j] * fl1_fx + 0.75 * pb_zzz[j] * fl2_fx + 3.0 * pa2pb_xx_zzz[j] * fl1_fx + pa2pb_xxxx_zzz[j]);

                t_xxxx_zzz[j] += fl_r_0_0 * (-2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx + 9.0 * pb_z[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_xxxx_z[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xxxx_z[j] * fl1_fz * fl1_fx - 3.0 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_xx_zzz[j] * fl1_fz * fl1_fga + 7.5 * pb_zzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_xx_zzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xxxx_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_10_15(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_xy_xxx = pa2pbDistances.data(646 * idx + 85);

            auto pa2pb_xy_xxy = pa2pbDistances.data(646 * idx + 86);

            auto pa2pb_xy_xxz = pa2pbDistances.data(646 * idx + 87);

            auto pa2pb_xy_xyy = pa2pbDistances.data(646 * idx + 88);

            auto pa2pb_xy_xyz = pa2pbDistances.data(646 * idx + 89);

            auto pa2pb_xxx_xx = pa2pbDistances.data(646 * idx + 174);

            auto pa2pb_xxx_xy = pa2pbDistances.data(646 * idx + 175);

            auto pa2pb_xxx_xz = pa2pbDistances.data(646 * idx + 176);

            auto pa2pb_xxy_xx = pa2pbDistances.data(646 * idx + 193);

            auto pa2pb_xxy_xy = pa2pbDistances.data(646 * idx + 194);

            auto pa2pb_xxy_xz = pa2pbDistances.data(646 * idx + 195);

            auto pa2pb_xxy_yy = pa2pbDistances.data(646 * idx + 196);

            auto pa2pb_xxy_yz = pa2pbDistances.data(646 * idx + 197);

            auto pa2pb_xxxy_x = pa2pbDistances.data(646 * idx + 380);

            auto pa2pb_xxxy_y = pa2pbDistances.data(646 * idx + 381);

            auto pa2pb_xxxy_z = pa2pbDistances.data(646 * idx + 382);

            auto pa2pb_xxxy_xxx = pa2pbDistances.data(646 * idx + 389);

            auto pa2pb_xxxy_xxy = pa2pbDistances.data(646 * idx + 390);

            auto pa2pb_xxxy_xxz = pa2pbDistances.data(646 * idx + 391);

            auto pa2pb_xxxy_xyy = pa2pbDistances.data(646 * idx + 392);

            auto pa2pb_xxxy_xyz = pa2pbDistances.data(646 * idx + 393);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxy_xxx = primBuffer.data(150 * idx + 10);

            auto t_xxxy_xxy = primBuffer.data(150 * idx + 11);

            auto t_xxxy_xxz = primBuffer.data(150 * idx + 12);

            auto t_xxxy_xyy = primBuffer.data(150 * idx + 13);

            auto t_xxxy_xyz = primBuffer.data(150 * idx + 14);

            // Batch of Integrals (10,15)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_xx_x, pa2pb_xx_y, \
                                     pa2pb_xx_z, pa2pb_xxx_xx, pa2pb_xxx_xy, pa2pb_xxx_xz, pa2pb_xxxy_x, \
                                     pa2pb_xxxy_xxx, pa2pb_xxxy_xxy, pa2pb_xxxy_xxz, pa2pb_xxxy_xyy, pa2pb_xxxy_xyz, \
                                     pa2pb_xxxy_y, pa2pb_xxxy_z, pa2pb_xxy_xx, pa2pb_xxy_xy, pa2pb_xxy_xz, pa2pb_xxy_yy, \
                                     pa2pb_xxy_yz, pa2pb_xy_x, pa2pb_xy_xxx, pa2pb_xy_xxy, pa2pb_xy_xxz, pa2pb_xy_xyy, \
                                     pa2pb_xy_xyz, pa2pb_xy_y, pa2pb_xy_z, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, \
                                     pa2pb_y_yz, pa_x, pa_xxx, pa_xxy, pa_y, pb_x, pb_y, pb_z, r_0_0, s_0_0, t_xxxy_xxx, \
                                     t_xxxy_xxy, t_xxxy_xxz, t_xxxy_xyy, t_xxxy_xyz: VLX_ALIGN)
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

                t_xxxy_xxx[j] = fl_s_0_0 * (1.875 * pa_y[j] * fl3_fx + 2.25 * pa_xxy[j] * fl2_fx + 6.75 * pa2pb_xy_x[j] * fl2_fx + 2.25 * pa2pb_y_xx[j] * fl2_fx + 1.5 * pa2pb_xxxy_x[j] * fl1_fx + 4.5 * pa2pb_xxy_xx[j] * fl1_fx + 1.5 * pa2pb_xy_xxx[j] * fl1_fx + pa2pb_xxxy_xxx[j]);

                t_xxxy_xxx[j] += fl_r_0_0 * (15.0 * pa_y[j] * fl3_fx * fl1_fz - 2.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_xxy[j] * fl2_fx * fl1_fz + 67.5 * pa2pb_xy_x[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxy_x[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxxy_x[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_xxy_xx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_xxx[j] * fl1_fz);

                t_xxxy_xxy[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pb_x[j] * fl3_fx + 0.25 * pa_xxx[j] * fl2_fx + 1.5 * pa2pb_xx_x[j] * fl2_fx + 2.25 * pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pa2pb_y_xy[j] * fl2_fx + 0.5 * pa2pb_xxxy_y[j] * fl1_fx + 0.5 * pa2pb_xxx_xx[j] * fl1_fx + 3.0 * pa2pb_xxy_xy[j] * fl1_fx + 1.5 * pa2pb_xy_xxy[j] * fl1_fx + pa2pb_xxxy_xxy[j]);

                t_xxxy_xxy[j] += fl_r_0_0 * (9.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xx_x[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_xy_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxy_y[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxx_xx[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xxy_xy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_xxy[j] * fl1_fz);

                t_xxxy_xxz[j] = fl_s_0_0 * (2.25 * pa2pb_xy_z[j] * fl2_fx + 1.5 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pa2pb_xxxy_z[j] * fl1_fx + 3.0 * pa2pb_xxy_xz[j] * fl1_fx + 1.5 * pa2pb_xy_xxz[j] * fl1_fx + pa2pb_xxxy_xxz[j]);

                t_xxxy_xxz[j] += fl_r_0_0 * (22.5 * pa2pb_xy_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxy_z[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xxy_xz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_xxz[j] * fl1_fz);

                t_xxxy_xyy[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pb_y[j] * fl3_fx + 0.75 * pa_xxy[j] * fl2_fx + 1.5 * pa2pb_xx_y[j] * fl2_fx + 0.75 * pa2pb_xy_x[j] * fl2_fx + 1.5 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pa2pb_y_yy[j] * fl2_fx + 0.5 * pa2pb_xxxy_x[j] * fl1_fx + pa2pb_xxx_xy[j] * fl1_fx + 1.5 * pa2pb_xxy_yy[j] * fl1_fx + 1.5 * pa2pb_xy_xyy[j] * fl1_fx + pa2pb_xxxy_xyy[j]);

                t_xxxy_xyy[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 6.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa_xxy[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xx_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxy_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxy_yy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_xyy[j] * fl1_fz);

                t_xxxy_xyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_xx_z[j] * fl2_fx + 0.75 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pa2pb_y_yz[j] * fl2_fx + 0.5 * pa2pb_xxx_xz[j] * fl1_fx + 1.5 * pa2pb_xxy_yz[j] * fl1_fx + 1.5 * pa2pb_xy_xyz[j] * fl1_fx + pa2pb_xxxy_xyz[j]);

                t_xxxy_xyz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fga + 7.5 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxx_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxy_yz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_15_20(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_xy_xzz = pa2pbDistances.data(646 * idx + 90);

            auto pa2pb_xy_yyy = pa2pbDistances.data(646 * idx + 91);

            auto pa2pb_xy_yyz = pa2pbDistances.data(646 * idx + 92);

            auto pa2pb_xy_yzz = pa2pbDistances.data(646 * idx + 93);

            auto pa2pb_xy_zzz = pa2pbDistances.data(646 * idx + 94);

            auto pa2pb_xxx_yy = pa2pbDistances.data(646 * idx + 177);

            auto pa2pb_xxx_yz = pa2pbDistances.data(646 * idx + 178);

            auto pa2pb_xxx_zz = pa2pbDistances.data(646 * idx + 179);

            auto pa2pb_xxy_zz = pa2pbDistances.data(646 * idx + 198);

            auto pa2pb_xxxy_x = pa2pbDistances.data(646 * idx + 380);

            auto pa2pb_xxxy_y = pa2pbDistances.data(646 * idx + 381);

            auto pa2pb_xxxy_z = pa2pbDistances.data(646 * idx + 382);

            auto pa2pb_xxxy_xzz = pa2pbDistances.data(646 * idx + 394);

            auto pa2pb_xxxy_yyy = pa2pbDistances.data(646 * idx + 395);

            auto pa2pb_xxxy_yyz = pa2pbDistances.data(646 * idx + 396);

            auto pa2pb_xxxy_yzz = pa2pbDistances.data(646 * idx + 397);

            auto pa2pb_xxxy_zzz = pa2pbDistances.data(646 * idx + 398);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxy_xzz = primBuffer.data(150 * idx + 15);

            auto t_xxxy_yyy = primBuffer.data(150 * idx + 16);

            auto t_xxxy_yyz = primBuffer.data(150 * idx + 17);

            auto t_xxxy_yzz = primBuffer.data(150 * idx + 18);

            auto t_xxxy_zzz = primBuffer.data(150 * idx + 19);

            // Batch of Integrals (15,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xxx_yy, \
                                     pa2pb_xxx_yz, pa2pb_xxx_zz, pa2pb_xxxy_x, pa2pb_xxxy_xzz, pa2pb_xxxy_y, \
                                     pa2pb_xxxy_yyy, pa2pb_xxxy_yyz, pa2pb_xxxy_yzz, pa2pb_xxxy_z, pa2pb_xxxy_zzz, \
                                     pa2pb_xxy_zz, pa2pb_xy_x, pa2pb_xy_xzz, pa2pb_xy_y, pa2pb_xy_yyy, pa2pb_xy_yyz, \
                                     pa2pb_xy_yzz, pa2pb_xy_z, pa2pb_xy_zzz, pa2pb_y_zz, pa_x, pa_xxx, pa_xxy, pa_y, r_0_0, \
                                     s_0_0, t_xxxy_xzz, t_xxxy_yyy, t_xxxy_yyz, t_xxxy_yzz, t_xxxy_zzz: VLX_ALIGN)
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

                t_xxxy_xzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa_xxy[j] * fl2_fx + 0.75 * pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pa2pb_xxxy_x[j] * fl1_fx + 1.5 * pa2pb_xxy_zz[j] * fl1_fx + 1.5 * pa2pb_xy_xzz[j] * fl1_fx + pa2pb_xxxy_xzz[j]);

                t_xxxy_xzz[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 7.5 * pa_xxy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxy_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxy_x[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxy_zz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_xzz[j] * fl1_fz);

                t_xxxy_yyy[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pa_xxx[j] * fl2_fx + 2.25 * pa2pb_xy_y[j] * fl2_fx + 2.25 * pa2pb_x_yy[j] * fl2_fx + 1.5 * pa2pb_xxxy_y[j] * fl1_fx + 1.5 * pa2pb_xxx_yy[j] * fl1_fx + 1.5 * pa2pb_xy_yyy[j] * fl1_fx + pa2pb_xxxy_yyy[j]);

                t_xxxy_yyy[j] += fl_r_0_0 * (-2.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xxx[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxxy_y[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_xy_y[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxxy_y[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxx_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_yyy[j] * fl1_fz);

                t_xxxy_yyz[j] = fl_s_0_0 * (0.75 * pa2pb_xy_z[j] * fl2_fx + 1.5 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_xxxy_z[j] * fl1_fx + pa2pb_xxx_yz[j] * fl1_fx + 1.5 * pa2pb_xy_yyz[j] * fl1_fx + pa2pb_xxxy_yyz[j]);

                t_xxxy_yyz[j] += fl_r_0_0 * (-1.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxxy_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxy_z[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_yyz[j] * fl1_fz);

                t_xxxy_yzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xxx[j] * fl2_fx + 0.75 * pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_x_zz[j] * fl2_fx + 0.5 * pa2pb_xxxy_y[j] * fl1_fx + 0.5 * pa2pb_xxx_zz[j] * fl1_fx + 1.5 * pa2pb_xy_yzz[j] * fl1_fx + pa2pb_xxxy_yzz[j]);

                t_xxxy_yzz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxxy_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_y[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxy_y[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxx_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_yzz[j] * fl1_fz);

                t_xxxy_zzz[j] = fl_s_0_0 * (2.25 * pa2pb_xy_z[j] * fl2_fx + 1.5 * pa2pb_xxxy_z[j] * fl1_fx + 1.5 * pa2pb_xy_zzz[j] * fl1_fx + pa2pb_xxxy_zzz[j]);

                t_xxxy_zzz[j] += fl_r_0_0 * (-4.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxxy_z[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_xy_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxxy_z[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_zzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxy_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_20_25(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_xz_xxx = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xz_xxy = pa2pbDistances.data(646 * idx + 105);

            auto pa2pb_xz_xxz = pa2pbDistances.data(646 * idx + 106);

            auto pa2pb_xz_xyy = pa2pbDistances.data(646 * idx + 107);

            auto pa2pb_xz_xyz = pa2pbDistances.data(646 * idx + 108);

            auto pa2pb_xxx_xx = pa2pbDistances.data(646 * idx + 174);

            auto pa2pb_xxx_xy = pa2pbDistances.data(646 * idx + 175);

            auto pa2pb_xxz_xx = pa2pbDistances.data(646 * idx + 212);

            auto pa2pb_xxz_xy = pa2pbDistances.data(646 * idx + 213);

            auto pa2pb_xxz_xz = pa2pbDistances.data(646 * idx + 214);

            auto pa2pb_xxz_yy = pa2pbDistances.data(646 * idx + 215);

            auto pa2pb_xxz_yz = pa2pbDistances.data(646 * idx + 216);

            auto pa2pb_xxxz_x = pa2pbDistances.data(646 * idx + 399);

            auto pa2pb_xxxz_y = pa2pbDistances.data(646 * idx + 400);

            auto pa2pb_xxxz_z = pa2pbDistances.data(646 * idx + 401);

            auto pa2pb_xxxz_xxx = pa2pbDistances.data(646 * idx + 408);

            auto pa2pb_xxxz_xxy = pa2pbDistances.data(646 * idx + 409);

            auto pa2pb_xxxz_xxz = pa2pbDistances.data(646 * idx + 410);

            auto pa2pb_xxxz_xyy = pa2pbDistances.data(646 * idx + 411);

            auto pa2pb_xxxz_xyz = pa2pbDistances.data(646 * idx + 412);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxz_xxx = primBuffer.data(150 * idx + 20);

            auto t_xxxz_xxy = primBuffer.data(150 * idx + 21);

            auto t_xxxz_xxz = primBuffer.data(150 * idx + 22);

            auto t_xxxz_xyy = primBuffer.data(150 * idx + 23);

            auto t_xxxz_xyz = primBuffer.data(150 * idx + 24);

            // Batch of Integrals (20,25)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_xx_x, pa2pb_xx_y, \
                                     pa2pb_xxx_xx, pa2pb_xxx_xy, pa2pb_xxxz_x, pa2pb_xxxz_xxx, pa2pb_xxxz_xxy, \
                                     pa2pb_xxxz_xxz, pa2pb_xxxz_xyy, pa2pb_xxxz_xyz, pa2pb_xxxz_y, pa2pb_xxxz_z, \
                                     pa2pb_xxz_xx, pa2pb_xxz_xy, pa2pb_xxz_xz, pa2pb_xxz_yy, pa2pb_xxz_yz, pa2pb_xz_x, \
                                     pa2pb_xz_xxx, pa2pb_xz_xxy, pa2pb_xz_xxz, pa2pb_xz_xyy, pa2pb_xz_xyz, pa2pb_xz_y, \
                                     pa2pb_xz_z, pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, pa_x, pa_xxx, \
                                     pa_xxz, pa_z, pb_x, pb_y, r_0_0, s_0_0, t_xxxz_xxx, t_xxxz_xxy, t_xxxz_xxz, \
                                     t_xxxz_xyy, t_xxxz_xyz: VLX_ALIGN)
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

                t_xxxz_xxx[j] = fl_s_0_0 * (1.875 * pa_z[j] * fl3_fx + 2.25 * pa_xxz[j] * fl2_fx + 6.75 * pa2pb_xz_x[j] * fl2_fx + 2.25 * pa2pb_z_xx[j] * fl2_fx + 1.5 * pa2pb_xxxz_x[j] * fl1_fx + 4.5 * pa2pb_xxz_xx[j] * fl1_fx + 1.5 * pa2pb_xz_xxx[j] * fl1_fx + pa2pb_xxxz_xxx[j]);

                t_xxxz_xxx[j] += fl_r_0_0 * (15.0 * pa_z[j] * fl3_fx * fl1_fz - 2.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_xxz[j] * fl2_fx * fl1_fz + 67.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xxxz_x[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxxz_x[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_xxz_xx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_xxx[j] * fl1_fz);

                t_xxxz_xxy[j] = fl_s_0_0 * (2.25 * pa2pb_xz_y[j] * fl2_fx + 1.5 * pa2pb_z_xy[j] * fl2_fx + 0.5 * pa2pb_xxxz_y[j] * fl1_fx + 3.0 * pa2pb_xxz_xy[j] * fl1_fx + 1.5 * pa2pb_xz_xxy[j] * fl1_fx + pa2pb_xxxz_xxy[j]);

                t_xxxz_xxy[j] += fl_r_0_0 * (22.5 * pa2pb_xz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxz_y[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xxz_xy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_xxy[j] * fl1_fz);

                t_xxxz_xxz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pb_x[j] * fl3_fx + 0.25 * pa_xxx[j] * fl2_fx + 1.5 * pa2pb_xx_x[j] * fl2_fx + 2.25 * pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_xxxz_z[j] * fl1_fx + 0.5 * pa2pb_xxx_xx[j] * fl1_fx + 3.0 * pa2pb_xxz_xz[j] * fl1_fx + 1.5 * pa2pb_xz_xxz[j] * fl1_fx + pa2pb_xxxz_xxz[j]);

                t_xxxz_xxz[j] += fl_r_0_0 * (9.0 * pa_x[j] * fl3_fx * fl1_fz - 0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_xx_x[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxz_z[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxx_xx[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xxz_xz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_xxz[j] * fl1_fz);

                t_xxxz_xyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa_xxz[j] * fl2_fx + 0.75 * pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pa2pb_xxxz_x[j] * fl1_fx + 1.5 * pa2pb_xxz_yy[j] * fl1_fx + 1.5 * pa2pb_xz_xyy[j] * fl1_fx + pa2pb_xxxz_xyy[j]);

                t_xxxz_xyy[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 7.5 * pa_xxz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxz_x[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxz_yy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_xyy[j] * fl1_fz);

                t_xxxz_xyz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_xx_y[j] * fl2_fx + 0.75 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_xxx_xy[j] * fl1_fx + 1.5 * pa2pb_xxz_yz[j] * fl1_fx + 1.5 * pa2pb_xz_xyz[j] * fl1_fx + pa2pb_xxxz_xyz[j]);

                t_xxxz_xyz[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fga + 7.5 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxx_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxz_yz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_25_30(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_xz_xzz = pa2pbDistances.data(646 * idx + 109);

            auto pa2pb_xz_yyy = pa2pbDistances.data(646 * idx + 110);

            auto pa2pb_xz_yyz = pa2pbDistances.data(646 * idx + 111);

            auto pa2pb_xz_yzz = pa2pbDistances.data(646 * idx + 112);

            auto pa2pb_xz_zzz = pa2pbDistances.data(646 * idx + 113);

            auto pa2pb_xxx_xz = pa2pbDistances.data(646 * idx + 176);

            auto pa2pb_xxx_yy = pa2pbDistances.data(646 * idx + 177);

            auto pa2pb_xxx_yz = pa2pbDistances.data(646 * idx + 178);

            auto pa2pb_xxx_zz = pa2pbDistances.data(646 * idx + 179);

            auto pa2pb_xxz_zz = pa2pbDistances.data(646 * idx + 217);

            auto pa2pb_xxxz_x = pa2pbDistances.data(646 * idx + 399);

            auto pa2pb_xxxz_y = pa2pbDistances.data(646 * idx + 400);

            auto pa2pb_xxxz_z = pa2pbDistances.data(646 * idx + 401);

            auto pa2pb_xxxz_xzz = pa2pbDistances.data(646 * idx + 413);

            auto pa2pb_xxxz_yyy = pa2pbDistances.data(646 * idx + 414);

            auto pa2pb_xxxz_yyz = pa2pbDistances.data(646 * idx + 415);

            auto pa2pb_xxxz_yzz = pa2pbDistances.data(646 * idx + 416);

            auto pa2pb_xxxz_zzz = pa2pbDistances.data(646 * idx + 417);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxxz_xzz = primBuffer.data(150 * idx + 25);

            auto t_xxxz_yyy = primBuffer.data(150 * idx + 26);

            auto t_xxxz_yyz = primBuffer.data(150 * idx + 27);

            auto t_xxxz_yzz = primBuffer.data(150 * idx + 28);

            auto t_xxxz_zzz = primBuffer.data(150 * idx + 29);

            // Batch of Integrals (25,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xx_z, \
                                     pa2pb_xxx_xz, pa2pb_xxx_yy, pa2pb_xxx_yz, pa2pb_xxx_zz, pa2pb_xxxz_x, \
                                     pa2pb_xxxz_xzz, pa2pb_xxxz_y, pa2pb_xxxz_yyy, pa2pb_xxxz_yyz, pa2pb_xxxz_yzz, \
                                     pa2pb_xxxz_z, pa2pb_xxxz_zzz, pa2pb_xxz_zz, pa2pb_xz_x, pa2pb_xz_xzz, pa2pb_xz_y, \
                                     pa2pb_xz_yyy, pa2pb_xz_yyz, pa2pb_xz_yzz, pa2pb_xz_z, pa2pb_xz_zzz, pa2pb_z_zz, pa_x, \
                                     pa_xxx, pa_xxz, pa_z, pb_z, r_0_0, s_0_0, t_xxxz_xzz, t_xxxz_yyy, t_xxxz_yyz, \
                                     t_xxxz_yzz, t_xxxz_zzz: VLX_ALIGN)
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

                t_xxxz_xzz[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pb_z[j] * fl3_fx + 0.75 * pa_xxz[j] * fl2_fx + 1.5 * pa2pb_xx_z[j] * fl2_fx + 0.75 * pa2pb_xz_x[j] * fl2_fx + 1.5 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pa2pb_xxxz_x[j] * fl1_fx + pa2pb_xxx_xz[j] * fl1_fx + 1.5 * pa2pb_xxz_zz[j] * fl1_fx + 1.5 * pa2pb_xz_xzz[j] * fl1_fx + pa2pb_xxxz_xzz[j]);

                t_xxxz_xzz[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 6.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa_xxz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_xx_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xxxz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxz_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_xz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxz_zz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_xzz[j] * fl1_fz);

                t_xxxz_yyy[j] = fl_s_0_0 * (2.25 * pa2pb_xz_y[j] * fl2_fx + 1.5 * pa2pb_xxxz_y[j] * fl1_fx + 1.5 * pa2pb_xz_yyy[j] * fl1_fx + pa2pb_xxxz_yyy[j]);

                t_xxxz_yyy[j] += fl_r_0_0 * (-4.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxxz_y[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_xz_y[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxxz_y[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_yyy[j] * fl1_fz);

                t_xxxz_yyz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xxx[j] * fl2_fx + 0.75 * pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_yy[j] * fl2_fx + 0.5 * pa2pb_xxxz_z[j] * fl1_fx + 0.5 * pa2pb_xxx_yy[j] * fl1_fx + 1.5 * pa2pb_xz_yyz[j] * fl1_fx + pa2pb_xxxz_yyz[j]);

                t_xxxz_yyz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xxx[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxxz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxz_z[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxx_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_yyz[j] * fl1_fz);

                t_xxxz_yzz[j] = fl_s_0_0 * (0.75 * pa2pb_xz_y[j] * fl2_fx + 1.5 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_xxxz_y[j] * fl1_fx + pa2pb_xxx_yz[j] * fl1_fx + 1.5 * pa2pb_xz_yzz[j] * fl1_fx + pa2pb_xxxz_yzz[j]);

                t_xxxz_yzz[j] += fl_r_0_0 * (-1.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxxz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_y[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxxz_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxx_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_yzz[j] * fl1_fz);

                t_xxxz_zzz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx + 0.75 * pa_xxx[j] * fl2_fx + 2.25 * pa2pb_xz_z[j] * fl2_fx + 2.25 * pa2pb_x_zz[j] * fl2_fx + 1.5 * pa2pb_xxxz_z[j] * fl1_fx + 1.5 * pa2pb_xxx_zz[j] * fl1_fx + 1.5 * pa2pb_xz_zzz[j] * fl1_fx + pa2pb_xxxz_zzz[j]);

                t_xxxz_zzz[j] += fl_r_0_0 * (-2.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xxx[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xxx[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxxz_z[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxxz_z[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxx_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_zzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxxz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_30_35(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xx_xxx = pa2pbDistances.data(646 * idx + 66);

            auto pa2pb_xx_xxy = pa2pbDistances.data(646 * idx + 67);

            auto pa2pb_xx_xxz = pa2pbDistances.data(646 * idx + 68);

            auto pa2pb_xx_xyy = pa2pbDistances.data(646 * idx + 69);

            auto pa2pb_xx_xyz = pa2pbDistances.data(646 * idx + 70);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_yy_xxx = pa2pbDistances.data(646 * idx + 123);

            auto pa2pb_yy_xxy = pa2pbDistances.data(646 * idx + 124);

            auto pa2pb_yy_xxz = pa2pbDistances.data(646 * idx + 125);

            auto pa2pb_yy_xyy = pa2pbDistances.data(646 * idx + 126);

            auto pa2pb_yy_xyz = pa2pbDistances.data(646 * idx + 127);

            auto pa2pb_xxy_xx = pa2pbDistances.data(646 * idx + 193);

            auto pa2pb_xxy_xy = pa2pbDistances.data(646 * idx + 194);

            auto pa2pb_xxy_xz = pa2pbDistances.data(646 * idx + 195);

            auto pa2pb_xyy_xx = pa2pbDistances.data(646 * idx + 231);

            auto pa2pb_xyy_xy = pa2pbDistances.data(646 * idx + 232);

            auto pa2pb_xyy_xz = pa2pbDistances.data(646 * idx + 233);

            auto pa2pb_xyy_yy = pa2pbDistances.data(646 * idx + 234);

            auto pa2pb_xyy_yz = pa2pbDistances.data(646 * idx + 235);

            auto pa2pb_xxyy_x = pa2pbDistances.data(646 * idx + 418);

            auto pa2pb_xxyy_y = pa2pbDistances.data(646 * idx + 419);

            auto pa2pb_xxyy_z = pa2pbDistances.data(646 * idx + 420);

            auto pa2pb_xxyy_xxx = pa2pbDistances.data(646 * idx + 427);

            auto pa2pb_xxyy_xxy = pa2pbDistances.data(646 * idx + 428);

            auto pa2pb_xxyy_xxz = pa2pbDistances.data(646 * idx + 429);

            auto pa2pb_xxyy_xyy = pa2pbDistances.data(646 * idx + 430);

            auto pa2pb_xxyy_xyz = pa2pbDistances.data(646 * idx + 431);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyy_xxx = primBuffer.data(150 * idx + 30);

            auto t_xxyy_xxy = primBuffer.data(150 * idx + 31);

            auto t_xxyy_xxz = primBuffer.data(150 * idx + 32);

            auto t_xxyy_xyy = primBuffer.data(150 * idx + 33);

            auto t_xxyy_xyz = primBuffer.data(150 * idx + 34);

            // Batch of Integrals (30,35)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_xx_x, pa2pb_xx_xxx, pa2pb_xx_xxy, pa2pb_xx_xxz, pa2pb_xx_xyy, pa2pb_xx_xyz, \
                                     pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxy_xx, pa2pb_xxy_xy, pa2pb_xxy_xz, pa2pb_xxyy_x, \
                                     pa2pb_xxyy_xxx, pa2pb_xxyy_xxy, pa2pb_xxyy_xxz, pa2pb_xxyy_xyy, pa2pb_xxyy_xyz, \
                                     pa2pb_xxyy_y, pa2pb_xxyy_z, pa2pb_xy_x, pa2pb_xy_y, pa2pb_xy_z, pa2pb_xyy_xx, \
                                     pa2pb_xyy_xy, pa2pb_xyy_xz, pa2pb_xyy_yy, pa2pb_xyy_yz, pa2pb_y_xx, pa2pb_y_xy, \
                                     pa2pb_y_xz, pa2pb_yy_x, pa2pb_yy_xxx, pa2pb_yy_xxy, pa2pb_yy_xxz, pa2pb_yy_xyy, \
                                     pa2pb_yy_xyz, pa2pb_yy_y, pa2pb_yy_z, pa_x, pa_xxy, pa_xyy, pa_y, pb_x, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xyy, pb_xyz, pb_y, pb_z, r_0_0, s_0_0, t_xxyy_xxx, t_xxyy_xxy, t_xxyy_xxz, \
                                     t_xxyy_xyy, t_xxyy_xyz: VLX_ALIGN)
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

                t_xxyy_xxx[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 1.125 * pb_x[j] * fl3_fx + 1.5 * pa_xyy[j] * fl2_fx + 2.25 * pa2pb_yy_x[j] * fl2_fx + 0.75 * pa2pb_xx_x[j] * fl2_fx + 1.5 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pa2pb_xxyy_x[j] * fl1_fx + 3.0 * pa2pb_xyy_xx[j] * fl1_fx + 0.25 * pb_xxx[j] * fl2_fx + 0.5 * pa2pb_xx_xxx[j] * fl1_fx + 0.5 * pa2pb_yy_xxx[j] * fl1_fx + pa2pb_xxyy_xxx[j]);

                t_xxyy_xxx[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz + 9.0 * pb_x[j] * fl3_fx * fl1_fz + 15.0 * pa_xyy[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_yy_x[j] * fl2_fx * fl1_fz - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxyy_x[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xyy_xx[j] * fl1_fx * fl1_fz - pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxx[j] * fl1_fz * fl1_fga + 2.5 * pb_xxx[j] * fl2_fx * fl1_fz - pa2pb_yy_xxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xxx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_xxx[j] * fl1_fz);

                t_xxyy_xxy[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 0.375 * pb_y[j] * fl3_fx + 0.5 * pa_xxy[j] * fl2_fx + 2.0 * pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_yy_y[j] * fl2_fx + 0.25 * pa2pb_xx_y[j] * fl2_fx + pa2pb_x_xy[j] * fl2_fx + 0.5 * pa2pb_y_xx[j] * fl2_fx + 0.5 * pa2pb_xxyy_y[j] * fl1_fx + pa2pb_xxy_xx[j] * fl1_fx + 2.0 * pa2pb_xyy_xy[j] * fl1_fx + 0.25 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_xx_xxy[j] * fl1_fx + 0.5 * pa2pb_yy_xxy[j] * fl1_fx + pa2pb_xxyy_xxy[j]);

                t_xxyy_xxy[j] += fl_r_0_0 * (6.0 * pa_y[j] * fl3_fx * fl1_fz - 0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 5.0 * pa_xxy[j] * fl1_fz * fl2_fx + 20.0 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_yy_y[j] * fl2_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyy_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxy_xx[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xyy_xy[j] * fl1_fx * fl1_fz - pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxy[j] * fl1_fz * fl1_fga + 2.5 * pb_xxy[j] * fl2_fx * fl1_fz - pa2pb_yy_xxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xxy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_xxy[j] * fl1_fz);

                t_xxyy_xxz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_yy_z[j] * fl2_fx + 0.25 * pa2pb_xx_z[j] * fl2_fx + pa2pb_x_xz[j] * fl2_fx + 0.5 * pa2pb_xxyy_z[j] * fl1_fx + 2.0 * pa2pb_xyy_xz[j] * fl1_fx + 0.25 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_xx_xxz[j] * fl1_fx + 0.5 * pa2pb_yy_xxz[j] * fl1_fx + pa2pb_xxyy_xxz[j]);

                t_xxyy_xxz[j] += fl_r_0_0 * (-pb_z[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_z[j] * fl2_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyy_z[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xyy_xz[j] * fl1_fx * fl1_fz - pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxz[j] * fl1_fz * fl1_fga + 2.5 * pb_xxz[j] * fl2_fx * fl1_fz - pa2pb_yy_xxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xxz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_xxz[j] * fl1_fz);

                t_xxyy_xyy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_xx_x[j] * fl2_fx + 0.5 * pa_xyy[j] * fl2_fx + 2.0 * pa2pb_xy_y[j] * fl2_fx + 0.5 * pa2pb_x_yy[j] * fl2_fx + 0.25 * pa2pb_yy_x[j] * fl2_fx + pa2pb_y_xy[j] * fl2_fx + 0.5 * pa2pb_xxyy_x[j] * fl1_fx + 2.0 * pa2pb_xxy_xy[j] * fl1_fx + pa2pb_xyy_yy[j] * fl1_fx + 0.25 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_xx_xyy[j] * fl1_fx + 0.5 * pa2pb_yy_xyy[j] * fl1_fx + pa2pb_xxyy_xyy[j]);

                t_xxyy_xyy[j] += fl_r_0_0 * (6.0 * pa_x[j] * fl3_fx * fl1_fz - 0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - pb_x[j] * fl1_fz * fl1_fga * fl2_fx - pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_x[j] * fl2_fx * fl1_fz + 5.0 * pa_xyy[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_xy_y[j] * fl2_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_x[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_yy_x[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyy_x[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xxy_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_yy[j] * fl1_fx * fl1_fz - pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyy[j] * fl1_fz * fl1_fga + 2.5 * pb_xyy[j] * fl2_fx * fl1_fz - pa2pb_yy_xyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xyy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_xyy[j] * fl1_fz);

                t_xxyy_xyz[j] = fl_s_0_0 * (pa2pb_xy_z[j] * fl2_fx + 0.5 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_y_xz[j] * fl2_fx + pa2pb_xxy_xz[j] * fl1_fx + pa2pb_xyy_yz[j] * fl1_fx + 0.25 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_xx_xyz[j] * fl1_fx + 0.5 * pa2pb_yy_xyz[j] * fl1_fx + pa2pb_xxyy_xyz[j]);

                t_xxyy_xyz[j] += fl_r_0_0 * (10.0 * pa2pb_xy_z[j] * fl2_fx * fl1_fz - pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx + 5.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xxy_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_yz[j] * fl1_fx * fl1_fz - pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyz[j] * fl1_fz * fl1_fga + 2.5 * pb_xyz[j] * fl2_fx * fl1_fz - pa2pb_yy_xyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_35_40(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xx_xzz = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_xx_yyy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_xx_yyz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_xx_yzz = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_xx_zzz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_yy_xzz = pa2pbDistances.data(646 * idx + 128);

            auto pa2pb_yy_yyy = pa2pbDistances.data(646 * idx + 129);

            auto pa2pb_yy_yyz = pa2pbDistances.data(646 * idx + 130);

            auto pa2pb_yy_yzz = pa2pbDistances.data(646 * idx + 131);

            auto pa2pb_yy_zzz = pa2pbDistances.data(646 * idx + 132);

            auto pa2pb_xxy_yy = pa2pbDistances.data(646 * idx + 196);

            auto pa2pb_xxy_yz = pa2pbDistances.data(646 * idx + 197);

            auto pa2pb_xxy_zz = pa2pbDistances.data(646 * idx + 198);

            auto pa2pb_xyy_zz = pa2pbDistances.data(646 * idx + 236);

            auto pa2pb_xxyy_x = pa2pbDistances.data(646 * idx + 418);

            auto pa2pb_xxyy_y = pa2pbDistances.data(646 * idx + 419);

            auto pa2pb_xxyy_z = pa2pbDistances.data(646 * idx + 420);

            auto pa2pb_xxyy_xzz = pa2pbDistances.data(646 * idx + 432);

            auto pa2pb_xxyy_yyy = pa2pbDistances.data(646 * idx + 433);

            auto pa2pb_xxyy_yyz = pa2pbDistances.data(646 * idx + 434);

            auto pa2pb_xxyy_yzz = pa2pbDistances.data(646 * idx + 435);

            auto pa2pb_xxyy_zzz = pa2pbDistances.data(646 * idx + 436);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyy_xzz = primBuffer.data(150 * idx + 35);

            auto t_xxyy_yyy = primBuffer.data(150 * idx + 36);

            auto t_xxyy_yyz = primBuffer.data(150 * idx + 37);

            auto t_xxyy_yzz = primBuffer.data(150 * idx + 38);

            auto t_xxyy_zzz = primBuffer.data(150 * idx + 39);

            // Batch of Integrals (35,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_zz, pa2pb_xx_x, pa2pb_xx_xzz, pa2pb_xx_y, \
                                     pa2pb_xx_yyy, pa2pb_xx_yyz, pa2pb_xx_yzz, pa2pb_xx_z, pa2pb_xx_zzz, pa2pb_xxy_yy, \
                                     pa2pb_xxy_yz, pa2pb_xxy_zz, pa2pb_xxyy_x, pa2pb_xxyy_xzz, pa2pb_xxyy_y, \
                                     pa2pb_xxyy_yyy, pa2pb_xxyy_yyz, pa2pb_xxyy_yzz, pa2pb_xxyy_z, pa2pb_xxyy_zzz, \
                                     pa2pb_xyy_zz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yy_x, pa2pb_yy_xzz, \
                                     pa2pb_yy_y, pa2pb_yy_yyy, pa2pb_yy_yyz, pa2pb_yy_yzz, pa2pb_yy_z, pa2pb_yy_zzz, \
                                     pa_x, pa_xxy, pa_xyy, pa_y, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, \
                                     r_0_0, s_0_0, t_xxyy_xzz, t_xxyy_yyy, t_xxyy_yyz, t_xxyy_yzz, t_xxyy_zzz: VLX_ALIGN)
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

                t_xxyy_xzz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl3_fx + 0.5 * pa_xyy[j] * fl2_fx + 0.125 * pb_x[j] * fl3_fx + 0.25 * pa2pb_xx_x[j] * fl2_fx + 0.5 * pa2pb_x_zz[j] * fl2_fx + 0.25 * pa2pb_yy_x[j] * fl2_fx + 0.5 * pa2pb_xxyy_x[j] * fl1_fx + pa2pb_xyy_zz[j] * fl1_fx + 0.25 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_xx_xzz[j] * fl1_fx + 0.5 * pa2pb_yy_xzz[j] * fl1_fx + pa2pb_xxyy_xzz[j]);

                t_xxyy_xzz[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl3_fx * fl1_fz + 5.0 * pa_xyy[j] * fl2_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb + pb_x[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_yy_x[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyy_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_zz[j] * fl1_fx * fl1_fz - pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xzz[j] * fl1_fz * fl1_fga + 2.5 * pb_xzz[j] * fl2_fx * fl1_fz - pa2pb_yy_xzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_xzz[j] * fl1_fz);

                t_xxyy_yyy[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 1.125 * pb_y[j] * fl3_fx + 1.5 * pa_xxy[j] * fl2_fx + 2.25 * pa2pb_xx_y[j] * fl2_fx + 0.75 * pa2pb_yy_y[j] * fl2_fx + 1.5 * pa2pb_y_yy[j] * fl2_fx + 1.5 * pa2pb_xxyy_y[j] * fl1_fx + 3.0 * pa2pb_xxy_yy[j] * fl1_fx + 0.25 * pb_yyy[j] * fl2_fx + 0.5 * pa2pb_xx_yyy[j] * fl1_fx + 0.5 * pa2pb_yy_yyy[j] * fl1_fx + pa2pb_xxyy_yyy[j]);

                t_xxyy_yyy[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl3_fx * fl1_fz + 9.0 * pb_y[j] * fl3_fx * fl1_fz + 15.0 * pa_xxy[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_xx_y[j] * fl2_fx * fl1_fz - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yy_y[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxyy_y[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xxy_yy[j] * fl1_fz * fl1_fx - pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyy[j] * fl1_fz * fl1_fga + 2.5 * pb_yyy[j] * fl2_fx * fl1_fz - pa2pb_yy_yyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_yyy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_yyy[j] * fl1_fz);

                t_xxyy_yyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_xx_z[j] * fl2_fx + 0.25 * pa2pb_yy_z[j] * fl2_fx + pa2pb_y_yz[j] * fl2_fx + 0.5 * pa2pb_xxyy_z[j] * fl1_fx + 2.0 * pa2pb_xxy_yz[j] * fl1_fx + 0.25 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_xx_yyz[j] * fl1_fx + 0.5 * pa2pb_yy_yyz[j] * fl1_fx + pa2pb_xxyy_yyz[j]);

                t_xxyy_yyz[j] += fl_r_0_0 * (-pb_z[j] * fl1_fz * fl1_fga * fl2_fx + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_z[j] * fl2_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yy_z[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyy_z[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xxy_yz[j] * fl1_fz * fl1_fx - pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyz[j] * fl1_fz * fl1_fga + 2.5 * pb_yyz[j] * fl2_fx * fl1_fz - pa2pb_yy_yyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_yyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_yyz[j] * fl1_fz);

                t_xxyy_yzz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl3_fx + 0.5 * pa_xxy[j] * fl2_fx + 0.125 * pb_y[j] * fl3_fx + 0.25 * pa2pb_xx_y[j] * fl2_fx + 0.25 * pa2pb_yy_y[j] * fl2_fx + 0.5 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pa2pb_xxyy_y[j] * fl1_fx + pa2pb_xxy_zz[j] * fl1_fx + 0.25 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_xx_yzz[j] * fl1_fx + 0.5 * pa2pb_yy_yzz[j] * fl1_fx + pa2pb_xxyy_yzz[j]);

                t_xxyy_yzz[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_y[j] * fl3_fx * fl1_fz + 5.0 * pa_xxy[j] * fl1_fz * fl2_fx - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb + pb_y[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyy_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_yy_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyy_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxy_zz[j] * fl1_fz * fl1_fx - pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yzz[j] * fl1_fz * fl1_fga + 2.5 * pb_yzz[j] * fl2_fx * fl1_fz - pa2pb_yy_yzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_yzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_yzz[j] * fl1_fz);

                t_xxyy_zzz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_xx_z[j] * fl2_fx + 0.75 * pa2pb_yy_z[j] * fl2_fx + 1.5 * pa2pb_xxyy_z[j] * fl1_fx + 0.25 * pb_zzz[j] * fl2_fx + 0.5 * pa2pb_xx_zzz[j] * fl1_fx + 0.5 * pa2pb_yy_zzz[j] * fl1_fx + pa2pb_xxyy_zzz[j]);

                t_xxyy_zzz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyy_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_yy_z[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxyy_z[j] * fl1_fz * fl1_fx - pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_zzz[j] * fl1_fz * fl1_fga + 2.5 * pb_zzz[j] * fl2_fx * fl1_fz - pa2pb_yy_zzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_zzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yy_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyy_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_40_45(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_yz_xxx = pa2pbDistances.data(646 * idx + 142);

            auto pa2pb_yz_xxy = pa2pbDistances.data(646 * idx + 143);

            auto pa2pb_yz_xxz = pa2pbDistances.data(646 * idx + 144);

            auto pa2pb_yz_xyy = pa2pbDistances.data(646 * idx + 145);

            auto pa2pb_yz_xyz = pa2pbDistances.data(646 * idx + 146);

            auto pa2pb_xxy_xx = pa2pbDistances.data(646 * idx + 193);

            auto pa2pb_xxy_xy = pa2pbDistances.data(646 * idx + 194);

            auto pa2pb_xxz_xx = pa2pbDistances.data(646 * idx + 212);

            auto pa2pb_xxz_xy = pa2pbDistances.data(646 * idx + 213);

            auto pa2pb_xxz_xz = pa2pbDistances.data(646 * idx + 214);

            auto pa2pb_xyz_xx = pa2pbDistances.data(646 * idx + 250);

            auto pa2pb_xyz_xy = pa2pbDistances.data(646 * idx + 251);

            auto pa2pb_xyz_xz = pa2pbDistances.data(646 * idx + 252);

            auto pa2pb_xyz_yy = pa2pbDistances.data(646 * idx + 253);

            auto pa2pb_xyz_yz = pa2pbDistances.data(646 * idx + 254);

            auto pa2pb_xxyz_x = pa2pbDistances.data(646 * idx + 437);

            auto pa2pb_xxyz_y = pa2pbDistances.data(646 * idx + 438);

            auto pa2pb_xxyz_z = pa2pbDistances.data(646 * idx + 439);

            auto pa2pb_xxyz_xxx = pa2pbDistances.data(646 * idx + 446);

            auto pa2pb_xxyz_xxy = pa2pbDistances.data(646 * idx + 447);

            auto pa2pb_xxyz_xxz = pa2pbDistances.data(646 * idx + 448);

            auto pa2pb_xxyz_xyy = pa2pbDistances.data(646 * idx + 449);

            auto pa2pb_xxyz_xyz = pa2pbDistances.data(646 * idx + 450);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyz_xxx = primBuffer.data(150 * idx + 40);

            auto t_xxyz_xxy = primBuffer.data(150 * idx + 41);

            auto t_xxyz_xxz = primBuffer.data(150 * idx + 42);

            auto t_xxyz_xyy = primBuffer.data(150 * idx + 43);

            auto t_xxyz_xyz = primBuffer.data(150 * idx + 44);

            // Batch of Integrals (40,45)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_xx_x, pa2pb_xxy_xx, pa2pb_xxy_xy, pa2pb_xxyz_x, \
                                     pa2pb_xxyz_xxx, pa2pb_xxyz_xxy, pa2pb_xxyz_xxz, pa2pb_xxyz_xyy, pa2pb_xxyz_xyz, \
                                     pa2pb_xxyz_y, pa2pb_xxyz_z, pa2pb_xxz_xx, pa2pb_xxz_xy, pa2pb_xxz_xz, pa2pb_xy_x, \
                                     pa2pb_xy_y, pa2pb_xyz_xx, pa2pb_xyz_xy, pa2pb_xyz_xz, pa2pb_xyz_yy, pa2pb_xyz_yz, \
                                     pa2pb_xz_x, pa2pb_xz_y, pa2pb_xz_z, pa2pb_y_xx, pa2pb_y_xy, pa2pb_yz_x, \
                                     pa2pb_yz_xxx, pa2pb_yz_xxy, pa2pb_yz_xxz, pa2pb_yz_xyy, pa2pb_yz_xyz, pa2pb_yz_y, \
                                     pa2pb_yz_z, pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa_x, pa_xxy, pa_xxz, pa_xyz, pa_y, pa_z, \
                                     pb_x, r_0_0, s_0_0, t_xxyz_xxx, t_xxyz_xxy, t_xxyz_xxz, t_xxyz_xyy, t_xxyz_xyz: VLX_ALIGN)
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

                t_xxyz_xxx[j] = fl_s_0_0 * (1.5 * pa_xyz[j] * fl2_fx + 2.25 * pa2pb_yz_x[j] * fl2_fx + 1.5 * pa2pb_xxyz_x[j] * fl1_fx + 3.0 * pa2pb_xyz_xx[j] * fl1_fx + 0.5 * pa2pb_yz_xxx[j] * fl1_fx + pa2pb_xxyz_xxx[j]);

                t_xxyz_xxx[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa_xyz[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_yz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_x[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xxyz_x[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xyz_xx[j] * fl1_fx * fl1_fz - pa2pb_yz_xxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_xxx[j] * fl1_fz);

                t_xxyz_xxy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.25 * pa_xxz[j] * fl2_fx + pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_yz_y[j] * fl2_fx + 0.25 * pa2pb_z_xx[j] * fl2_fx + 0.5 * pa2pb_xxyz_y[j] * fl1_fx + 0.5 * pa2pb_xxz_xx[j] * fl1_fx + 2.0 * pa2pb_xyz_xy[j] * fl1_fx + 0.5 * pa2pb_yz_xxy[j] * fl1_fx + pa2pb_xxyz_xxy[j]);

                t_xxyz_xxy[j] += fl_r_0_0 * (3.0 * pa_z[j] * fl3_fx * fl1_fz - 0.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xxz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyz_y[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxz_xx[j] * fl1_fx * fl1_fz + 24.0 * pa2pb_xyz_xy[j] * fl1_fx * fl1_fz - pa2pb_yz_xxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_xxy[j] * fl1_fz);

                t_xxyz_xxz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.25 * pa_xxy[j] * fl2_fx + pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_yz_z[j] * fl2_fx + 0.25 * pa2pb_y_xx[j] * fl2_fx + 0.5 * pa2pb_xxyz_z[j] * fl1_fx + 0.5 * pa2pb_xxy_xx[j] * fl1_fx + 2.0 * pa2pb_xyz_xz[j] * fl1_fx + 0.5 * pa2pb_yz_xxz[j] * fl1_fx + pa2pb_xxyz_xxz[j]);

                t_xxyz_xxz[j] += fl_r_0_0 * (3.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xxy[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyz_z[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxy_xx[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xyz_xz[j] * fl1_fx * fl1_fz - pa2pb_yz_xxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_xxz[j] * fl1_fz);

                t_xxyz_xyy[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl2_fx + pa2pb_xz_y[j] * fl2_fx + 0.25 * pa2pb_yz_x[j] * fl2_fx + 0.5 * pa2pb_z_xy[j] * fl2_fx + 0.5 * pa2pb_xxyz_x[j] * fl1_fx + pa2pb_xxz_xy[j] * fl1_fx + pa2pb_xyz_yy[j] * fl1_fx + 0.5 * pa2pb_yz_xyy[j] * fl1_fx + pa2pb_xxyz_xyy[j]);

                t_xxyz_xyy[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yz_x[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyz_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxz_xy[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_yy[j] * fl1_fx * fl1_fz - pa2pb_yz_xyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_xyy[j] * fl1_fz);

                t_xxyz_xyz[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl3_fx + 0.125 * pb_x[j] * fl3_fx + 0.25 * pa2pb_xx_x[j] * fl2_fx + 0.5 * pa2pb_xy_y[j] * fl2_fx + 0.5 * pa2pb_xz_z[j] * fl2_fx + 0.25 * pa2pb_y_xy[j] * fl2_fx + 0.25 * pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_xxy_xy[j] * fl1_fx + 0.5 * pa2pb_xxz_xz[j] * fl1_fx + pa2pb_xyz_yz[j] * fl1_fx + 0.5 * pa2pb_yz_xyz[j] * fl1_fx + pa2pb_xxyz_xyz[j]);

                t_xxyz_xyz[j] += fl_r_0_0 * (2.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pb_x[j] * fl1_fz * fl1_fga * fl2_fx + pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa2pb_xx_x[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xy_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xz_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx + 2.5 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxy_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxz_xz[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_xyz_yz[j] * fl1_fx * fl1_fz - pa2pb_yz_xyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_45_50(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_yz_xzz = pa2pbDistances.data(646 * idx + 147);

            auto pa2pb_yz_yyy = pa2pbDistances.data(646 * idx + 148);

            auto pa2pb_yz_yyz = pa2pbDistances.data(646 * idx + 149);

            auto pa2pb_yz_yzz = pa2pbDistances.data(646 * idx + 150);

            auto pa2pb_yz_zzz = pa2pbDistances.data(646 * idx + 151);

            auto pa2pb_xxy_xz = pa2pbDistances.data(646 * idx + 195);

            auto pa2pb_xxy_yy = pa2pbDistances.data(646 * idx + 196);

            auto pa2pb_xxy_yz = pa2pbDistances.data(646 * idx + 197);

            auto pa2pb_xxy_zz = pa2pbDistances.data(646 * idx + 198);

            auto pa2pb_xxz_yy = pa2pbDistances.data(646 * idx + 215);

            auto pa2pb_xxz_yz = pa2pbDistances.data(646 * idx + 216);

            auto pa2pb_xxz_zz = pa2pbDistances.data(646 * idx + 217);

            auto pa2pb_xyz_zz = pa2pbDistances.data(646 * idx + 255);

            auto pa2pb_xxyz_x = pa2pbDistances.data(646 * idx + 437);

            auto pa2pb_xxyz_y = pa2pbDistances.data(646 * idx + 438);

            auto pa2pb_xxyz_z = pa2pbDistances.data(646 * idx + 439);

            auto pa2pb_xxyz_xzz = pa2pbDistances.data(646 * idx + 451);

            auto pa2pb_xxyz_yyy = pa2pbDistances.data(646 * idx + 452);

            auto pa2pb_xxyz_yyz = pa2pbDistances.data(646 * idx + 453);

            auto pa2pb_xxyz_yzz = pa2pbDistances.data(646 * idx + 454);

            auto pa2pb_xxyz_zzz = pa2pbDistances.data(646 * idx + 455);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxyz_xzz = primBuffer.data(150 * idx + 45);

            auto t_xxyz_yyy = primBuffer.data(150 * idx + 46);

            auto t_xxyz_yyz = primBuffer.data(150 * idx + 47);

            auto t_xxyz_yzz = primBuffer.data(150 * idx + 48);

            auto t_xxyz_zzz = primBuffer.data(150 * idx + 49);

            // Batch of Integrals (45,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxy_xz, pa2pb_xxy_yy, \
                                     pa2pb_xxy_yz, pa2pb_xxy_zz, pa2pb_xxyz_x, pa2pb_xxyz_xzz, pa2pb_xxyz_y, \
                                     pa2pb_xxyz_yyy, pa2pb_xxyz_yyz, pa2pb_xxyz_yzz, pa2pb_xxyz_z, pa2pb_xxyz_zzz, \
                                     pa2pb_xxz_yy, pa2pb_xxz_yz, pa2pb_xxz_zz, pa2pb_xy_z, pa2pb_xyz_zz, pa2pb_y_xz, \
                                     pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yz_x, pa2pb_yz_xzz, pa2pb_yz_y, \
                                     pa2pb_yz_yyy, pa2pb_yz_yyz, pa2pb_yz_yzz, pa2pb_yz_z, pa2pb_yz_zzz, pa2pb_z_yy, \
                                     pa2pb_z_yz, pa2pb_z_zz, pa_xxy, pa_xxz, pa_xyz, pa_y, pa_z, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_xxyz_xzz, t_xxyz_yyy, t_xxyz_yyz, t_xxyz_yzz, t_xxyz_zzz: VLX_ALIGN)
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

                t_xxyz_xzz[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl2_fx + pa2pb_xy_z[j] * fl2_fx + 0.25 * pa2pb_yz_x[j] * fl2_fx + 0.5 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pa2pb_xxyz_x[j] * fl1_fx + pa2pb_xxy_xz[j] * fl1_fx + pa2pb_xyz_zz[j] * fl1_fx + 0.5 * pa2pb_yz_xzz[j] * fl1_fx + pa2pb_xxyz_xzz[j]);

                t_xxyz_xzz[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_xy_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yz_x[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyz_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxy_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyz_zz[j] * fl1_fx * fl1_fz - pa2pb_yz_xzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_xzz[j] * fl1_fz);

                t_xxyz_yyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa_xxz[j] * fl2_fx + 0.75 * pa2pb_yz_y[j] * fl2_fx + 0.75 * pa2pb_z_yy[j] * fl2_fx + 1.5 * pa2pb_xxyz_y[j] * fl1_fx + 1.5 * pa2pb_xxz_yy[j] * fl1_fx + 0.5 * pa2pb_yz_yyy[j] * fl1_fx + pa2pb_xxyz_yyy[j]);

                t_xxyz_yyy[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 7.5 * pa_xxz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxyz_y[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxz_yy[j] * fl1_fx * fl1_fz - pa2pb_yz_yyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_yyy[j] * fl1_fz);

                t_xxyz_yyz[j] = fl_s_0_0 * (0.125 * pa_y[j] * fl3_fx + 0.25 * pb_y[j] * fl3_fx + 0.25 * pa_xxy[j] * fl2_fx + 0.5 * pa2pb_xx_y[j] * fl2_fx + 0.25 * pa2pb_yz_z[j] * fl2_fx + 0.25 * pa2pb_y_yy[j] * fl2_fx + 0.5 * pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_xxyz_z[j] * fl1_fx + 0.5 * pa2pb_xxy_yy[j] * fl1_fx + pa2pb_xxz_yz[j] * fl1_fx + 0.5 * pa2pb_yz_yyz[j] * fl1_fx + pa2pb_xxyz_yyz[j]);

                t_xxyz_yyz[j] += fl_r_0_0 * (-0.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + pa_y[j] * fl3_fx * fl1_fz + 2.0 * pb_y[j] * fl3_fx * fl1_fz + 2.5 * pa_xxy[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_xx_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyz_z[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxy_yy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxz_yz[j] * fl1_fx * fl1_fz - pa2pb_yz_yyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_yyz[j] * fl1_fz);

                t_xxyz_yzz[j] = fl_s_0_0 * (0.125 * pa_z[j] * fl3_fx + 0.25 * pb_z[j] * fl3_fx + 0.25 * pa_xxz[j] * fl2_fx + 0.5 * pa2pb_xx_z[j] * fl2_fx + 0.25 * pa2pb_yz_y[j] * fl2_fx + 0.5 * pa2pb_y_yz[j] * fl2_fx + 0.25 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pa2pb_xxyz_y[j] * fl1_fx + pa2pb_xxy_yz[j] * fl1_fx + 0.5 * pa2pb_xxz_zz[j] * fl1_fx + 0.5 * pa2pb_yz_yzz[j] * fl1_fx + pa2pb_xxyz_yzz[j]);

                t_xxyz_yzz[j] += fl_r_0_0 * (-0.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + pa_z[j] * fl3_fx * fl1_fz + 2.0 * pb_z[j] * fl3_fx * fl1_fz + 2.5 * pa_xxz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_xx_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxyz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxyz_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxy_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xxz_zz[j] * fl1_fx * fl1_fz - pa2pb_yz_yzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_yzz[j] * fl1_fz);

                t_xxyz_zzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa_xxy[j] * fl2_fx + 0.75 * pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_zz[j] * fl2_fx + 1.5 * pa2pb_xxyz_z[j] * fl1_fx + 1.5 * pa2pb_xxy_zz[j] * fl1_fx + 0.5 * pa2pb_yz_zzz[j] * fl1_fx + pa2pb_xxyz_zzz[j]);

                t_xxyz_zzz[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xxy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 7.5 * pa_xxy[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxyz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxyz_z[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xxy_zz[j] * fl1_fz * fl1_fx - pa2pb_yz_zzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yz_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxyz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_50_55(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xx_xxx = pa2pbDistances.data(646 * idx + 66);

            auto pa2pb_xx_xxy = pa2pbDistances.data(646 * idx + 67);

            auto pa2pb_xx_xxz = pa2pbDistances.data(646 * idx + 68);

            auto pa2pb_xx_xyy = pa2pbDistances.data(646 * idx + 69);

            auto pa2pb_xx_xyz = pa2pbDistances.data(646 * idx + 70);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_zz_xxx = pa2pbDistances.data(646 * idx + 161);

            auto pa2pb_zz_xxy = pa2pbDistances.data(646 * idx + 162);

            auto pa2pb_zz_xxz = pa2pbDistances.data(646 * idx + 163);

            auto pa2pb_zz_xyy = pa2pbDistances.data(646 * idx + 164);

            auto pa2pb_zz_xyz = pa2pbDistances.data(646 * idx + 165);

            auto pa2pb_xxz_xx = pa2pbDistances.data(646 * idx + 212);

            auto pa2pb_xxz_xy = pa2pbDistances.data(646 * idx + 213);

            auto pa2pb_xzz_xx = pa2pbDistances.data(646 * idx + 269);

            auto pa2pb_xzz_xy = pa2pbDistances.data(646 * idx + 270);

            auto pa2pb_xzz_xz = pa2pbDistances.data(646 * idx + 271);

            auto pa2pb_xzz_yy = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_xzz_yz = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_xxzz_x = pa2pbDistances.data(646 * idx + 456);

            auto pa2pb_xxzz_y = pa2pbDistances.data(646 * idx + 457);

            auto pa2pb_xxzz_z = pa2pbDistances.data(646 * idx + 458);

            auto pa2pb_xxzz_xxx = pa2pbDistances.data(646 * idx + 465);

            auto pa2pb_xxzz_xxy = pa2pbDistances.data(646 * idx + 466);

            auto pa2pb_xxzz_xxz = pa2pbDistances.data(646 * idx + 467);

            auto pa2pb_xxzz_xyy = pa2pbDistances.data(646 * idx + 468);

            auto pa2pb_xxzz_xyz = pa2pbDistances.data(646 * idx + 469);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxzz_xxx = primBuffer.data(150 * idx + 50);

            auto t_xxzz_xxy = primBuffer.data(150 * idx + 51);

            auto t_xxzz_xxz = primBuffer.data(150 * idx + 52);

            auto t_xxzz_xyy = primBuffer.data(150 * idx + 53);

            auto t_xxzz_xyz = primBuffer.data(150 * idx + 54);

            // Batch of Integrals (50,55)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, \
                                     pa2pb_xx_x, pa2pb_xx_xxx, pa2pb_xx_xxy, pa2pb_xx_xxz, pa2pb_xx_xyy, pa2pb_xx_xyz, \
                                     pa2pb_xx_y, pa2pb_xx_z, pa2pb_xxz_xx, pa2pb_xxz_xy, pa2pb_xxzz_x, pa2pb_xxzz_xxx, \
                                     pa2pb_xxzz_xxy, pa2pb_xxzz_xxz, pa2pb_xxzz_xyy, pa2pb_xxzz_xyz, pa2pb_xxzz_y, \
                                     pa2pb_xxzz_z, pa2pb_xz_x, pa2pb_xz_y, pa2pb_xzz_xx, pa2pb_xzz_xy, pa2pb_xzz_xz, \
                                     pa2pb_xzz_yy, pa2pb_xzz_yz, pa2pb_z_xx, pa2pb_z_xy, pa2pb_zz_x, pa2pb_zz_xxx, \
                                     pa2pb_zz_xxy, pa2pb_zz_xxz, pa2pb_zz_xyy, pa2pb_zz_xyz, pa2pb_zz_y, pa2pb_zz_z, pa_x, \
                                     pa_xxz, pa_xzz, pa_z, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_y, pb_z, r_0_0, \
                                     s_0_0, t_xxzz_xxx, t_xxzz_xxy, t_xxzz_xxz, t_xxzz_xyy, t_xxzz_xyz: VLX_ALIGN)
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

                t_xxzz_xxx[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 1.125 * pb_x[j] * fl3_fx + 1.5 * pa_xzz[j] * fl2_fx + 2.25 * pa2pb_zz_x[j] * fl2_fx + 0.75 * pa2pb_xx_x[j] * fl2_fx + 1.5 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pa2pb_xxzz_x[j] * fl1_fx + 3.0 * pa2pb_xzz_xx[j] * fl1_fx + 0.25 * pb_xxx[j] * fl2_fx + 0.5 * pa2pb_xx_xxx[j] * fl1_fx + 0.5 * pa2pb_zz_xxx[j] * fl1_fx + pa2pb_xxzz_xxx[j]);

                t_xxzz_xxx[j] += fl_r_0_0 * (-1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_x[j] * fl3_fx * fl1_fz + 9.0 * pb_x[j] * fl3_fx * fl1_fz + 15.0 * pa_xzz[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxzz_x[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xzz_xx[j] * fl1_fx * fl1_fz - pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxx[j] * fl1_fz * fl1_fga + 2.5 * pb_xxx[j] * fl2_fx * fl1_fz - pa2pb_zz_xxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xxx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_xxx[j] * fl1_fz);

                t_xxzz_xxy[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_zz_y[j] * fl2_fx + 0.25 * pa2pb_xx_y[j] * fl2_fx + pa2pb_x_xy[j] * fl2_fx + 0.5 * pa2pb_xxzz_y[j] * fl1_fx + 2.0 * pa2pb_xzz_xy[j] * fl1_fx + 0.25 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_xx_xxy[j] * fl1_fx + 0.5 * pa2pb_zz_xxy[j] * fl1_fx + pa2pb_xxzz_xxy[j]);

                t_xxzz_xxy[j] += fl_r_0_0 * (-pb_y[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxzz_y[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xzz_xy[j] * fl1_fx * fl1_fz - pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxy[j] * fl1_fz * fl1_fga + 2.5 * pb_xxy[j] * fl2_fx * fl1_fz - pa2pb_zz_xxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xxy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_xxy[j] * fl1_fz);

                t_xxzz_xxz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl3_fx + 0.375 * pb_z[j] * fl3_fx + 0.5 * pa_xxz[j] * fl2_fx + 2.0 * pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_zz_z[j] * fl2_fx + 0.25 * pa2pb_xx_z[j] * fl2_fx + pa2pb_x_xz[j] * fl2_fx + 0.5 * pa2pb_z_xx[j] * fl2_fx + 0.5 * pa2pb_xxzz_z[j] * fl1_fx + pa2pb_xxz_xx[j] * fl1_fx + 2.0 * pa2pb_xzz_xz[j] * fl1_fx + 0.25 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_xx_xxz[j] * fl1_fx + 0.5 * pa2pb_zz_xxz[j] * fl1_fx + pa2pb_xxzz_xxz[j]);

                t_xxzz_xxz[j] += fl_r_0_0 * (6.0 * pa_z[j] * fl3_fx * fl1_fz - 0.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 5.0 * pa_xxz[j] * fl1_fz * fl2_fx + 20.0 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxzz_z[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxz_xx[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xzz_xz[j] * fl1_fx * fl1_fz - pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xxz[j] * fl1_fz * fl1_fga + 2.5 * pb_xxz[j] * fl2_fx * fl1_fz - pa2pb_zz_xxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xxz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_xxz[j] * fl1_fz);

                t_xxzz_xyy[j] = fl_s_0_0 * (0.25 * pa_x[j] * fl3_fx + 0.5 * pa_xzz[j] * fl2_fx + 0.125 * pb_x[j] * fl3_fx + 0.25 * pa2pb_xx_x[j] * fl2_fx + 0.5 * pa2pb_x_yy[j] * fl2_fx + 0.25 * pa2pb_zz_x[j] * fl2_fx + 0.5 * pa2pb_xxzz_x[j] * fl1_fx + pa2pb_xzz_yy[j] * fl1_fx + 0.25 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_xx_xyy[j] * fl1_fx + 0.5 * pa2pb_zz_xyy[j] * fl1_fx + pa2pb_xxzz_xyy[j]);

                t_xxzz_xyy[j] += fl_r_0_0 * (-0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_x[j] * fl3_fx * fl1_fz + 5.0 * pa_xzz[j] * fl2_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + pb_x[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_x[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxzz_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_yy[j] * fl1_fx * fl1_fz - pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyy[j] * fl1_fz * fl1_fga + 2.5 * pb_xyy[j] * fl2_fx * fl1_fz - pa2pb_zz_xyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xyy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_xyy[j] * fl1_fz);

                t_xxzz_xyz[j] = fl_s_0_0 * (pa2pb_xz_y[j] * fl2_fx + 0.5 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_z_xy[j] * fl2_fx + pa2pb_xxz_xy[j] * fl1_fx + pa2pb_xzz_yz[j] * fl1_fx + 0.25 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_xx_xyz[j] * fl1_fx + 0.5 * pa2pb_zz_xyz[j] * fl1_fx + pa2pb_xxzz_xyz[j]);

                t_xxzz_xyz[j] += fl_r_0_0 * (10.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz - pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx + 5.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xxz_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_yz[j] * fl1_fx * fl1_fz - pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xyz[j] * fl1_fz * fl1_fga + 2.5 * pb_xyz[j] * fl2_fx * fl1_fz - pa2pb_zz_xyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_55_60(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_xx_x = pa2pbDistances.data(646 * idx + 57);

            auto pa2pb_xx_y = pa2pbDistances.data(646 * idx + 58);

            auto pa2pb_xx_z = pa2pbDistances.data(646 * idx + 59);

            auto pa2pb_xx_xzz = pa2pbDistances.data(646 * idx + 71);

            auto pa2pb_xx_yyy = pa2pbDistances.data(646 * idx + 72);

            auto pa2pb_xx_yyz = pa2pbDistances.data(646 * idx + 73);

            auto pa2pb_xx_yzz = pa2pbDistances.data(646 * idx + 74);

            auto pa2pb_xx_zzz = pa2pbDistances.data(646 * idx + 75);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_zz_xzz = pa2pbDistances.data(646 * idx + 166);

            auto pa2pb_zz_yyy = pa2pbDistances.data(646 * idx + 167);

            auto pa2pb_zz_yyz = pa2pbDistances.data(646 * idx + 168);

            auto pa2pb_zz_yzz = pa2pbDistances.data(646 * idx + 169);

            auto pa2pb_zz_zzz = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_xxz_xz = pa2pbDistances.data(646 * idx + 214);

            auto pa2pb_xxz_yy = pa2pbDistances.data(646 * idx + 215);

            auto pa2pb_xxz_yz = pa2pbDistances.data(646 * idx + 216);

            auto pa2pb_xxz_zz = pa2pbDistances.data(646 * idx + 217);

            auto pa2pb_xzz_zz = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_xxzz_x = pa2pbDistances.data(646 * idx + 456);

            auto pa2pb_xxzz_y = pa2pbDistances.data(646 * idx + 457);

            auto pa2pb_xxzz_z = pa2pbDistances.data(646 * idx + 458);

            auto pa2pb_xxzz_xzz = pa2pbDistances.data(646 * idx + 470);

            auto pa2pb_xxzz_yyy = pa2pbDistances.data(646 * idx + 471);

            auto pa2pb_xxzz_yyz = pa2pbDistances.data(646 * idx + 472);

            auto pa2pb_xxzz_yzz = pa2pbDistances.data(646 * idx + 473);

            auto pa2pb_xxzz_zzz = pa2pbDistances.data(646 * idx + 474);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xxzz_xzz = primBuffer.data(150 * idx + 55);

            auto t_xxzz_yyy = primBuffer.data(150 * idx + 56);

            auto t_xxzz_yyz = primBuffer.data(150 * idx + 57);

            auto t_xxzz_yzz = primBuffer.data(150 * idx + 58);

            auto t_xxzz_zzz = primBuffer.data(150 * idx + 59);

            // Batch of Integrals (55,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_zz, pa2pb_xx_x, pa2pb_xx_xzz, pa2pb_xx_y, \
                                     pa2pb_xx_yyy, pa2pb_xx_yyz, pa2pb_xx_yzz, pa2pb_xx_z, pa2pb_xx_zzz, pa2pb_xxz_xz, \
                                     pa2pb_xxz_yy, pa2pb_xxz_yz, pa2pb_xxz_zz, pa2pb_xxzz_x, pa2pb_xxzz_xzz, \
                                     pa2pb_xxzz_y, pa2pb_xxzz_yyy, pa2pb_xxzz_yyz, pa2pb_xxzz_yzz, pa2pb_xxzz_z, \
                                     pa2pb_xxzz_zzz, pa2pb_xz_z, pa2pb_xzz_zz, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, \
                                     pa2pb_z_zz, pa2pb_zz_x, pa2pb_zz_xzz, pa2pb_zz_y, pa2pb_zz_yyy, pa2pb_zz_yyz, \
                                     pa2pb_zz_yzz, pa2pb_zz_z, pa2pb_zz_zzz, pa_x, pa_xxz, pa_xzz, pa_z, pb_x, pb_xzz, pb_y, \
                                     pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zzz, r_0_0, s_0_0, t_xxzz_xzz, t_xxzz_yyy, \
                                     t_xxzz_yyz, t_xxzz_yzz, t_xxzz_zzz: VLX_ALIGN)
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

                t_xxzz_xzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx + 0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_xx_x[j] * fl2_fx + 0.5 * pa_xzz[j] * fl2_fx + 2.0 * pa2pb_xz_z[j] * fl2_fx + 0.5 * pa2pb_x_zz[j] * fl2_fx + 0.25 * pa2pb_zz_x[j] * fl2_fx + pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_xxzz_x[j] * fl1_fx + 2.0 * pa2pb_xxz_xz[j] * fl1_fx + pa2pb_xzz_zz[j] * fl1_fx + 0.25 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_xx_xzz[j] * fl1_fx + 0.5 * pa2pb_zz_xzz[j] * fl1_fx + pa2pb_xxzz_xzz[j]);

                t_xxzz_xzz[j] += fl_r_0_0 * (6.0 * pa_x[j] * fl3_fx * fl1_fz - 0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - pb_x[j] * fl1_fz * fl1_fga * fl2_fx - pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_x[j] * fl2_fx * fl1_fz + 5.0 * pa_xzz[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_xz_z[j] * fl2_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_x[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxzz_x[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xxz_xz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_zz[j] * fl1_fx * fl1_fz - pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_xzz[j] * fl1_fz * fl1_fga + 2.5 * pb_xzz[j] * fl2_fx * fl1_fz - pa2pb_zz_xzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_xzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_xzz[j] * fl1_fz);

                t_xxzz_yyy[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_xx_y[j] * fl2_fx + 0.75 * pa2pb_zz_y[j] * fl2_fx + 1.5 * pa2pb_xxzz_y[j] * fl1_fx + 0.25 * pb_yyy[j] * fl2_fx + 0.5 * pa2pb_xx_yyy[j] * fl1_fx + 0.5 * pa2pb_zz_yyy[j] * fl1_fx + pa2pb_xxzz_yyy[j]);

                t_xxzz_yyy[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xx_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxzz_y[j] * fl1_fz * fl1_fx - pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyy[j] * fl1_fz * fl1_fga + 2.5 * pb_yyy[j] * fl2_fx * fl1_fz - pa2pb_zz_yyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_yyy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_yyy[j] * fl1_fz);

                t_xxzz_yyz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl3_fx + 0.5 * pa_xxz[j] * fl2_fx + 0.125 * pb_z[j] * fl3_fx + 0.25 * pa2pb_xx_z[j] * fl2_fx + 0.25 * pa2pb_zz_z[j] * fl2_fx + 0.5 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pa2pb_xxzz_z[j] * fl1_fx + pa2pb_xxz_yy[j] * fl1_fx + 0.25 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_xx_yyz[j] * fl1_fx + 0.5 * pa2pb_zz_yyz[j] * fl1_fx + pa2pb_xxzz_yyz[j]);

                t_xxzz_yyz[j] += fl_r_0_0 * (-0.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_z[j] * fl3_fx * fl1_fz + 5.0 * pa_xxz[j] * fl1_fz * fl2_fx - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + pb_z[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xx_z[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxzz_z[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xxz_yy[j] * fl1_fz * fl1_fx - pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yyz[j] * fl1_fz * fl1_fga + 2.5 * pb_yyz[j] * fl2_fx * fl1_fz - pa2pb_zz_yyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_yyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_yyz[j] * fl1_fz);

                t_xxzz_yzz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_xx_y[j] * fl2_fx + 0.25 * pa2pb_zz_y[j] * fl2_fx + pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_xxzz_y[j] * fl1_fx + 2.0 * pa2pb_xxz_yz[j] * fl1_fx + 0.25 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_xx_yzz[j] * fl1_fx + 0.5 * pa2pb_zz_yzz[j] * fl1_fx + pa2pb_xxzz_yzz[j]);

                t_xxzz_yzz[j] += fl_r_0_0 * (-pb_y[j] * fl1_fz * fl1_fga * fl2_fx + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_xx_y[j] * fl2_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xx_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xxzz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xxzz_y[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xxz_yz[j] * fl1_fz * fl1_fx - pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_yzz[j] * fl1_fz * fl1_fga + 2.5 * pb_yzz[j] * fl2_fx * fl1_fz - pa2pb_zz_yzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_yzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_yzz[j] * fl1_fz);

                t_xxzz_zzz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl3_fx + 1.125 * pb_z[j] * fl3_fx + 1.5 * pa_xxz[j] * fl2_fx + 2.25 * pa2pb_xx_z[j] * fl2_fx + 0.75 * pa2pb_zz_z[j] * fl2_fx + 1.5 * pa2pb_z_zz[j] * fl2_fx + 1.5 * pa2pb_xxzz_z[j] * fl1_fx + 3.0 * pa2pb_xxz_zz[j] * fl1_fx + 0.25 * pb_zzz[j] * fl2_fx + 0.5 * pa2pb_xx_zzz[j] * fl1_fx + 0.5 * pa2pb_zz_zzz[j] * fl1_fx + pa2pb_xxzz_zzz[j]);

                t_xxzz_zzz[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_xxz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_z[j] * fl3_fx * fl1_fz + 9.0 * pb_z[j] * fl3_fx * fl1_fz + 15.0 * pa_xxz[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_xx_z[j] * fl2_fx * fl1_fz - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xx_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xxzz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xxzz_z[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xxz_zz[j] * fl1_fz * fl1_fx - pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xx_zzz[j] * fl1_fz * fl1_fga + 2.5 * pb_zzz[j] * fl2_fx * fl1_fz - pa2pb_zz_zzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xx_zzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xxzz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_60_65(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_xy_xxx = pa2pbDistances.data(646 * idx + 85);

            auto pa2pb_xy_xxy = pa2pbDistances.data(646 * idx + 86);

            auto pa2pb_xy_xxz = pa2pbDistances.data(646 * idx + 87);

            auto pa2pb_xy_xyy = pa2pbDistances.data(646 * idx + 88);

            auto pa2pb_xy_xyz = pa2pbDistances.data(646 * idx + 89);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_xyy_xx = pa2pbDistances.data(646 * idx + 231);

            auto pa2pb_xyy_xy = pa2pbDistances.data(646 * idx + 232);

            auto pa2pb_xyy_xz = pa2pbDistances.data(646 * idx + 233);

            auto pa2pb_yyy_xx = pa2pbDistances.data(646 * idx + 288);

            auto pa2pb_yyy_xy = pa2pbDistances.data(646 * idx + 289);

            auto pa2pb_yyy_xz = pa2pbDistances.data(646 * idx + 290);

            auto pa2pb_yyy_yy = pa2pbDistances.data(646 * idx + 291);

            auto pa2pb_yyy_yz = pa2pbDistances.data(646 * idx + 292);

            auto pa2pb_xyyy_x = pa2pbDistances.data(646 * idx + 475);

            auto pa2pb_xyyy_y = pa2pbDistances.data(646 * idx + 476);

            auto pa2pb_xyyy_z = pa2pbDistances.data(646 * idx + 477);

            auto pa2pb_xyyy_xxx = pa2pbDistances.data(646 * idx + 484);

            auto pa2pb_xyyy_xxy = pa2pbDistances.data(646 * idx + 485);

            auto pa2pb_xyyy_xxz = pa2pbDistances.data(646 * idx + 486);

            auto pa2pb_xyyy_xyy = pa2pbDistances.data(646 * idx + 487);

            auto pa2pb_xyyy_xyz = pa2pbDistances.data(646 * idx + 488);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyy_xxx = primBuffer.data(150 * idx + 60);

            auto t_xyyy_xxy = primBuffer.data(150 * idx + 61);

            auto t_xyyy_xxz = primBuffer.data(150 * idx + 62);

            auto t_xyyy_xyy = primBuffer.data(150 * idx + 63);

            auto t_xyyy_xyz = primBuffer.data(150 * idx + 64);

            // Batch of Integrals (60,65)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_xy_x, \
                                     pa2pb_xy_xxx, pa2pb_xy_xxy, pa2pb_xy_xxz, pa2pb_xy_xyy, pa2pb_xy_xyz, pa2pb_xy_y, \
                                     pa2pb_xy_z, pa2pb_xyy_xx, pa2pb_xyy_xy, pa2pb_xyy_xz, pa2pb_xyyy_x, \
                                     pa2pb_xyyy_xxx, pa2pb_xyyy_xxy, pa2pb_xyyy_xxz, pa2pb_xyyy_xyy, pa2pb_xyyy_xyz, \
                                     pa2pb_xyyy_y, pa2pb_xyyy_z, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, \
                                     pa2pb_y_yz, pa2pb_yy_x, pa2pb_yy_y, pa2pb_yy_z, pa2pb_yyy_xx, pa2pb_yyy_xy, \
                                     pa2pb_yyy_xz, pa2pb_yyy_yy, pa2pb_yyy_yz, pa_x, pa_xyy, pa_y, pa_yyy, pb_x, pb_y, pb_z, \
                                     r_0_0, s_0_0, t_xyyy_xxx, t_xyyy_xxy, t_xyyy_xxz, t_xyyy_xyy, t_xyyy_xyz: VLX_ALIGN)
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

                t_xyyy_xxx[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pa_yyy[j] * fl2_fx + 2.25 * pa2pb_xy_x[j] * fl2_fx + 2.25 * pa2pb_y_xx[j] * fl2_fx + 1.5 * pa2pb_xyyy_x[j] * fl1_fx + 1.5 * pa2pb_yyy_xx[j] * fl1_fx + 1.5 * pa2pb_xy_xxx[j] * fl1_fx + pa2pb_xyyy_xxx[j]);

                t_xyyy_xxx[j] += fl_r_0_0 * (-2.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl3_fx * fl1_fz + 7.5 * pa_yyy[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyy_x[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyyy_x[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yyy_xx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_xxx[j] * fl1_fz);

                t_xyyy_xxy[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pb_x[j] * fl3_fx + 0.75 * pa_xyy[j] * fl2_fx + 1.5 * pa2pb_yy_x[j] * fl2_fx + 0.75 * pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pa2pb_y_xy[j] * fl2_fx + 0.5 * pa2pb_xyyy_y[j] * fl1_fx + 1.5 * pa2pb_xyy_xx[j] * fl1_fx + pa2pb_yyy_xy[j] * fl1_fx + 1.5 * pa2pb_xy_xxy[j] * fl1_fx + pa2pb_xyyy_xxy[j]);

                t_xyyy_xxy[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 6.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yy_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyy_y[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xyy_xx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_xxy[j] * fl1_fz);

                t_xyyy_xxz[j] = fl_s_0_0 * (0.75 * pa2pb_xy_z[j] * fl2_fx + 1.5 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pa2pb_xyyy_z[j] * fl1_fx + pa2pb_yyy_xz[j] * fl1_fx + 1.5 * pa2pb_xy_xxz[j] * fl1_fx + pa2pb_xyyy_xxz[j]);

                t_xyyy_xxz[j] += fl_r_0_0 * (-1.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_z[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyy_z[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_xxz[j] * fl1_fz);

                t_xyyy_xyy[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pb_y[j] * fl3_fx + 2.25 * pa2pb_xy_x[j] * fl2_fx + 0.25 * pa_yyy[j] * fl2_fx + 1.5 * pa2pb_yy_y[j] * fl2_fx + 1.5 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pa2pb_y_yy[j] * fl2_fx + 0.5 * pa2pb_xyyy_x[j] * fl1_fx + 3.0 * pa2pb_xyy_xy[j] * fl1_fx + 0.5 * pa2pb_yyy_yy[j] * fl1_fx + 1.5 * pa2pb_xy_xyy[j] * fl1_fx + pa2pb_xyyy_xyy[j]);

                t_xyyy_xyy[j] += fl_r_0_0 * (9.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_y[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 2.5 * pa_yyy[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_yy_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyy_x[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xyy_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyy_yy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_xyy[j] * fl1_fz);

                t_xyyy_xyz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 0.75 * pa2pb_yy_z[j] * fl2_fx + 0.75 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pa2pb_y_yz[j] * fl2_fx + 1.5 * pa2pb_xyy_xz[j] * fl1_fx + 0.5 * pa2pb_yyy_yz[j] * fl1_fx + 1.5 * pa2pb_xy_xyz[j] * fl1_fx + pa2pb_xyyy_xyz[j]);

                t_xyyy_xyz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fga + 7.5 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyy_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyy_yz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_65_70(      CMemBlock2D<double>& primBuffer,
                                 const CMemBlock2D<double>& auxBuffer,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_xy_xzz = pa2pbDistances.data(646 * idx + 90);

            auto pa2pb_xy_yyy = pa2pbDistances.data(646 * idx + 91);

            auto pa2pb_xy_yyz = pa2pbDistances.data(646 * idx + 92);

            auto pa2pb_xy_yzz = pa2pbDistances.data(646 * idx + 93);

            auto pa2pb_xy_zzz = pa2pbDistances.data(646 * idx + 94);

            auto pa2pb_xyy_yy = pa2pbDistances.data(646 * idx + 234);

            auto pa2pb_xyy_yz = pa2pbDistances.data(646 * idx + 235);

            auto pa2pb_xyy_zz = pa2pbDistances.data(646 * idx + 236);

            auto pa2pb_yyy_zz = pa2pbDistances.data(646 * idx + 293);

            auto pa2pb_xyyy_x = pa2pbDistances.data(646 * idx + 475);

            auto pa2pb_xyyy_y = pa2pbDistances.data(646 * idx + 476);

            auto pa2pb_xyyy_z = pa2pbDistances.data(646 * idx + 477);

            auto pa2pb_xyyy_xzz = pa2pbDistances.data(646 * idx + 489);

            auto pa2pb_xyyy_yyy = pa2pbDistances.data(646 * idx + 490);

            auto pa2pb_xyyy_yyz = pa2pbDistances.data(646 * idx + 491);

            auto pa2pb_xyyy_yzz = pa2pbDistances.data(646 * idx + 492);

            auto pa2pb_xyyy_zzz = pa2pbDistances.data(646 * idx + 493);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyy_xzz = primBuffer.data(150 * idx + 65);

            auto t_xyyy_yyy = primBuffer.data(150 * idx + 66);

            auto t_xyyy_yyz = primBuffer.data(150 * idx + 67);

            auto t_xyyy_yzz = primBuffer.data(150 * idx + 68);

            auto t_xyyy_zzz = primBuffer.data(150 * idx + 69);

            // Batch of Integrals (65,70)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xy_x, \
                                     pa2pb_xy_xzz, pa2pb_xy_y, pa2pb_xy_yyy, pa2pb_xy_yyz, pa2pb_xy_yzz, pa2pb_xy_z, \
                                     pa2pb_xy_zzz, pa2pb_xyy_yy, pa2pb_xyy_yz, pa2pb_xyy_zz, pa2pb_xyyy_x, \
                                     pa2pb_xyyy_xzz, pa2pb_xyyy_y, pa2pb_xyyy_yyy, pa2pb_xyyy_yyz, pa2pb_xyyy_yzz, \
                                     pa2pb_xyyy_z, pa2pb_xyyy_zzz, pa2pb_y_zz, pa2pb_yyy_zz, pa_x, pa_xyy, pa_y, pa_yyy, \
                                     r_0_0, s_0_0, t_xyyy_xzz, t_xyyy_yyy, t_xyyy_yyz, t_xyyy_yzz, t_xyyy_zzz: VLX_ALIGN)
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

                t_xyyy_xzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.25 * pa_yyy[j] * fl2_fx + 0.75 * pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pa2pb_xyyy_x[j] * fl1_fx + 0.5 * pa2pb_yyy_zz[j] * fl1_fx + 1.5 * pa2pb_xy_xzz[j] * fl1_fx + pa2pb_xyyy_xzz[j]);

                t_xyyy_xzz[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 2.5 * pa_yyy[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyy_x[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyy_zz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xy_xzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_xzz[j] * fl1_fz);

                t_xyyy_yyy[j] = fl_s_0_0 * (1.875 * pa_x[j] * fl3_fx + 2.25 * pa_xyy[j] * fl2_fx + 6.75 * pa2pb_xy_y[j] * fl2_fx + 2.25 * pa2pb_x_yy[j] * fl2_fx + 1.5 * pa2pb_xyyy_y[j] * fl1_fx + 4.5 * pa2pb_xyy_yy[j] * fl1_fx + 1.5 * pa2pb_xy_yyy[j] * fl1_fx + pa2pb_xyyy_yyy[j]);

                t_xyyy_yyy[j] += fl_r_0_0 * (15.0 * pa_x[j] * fl3_fx * fl1_fz - 2.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_xyy[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_xy_y[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyy_y[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyyy_y[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_xyy_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_yyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_yyy[j] * fl1_fz);

                t_xyyy_yyz[j] = fl_s_0_0 * (2.25 * pa2pb_xy_z[j] * fl2_fx + 1.5 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_xyyy_z[j] * fl1_fx + 3.0 * pa2pb_xyy_yz[j] * fl1_fx + 1.5 * pa2pb_xy_yyz[j] * fl1_fx + pa2pb_xyyy_yyz[j]);

                t_xyyy_yyz[j] += fl_r_0_0 * (22.5 * pa2pb_xy_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyy_z[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xyy_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_yyz[j] * fl1_fz);

                t_xyyy_yzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xyy[j] * fl2_fx + 0.75 * pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_x_zz[j] * fl2_fx + 0.5 * pa2pb_xyyy_y[j] * fl1_fx + 1.5 * pa2pb_xyy_zz[j] * fl1_fx + 1.5 * pa2pb_xy_yzz[j] * fl1_fx + pa2pb_xyyy_yzz[j]);

                t_xyyy_yzz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyy_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyy_y[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xyy_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_yzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_yzz[j] * fl1_fz);

                t_xyyy_zzz[j] = fl_s_0_0 * (2.25 * pa2pb_xy_z[j] * fl2_fx + 1.5 * pa2pb_xyyy_z[j] * fl1_fx + 1.5 * pa2pb_xy_zzz[j] * fl1_fx + pa2pb_xyyy_zzz[j]);

                t_xyyy_zzz[j] += fl_r_0_0 * (-4.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyyy_z[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_xy_z[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xyyy_z[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xy_zzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xy_zzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyyy_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_70_75(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_xz_xxx = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xz_xxy = pa2pbDistances.data(646 * idx + 105);

            auto pa2pb_xz_xxz = pa2pbDistances.data(646 * idx + 106);

            auto pa2pb_xz_xyy = pa2pbDistances.data(646 * idx + 107);

            auto pa2pb_xz_xyz = pa2pbDistances.data(646 * idx + 108);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_xyy_xx = pa2pbDistances.data(646 * idx + 231);

            auto pa2pb_xyy_xy = pa2pbDistances.data(646 * idx + 232);

            auto pa2pb_xyz_xx = pa2pbDistances.data(646 * idx + 250);

            auto pa2pb_xyz_xy = pa2pbDistances.data(646 * idx + 251);

            auto pa2pb_xyz_xz = pa2pbDistances.data(646 * idx + 252);

            auto pa2pb_yyz_xx = pa2pbDistances.data(646 * idx + 307);

            auto pa2pb_yyz_xy = pa2pbDistances.data(646 * idx + 308);

            auto pa2pb_yyz_xz = pa2pbDistances.data(646 * idx + 309);

            auto pa2pb_yyz_yy = pa2pbDistances.data(646 * idx + 310);

            auto pa2pb_yyz_yz = pa2pbDistances.data(646 * idx + 311);

            auto pa2pb_xyyz_x = pa2pbDistances.data(646 * idx + 494);

            auto pa2pb_xyyz_y = pa2pbDistances.data(646 * idx + 495);

            auto pa2pb_xyyz_z = pa2pbDistances.data(646 * idx + 496);

            auto pa2pb_xyyz_xxx = pa2pbDistances.data(646 * idx + 503);

            auto pa2pb_xyyz_xxy = pa2pbDistances.data(646 * idx + 504);

            auto pa2pb_xyyz_xxz = pa2pbDistances.data(646 * idx + 505);

            auto pa2pb_xyyz_xyy = pa2pbDistances.data(646 * idx + 506);

            auto pa2pb_xyyz_xyz = pa2pbDistances.data(646 * idx + 507);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyz_xxx = primBuffer.data(150 * idx + 70);

            auto t_xyyz_xxy = primBuffer.data(150 * idx + 71);

            auto t_xyyz_xxz = primBuffer.data(150 * idx + 72);

            auto t_xyyz_xyy = primBuffer.data(150 * idx + 73);

            auto t_xyyz_xyz = primBuffer.data(150 * idx + 74);

            // Batch of Integrals (70,75)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_xy_x, pa2pb_xyy_xx, \
                                     pa2pb_xyy_xy, pa2pb_xyyz_x, pa2pb_xyyz_xxx, pa2pb_xyyz_xxy, pa2pb_xyyz_xxz, \
                                     pa2pb_xyyz_xyy, pa2pb_xyyz_xyz, pa2pb_xyyz_y, pa2pb_xyyz_z, pa2pb_xyz_xx, \
                                     pa2pb_xyz_xy, pa2pb_xyz_xz, pa2pb_xz_x, pa2pb_xz_xxx, pa2pb_xz_xxy, pa2pb_xz_xxz, \
                                     pa2pb_xz_xyy, pa2pb_xz_xyz, pa2pb_xz_y, pa2pb_xz_z, pa2pb_yy_x, pa2pb_yy_y, \
                                     pa2pb_yyz_xx, pa2pb_yyz_xy, pa2pb_yyz_xz, pa2pb_yyz_yy, pa2pb_yyz_yz, pa2pb_yz_x, \
                                     pa2pb_yz_y, pa2pb_yz_z, pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, \
                                     pa_x, pa_xyy, pa_xyz, pa_y, pa_yyz, pa_z, pb_x, pb_y, r_0_0, s_0_0, t_xyyz_xxx, \
                                     t_xyyz_xxy, t_xyyz_xxz, t_xyyz_xyy, t_xyyz_xyz: VLX_ALIGN)
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

                t_xyyz_xxx[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa_yyz[j] * fl2_fx + 0.75 * pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_z_xx[j] * fl2_fx + 1.5 * pa2pb_xyyz_x[j] * fl1_fx + 1.5 * pa2pb_yyz_xx[j] * fl1_fx + 0.5 * pa2pb_xz_xxx[j] * fl1_fx + pa2pb_xyyz_xxx[j]);

                t_xyyz_xxx[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 7.5 * pa_yyz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyyz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyyz_x[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yyz_xx[j] * fl1_fx * fl1_fz - pa2pb_xz_xxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_xxx[j] * fl1_fz);

                t_xyyz_xxy[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl2_fx + pa2pb_yz_x[j] * fl2_fx + 0.25 * pa2pb_xz_y[j] * fl2_fx + 0.5 * pa2pb_z_xy[j] * fl2_fx + 0.5 * pa2pb_xyyz_y[j] * fl1_fx + pa2pb_xyz_xx[j] * fl1_fx + pa2pb_yyz_xy[j] * fl1_fx + 0.5 * pa2pb_xz_xxy[j] * fl1_fx + pa2pb_xyyz_xxy[j]);

                t_xyyz_xxy[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xz_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyz_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyz_xx[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yyz_xy[j] * fl1_fx * fl1_fz - pa2pb_xz_xxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_xxy[j] * fl1_fz);

                t_xyyz_xxz[j] = fl_s_0_0 * (0.125 * pa_x[j] * fl3_fx + 0.25 * pb_x[j] * fl3_fx + 0.25 * pa_xyy[j] * fl2_fx + 0.5 * pa2pb_yy_x[j] * fl2_fx + 0.25 * pa2pb_xz_z[j] * fl2_fx + 0.25 * pa2pb_x_xx[j] * fl2_fx + 0.5 * pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_xyyz_z[j] * fl1_fx + 0.5 * pa2pb_xyy_xx[j] * fl1_fx + pa2pb_yyz_xz[j] * fl1_fx + 0.5 * pa2pb_xz_xxz[j] * fl1_fx + pa2pb_xyyz_xxz[j]);

                t_xyyz_xxz[j] += fl_r_0_0 * (-0.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + pa_x[j] * fl3_fx * fl1_fz + 2.0 * pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xyy[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_yy_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyz_z[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xyy_xx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyz_xz[j] * fl1_fx * fl1_fz - pa2pb_xz_xxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_xxz[j] * fl1_fz);

                t_xyyz_xyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa2pb_xz_x[j] * fl2_fx + 0.25 * pa_yyz[j] * fl2_fx + pa2pb_yz_y[j] * fl2_fx + 0.25 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pa2pb_xyyz_x[j] * fl1_fx + 2.0 * pa2pb_xyz_xy[j] * fl1_fx + 0.5 * pa2pb_yyz_yy[j] * fl1_fx + 0.5 * pa2pb_xz_xyy[j] * fl1_fx + pa2pb_xyyz_xyy[j]);

                t_xyyz_xyy[j] += fl_r_0_0 * (3.0 * pa_z[j] * fl3_fx * fl1_fz - 0.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 2.5 * pa_yyz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_yz_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyz_x[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xyz_xy[j] * fl1_fx * fl1_fz + 6.0 * pa2pb_yyz_yy[j] * fl1_fx * fl1_fz - pa2pb_xz_xyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_xyy[j] * fl1_fz);

                t_xyyz_xyz[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl3_fx + 0.125 * pb_y[j] * fl3_fx + 0.5 * pa2pb_xy_x[j] * fl2_fx + 0.25 * pa2pb_yy_y[j] * fl2_fx + 0.5 * pa2pb_yz_z[j] * fl2_fx + 0.25 * pa2pb_x_xy[j] * fl2_fx + 0.25 * pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_xyy_xy[j] * fl1_fx + pa2pb_xyz_xz[j] * fl1_fx + 0.5 * pa2pb_yyz_yz[j] * fl1_fx + 0.5 * pa2pb_xz_xyz[j] * fl1_fx + pa2pb_xyyz_xyz[j]);

                t_xyyz_xyz[j] += fl_r_0_0 * (2.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga + pb_y[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_yy_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yz_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xy[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fga + 2.5 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyy_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyz_xz[j] * fl1_fx * fl1_fz + 6.0 * pa2pb_yyz_yz[j] * fl1_fx * fl1_fz - pa2pb_xz_xyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_75_80(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_xz_xzz = pa2pbDistances.data(646 * idx + 109);

            auto pa2pb_xz_yyy = pa2pbDistances.data(646 * idx + 110);

            auto pa2pb_xz_yyz = pa2pbDistances.data(646 * idx + 111);

            auto pa2pb_xz_yzz = pa2pbDistances.data(646 * idx + 112);

            auto pa2pb_xz_zzz = pa2pbDistances.data(646 * idx + 113);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_xyy_xz = pa2pbDistances.data(646 * idx + 233);

            auto pa2pb_xyy_yy = pa2pbDistances.data(646 * idx + 234);

            auto pa2pb_xyy_yz = pa2pbDistances.data(646 * idx + 235);

            auto pa2pb_xyy_zz = pa2pbDistances.data(646 * idx + 236);

            auto pa2pb_xyz_yy = pa2pbDistances.data(646 * idx + 253);

            auto pa2pb_xyz_yz = pa2pbDistances.data(646 * idx + 254);

            auto pa2pb_xyz_zz = pa2pbDistances.data(646 * idx + 255);

            auto pa2pb_yyz_zz = pa2pbDistances.data(646 * idx + 312);

            auto pa2pb_xyyz_x = pa2pbDistances.data(646 * idx + 494);

            auto pa2pb_xyyz_y = pa2pbDistances.data(646 * idx + 495);

            auto pa2pb_xyyz_z = pa2pbDistances.data(646 * idx + 496);

            auto pa2pb_xyyz_xzz = pa2pbDistances.data(646 * idx + 508);

            auto pa2pb_xyyz_yyy = pa2pbDistances.data(646 * idx + 509);

            auto pa2pb_xyyz_yyz = pa2pbDistances.data(646 * idx + 510);

            auto pa2pb_xyyz_yzz = pa2pbDistances.data(646 * idx + 511);

            auto pa2pb_xyyz_zzz = pa2pbDistances.data(646 * idx + 512);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyyz_xzz = primBuffer.data(150 * idx + 75);

            auto t_xyyz_yyy = primBuffer.data(150 * idx + 76);

            auto t_xyyz_yyz = primBuffer.data(150 * idx + 77);

            auto t_xyyz_yzz = primBuffer.data(150 * idx + 78);

            auto t_xyyz_zzz = primBuffer.data(150 * idx + 79);

            // Batch of Integrals (75,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xy_y, \
                                     pa2pb_xy_z, pa2pb_xyy_xz, pa2pb_xyy_yy, pa2pb_xyy_yz, pa2pb_xyy_zz, pa2pb_xyyz_x, \
                                     pa2pb_xyyz_xzz, pa2pb_xyyz_y, pa2pb_xyyz_yyy, pa2pb_xyyz_yyz, pa2pb_xyyz_yzz, \
                                     pa2pb_xyyz_z, pa2pb_xyyz_zzz, pa2pb_xyz_yy, pa2pb_xyz_yz, pa2pb_xyz_zz, pa2pb_xz_x, \
                                     pa2pb_xz_xzz, pa2pb_xz_y, pa2pb_xz_yyy, pa2pb_xz_yyz, pa2pb_xz_yzz, pa2pb_xz_z, \
                                     pa2pb_xz_zzz, pa2pb_yy_z, pa2pb_yyz_zz, pa2pb_z_zz, pa_x, pa_xyy, pa_xyz, pa_yyz, pa_z, \
                                     pb_z, r_0_0, s_0_0, t_xyyz_xzz, t_xyyz_yyy, t_xyyz_yyz, t_xyyz_yzz, t_xyyz_zzz: VLX_ALIGN)
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

                t_xyyz_xzz[j] = fl_s_0_0 * (0.125 * pa_z[j] * fl3_fx + 0.25 * pb_z[j] * fl3_fx + 0.25 * pa_yyz[j] * fl2_fx + 0.5 * pa2pb_yy_z[j] * fl2_fx + 0.25 * pa2pb_xz_x[j] * fl2_fx + 0.5 * pa2pb_x_xz[j] * fl2_fx + 0.25 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pa2pb_xyyz_x[j] * fl1_fx + pa2pb_xyy_xz[j] * fl1_fx + 0.5 * pa2pb_yyz_zz[j] * fl1_fx + 0.5 * pa2pb_xz_xzz[j] * fl1_fx + pa2pb_xyyz_xzz[j]);

                t_xyyz_xzz[j] += fl_r_0_0 * (-0.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + pa_z[j] * fl3_fx * fl1_fz + 2.0 * pb_z[j] * fl3_fx * fl1_fz + 2.5 * pa_yyz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yy_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_xz[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyyz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyz_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyz_zz[j] * fl1_fx * fl1_fz - pa2pb_xz_xzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_xzz[j] * fl1_fz);

                t_xyyz_yyy[j] = fl_s_0_0 * (1.5 * pa_xyz[j] * fl2_fx + 2.25 * pa2pb_xz_y[j] * fl2_fx + 1.5 * pa2pb_xyyz_y[j] * fl1_fx + 3.0 * pa2pb_xyz_yy[j] * fl1_fx + 0.5 * pa2pb_xz_yyy[j] * fl1_fx + pa2pb_xyyz_yyy[j]);

                t_xyyz_yyy[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa_xyz[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_xz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyyz_y[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyyz_y[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xyz_yy[j] * fl1_fx * fl1_fz - pa2pb_xz_yyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_yyy[j] * fl1_fz);

                t_xyyz_yyz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.25 * pa_xyy[j] * fl2_fx + pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_xz_z[j] * fl2_fx + 0.25 * pa2pb_x_yy[j] * fl2_fx + 0.5 * pa2pb_xyyz_z[j] * fl1_fx + 0.5 * pa2pb_xyy_yy[j] * fl1_fx + 2.0 * pa2pb_xyz_yz[j] * fl1_fx + 0.5 * pa2pb_xz_yyz[j] * fl1_fx + pa2pb_xyyz_yyz[j]);

                t_xyyz_yyz[j] += fl_r_0_0 * (3.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 2.5 * pa_xyy[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xy_y[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xyyz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyz_z[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xyy_yy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xyz_yz[j] * fl1_fx * fl1_fz - pa2pb_xz_yyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_yyz[j] * fl1_fz);

                t_xyyz_yzz[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl2_fx + pa2pb_xy_z[j] * fl2_fx + 0.25 * pa2pb_xz_y[j] * fl2_fx + 0.5 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_xyyz_y[j] * fl1_fx + pa2pb_xyy_yz[j] * fl1_fx + pa2pb_xyz_zz[j] * fl1_fx + 0.5 * pa2pb_xz_yzz[j] * fl1_fx + pa2pb_xyyz_yzz[j]);

                t_xyyz_yzz[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_xy_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_xyyz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xz_y[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyyz_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyy_yz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyz_zz[j] * fl1_fx * fl1_fz - pa2pb_xz_yzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_yzz[j] * fl1_fz);

                t_xyyz_zzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xyy[j] * fl2_fx + 0.75 * pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_zz[j] * fl2_fx + 1.5 * pa2pb_xyyz_z[j] * fl1_fx + 1.5 * pa2pb_xyy_zz[j] * fl1_fx + 0.5 * pa2pb_xz_zzz[j] * fl1_fx + pa2pb_xyyz_zzz[j]);

                t_xyyz_zzz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_xyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xyy[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyyz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyyz_z[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xyy_zz[j] * fl1_fz * fl1_fx - pa2pb_xz_zzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xz_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_xyyz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_80_85(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_xy_xxx = pa2pbDistances.data(646 * idx + 85);

            auto pa2pb_xy_xxy = pa2pbDistances.data(646 * idx + 86);

            auto pa2pb_xy_xxz = pa2pbDistances.data(646 * idx + 87);

            auto pa2pb_xy_xyy = pa2pbDistances.data(646 * idx + 88);

            auto pa2pb_xy_xyz = pa2pbDistances.data(646 * idx + 89);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_xyz_xx = pa2pbDistances.data(646 * idx + 250);

            auto pa2pb_xyz_xy = pa2pbDistances.data(646 * idx + 251);

            auto pa2pb_xzz_xx = pa2pbDistances.data(646 * idx + 269);

            auto pa2pb_xzz_xy = pa2pbDistances.data(646 * idx + 270);

            auto pa2pb_xzz_xz = pa2pbDistances.data(646 * idx + 271);

            auto pa2pb_yzz_xx = pa2pbDistances.data(646 * idx + 326);

            auto pa2pb_yzz_xy = pa2pbDistances.data(646 * idx + 327);

            auto pa2pb_yzz_xz = pa2pbDistances.data(646 * idx + 328);

            auto pa2pb_yzz_yy = pa2pbDistances.data(646 * idx + 329);

            auto pa2pb_yzz_yz = pa2pbDistances.data(646 * idx + 330);

            auto pa2pb_xyzz_x = pa2pbDistances.data(646 * idx + 513);

            auto pa2pb_xyzz_y = pa2pbDistances.data(646 * idx + 514);

            auto pa2pb_xyzz_z = pa2pbDistances.data(646 * idx + 515);

            auto pa2pb_xyzz_xxx = pa2pbDistances.data(646 * idx + 522);

            auto pa2pb_xyzz_xxy = pa2pbDistances.data(646 * idx + 523);

            auto pa2pb_xyzz_xxz = pa2pbDistances.data(646 * idx + 524);

            auto pa2pb_xyzz_xyy = pa2pbDistances.data(646 * idx + 525);

            auto pa2pb_xyzz_xyz = pa2pbDistances.data(646 * idx + 526);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyzz_xxx = primBuffer.data(150 * idx + 80);

            auto t_xyzz_xxy = primBuffer.data(150 * idx + 81);

            auto t_xyzz_xxz = primBuffer.data(150 * idx + 82);

            auto t_xyzz_xyy = primBuffer.data(150 * idx + 83);

            auto t_xyzz_xyz = primBuffer.data(150 * idx + 84);

            // Batch of Integrals (80,85)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_x_xz, pa2pb_xy_x, \
                                     pa2pb_xy_xxx, pa2pb_xy_xxy, pa2pb_xy_xxz, pa2pb_xy_xyy, pa2pb_xy_xyz, pa2pb_xy_y, \
                                     pa2pb_xy_z, pa2pb_xyz_xx, pa2pb_xyz_xy, pa2pb_xyzz_x, pa2pb_xyzz_xxx, \
                                     pa2pb_xyzz_xxy, pa2pb_xyzz_xxz, pa2pb_xyzz_xyy, pa2pb_xyzz_xyz, pa2pb_xyzz_y, \
                                     pa2pb_xyzz_z, pa2pb_xz_x, pa2pb_xzz_xx, pa2pb_xzz_xy, pa2pb_xzz_xz, pa2pb_y_xx, \
                                     pa2pb_y_xy, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_yz_x, pa2pb_yz_y, \
                                     pa2pb_yzz_xx, pa2pb_yzz_xy, pa2pb_yzz_xz, pa2pb_yzz_yy, pa2pb_yzz_yz, pa2pb_zz_x, \
                                     pa2pb_zz_y, pa2pb_zz_z, pa_x, pa_xyz, pa_xzz, pa_y, pa_yzz, pa_z, pb_x, pb_y, pb_z, r_0_0, \
                                     s_0_0, t_xyzz_xxx, t_xyzz_xxy, t_xyzz_xxz, t_xyzz_xyy, t_xyzz_xyz: VLX_ALIGN)
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

                t_xyzz_xxx[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa_yzz[j] * fl2_fx + 0.75 * pa2pb_xy_x[j] * fl2_fx + 0.75 * pa2pb_y_xx[j] * fl2_fx + 1.5 * pa2pb_xyzz_x[j] * fl1_fx + 1.5 * pa2pb_yzz_xx[j] * fl1_fx + 0.5 * pa2pb_xy_xxx[j] * fl1_fx + pa2pb_xyzz_xxx[j]);

                t_xyzz_xxx[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 7.5 * pa_yzz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyzz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyzz_x[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yzz_xx[j] * fl1_fx * fl1_fz - pa2pb_xy_xxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_xxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_xxx[j] * fl1_fz);

                t_xyzz_xxy[j] = fl_s_0_0 * (0.125 * pa_x[j] * fl3_fx + 0.25 * pb_x[j] * fl3_fx + 0.25 * pa_xzz[j] * fl2_fx + 0.5 * pa2pb_zz_x[j] * fl2_fx + 0.25 * pa2pb_xy_y[j] * fl2_fx + 0.25 * pa2pb_x_xx[j] * fl2_fx + 0.5 * pa2pb_y_xy[j] * fl2_fx + 0.5 * pa2pb_xyzz_y[j] * fl1_fx + 0.5 * pa2pb_xzz_xx[j] * fl1_fx + pa2pb_yzz_xy[j] * fl1_fx + 0.5 * pa2pb_xy_xxy[j] * fl1_fx + pa2pb_xyzz_xxy[j]);

                t_xyzz_xxy[j] += fl_r_0_0 * (-0.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + pa_x[j] * fl3_fx * fl1_fz + 2.0 * pb_x[j] * fl3_fx * fl1_fz + 2.5 * pa_xzz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyzz_y[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xzz_xx[j] * fl1_fx * fl1_fz + 12.0 * pa2pb_yzz_xy[j] * fl1_fx * fl1_fz - pa2pb_xy_xxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_xxy[j] * fl1_fz);

                t_xyzz_xxz[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl2_fx + pa2pb_yz_x[j] * fl2_fx + 0.25 * pa2pb_xy_z[j] * fl2_fx + 0.5 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pa2pb_xyzz_z[j] * fl1_fx + pa2pb_xyz_xx[j] * fl1_fx + pa2pb_yzz_xz[j] * fl1_fx + 0.5 * pa2pb_xy_xxz[j] * fl1_fx + pa2pb_xyzz_xxz[j]);

                t_xyzz_xxz[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xy_z[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyzz_z[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyz_xx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xz[j] * fl1_fx * fl1_fz - pa2pb_xy_xxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_xxz[j] * fl1_fz);

                t_xyzz_xyy[j] = fl_s_0_0 * (0.125 * pa_y[j] * fl3_fx + 0.25 * pb_y[j] * fl3_fx + 0.25 * pa_yzz[j] * fl2_fx + 0.5 * pa2pb_zz_y[j] * fl2_fx + 0.25 * pa2pb_xy_x[j] * fl2_fx + 0.5 * pa2pb_x_xy[j] * fl2_fx + 0.25 * pa2pb_y_yy[j] * fl2_fx + 0.5 * pa2pb_xyzz_x[j] * fl1_fx + pa2pb_xzz_xy[j] * fl1_fx + 0.5 * pa2pb_yzz_yy[j] * fl1_fx + 0.5 * pa2pb_xy_xyy[j] * fl1_fx + pa2pb_xyzz_xyy[j]);

                t_xyzz_xyy[j] += fl_r_0_0 * (-0.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + pa_y[j] * fl3_fx * fl1_fz + 2.0 * pb_y[j] * fl3_fx * fl1_fz + 2.5 * pa_yzz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xy_x[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyzz_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_xy[j] * fl1_fx * fl1_fz + 6.0 * pa2pb_yzz_yy[j] * fl1_fx * fl1_fz - pa2pb_xy_xyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_xyy[j] * fl1_fz);

                t_xyzz_xyz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl3_fx + 0.125 * pb_z[j] * fl3_fx + 0.5 * pa2pb_xz_x[j] * fl2_fx + 0.5 * pa2pb_yz_y[j] * fl2_fx + 0.25 * pa2pb_zz_z[j] * fl2_fx + 0.25 * pa2pb_x_xz[j] * fl2_fx + 0.25 * pa2pb_y_yz[j] * fl2_fx + pa2pb_xyz_xy[j] * fl1_fx + 0.5 * pa2pb_xzz_xz[j] * fl1_fx + 0.5 * pa2pb_yzz_yz[j] * fl1_fx + 0.5 * pa2pb_xy_xyz[j] * fl1_fx + pa2pb_xyzz_xyz[j]);

                t_xyzz_xyz[j] += fl_r_0_0 * (2.0 * pa_z[j] * fl3_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga + pb_z[j] * fl3_fx * fl1_fz + 5.0 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fga + 2.5 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_xyz_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xzz_xz[j] * fl1_fx * fl1_fz + 6.0 * pa2pb_yzz_yz[j] * fl1_fx * fl1_fz - pa2pb_xy_xyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_85_90(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_xy_x = pa2pbDistances.data(646 * idx + 76);

            auto pa2pb_xy_y = pa2pbDistances.data(646 * idx + 77);

            auto pa2pb_xy_z = pa2pbDistances.data(646 * idx + 78);

            auto pa2pb_xy_xzz = pa2pbDistances.data(646 * idx + 90);

            auto pa2pb_xy_yyy = pa2pbDistances.data(646 * idx + 91);

            auto pa2pb_xy_yyz = pa2pbDistances.data(646 * idx + 92);

            auto pa2pb_xy_yzz = pa2pbDistances.data(646 * idx + 93);

            auto pa2pb_xy_zzz = pa2pbDistances.data(646 * idx + 94);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_xyz_xz = pa2pbDistances.data(646 * idx + 252);

            auto pa2pb_xyz_yy = pa2pbDistances.data(646 * idx + 253);

            auto pa2pb_xyz_yz = pa2pbDistances.data(646 * idx + 254);

            auto pa2pb_xyz_zz = pa2pbDistances.data(646 * idx + 255);

            auto pa2pb_xzz_yy = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_xzz_yz = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_xzz_zz = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_yzz_zz = pa2pbDistances.data(646 * idx + 331);

            auto pa2pb_xyzz_x = pa2pbDistances.data(646 * idx + 513);

            auto pa2pb_xyzz_y = pa2pbDistances.data(646 * idx + 514);

            auto pa2pb_xyzz_z = pa2pbDistances.data(646 * idx + 515);

            auto pa2pb_xyzz_xzz = pa2pbDistances.data(646 * idx + 527);

            auto pa2pb_xyzz_yyy = pa2pbDistances.data(646 * idx + 528);

            auto pa2pb_xyzz_yyz = pa2pbDistances.data(646 * idx + 529);

            auto pa2pb_xyzz_yzz = pa2pbDistances.data(646 * idx + 530);

            auto pa2pb_xyzz_zzz = pa2pbDistances.data(646 * idx + 531);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xyzz_xzz = primBuffer.data(150 * idx + 85);

            auto t_xyzz_yyy = primBuffer.data(150 * idx + 86);

            auto t_xyzz_yyz = primBuffer.data(150 * idx + 87);

            auto t_xyzz_yzz = primBuffer.data(150 * idx + 88);

            auto t_xyzz_zzz = primBuffer.data(150 * idx + 89);

            // Batch of Integrals (85,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xy_x, \
                                     pa2pb_xy_xzz, pa2pb_xy_y, pa2pb_xy_yyy, pa2pb_xy_yyz, pa2pb_xy_yzz, pa2pb_xy_z, \
                                     pa2pb_xy_zzz, pa2pb_xyz_xz, pa2pb_xyz_yy, pa2pb_xyz_yz, pa2pb_xyz_zz, pa2pb_xyzz_x, \
                                     pa2pb_xyzz_xzz, pa2pb_xyzz_y, pa2pb_xyzz_yyy, pa2pb_xyzz_yyz, pa2pb_xyzz_yzz, \
                                     pa2pb_xyzz_z, pa2pb_xyzz_zzz, pa2pb_xz_y, pa2pb_xz_z, pa2pb_xzz_yy, pa2pb_xzz_yz, \
                                     pa2pb_xzz_zz, pa2pb_y_zz, pa2pb_yz_z, pa2pb_yzz_zz, pa_x, pa_xyz, pa_xzz, pa_y, pa_yzz, \
                                     r_0_0, s_0_0, t_xyzz_xzz, t_xyzz_yyy, t_xyzz_yyz, t_xyzz_yzz, t_xyzz_zzz: VLX_ALIGN)
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

                t_xyzz_xzz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa2pb_xy_x[j] * fl2_fx + 0.25 * pa_yzz[j] * fl2_fx + pa2pb_yz_z[j] * fl2_fx + 0.25 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pa2pb_xyzz_x[j] * fl1_fx + 2.0 * pa2pb_xyz_xz[j] * fl1_fx + 0.5 * pa2pb_yzz_zz[j] * fl1_fx + 0.5 * pa2pb_xy_xzz[j] * fl1_fx + pa2pb_xyzz_xzz[j]);

                t_xyzz_xzz[j] += fl_r_0_0 * (3.0 * pa_y[j] * fl3_fx * fl1_fz - 0.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_x[j] * fl2_fx * fl1_fz + 2.5 * pa_yzz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_yz_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_x[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyzz_x[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xyz_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yzz_zz[j] * fl1_fx * fl1_fz - pa2pb_xy_xzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_xzz[j] * fl1_fz);

                t_xyzz_yyy[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xzz[j] * fl2_fx + 0.75 * pa2pb_xy_y[j] * fl2_fx + 0.75 * pa2pb_x_yy[j] * fl2_fx + 1.5 * pa2pb_xyzz_y[j] * fl1_fx + 1.5 * pa2pb_xzz_yy[j] * fl1_fx + 0.5 * pa2pb_xy_yyy[j] * fl1_fx + pa2pb_xyzz_yyy[j]);

                t_xyzz_yyy[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xzz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xyzz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xyzz_y[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xzz_yy[j] * fl1_fx * fl1_fz - pa2pb_xy_yyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_yyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_yyy[j] * fl1_fz);

                t_xyzz_yyz[j] = fl_s_0_0 * (0.5 * pa_xyz[j] * fl2_fx + pa2pb_xz_y[j] * fl2_fx + 0.25 * pa2pb_xy_z[j] * fl2_fx + 0.5 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_xyzz_z[j] * fl1_fx + pa2pb_xyz_yy[j] * fl1_fx + pa2pb_xzz_yz[j] * fl1_fx + 0.5 * pa2pb_xy_yyz[j] * fl1_fx + pa2pb_xyzz_yyz[j]);

                t_xyzz_yyz[j] += fl_r_0_0 * (-pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 5.0 * pa_xyz[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_xz_y[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_xy_z[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyzz_z[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xyz_yy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_xzz_yz[j] * fl1_fx * fl1_fz - pa2pb_xy_yyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_yyz[j] * fl1_fz);

                t_xyzz_yzz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa2pb_xy_y[j] * fl2_fx + 0.25 * pa_xzz[j] * fl2_fx + pa2pb_xz_z[j] * fl2_fx + 0.25 * pa2pb_x_zz[j] * fl2_fx + 0.5 * pa2pb_xyzz_y[j] * fl1_fx + 2.0 * pa2pb_xyz_yz[j] * fl1_fx + 0.5 * pa2pb_xzz_zz[j] * fl1_fx + 0.5 * pa2pb_xy_yzz[j] * fl1_fx + pa2pb_xyzz_yzz[j]);

                t_xyzz_yzz[j] += fl_r_0_0 * (3.0 * pa_x[j] * fl3_fx * fl1_fz - 0.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 7.5 * pa2pb_xy_y[j] * fl2_fx * fl1_fz + 2.5 * pa_xzz[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_xz_z[j] * fl2_fx * fl1_fz - 0.5 * pa2pb_xy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_xy_y[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xyzz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xyzz_y[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_xyz_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_xzz_zz[j] * fl1_fx * fl1_fz - pa2pb_xy_yzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_yzz[j] * fl1_fz);

                t_xyzz_zzz[j] = fl_s_0_0 * (1.5 * pa_xyz[j] * fl2_fx + 2.25 * pa2pb_xy_z[j] * fl2_fx + 1.5 * pa2pb_xyzz_z[j] * fl1_fx + 3.0 * pa2pb_xyz_zz[j] * fl1_fx + 0.5 * pa2pb_xy_zzz[j] * fl1_fx + pa2pb_xyzz_zzz[j]);

                t_xyzz_zzz[j] += fl_r_0_0 * (-3.0 * pa_xyz[j] * fl1_fx * fl1_fz * fl1_fgb + 15.0 * pa_xyz[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_xy_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xy_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xyzz_z[j] * fl1_fz * fl1_fgb + 18.0 * pa2pb_xyzz_z[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xyz_zz[j] * fl1_fz * fl1_fx - pa2pb_xy_zzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_xy_zzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xyzz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_90_95(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xx = pa2pbDistances.data(646 * idx + 3);

            auto pa2pb_x_xy = pa2pbDistances.data(646 * idx + 4);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_xz_xxx = pa2pbDistances.data(646 * idx + 104);

            auto pa2pb_xz_xxy = pa2pbDistances.data(646 * idx + 105);

            auto pa2pb_xz_xxz = pa2pbDistances.data(646 * idx + 106);

            auto pa2pb_xz_xyy = pa2pbDistances.data(646 * idx + 107);

            auto pa2pb_xz_xyz = pa2pbDistances.data(646 * idx + 108);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_xzz_xx = pa2pbDistances.data(646 * idx + 269);

            auto pa2pb_xzz_xy = pa2pbDistances.data(646 * idx + 270);

            auto pa2pb_zzz_xx = pa2pbDistances.data(646 * idx + 345);

            auto pa2pb_zzz_xy = pa2pbDistances.data(646 * idx + 346);

            auto pa2pb_zzz_xz = pa2pbDistances.data(646 * idx + 347);

            auto pa2pb_zzz_yy = pa2pbDistances.data(646 * idx + 348);

            auto pa2pb_zzz_yz = pa2pbDistances.data(646 * idx + 349);

            auto pa2pb_xzzz_x = pa2pbDistances.data(646 * idx + 532);

            auto pa2pb_xzzz_y = pa2pbDistances.data(646 * idx + 533);

            auto pa2pb_xzzz_z = pa2pbDistances.data(646 * idx + 534);

            auto pa2pb_xzzz_xxx = pa2pbDistances.data(646 * idx + 541);

            auto pa2pb_xzzz_xxy = pa2pbDistances.data(646 * idx + 542);

            auto pa2pb_xzzz_xxz = pa2pbDistances.data(646 * idx + 543);

            auto pa2pb_xzzz_xyy = pa2pbDistances.data(646 * idx + 544);

            auto pa2pb_xzzz_xyz = pa2pbDistances.data(646 * idx + 545);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzzz_xxx = primBuffer.data(150 * idx + 90);

            auto t_xzzz_xxy = primBuffer.data(150 * idx + 91);

            auto t_xzzz_xxz = primBuffer.data(150 * idx + 92);

            auto t_xzzz_xyy = primBuffer.data(150 * idx + 93);

            auto t_xzzz_xyz = primBuffer.data(150 * idx + 94);

            // Batch of Integrals (90,95)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xx, pa2pb_x_xy, pa2pb_xz_x, pa2pb_xz_xxx, \
                                     pa2pb_xz_xxy, pa2pb_xz_xxz, pa2pb_xz_xyy, pa2pb_xz_xyz, pa2pb_xz_y, pa2pb_xz_z, \
                                     pa2pb_xzz_xx, pa2pb_xzz_xy, pa2pb_xzzz_x, pa2pb_xzzz_xxx, pa2pb_xzzz_xxy, \
                                     pa2pb_xzzz_xxz, pa2pb_xzzz_xyy, pa2pb_xzzz_xyz, pa2pb_xzzz_y, pa2pb_xzzz_z, \
                                     pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_zz_x, pa2pb_zz_y, \
                                     pa2pb_zzz_xx, pa2pb_zzz_xy, pa2pb_zzz_xz, pa2pb_zzz_yy, pa2pb_zzz_yz, pa_x, pa_xzz, \
                                     pa_z, pa_zzz, pb_x, pb_y, r_0_0, s_0_0, t_xzzz_xxx, t_xzzz_xxy, t_xzzz_xxz, \
                                     t_xzzz_xyy, t_xzzz_xyz: VLX_ALIGN)
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

                t_xzzz_xxx[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pa_zzz[j] * fl2_fx + 2.25 * pa2pb_xz_x[j] * fl2_fx + 2.25 * pa2pb_z_xx[j] * fl2_fx + 1.5 * pa2pb_xzzz_x[j] * fl1_fx + 1.5 * pa2pb_zzz_xx[j] * fl1_fx + 1.5 * pa2pb_xz_xxx[j] * fl1_fx + pa2pb_xzzz_xxx[j]);

                t_xzzz_xxx[j] += fl_r_0_0 * (-2.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl3_fx * fl1_fz + 7.5 * pa_zzz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzzz_x[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_xz_x[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xzzz_x[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zzz_xx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_xxx[j] * fl1_fz);

                t_xzzz_xxy[j] = fl_s_0_0 * (0.75 * pa2pb_xz_y[j] * fl2_fx + 1.5 * pa2pb_z_xy[j] * fl2_fx + 0.5 * pa2pb_xzzz_y[j] * fl1_fx + pa2pb_zzz_xy[j] * fl1_fx + 1.5 * pa2pb_xz_xxy[j] * fl1_fx + pa2pb_xzzz_xxy[j]);

                t_xzzz_xxy[j] += fl_r_0_0 * (-1.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_y[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzzz_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_xxy[j] * fl1_fz);

                t_xzzz_xxz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pb_x[j] * fl3_fx + 0.75 * pa_xzz[j] * fl2_fx + 1.5 * pa2pb_zz_x[j] * fl2_fx + 0.75 * pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_xx[j] * fl2_fx + 1.5 * pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_xzzz_z[j] * fl1_fx + 1.5 * pa2pb_xzz_xx[j] * fl1_fx + pa2pb_zzz_xz[j] * fl1_fx + 1.5 * pa2pb_xz_xxz[j] * fl1_fx + pa2pb_xzzz_xxz[j]);

                t_xzzz_xxz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 6.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xzz[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_xx[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_xx[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzzz_z[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xzz_xx[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_xxz[j] * fl1_fz);

                t_xzzz_xyy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.25 * pa_zzz[j] * fl2_fx + 0.75 * pa2pb_xz_x[j] * fl2_fx + 0.75 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pa2pb_xzzz_x[j] * fl1_fx + 0.5 * pa2pb_zzz_yy[j] * fl1_fx + 1.5 * pa2pb_xz_xyy[j] * fl1_fx + pa2pb_xzzz_xyy[j]);

                t_xzzz_xyy[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 2.5 * pa_zzz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_x[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzzz_x[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_yy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_xyy[j] * fl1_fz);

                t_xzzz_xyz[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_zz_y[j] * fl2_fx + 0.75 * pa2pb_x_xy[j] * fl2_fx + 0.75 * pa2pb_z_yz[j] * fl2_fx + 1.5 * pa2pb_xzz_xy[j] * fl1_fx + 0.5 * pa2pb_zzz_yz[j] * fl1_fx + 1.5 * pa2pb_xz_xyz[j] * fl1_fx + pa2pb_xzzz_xyz[j]);

                t_xzzz_xyz[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_x_xy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fga + 7.5 * pa2pb_x_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xzz_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_yz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_95_100(      CMemBlock2D<double>& primBuffer,
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

            auto pa_x = paDistances.data(34 * idx);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_x_xz = pa2pbDistances.data(646 * idx + 5);

            auto pa2pb_x_yy = pa2pbDistances.data(646 * idx + 6);

            auto pa2pb_x_yz = pa2pbDistances.data(646 * idx + 7);

            auto pa2pb_x_zz = pa2pbDistances.data(646 * idx + 8);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_xz_x = pa2pbDistances.data(646 * idx + 95);

            auto pa2pb_xz_y = pa2pbDistances.data(646 * idx + 96);

            auto pa2pb_xz_z = pa2pbDistances.data(646 * idx + 97);

            auto pa2pb_xz_xzz = pa2pbDistances.data(646 * idx + 109);

            auto pa2pb_xz_yyy = pa2pbDistances.data(646 * idx + 110);

            auto pa2pb_xz_yyz = pa2pbDistances.data(646 * idx + 111);

            auto pa2pb_xz_yzz = pa2pbDistances.data(646 * idx + 112);

            auto pa2pb_xz_zzz = pa2pbDistances.data(646 * idx + 113);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_xzz_xz = pa2pbDistances.data(646 * idx + 271);

            auto pa2pb_xzz_yy = pa2pbDistances.data(646 * idx + 272);

            auto pa2pb_xzz_yz = pa2pbDistances.data(646 * idx + 273);

            auto pa2pb_xzz_zz = pa2pbDistances.data(646 * idx + 274);

            auto pa2pb_zzz_zz = pa2pbDistances.data(646 * idx + 350);

            auto pa2pb_xzzz_x = pa2pbDistances.data(646 * idx + 532);

            auto pa2pb_xzzz_y = pa2pbDistances.data(646 * idx + 533);

            auto pa2pb_xzzz_z = pa2pbDistances.data(646 * idx + 534);

            auto pa2pb_xzzz_xzz = pa2pbDistances.data(646 * idx + 546);

            auto pa2pb_xzzz_yyy = pa2pbDistances.data(646 * idx + 547);

            auto pa2pb_xzzz_yyz = pa2pbDistances.data(646 * idx + 548);

            auto pa2pb_xzzz_yzz = pa2pbDistances.data(646 * idx + 549);

            auto pa2pb_xzzz_zzz = pa2pbDistances.data(646 * idx + 550);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_xzzz_xzz = primBuffer.data(150 * idx + 95);

            auto t_xzzz_yyy = primBuffer.data(150 * idx + 96);

            auto t_xzzz_yyz = primBuffer.data(150 * idx + 97);

            auto t_xzzz_yzz = primBuffer.data(150 * idx + 98);

            auto t_xzzz_zzz = primBuffer.data(150 * idx + 99);

            // Batch of Integrals (95,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_x_xz, pa2pb_x_yy, pa2pb_x_yz, pa2pb_x_zz, pa2pb_xz_x, \
                                     pa2pb_xz_xzz, pa2pb_xz_y, pa2pb_xz_yyy, pa2pb_xz_yyz, pa2pb_xz_yzz, pa2pb_xz_z, \
                                     pa2pb_xz_zzz, pa2pb_xzz_xz, pa2pb_xzz_yy, pa2pb_xzz_yz, pa2pb_xzz_zz, pa2pb_xzzz_x, \
                                     pa2pb_xzzz_xzz, pa2pb_xzzz_y, pa2pb_xzzz_yyy, pa2pb_xzzz_yyz, pa2pb_xzzz_yzz, \
                                     pa2pb_xzzz_z, pa2pb_xzzz_zzz, pa2pb_z_zz, pa2pb_zz_z, pa2pb_zzz_zz, pa_x, pa_xzz, pa_z, \
                                     pa_zzz, pb_z, r_0_0, s_0_0, t_xzzz_xzz, t_xzzz_yyy, t_xzzz_yyz, t_xzzz_yzz, \
                                     t_xzzz_zzz: VLX_ALIGN)
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

                t_xzzz_xzz[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pb_z[j] * fl3_fx + 2.25 * pa2pb_xz_x[j] * fl2_fx + 0.25 * pa_zzz[j] * fl2_fx + 1.5 * pa2pb_zz_z[j] * fl2_fx + 1.5 * pa2pb_x_xz[j] * fl2_fx + 0.75 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pa2pb_xzzz_x[j] * fl1_fx + 3.0 * pa2pb_xzz_xz[j] * fl1_fx + 0.5 * pa2pb_zzz_zz[j] * fl1_fx + 1.5 * pa2pb_xz_xzz[j] * fl1_fx + pa2pb_xzzz_xzz[j]);

                t_xzzz_xzz[j] += fl_r_0_0 * (9.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_z[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_xz_x[j] * fl2_fx * fl1_fz + 2.5 * pa_zzz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_zz_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_xz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_xz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzzz_x[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xzz_xz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_zz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_xz_xzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_xzz[j] * fl1_fz);

                t_xzzz_yyy[j] = fl_s_0_0 * (2.25 * pa2pb_xz_y[j] * fl2_fx + 1.5 * pa2pb_xzzz_y[j] * fl1_fx + 1.5 * pa2pb_xz_yyy[j] * fl1_fx + pa2pb_xzzz_yyy[j]);

                t_xzzz_yyy[j] += fl_r_0_0 * (-4.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_xzzz_y[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_xz_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_xzzz_y[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_yyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_yyy[j] * fl1_fz);

                t_xzzz_yyz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx + 0.75 * pa_xzz[j] * fl2_fx + 0.75 * pa2pb_xz_z[j] * fl2_fx + 0.75 * pa2pb_x_yy[j] * fl2_fx + 0.5 * pa2pb_xzzz_z[j] * fl1_fx + 1.5 * pa2pb_xzz_yy[j] * fl1_fx + 1.5 * pa2pb_xz_yyz[j] * fl1_fx + pa2pb_xzzz_yyz[j]);

                t_xzzz_yyz[j] += fl_r_0_0 * (-0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_x[j] * fl3_fx * fl1_fz + 7.5 * pa_xzz[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_x_yy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_xz_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_x_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzzz_z[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_xzz_yy[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_yyz[j] * fl1_fz);

                t_xzzz_yzz[j] = fl_s_0_0 * (2.25 * pa2pb_xz_y[j] * fl2_fx + 1.5 * pa2pb_x_yz[j] * fl2_fx + 0.5 * pa2pb_xzzz_y[j] * fl1_fx + 3.0 * pa2pb_xzz_yz[j] * fl1_fx + 1.5 * pa2pb_xz_yzz[j] * fl1_fx + pa2pb_xzzz_yzz[j]);

                t_xzzz_yzz[j] += fl_r_0_0 * (22.5 * pa2pb_xz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_xz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_xz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_x_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_xzzz_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_x_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_xzzz_y[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_xzz_yz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_yzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_yzz[j] * fl1_fz);

                t_xzzz_zzz[j] = fl_s_0_0 * (1.875 * pa_x[j] * fl3_fx + 2.25 * pa_xzz[j] * fl2_fx + 6.75 * pa2pb_xz_z[j] * fl2_fx + 2.25 * pa2pb_x_zz[j] * fl2_fx + 1.5 * pa2pb_xzzz_z[j] * fl1_fx + 4.5 * pa2pb_xzz_zz[j] * fl1_fx + 1.5 * pa2pb_xz_zzz[j] * fl1_fx + pa2pb_xzzz_zzz[j]);

                t_xzzz_zzz[j] += fl_r_0_0 * (15.0 * pa_x[j] * fl3_fx * fl1_fz - 2.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_x[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_xzz[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_xzz[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_xz_z[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_xz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_xz_z[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_x_zz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_xzzz_z[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_x_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_xzzz_z[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_xzz_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_xz_zzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_xz_zzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_xzzz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_100_105(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_yy_xxx = pa2pbDistances.data(646 * idx + 123);

            auto pa2pb_yy_xxy = pa2pbDistances.data(646 * idx + 124);

            auto pa2pb_yy_xxz = pa2pbDistances.data(646 * idx + 125);

            auto pa2pb_yy_xyy = pa2pbDistances.data(646 * idx + 126);

            auto pa2pb_yy_xyz = pa2pbDistances.data(646 * idx + 127);

            auto pa2pb_yyy_xx = pa2pbDistances.data(646 * idx + 288);

            auto pa2pb_yyy_xy = pa2pbDistances.data(646 * idx + 289);

            auto pa2pb_yyy_xz = pa2pbDistances.data(646 * idx + 290);

            auto pa2pb_yyyy_x = pa2pbDistances.data(646 * idx + 551);

            auto pa2pb_yyyy_y = pa2pbDistances.data(646 * idx + 552);

            auto pa2pb_yyyy_z = pa2pbDistances.data(646 * idx + 553);

            auto pa2pb_yyyy_xxx = pa2pbDistances.data(646 * idx + 560);

            auto pa2pb_yyyy_xxy = pa2pbDistances.data(646 * idx + 561);

            auto pa2pb_yyyy_xxz = pa2pbDistances.data(646 * idx + 562);

            auto pa2pb_yyyy_xyy = pa2pbDistances.data(646 * idx + 563);

            auto pa2pb_yyyy_xyz = pa2pbDistances.data(646 * idx + 564);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyy_xxx = primBuffer.data(150 * idx + 100);

            auto t_yyyy_xxy = primBuffer.data(150 * idx + 101);

            auto t_yyyy_xxz = primBuffer.data(150 * idx + 102);

            auto t_yyyy_xyy = primBuffer.data(150 * idx + 103);

            auto t_yyyy_xyz = primBuffer.data(150 * idx + 104);

            // Batch of Integrals (100,105)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_xz, pa2pb_yy_x, \
                                     pa2pb_yy_xxx, pa2pb_yy_xxy, pa2pb_yy_xxz, pa2pb_yy_xyy, pa2pb_yy_xyz, pa2pb_yy_y, \
                                     pa2pb_yy_z, pa2pb_yyy_xx, pa2pb_yyy_xy, pa2pb_yyy_xz, pa2pb_yyyy_x, \
                                     pa2pb_yyyy_xxx, pa2pb_yyyy_xxy, pa2pb_yyyy_xxz, pa2pb_yyyy_xyy, pa2pb_yyyy_xyz, \
                                     pa2pb_yyyy_y, pa2pb_yyyy_z, pa_y, pa_yyy, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_y, \
                                     pb_z, r_0_0, s_0_0, t_yyyy_xxx, t_yyyy_xxy, t_yyyy_xxz, t_yyyy_xyy, t_yyyy_xyz: VLX_ALIGN)
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

                t_yyyy_xxx[j] = fl_s_0_0 * (1.125 * pb_x[j] * fl3_fx + 4.5 * pa2pb_yy_x[j] * fl2_fx + 1.5 * pa2pb_yyyy_x[j] * fl1_fx + 0.75 * pb_xxx[j] * fl2_fx + 3.0 * pa2pb_yy_xxx[j] * fl1_fx + pa2pb_yyyy_xxx[j]);

                t_yyyy_xxx[j] += fl_r_0_0 * (-2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx + 9.0 * pb_x[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_yyyy_x[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_yy_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yyyy_x[j] * fl1_fz * fl1_fx - 3.0 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxx[j] * fl1_fz * fl1_fga + 7.5 * pb_xxx[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_xxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_xxx[j] * fl1_fz);

                t_yyyy_xxy[j] = fl_s_0_0 * (1.5 * pa_y[j] * fl3_fx + pa_yyy[j] * fl2_fx + 0.375 * pb_y[j] * fl3_fx + 1.5 * pa2pb_yy_y[j] * fl2_fx + 3.0 * pa2pb_y_xx[j] * fl2_fx + 0.5 * pa2pb_yyyy_y[j] * fl1_fx + 2.0 * pa2pb_yyy_xx[j] * fl1_fx + 0.75 * pb_xxy[j] * fl2_fx + 3.0 * pa2pb_yy_xxy[j] * fl1_fx + pa2pb_yyyy_xxy[j]);

                t_yyyy_xxy[j] += fl_r_0_0 * (-3.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_y[j] * fl3_fx * fl1_fz + 10.0 * pa_yyy[j] * fl1_fz * fl2_fx - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fga + 3.0 * pb_y[j] * fl3_fx * fl1_fz - pa2pb_yyyy_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyy_y[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yyy_xx[j] * fl1_fz * fl1_fx - 3.0 * pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxy[j] * fl1_fz * fl1_fga + 7.5 * pb_xxy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_xxy[j] * fl1_fz);

                t_yyyy_xxz[j] = fl_s_0_0 * (0.375 * pb_z[j] * fl3_fx + 1.5 * pa2pb_yy_z[j] * fl2_fx + 0.5 * pa2pb_yyyy_z[j] * fl1_fx + 0.75 * pb_xxz[j] * fl2_fx + 3.0 * pa2pb_yy_xxz[j] * fl1_fx + pa2pb_yyyy_xxz[j]);

                t_yyyy_xxz[j] += fl_r_0_0 * (-0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx + 3.0 * pb_z[j] * fl3_fx * fl1_fz - pa2pb_yyyy_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_yy_z[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yyyy_z[j] * fl1_fz * fl1_fx - 3.0 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xxz[j] * fl1_fz * fl1_fga + 7.5 * pb_xxz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_xxz[j] * fl1_fz);

                t_yyyy_xyy[j] = fl_s_0_0 * (1.875 * pb_x[j] * fl3_fx + 4.5 * pa2pb_yy_x[j] * fl2_fx + 6.0 * pa2pb_y_xy[j] * fl2_fx + 0.5 * pa2pb_yyyy_x[j] * fl1_fx + 4.0 * pa2pb_yyy_xy[j] * fl1_fx + 0.75 * pb_xyy[j] * fl2_fx + 3.0 * pa2pb_yy_xyy[j] * fl1_fx + pa2pb_yyyy_xyy[j]);

                t_yyyy_xyy[j] += fl_r_0_0 * (-4.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga + 15.0 * pb_x[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_yy_x[j] * fl2_fx * fl1_fz - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyy_x[j] * fl1_fz * fl1_fgb + 60.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyy_x[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_yyy_xy[j] * fl1_fz * fl1_fx - 3.0 * pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xyy[j] * fl1_fz * fl1_fga + 7.5 * pb_xyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_xyy[j] * fl1_fz);

                t_yyyy_xyz[j] = fl_s_0_0 * (3.0 * pa2pb_y_xz[j] * fl2_fx + 2.0 * pa2pb_yyy_xz[j] * fl1_fx + 0.75 * pb_xyz[j] * fl2_fx + 3.0 * pa2pb_yy_xyz[j] * fl1_fx + pa2pb_yyyy_xyz[j]);

                t_yyyy_xyz[j] += fl_r_0_0 * (-6.0 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fga + 30.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_yyy_xz[j] * fl1_fz * fl1_fx - 3.0 * pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xyz[j] * fl1_fz * fl1_fga + 7.5 * pb_xyz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_105_110(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_yy_xzz = pa2pbDistances.data(646 * idx + 128);

            auto pa2pb_yy_yyy = pa2pbDistances.data(646 * idx + 129);

            auto pa2pb_yy_yyz = pa2pbDistances.data(646 * idx + 130);

            auto pa2pb_yy_yzz = pa2pbDistances.data(646 * idx + 131);

            auto pa2pb_yy_zzz = pa2pbDistances.data(646 * idx + 132);

            auto pa2pb_yyy_yy = pa2pbDistances.data(646 * idx + 291);

            auto pa2pb_yyy_yz = pa2pbDistances.data(646 * idx + 292);

            auto pa2pb_yyy_zz = pa2pbDistances.data(646 * idx + 293);

            auto pa2pb_yyyy_x = pa2pbDistances.data(646 * idx + 551);

            auto pa2pb_yyyy_y = pa2pbDistances.data(646 * idx + 552);

            auto pa2pb_yyyy_z = pa2pbDistances.data(646 * idx + 553);

            auto pa2pb_yyyy_xzz = pa2pbDistances.data(646 * idx + 565);

            auto pa2pb_yyyy_yyy = pa2pbDistances.data(646 * idx + 566);

            auto pa2pb_yyyy_yyz = pa2pbDistances.data(646 * idx + 567);

            auto pa2pb_yyyy_yzz = pa2pbDistances.data(646 * idx + 568);

            auto pa2pb_yyyy_zzz = pa2pbDistances.data(646 * idx + 569);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyy_xzz = primBuffer.data(150 * idx + 105);

            auto t_yyyy_yyy = primBuffer.data(150 * idx + 106);

            auto t_yyyy_yyz = primBuffer.data(150 * idx + 107);

            auto t_yyyy_yzz = primBuffer.data(150 * idx + 108);

            auto t_yyyy_zzz = primBuffer.data(150 * idx + 109);

            // Batch of Integrals (105,110)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yy_x, \
                                     pa2pb_yy_xzz, pa2pb_yy_y, pa2pb_yy_yyy, pa2pb_yy_yyz, pa2pb_yy_yzz, pa2pb_yy_z, \
                                     pa2pb_yy_zzz, pa2pb_yyy_yy, pa2pb_yyy_yz, pa2pb_yyy_zz, pa2pb_yyyy_x, \
                                     pa2pb_yyyy_xzz, pa2pb_yyyy_y, pa2pb_yyyy_yyy, pa2pb_yyyy_yyz, pa2pb_yyyy_yzz, \
                                     pa2pb_yyyy_z, pa2pb_yyyy_zzz, pa_y, pa_yyy, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, \
                                     pb_zzz, r_0_0, s_0_0, t_yyyy_xzz, t_yyyy_yyy, t_yyyy_yyz, t_yyyy_yzz, t_yyyy_zzz: VLX_ALIGN)
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

                t_yyyy_xzz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 1.5 * pa2pb_yy_x[j] * fl2_fx + 0.5 * pa2pb_yyyy_x[j] * fl1_fx + 0.75 * pb_xzz[j] * fl2_fx + 3.0 * pa2pb_yy_xzz[j] * fl1_fx + pa2pb_yyyy_xzz[j]);

                t_yyyy_xzz[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx + 3.0 * pb_x[j] * fl3_fx * fl1_fz - pa2pb_yyyy_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_yy_x[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_yyyy_x[j] * fl1_fz * fl1_fx - 3.0 * pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_xzz[j] * fl1_fz * fl1_fga + 7.5 * pb_xzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_xzz[j] * fl1_fz);

                t_yyyy_yyy[j] = fl_s_0_0 * (7.5 * pa_y[j] * fl3_fx + 5.625 * pb_y[j] * fl3_fx + 3.0 * pa_yyy[j] * fl2_fx + 13.5 * pa2pb_yy_y[j] * fl2_fx + 9.0 * pa2pb_y_yy[j] * fl2_fx + 1.5 * pa2pb_yyyy_y[j] * fl1_fx + 6.0 * pa2pb_yyy_yy[j] * fl1_fx + 0.75 * pb_yyy[j] * fl2_fx + 3.0 * pa2pb_yy_yyy[j] * fl1_fx + pa2pb_yyyy_yyy[j]);

                t_yyyy_yyy[j] += fl_r_0_0 * (60.0 * pa_y[j] * fl3_fx * fl1_fz - 9.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pb_y[j] * fl3_fx * fl1_fz + 30.0 * pa_yyy[j] * fl1_fz * fl2_fx + 135.0 * pa2pb_yy_y[j] * fl2_fx * fl1_fz - 2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyyy_y[j] * fl1_fz * fl1_fgb + 90.0 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyyy_y[j] * fl1_fz * fl1_fx + 72.0 * pa2pb_yyy_yy[j] * fl1_fz * fl1_fx - 3.0 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_yyy[j] * fl1_fz * fl1_fga + 7.5 * pb_yyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_yyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_yyy[j] * fl1_fz);

                t_yyyy_yyz[j] = fl_s_0_0 * (1.875 * pb_z[j] * fl3_fx + 4.5 * pa2pb_yy_z[j] * fl2_fx + 6.0 * pa2pb_y_yz[j] * fl2_fx + 0.5 * pa2pb_yyyy_z[j] * fl1_fx + 4.0 * pa2pb_yyy_yz[j] * fl1_fx + 0.75 * pb_yyz[j] * fl2_fx + 3.0 * pa2pb_yy_yyz[j] * fl1_fx + pa2pb_yyyy_yyz[j]);

                t_yyyy_yyz[j] += fl_r_0_0 * (-4.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga + 15.0 * pb_z[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_yy_z[j] * fl2_fx * fl1_fz - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyy_z[j] * fl1_fz * fl1_fgb + 60.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyy_z[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_yyy_yz[j] * fl1_fz * fl1_fx - 3.0 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_yyz[j] * fl1_fz * fl1_fga + 7.5 * pb_yyz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_yyz[j] * fl1_fz);

                t_yyyy_yzz[j] = fl_s_0_0 * (1.5 * pa_y[j] * fl3_fx + pa_yyy[j] * fl2_fx + 0.375 * pb_y[j] * fl3_fx + 1.5 * pa2pb_yy_y[j] * fl2_fx + 3.0 * pa2pb_y_zz[j] * fl2_fx + 0.5 * pa2pb_yyyy_y[j] * fl1_fx + 2.0 * pa2pb_yyy_zz[j] * fl1_fx + 0.75 * pb_yzz[j] * fl2_fx + 3.0 * pa2pb_yy_yzz[j] * fl1_fx + pa2pb_yyyy_yzz[j]);

                t_yyyy_yzz[j] += fl_r_0_0 * (-3.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_y[j] * fl3_fx * fl1_fz + 10.0 * pa_yyy[j] * fl1_fz * fl2_fx - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fga + 3.0 * pb_y[j] * fl3_fx * fl1_fz - pa2pb_yyyy_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyy_y[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yyy_zz[j] * fl1_fz * fl1_fx - 3.0 * pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_yzz[j] * fl1_fz * fl1_fga + 7.5 * pb_yzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_yzz[j] * fl1_fz);

                t_yyyy_zzz[j] = fl_s_0_0 * (1.125 * pb_z[j] * fl3_fx + 4.5 * pa2pb_yy_z[j] * fl2_fx + 1.5 * pa2pb_yyyy_z[j] * fl1_fx + 0.75 * pb_zzz[j] * fl2_fx + 3.0 * pa2pb_yy_zzz[j] * fl1_fx + pa2pb_yyyy_zzz[j]);

                t_yyyy_zzz[j] += fl_r_0_0 * (-2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx + 9.0 * pb_z[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_yyyy_z[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_yy_z[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yyyy_z[j] * fl1_fz * fl1_fx - 3.0 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_yy_zzz[j] * fl1_fz * fl1_fga + 7.5 * pb_zzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_yy_zzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yyyy_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_110_115(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_yz_xxx = pa2pbDistances.data(646 * idx + 142);

            auto pa2pb_yz_xxy = pa2pbDistances.data(646 * idx + 143);

            auto pa2pb_yz_xxz = pa2pbDistances.data(646 * idx + 144);

            auto pa2pb_yz_xyy = pa2pbDistances.data(646 * idx + 145);

            auto pa2pb_yz_xyz = pa2pbDistances.data(646 * idx + 146);

            auto pa2pb_yyy_xx = pa2pbDistances.data(646 * idx + 288);

            auto pa2pb_yyy_xy = pa2pbDistances.data(646 * idx + 289);

            auto pa2pb_yyz_xx = pa2pbDistances.data(646 * idx + 307);

            auto pa2pb_yyz_xy = pa2pbDistances.data(646 * idx + 308);

            auto pa2pb_yyz_xz = pa2pbDistances.data(646 * idx + 309);

            auto pa2pb_yyyz_x = pa2pbDistances.data(646 * idx + 570);

            auto pa2pb_yyyz_y = pa2pbDistances.data(646 * idx + 571);

            auto pa2pb_yyyz_z = pa2pbDistances.data(646 * idx + 572);

            auto pa2pb_yyyz_xxx = pa2pbDistances.data(646 * idx + 579);

            auto pa2pb_yyyz_xxy = pa2pbDistances.data(646 * idx + 580);

            auto pa2pb_yyyz_xxz = pa2pbDistances.data(646 * idx + 581);

            auto pa2pb_yyyz_xyy = pa2pbDistances.data(646 * idx + 582);

            auto pa2pb_yyyz_xyz = pa2pbDistances.data(646 * idx + 583);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyz_xxx = primBuffer.data(150 * idx + 110);

            auto t_yyyz_xxy = primBuffer.data(150 * idx + 111);

            auto t_yyyz_xxz = primBuffer.data(150 * idx + 112);

            auto t_yyyz_xyy = primBuffer.data(150 * idx + 113);

            auto t_yyyz_xyz = primBuffer.data(150 * idx + 114);

            // Batch of Integrals (110,115)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_y_xy, pa2pb_yy_x, pa2pb_yyy_xx, \
                                     pa2pb_yyy_xy, pa2pb_yyyz_x, pa2pb_yyyz_xxx, pa2pb_yyyz_xxy, pa2pb_yyyz_xxz, \
                                     pa2pb_yyyz_xyy, pa2pb_yyyz_xyz, pa2pb_yyyz_y, pa2pb_yyyz_z, pa2pb_yyz_xx, \
                                     pa2pb_yyz_xy, pa2pb_yyz_xz, pa2pb_yz_x, pa2pb_yz_xxx, pa2pb_yz_xxy, pa2pb_yz_xxz, \
                                     pa2pb_yz_xyy, pa2pb_yz_xyz, pa2pb_yz_y, pa2pb_yz_z, pa2pb_z_xx, pa2pb_z_xy, \
                                     pa2pb_z_xz, pa_y, pa_yyy, pa_yyz, pa_z, pb_x, r_0_0, s_0_0, t_yyyz_xxx, t_yyyz_xxy, \
                                     t_yyyz_xxz, t_yyyz_xyy, t_yyyz_xyz: VLX_ALIGN)
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

                t_yyyz_xxx[j] = fl_s_0_0 * (2.25 * pa2pb_yz_x[j] * fl2_fx + 1.5 * pa2pb_yyyz_x[j] * fl1_fx + 1.5 * pa2pb_yz_xxx[j] * fl1_fx + pa2pb_yyyz_xxx[j]);

                t_yyyz_xxx[j] += fl_r_0_0 * (-4.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyyz_x[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_yz_x[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyyz_x[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_xxx[j] * fl1_fz);

                t_yyyz_xxy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pa_yyz[j] * fl2_fx + 0.75 * pa2pb_yz_y[j] * fl2_fx + 0.75 * pa2pb_z_xx[j] * fl2_fx + 0.5 * pa2pb_yyyz_y[j] * fl1_fx + 1.5 * pa2pb_yyz_xx[j] * fl1_fx + 1.5 * pa2pb_yz_xxy[j] * fl1_fx + pa2pb_yyyz_xxy[j]);

                t_yyyz_xxy[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 7.5 * pa_yyz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyz_y[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yyz_xx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_xxy[j] * fl1_fz);

                t_yyyz_xxz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.25 * pa_yyy[j] * fl2_fx + 0.75 * pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_xx[j] * fl2_fx + 0.5 * pa2pb_yyyz_z[j] * fl1_fx + 0.5 * pa2pb_yyy_xx[j] * fl1_fx + 1.5 * pa2pb_yz_xxz[j] * fl1_fx + pa2pb_yyyz_xxz[j]);

                t_yyyz_xxz[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 0.5 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 2.5 * pa_yyy[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyyz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyz_z[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyy_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_xxz[j] * fl1_fz);

                t_yyyz_xyy[j] = fl_s_0_0 * (2.25 * pa2pb_yz_x[j] * fl2_fx + 1.5 * pa2pb_z_xy[j] * fl2_fx + 0.5 * pa2pb_yyyz_x[j] * fl1_fx + 3.0 * pa2pb_yyz_xy[j] * fl1_fx + 1.5 * pa2pb_yz_xyy[j] * fl1_fx + pa2pb_yyyz_xyy[j]);

                t_yyyz_xyy[j] += fl_r_0_0 * (22.5 * pa2pb_yz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyz_x[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yyz_xy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_xyy[j] * fl1_fz);

                t_yyyz_xyz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_yy_x[j] * fl2_fx + 0.75 * pa2pb_y_xy[j] * fl2_fx + 0.75 * pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_yyy_xy[j] * fl1_fx + 1.5 * pa2pb_yyz_xz[j] * fl1_fx + 1.5 * pa2pb_yz_xyz[j] * fl1_fx + pa2pb_yyyz_xyz[j]);

                t_yyyz_xyz[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xy[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fga + 7.5 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyy_xy[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yyz_xz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_115_120(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_yz_xzz = pa2pbDistances.data(646 * idx + 147);

            auto pa2pb_yz_yyy = pa2pbDistances.data(646 * idx + 148);

            auto pa2pb_yz_yyz = pa2pbDistances.data(646 * idx + 149);

            auto pa2pb_yz_yzz = pa2pbDistances.data(646 * idx + 150);

            auto pa2pb_yz_zzz = pa2pbDistances.data(646 * idx + 151);

            auto pa2pb_yyy_xz = pa2pbDistances.data(646 * idx + 290);

            auto pa2pb_yyy_yy = pa2pbDistances.data(646 * idx + 291);

            auto pa2pb_yyy_yz = pa2pbDistances.data(646 * idx + 292);

            auto pa2pb_yyy_zz = pa2pbDistances.data(646 * idx + 293);

            auto pa2pb_yyz_yy = pa2pbDistances.data(646 * idx + 310);

            auto pa2pb_yyz_yz = pa2pbDistances.data(646 * idx + 311);

            auto pa2pb_yyz_zz = pa2pbDistances.data(646 * idx + 312);

            auto pa2pb_yyyz_x = pa2pbDistances.data(646 * idx + 570);

            auto pa2pb_yyyz_y = pa2pbDistances.data(646 * idx + 571);

            auto pa2pb_yyyz_z = pa2pbDistances.data(646 * idx + 572);

            auto pa2pb_yyyz_xzz = pa2pbDistances.data(646 * idx + 584);

            auto pa2pb_yyyz_yyy = pa2pbDistances.data(646 * idx + 585);

            auto pa2pb_yyyz_yyz = pa2pbDistances.data(646 * idx + 586);

            auto pa2pb_yyyz_yzz = pa2pbDistances.data(646 * idx + 587);

            auto pa2pb_yyyz_zzz = pa2pbDistances.data(646 * idx + 588);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyyz_xzz = primBuffer.data(150 * idx + 115);

            auto t_yyyz_yyy = primBuffer.data(150 * idx + 116);

            auto t_yyyz_yyz = primBuffer.data(150 * idx + 117);

            auto t_yyyz_yzz = primBuffer.data(150 * idx + 118);

            auto t_yyyz_zzz = primBuffer.data(150 * idx + 119);

            // Batch of Integrals (115,120)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yy_y, \
                                     pa2pb_yy_z, pa2pb_yyy_xz, pa2pb_yyy_yy, pa2pb_yyy_yz, pa2pb_yyy_zz, pa2pb_yyyz_x, \
                                     pa2pb_yyyz_xzz, pa2pb_yyyz_y, pa2pb_yyyz_yyy, pa2pb_yyyz_yyz, pa2pb_yyyz_yzz, \
                                     pa2pb_yyyz_z, pa2pb_yyyz_zzz, pa2pb_yyz_yy, pa2pb_yyz_yz, pa2pb_yyz_zz, pa2pb_yz_x, \
                                     pa2pb_yz_xzz, pa2pb_yz_y, pa2pb_yz_yyy, pa2pb_yz_yyz, pa2pb_yz_yzz, pa2pb_yz_z, \
                                     pa2pb_yz_zzz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa_y, pa_yyy, pa_yyz, pa_z, pb_y, pb_z, \
                                     r_0_0, s_0_0, t_yyyz_xzz, t_yyyz_yyy, t_yyyz_yyz, t_yyyz_yzz, t_yyyz_zzz: VLX_ALIGN)
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

                t_yyyz_xzz[j] = fl_s_0_0 * (0.75 * pa2pb_yz_x[j] * fl2_fx + 1.5 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pa2pb_yyyz_x[j] * fl1_fx + pa2pb_yyy_xz[j] * fl1_fx + 1.5 * pa2pb_yz_xzz[j] * fl1_fx + pa2pb_yyyz_xzz[j]);

                t_yyyz_xzz[j] += fl_r_0_0 * (-1.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyyz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_x[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyz_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_xz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_xzz[j] * fl1_fz);

                t_yyyz_yyy[j] = fl_s_0_0 * (1.875 * pa_z[j] * fl3_fx + 2.25 * pa_yyz[j] * fl2_fx + 6.75 * pa2pb_yz_y[j] * fl2_fx + 2.25 * pa2pb_z_yy[j] * fl2_fx + 1.5 * pa2pb_yyyz_y[j] * fl1_fx + 4.5 * pa2pb_yyz_yy[j] * fl1_fx + 1.5 * pa2pb_yz_yyy[j] * fl1_fx + pa2pb_yyyz_yyy[j]);

                t_yyyz_yyy[j] += fl_r_0_0 * (15.0 * pa_z[j] * fl3_fx * fl1_fz - 2.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_yyz[j] * fl2_fx * fl1_fz + 67.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yyyz_y[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyyz_y[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_yyz_yy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_yyy[j] * fl1_fz);

                t_yyyz_yyz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pb_y[j] * fl3_fx + 0.25 * pa_yyy[j] * fl2_fx + 1.5 * pa2pb_yy_y[j] * fl2_fx + 2.25 * pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_yy[j] * fl2_fx + 1.5 * pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_yyyz_z[j] * fl1_fx + 0.5 * pa2pb_yyy_yy[j] * fl1_fx + 3.0 * pa2pb_yyz_yz[j] * fl1_fx + 1.5 * pa2pb_yz_yyz[j] * fl1_fx + pa2pb_yyyz_yyz[j]);

                t_yyyz_yyz[j] += fl_r_0_0 * (9.0 * pa_y[j] * fl3_fx * fl1_fz - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_y[j] * fl3_fx * fl1_fz + 2.5 * pa_yyy[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_yy_y[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yy[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyz_z[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_yyy_yy[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yyz_yz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_yyz[j] * fl1_fz);

                t_yyyz_yzz[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.75 * pb_z[j] * fl3_fx + 0.75 * pa_yyz[j] * fl2_fx + 1.5 * pa2pb_yy_z[j] * fl2_fx + 0.75 * pa2pb_yz_y[j] * fl2_fx + 1.5 * pa2pb_y_yz[j] * fl2_fx + 0.75 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pa2pb_yyyz_y[j] * fl1_fx + pa2pb_yyy_yz[j] * fl1_fx + 1.5 * pa2pb_yyz_zz[j] * fl1_fx + 1.5 * pa2pb_yz_yzz[j] * fl1_fx + pa2pb_yyyz_yzz[j]);

                t_yyyz_yzz[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 6.0 * pb_z[j] * fl3_fx * fl1_fz + 7.5 * pa_yyz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_yy_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_yz[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yyyz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyyz_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyy_yz[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yyz_zz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_yzz[j] * fl1_fz);

                t_yyyz_zzz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx + 0.75 * pa_yyy[j] * fl2_fx + 2.25 * pa2pb_yz_z[j] * fl2_fx + 2.25 * pa2pb_y_zz[j] * fl2_fx + 1.5 * pa2pb_yyyz_z[j] * fl1_fx + 1.5 * pa2pb_yyy_zz[j] * fl1_fx + 1.5 * pa2pb_yz_zzz[j] * fl1_fx + pa2pb_yyyz_zzz[j]);

                t_yyyz_zzz[j] += fl_r_0_0 * (-2.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_y[j] * fl1_fz * fl1_fga * fl2_fx - 1.5 * pa_yyy[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_y[j] * fl3_fx * fl1_fz + 7.5 * pa_yyy[j] * fl1_fz * fl2_fx - 4.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyyz_z[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyyz_z[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yyy_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_zzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyyz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_120_125(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_yy_xxx = pa2pbDistances.data(646 * idx + 123);

            auto pa2pb_yy_xxy = pa2pbDistances.data(646 * idx + 124);

            auto pa2pb_yy_xxz = pa2pbDistances.data(646 * idx + 125);

            auto pa2pb_yy_xyy = pa2pbDistances.data(646 * idx + 126);

            auto pa2pb_yy_xyz = pa2pbDistances.data(646 * idx + 127);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_zz_xxx = pa2pbDistances.data(646 * idx + 161);

            auto pa2pb_zz_xxy = pa2pbDistances.data(646 * idx + 162);

            auto pa2pb_zz_xxz = pa2pbDistances.data(646 * idx + 163);

            auto pa2pb_zz_xyy = pa2pbDistances.data(646 * idx + 164);

            auto pa2pb_zz_xyz = pa2pbDistances.data(646 * idx + 165);

            auto pa2pb_yyz_xx = pa2pbDistances.data(646 * idx + 307);

            auto pa2pb_yyz_xy = pa2pbDistances.data(646 * idx + 308);

            auto pa2pb_yzz_xx = pa2pbDistances.data(646 * idx + 326);

            auto pa2pb_yzz_xy = pa2pbDistances.data(646 * idx + 327);

            auto pa2pb_yzz_xz = pa2pbDistances.data(646 * idx + 328);

            auto pa2pb_yyzz_x = pa2pbDistances.data(646 * idx + 589);

            auto pa2pb_yyzz_y = pa2pbDistances.data(646 * idx + 590);

            auto pa2pb_yyzz_z = pa2pbDistances.data(646 * idx + 591);

            auto pa2pb_yyzz_xxx = pa2pbDistances.data(646 * idx + 598);

            auto pa2pb_yyzz_xxy = pa2pbDistances.data(646 * idx + 599);

            auto pa2pb_yyzz_xxz = pa2pbDistances.data(646 * idx + 600);

            auto pa2pb_yyzz_xyy = pa2pbDistances.data(646 * idx + 601);

            auto pa2pb_yyzz_xyz = pa2pbDistances.data(646 * idx + 602);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyzz_xxx = primBuffer.data(150 * idx + 120);

            auto t_yyzz_xxy = primBuffer.data(150 * idx + 121);

            auto t_yyzz_xxz = primBuffer.data(150 * idx + 122);

            auto t_yyzz_xyy = primBuffer.data(150 * idx + 123);

            auto t_yyzz_xyz = primBuffer.data(150 * idx + 124);

            // Batch of Integrals (120,125)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_y_xy, pa2pb_y_xz, pa2pb_yy_x, \
                                     pa2pb_yy_xxx, pa2pb_yy_xxy, pa2pb_yy_xxz, pa2pb_yy_xyy, pa2pb_yy_xyz, pa2pb_yy_y, \
                                     pa2pb_yy_z, pa2pb_yyz_xx, pa2pb_yyz_xy, pa2pb_yyzz_x, pa2pb_yyzz_xxx, \
                                     pa2pb_yyzz_xxy, pa2pb_yyzz_xxz, pa2pb_yyzz_xyy, pa2pb_yyzz_xyz, pa2pb_yyzz_y, \
                                     pa2pb_yyzz_z, pa2pb_yz_x, pa2pb_yzz_xx, pa2pb_yzz_xy, pa2pb_yzz_xz, pa2pb_z_xx, \
                                     pa2pb_z_xy, pa2pb_zz_x, pa2pb_zz_xxx, pa2pb_zz_xxy, pa2pb_zz_xxz, pa2pb_zz_xyy, \
                                     pa2pb_zz_xyz, pa2pb_zz_y, pa2pb_zz_z, pa_y, pa_yyz, pa_yzz, pa_z, pb_x, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xyy, pb_xyz, pb_y, pb_z, r_0_0, s_0_0, t_yyzz_xxx, t_yyzz_xxy, t_yyzz_xxz, \
                                     t_yyzz_xyy, t_yyzz_xyz: VLX_ALIGN)
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

                t_yyzz_xxx[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_yy_x[j] * fl2_fx + 0.75 * pa2pb_zz_x[j] * fl2_fx + 1.5 * pa2pb_yyzz_x[j] * fl1_fx + 0.25 * pb_xxx[j] * fl2_fx + 0.5 * pa2pb_yy_xxx[j] * fl1_fx + 0.5 * pa2pb_zz_xxx[j] * fl1_fx + pa2pb_yyzz_xxx[j]);

                t_yyzz_xxx[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_x[j] * fl3_fx * fl1_fz - 1.5 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yy_x[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyzz_x[j] * fl1_fz * fl1_fx - pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxx[j] * fl1_fz * fl1_fga + 2.5 * pb_xxx[j] * fl2_fx * fl1_fz - pa2pb_zz_xxx[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_xxx[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xxx[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_xxx[j] * fl1_fz);

                t_yyzz_xxy[j] = fl_s_0_0 * (0.25 * pa_y[j] * fl3_fx + 0.5 * pa_yzz[j] * fl2_fx + 0.125 * pb_y[j] * fl3_fx + 0.25 * pa2pb_yy_y[j] * fl2_fx + 0.5 * pa2pb_y_xx[j] * fl2_fx + 0.25 * pa2pb_zz_y[j] * fl2_fx + 0.5 * pa2pb_yyzz_y[j] * fl1_fx + pa2pb_yzz_xx[j] * fl1_fx + 0.25 * pb_xxy[j] * fl2_fx + 0.5 * pa2pb_yy_xxy[j] * fl1_fx + 0.5 * pa2pb_zz_xxy[j] * fl1_fx + pa2pb_yyzz_xxy[j]);

                t_yyzz_xxy[j] += fl_r_0_0 * (-0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_y[j] * fl3_fx * fl1_fz + 5.0 * pa_yzz[j] * fl2_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb + pb_y[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_y[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 5.0 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyzz_y[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xx[j] * fl1_fx * fl1_fz - pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxy[j] * fl1_fz * fl1_fga + 2.5 * pb_xxy[j] * fl2_fx * fl1_fz - pa2pb_zz_xxy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_xxy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xxy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_xxy[j] * fl1_fz);

                t_yyzz_xxz[j] = fl_s_0_0 * (0.25 * pa_z[j] * fl3_fx + 0.5 * pa_yyz[j] * fl2_fx + 0.125 * pb_z[j] * fl3_fx + 0.25 * pa2pb_yy_z[j] * fl2_fx + 0.25 * pa2pb_zz_z[j] * fl2_fx + 0.5 * pa2pb_z_xx[j] * fl2_fx + 0.5 * pa2pb_yyzz_z[j] * fl1_fx + pa2pb_yyz_xx[j] * fl1_fx + 0.25 * pb_xxz[j] * fl2_fx + 0.5 * pa2pb_yy_xxz[j] * fl1_fx + 0.5 * pa2pb_zz_xxz[j] * fl1_fx + pa2pb_yyzz_xxz[j]);

                t_yyzz_xxz[j] += fl_r_0_0 * (-0.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 2.0 * pa_z[j] * fl3_fx * fl1_fz + 5.0 * pa_yyz[j] * fl1_fz * fl2_fx - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb + pb_z[j] * fl3_fx * fl1_fz - 0.5 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_xx[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yy_z[j] * fl1_fz * fl2_fx + 2.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyzz_z[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyz_xx[j] * fl1_fz * fl1_fx - pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xxz[j] * fl1_fz * fl1_fga + 2.5 * pb_xxz[j] * fl2_fx * fl1_fz - pa2pb_zz_xxz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_xxz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xxz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_xxz[j] * fl1_fz);

                t_yyzz_xyy[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_zz_x[j] * fl2_fx + 0.25 * pa2pb_yy_x[j] * fl2_fx + pa2pb_y_xy[j] * fl2_fx + 0.5 * pa2pb_yyzz_x[j] * fl1_fx + 2.0 * pa2pb_yzz_xy[j] * fl1_fx + 0.25 * pb_xyy[j] * fl2_fx + 0.5 * pa2pb_yy_xyy[j] * fl1_fx + 0.5 * pa2pb_zz_xyy[j] * fl1_fx + pa2pb_yyzz_xyy[j]);

                t_yyzz_xyy[j] += fl_r_0_0 * (-pb_x[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yy_x[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyzz_x[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yzz_xy[j] * fl1_fx * fl1_fz - pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xyy[j] * fl1_fz * fl1_fga + 2.5 * pb_xyy[j] * fl2_fx * fl1_fz - pa2pb_zz_xyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_xyy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_xyy[j] * fl1_fz);

                t_yyzz_xyz[j] = fl_s_0_0 * (pa2pb_yz_x[j] * fl2_fx + 0.5 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pa2pb_z_xy[j] * fl2_fx + pa2pb_yyz_xy[j] * fl1_fx + pa2pb_yzz_xz[j] * fl1_fx + 0.25 * pb_xyz[j] * fl2_fx + 0.5 * pa2pb_yy_xyz[j] * fl1_fx + 0.5 * pa2pb_zz_xyz[j] * fl1_fx + pa2pb_yyzz_xyz[j]);

                t_yyzz_xyz[j] += fl_r_0_0 * (10.0 * pa2pb_yz_x[j] * fl2_fx * fl1_fz - pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_z_xy[j] * fl1_fz * fl1_fga * fl1_fx + 5.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 12.0 * pa2pb_yyz_xy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_xz[j] * fl1_fx * fl1_fz - pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xyz[j] * fl1_fz * fl1_fga + 2.5 * pb_xyz[j] * fl2_fx * fl1_fz - pa2pb_zz_xyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_xyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_125_130(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_yy_x = pa2pbDistances.data(646 * idx + 114);

            auto pa2pb_yy_y = pa2pbDistances.data(646 * idx + 115);

            auto pa2pb_yy_z = pa2pbDistances.data(646 * idx + 116);

            auto pa2pb_yy_xzz = pa2pbDistances.data(646 * idx + 128);

            auto pa2pb_yy_yyy = pa2pbDistances.data(646 * idx + 129);

            auto pa2pb_yy_yyz = pa2pbDistances.data(646 * idx + 130);

            auto pa2pb_yy_yzz = pa2pbDistances.data(646 * idx + 131);

            auto pa2pb_yy_zzz = pa2pbDistances.data(646 * idx + 132);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_zz_xzz = pa2pbDistances.data(646 * idx + 166);

            auto pa2pb_zz_yyy = pa2pbDistances.data(646 * idx + 167);

            auto pa2pb_zz_yyz = pa2pbDistances.data(646 * idx + 168);

            auto pa2pb_zz_yzz = pa2pbDistances.data(646 * idx + 169);

            auto pa2pb_zz_zzz = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_yyz_xz = pa2pbDistances.data(646 * idx + 309);

            auto pa2pb_yyz_yy = pa2pbDistances.data(646 * idx + 310);

            auto pa2pb_yyz_yz = pa2pbDistances.data(646 * idx + 311);

            auto pa2pb_yyz_zz = pa2pbDistances.data(646 * idx + 312);

            auto pa2pb_yzz_yy = pa2pbDistances.data(646 * idx + 329);

            auto pa2pb_yzz_yz = pa2pbDistances.data(646 * idx + 330);

            auto pa2pb_yzz_zz = pa2pbDistances.data(646 * idx + 331);

            auto pa2pb_yyzz_x = pa2pbDistances.data(646 * idx + 589);

            auto pa2pb_yyzz_y = pa2pbDistances.data(646 * idx + 590);

            auto pa2pb_yyzz_z = pa2pbDistances.data(646 * idx + 591);

            auto pa2pb_yyzz_xzz = pa2pbDistances.data(646 * idx + 603);

            auto pa2pb_yyzz_yyy = pa2pbDistances.data(646 * idx + 604);

            auto pa2pb_yyzz_yyz = pa2pbDistances.data(646 * idx + 605);

            auto pa2pb_yyzz_yzz = pa2pbDistances.data(646 * idx + 606);

            auto pa2pb_yyzz_zzz = pa2pbDistances.data(646 * idx + 607);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yyzz_xzz = primBuffer.data(150 * idx + 125);

            auto t_yyzz_yyy = primBuffer.data(150 * idx + 126);

            auto t_yyzz_yyz = primBuffer.data(150 * idx + 127);

            auto t_yyzz_yzz = primBuffer.data(150 * idx + 128);

            auto t_yyzz_zzz = primBuffer.data(150 * idx + 129);

            // Batch of Integrals (125,130)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yy_x, \
                                     pa2pb_yy_xzz, pa2pb_yy_y, pa2pb_yy_yyy, pa2pb_yy_yyz, pa2pb_yy_yzz, pa2pb_yy_z, \
                                     pa2pb_yy_zzz, pa2pb_yyz_xz, pa2pb_yyz_yy, pa2pb_yyz_yz, pa2pb_yyz_zz, pa2pb_yyzz_x, \
                                     pa2pb_yyzz_xzz, pa2pb_yyzz_y, pa2pb_yyzz_yyy, pa2pb_yyzz_yyz, pa2pb_yyzz_yzz, \
                                     pa2pb_yyzz_z, pa2pb_yyzz_zzz, pa2pb_yz_y, pa2pb_yz_z, pa2pb_yzz_yy, pa2pb_yzz_yz, \
                                     pa2pb_yzz_zz, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa2pb_zz_x, \
                                     pa2pb_zz_xzz, pa2pb_zz_y, pa2pb_zz_yyy, pa2pb_zz_yyz, pa2pb_zz_yzz, pa2pb_zz_z, \
                                     pa2pb_zz_zzz, pa_y, pa_yyz, pa_yzz, pa_z, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, \
                                     pb_zzz, r_0_0, s_0_0, t_yyzz_xzz, t_yyzz_yyy, t_yyzz_yyz, t_yyzz_yzz, t_yyzz_zzz: VLX_ALIGN)
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

                t_yyzz_xzz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_yy_x[j] * fl2_fx + 0.25 * pa2pb_zz_x[j] * fl2_fx + pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_yyzz_x[j] * fl1_fx + 2.0 * pa2pb_yyz_xz[j] * fl1_fx + 0.25 * pb_xzz[j] * fl2_fx + 0.5 * pa2pb_yy_xzz[j] * fl1_fx + 0.5 * pa2pb_zz_xzz[j] * fl1_fx + pa2pb_yyzz_xzz[j]);

                t_yyzz_xzz[j] += fl_r_0_0 * (-pb_x[j] * fl1_fz * fl1_fga * fl2_fx + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_x[j] * fl2_fx * fl1_fz - 0.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_x[j] * fl1_fz * fl1_fga * fl1_fx - 0.5 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_xz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_x[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyzz_x[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yyz_xz[j] * fl1_fz * fl1_fx - pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_xzz[j] * fl1_fz * fl1_fga + 2.5 * pb_xzz[j] * fl2_fx * fl1_fz - pa2pb_zz_xzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_xzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_xzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_xzz[j] * fl1_fz);

                t_yyzz_yyy[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 1.125 * pb_y[j] * fl3_fx + 1.5 * pa_yzz[j] * fl2_fx + 2.25 * pa2pb_zz_y[j] * fl2_fx + 0.75 * pa2pb_yy_y[j] * fl2_fx + 1.5 * pa2pb_y_yy[j] * fl2_fx + 1.5 * pa2pb_yyzz_y[j] * fl1_fx + 3.0 * pa2pb_yzz_yy[j] * fl1_fx + 0.25 * pb_yyy[j] * fl2_fx + 0.5 * pa2pb_yy_yyy[j] * fl1_fx + 0.5 * pa2pb_zz_yyy[j] * fl1_fx + pa2pb_yyzz_yyy[j]);

                t_yyzz_yyy[j] += fl_r_0_0 * (-1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_y[j] * fl3_fx * fl1_fz + 9.0 * pb_y[j] * fl3_fx * fl1_fz + 15.0 * pa_yzz[j] * fl2_fx * fl1_fz + 22.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yy_y[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyzz_y[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yzz_yy[j] * fl1_fx * fl1_fz - pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_yyy[j] * fl1_fz * fl1_fga + 2.5 * pb_yyy[j] * fl2_fx * fl1_fz - pa2pb_zz_yyy[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_yyy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yyy[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_yyy[j] * fl1_fz);

                t_yyzz_yyz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl3_fx + 0.375 * pb_z[j] * fl3_fx + 0.5 * pa_yyz[j] * fl2_fx + 2.0 * pa2pb_yz_y[j] * fl2_fx + 0.75 * pa2pb_zz_z[j] * fl2_fx + 0.25 * pa2pb_yy_z[j] * fl2_fx + pa2pb_y_yz[j] * fl2_fx + 0.5 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pa2pb_yyzz_z[j] * fl1_fx + pa2pb_yyz_yy[j] * fl1_fx + 2.0 * pa2pb_yzz_yz[j] * fl1_fx + 0.25 * pb_yyz[j] * fl2_fx + 0.5 * pa2pb_yy_yyz[j] * fl1_fx + 0.5 * pa2pb_zz_yyz[j] * fl1_fx + pa2pb_yyzz_yyz[j]);

                t_yyzz_yyz[j] += fl_r_0_0 * (6.0 * pa_z[j] * fl3_fx * fl1_fz - 0.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_z[j] * fl3_fx * fl1_fz + 5.0 * pa_yyz[j] * fl1_fz * fl2_fx + 20.0 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz - 0.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_z_yy[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_z[j] * fl1_fz * fl1_fgb + 2.5 * pa2pb_yy_z[j] * fl1_fz * fl2_fx + 10.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 5.0 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyzz_z[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yyz_yy[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yzz_yz[j] * fl1_fx * fl1_fz - pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_yyz[j] * fl1_fz * fl1_fga + 2.5 * pb_yyz[j] * fl2_fx * fl1_fz - pa2pb_zz_yyz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_yyz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yyz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_yyz[j] * fl1_fz);

                t_yyzz_yzz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx + 0.375 * pb_y[j] * fl3_fx + 0.75 * pa2pb_yy_y[j] * fl2_fx + 0.5 * pa_yzz[j] * fl2_fx + 2.0 * pa2pb_yz_z[j] * fl2_fx + 0.5 * pa2pb_y_zz[j] * fl2_fx + 0.25 * pa2pb_zz_y[j] * fl2_fx + pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_yyzz_y[j] * fl1_fx + 2.0 * pa2pb_yyz_yz[j] * fl1_fx + pa2pb_yzz_zz[j] * fl1_fx + 0.25 * pb_yzz[j] * fl2_fx + 0.5 * pa2pb_yy_yzz[j] * fl1_fx + 0.5 * pa2pb_zz_yzz[j] * fl1_fx + pa2pb_yyzz_yzz[j]);

                t_yyzz_yzz[j] += fl_r_0_0 * (6.0 * pa_y[j] * fl3_fx * fl1_fz - 0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - pb_y[j] * fl1_fz * fl1_fga * fl2_fx - pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_yy_y[j] * fl2_fx * fl1_fz + 5.0 * pa_yzz[j] * fl2_fx * fl1_fz + 20.0 * pa2pb_yz_z[j] * fl2_fx * fl1_fz - 0.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_yy_y[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fga - 0.5 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 0.5 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx - 2.0 * pa2pb_z_yz[j] * fl1_fz * fl1_fga * fl1_fx - pa2pb_yyzz_y[j] * fl1_fz * fl1_fgb + 5.0 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 2.5 * pa2pb_zz_y[j] * fl2_fx * fl1_fz + 10.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yyzz_y[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_yyz_yz[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_yzz_zz[j] * fl1_fx * fl1_fz - pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_yzz[j] * fl1_fz * fl1_fga + 2.5 * pb_yzz[j] * fl2_fx * fl1_fz - pa2pb_zz_yzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_yzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_yzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_yzz[j] * fl1_fz);

                t_yyzz_zzz[j] = fl_s_0_0 * (0.75 * pa_z[j] * fl3_fx + 1.125 * pb_z[j] * fl3_fx + 1.5 * pa_yyz[j] * fl2_fx + 2.25 * pa2pb_yy_z[j] * fl2_fx + 0.75 * pa2pb_zz_z[j] * fl2_fx + 1.5 * pa2pb_z_zz[j] * fl2_fx + 1.5 * pa2pb_yyzz_z[j] * fl1_fx + 3.0 * pa2pb_yyz_zz[j] * fl1_fx + 0.25 * pb_zzz[j] * fl2_fx + 0.5 * pa2pb_yy_zzz[j] * fl1_fx + 0.5 * pa2pb_zz_zzz[j] * fl1_fx + pa2pb_yyzz_zzz[j]);

                t_yyzz_zzz[j] += fl_r_0_0 * (-1.5 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pb_z[j] * fl1_fz * fl1_fga * fl2_fx - 3.0 * pa_yyz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pa_z[j] * fl3_fx * fl1_fz + 9.0 * pb_z[j] * fl3_fx * fl1_fz + 15.0 * pa_yyz[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_yy_z[j] * fl2_fx * fl1_fz - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yy_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_zz[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yyzz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_zz_z[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yyzz_z[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yyz_zz[j] * fl1_fz * fl1_fx - pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yy_zzz[j] * fl1_fz * fl1_fga + 2.5 * pb_zzz[j] * fl2_fx * fl1_fz - pa2pb_zz_zzz[j] * fl1_fz * fl1_fga + 6.0 * pa2pb_yy_zzz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zz_zzz[j] * fl1_fx * fl1_fz + 14.0 * pa2pb_yyzz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_130_135(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xx = pa2pbDistances.data(646 * idx + 22);

            auto pa2pb_y_xy = pa2pbDistances.data(646 * idx + 23);

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_yz_xxx = pa2pbDistances.data(646 * idx + 142);

            auto pa2pb_yz_xxy = pa2pbDistances.data(646 * idx + 143);

            auto pa2pb_yz_xxz = pa2pbDistances.data(646 * idx + 144);

            auto pa2pb_yz_xyy = pa2pbDistances.data(646 * idx + 145);

            auto pa2pb_yz_xyz = pa2pbDistances.data(646 * idx + 146);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_yzz_xx = pa2pbDistances.data(646 * idx + 326);

            auto pa2pb_yzz_xy = pa2pbDistances.data(646 * idx + 327);

            auto pa2pb_zzz_xx = pa2pbDistances.data(646 * idx + 345);

            auto pa2pb_zzz_xy = pa2pbDistances.data(646 * idx + 346);

            auto pa2pb_zzz_xz = pa2pbDistances.data(646 * idx + 347);

            auto pa2pb_yzzz_x = pa2pbDistances.data(646 * idx + 608);

            auto pa2pb_yzzz_y = pa2pbDistances.data(646 * idx + 609);

            auto pa2pb_yzzz_z = pa2pbDistances.data(646 * idx + 610);

            auto pa2pb_yzzz_xxx = pa2pbDistances.data(646 * idx + 617);

            auto pa2pb_yzzz_xxy = pa2pbDistances.data(646 * idx + 618);

            auto pa2pb_yzzz_xxz = pa2pbDistances.data(646 * idx + 619);

            auto pa2pb_yzzz_xyy = pa2pbDistances.data(646 * idx + 620);

            auto pa2pb_yzzz_xyz = pa2pbDistances.data(646 * idx + 621);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzzz_xxx = primBuffer.data(150 * idx + 130);

            auto t_yzzz_xxy = primBuffer.data(150 * idx + 131);

            auto t_yzzz_xxz = primBuffer.data(150 * idx + 132);

            auto t_yzzz_xyy = primBuffer.data(150 * idx + 133);

            auto t_yzzz_xyz = primBuffer.data(150 * idx + 134);

            // Batch of Integrals (130,135)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xx, pa2pb_y_xy, pa2pb_yz_x, pa2pb_yz_xxx, \
                                     pa2pb_yz_xxy, pa2pb_yz_xxz, pa2pb_yz_xyy, pa2pb_yz_xyz, pa2pb_yz_y, pa2pb_yz_z, \
                                     pa2pb_yzz_xx, pa2pb_yzz_xy, pa2pb_yzzz_x, pa2pb_yzzz_xxx, pa2pb_yzzz_xxy, \
                                     pa2pb_yzzz_xxz, pa2pb_yzzz_xyy, pa2pb_yzzz_xyz, pa2pb_yzzz_y, pa2pb_yzzz_z, \
                                     pa2pb_z_xx, pa2pb_z_xy, pa2pb_z_xz, pa2pb_zz_x, pa2pb_zzz_xx, pa2pb_zzz_xy, \
                                     pa2pb_zzz_xz, pa_y, pa_yzz, pa_z, pa_zzz, pb_x, r_0_0, s_0_0, t_yzzz_xxx, t_yzzz_xxy, \
                                     t_yzzz_xxz, t_yzzz_xyy, t_yzzz_xyz: VLX_ALIGN)
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

                t_yzzz_xxx[j] = fl_s_0_0 * (2.25 * pa2pb_yz_x[j] * fl2_fx + 1.5 * pa2pb_yzzz_x[j] * fl1_fx + 1.5 * pa2pb_yz_xxx[j] * fl1_fx + pa2pb_yzzz_xxx[j]);

                t_yzzz_xxx[j] += fl_r_0_0 * (-4.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_yzzz_x[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_yz_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_yzzz_x[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxx[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_xxx[j] * fl1_fz);

                t_yzzz_xxy[j] = fl_s_0_0 * (0.375 * pa_z[j] * fl3_fx + 0.25 * pa_zzz[j] * fl2_fx + 0.75 * pa2pb_yz_y[j] * fl2_fx + 0.75 * pa2pb_z_xx[j] * fl2_fx + 0.5 * pa2pb_yzzz_y[j] * fl1_fx + 0.5 * pa2pb_zzz_xx[j] * fl1_fx + 1.5 * pa2pb_yz_xxy[j] * fl1_fx + pa2pb_yzzz_xxy[j]);

                t_yzzz_xxy[j] += fl_r_0_0 * (-0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_z[j] * fl3_fx * fl1_fz + 2.5 * pa_zzz[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_y[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_y[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzzz_y[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_xx[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xxy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_xxy[j] * fl1_fz);

                t_yzzz_xxz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pa_yzz[j] * fl2_fx + 0.75 * pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_xx[j] * fl2_fx + 0.5 * pa2pb_yzzz_z[j] * fl1_fx + 1.5 * pa2pb_yzz_xx[j] * fl1_fx + 1.5 * pa2pb_yz_xxz[j] * fl1_fx + pa2pb_yzzz_xxz[j]);

                t_yzzz_xxz[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 7.5 * pa_yzz[j] * fl1_fz * fl2_fx - 1.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_xx[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_y_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzzz_z[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yzz_xx[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xxz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_xxz[j] * fl1_fz);

                t_yzzz_xyy[j] = fl_s_0_0 * (0.75 * pa2pb_yz_x[j] * fl2_fx + 1.5 * pa2pb_z_xy[j] * fl2_fx + 0.5 * pa2pb_yzzz_x[j] * fl1_fx + pa2pb_zzz_xy[j] * fl1_fx + 1.5 * pa2pb_yz_xyy[j] * fl1_fx + pa2pb_yzzz_xyy[j]);

                t_yzzz_xyy[j] += fl_r_0_0 * (-1.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_x[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_x[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzzz_x[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_xy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_xyy[j] * fl1_fz);

                t_yzzz_xyz[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 0.75 * pa2pb_zz_x[j] * fl2_fx + 0.75 * pa2pb_y_xy[j] * fl2_fx + 0.75 * pa2pb_z_xz[j] * fl2_fx + 1.5 * pa2pb_yzz_xy[j] * fl1_fx + 0.5 * pa2pb_zzz_xz[j] * fl1_fx + 1.5 * pa2pb_yz_xyz[j] * fl1_fx + pa2pb_yzzz_xyz[j]);

                t_yzzz_xyz[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga + 3.0 * pb_x[j] * fl3_fx * fl1_fz + 7.5 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_y_xy[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fga + 7.5 * pa2pb_y_xy[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yzz_xy[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_xz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_xyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_135_140(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_y_xz = pa2pbDistances.data(646 * idx + 24);

            auto pa2pb_y_yy = pa2pbDistances.data(646 * idx + 25);

            auto pa2pb_y_yz = pa2pbDistances.data(646 * idx + 26);

            auto pa2pb_y_zz = pa2pbDistances.data(646 * idx + 27);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_yz_x = pa2pbDistances.data(646 * idx + 133);

            auto pa2pb_yz_y = pa2pbDistances.data(646 * idx + 134);

            auto pa2pb_yz_z = pa2pbDistances.data(646 * idx + 135);

            auto pa2pb_yz_xzz = pa2pbDistances.data(646 * idx + 147);

            auto pa2pb_yz_yyy = pa2pbDistances.data(646 * idx + 148);

            auto pa2pb_yz_yyz = pa2pbDistances.data(646 * idx + 149);

            auto pa2pb_yz_yzz = pa2pbDistances.data(646 * idx + 150);

            auto pa2pb_yz_zzz = pa2pbDistances.data(646 * idx + 151);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_yzz_xz = pa2pbDistances.data(646 * idx + 328);

            auto pa2pb_yzz_yy = pa2pbDistances.data(646 * idx + 329);

            auto pa2pb_yzz_yz = pa2pbDistances.data(646 * idx + 330);

            auto pa2pb_yzz_zz = pa2pbDistances.data(646 * idx + 331);

            auto pa2pb_zzz_yy = pa2pbDistances.data(646 * idx + 348);

            auto pa2pb_zzz_yz = pa2pbDistances.data(646 * idx + 349);

            auto pa2pb_zzz_zz = pa2pbDistances.data(646 * idx + 350);

            auto pa2pb_yzzz_x = pa2pbDistances.data(646 * idx + 608);

            auto pa2pb_yzzz_y = pa2pbDistances.data(646 * idx + 609);

            auto pa2pb_yzzz_z = pa2pbDistances.data(646 * idx + 610);

            auto pa2pb_yzzz_xzz = pa2pbDistances.data(646 * idx + 622);

            auto pa2pb_yzzz_yyy = pa2pbDistances.data(646 * idx + 623);

            auto pa2pb_yzzz_yyz = pa2pbDistances.data(646 * idx + 624);

            auto pa2pb_yzzz_yzz = pa2pbDistances.data(646 * idx + 625);

            auto pa2pb_yzzz_zzz = pa2pbDistances.data(646 * idx + 626);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_yzzz_xzz = primBuffer.data(150 * idx + 135);

            auto t_yzzz_yyy = primBuffer.data(150 * idx + 136);

            auto t_yzzz_yyz = primBuffer.data(150 * idx + 137);

            auto t_yzzz_yzz = primBuffer.data(150 * idx + 138);

            auto t_yzzz_zzz = primBuffer.data(150 * idx + 139);

            // Batch of Integrals (135,140)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_y_xz, pa2pb_y_yy, pa2pb_y_yz, pa2pb_y_zz, pa2pb_yz_x, \
                                     pa2pb_yz_xzz, pa2pb_yz_y, pa2pb_yz_yyy, pa2pb_yz_yyz, pa2pb_yz_yzz, pa2pb_yz_z, \
                                     pa2pb_yz_zzz, pa2pb_yzz_xz, pa2pb_yzz_yy, pa2pb_yzz_yz, pa2pb_yzz_zz, pa2pb_yzzz_x, \
                                     pa2pb_yzzz_xzz, pa2pb_yzzz_y, pa2pb_yzzz_yyy, pa2pb_yzzz_yyz, pa2pb_yzzz_yzz, \
                                     pa2pb_yzzz_z, pa2pb_yzzz_zzz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa2pb_zz_y, \
                                     pa2pb_zz_z, pa2pb_zzz_yy, pa2pb_zzz_yz, pa2pb_zzz_zz, pa_y, pa_yzz, pa_z, pa_zzz, pb_y, \
                                     pb_z, r_0_0, s_0_0, t_yzzz_xzz, t_yzzz_yyy, t_yzzz_yyz, t_yzzz_yzz, t_yzzz_zzz: VLX_ALIGN)
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

                t_yzzz_xzz[j] = fl_s_0_0 * (2.25 * pa2pb_yz_x[j] * fl2_fx + 1.5 * pa2pb_y_xz[j] * fl2_fx + 0.5 * pa2pb_yzzz_x[j] * fl1_fx + 3.0 * pa2pb_yzz_xz[j] * fl1_fx + 1.5 * pa2pb_yz_xzz[j] * fl1_fx + pa2pb_yzzz_xzz[j]);

                t_yzzz_xzz[j] += fl_r_0_0 * (22.5 * pa2pb_yz_x[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_x[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzzz_x[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yzz_xz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_xzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_xzz[j] * fl1_fz);

                t_yzzz_yyy[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pa_zzz[j] * fl2_fx + 2.25 * pa2pb_yz_y[j] * fl2_fx + 2.25 * pa2pb_z_yy[j] * fl2_fx + 1.5 * pa2pb_yzzz_y[j] * fl1_fx + 1.5 * pa2pb_zzz_yy[j] * fl1_fx + 1.5 * pa2pb_yz_yyy[j] * fl1_fx + pa2pb_yzzz_yyy[j]);

                t_yzzz_yyy[j] += fl_r_0_0 * (-2.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 9.0 * pa_z[j] * fl3_fx * fl1_fz + 7.5 * pa_zzz[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzzz_y[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_yz_y[j] * fl1_fz * fl2_fx + 22.5 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yzzz_y[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_zzz_yy[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyy[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_yyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_yyy[j] * fl1_fz);

                t_yzzz_yyz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx + 0.75 * pb_y[j] * fl3_fx + 0.75 * pa_yzz[j] * fl2_fx + 1.5 * pa2pb_zz_y[j] * fl2_fx + 0.75 * pa2pb_yz_z[j] * fl2_fx + 0.75 * pa2pb_y_yy[j] * fl2_fx + 1.5 * pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_yzzz_z[j] * fl1_fx + 1.5 * pa2pb_yzz_yy[j] * fl1_fx + pa2pb_zzz_yz[j] * fl1_fx + 1.5 * pa2pb_yz_yyz[j] * fl1_fx + pa2pb_yzzz_yyz[j]);

                t_yzzz_yyz[j] += fl_r_0_0 * (-0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 3.0 * pa_y[j] * fl3_fx * fl1_fz + 6.0 * pb_y[j] * fl3_fx * fl1_fz + 7.5 * pa_yzz[j] * fl1_fz * fl2_fx + 15.0 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 1.5 * pa2pb_y_yy[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_z[j] * fl1_fz * fl1_fgb + 7.5 * pa2pb_yz_z[j] * fl1_fz * fl2_fx + 7.5 * pa2pb_y_yy[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzzz_z[j] * fl1_fz * fl1_fx + 18.0 * pa2pb_yzz_yy[j] * fl1_fz * fl1_fx + 12.0 * pa2pb_zzz_yz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yyz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_yyz[j] * fl1_fz);

                t_yzzz_yzz[j] = fl_s_0_0 * (1.125 * pa_z[j] * fl3_fx + 0.75 * pb_z[j] * fl3_fx + 2.25 * pa2pb_yz_y[j] * fl2_fx + 0.25 * pa_zzz[j] * fl2_fx + 1.5 * pa2pb_zz_z[j] * fl2_fx + 1.5 * pa2pb_y_yz[j] * fl2_fx + 0.75 * pa2pb_z_zz[j] * fl2_fx + 0.5 * pa2pb_yzzz_y[j] * fl1_fx + 3.0 * pa2pb_yzz_yz[j] * fl1_fx + 0.5 * pa2pb_zzz_zz[j] * fl1_fx + 1.5 * pa2pb_yz_yzz[j] * fl1_fx + pa2pb_yzzz_yzz[j]);

                t_yzzz_yzz[j] += fl_r_0_0 * (9.0 * pa_z[j] * fl3_fx * fl1_fz - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 0.75 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 0.5 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 6.0 * pb_z[j] * fl3_fx * fl1_fz + 22.5 * pa2pb_yz_y[j] * fl2_fx * fl1_fz + 2.5 * pa_zzz[j] * fl2_fx * fl1_fz + 15.0 * pa2pb_zz_z[j] * fl2_fx * fl1_fz - 1.5 * pa2pb_yz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 1.5 * pa2pb_yz_y[j] * fl1_fz * fl1_fga * fl1_fx - 3.0 * pa2pb_y_yz[j] * fl1_fx * fl1_fz * fl1_fga - 1.5 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_yzzz_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_y_yz[j] * fl2_fx * fl1_fz + 7.5 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_yzzz_y[j] * fl1_fz * fl1_fx + 36.0 * pa2pb_yzz_yz[j] * fl1_fz * fl1_fx + 6.0 * pa2pb_zzz_zz[j] * fl1_fx * fl1_fz - 3.0 * pa2pb_yz_yzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_yzz[j] * fl1_fz);

                t_yzzz_zzz[j] = fl_s_0_0 * (1.875 * pa_y[j] * fl3_fx + 2.25 * pa_yzz[j] * fl2_fx + 6.75 * pa2pb_yz_z[j] * fl2_fx + 2.25 * pa2pb_y_zz[j] * fl2_fx + 1.5 * pa2pb_yzzz_z[j] * fl1_fx + 4.5 * pa2pb_yzz_zz[j] * fl1_fx + 1.5 * pa2pb_yz_zzz[j] * fl1_fx + pa2pb_yzzz_zzz[j]);

                t_yzzz_zzz[j] += fl_r_0_0 * (15.0 * pa_y[j] * fl3_fx * fl1_fz - 2.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fgb - 2.25 * pa_y[j] * fl2_fx * fl1_fz * fl1_fga - 4.5 * pa_yzz[j] * fl1_fx * fl1_fz * fl1_fgb + 22.5 * pa_yzz[j] * fl1_fz * fl2_fx + 67.5 * pa2pb_yz_z[j] * fl2_fx * fl1_fz - 4.5 * pa2pb_yz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 4.5 * pa2pb_yz_z[j] * fl1_fz * fl1_fga * fl1_fx - 4.5 * pa2pb_y_zz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_yzzz_z[j] * fl1_fz * fl1_fgb + 22.5 * pa2pb_y_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_yzzz_z[j] * fl1_fz * fl1_fx + 54.0 * pa2pb_yzz_zz[j] * fl1_fz * fl1_fx - 3.0 * pa2pb_yz_zzz[j] * fl1_fz * fl1_fga + 18.0 * pa2pb_yz_zzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_yzzz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_140_145(      CMemBlock2D<double>& primBuffer,
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

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_xx = pa2pbDistances.data(646 * idx + 41);

            auto pa2pb_z_xy = pa2pbDistances.data(646 * idx + 42);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_zz_xxx = pa2pbDistances.data(646 * idx + 161);

            auto pa2pb_zz_xxy = pa2pbDistances.data(646 * idx + 162);

            auto pa2pb_zz_xxz = pa2pbDistances.data(646 * idx + 163);

            auto pa2pb_zz_xyy = pa2pbDistances.data(646 * idx + 164);

            auto pa2pb_zz_xyz = pa2pbDistances.data(646 * idx + 165);

            auto pa2pb_zzz_xx = pa2pbDistances.data(646 * idx + 345);

            auto pa2pb_zzz_xy = pa2pbDistances.data(646 * idx + 346);

            auto pa2pb_zzzz_x = pa2pbDistances.data(646 * idx + 627);

            auto pa2pb_zzzz_y = pa2pbDistances.data(646 * idx + 628);

            auto pa2pb_zzzz_z = pa2pbDistances.data(646 * idx + 629);

            auto pa2pb_zzzz_xxx = pa2pbDistances.data(646 * idx + 636);

            auto pa2pb_zzzz_xxy = pa2pbDistances.data(646 * idx + 637);

            auto pa2pb_zzzz_xxz = pa2pbDistances.data(646 * idx + 638);

            auto pa2pb_zzzz_xyy = pa2pbDistances.data(646 * idx + 639);

            auto pa2pb_zzzz_xyz = pa2pbDistances.data(646 * idx + 640);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzzz_xxx = primBuffer.data(150 * idx + 140);

            auto t_zzzz_xxy = primBuffer.data(150 * idx + 141);

            auto t_zzzz_xxz = primBuffer.data(150 * idx + 142);

            auto t_zzzz_xyy = primBuffer.data(150 * idx + 143);

            auto t_zzzz_xyz = primBuffer.data(150 * idx + 144);

            // Batch of Integrals (140,145)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_xx, pa2pb_z_xy, pa2pb_zz_x, pa2pb_zz_xxx, \
                                     pa2pb_zz_xxy, pa2pb_zz_xxz, pa2pb_zz_xyy, pa2pb_zz_xyz, pa2pb_zz_y, pa2pb_zz_z, \
                                     pa2pb_zzz_xx, pa2pb_zzz_xy, pa2pb_zzzz_x, pa2pb_zzzz_xxx, pa2pb_zzzz_xxy, \
                                     pa2pb_zzzz_xxz, pa2pb_zzzz_xyy, pa2pb_zzzz_xyz, pa2pb_zzzz_y, pa2pb_zzzz_z, pa_z, \
                                     pa_zzz, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_zzzz_xxx, t_zzzz_xxy, t_zzzz_xxz, t_zzzz_xyy, t_zzzz_xyz: VLX_ALIGN)
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

                t_zzzz_xxx[j] = fl_s_0_0 * (1.125 * pb_x[j] * fl3_fx + 4.5 * pa2pb_zz_x[j] * fl2_fx + 1.5 * pa2pb_zzzz_x[j] * fl1_fx + 0.75 * pb_xxx[j] * fl2_fx + 3.0 * pa2pb_zz_xxx[j] * fl1_fx + pa2pb_zzzz_xxx[j]);

                t_zzzz_xxx[j] += fl_r_0_0 * (-2.25 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx + 9.0 * pb_x[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_zzzz_x[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_zz_x[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zzzz_x[j] * fl1_fz * fl1_fx - 3.0 * pb_xxx[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxx[j] * fl1_fz * fl1_fga + 7.5 * pb_xxx[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_xxx[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_xxx[j] * fl1_fz);

                t_zzzz_xxy[j] = fl_s_0_0 * (0.375 * pb_y[j] * fl3_fx + 1.5 * pa2pb_zz_y[j] * fl2_fx + 0.5 * pa2pb_zzzz_y[j] * fl1_fx + 0.75 * pb_xxy[j] * fl2_fx + 3.0 * pa2pb_zz_xxy[j] * fl1_fx + pa2pb_zzzz_xxy[j]);

                t_zzzz_xxy[j] += fl_r_0_0 * (-0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx + 3.0 * pb_y[j] * fl3_fx * fl1_fz - pa2pb_zzzz_y[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_zz_y[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_zzzz_y[j] * fl1_fz * fl1_fx - 3.0 * pb_xxy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxy[j] * fl1_fz * fl1_fga + 7.5 * pb_xxy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_xxy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_xxy[j] * fl1_fz);

                t_zzzz_xxz[j] = fl_s_0_0 * (1.5 * pa_z[j] * fl3_fx + pa_zzz[j] * fl2_fx + 0.375 * pb_z[j] * fl3_fx + 1.5 * pa2pb_zz_z[j] * fl2_fx + 3.0 * pa2pb_z_xx[j] * fl2_fx + 0.5 * pa2pb_zzzz_z[j] * fl1_fx + 2.0 * pa2pb_zzz_xx[j] * fl1_fx + 0.75 * pb_xxz[j] * fl2_fx + 3.0 * pa2pb_zz_xxz[j] * fl1_fx + pa2pb_zzzz_xxz[j]);

                t_zzzz_xxz[j] += fl_r_0_0 * (-3.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_z[j] * fl3_fx * fl1_fz + 10.0 * pa_zzz[j] * fl1_fz * fl2_fx - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_xx[j] * fl1_fx * fl1_fz * fl1_fga + 3.0 * pb_z[j] * fl3_fx * fl1_fz - pa2pb_zzzz_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_zz_z[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_z_xx[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzzz_z[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_zzz_xx[j] * fl1_fz * fl1_fx - 3.0 * pb_xxz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xxz[j] * fl1_fz * fl1_fga + 7.5 * pb_xxz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_xxz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_xxz[j] * fl1_fz);

                t_zzzz_xyy[j] = fl_s_0_0 * (0.375 * pb_x[j] * fl3_fx + 1.5 * pa2pb_zz_x[j] * fl2_fx + 0.5 * pa2pb_zzzz_x[j] * fl1_fx + 0.75 * pb_xyy[j] * fl2_fx + 3.0 * pa2pb_zz_xyy[j] * fl1_fx + pa2pb_zzzz_xyy[j]);

                t_zzzz_xyy[j] += fl_r_0_0 * (-0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx + 3.0 * pb_x[j] * fl3_fx * fl1_fz - pa2pb_zzzz_x[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_zz_x[j] * fl1_fz * fl2_fx + 6.0 * pa2pb_zzzz_x[j] * fl1_fz * fl1_fx - 3.0 * pb_xyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xyy[j] * fl1_fz * fl1_fga + 7.5 * pb_xyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_xyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_xyy[j] * fl1_fz);

                t_zzzz_xyz[j] = fl_s_0_0 * (3.0 * pa2pb_z_xy[j] * fl2_fx + 2.0 * pa2pb_zzz_xy[j] * fl1_fx + 0.75 * pb_xyz[j] * fl2_fx + 3.0 * pa2pb_zz_xyz[j] * fl1_fx + pa2pb_zzzz_xyz[j]);

                t_zzzz_xyz[j] += fl_r_0_0 * (-6.0 * pa2pb_z_xy[j] * fl1_fx * fl1_fz * fl1_fga + 30.0 * pa2pb_z_xy[j] * fl2_fx * fl1_fz + 24.0 * pa2pb_zzz_xy[j] * fl1_fz * fl1_fx - 3.0 * pb_xyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xyz[j] * fl1_fz * fl1_fga + 7.5 * pb_xyz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_xyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_xyz[j] * fl1_fz);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGF_145_150(      CMemBlock2D<double>& primBuffer,
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

            auto pa_z = paDistances.data(34 * idx + 2);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to tensors product of distances R(PA)xR(PB)

            auto pa2pb_z_xz = pa2pbDistances.data(646 * idx + 43);

            auto pa2pb_z_yy = pa2pbDistances.data(646 * idx + 44);

            auto pa2pb_z_yz = pa2pbDistances.data(646 * idx + 45);

            auto pa2pb_z_zz = pa2pbDistances.data(646 * idx + 46);

            auto pa2pb_zz_x = pa2pbDistances.data(646 * idx + 152);

            auto pa2pb_zz_y = pa2pbDistances.data(646 * idx + 153);

            auto pa2pb_zz_z = pa2pbDistances.data(646 * idx + 154);

            auto pa2pb_zz_xzz = pa2pbDistances.data(646 * idx + 166);

            auto pa2pb_zz_yyy = pa2pbDistances.data(646 * idx + 167);

            auto pa2pb_zz_yyz = pa2pbDistances.data(646 * idx + 168);

            auto pa2pb_zz_yzz = pa2pbDistances.data(646 * idx + 169);

            auto pa2pb_zz_zzz = pa2pbDistances.data(646 * idx + 170);

            auto pa2pb_zzz_xz = pa2pbDistances.data(646 * idx + 347);

            auto pa2pb_zzz_yy = pa2pbDistances.data(646 * idx + 348);

            auto pa2pb_zzz_yz = pa2pbDistances.data(646 * idx + 349);

            auto pa2pb_zzz_zz = pa2pbDistances.data(646 * idx + 350);

            auto pa2pb_zzzz_x = pa2pbDistances.data(646 * idx + 627);

            auto pa2pb_zzzz_y = pa2pbDistances.data(646 * idx + 628);

            auto pa2pb_zzzz_z = pa2pbDistances.data(646 * idx + 629);

            auto pa2pb_zzzz_xzz = pa2pbDistances.data(646 * idx + 641);

            auto pa2pb_zzzz_yyy = pa2pbDistances.data(646 * idx + 642);

            auto pa2pb_zzzz_yyz = pa2pbDistances.data(646 * idx + 643);

            auto pa2pb_zzzz_yzz = pa2pbDistances.data(646 * idx + 644);

            auto pa2pb_zzzz_zzz = pa2pbDistances.data(646 * idx + 645);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

            auto t_zzzz_xzz = primBuffer.data(150 * idx + 145);

            auto t_zzzz_yyy = primBuffer.data(150 * idx + 146);

            auto t_zzzz_yyz = primBuffer.data(150 * idx + 147);

            auto t_zzzz_yzz = primBuffer.data(150 * idx + 148);

            auto t_zzzz_zzz = primBuffer.data(150 * idx + 149);

            // Batch of Integrals (145,150)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa2pb_z_xz, pa2pb_z_yy, pa2pb_z_yz, pa2pb_z_zz, pa2pb_zz_x, \
                                     pa2pb_zz_xzz, pa2pb_zz_y, pa2pb_zz_yyy, pa2pb_zz_yyz, pa2pb_zz_yzz, pa2pb_zz_z, \
                                     pa2pb_zz_zzz, pa2pb_zzz_xz, pa2pb_zzz_yy, pa2pb_zzz_yz, pa2pb_zzz_zz, pa2pb_zzzz_x, \
                                     pa2pb_zzzz_xzz, pa2pb_zzzz_y, pa2pb_zzzz_yyy, pa2pb_zzzz_yyz, pa2pb_zzzz_yzz, \
                                     pa2pb_zzzz_z, pa2pb_zzzz_zzz, pa_z, pa_zzz, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, \
                                     pb_zzz, r_0_0, s_0_0, t_zzzz_xzz, t_zzzz_yyy, t_zzzz_yyz, t_zzzz_yzz, t_zzzz_zzz: VLX_ALIGN)
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

                t_zzzz_xzz[j] = fl_s_0_0 * (1.875 * pb_x[j] * fl3_fx + 4.5 * pa2pb_zz_x[j] * fl2_fx + 6.0 * pa2pb_z_xz[j] * fl2_fx + 0.5 * pa2pb_zzzz_x[j] * fl1_fx + 4.0 * pa2pb_zzz_xz[j] * fl1_fx + 0.75 * pb_xzz[j] * fl2_fx + 3.0 * pa2pb_zz_xzz[j] * fl1_fx + pa2pb_zzzz_xzz[j]);

                t_zzzz_xzz[j] += fl_r_0_0 * (-4.5 * pb_x[j] * fl2_fx * fl1_fz * fl1_fga + 15.0 * pb_x[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_zz_x[j] * fl2_fx * fl1_fz - 0.75 * pb_x[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_x[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_x[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_z_xz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzzz_x[j] * fl1_fz * fl1_fgb + 60.0 * pa2pb_z_xz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzzz_x[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_zzz_xz[j] * fl1_fz * fl1_fx - 3.0 * pb_xzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_xzz[j] * fl1_fz * fl1_fga + 7.5 * pb_xzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_xzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_xzz[j] * fl1_fz);

                t_zzzz_yyy[j] = fl_s_0_0 * (1.125 * pb_y[j] * fl3_fx + 4.5 * pa2pb_zz_y[j] * fl2_fx + 1.5 * pa2pb_zzzz_y[j] * fl1_fx + 0.75 * pb_yyy[j] * fl2_fx + 3.0 * pa2pb_zz_yyy[j] * fl1_fx + pa2pb_zzzz_yyy[j]);

                t_zzzz_yyy[j] += fl_r_0_0 * (-2.25 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 4.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga - 9.0 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx + 9.0 * pb_y[j] * fl3_fx * fl1_fz - 3.0 * pa2pb_zzzz_y[j] * fl1_fz * fl1_fgb + 45.0 * pa2pb_zz_y[j] * fl1_fz * fl2_fx + 18.0 * pa2pb_zzzz_y[j] * fl1_fz * fl1_fx - 3.0 * pb_yyy[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_yyy[j] * fl1_fz * fl1_fga + 7.5 * pb_yyy[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_yyy[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_yyy[j] * fl1_fz);

                t_zzzz_yyz[j] = fl_s_0_0 * (1.5 * pa_z[j] * fl3_fx + pa_zzz[j] * fl2_fx + 0.375 * pb_z[j] * fl3_fx + 1.5 * pa2pb_zz_z[j] * fl2_fx + 3.0 * pa2pb_z_yy[j] * fl2_fx + 0.5 * pa2pb_zzzz_z[j] * fl1_fx + 2.0 * pa2pb_zzz_yy[j] * fl1_fx + 0.75 * pb_yyz[j] * fl2_fx + 3.0 * pa2pb_zz_yyz[j] * fl1_fx + pa2pb_zzzz_yyz[j]);

                t_zzzz_yyz[j] += fl_r_0_0 * (-3.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 2.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 12.0 * pa_z[j] * fl3_fx * fl1_fz + 10.0 * pa_zzz[j] * fl1_fz * fl2_fx - 0.75 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 1.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - 6.0 * pa2pb_z_yy[j] * fl1_fx * fl1_fz * fl1_fga + 3.0 * pb_z[j] * fl3_fx * fl1_fz - pa2pb_zzzz_z[j] * fl1_fz * fl1_fgb + 15.0 * pa2pb_zz_z[j] * fl1_fz * fl2_fx + 30.0 * pa2pb_z_yy[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzzz_z[j] * fl1_fz * fl1_fx + 24.0 * pa2pb_zzz_yy[j] * fl1_fz * fl1_fx - 3.0 * pb_yyz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_yyz[j] * fl1_fz * fl1_fga + 7.5 * pb_yyz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_yyz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_yyz[j] * fl1_fz);

                t_zzzz_yzz[j] = fl_s_0_0 * (1.875 * pb_y[j] * fl3_fx + 4.5 * pa2pb_zz_y[j] * fl2_fx + 6.0 * pa2pb_z_yz[j] * fl2_fx + 0.5 * pa2pb_zzzz_y[j] * fl1_fx + 4.0 * pa2pb_zzz_yz[j] * fl1_fx + 0.75 * pb_yzz[j] * fl2_fx + 3.0 * pa2pb_zz_yzz[j] * fl1_fx + pa2pb_zzzz_yzz[j]);

                t_zzzz_yzz[j] += fl_r_0_0 * (-4.5 * pb_y[j] * fl2_fx * fl1_fz * fl1_fga + 15.0 * pb_y[j] * fl3_fx * fl1_fz + 45.0 * pa2pb_zz_y[j] * fl2_fx * fl1_fz - 0.75 * pb_y[j] * fl2_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_y[j] * fl1_fx * fl1_fz * fl1_fgb - 3.0 * pa2pb_zz_y[j] * fl1_fz * fl1_fga * fl1_fx - 12.0 * pa2pb_z_yz[j] * fl1_fx * fl1_fz * fl1_fga - pa2pb_zzzz_y[j] * fl1_fz * fl1_fgb + 60.0 * pa2pb_z_yz[j] * fl2_fx * fl1_fz + 6.0 * pa2pb_zzzz_y[j] * fl1_fz * fl1_fx + 48.0 * pa2pb_zzz_yz[j] * fl1_fz * fl1_fx - 3.0 * pb_yzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_yzz[j] * fl1_fz * fl1_fga + 7.5 * pb_yzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_yzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_yzz[j] * fl1_fz);

                t_zzzz_zzz[j] = fl_s_0_0 * (7.5 * pa_z[j] * fl3_fx + 5.625 * pb_z[j] * fl3_fx + 3.0 * pa_zzz[j] * fl2_fx + 13.5 * pa2pb_zz_z[j] * fl2_fx + 9.0 * pa2pb_z_zz[j] * fl2_fx + 1.5 * pa2pb_zzzz_z[j] * fl1_fx + 6.0 * pa2pb_zzz_zz[j] * fl1_fx + 0.75 * pb_zzz[j] * fl2_fx + 3.0 * pa2pb_zz_zzz[j] * fl1_fx + pa2pb_zzzz_zzz[j]);

                t_zzzz_zzz[j] += fl_r_0_0 * (60.0 * pa_z[j] * fl3_fx * fl1_fz - 9.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa_z[j] * fl2_fx * fl1_fz * fl1_fga - 13.5 * pb_z[j] * fl2_fx * fl1_fz * fl1_fga - 6.0 * pa_zzz[j] * fl1_fx * fl1_fz * fl1_fgb + 45.0 * pb_z[j] * fl3_fx * fl1_fz + 30.0 * pa_zzz[j] * fl1_fz * fl2_fx + 135.0 * pa2pb_zz_z[j] * fl2_fx * fl1_fz - 2.25 * pb_z[j] * fl2_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_z[j] * fl1_fx * fl1_fz * fl1_fgb - 9.0 * pa2pb_zz_z[j] * fl1_fz * fl1_fga * fl1_fx - 18.0 * pa2pb_z_zz[j] * fl1_fx * fl1_fz * fl1_fga - 3.0 * pa2pb_zzzz_z[j] * fl1_fz * fl1_fgb + 90.0 * pa2pb_z_zz[j] * fl2_fx * fl1_fz + 18.0 * pa2pb_zzzz_z[j] * fl1_fz * fl1_fx + 72.0 * pa2pb_zzz_zz[j] * fl1_fz * fl1_fx - 3.0 * pb_zzz[j] * fl1_fx * fl1_fz * fl1_fga - 6.0 * pa2pb_zz_zzz[j] * fl1_fz * fl1_fga + 7.5 * pb_zzz[j] * fl2_fx * fl1_fz + 36.0 * pa2pb_zz_zzz[j] * fl1_fz * fl1_fx + 14.0 * pa2pb_zzzz_zzz[j] * fl1_fz);
            }

            idx++;
        }
    }


} // kinrecfunc namespace

