//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionBRRRecFuncForDXYY.hpp"

#include "AngularMomentum.hpp"

namespace eribrrfunc { // eribrrfunc namespace

    void
    compElectronRepulsionForDDXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 2, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_2_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 2, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 3, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_x_xx_0 = braBuffer.data(pidx_g_1_2_m0 + i); 

            auto tg_x_xy_0 = braBuffer.data(pidx_g_1_2_m0 + kcomp + i); 

            auto tg_x_xz_0 = braBuffer.data(pidx_g_1_2_m0 + 2 * kcomp + i); 

            auto tg_x_yy_0 = braBuffer.data(pidx_g_1_2_m0 + 3 * kcomp + i); 

            auto tg_x_yz_0 = braBuffer.data(pidx_g_1_2_m0 + 4 * kcomp + i); 

            auto tg_x_zz_0 = braBuffer.data(pidx_g_1_2_m0 + 5 * kcomp + i); 

            auto tg_y_xx_0 = braBuffer.data(pidx_g_1_2_m0 + 6 * kcomp + i); 

            auto tg_y_xy_0 = braBuffer.data(pidx_g_1_2_m0 + 7 * kcomp + i); 

            auto tg_y_xz_0 = braBuffer.data(pidx_g_1_2_m0 + 8 * kcomp + i); 

            auto tg_y_yy_0 = braBuffer.data(pidx_g_1_2_m0 + 9 * kcomp + i); 

            auto tg_y_yz_0 = braBuffer.data(pidx_g_1_2_m0 + 10 * kcomp + i); 

            auto tg_y_zz_0 = braBuffer.data(pidx_g_1_2_m0 + 11 * kcomp + i); 

            auto tg_z_xx_0 = braBuffer.data(pidx_g_1_2_m0 + 12 * kcomp + i); 

            auto tg_z_xy_0 = braBuffer.data(pidx_g_1_2_m0 + 13 * kcomp + i); 

            auto tg_z_xz_0 = braBuffer.data(pidx_g_1_2_m0 + 14 * kcomp + i); 

            auto tg_z_yy_0 = braBuffer.data(pidx_g_1_2_m0 + 15 * kcomp + i); 

            auto tg_z_yz_0 = braBuffer.data(pidx_g_1_2_m0 + 16 * kcomp + i); 

            auto tg_z_zz_0 = braBuffer.data(pidx_g_1_2_m0 + 17 * kcomp + i); 

            auto tg_x_xxx_0 = braBuffer.data(pidx_g_1_3_m0 + i); 

            auto tg_x_xxy_0 = braBuffer.data(pidx_g_1_3_m0 + kcomp + i); 

            auto tg_x_xxz_0 = braBuffer.data(pidx_g_1_3_m0 + 2 * kcomp + i); 

            auto tg_x_xyy_0 = braBuffer.data(pidx_g_1_3_m0 + 3 * kcomp + i); 

            auto tg_x_xyz_0 = braBuffer.data(pidx_g_1_3_m0 + 4 * kcomp + i); 

            auto tg_x_xzz_0 = braBuffer.data(pidx_g_1_3_m0 + 5 * kcomp + i); 

            auto tg_y_xxx_0 = braBuffer.data(pidx_g_1_3_m0 + 10 * kcomp + i); 

            auto tg_y_xxy_0 = braBuffer.data(pidx_g_1_3_m0 + 11 * kcomp + i); 

            auto tg_y_xxz_0 = braBuffer.data(pidx_g_1_3_m0 + 12 * kcomp + i); 

            auto tg_y_xyy_0 = braBuffer.data(pidx_g_1_3_m0 + 13 * kcomp + i); 

            auto tg_y_xyz_0 = braBuffer.data(pidx_g_1_3_m0 + 14 * kcomp + i); 

            auto tg_y_xzz_0 = braBuffer.data(pidx_g_1_3_m0 + 15 * kcomp + i); 

            auto tg_y_yyy_0 = braBuffer.data(pidx_g_1_3_m0 + 16 * kcomp + i); 

            auto tg_y_yyz_0 = braBuffer.data(pidx_g_1_3_m0 + 17 * kcomp + i); 

            auto tg_y_yzz_0 = braBuffer.data(pidx_g_1_3_m0 + 18 * kcomp + i); 

            auto tg_z_xxx_0 = braBuffer.data(pidx_g_1_3_m0 + 20 * kcomp + i); 

            auto tg_z_xxy_0 = braBuffer.data(pidx_g_1_3_m0 + 21 * kcomp + i); 

            auto tg_z_xxz_0 = braBuffer.data(pidx_g_1_3_m0 + 22 * kcomp + i); 

            auto tg_z_xyy_0 = braBuffer.data(pidx_g_1_3_m0 + 23 * kcomp + i); 

            auto tg_z_xyz_0 = braBuffer.data(pidx_g_1_3_m0 + 24 * kcomp + i); 

            auto tg_z_xzz_0 = braBuffer.data(pidx_g_1_3_m0 + 25 * kcomp + i); 

            auto tg_z_yyy_0 = braBuffer.data(pidx_g_1_3_m0 + 26 * kcomp + i); 

            auto tg_z_yyz_0 = braBuffer.data(pidx_g_1_3_m0 + 27 * kcomp + i); 

            auto tg_z_yzz_0 = braBuffer.data(pidx_g_1_3_m0 + 28 * kcomp + i); 

            auto tg_z_zzz_0 = braBuffer.data(pidx_g_1_3_m0 + 29 * kcomp + i); 

            // set up pointers to integrals

            auto tg_xx_xx_0 = braBuffer.data(pidx_g_2_2_m0 + i); 

            auto tg_xx_xy_0 = braBuffer.data(pidx_g_2_2_m0 + kcomp + i); 

            auto tg_xx_xz_0 = braBuffer.data(pidx_g_2_2_m0 + 2 * kcomp + i); 

            auto tg_xx_yy_0 = braBuffer.data(pidx_g_2_2_m0 + 3 * kcomp + i); 

            auto tg_xx_yz_0 = braBuffer.data(pidx_g_2_2_m0 + 4 * kcomp + i); 

            auto tg_xx_zz_0 = braBuffer.data(pidx_g_2_2_m0 + 5 * kcomp + i); 

            auto tg_xy_xx_0 = braBuffer.data(pidx_g_2_2_m0 + 6 * kcomp + i); 

            auto tg_xy_xy_0 = braBuffer.data(pidx_g_2_2_m0 + 7 * kcomp + i); 

            auto tg_xy_xz_0 = braBuffer.data(pidx_g_2_2_m0 + 8 * kcomp + i); 

            auto tg_xy_yy_0 = braBuffer.data(pidx_g_2_2_m0 + 9 * kcomp + i); 

            auto tg_xy_yz_0 = braBuffer.data(pidx_g_2_2_m0 + 10 * kcomp + i); 

            auto tg_xy_zz_0 = braBuffer.data(pidx_g_2_2_m0 + 11 * kcomp + i); 

            auto tg_xz_xx_0 = braBuffer.data(pidx_g_2_2_m0 + 12 * kcomp + i); 

            auto tg_xz_xy_0 = braBuffer.data(pidx_g_2_2_m0 + 13 * kcomp + i); 

            auto tg_xz_xz_0 = braBuffer.data(pidx_g_2_2_m0 + 14 * kcomp + i); 

            auto tg_xz_yy_0 = braBuffer.data(pidx_g_2_2_m0 + 15 * kcomp + i); 

            auto tg_xz_yz_0 = braBuffer.data(pidx_g_2_2_m0 + 16 * kcomp + i); 

            auto tg_xz_zz_0 = braBuffer.data(pidx_g_2_2_m0 + 17 * kcomp + i); 

            auto tg_yy_xx_0 = braBuffer.data(pidx_g_2_2_m0 + 18 * kcomp + i); 

            auto tg_yy_xy_0 = braBuffer.data(pidx_g_2_2_m0 + 19 * kcomp + i); 

            auto tg_yy_xz_0 = braBuffer.data(pidx_g_2_2_m0 + 20 * kcomp + i); 

            auto tg_yy_yy_0 = braBuffer.data(pidx_g_2_2_m0 + 21 * kcomp + i); 

            auto tg_yy_yz_0 = braBuffer.data(pidx_g_2_2_m0 + 22 * kcomp + i); 

            auto tg_yy_zz_0 = braBuffer.data(pidx_g_2_2_m0 + 23 * kcomp + i); 

            auto tg_yz_xx_0 = braBuffer.data(pidx_g_2_2_m0 + 24 * kcomp + i); 

            auto tg_yz_xy_0 = braBuffer.data(pidx_g_2_2_m0 + 25 * kcomp + i); 

            auto tg_yz_xz_0 = braBuffer.data(pidx_g_2_2_m0 + 26 * kcomp + i); 

            auto tg_yz_yy_0 = braBuffer.data(pidx_g_2_2_m0 + 27 * kcomp + i); 

            auto tg_yz_yz_0 = braBuffer.data(pidx_g_2_2_m0 + 28 * kcomp + i); 

            auto tg_yz_zz_0 = braBuffer.data(pidx_g_2_2_m0 + 29 * kcomp + i); 

            auto tg_zz_xx_0 = braBuffer.data(pidx_g_2_2_m0 + 30 * kcomp + i); 

            auto tg_zz_xy_0 = braBuffer.data(pidx_g_2_2_m0 + 31 * kcomp + i); 

            auto tg_zz_xz_0 = braBuffer.data(pidx_g_2_2_m0 + 32 * kcomp + i); 

            auto tg_zz_yy_0 = braBuffer.data(pidx_g_2_2_m0 + 33 * kcomp + i); 

            auto tg_zz_yz_0 = braBuffer.data(pidx_g_2_2_m0 + 34 * kcomp + i); 

            auto tg_zz_zz_0 = braBuffer.data(pidx_g_2_2_m0 + 35 * kcomp + i); 

            #pragma omp simd aligned(tg_x_xx_0, tg_x_xxx_0, tg_x_xxy_0, tg_x_xxz_0, tg_x_xy_0, tg_x_xyy_0, \
                                     tg_x_xyz_0, tg_x_xz_0, tg_x_xzz_0, tg_x_yy_0, tg_x_yz_0, tg_x_zz_0, tg_xx_xx_0, \
                                     tg_xx_xy_0, tg_xx_xz_0, tg_xx_yy_0, tg_xx_yz_0, tg_xx_zz_0, tg_xy_xx_0, tg_xy_xy_0, \
                                     tg_xy_xz_0, tg_xy_yy_0, tg_xy_yz_0, tg_xy_zz_0, tg_xz_xx_0, tg_xz_xy_0, tg_xz_xz_0, \
                                     tg_xz_yy_0, tg_xz_yz_0, tg_xz_zz_0, tg_y_xx_0, tg_y_xxx_0, tg_y_xxy_0, tg_y_xxz_0, \
                                     tg_y_xy_0, tg_y_xyy_0, tg_y_xyz_0, tg_y_xz_0, tg_y_xzz_0, tg_y_yy_0, tg_y_yyy_0, \
                                     tg_y_yyz_0, tg_y_yz_0, tg_y_yzz_0, tg_y_zz_0, tg_yy_xx_0, tg_yy_xy_0, tg_yy_xz_0, \
                                     tg_yy_yy_0, tg_yy_yz_0, tg_yy_zz_0, tg_yz_xx_0, tg_yz_xy_0, tg_yz_xz_0, tg_yz_yy_0, \
                                     tg_yz_yz_0, tg_yz_zz_0, tg_z_xx_0, tg_z_xxx_0, tg_z_xxy_0, tg_z_xxz_0, tg_z_xy_0, \
                                     tg_z_xyy_0, tg_z_xyz_0, tg_z_xz_0, tg_z_xzz_0, tg_z_yy_0, tg_z_yyy_0, tg_z_yyz_0, \
                                     tg_z_yz_0, tg_z_yzz_0, tg_z_zz_0, tg_z_zzz_0, tg_zz_xx_0, tg_zz_xy_0, tg_zz_xz_0, \
                                     tg_zz_yy_0, tg_zz_yz_0, tg_zz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xx_xx_0[j] = -ab_x * tg_x_xx_0[j] + tg_x_xxx_0[j];

                tg_xx_xy_0[j] = -ab_x * tg_x_xy_0[j] + tg_x_xxy_0[j];

                tg_xx_xz_0[j] = -ab_x * tg_x_xz_0[j] + tg_x_xxz_0[j];

                tg_xx_yy_0[j] = -ab_x * tg_x_yy_0[j] + tg_x_xyy_0[j];

                tg_xx_yz_0[j] = -ab_x * tg_x_yz_0[j] + tg_x_xyz_0[j];

                tg_xx_zz_0[j] = -ab_x * tg_x_zz_0[j] + tg_x_xzz_0[j];

                tg_xy_xx_0[j] = -ab_x * tg_y_xx_0[j] + tg_y_xxx_0[j];

                tg_xy_xy_0[j] = -ab_x * tg_y_xy_0[j] + tg_y_xxy_0[j];

                tg_xy_xz_0[j] = -ab_x * tg_y_xz_0[j] + tg_y_xxz_0[j];

                tg_xy_yy_0[j] = -ab_x * tg_y_yy_0[j] + tg_y_xyy_0[j];

                tg_xy_yz_0[j] = -ab_x * tg_y_yz_0[j] + tg_y_xyz_0[j];

                tg_xy_zz_0[j] = -ab_x * tg_y_zz_0[j] + tg_y_xzz_0[j];

                tg_xz_xx_0[j] = -ab_x * tg_z_xx_0[j] + tg_z_xxx_0[j];

                tg_xz_xy_0[j] = -ab_x * tg_z_xy_0[j] + tg_z_xxy_0[j];

                tg_xz_xz_0[j] = -ab_x * tg_z_xz_0[j] + tg_z_xxz_0[j];

                tg_xz_yy_0[j] = -ab_x * tg_z_yy_0[j] + tg_z_xyy_0[j];

                tg_xz_yz_0[j] = -ab_x * tg_z_yz_0[j] + tg_z_xyz_0[j];

                tg_xz_zz_0[j] = -ab_x * tg_z_zz_0[j] + tg_z_xzz_0[j];

                tg_yy_xx_0[j] = -ab_y * tg_y_xx_0[j] + tg_y_xxy_0[j];

                tg_yy_xy_0[j] = -ab_y * tg_y_xy_0[j] + tg_y_xyy_0[j];

                tg_yy_xz_0[j] = -ab_y * tg_y_xz_0[j] + tg_y_xyz_0[j];

                tg_yy_yy_0[j] = -ab_y * tg_y_yy_0[j] + tg_y_yyy_0[j];

                tg_yy_yz_0[j] = -ab_y * tg_y_yz_0[j] + tg_y_yyz_0[j];

                tg_yy_zz_0[j] = -ab_y * tg_y_zz_0[j] + tg_y_yzz_0[j];

                tg_yz_xx_0[j] = -ab_y * tg_z_xx_0[j] + tg_z_xxy_0[j];

                tg_yz_xy_0[j] = -ab_y * tg_z_xy_0[j] + tg_z_xyy_0[j];

                tg_yz_xz_0[j] = -ab_y * tg_z_xz_0[j] + tg_z_xyz_0[j];

                tg_yz_yy_0[j] = -ab_y * tg_z_yy_0[j] + tg_z_yyy_0[j];

                tg_yz_yz_0[j] = -ab_y * tg_z_yz_0[j] + tg_z_yyz_0[j];

                tg_yz_zz_0[j] = -ab_y * tg_z_zz_0[j] + tg_z_yzz_0[j];

                tg_zz_xx_0[j] = -ab_z * tg_z_xx_0[j] + tg_z_xxz_0[j];

                tg_zz_xy_0[j] = -ab_z * tg_z_xy_0[j] + tg_z_xyz_0[j];

                tg_zz_xz_0[j] = -ab_z * tg_z_xz_0[j] + tg_z_xzz_0[j];

                tg_zz_yy_0[j] = -ab_z * tg_z_yy_0[j] + tg_z_yyz_0[j];

                tg_zz_yz_0[j] = -ab_z * tg_z_yz_0[j] + tg_z_yzz_0[j];

                tg_zz_zz_0[j] = -ab_z * tg_z_zz_0[j] + tg_z_zzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForDFXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 3, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 3, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_x_xxx_0 = braBuffer.data(pidx_g_1_3_m0 + i); 

            auto tg_x_xxy_0 = braBuffer.data(pidx_g_1_3_m0 + kcomp + i); 

            auto tg_x_xxz_0 = braBuffer.data(pidx_g_1_3_m0 + 2 * kcomp + i); 

            auto tg_x_xyy_0 = braBuffer.data(pidx_g_1_3_m0 + 3 * kcomp + i); 

            auto tg_x_xyz_0 = braBuffer.data(pidx_g_1_3_m0 + 4 * kcomp + i); 

            auto tg_x_xzz_0 = braBuffer.data(pidx_g_1_3_m0 + 5 * kcomp + i); 

            auto tg_x_yyy_0 = braBuffer.data(pidx_g_1_3_m0 + 6 * kcomp + i); 

            auto tg_x_yyz_0 = braBuffer.data(pidx_g_1_3_m0 + 7 * kcomp + i); 

            auto tg_x_yzz_0 = braBuffer.data(pidx_g_1_3_m0 + 8 * kcomp + i); 

            auto tg_x_zzz_0 = braBuffer.data(pidx_g_1_3_m0 + 9 * kcomp + i); 

            auto tg_y_xxx_0 = braBuffer.data(pidx_g_1_3_m0 + 10 * kcomp + i); 

            auto tg_y_xxy_0 = braBuffer.data(pidx_g_1_3_m0 + 11 * kcomp + i); 

            auto tg_y_xxz_0 = braBuffer.data(pidx_g_1_3_m0 + 12 * kcomp + i); 

            auto tg_y_xyy_0 = braBuffer.data(pidx_g_1_3_m0 + 13 * kcomp + i); 

            auto tg_y_xyz_0 = braBuffer.data(pidx_g_1_3_m0 + 14 * kcomp + i); 

            auto tg_y_xzz_0 = braBuffer.data(pidx_g_1_3_m0 + 15 * kcomp + i); 

            auto tg_y_yyy_0 = braBuffer.data(pidx_g_1_3_m0 + 16 * kcomp + i); 

            auto tg_y_yyz_0 = braBuffer.data(pidx_g_1_3_m0 + 17 * kcomp + i); 

            auto tg_y_yzz_0 = braBuffer.data(pidx_g_1_3_m0 + 18 * kcomp + i); 

            auto tg_y_zzz_0 = braBuffer.data(pidx_g_1_3_m0 + 19 * kcomp + i); 

            auto tg_z_xxx_0 = braBuffer.data(pidx_g_1_3_m0 + 20 * kcomp + i); 

            auto tg_z_xxy_0 = braBuffer.data(pidx_g_1_3_m0 + 21 * kcomp + i); 

            auto tg_z_xxz_0 = braBuffer.data(pidx_g_1_3_m0 + 22 * kcomp + i); 

            auto tg_z_xyy_0 = braBuffer.data(pidx_g_1_3_m0 + 23 * kcomp + i); 

            auto tg_z_xyz_0 = braBuffer.data(pidx_g_1_3_m0 + 24 * kcomp + i); 

            auto tg_z_xzz_0 = braBuffer.data(pidx_g_1_3_m0 + 25 * kcomp + i); 

            auto tg_z_yyy_0 = braBuffer.data(pidx_g_1_3_m0 + 26 * kcomp + i); 

            auto tg_z_yyz_0 = braBuffer.data(pidx_g_1_3_m0 + 27 * kcomp + i); 

            auto tg_z_yzz_0 = braBuffer.data(pidx_g_1_3_m0 + 28 * kcomp + i); 

            auto tg_z_zzz_0 = braBuffer.data(pidx_g_1_3_m0 + 29 * kcomp + i); 

            auto tg_x_xxxx_0 = braBuffer.data(pidx_g_1_4_m0 + i); 

            auto tg_x_xxxy_0 = braBuffer.data(pidx_g_1_4_m0 + kcomp + i); 

            auto tg_x_xxxz_0 = braBuffer.data(pidx_g_1_4_m0 + 2 * kcomp + i); 

            auto tg_x_xxyy_0 = braBuffer.data(pidx_g_1_4_m0 + 3 * kcomp + i); 

            auto tg_x_xxyz_0 = braBuffer.data(pidx_g_1_4_m0 + 4 * kcomp + i); 

            auto tg_x_xxzz_0 = braBuffer.data(pidx_g_1_4_m0 + 5 * kcomp + i); 

            auto tg_x_xyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 6 * kcomp + i); 

            auto tg_x_xyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 7 * kcomp + i); 

            auto tg_x_xyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 8 * kcomp + i); 

            auto tg_x_xzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 9 * kcomp + i); 

            auto tg_y_xxxx_0 = braBuffer.data(pidx_g_1_4_m0 + 15 * kcomp + i); 

            auto tg_y_xxxy_0 = braBuffer.data(pidx_g_1_4_m0 + 16 * kcomp + i); 

            auto tg_y_xxxz_0 = braBuffer.data(pidx_g_1_4_m0 + 17 * kcomp + i); 

            auto tg_y_xxyy_0 = braBuffer.data(pidx_g_1_4_m0 + 18 * kcomp + i); 

            auto tg_y_xxyz_0 = braBuffer.data(pidx_g_1_4_m0 + 19 * kcomp + i); 

            auto tg_y_xxzz_0 = braBuffer.data(pidx_g_1_4_m0 + 20 * kcomp + i); 

            auto tg_y_xyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 21 * kcomp + i); 

            auto tg_y_xyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 22 * kcomp + i); 

            auto tg_y_xyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 23 * kcomp + i); 

            auto tg_y_xzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 24 * kcomp + i); 

            auto tg_y_yyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 25 * kcomp + i); 

            auto tg_y_yyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 26 * kcomp + i); 

            auto tg_y_yyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 27 * kcomp + i); 

            auto tg_y_yzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 28 * kcomp + i); 

            auto tg_z_xxxx_0 = braBuffer.data(pidx_g_1_4_m0 + 30 * kcomp + i); 

            auto tg_z_xxxy_0 = braBuffer.data(pidx_g_1_4_m0 + 31 * kcomp + i); 

            auto tg_z_xxxz_0 = braBuffer.data(pidx_g_1_4_m0 + 32 * kcomp + i); 

            auto tg_z_xxyy_0 = braBuffer.data(pidx_g_1_4_m0 + 33 * kcomp + i); 

            auto tg_z_xxyz_0 = braBuffer.data(pidx_g_1_4_m0 + 34 * kcomp + i); 

            auto tg_z_xxzz_0 = braBuffer.data(pidx_g_1_4_m0 + 35 * kcomp + i); 

            auto tg_z_xyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 36 * kcomp + i); 

            auto tg_z_xyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 37 * kcomp + i); 

            auto tg_z_xyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 38 * kcomp + i); 

            auto tg_z_xzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 39 * kcomp + i); 

            auto tg_z_yyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 40 * kcomp + i); 

            auto tg_z_yyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 41 * kcomp + i); 

            auto tg_z_yyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 42 * kcomp + i); 

            auto tg_z_yzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 43 * kcomp + i); 

            auto tg_z_zzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 44 * kcomp + i); 

            // set up pointers to integrals

            auto tg_xx_xxx_0 = braBuffer.data(pidx_g_2_3_m0 + i); 

            auto tg_xx_xxy_0 = braBuffer.data(pidx_g_2_3_m0 + kcomp + i); 

            auto tg_xx_xxz_0 = braBuffer.data(pidx_g_2_3_m0 + 2 * kcomp + i); 

            auto tg_xx_xyy_0 = braBuffer.data(pidx_g_2_3_m0 + 3 * kcomp + i); 

            auto tg_xx_xyz_0 = braBuffer.data(pidx_g_2_3_m0 + 4 * kcomp + i); 

            auto tg_xx_xzz_0 = braBuffer.data(pidx_g_2_3_m0 + 5 * kcomp + i); 

            auto tg_xx_yyy_0 = braBuffer.data(pidx_g_2_3_m0 + 6 * kcomp + i); 

            auto tg_xx_yyz_0 = braBuffer.data(pidx_g_2_3_m0 + 7 * kcomp + i); 

            auto tg_xx_yzz_0 = braBuffer.data(pidx_g_2_3_m0 + 8 * kcomp + i); 

            auto tg_xx_zzz_0 = braBuffer.data(pidx_g_2_3_m0 + 9 * kcomp + i); 

            auto tg_xy_xxx_0 = braBuffer.data(pidx_g_2_3_m0 + 10 * kcomp + i); 

            auto tg_xy_xxy_0 = braBuffer.data(pidx_g_2_3_m0 + 11 * kcomp + i); 

            auto tg_xy_xxz_0 = braBuffer.data(pidx_g_2_3_m0 + 12 * kcomp + i); 

            auto tg_xy_xyy_0 = braBuffer.data(pidx_g_2_3_m0 + 13 * kcomp + i); 

            auto tg_xy_xyz_0 = braBuffer.data(pidx_g_2_3_m0 + 14 * kcomp + i); 

            auto tg_xy_xzz_0 = braBuffer.data(pidx_g_2_3_m0 + 15 * kcomp + i); 

            auto tg_xy_yyy_0 = braBuffer.data(pidx_g_2_3_m0 + 16 * kcomp + i); 

            auto tg_xy_yyz_0 = braBuffer.data(pidx_g_2_3_m0 + 17 * kcomp + i); 

            auto tg_xy_yzz_0 = braBuffer.data(pidx_g_2_3_m0 + 18 * kcomp + i); 

            auto tg_xy_zzz_0 = braBuffer.data(pidx_g_2_3_m0 + 19 * kcomp + i); 

            auto tg_xz_xxx_0 = braBuffer.data(pidx_g_2_3_m0 + 20 * kcomp + i); 

            auto tg_xz_xxy_0 = braBuffer.data(pidx_g_2_3_m0 + 21 * kcomp + i); 

            auto tg_xz_xxz_0 = braBuffer.data(pidx_g_2_3_m0 + 22 * kcomp + i); 

            auto tg_xz_xyy_0 = braBuffer.data(pidx_g_2_3_m0 + 23 * kcomp + i); 

            auto tg_xz_xyz_0 = braBuffer.data(pidx_g_2_3_m0 + 24 * kcomp + i); 

            auto tg_xz_xzz_0 = braBuffer.data(pidx_g_2_3_m0 + 25 * kcomp + i); 

            auto tg_xz_yyy_0 = braBuffer.data(pidx_g_2_3_m0 + 26 * kcomp + i); 

            auto tg_xz_yyz_0 = braBuffer.data(pidx_g_2_3_m0 + 27 * kcomp + i); 

            auto tg_xz_yzz_0 = braBuffer.data(pidx_g_2_3_m0 + 28 * kcomp + i); 

            auto tg_xz_zzz_0 = braBuffer.data(pidx_g_2_3_m0 + 29 * kcomp + i); 

            auto tg_yy_xxx_0 = braBuffer.data(pidx_g_2_3_m0 + 30 * kcomp + i); 

            auto tg_yy_xxy_0 = braBuffer.data(pidx_g_2_3_m0 + 31 * kcomp + i); 

            auto tg_yy_xxz_0 = braBuffer.data(pidx_g_2_3_m0 + 32 * kcomp + i); 

            auto tg_yy_xyy_0 = braBuffer.data(pidx_g_2_3_m0 + 33 * kcomp + i); 

            auto tg_yy_xyz_0 = braBuffer.data(pidx_g_2_3_m0 + 34 * kcomp + i); 

            auto tg_yy_xzz_0 = braBuffer.data(pidx_g_2_3_m0 + 35 * kcomp + i); 

            auto tg_yy_yyy_0 = braBuffer.data(pidx_g_2_3_m0 + 36 * kcomp + i); 

            auto tg_yy_yyz_0 = braBuffer.data(pidx_g_2_3_m0 + 37 * kcomp + i); 

            auto tg_yy_yzz_0 = braBuffer.data(pidx_g_2_3_m0 + 38 * kcomp + i); 

            auto tg_yy_zzz_0 = braBuffer.data(pidx_g_2_3_m0 + 39 * kcomp + i); 

            auto tg_yz_xxx_0 = braBuffer.data(pidx_g_2_3_m0 + 40 * kcomp + i); 

            auto tg_yz_xxy_0 = braBuffer.data(pidx_g_2_3_m0 + 41 * kcomp + i); 

            auto tg_yz_xxz_0 = braBuffer.data(pidx_g_2_3_m0 + 42 * kcomp + i); 

            auto tg_yz_xyy_0 = braBuffer.data(pidx_g_2_3_m0 + 43 * kcomp + i); 

            auto tg_yz_xyz_0 = braBuffer.data(pidx_g_2_3_m0 + 44 * kcomp + i); 

            auto tg_yz_xzz_0 = braBuffer.data(pidx_g_2_3_m0 + 45 * kcomp + i); 

            auto tg_yz_yyy_0 = braBuffer.data(pidx_g_2_3_m0 + 46 * kcomp + i); 

            auto tg_yz_yyz_0 = braBuffer.data(pidx_g_2_3_m0 + 47 * kcomp + i); 

            auto tg_yz_yzz_0 = braBuffer.data(pidx_g_2_3_m0 + 48 * kcomp + i); 

            auto tg_yz_zzz_0 = braBuffer.data(pidx_g_2_3_m0 + 49 * kcomp + i); 

            auto tg_zz_xxx_0 = braBuffer.data(pidx_g_2_3_m0 + 50 * kcomp + i); 

            auto tg_zz_xxy_0 = braBuffer.data(pidx_g_2_3_m0 + 51 * kcomp + i); 

            auto tg_zz_xxz_0 = braBuffer.data(pidx_g_2_3_m0 + 52 * kcomp + i); 

            auto tg_zz_xyy_0 = braBuffer.data(pidx_g_2_3_m0 + 53 * kcomp + i); 

            auto tg_zz_xyz_0 = braBuffer.data(pidx_g_2_3_m0 + 54 * kcomp + i); 

            auto tg_zz_xzz_0 = braBuffer.data(pidx_g_2_3_m0 + 55 * kcomp + i); 

            auto tg_zz_yyy_0 = braBuffer.data(pidx_g_2_3_m0 + 56 * kcomp + i); 

            auto tg_zz_yyz_0 = braBuffer.data(pidx_g_2_3_m0 + 57 * kcomp + i); 

            auto tg_zz_yzz_0 = braBuffer.data(pidx_g_2_3_m0 + 58 * kcomp + i); 

            auto tg_zz_zzz_0 = braBuffer.data(pidx_g_2_3_m0 + 59 * kcomp + i); 

            #pragma omp simd aligned(tg_x_xxx_0, tg_x_xxxx_0, tg_x_xxxy_0, tg_x_xxxz_0, tg_x_xxy_0, \
                                     tg_x_xxyy_0, tg_x_xxyz_0, tg_x_xxz_0, tg_x_xxzz_0, tg_x_xyy_0, tg_x_xyyy_0, \
                                     tg_x_xyyz_0, tg_x_xyz_0, tg_x_xyzz_0, tg_x_xzz_0, tg_x_xzzz_0, tg_x_yyy_0, \
                                     tg_x_yyz_0, tg_x_yzz_0, tg_x_zzz_0, tg_xx_xxx_0, tg_xx_xxy_0, tg_xx_xxz_0, \
                                     tg_xx_xyy_0, tg_xx_xyz_0, tg_xx_xzz_0, tg_xx_yyy_0, tg_xx_yyz_0, tg_xx_yzz_0, \
                                     tg_xx_zzz_0, tg_xy_xxx_0, tg_xy_xxy_0, tg_xy_xxz_0, tg_xy_xyy_0, tg_xy_xyz_0, \
                                     tg_xy_xzz_0, tg_xy_yyy_0, tg_xy_yyz_0, tg_xy_yzz_0, tg_xy_zzz_0, tg_xz_xxx_0, \
                                     tg_xz_xxy_0, tg_xz_xxz_0, tg_xz_xyy_0, tg_xz_xyz_0, tg_xz_xzz_0, tg_xz_yyy_0, \
                                     tg_xz_yyz_0, tg_xz_yzz_0, tg_xz_zzz_0, tg_y_xxx_0, tg_y_xxxx_0, tg_y_xxxy_0, \
                                     tg_y_xxxz_0, tg_y_xxy_0, tg_y_xxyy_0, tg_y_xxyz_0, tg_y_xxz_0, tg_y_xxzz_0, \
                                     tg_y_xyy_0, tg_y_xyyy_0, tg_y_xyyz_0, tg_y_xyz_0, tg_y_xyzz_0, tg_y_xzz_0, \
                                     tg_y_xzzz_0, tg_y_yyy_0, tg_y_yyyy_0, tg_y_yyyz_0, tg_y_yyz_0, tg_y_yyzz_0, \
                                     tg_y_yzz_0, tg_y_yzzz_0, tg_y_zzz_0, tg_yy_xxx_0, tg_yy_xxy_0, tg_yy_xxz_0, \
                                     tg_yy_xyy_0, tg_yy_xyz_0, tg_yy_xzz_0, tg_yy_yyy_0, tg_yy_yyz_0, tg_yy_yzz_0, \
                                     tg_yy_zzz_0, tg_yz_xxx_0, tg_yz_xxy_0, tg_yz_xxz_0, tg_yz_xyy_0, tg_yz_xyz_0, \
                                     tg_yz_xzz_0, tg_yz_yyy_0, tg_yz_yyz_0, tg_yz_yzz_0, tg_yz_zzz_0, tg_z_xxx_0, \
                                     tg_z_xxxx_0, tg_z_xxxy_0, tg_z_xxxz_0, tg_z_xxy_0, tg_z_xxyy_0, tg_z_xxyz_0, \
                                     tg_z_xxz_0, tg_z_xxzz_0, tg_z_xyy_0, tg_z_xyyy_0, tg_z_xyyz_0, tg_z_xyz_0, \
                                     tg_z_xyzz_0, tg_z_xzz_0, tg_z_xzzz_0, tg_z_yyy_0, tg_z_yyyy_0, tg_z_yyyz_0, \
                                     tg_z_yyz_0, tg_z_yyzz_0, tg_z_yzz_0, tg_z_yzzz_0, tg_z_zzz_0, tg_z_zzzz_0, \
                                     tg_zz_xxx_0, tg_zz_xxy_0, tg_zz_xxz_0, tg_zz_xyy_0, tg_zz_xyz_0, tg_zz_xzz_0, \
                                     tg_zz_yyy_0, tg_zz_yyz_0, tg_zz_yzz_0, tg_zz_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xx_xxx_0[j] = -ab_x * tg_x_xxx_0[j] + tg_x_xxxx_0[j];

                tg_xx_xxy_0[j] = -ab_x * tg_x_xxy_0[j] + tg_x_xxxy_0[j];

                tg_xx_xxz_0[j] = -ab_x * tg_x_xxz_0[j] + tg_x_xxxz_0[j];

                tg_xx_xyy_0[j] = -ab_x * tg_x_xyy_0[j] + tg_x_xxyy_0[j];

                tg_xx_xyz_0[j] = -ab_x * tg_x_xyz_0[j] + tg_x_xxyz_0[j];

                tg_xx_xzz_0[j] = -ab_x * tg_x_xzz_0[j] + tg_x_xxzz_0[j];

                tg_xx_yyy_0[j] = -ab_x * tg_x_yyy_0[j] + tg_x_xyyy_0[j];

                tg_xx_yyz_0[j] = -ab_x * tg_x_yyz_0[j] + tg_x_xyyz_0[j];

                tg_xx_yzz_0[j] = -ab_x * tg_x_yzz_0[j] + tg_x_xyzz_0[j];

                tg_xx_zzz_0[j] = -ab_x * tg_x_zzz_0[j] + tg_x_xzzz_0[j];

                tg_xy_xxx_0[j] = -ab_x * tg_y_xxx_0[j] + tg_y_xxxx_0[j];

                tg_xy_xxy_0[j] = -ab_x * tg_y_xxy_0[j] + tg_y_xxxy_0[j];

                tg_xy_xxz_0[j] = -ab_x * tg_y_xxz_0[j] + tg_y_xxxz_0[j];

                tg_xy_xyy_0[j] = -ab_x * tg_y_xyy_0[j] + tg_y_xxyy_0[j];

                tg_xy_xyz_0[j] = -ab_x * tg_y_xyz_0[j] + tg_y_xxyz_0[j];

                tg_xy_xzz_0[j] = -ab_x * tg_y_xzz_0[j] + tg_y_xxzz_0[j];

                tg_xy_yyy_0[j] = -ab_x * tg_y_yyy_0[j] + tg_y_xyyy_0[j];

                tg_xy_yyz_0[j] = -ab_x * tg_y_yyz_0[j] + tg_y_xyyz_0[j];

                tg_xy_yzz_0[j] = -ab_x * tg_y_yzz_0[j] + tg_y_xyzz_0[j];

                tg_xy_zzz_0[j] = -ab_x * tg_y_zzz_0[j] + tg_y_xzzz_0[j];

                tg_xz_xxx_0[j] = -ab_x * tg_z_xxx_0[j] + tg_z_xxxx_0[j];

                tg_xz_xxy_0[j] = -ab_x * tg_z_xxy_0[j] + tg_z_xxxy_0[j];

                tg_xz_xxz_0[j] = -ab_x * tg_z_xxz_0[j] + tg_z_xxxz_0[j];

                tg_xz_xyy_0[j] = -ab_x * tg_z_xyy_0[j] + tg_z_xxyy_0[j];

                tg_xz_xyz_0[j] = -ab_x * tg_z_xyz_0[j] + tg_z_xxyz_0[j];

                tg_xz_xzz_0[j] = -ab_x * tg_z_xzz_0[j] + tg_z_xxzz_0[j];

                tg_xz_yyy_0[j] = -ab_x * tg_z_yyy_0[j] + tg_z_xyyy_0[j];

                tg_xz_yyz_0[j] = -ab_x * tg_z_yyz_0[j] + tg_z_xyyz_0[j];

                tg_xz_yzz_0[j] = -ab_x * tg_z_yzz_0[j] + tg_z_xyzz_0[j];

                tg_xz_zzz_0[j] = -ab_x * tg_z_zzz_0[j] + tg_z_xzzz_0[j];

                tg_yy_xxx_0[j] = -ab_y * tg_y_xxx_0[j] + tg_y_xxxy_0[j];

                tg_yy_xxy_0[j] = -ab_y * tg_y_xxy_0[j] + tg_y_xxyy_0[j];

                tg_yy_xxz_0[j] = -ab_y * tg_y_xxz_0[j] + tg_y_xxyz_0[j];

                tg_yy_xyy_0[j] = -ab_y * tg_y_xyy_0[j] + tg_y_xyyy_0[j];

                tg_yy_xyz_0[j] = -ab_y * tg_y_xyz_0[j] + tg_y_xyyz_0[j];

                tg_yy_xzz_0[j] = -ab_y * tg_y_xzz_0[j] + tg_y_xyzz_0[j];

                tg_yy_yyy_0[j] = -ab_y * tg_y_yyy_0[j] + tg_y_yyyy_0[j];

                tg_yy_yyz_0[j] = -ab_y * tg_y_yyz_0[j] + tg_y_yyyz_0[j];

                tg_yy_yzz_0[j] = -ab_y * tg_y_yzz_0[j] + tg_y_yyzz_0[j];

                tg_yy_zzz_0[j] = -ab_y * tg_y_zzz_0[j] + tg_y_yzzz_0[j];

                tg_yz_xxx_0[j] = -ab_y * tg_z_xxx_0[j] + tg_z_xxxy_0[j];

                tg_yz_xxy_0[j] = -ab_y * tg_z_xxy_0[j] + tg_z_xxyy_0[j];

                tg_yz_xxz_0[j] = -ab_y * tg_z_xxz_0[j] + tg_z_xxyz_0[j];

                tg_yz_xyy_0[j] = -ab_y * tg_z_xyy_0[j] + tg_z_xyyy_0[j];

                tg_yz_xyz_0[j] = -ab_y * tg_z_xyz_0[j] + tg_z_xyyz_0[j];

                tg_yz_xzz_0[j] = -ab_y * tg_z_xzz_0[j] + tg_z_xyzz_0[j];

                tg_yz_yyy_0[j] = -ab_y * tg_z_yyy_0[j] + tg_z_yyyy_0[j];

                tg_yz_yyz_0[j] = -ab_y * tg_z_yyz_0[j] + tg_z_yyyz_0[j];

                tg_yz_yzz_0[j] = -ab_y * tg_z_yzz_0[j] + tg_z_yyzz_0[j];

                tg_yz_zzz_0[j] = -ab_y * tg_z_zzz_0[j] + tg_z_yzzz_0[j];

                tg_zz_xxx_0[j] = -ab_z * tg_z_xxx_0[j] + tg_z_xxxz_0[j];

                tg_zz_xxy_0[j] = -ab_z * tg_z_xxy_0[j] + tg_z_xxyz_0[j];

                tg_zz_xxz_0[j] = -ab_z * tg_z_xxz_0[j] + tg_z_xxzz_0[j];

                tg_zz_xyy_0[j] = -ab_z * tg_z_xyy_0[j] + tg_z_xyyz_0[j];

                tg_zz_xyz_0[j] = -ab_z * tg_z_xyz_0[j] + tg_z_xyzz_0[j];

                tg_zz_xzz_0[j] = -ab_z * tg_z_xzz_0[j] + tg_z_xzzz_0[j];

                tg_zz_yyy_0[j] = -ab_z * tg_z_yyy_0[j] + tg_z_yyyz_0[j];

                tg_zz_yyz_0[j] = -ab_z * tg_z_yyz_0[j] + tg_z_yyzz_0[j];

                tg_zz_yzz_0[j] = -ab_z * tg_z_yzz_0[j] + tg_z_yzzz_0[j];

                tg_zz_zzz_0[j] = -ab_z * tg_z_zzz_0[j] + tg_z_zzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForDGXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_1_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_x_xxxx_0 = braBuffer.data(pidx_g_1_4_m0 + i); 

            auto tg_x_xxxy_0 = braBuffer.data(pidx_g_1_4_m0 + kcomp + i); 

            auto tg_x_xxxz_0 = braBuffer.data(pidx_g_1_4_m0 + 2 * kcomp + i); 

            auto tg_x_xxyy_0 = braBuffer.data(pidx_g_1_4_m0 + 3 * kcomp + i); 

            auto tg_x_xxyz_0 = braBuffer.data(pidx_g_1_4_m0 + 4 * kcomp + i); 

            auto tg_x_xxzz_0 = braBuffer.data(pidx_g_1_4_m0 + 5 * kcomp + i); 

            auto tg_x_xyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 6 * kcomp + i); 

            auto tg_x_xyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 7 * kcomp + i); 

            auto tg_x_xyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 8 * kcomp + i); 

            auto tg_x_xzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 9 * kcomp + i); 

            auto tg_x_yyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 10 * kcomp + i); 

            auto tg_x_yyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 11 * kcomp + i); 

            auto tg_x_yyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 12 * kcomp + i); 

            auto tg_x_yzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 13 * kcomp + i); 

            auto tg_x_zzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 14 * kcomp + i); 

            auto tg_y_xxxx_0 = braBuffer.data(pidx_g_1_4_m0 + 15 * kcomp + i); 

            auto tg_y_xxxy_0 = braBuffer.data(pidx_g_1_4_m0 + 16 * kcomp + i); 

            auto tg_y_xxxz_0 = braBuffer.data(pidx_g_1_4_m0 + 17 * kcomp + i); 

            auto tg_y_xxyy_0 = braBuffer.data(pidx_g_1_4_m0 + 18 * kcomp + i); 

            auto tg_y_xxyz_0 = braBuffer.data(pidx_g_1_4_m0 + 19 * kcomp + i); 

            auto tg_y_xxzz_0 = braBuffer.data(pidx_g_1_4_m0 + 20 * kcomp + i); 

            auto tg_y_xyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 21 * kcomp + i); 

            auto tg_y_xyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 22 * kcomp + i); 

            auto tg_y_xyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 23 * kcomp + i); 

            auto tg_y_xzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 24 * kcomp + i); 

            auto tg_y_yyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 25 * kcomp + i); 

            auto tg_y_yyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 26 * kcomp + i); 

            auto tg_y_yyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 27 * kcomp + i); 

            auto tg_y_yzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 28 * kcomp + i); 

            auto tg_y_zzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 29 * kcomp + i); 

            auto tg_z_xxxx_0 = braBuffer.data(pidx_g_1_4_m0 + 30 * kcomp + i); 

            auto tg_z_xxxy_0 = braBuffer.data(pidx_g_1_4_m0 + 31 * kcomp + i); 

            auto tg_z_xxxz_0 = braBuffer.data(pidx_g_1_4_m0 + 32 * kcomp + i); 

            auto tg_z_xxyy_0 = braBuffer.data(pidx_g_1_4_m0 + 33 * kcomp + i); 

            auto tg_z_xxyz_0 = braBuffer.data(pidx_g_1_4_m0 + 34 * kcomp + i); 

            auto tg_z_xxzz_0 = braBuffer.data(pidx_g_1_4_m0 + 35 * kcomp + i); 

            auto tg_z_xyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 36 * kcomp + i); 

            auto tg_z_xyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 37 * kcomp + i); 

            auto tg_z_xyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 38 * kcomp + i); 

            auto tg_z_xzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 39 * kcomp + i); 

            auto tg_z_yyyy_0 = braBuffer.data(pidx_g_1_4_m0 + 40 * kcomp + i); 

            auto tg_z_yyyz_0 = braBuffer.data(pidx_g_1_4_m0 + 41 * kcomp + i); 

            auto tg_z_yyzz_0 = braBuffer.data(pidx_g_1_4_m0 + 42 * kcomp + i); 

            auto tg_z_yzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 43 * kcomp + i); 

            auto tg_z_zzzz_0 = braBuffer.data(pidx_g_1_4_m0 + 44 * kcomp + i); 

            auto tg_x_xxxxx_0 = braBuffer.data(pidx_g_1_5_m0 + i); 

            auto tg_x_xxxxy_0 = braBuffer.data(pidx_g_1_5_m0 + kcomp + i); 

            auto tg_x_xxxxz_0 = braBuffer.data(pidx_g_1_5_m0 + 2 * kcomp + i); 

            auto tg_x_xxxyy_0 = braBuffer.data(pidx_g_1_5_m0 + 3 * kcomp + i); 

            auto tg_x_xxxyz_0 = braBuffer.data(pidx_g_1_5_m0 + 4 * kcomp + i); 

            auto tg_x_xxxzz_0 = braBuffer.data(pidx_g_1_5_m0 + 5 * kcomp + i); 

            auto tg_x_xxyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 6 * kcomp + i); 

            auto tg_x_xxyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 7 * kcomp + i); 

            auto tg_x_xxyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 8 * kcomp + i); 

            auto tg_x_xxzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 9 * kcomp + i); 

            auto tg_x_xyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 10 * kcomp + i); 

            auto tg_x_xyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 11 * kcomp + i); 

            auto tg_x_xyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 12 * kcomp + i); 

            auto tg_x_xyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 13 * kcomp + i); 

            auto tg_x_xzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 14 * kcomp + i); 

            auto tg_y_xxxxx_0 = braBuffer.data(pidx_g_1_5_m0 + 21 * kcomp + i); 

            auto tg_y_xxxxy_0 = braBuffer.data(pidx_g_1_5_m0 + 22 * kcomp + i); 

            auto tg_y_xxxxz_0 = braBuffer.data(pidx_g_1_5_m0 + 23 * kcomp + i); 

            auto tg_y_xxxyy_0 = braBuffer.data(pidx_g_1_5_m0 + 24 * kcomp + i); 

            auto tg_y_xxxyz_0 = braBuffer.data(pidx_g_1_5_m0 + 25 * kcomp + i); 

            auto tg_y_xxxzz_0 = braBuffer.data(pidx_g_1_5_m0 + 26 * kcomp + i); 

            auto tg_y_xxyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 27 * kcomp + i); 

            auto tg_y_xxyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 28 * kcomp + i); 

            auto tg_y_xxyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 29 * kcomp + i); 

            auto tg_y_xxzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 30 * kcomp + i); 

            auto tg_y_xyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 31 * kcomp + i); 

            auto tg_y_xyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 32 * kcomp + i); 

            auto tg_y_xyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 33 * kcomp + i); 

            auto tg_y_xyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 34 * kcomp + i); 

            auto tg_y_xzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 35 * kcomp + i); 

            auto tg_y_yyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 36 * kcomp + i); 

            auto tg_y_yyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 37 * kcomp + i); 

            auto tg_y_yyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 38 * kcomp + i); 

            auto tg_y_yyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 39 * kcomp + i); 

            auto tg_y_yzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 40 * kcomp + i); 

            auto tg_z_xxxxx_0 = braBuffer.data(pidx_g_1_5_m0 + 42 * kcomp + i); 

            auto tg_z_xxxxy_0 = braBuffer.data(pidx_g_1_5_m0 + 43 * kcomp + i); 

            auto tg_z_xxxxz_0 = braBuffer.data(pidx_g_1_5_m0 + 44 * kcomp + i); 

            auto tg_z_xxxyy_0 = braBuffer.data(pidx_g_1_5_m0 + 45 * kcomp + i); 

            auto tg_z_xxxyz_0 = braBuffer.data(pidx_g_1_5_m0 + 46 * kcomp + i); 

            auto tg_z_xxxzz_0 = braBuffer.data(pidx_g_1_5_m0 + 47 * kcomp + i); 

            auto tg_z_xxyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 48 * kcomp + i); 

            auto tg_z_xxyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 49 * kcomp + i); 

            auto tg_z_xxyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 50 * kcomp + i); 

            auto tg_z_xxzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 51 * kcomp + i); 

            auto tg_z_xyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 52 * kcomp + i); 

            auto tg_z_xyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 53 * kcomp + i); 

            auto tg_z_xyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 54 * kcomp + i); 

            auto tg_z_xyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 55 * kcomp + i); 

            auto tg_z_xzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 56 * kcomp + i); 

            auto tg_z_yyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 57 * kcomp + i); 

            auto tg_z_yyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 58 * kcomp + i); 

            auto tg_z_yyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 59 * kcomp + i); 

            auto tg_z_yyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 60 * kcomp + i); 

            auto tg_z_yzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 61 * kcomp + i); 

            auto tg_z_zzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 62 * kcomp + i); 

            // set up pointers to integrals

            auto tg_xx_xxxx_0 = braBuffer.data(pidx_g_2_4_m0 + i); 

            auto tg_xx_xxxy_0 = braBuffer.data(pidx_g_2_4_m0 + kcomp + i); 

            auto tg_xx_xxxz_0 = braBuffer.data(pidx_g_2_4_m0 + 2 * kcomp + i); 

            auto tg_xx_xxyy_0 = braBuffer.data(pidx_g_2_4_m0 + 3 * kcomp + i); 

            auto tg_xx_xxyz_0 = braBuffer.data(pidx_g_2_4_m0 + 4 * kcomp + i); 

            auto tg_xx_xxzz_0 = braBuffer.data(pidx_g_2_4_m0 + 5 * kcomp + i); 

            auto tg_xx_xyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 6 * kcomp + i); 

            auto tg_xx_xyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 7 * kcomp + i); 

            auto tg_xx_xyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 8 * kcomp + i); 

            auto tg_xx_xzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 9 * kcomp + i); 

            auto tg_xx_yyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 10 * kcomp + i); 

            auto tg_xx_yyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 11 * kcomp + i); 

            auto tg_xx_yyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 12 * kcomp + i); 

            auto tg_xx_yzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 13 * kcomp + i); 

            auto tg_xx_zzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 14 * kcomp + i); 

            auto tg_xy_xxxx_0 = braBuffer.data(pidx_g_2_4_m0 + 15 * kcomp + i); 

            auto tg_xy_xxxy_0 = braBuffer.data(pidx_g_2_4_m0 + 16 * kcomp + i); 

            auto tg_xy_xxxz_0 = braBuffer.data(pidx_g_2_4_m0 + 17 * kcomp + i); 

            auto tg_xy_xxyy_0 = braBuffer.data(pidx_g_2_4_m0 + 18 * kcomp + i); 

            auto tg_xy_xxyz_0 = braBuffer.data(pidx_g_2_4_m0 + 19 * kcomp + i); 

            auto tg_xy_xxzz_0 = braBuffer.data(pidx_g_2_4_m0 + 20 * kcomp + i); 

            auto tg_xy_xyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 21 * kcomp + i); 

            auto tg_xy_xyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 22 * kcomp + i); 

            auto tg_xy_xyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 23 * kcomp + i); 

            auto tg_xy_xzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 24 * kcomp + i); 

            auto tg_xy_yyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 25 * kcomp + i); 

            auto tg_xy_yyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 26 * kcomp + i); 

            auto tg_xy_yyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 27 * kcomp + i); 

            auto tg_xy_yzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 28 * kcomp + i); 

            auto tg_xy_zzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 29 * kcomp + i); 

            auto tg_xz_xxxx_0 = braBuffer.data(pidx_g_2_4_m0 + 30 * kcomp + i); 

            auto tg_xz_xxxy_0 = braBuffer.data(pidx_g_2_4_m0 + 31 * kcomp + i); 

            auto tg_xz_xxxz_0 = braBuffer.data(pidx_g_2_4_m0 + 32 * kcomp + i); 

            auto tg_xz_xxyy_0 = braBuffer.data(pidx_g_2_4_m0 + 33 * kcomp + i); 

            auto tg_xz_xxyz_0 = braBuffer.data(pidx_g_2_4_m0 + 34 * kcomp + i); 

            auto tg_xz_xxzz_0 = braBuffer.data(pidx_g_2_4_m0 + 35 * kcomp + i); 

            auto tg_xz_xyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 36 * kcomp + i); 

            auto tg_xz_xyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 37 * kcomp + i); 

            auto tg_xz_xyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 38 * kcomp + i); 

            auto tg_xz_xzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 39 * kcomp + i); 

            auto tg_xz_yyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 40 * kcomp + i); 

            auto tg_xz_yyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 41 * kcomp + i); 

            auto tg_xz_yyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 42 * kcomp + i); 

            auto tg_xz_yzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 43 * kcomp + i); 

            auto tg_xz_zzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 44 * kcomp + i); 

            auto tg_yy_xxxx_0 = braBuffer.data(pidx_g_2_4_m0 + 45 * kcomp + i); 

            auto tg_yy_xxxy_0 = braBuffer.data(pidx_g_2_4_m0 + 46 * kcomp + i); 

            auto tg_yy_xxxz_0 = braBuffer.data(pidx_g_2_4_m0 + 47 * kcomp + i); 

            auto tg_yy_xxyy_0 = braBuffer.data(pidx_g_2_4_m0 + 48 * kcomp + i); 

            auto tg_yy_xxyz_0 = braBuffer.data(pidx_g_2_4_m0 + 49 * kcomp + i); 

            auto tg_yy_xxzz_0 = braBuffer.data(pidx_g_2_4_m0 + 50 * kcomp + i); 

            auto tg_yy_xyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 51 * kcomp + i); 

            auto tg_yy_xyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 52 * kcomp + i); 

            auto tg_yy_xyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 53 * kcomp + i); 

            auto tg_yy_xzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 54 * kcomp + i); 

            auto tg_yy_yyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 55 * kcomp + i); 

            auto tg_yy_yyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 56 * kcomp + i); 

            auto tg_yy_yyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 57 * kcomp + i); 

            auto tg_yy_yzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 58 * kcomp + i); 

            auto tg_yy_zzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 59 * kcomp + i); 

            auto tg_yz_xxxx_0 = braBuffer.data(pidx_g_2_4_m0 + 60 * kcomp + i); 

            auto tg_yz_xxxy_0 = braBuffer.data(pidx_g_2_4_m0 + 61 * kcomp + i); 

            auto tg_yz_xxxz_0 = braBuffer.data(pidx_g_2_4_m0 + 62 * kcomp + i); 

            auto tg_yz_xxyy_0 = braBuffer.data(pidx_g_2_4_m0 + 63 * kcomp + i); 

            auto tg_yz_xxyz_0 = braBuffer.data(pidx_g_2_4_m0 + 64 * kcomp + i); 

            auto tg_yz_xxzz_0 = braBuffer.data(pidx_g_2_4_m0 + 65 * kcomp + i); 

            auto tg_yz_xyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 66 * kcomp + i); 

            auto tg_yz_xyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 67 * kcomp + i); 

            auto tg_yz_xyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 68 * kcomp + i); 

            auto tg_yz_xzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 69 * kcomp + i); 

            auto tg_yz_yyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 70 * kcomp + i); 

            auto tg_yz_yyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 71 * kcomp + i); 

            auto tg_yz_yyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 72 * kcomp + i); 

            auto tg_yz_yzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 73 * kcomp + i); 

            auto tg_yz_zzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 74 * kcomp + i); 

            auto tg_zz_xxxx_0 = braBuffer.data(pidx_g_2_4_m0 + 75 * kcomp + i); 

            auto tg_zz_xxxy_0 = braBuffer.data(pidx_g_2_4_m0 + 76 * kcomp + i); 

            auto tg_zz_xxxz_0 = braBuffer.data(pidx_g_2_4_m0 + 77 * kcomp + i); 

            auto tg_zz_xxyy_0 = braBuffer.data(pidx_g_2_4_m0 + 78 * kcomp + i); 

            auto tg_zz_xxyz_0 = braBuffer.data(pidx_g_2_4_m0 + 79 * kcomp + i); 

            auto tg_zz_xxzz_0 = braBuffer.data(pidx_g_2_4_m0 + 80 * kcomp + i); 

            auto tg_zz_xyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 81 * kcomp + i); 

            auto tg_zz_xyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 82 * kcomp + i); 

            auto tg_zz_xyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 83 * kcomp + i); 

            auto tg_zz_xzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 84 * kcomp + i); 

            auto tg_zz_yyyy_0 = braBuffer.data(pidx_g_2_4_m0 + 85 * kcomp + i); 

            auto tg_zz_yyyz_0 = braBuffer.data(pidx_g_2_4_m0 + 86 * kcomp + i); 

            auto tg_zz_yyzz_0 = braBuffer.data(pidx_g_2_4_m0 + 87 * kcomp + i); 

            auto tg_zz_yzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 88 * kcomp + i); 

            auto tg_zz_zzzz_0 = braBuffer.data(pidx_g_2_4_m0 + 89 * kcomp + i); 

            #pragma omp simd aligned(tg_x_xxxx_0, tg_x_xxxxx_0, tg_x_xxxxy_0, tg_x_xxxxz_0, tg_x_xxxy_0, \
                                     tg_x_xxxyy_0, tg_x_xxxyz_0, tg_x_xxxz_0, tg_x_xxxzz_0, tg_x_xxyy_0, tg_x_xxyyy_0, \
                                     tg_x_xxyyz_0, tg_x_xxyz_0, tg_x_xxyzz_0, tg_x_xxzz_0, tg_x_xxzzz_0, tg_x_xyyy_0, \
                                     tg_x_xyyyy_0, tg_x_xyyyz_0, tg_x_xyyz_0, tg_x_xyyzz_0, tg_x_xyzz_0, tg_x_xyzzz_0, \
                                     tg_x_xzzz_0, tg_x_xzzzz_0, tg_x_yyyy_0, tg_x_yyyz_0, tg_x_yyzz_0, tg_x_yzzz_0, \
                                     tg_x_zzzz_0, tg_xx_xxxx_0, tg_xx_xxxy_0, tg_xx_xxxz_0, tg_xx_xxyy_0, tg_xx_xxyz_0, \
                                     tg_xx_xxzz_0, tg_xx_xyyy_0, tg_xx_xyyz_0, tg_xx_xyzz_0, tg_xx_xzzz_0, tg_xx_yyyy_0, \
                                     tg_xx_yyyz_0, tg_xx_yyzz_0, tg_xx_yzzz_0, tg_xx_zzzz_0, tg_xy_xxxx_0, tg_xy_xxxy_0, \
                                     tg_xy_xxxz_0, tg_xy_xxyy_0, tg_xy_xxyz_0, tg_xy_xxzz_0, tg_xy_xyyy_0, tg_xy_xyyz_0, \
                                     tg_xy_xyzz_0, tg_xy_xzzz_0, tg_xy_yyyy_0, tg_xy_yyyz_0, tg_xy_yyzz_0, tg_xy_yzzz_0, \
                                     tg_xy_zzzz_0, tg_xz_xxxx_0, tg_xz_xxxy_0, tg_xz_xxxz_0, tg_xz_xxyy_0, tg_xz_xxyz_0, \
                                     tg_xz_xxzz_0, tg_xz_xyyy_0, tg_xz_xyyz_0, tg_xz_xyzz_0, tg_xz_xzzz_0, tg_xz_yyyy_0, \
                                     tg_xz_yyyz_0, tg_xz_yyzz_0, tg_xz_yzzz_0, tg_xz_zzzz_0, tg_y_xxxx_0, tg_y_xxxxx_0, \
                                     tg_y_xxxxy_0, tg_y_xxxxz_0, tg_y_xxxy_0, tg_y_xxxyy_0, tg_y_xxxyz_0, tg_y_xxxz_0, \
                                     tg_y_xxxzz_0, tg_y_xxyy_0, tg_y_xxyyy_0, tg_y_xxyyz_0, tg_y_xxyz_0, tg_y_xxyzz_0, \
                                     tg_y_xxzz_0, tg_y_xxzzz_0, tg_y_xyyy_0, tg_y_xyyyy_0, tg_y_xyyyz_0, tg_y_xyyz_0, \
                                     tg_y_xyyzz_0, tg_y_xyzz_0, tg_y_xyzzz_0, tg_y_xzzz_0, tg_y_xzzzz_0, tg_y_yyyy_0, \
                                     tg_y_yyyyy_0, tg_y_yyyyz_0, tg_y_yyyz_0, tg_y_yyyzz_0, tg_y_yyzz_0, tg_y_yyzzz_0, \
                                     tg_y_yzzz_0, tg_y_yzzzz_0, tg_y_zzzz_0, tg_yy_xxxx_0, tg_yy_xxxy_0, tg_yy_xxxz_0, \
                                     tg_yy_xxyy_0, tg_yy_xxyz_0, tg_yy_xxzz_0, tg_yy_xyyy_0, tg_yy_xyyz_0, tg_yy_xyzz_0, \
                                     tg_yy_xzzz_0, tg_yy_yyyy_0, tg_yy_yyyz_0, tg_yy_yyzz_0, tg_yy_yzzz_0, tg_yy_zzzz_0, \
                                     tg_yz_xxxx_0, tg_yz_xxxy_0, tg_yz_xxxz_0, tg_yz_xxyy_0, tg_yz_xxyz_0, tg_yz_xxzz_0, \
                                     tg_yz_xyyy_0, tg_yz_xyyz_0, tg_yz_xyzz_0, tg_yz_xzzz_0, tg_yz_yyyy_0, tg_yz_yyyz_0, \
                                     tg_yz_yyzz_0, tg_yz_yzzz_0, tg_yz_zzzz_0, tg_z_xxxx_0, tg_z_xxxxx_0, tg_z_xxxxy_0, \
                                     tg_z_xxxxz_0, tg_z_xxxy_0, tg_z_xxxyy_0, tg_z_xxxyz_0, tg_z_xxxz_0, tg_z_xxxzz_0, \
                                     tg_z_xxyy_0, tg_z_xxyyy_0, tg_z_xxyyz_0, tg_z_xxyz_0, tg_z_xxyzz_0, tg_z_xxzz_0, \
                                     tg_z_xxzzz_0, tg_z_xyyy_0, tg_z_xyyyy_0, tg_z_xyyyz_0, tg_z_xyyz_0, tg_z_xyyzz_0, \
                                     tg_z_xyzz_0, tg_z_xyzzz_0, tg_z_xzzz_0, tg_z_xzzzz_0, tg_z_yyyy_0, tg_z_yyyyy_0, \
                                     tg_z_yyyyz_0, tg_z_yyyz_0, tg_z_yyyzz_0, tg_z_yyzz_0, tg_z_yyzzz_0, tg_z_yzzz_0, \
                                     tg_z_yzzzz_0, tg_z_zzzz_0, tg_z_zzzzz_0, tg_zz_xxxx_0, tg_zz_xxxy_0, tg_zz_xxxz_0, \
                                     tg_zz_xxyy_0, tg_zz_xxyz_0, tg_zz_xxzz_0, tg_zz_xyyy_0, tg_zz_xyyz_0, tg_zz_xyzz_0, \
                                     tg_zz_xzzz_0, tg_zz_yyyy_0, tg_zz_yyyz_0, tg_zz_yyzz_0, tg_zz_yzzz_0, tg_zz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xx_xxxx_0[j] = -ab_x * tg_x_xxxx_0[j] + tg_x_xxxxx_0[j];

                tg_xx_xxxy_0[j] = -ab_x * tg_x_xxxy_0[j] + tg_x_xxxxy_0[j];

                tg_xx_xxxz_0[j] = -ab_x * tg_x_xxxz_0[j] + tg_x_xxxxz_0[j];

                tg_xx_xxyy_0[j] = -ab_x * tg_x_xxyy_0[j] + tg_x_xxxyy_0[j];

                tg_xx_xxyz_0[j] = -ab_x * tg_x_xxyz_0[j] + tg_x_xxxyz_0[j];

                tg_xx_xxzz_0[j] = -ab_x * tg_x_xxzz_0[j] + tg_x_xxxzz_0[j];

                tg_xx_xyyy_0[j] = -ab_x * tg_x_xyyy_0[j] + tg_x_xxyyy_0[j];

                tg_xx_xyyz_0[j] = -ab_x * tg_x_xyyz_0[j] + tg_x_xxyyz_0[j];

                tg_xx_xyzz_0[j] = -ab_x * tg_x_xyzz_0[j] + tg_x_xxyzz_0[j];

                tg_xx_xzzz_0[j] = -ab_x * tg_x_xzzz_0[j] + tg_x_xxzzz_0[j];

                tg_xx_yyyy_0[j] = -ab_x * tg_x_yyyy_0[j] + tg_x_xyyyy_0[j];

                tg_xx_yyyz_0[j] = -ab_x * tg_x_yyyz_0[j] + tg_x_xyyyz_0[j];

                tg_xx_yyzz_0[j] = -ab_x * tg_x_yyzz_0[j] + tg_x_xyyzz_0[j];

                tg_xx_yzzz_0[j] = -ab_x * tg_x_yzzz_0[j] + tg_x_xyzzz_0[j];

                tg_xx_zzzz_0[j] = -ab_x * tg_x_zzzz_0[j] + tg_x_xzzzz_0[j];

                tg_xy_xxxx_0[j] = -ab_x * tg_y_xxxx_0[j] + tg_y_xxxxx_0[j];

                tg_xy_xxxy_0[j] = -ab_x * tg_y_xxxy_0[j] + tg_y_xxxxy_0[j];

                tg_xy_xxxz_0[j] = -ab_x * tg_y_xxxz_0[j] + tg_y_xxxxz_0[j];

                tg_xy_xxyy_0[j] = -ab_x * tg_y_xxyy_0[j] + tg_y_xxxyy_0[j];

                tg_xy_xxyz_0[j] = -ab_x * tg_y_xxyz_0[j] + tg_y_xxxyz_0[j];

                tg_xy_xxzz_0[j] = -ab_x * tg_y_xxzz_0[j] + tg_y_xxxzz_0[j];

                tg_xy_xyyy_0[j] = -ab_x * tg_y_xyyy_0[j] + tg_y_xxyyy_0[j];

                tg_xy_xyyz_0[j] = -ab_x * tg_y_xyyz_0[j] + tg_y_xxyyz_0[j];

                tg_xy_xyzz_0[j] = -ab_x * tg_y_xyzz_0[j] + tg_y_xxyzz_0[j];

                tg_xy_xzzz_0[j] = -ab_x * tg_y_xzzz_0[j] + tg_y_xxzzz_0[j];

                tg_xy_yyyy_0[j] = -ab_x * tg_y_yyyy_0[j] + tg_y_xyyyy_0[j];

                tg_xy_yyyz_0[j] = -ab_x * tg_y_yyyz_0[j] + tg_y_xyyyz_0[j];

                tg_xy_yyzz_0[j] = -ab_x * tg_y_yyzz_0[j] + tg_y_xyyzz_0[j];

                tg_xy_yzzz_0[j] = -ab_x * tg_y_yzzz_0[j] + tg_y_xyzzz_0[j];

                tg_xy_zzzz_0[j] = -ab_x * tg_y_zzzz_0[j] + tg_y_xzzzz_0[j];

                tg_xz_xxxx_0[j] = -ab_x * tg_z_xxxx_0[j] + tg_z_xxxxx_0[j];

                tg_xz_xxxy_0[j] = -ab_x * tg_z_xxxy_0[j] + tg_z_xxxxy_0[j];

                tg_xz_xxxz_0[j] = -ab_x * tg_z_xxxz_0[j] + tg_z_xxxxz_0[j];

                tg_xz_xxyy_0[j] = -ab_x * tg_z_xxyy_0[j] + tg_z_xxxyy_0[j];

                tg_xz_xxyz_0[j] = -ab_x * tg_z_xxyz_0[j] + tg_z_xxxyz_0[j];

                tg_xz_xxzz_0[j] = -ab_x * tg_z_xxzz_0[j] + tg_z_xxxzz_0[j];

                tg_xz_xyyy_0[j] = -ab_x * tg_z_xyyy_0[j] + tg_z_xxyyy_0[j];

                tg_xz_xyyz_0[j] = -ab_x * tg_z_xyyz_0[j] + tg_z_xxyyz_0[j];

                tg_xz_xyzz_0[j] = -ab_x * tg_z_xyzz_0[j] + tg_z_xxyzz_0[j];

                tg_xz_xzzz_0[j] = -ab_x * tg_z_xzzz_0[j] + tg_z_xxzzz_0[j];

                tg_xz_yyyy_0[j] = -ab_x * tg_z_yyyy_0[j] + tg_z_xyyyy_0[j];

                tg_xz_yyyz_0[j] = -ab_x * tg_z_yyyz_0[j] + tg_z_xyyyz_0[j];

                tg_xz_yyzz_0[j] = -ab_x * tg_z_yyzz_0[j] + tg_z_xyyzz_0[j];

                tg_xz_yzzz_0[j] = -ab_x * tg_z_yzzz_0[j] + tg_z_xyzzz_0[j];

                tg_xz_zzzz_0[j] = -ab_x * tg_z_zzzz_0[j] + tg_z_xzzzz_0[j];

                tg_yy_xxxx_0[j] = -ab_y * tg_y_xxxx_0[j] + tg_y_xxxxy_0[j];

                tg_yy_xxxy_0[j] = -ab_y * tg_y_xxxy_0[j] + tg_y_xxxyy_0[j];

                tg_yy_xxxz_0[j] = -ab_y * tg_y_xxxz_0[j] + tg_y_xxxyz_0[j];

                tg_yy_xxyy_0[j] = -ab_y * tg_y_xxyy_0[j] + tg_y_xxyyy_0[j];

                tg_yy_xxyz_0[j] = -ab_y * tg_y_xxyz_0[j] + tg_y_xxyyz_0[j];

                tg_yy_xxzz_0[j] = -ab_y * tg_y_xxzz_0[j] + tg_y_xxyzz_0[j];

                tg_yy_xyyy_0[j] = -ab_y * tg_y_xyyy_0[j] + tg_y_xyyyy_0[j];

                tg_yy_xyyz_0[j] = -ab_y * tg_y_xyyz_0[j] + tg_y_xyyyz_0[j];

                tg_yy_xyzz_0[j] = -ab_y * tg_y_xyzz_0[j] + tg_y_xyyzz_0[j];

                tg_yy_xzzz_0[j] = -ab_y * tg_y_xzzz_0[j] + tg_y_xyzzz_0[j];

                tg_yy_yyyy_0[j] = -ab_y * tg_y_yyyy_0[j] + tg_y_yyyyy_0[j];

                tg_yy_yyyz_0[j] = -ab_y * tg_y_yyyz_0[j] + tg_y_yyyyz_0[j];

                tg_yy_yyzz_0[j] = -ab_y * tg_y_yyzz_0[j] + tg_y_yyyzz_0[j];

                tg_yy_yzzz_0[j] = -ab_y * tg_y_yzzz_0[j] + tg_y_yyzzz_0[j];

                tg_yy_zzzz_0[j] = -ab_y * tg_y_zzzz_0[j] + tg_y_yzzzz_0[j];

                tg_yz_xxxx_0[j] = -ab_y * tg_z_xxxx_0[j] + tg_z_xxxxy_0[j];

                tg_yz_xxxy_0[j] = -ab_y * tg_z_xxxy_0[j] + tg_z_xxxyy_0[j];

                tg_yz_xxxz_0[j] = -ab_y * tg_z_xxxz_0[j] + tg_z_xxxyz_0[j];

                tg_yz_xxyy_0[j] = -ab_y * tg_z_xxyy_0[j] + tg_z_xxyyy_0[j];

                tg_yz_xxyz_0[j] = -ab_y * tg_z_xxyz_0[j] + tg_z_xxyyz_0[j];

                tg_yz_xxzz_0[j] = -ab_y * tg_z_xxzz_0[j] + tg_z_xxyzz_0[j];

                tg_yz_xyyy_0[j] = -ab_y * tg_z_xyyy_0[j] + tg_z_xyyyy_0[j];

                tg_yz_xyyz_0[j] = -ab_y * tg_z_xyyz_0[j] + tg_z_xyyyz_0[j];

                tg_yz_xyzz_0[j] = -ab_y * tg_z_xyzz_0[j] + tg_z_xyyzz_0[j];

                tg_yz_xzzz_0[j] = -ab_y * tg_z_xzzz_0[j] + tg_z_xyzzz_0[j];

                tg_yz_yyyy_0[j] = -ab_y * tg_z_yyyy_0[j] + tg_z_yyyyy_0[j];

                tg_yz_yyyz_0[j] = -ab_y * tg_z_yyyz_0[j] + tg_z_yyyyz_0[j];

                tg_yz_yyzz_0[j] = -ab_y * tg_z_yyzz_0[j] + tg_z_yyyzz_0[j];

                tg_yz_yzzz_0[j] = -ab_y * tg_z_yzzz_0[j] + tg_z_yyzzz_0[j];

                tg_yz_zzzz_0[j] = -ab_y * tg_z_zzzz_0[j] + tg_z_yzzzz_0[j];

                tg_zz_xxxx_0[j] = -ab_z * tg_z_xxxx_0[j] + tg_z_xxxxz_0[j];

                tg_zz_xxxy_0[j] = -ab_z * tg_z_xxxy_0[j] + tg_z_xxxyz_0[j];

                tg_zz_xxxz_0[j] = -ab_z * tg_z_xxxz_0[j] + tg_z_xxxzz_0[j];

                tg_zz_xxyy_0[j] = -ab_z * tg_z_xxyy_0[j] + tg_z_xxyyz_0[j];

                tg_zz_xxyz_0[j] = -ab_z * tg_z_xxyz_0[j] + tg_z_xxyzz_0[j];

                tg_zz_xxzz_0[j] = -ab_z * tg_z_xxzz_0[j] + tg_z_xxzzz_0[j];

                tg_zz_xyyy_0[j] = -ab_z * tg_z_xyyy_0[j] + tg_z_xyyyz_0[j];

                tg_zz_xyyz_0[j] = -ab_z * tg_z_xyyz_0[j] + tg_z_xyyzz_0[j];

                tg_zz_xyzz_0[j] = -ab_z * tg_z_xyzz_0[j] + tg_z_xyzzz_0[j];

                tg_zz_xzzz_0[j] = -ab_z * tg_z_xzzz_0[j] + tg_z_xzzzz_0[j];

                tg_zz_yyyy_0[j] = -ab_z * tg_z_yyyy_0[j] + tg_z_yyyyz_0[j];

                tg_zz_yyyz_0[j] = -ab_z * tg_z_yyyz_0[j] + tg_z_yyyzz_0[j];

                tg_zz_yyzz_0[j] = -ab_z * tg_z_yyzz_0[j] + tg_z_yyzzz_0[j];

                tg_zz_yzzz_0[j] = -ab_z * tg_z_yzzz_0[j] + tg_z_yzzzz_0[j];

                tg_zz_zzzz_0[j] = -ab_z * tg_z_zzzz_0[j] + tg_z_zzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForDHXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        eribrrfunc::compElectronRepulsionForDHXY_0_63(braBuffer,
                                                      recursionMap,
                                                      abDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        eribrrfunc::compElectronRepulsionForDHXY_63_126(braBuffer,
                                                        recursionMap,
                                                        abDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForDHXY_0_63(      CMemBlock2D<double>& braBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& abDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,63)

        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_2_5_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_1_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_1_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_x_xxxxx_0 = braBuffer.data(pidx_g_1_5_m0 + i); 

            auto tg_x_xxxxy_0 = braBuffer.data(pidx_g_1_5_m0 + kcomp + i); 

            auto tg_x_xxxxz_0 = braBuffer.data(pidx_g_1_5_m0 + 2 * kcomp + i); 

            auto tg_x_xxxyy_0 = braBuffer.data(pidx_g_1_5_m0 + 3 * kcomp + i); 

            auto tg_x_xxxyz_0 = braBuffer.data(pidx_g_1_5_m0 + 4 * kcomp + i); 

            auto tg_x_xxxzz_0 = braBuffer.data(pidx_g_1_5_m0 + 5 * kcomp + i); 

            auto tg_x_xxyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 6 * kcomp + i); 

            auto tg_x_xxyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 7 * kcomp + i); 

            auto tg_x_xxyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 8 * kcomp + i); 

            auto tg_x_xxzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 9 * kcomp + i); 

            auto tg_x_xyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 10 * kcomp + i); 

            auto tg_x_xyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 11 * kcomp + i); 

            auto tg_x_xyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 12 * kcomp + i); 

            auto tg_x_xyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 13 * kcomp + i); 

            auto tg_x_xzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 14 * kcomp + i); 

            auto tg_x_yyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 15 * kcomp + i); 

            auto tg_x_yyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 16 * kcomp + i); 

            auto tg_x_yyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 17 * kcomp + i); 

            auto tg_x_yyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 18 * kcomp + i); 

            auto tg_x_yzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 19 * kcomp + i); 

            auto tg_x_zzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 20 * kcomp + i); 

            auto tg_y_xxxxx_0 = braBuffer.data(pidx_g_1_5_m0 + 21 * kcomp + i); 

            auto tg_y_xxxxy_0 = braBuffer.data(pidx_g_1_5_m0 + 22 * kcomp + i); 

            auto tg_y_xxxxz_0 = braBuffer.data(pidx_g_1_5_m0 + 23 * kcomp + i); 

            auto tg_y_xxxyy_0 = braBuffer.data(pidx_g_1_5_m0 + 24 * kcomp + i); 

            auto tg_y_xxxyz_0 = braBuffer.data(pidx_g_1_5_m0 + 25 * kcomp + i); 

            auto tg_y_xxxzz_0 = braBuffer.data(pidx_g_1_5_m0 + 26 * kcomp + i); 

            auto tg_y_xxyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 27 * kcomp + i); 

            auto tg_y_xxyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 28 * kcomp + i); 

            auto tg_y_xxyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 29 * kcomp + i); 

            auto tg_y_xxzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 30 * kcomp + i); 

            auto tg_y_xyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 31 * kcomp + i); 

            auto tg_y_xyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 32 * kcomp + i); 

            auto tg_y_xyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 33 * kcomp + i); 

            auto tg_y_xyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 34 * kcomp + i); 

            auto tg_y_xzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 35 * kcomp + i); 

            auto tg_y_yyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 36 * kcomp + i); 

            auto tg_y_yyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 37 * kcomp + i); 

            auto tg_y_yyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 38 * kcomp + i); 

            auto tg_y_yyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 39 * kcomp + i); 

            auto tg_y_yzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 40 * kcomp + i); 

            auto tg_y_zzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 41 * kcomp + i); 

            auto tg_z_xxxxx_0 = braBuffer.data(pidx_g_1_5_m0 + 42 * kcomp + i); 

            auto tg_z_xxxxy_0 = braBuffer.data(pidx_g_1_5_m0 + 43 * kcomp + i); 

            auto tg_z_xxxxz_0 = braBuffer.data(pidx_g_1_5_m0 + 44 * kcomp + i); 

            auto tg_z_xxxyy_0 = braBuffer.data(pidx_g_1_5_m0 + 45 * kcomp + i); 

            auto tg_z_xxxyz_0 = braBuffer.data(pidx_g_1_5_m0 + 46 * kcomp + i); 

            auto tg_z_xxxzz_0 = braBuffer.data(pidx_g_1_5_m0 + 47 * kcomp + i); 

            auto tg_z_xxyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 48 * kcomp + i); 

            auto tg_z_xxyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 49 * kcomp + i); 

            auto tg_z_xxyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 50 * kcomp + i); 

            auto tg_z_xxzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 51 * kcomp + i); 

            auto tg_z_xyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 52 * kcomp + i); 

            auto tg_z_xyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 53 * kcomp + i); 

            auto tg_z_xyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 54 * kcomp + i); 

            auto tg_z_xyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 55 * kcomp + i); 

            auto tg_z_xzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 56 * kcomp + i); 

            auto tg_z_yyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 57 * kcomp + i); 

            auto tg_z_yyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 58 * kcomp + i); 

            auto tg_z_yyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 59 * kcomp + i); 

            auto tg_z_yyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 60 * kcomp + i); 

            auto tg_z_yzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 61 * kcomp + i); 

            auto tg_z_zzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 62 * kcomp + i); 

            auto tg_x_xxxxxx_0 = braBuffer.data(pidx_g_1_6_m0 + i); 

            auto tg_x_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + kcomp + i); 

            auto tg_x_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 2 * kcomp + i); 

            auto tg_x_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 3 * kcomp + i); 

            auto tg_x_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 4 * kcomp + i); 

            auto tg_x_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 5 * kcomp + i); 

            auto tg_x_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 6 * kcomp + i); 

            auto tg_x_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 7 * kcomp + i); 

            auto tg_x_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 8 * kcomp + i); 

            auto tg_x_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 9 * kcomp + i); 

            auto tg_x_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 10 * kcomp + i); 

            auto tg_x_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 11 * kcomp + i); 

            auto tg_x_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 12 * kcomp + i); 

            auto tg_x_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 13 * kcomp + i); 

            auto tg_x_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 14 * kcomp + i); 

            auto tg_x_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 15 * kcomp + i); 

            auto tg_x_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 16 * kcomp + i); 

            auto tg_x_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 17 * kcomp + i); 

            auto tg_x_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 18 * kcomp + i); 

            auto tg_x_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 19 * kcomp + i); 

            auto tg_x_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 20 * kcomp + i); 

            auto tg_y_xxxxxx_0 = braBuffer.data(pidx_g_1_6_m0 + 28 * kcomp + i); 

            auto tg_y_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + 29 * kcomp + i); 

            auto tg_y_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 30 * kcomp + i); 

            auto tg_y_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 31 * kcomp + i); 

            auto tg_y_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 32 * kcomp + i); 

            auto tg_y_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 33 * kcomp + i); 

            auto tg_y_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 34 * kcomp + i); 

            auto tg_y_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 35 * kcomp + i); 

            auto tg_y_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 36 * kcomp + i); 

            auto tg_y_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 37 * kcomp + i); 

            auto tg_y_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 38 * kcomp + i); 

            auto tg_y_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 39 * kcomp + i); 

            auto tg_y_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 40 * kcomp + i); 

            auto tg_y_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 41 * kcomp + i); 

            auto tg_y_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 42 * kcomp + i); 

            auto tg_y_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 43 * kcomp + i); 

            auto tg_y_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 44 * kcomp + i); 

            auto tg_y_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 45 * kcomp + i); 

            auto tg_y_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 46 * kcomp + i); 

            auto tg_y_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 47 * kcomp + i); 

            auto tg_y_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 48 * kcomp + i); 

            auto tg_z_xxxxxx_0 = braBuffer.data(pidx_g_1_6_m0 + 56 * kcomp + i); 

            auto tg_z_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + 57 * kcomp + i); 

            auto tg_z_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 58 * kcomp + i); 

            auto tg_z_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 59 * kcomp + i); 

            auto tg_z_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 60 * kcomp + i); 

            auto tg_z_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 61 * kcomp + i); 

            auto tg_z_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 62 * kcomp + i); 

            auto tg_z_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 63 * kcomp + i); 

            auto tg_z_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 64 * kcomp + i); 

            auto tg_z_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 65 * kcomp + i); 

            auto tg_z_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 66 * kcomp + i); 

            auto tg_z_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 67 * kcomp + i); 

            auto tg_z_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 68 * kcomp + i); 

            auto tg_z_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 69 * kcomp + i); 

            auto tg_z_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 70 * kcomp + i); 

            auto tg_z_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 71 * kcomp + i); 

            auto tg_z_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 72 * kcomp + i); 

            auto tg_z_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 73 * kcomp + i); 

            auto tg_z_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 74 * kcomp + i); 

            auto tg_z_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 75 * kcomp + i); 

            auto tg_z_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 76 * kcomp + i); 

            // set up pointers to integrals

            auto tg_xx_xxxxx_0 = braBuffer.data(pidx_g_2_5_m0 + i); 

            auto tg_xx_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + kcomp + i); 

            auto tg_xx_xxxxz_0 = braBuffer.data(pidx_g_2_5_m0 + 2 * kcomp + i); 

            auto tg_xx_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 3 * kcomp + i); 

            auto tg_xx_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 4 * kcomp + i); 

            auto tg_xx_xxxzz_0 = braBuffer.data(pidx_g_2_5_m0 + 5 * kcomp + i); 

            auto tg_xx_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 6 * kcomp + i); 

            auto tg_xx_xxyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 7 * kcomp + i); 

            auto tg_xx_xxyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 8 * kcomp + i); 

            auto tg_xx_xxzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 9 * kcomp + i); 

            auto tg_xx_xyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 10 * kcomp + i); 

            auto tg_xx_xyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 11 * kcomp + i); 

            auto tg_xx_xyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 12 * kcomp + i); 

            auto tg_xx_xyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 13 * kcomp + i); 

            auto tg_xx_xzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 14 * kcomp + i); 

            auto tg_xx_yyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 15 * kcomp + i); 

            auto tg_xx_yyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 16 * kcomp + i); 

            auto tg_xx_yyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 17 * kcomp + i); 

            auto tg_xx_yyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 18 * kcomp + i); 

            auto tg_xx_yzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 19 * kcomp + i); 

            auto tg_xx_zzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 20 * kcomp + i); 

            auto tg_xy_xxxxx_0 = braBuffer.data(pidx_g_2_5_m0 + 21 * kcomp + i); 

            auto tg_xy_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + 22 * kcomp + i); 

            auto tg_xy_xxxxz_0 = braBuffer.data(pidx_g_2_5_m0 + 23 * kcomp + i); 

            auto tg_xy_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 24 * kcomp + i); 

            auto tg_xy_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 25 * kcomp + i); 

            auto tg_xy_xxxzz_0 = braBuffer.data(pidx_g_2_5_m0 + 26 * kcomp + i); 

            auto tg_xy_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 27 * kcomp + i); 

            auto tg_xy_xxyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 28 * kcomp + i); 

            auto tg_xy_xxyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 29 * kcomp + i); 

            auto tg_xy_xxzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 30 * kcomp + i); 

            auto tg_xy_xyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 31 * kcomp + i); 

            auto tg_xy_xyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 32 * kcomp + i); 

            auto tg_xy_xyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 33 * kcomp + i); 

            auto tg_xy_xyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 34 * kcomp + i); 

            auto tg_xy_xzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 35 * kcomp + i); 

            auto tg_xy_yyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 36 * kcomp + i); 

            auto tg_xy_yyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 37 * kcomp + i); 

            auto tg_xy_yyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 38 * kcomp + i); 

            auto tg_xy_yyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 39 * kcomp + i); 

            auto tg_xy_yzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 40 * kcomp + i); 

            auto tg_xy_zzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 41 * kcomp + i); 

            auto tg_xz_xxxxx_0 = braBuffer.data(pidx_g_2_5_m0 + 42 * kcomp + i); 

            auto tg_xz_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + 43 * kcomp + i); 

            auto tg_xz_xxxxz_0 = braBuffer.data(pidx_g_2_5_m0 + 44 * kcomp + i); 

            auto tg_xz_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 45 * kcomp + i); 

            auto tg_xz_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 46 * kcomp + i); 

            auto tg_xz_xxxzz_0 = braBuffer.data(pidx_g_2_5_m0 + 47 * kcomp + i); 

            auto tg_xz_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 48 * kcomp + i); 

            auto tg_xz_xxyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 49 * kcomp + i); 

            auto tg_xz_xxyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 50 * kcomp + i); 

            auto tg_xz_xxzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 51 * kcomp + i); 

            auto tg_xz_xyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 52 * kcomp + i); 

            auto tg_xz_xyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 53 * kcomp + i); 

            auto tg_xz_xyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 54 * kcomp + i); 

            auto tg_xz_xyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 55 * kcomp + i); 

            auto tg_xz_xzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 56 * kcomp + i); 

            auto tg_xz_yyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 57 * kcomp + i); 

            auto tg_xz_yyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 58 * kcomp + i); 

            auto tg_xz_yyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 59 * kcomp + i); 

            auto tg_xz_yyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 60 * kcomp + i); 

            auto tg_xz_yzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 61 * kcomp + i); 

            auto tg_xz_zzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 62 * kcomp + i); 

            // Batch of Integrals (0,63)

            #pragma omp simd aligned(tg_x_xxxxx_0, tg_x_xxxxxx_0, tg_x_xxxxxy_0, tg_x_xxxxxz_0, \
                                     tg_x_xxxxy_0, tg_x_xxxxyy_0, tg_x_xxxxyz_0, tg_x_xxxxz_0, tg_x_xxxxzz_0, \
                                     tg_x_xxxyy_0, tg_x_xxxyyy_0, tg_x_xxxyyz_0, tg_x_xxxyz_0, tg_x_xxxyzz_0, \
                                     tg_x_xxxzz_0, tg_x_xxxzzz_0, tg_x_xxyyy_0, tg_x_xxyyyy_0, tg_x_xxyyyz_0, \
                                     tg_x_xxyyz_0, tg_x_xxyyzz_0, tg_x_xxyzz_0, tg_x_xxyzzz_0, tg_x_xxzzz_0, \
                                     tg_x_xxzzzz_0, tg_x_xyyyy_0, tg_x_xyyyyy_0, tg_x_xyyyyz_0, tg_x_xyyyz_0, \
                                     tg_x_xyyyzz_0, tg_x_xyyzz_0, tg_x_xyyzzz_0, tg_x_xyzzz_0, tg_x_xyzzzz_0, \
                                     tg_x_xzzzz_0, tg_x_xzzzzz_0, tg_x_yyyyy_0, tg_x_yyyyz_0, tg_x_yyyzz_0, \
                                     tg_x_yyzzz_0, tg_x_yzzzz_0, tg_x_zzzzz_0, tg_xx_xxxxx_0, tg_xx_xxxxy_0, \
                                     tg_xx_xxxxz_0, tg_xx_xxxyy_0, tg_xx_xxxyz_0, tg_xx_xxxzz_0, tg_xx_xxyyy_0, \
                                     tg_xx_xxyyz_0, tg_xx_xxyzz_0, tg_xx_xxzzz_0, tg_xx_xyyyy_0, tg_xx_xyyyz_0, \
                                     tg_xx_xyyzz_0, tg_xx_xyzzz_0, tg_xx_xzzzz_0, tg_xx_yyyyy_0, tg_xx_yyyyz_0, \
                                     tg_xx_yyyzz_0, tg_xx_yyzzz_0, tg_xx_yzzzz_0, tg_xx_zzzzz_0, tg_xy_xxxxx_0, \
                                     tg_xy_xxxxy_0, tg_xy_xxxxz_0, tg_xy_xxxyy_0, tg_xy_xxxyz_0, tg_xy_xxxzz_0, \
                                     tg_xy_xxyyy_0, tg_xy_xxyyz_0, tg_xy_xxyzz_0, tg_xy_xxzzz_0, tg_xy_xyyyy_0, \
                                     tg_xy_xyyyz_0, tg_xy_xyyzz_0, tg_xy_xyzzz_0, tg_xy_xzzzz_0, tg_xy_yyyyy_0, \
                                     tg_xy_yyyyz_0, tg_xy_yyyzz_0, tg_xy_yyzzz_0, tg_xy_yzzzz_0, tg_xy_zzzzz_0, \
                                     tg_xz_xxxxx_0, tg_xz_xxxxy_0, tg_xz_xxxxz_0, tg_xz_xxxyy_0, tg_xz_xxxyz_0, \
                                     tg_xz_xxxzz_0, tg_xz_xxyyy_0, tg_xz_xxyyz_0, tg_xz_xxyzz_0, tg_xz_xxzzz_0, \
                                     tg_xz_xyyyy_0, tg_xz_xyyyz_0, tg_xz_xyyzz_0, tg_xz_xyzzz_0, tg_xz_xzzzz_0, \
                                     tg_xz_yyyyy_0, tg_xz_yyyyz_0, tg_xz_yyyzz_0, tg_xz_yyzzz_0, tg_xz_yzzzz_0, \
                                     tg_xz_zzzzz_0, tg_y_xxxxx_0, tg_y_xxxxxx_0, tg_y_xxxxxy_0, tg_y_xxxxxz_0, \
                                     tg_y_xxxxy_0, tg_y_xxxxyy_0, tg_y_xxxxyz_0, tg_y_xxxxz_0, tg_y_xxxxzz_0, \
                                     tg_y_xxxyy_0, tg_y_xxxyyy_0, tg_y_xxxyyz_0, tg_y_xxxyz_0, tg_y_xxxyzz_0, \
                                     tg_y_xxxzz_0, tg_y_xxxzzz_0, tg_y_xxyyy_0, tg_y_xxyyyy_0, tg_y_xxyyyz_0, \
                                     tg_y_xxyyz_0, tg_y_xxyyzz_0, tg_y_xxyzz_0, tg_y_xxyzzz_0, tg_y_xxzzz_0, \
                                     tg_y_xxzzzz_0, tg_y_xyyyy_0, tg_y_xyyyyy_0, tg_y_xyyyyz_0, tg_y_xyyyz_0, \
                                     tg_y_xyyyzz_0, tg_y_xyyzz_0, tg_y_xyyzzz_0, tg_y_xyzzz_0, tg_y_xyzzzz_0, \
                                     tg_y_xzzzz_0, tg_y_xzzzzz_0, tg_y_yyyyy_0, tg_y_yyyyz_0, tg_y_yyyzz_0, \
                                     tg_y_yyzzz_0, tg_y_yzzzz_0, tg_y_zzzzz_0, tg_z_xxxxx_0, tg_z_xxxxxx_0, \
                                     tg_z_xxxxxy_0, tg_z_xxxxxz_0, tg_z_xxxxy_0, tg_z_xxxxyy_0, tg_z_xxxxyz_0, \
                                     tg_z_xxxxz_0, tg_z_xxxxzz_0, tg_z_xxxyy_0, tg_z_xxxyyy_0, tg_z_xxxyyz_0, \
                                     tg_z_xxxyz_0, tg_z_xxxyzz_0, tg_z_xxxzz_0, tg_z_xxxzzz_0, tg_z_xxyyy_0, \
                                     tg_z_xxyyyy_0, tg_z_xxyyyz_0, tg_z_xxyyz_0, tg_z_xxyyzz_0, tg_z_xxyzz_0, \
                                     tg_z_xxyzzz_0, tg_z_xxzzz_0, tg_z_xxzzzz_0, tg_z_xyyyy_0, tg_z_xyyyyy_0, \
                                     tg_z_xyyyyz_0, tg_z_xyyyz_0, tg_z_xyyyzz_0, tg_z_xyyzz_0, tg_z_xyyzzz_0, \
                                     tg_z_xyzzz_0, tg_z_xyzzzz_0, tg_z_xzzzz_0, tg_z_xzzzzz_0, tg_z_yyyyy_0, \
                                     tg_z_yyyyz_0, tg_z_yyyzz_0, tg_z_yyzzz_0, tg_z_yzzzz_0, tg_z_zzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xx_xxxxx_0[j] = -ab_x * tg_x_xxxxx_0[j] + tg_x_xxxxxx_0[j];

                tg_xx_xxxxy_0[j] = -ab_x * tg_x_xxxxy_0[j] + tg_x_xxxxxy_0[j];

                tg_xx_xxxxz_0[j] = -ab_x * tg_x_xxxxz_0[j] + tg_x_xxxxxz_0[j];

                tg_xx_xxxyy_0[j] = -ab_x * tg_x_xxxyy_0[j] + tg_x_xxxxyy_0[j];

                tg_xx_xxxyz_0[j] = -ab_x * tg_x_xxxyz_0[j] + tg_x_xxxxyz_0[j];

                tg_xx_xxxzz_0[j] = -ab_x * tg_x_xxxzz_0[j] + tg_x_xxxxzz_0[j];

                tg_xx_xxyyy_0[j] = -ab_x * tg_x_xxyyy_0[j] + tg_x_xxxyyy_0[j];

                tg_xx_xxyyz_0[j] = -ab_x * tg_x_xxyyz_0[j] + tg_x_xxxyyz_0[j];

                tg_xx_xxyzz_0[j] = -ab_x * tg_x_xxyzz_0[j] + tg_x_xxxyzz_0[j];

                tg_xx_xxzzz_0[j] = -ab_x * tg_x_xxzzz_0[j] + tg_x_xxxzzz_0[j];

                tg_xx_xyyyy_0[j] = -ab_x * tg_x_xyyyy_0[j] + tg_x_xxyyyy_0[j];

                tg_xx_xyyyz_0[j] = -ab_x * tg_x_xyyyz_0[j] + tg_x_xxyyyz_0[j];

                tg_xx_xyyzz_0[j] = -ab_x * tg_x_xyyzz_0[j] + tg_x_xxyyzz_0[j];

                tg_xx_xyzzz_0[j] = -ab_x * tg_x_xyzzz_0[j] + tg_x_xxyzzz_0[j];

                tg_xx_xzzzz_0[j] = -ab_x * tg_x_xzzzz_0[j] + tg_x_xxzzzz_0[j];

                tg_xx_yyyyy_0[j] = -ab_x * tg_x_yyyyy_0[j] + tg_x_xyyyyy_0[j];

                tg_xx_yyyyz_0[j] = -ab_x * tg_x_yyyyz_0[j] + tg_x_xyyyyz_0[j];

                tg_xx_yyyzz_0[j] = -ab_x * tg_x_yyyzz_0[j] + tg_x_xyyyzz_0[j];

                tg_xx_yyzzz_0[j] = -ab_x * tg_x_yyzzz_0[j] + tg_x_xyyzzz_0[j];

                tg_xx_yzzzz_0[j] = -ab_x * tg_x_yzzzz_0[j] + tg_x_xyzzzz_0[j];

                tg_xx_zzzzz_0[j] = -ab_x * tg_x_zzzzz_0[j] + tg_x_xzzzzz_0[j];

                tg_xy_xxxxx_0[j] = -ab_x * tg_y_xxxxx_0[j] + tg_y_xxxxxx_0[j];

                tg_xy_xxxxy_0[j] = -ab_x * tg_y_xxxxy_0[j] + tg_y_xxxxxy_0[j];

                tg_xy_xxxxz_0[j] = -ab_x * tg_y_xxxxz_0[j] + tg_y_xxxxxz_0[j];

                tg_xy_xxxyy_0[j] = -ab_x * tg_y_xxxyy_0[j] + tg_y_xxxxyy_0[j];

                tg_xy_xxxyz_0[j] = -ab_x * tg_y_xxxyz_0[j] + tg_y_xxxxyz_0[j];

                tg_xy_xxxzz_0[j] = -ab_x * tg_y_xxxzz_0[j] + tg_y_xxxxzz_0[j];

                tg_xy_xxyyy_0[j] = -ab_x * tg_y_xxyyy_0[j] + tg_y_xxxyyy_0[j];

                tg_xy_xxyyz_0[j] = -ab_x * tg_y_xxyyz_0[j] + tg_y_xxxyyz_0[j];

                tg_xy_xxyzz_0[j] = -ab_x * tg_y_xxyzz_0[j] + tg_y_xxxyzz_0[j];

                tg_xy_xxzzz_0[j] = -ab_x * tg_y_xxzzz_0[j] + tg_y_xxxzzz_0[j];

                tg_xy_xyyyy_0[j] = -ab_x * tg_y_xyyyy_0[j] + tg_y_xxyyyy_0[j];

                tg_xy_xyyyz_0[j] = -ab_x * tg_y_xyyyz_0[j] + tg_y_xxyyyz_0[j];

                tg_xy_xyyzz_0[j] = -ab_x * tg_y_xyyzz_0[j] + tg_y_xxyyzz_0[j];

                tg_xy_xyzzz_0[j] = -ab_x * tg_y_xyzzz_0[j] + tg_y_xxyzzz_0[j];

                tg_xy_xzzzz_0[j] = -ab_x * tg_y_xzzzz_0[j] + tg_y_xxzzzz_0[j];

                tg_xy_yyyyy_0[j] = -ab_x * tg_y_yyyyy_0[j] + tg_y_xyyyyy_0[j];

                tg_xy_yyyyz_0[j] = -ab_x * tg_y_yyyyz_0[j] + tg_y_xyyyyz_0[j];

                tg_xy_yyyzz_0[j] = -ab_x * tg_y_yyyzz_0[j] + tg_y_xyyyzz_0[j];

                tg_xy_yyzzz_0[j] = -ab_x * tg_y_yyzzz_0[j] + tg_y_xyyzzz_0[j];

                tg_xy_yzzzz_0[j] = -ab_x * tg_y_yzzzz_0[j] + tg_y_xyzzzz_0[j];

                tg_xy_zzzzz_0[j] = -ab_x * tg_y_zzzzz_0[j] + tg_y_xzzzzz_0[j];

                tg_xz_xxxxx_0[j] = -ab_x * tg_z_xxxxx_0[j] + tg_z_xxxxxx_0[j];

                tg_xz_xxxxy_0[j] = -ab_x * tg_z_xxxxy_0[j] + tg_z_xxxxxy_0[j];

                tg_xz_xxxxz_0[j] = -ab_x * tg_z_xxxxz_0[j] + tg_z_xxxxxz_0[j];

                tg_xz_xxxyy_0[j] = -ab_x * tg_z_xxxyy_0[j] + tg_z_xxxxyy_0[j];

                tg_xz_xxxyz_0[j] = -ab_x * tg_z_xxxyz_0[j] + tg_z_xxxxyz_0[j];

                tg_xz_xxxzz_0[j] = -ab_x * tg_z_xxxzz_0[j] + tg_z_xxxxzz_0[j];

                tg_xz_xxyyy_0[j] = -ab_x * tg_z_xxyyy_0[j] + tg_z_xxxyyy_0[j];

                tg_xz_xxyyz_0[j] = -ab_x * tg_z_xxyyz_0[j] + tg_z_xxxyyz_0[j];

                tg_xz_xxyzz_0[j] = -ab_x * tg_z_xxyzz_0[j] + tg_z_xxxyzz_0[j];

                tg_xz_xxzzz_0[j] = -ab_x * tg_z_xxzzz_0[j] + tg_z_xxxzzz_0[j];

                tg_xz_xyyyy_0[j] = -ab_x * tg_z_xyyyy_0[j] + tg_z_xxyyyy_0[j];

                tg_xz_xyyyz_0[j] = -ab_x * tg_z_xyyyz_0[j] + tg_z_xxyyyz_0[j];

                tg_xz_xyyzz_0[j] = -ab_x * tg_z_xyyzz_0[j] + tg_z_xxyyzz_0[j];

                tg_xz_xyzzz_0[j] = -ab_x * tg_z_xyzzz_0[j] + tg_z_xxyzzz_0[j];

                tg_xz_xzzzz_0[j] = -ab_x * tg_z_xzzzz_0[j] + tg_z_xxzzzz_0[j];

                tg_xz_yyyyy_0[j] = -ab_x * tg_z_yyyyy_0[j] + tg_z_xyyyyy_0[j];

                tg_xz_yyyyz_0[j] = -ab_x * tg_z_yyyyz_0[j] + tg_z_xyyyyz_0[j];

                tg_xz_yyyzz_0[j] = -ab_x * tg_z_yyyzz_0[j] + tg_z_xyyyzz_0[j];

                tg_xz_yyzzz_0[j] = -ab_x * tg_z_yyzzz_0[j] + tg_z_xyyzzz_0[j];

                tg_xz_yzzzz_0[j] = -ab_x * tg_z_yzzzz_0[j] + tg_z_xyzzzz_0[j];

                tg_xz_zzzzz_0[j] = -ab_x * tg_z_zzzzz_0[j] + tg_z_xzzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForDHXY_63_126(      CMemBlock2D<double>& braBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& abDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetContrPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (63,126)

        // set up distances R(AB) = A - B

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_2_5_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_1_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_1_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_y_xxxxx_0 = braBuffer.data(pidx_g_1_5_m0 + 21 * kcomp + i); 

            auto tg_y_xxxxy_0 = braBuffer.data(pidx_g_1_5_m0 + 22 * kcomp + i); 

            auto tg_y_xxxxz_0 = braBuffer.data(pidx_g_1_5_m0 + 23 * kcomp + i); 

            auto tg_y_xxxyy_0 = braBuffer.data(pidx_g_1_5_m0 + 24 * kcomp + i); 

            auto tg_y_xxxyz_0 = braBuffer.data(pidx_g_1_5_m0 + 25 * kcomp + i); 

            auto tg_y_xxxzz_0 = braBuffer.data(pidx_g_1_5_m0 + 26 * kcomp + i); 

            auto tg_y_xxyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 27 * kcomp + i); 

            auto tg_y_xxyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 28 * kcomp + i); 

            auto tg_y_xxyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 29 * kcomp + i); 

            auto tg_y_xxzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 30 * kcomp + i); 

            auto tg_y_xyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 31 * kcomp + i); 

            auto tg_y_xyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 32 * kcomp + i); 

            auto tg_y_xyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 33 * kcomp + i); 

            auto tg_y_xyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 34 * kcomp + i); 

            auto tg_y_xzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 35 * kcomp + i); 

            auto tg_y_yyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 36 * kcomp + i); 

            auto tg_y_yyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 37 * kcomp + i); 

            auto tg_y_yyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 38 * kcomp + i); 

            auto tg_y_yyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 39 * kcomp + i); 

            auto tg_y_yzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 40 * kcomp + i); 

            auto tg_y_zzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 41 * kcomp + i); 

            auto tg_z_xxxxx_0 = braBuffer.data(pidx_g_1_5_m0 + 42 * kcomp + i); 

            auto tg_z_xxxxy_0 = braBuffer.data(pidx_g_1_5_m0 + 43 * kcomp + i); 

            auto tg_z_xxxxz_0 = braBuffer.data(pidx_g_1_5_m0 + 44 * kcomp + i); 

            auto tg_z_xxxyy_0 = braBuffer.data(pidx_g_1_5_m0 + 45 * kcomp + i); 

            auto tg_z_xxxyz_0 = braBuffer.data(pidx_g_1_5_m0 + 46 * kcomp + i); 

            auto tg_z_xxxzz_0 = braBuffer.data(pidx_g_1_5_m0 + 47 * kcomp + i); 

            auto tg_z_xxyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 48 * kcomp + i); 

            auto tg_z_xxyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 49 * kcomp + i); 

            auto tg_z_xxyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 50 * kcomp + i); 

            auto tg_z_xxzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 51 * kcomp + i); 

            auto tg_z_xyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 52 * kcomp + i); 

            auto tg_z_xyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 53 * kcomp + i); 

            auto tg_z_xyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 54 * kcomp + i); 

            auto tg_z_xyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 55 * kcomp + i); 

            auto tg_z_xzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 56 * kcomp + i); 

            auto tg_z_yyyyy_0 = braBuffer.data(pidx_g_1_5_m0 + 57 * kcomp + i); 

            auto tg_z_yyyyz_0 = braBuffer.data(pidx_g_1_5_m0 + 58 * kcomp + i); 

            auto tg_z_yyyzz_0 = braBuffer.data(pidx_g_1_5_m0 + 59 * kcomp + i); 

            auto tg_z_yyzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 60 * kcomp + i); 

            auto tg_z_yzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 61 * kcomp + i); 

            auto tg_z_zzzzz_0 = braBuffer.data(pidx_g_1_5_m0 + 62 * kcomp + i); 

            auto tg_y_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + 29 * kcomp + i); 

            auto tg_y_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 31 * kcomp + i); 

            auto tg_y_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 32 * kcomp + i); 

            auto tg_y_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 34 * kcomp + i); 

            auto tg_y_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 35 * kcomp + i); 

            auto tg_y_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 36 * kcomp + i); 

            auto tg_y_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 38 * kcomp + i); 

            auto tg_y_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 39 * kcomp + i); 

            auto tg_y_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 40 * kcomp + i); 

            auto tg_y_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 41 * kcomp + i); 

            auto tg_y_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 43 * kcomp + i); 

            auto tg_y_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 44 * kcomp + i); 

            auto tg_y_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 45 * kcomp + i); 

            auto tg_y_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 46 * kcomp + i); 

            auto tg_y_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 47 * kcomp + i); 

            auto tg_y_yyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 49 * kcomp + i); 

            auto tg_y_yyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 50 * kcomp + i); 

            auto tg_y_yyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 51 * kcomp + i); 

            auto tg_y_yyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 52 * kcomp + i); 

            auto tg_y_yyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 53 * kcomp + i); 

            auto tg_y_yzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 54 * kcomp + i); 

            auto tg_z_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + 57 * kcomp + i); 

            auto tg_z_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 58 * kcomp + i); 

            auto tg_z_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 59 * kcomp + i); 

            auto tg_z_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 60 * kcomp + i); 

            auto tg_z_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 61 * kcomp + i); 

            auto tg_z_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 62 * kcomp + i); 

            auto tg_z_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 63 * kcomp + i); 

            auto tg_z_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 64 * kcomp + i); 

            auto tg_z_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 65 * kcomp + i); 

            auto tg_z_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 66 * kcomp + i); 

            auto tg_z_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 67 * kcomp + i); 

            auto tg_z_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 68 * kcomp + i); 

            auto tg_z_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 69 * kcomp + i); 

            auto tg_z_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 70 * kcomp + i); 

            auto tg_z_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 71 * kcomp + i); 

            auto tg_z_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 72 * kcomp + i); 

            auto tg_z_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 73 * kcomp + i); 

            auto tg_z_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 74 * kcomp + i); 

            auto tg_z_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 75 * kcomp + i); 

            auto tg_z_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 76 * kcomp + i); 

            auto tg_z_yyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 77 * kcomp + i); 

            auto tg_z_yyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 78 * kcomp + i); 

            auto tg_z_yyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 79 * kcomp + i); 

            auto tg_z_yyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 80 * kcomp + i); 

            auto tg_z_yyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 81 * kcomp + i); 

            auto tg_z_yzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 82 * kcomp + i); 

            auto tg_z_zzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 83 * kcomp + i); 

            // set up pointers to integrals

            auto tg_yy_xxxxx_0 = braBuffer.data(pidx_g_2_5_m0 + 63 * kcomp + i); 

            auto tg_yy_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + 64 * kcomp + i); 

            auto tg_yy_xxxxz_0 = braBuffer.data(pidx_g_2_5_m0 + 65 * kcomp + i); 

            auto tg_yy_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 66 * kcomp + i); 

            auto tg_yy_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 67 * kcomp + i); 

            auto tg_yy_xxxzz_0 = braBuffer.data(pidx_g_2_5_m0 + 68 * kcomp + i); 

            auto tg_yy_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 69 * kcomp + i); 

            auto tg_yy_xxyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 70 * kcomp + i); 

            auto tg_yy_xxyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 71 * kcomp + i); 

            auto tg_yy_xxzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 72 * kcomp + i); 

            auto tg_yy_xyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 73 * kcomp + i); 

            auto tg_yy_xyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 74 * kcomp + i); 

            auto tg_yy_xyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 75 * kcomp + i); 

            auto tg_yy_xyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 76 * kcomp + i); 

            auto tg_yy_xzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 77 * kcomp + i); 

            auto tg_yy_yyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 78 * kcomp + i); 

            auto tg_yy_yyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 79 * kcomp + i); 

            auto tg_yy_yyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 80 * kcomp + i); 

            auto tg_yy_yyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 81 * kcomp + i); 

            auto tg_yy_yzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 82 * kcomp + i); 

            auto tg_yy_zzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 83 * kcomp + i); 

            auto tg_yz_xxxxx_0 = braBuffer.data(pidx_g_2_5_m0 + 84 * kcomp + i); 

            auto tg_yz_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + 85 * kcomp + i); 

            auto tg_yz_xxxxz_0 = braBuffer.data(pidx_g_2_5_m0 + 86 * kcomp + i); 

            auto tg_yz_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 87 * kcomp + i); 

            auto tg_yz_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 88 * kcomp + i); 

            auto tg_yz_xxxzz_0 = braBuffer.data(pidx_g_2_5_m0 + 89 * kcomp + i); 

            auto tg_yz_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 90 * kcomp + i); 

            auto tg_yz_xxyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 91 * kcomp + i); 

            auto tg_yz_xxyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 92 * kcomp + i); 

            auto tg_yz_xxzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 93 * kcomp + i); 

            auto tg_yz_xyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 94 * kcomp + i); 

            auto tg_yz_xyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 95 * kcomp + i); 

            auto tg_yz_xyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 96 * kcomp + i); 

            auto tg_yz_xyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 97 * kcomp + i); 

            auto tg_yz_xzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 98 * kcomp + i); 

            auto tg_yz_yyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 99 * kcomp + i); 

            auto tg_yz_yyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 100 * kcomp + i); 

            auto tg_yz_yyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 101 * kcomp + i); 

            auto tg_yz_yyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 102 * kcomp + i); 

            auto tg_yz_yzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 103 * kcomp + i); 

            auto tg_yz_zzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 104 * kcomp + i); 

            auto tg_zz_xxxxx_0 = braBuffer.data(pidx_g_2_5_m0 + 105 * kcomp + i); 

            auto tg_zz_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + 106 * kcomp + i); 

            auto tg_zz_xxxxz_0 = braBuffer.data(pidx_g_2_5_m0 + 107 * kcomp + i); 

            auto tg_zz_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 108 * kcomp + i); 

            auto tg_zz_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 109 * kcomp + i); 

            auto tg_zz_xxxzz_0 = braBuffer.data(pidx_g_2_5_m0 + 110 * kcomp + i); 

            auto tg_zz_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 111 * kcomp + i); 

            auto tg_zz_xxyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 112 * kcomp + i); 

            auto tg_zz_xxyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 113 * kcomp + i); 

            auto tg_zz_xxzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 114 * kcomp + i); 

            auto tg_zz_xyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 115 * kcomp + i); 

            auto tg_zz_xyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 116 * kcomp + i); 

            auto tg_zz_xyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 117 * kcomp + i); 

            auto tg_zz_xyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 118 * kcomp + i); 

            auto tg_zz_xzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 119 * kcomp + i); 

            auto tg_zz_yyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 120 * kcomp + i); 

            auto tg_zz_yyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 121 * kcomp + i); 

            auto tg_zz_yyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 122 * kcomp + i); 

            auto tg_zz_yyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 123 * kcomp + i); 

            auto tg_zz_yzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 124 * kcomp + i); 

            auto tg_zz_zzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 125 * kcomp + i); 

            // Batch of Integrals (63,126)

            #pragma omp simd aligned(tg_y_xxxxx_0, tg_y_xxxxxy_0, tg_y_xxxxy_0, tg_y_xxxxyy_0, \
                                     tg_y_xxxxyz_0, tg_y_xxxxz_0, tg_y_xxxyy_0, tg_y_xxxyyy_0, tg_y_xxxyyz_0, \
                                     tg_y_xxxyz_0, tg_y_xxxyzz_0, tg_y_xxxzz_0, tg_y_xxyyy_0, tg_y_xxyyyy_0, \
                                     tg_y_xxyyyz_0, tg_y_xxyyz_0, tg_y_xxyyzz_0, tg_y_xxyzz_0, tg_y_xxyzzz_0, \
                                     tg_y_xxzzz_0, tg_y_xyyyy_0, tg_y_xyyyyy_0, tg_y_xyyyyz_0, tg_y_xyyyz_0, \
                                     tg_y_xyyyzz_0, tg_y_xyyzz_0, tg_y_xyyzzz_0, tg_y_xyzzz_0, tg_y_xyzzzz_0, \
                                     tg_y_xzzzz_0, tg_y_yyyyy_0, tg_y_yyyyyy_0, tg_y_yyyyyz_0, tg_y_yyyyz_0, \
                                     tg_y_yyyyzz_0, tg_y_yyyzz_0, tg_y_yyyzzz_0, tg_y_yyzzz_0, tg_y_yyzzzz_0, \
                                     tg_y_yzzzz_0, tg_y_yzzzzz_0, tg_y_zzzzz_0, tg_yy_xxxxx_0, tg_yy_xxxxy_0, \
                                     tg_yy_xxxxz_0, tg_yy_xxxyy_0, tg_yy_xxxyz_0, tg_yy_xxxzz_0, tg_yy_xxyyy_0, \
                                     tg_yy_xxyyz_0, tg_yy_xxyzz_0, tg_yy_xxzzz_0, tg_yy_xyyyy_0, tg_yy_xyyyz_0, \
                                     tg_yy_xyyzz_0, tg_yy_xyzzz_0, tg_yy_xzzzz_0, tg_yy_yyyyy_0, tg_yy_yyyyz_0, \
                                     tg_yy_yyyzz_0, tg_yy_yyzzz_0, tg_yy_yzzzz_0, tg_yy_zzzzz_0, tg_yz_xxxxx_0, \
                                     tg_yz_xxxxy_0, tg_yz_xxxxz_0, tg_yz_xxxyy_0, tg_yz_xxxyz_0, tg_yz_xxxzz_0, \
                                     tg_yz_xxyyy_0, tg_yz_xxyyz_0, tg_yz_xxyzz_0, tg_yz_xxzzz_0, tg_yz_xyyyy_0, \
                                     tg_yz_xyyyz_0, tg_yz_xyyzz_0, tg_yz_xyzzz_0, tg_yz_xzzzz_0, tg_yz_yyyyy_0, \
                                     tg_yz_yyyyz_0, tg_yz_yyyzz_0, tg_yz_yyzzz_0, tg_yz_yzzzz_0, tg_yz_zzzzz_0, \
                                     tg_z_xxxxx_0, tg_z_xxxxxy_0, tg_z_xxxxxz_0, tg_z_xxxxy_0, tg_z_xxxxyy_0, \
                                     tg_z_xxxxyz_0, tg_z_xxxxz_0, tg_z_xxxxzz_0, tg_z_xxxyy_0, tg_z_xxxyyy_0, \
                                     tg_z_xxxyyz_0, tg_z_xxxyz_0, tg_z_xxxyzz_0, tg_z_xxxzz_0, tg_z_xxxzzz_0, \
                                     tg_z_xxyyy_0, tg_z_xxyyyy_0, tg_z_xxyyyz_0, tg_z_xxyyz_0, tg_z_xxyyzz_0, \
                                     tg_z_xxyzz_0, tg_z_xxyzzz_0, tg_z_xxzzz_0, tg_z_xxzzzz_0, tg_z_xyyyy_0, \
                                     tg_z_xyyyyy_0, tg_z_xyyyyz_0, tg_z_xyyyz_0, tg_z_xyyyzz_0, tg_z_xyyzz_0, \
                                     tg_z_xyyzzz_0, tg_z_xyzzz_0, tg_z_xyzzzz_0, tg_z_xzzzz_0, tg_z_xzzzzz_0, \
                                     tg_z_yyyyy_0, tg_z_yyyyyy_0, tg_z_yyyyyz_0, tg_z_yyyyz_0, tg_z_yyyyzz_0, \
                                     tg_z_yyyzz_0, tg_z_yyyzzz_0, tg_z_yyzzz_0, tg_z_yyzzzz_0, tg_z_yzzzz_0, \
                                     tg_z_yzzzzz_0, tg_z_zzzzz_0, tg_z_zzzzzz_0, tg_zz_xxxxx_0, tg_zz_xxxxy_0, \
                                     tg_zz_xxxxz_0, tg_zz_xxxyy_0, tg_zz_xxxyz_0, tg_zz_xxxzz_0, tg_zz_xxyyy_0, \
                                     tg_zz_xxyyz_0, tg_zz_xxyzz_0, tg_zz_xxzzz_0, tg_zz_xyyyy_0, tg_zz_xyyyz_0, \
                                     tg_zz_xyyzz_0, tg_zz_xyzzz_0, tg_zz_xzzzz_0, tg_zz_yyyyy_0, tg_zz_yyyyz_0, \
                                     tg_zz_yyyzz_0, tg_zz_yyzzz_0, tg_zz_yzzzz_0, tg_zz_zzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_yy_xxxxx_0[j] = -ab_y * tg_y_xxxxx_0[j] + tg_y_xxxxxy_0[j];

                tg_yy_xxxxy_0[j] = -ab_y * tg_y_xxxxy_0[j] + tg_y_xxxxyy_0[j];

                tg_yy_xxxxz_0[j] = -ab_y * tg_y_xxxxz_0[j] + tg_y_xxxxyz_0[j];

                tg_yy_xxxyy_0[j] = -ab_y * tg_y_xxxyy_0[j] + tg_y_xxxyyy_0[j];

                tg_yy_xxxyz_0[j] = -ab_y * tg_y_xxxyz_0[j] + tg_y_xxxyyz_0[j];

                tg_yy_xxxzz_0[j] = -ab_y * tg_y_xxxzz_0[j] + tg_y_xxxyzz_0[j];

                tg_yy_xxyyy_0[j] = -ab_y * tg_y_xxyyy_0[j] + tg_y_xxyyyy_0[j];

                tg_yy_xxyyz_0[j] = -ab_y * tg_y_xxyyz_0[j] + tg_y_xxyyyz_0[j];

                tg_yy_xxyzz_0[j] = -ab_y * tg_y_xxyzz_0[j] + tg_y_xxyyzz_0[j];

                tg_yy_xxzzz_0[j] = -ab_y * tg_y_xxzzz_0[j] + tg_y_xxyzzz_0[j];

                tg_yy_xyyyy_0[j] = -ab_y * tg_y_xyyyy_0[j] + tg_y_xyyyyy_0[j];

                tg_yy_xyyyz_0[j] = -ab_y * tg_y_xyyyz_0[j] + tg_y_xyyyyz_0[j];

                tg_yy_xyyzz_0[j] = -ab_y * tg_y_xyyzz_0[j] + tg_y_xyyyzz_0[j];

                tg_yy_xyzzz_0[j] = -ab_y * tg_y_xyzzz_0[j] + tg_y_xyyzzz_0[j];

                tg_yy_xzzzz_0[j] = -ab_y * tg_y_xzzzz_0[j] + tg_y_xyzzzz_0[j];

                tg_yy_yyyyy_0[j] = -ab_y * tg_y_yyyyy_0[j] + tg_y_yyyyyy_0[j];

                tg_yy_yyyyz_0[j] = -ab_y * tg_y_yyyyz_0[j] + tg_y_yyyyyz_0[j];

                tg_yy_yyyzz_0[j] = -ab_y * tg_y_yyyzz_0[j] + tg_y_yyyyzz_0[j];

                tg_yy_yyzzz_0[j] = -ab_y * tg_y_yyzzz_0[j] + tg_y_yyyzzz_0[j];

                tg_yy_yzzzz_0[j] = -ab_y * tg_y_yzzzz_0[j] + tg_y_yyzzzz_0[j];

                tg_yy_zzzzz_0[j] = -ab_y * tg_y_zzzzz_0[j] + tg_y_yzzzzz_0[j];

                tg_yz_xxxxx_0[j] = -ab_y * tg_z_xxxxx_0[j] + tg_z_xxxxxy_0[j];

                tg_yz_xxxxy_0[j] = -ab_y * tg_z_xxxxy_0[j] + tg_z_xxxxyy_0[j];

                tg_yz_xxxxz_0[j] = -ab_y * tg_z_xxxxz_0[j] + tg_z_xxxxyz_0[j];

                tg_yz_xxxyy_0[j] = -ab_y * tg_z_xxxyy_0[j] + tg_z_xxxyyy_0[j];

                tg_yz_xxxyz_0[j] = -ab_y * tg_z_xxxyz_0[j] + tg_z_xxxyyz_0[j];

                tg_yz_xxxzz_0[j] = -ab_y * tg_z_xxxzz_0[j] + tg_z_xxxyzz_0[j];

                tg_yz_xxyyy_0[j] = -ab_y * tg_z_xxyyy_0[j] + tg_z_xxyyyy_0[j];

                tg_yz_xxyyz_0[j] = -ab_y * tg_z_xxyyz_0[j] + tg_z_xxyyyz_0[j];

                tg_yz_xxyzz_0[j] = -ab_y * tg_z_xxyzz_0[j] + tg_z_xxyyzz_0[j];

                tg_yz_xxzzz_0[j] = -ab_y * tg_z_xxzzz_0[j] + tg_z_xxyzzz_0[j];

                tg_yz_xyyyy_0[j] = -ab_y * tg_z_xyyyy_0[j] + tg_z_xyyyyy_0[j];

                tg_yz_xyyyz_0[j] = -ab_y * tg_z_xyyyz_0[j] + tg_z_xyyyyz_0[j];

                tg_yz_xyyzz_0[j] = -ab_y * tg_z_xyyzz_0[j] + tg_z_xyyyzz_0[j];

                tg_yz_xyzzz_0[j] = -ab_y * tg_z_xyzzz_0[j] + tg_z_xyyzzz_0[j];

                tg_yz_xzzzz_0[j] = -ab_y * tg_z_xzzzz_0[j] + tg_z_xyzzzz_0[j];

                tg_yz_yyyyy_0[j] = -ab_y * tg_z_yyyyy_0[j] + tg_z_yyyyyy_0[j];

                tg_yz_yyyyz_0[j] = -ab_y * tg_z_yyyyz_0[j] + tg_z_yyyyyz_0[j];

                tg_yz_yyyzz_0[j] = -ab_y * tg_z_yyyzz_0[j] + tg_z_yyyyzz_0[j];

                tg_yz_yyzzz_0[j] = -ab_y * tg_z_yyzzz_0[j] + tg_z_yyyzzz_0[j];

                tg_yz_yzzzz_0[j] = -ab_y * tg_z_yzzzz_0[j] + tg_z_yyzzzz_0[j];

                tg_yz_zzzzz_0[j] = -ab_y * tg_z_zzzzz_0[j] + tg_z_yzzzzz_0[j];

                tg_zz_xxxxx_0[j] = -ab_z * tg_z_xxxxx_0[j] + tg_z_xxxxxz_0[j];

                tg_zz_xxxxy_0[j] = -ab_z * tg_z_xxxxy_0[j] + tg_z_xxxxyz_0[j];

                tg_zz_xxxxz_0[j] = -ab_z * tg_z_xxxxz_0[j] + tg_z_xxxxzz_0[j];

                tg_zz_xxxyy_0[j] = -ab_z * tg_z_xxxyy_0[j] + tg_z_xxxyyz_0[j];

                tg_zz_xxxyz_0[j] = -ab_z * tg_z_xxxyz_0[j] + tg_z_xxxyzz_0[j];

                tg_zz_xxxzz_0[j] = -ab_z * tg_z_xxxzz_0[j] + tg_z_xxxzzz_0[j];

                tg_zz_xxyyy_0[j] = -ab_z * tg_z_xxyyy_0[j] + tg_z_xxyyyz_0[j];

                tg_zz_xxyyz_0[j] = -ab_z * tg_z_xxyyz_0[j] + tg_z_xxyyzz_0[j];

                tg_zz_xxyzz_0[j] = -ab_z * tg_z_xxyzz_0[j] + tg_z_xxyzzz_0[j];

                tg_zz_xxzzz_0[j] = -ab_z * tg_z_xxzzz_0[j] + tg_z_xxzzzz_0[j];

                tg_zz_xyyyy_0[j] = -ab_z * tg_z_xyyyy_0[j] + tg_z_xyyyyz_0[j];

                tg_zz_xyyyz_0[j] = -ab_z * tg_z_xyyyz_0[j] + tg_z_xyyyzz_0[j];

                tg_zz_xyyzz_0[j] = -ab_z * tg_z_xyyzz_0[j] + tg_z_xyyzzz_0[j];

                tg_zz_xyzzz_0[j] = -ab_z * tg_z_xyzzz_0[j] + tg_z_xyzzzz_0[j];

                tg_zz_xzzzz_0[j] = -ab_z * tg_z_xzzzz_0[j] + tg_z_xzzzzz_0[j];

                tg_zz_yyyyy_0[j] = -ab_z * tg_z_yyyyy_0[j] + tg_z_yyyyyz_0[j];

                tg_zz_yyyyz_0[j] = -ab_z * tg_z_yyyyz_0[j] + tg_z_yyyyzz_0[j];

                tg_zz_yyyzz_0[j] = -ab_z * tg_z_yyyzz_0[j] + tg_z_yyyzzz_0[j];

                tg_zz_yyzzz_0[j] = -ab_z * tg_z_yyzzz_0[j] + tg_z_yyzzzz_0[j];

                tg_zz_yzzzz_0[j] = -ab_z * tg_z_yzzzz_0[j] + tg_z_yzzzzz_0[j];

                tg_zz_zzzzz_0[j] = -ab_z * tg_z_zzzzz_0[j] + tg_z_zzzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForDIXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        eribrrfunc::compElectronRepulsionForDIXY_0_84(braBuffer,
                                                      recursionMap,
                                                      abDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        eribrrfunc::compElectronRepulsionForDIXY_84_168(braBuffer,
                                                        recursionMap,
                                                        abDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForDIXY_0_84(      CMemBlock2D<double>& braBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& abDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,84)

        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_2_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_2_6_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_1_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 7, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_x_xxxxxx_0 = braBuffer.data(pidx_g_1_6_m0 + i); 

            auto tg_x_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + kcomp + i); 

            auto tg_x_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 2 * kcomp + i); 

            auto tg_x_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 3 * kcomp + i); 

            auto tg_x_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 4 * kcomp + i); 

            auto tg_x_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 5 * kcomp + i); 

            auto tg_x_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 6 * kcomp + i); 

            auto tg_x_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 7 * kcomp + i); 

            auto tg_x_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 8 * kcomp + i); 

            auto tg_x_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 9 * kcomp + i); 

            auto tg_x_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 10 * kcomp + i); 

            auto tg_x_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 11 * kcomp + i); 

            auto tg_x_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 12 * kcomp + i); 

            auto tg_x_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 13 * kcomp + i); 

            auto tg_x_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 14 * kcomp + i); 

            auto tg_x_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 15 * kcomp + i); 

            auto tg_x_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 16 * kcomp + i); 

            auto tg_x_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 17 * kcomp + i); 

            auto tg_x_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 18 * kcomp + i); 

            auto tg_x_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 19 * kcomp + i); 

            auto tg_x_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 20 * kcomp + i); 

            auto tg_x_yyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 21 * kcomp + i); 

            auto tg_x_yyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 22 * kcomp + i); 

            auto tg_x_yyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 23 * kcomp + i); 

            auto tg_x_yyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 24 * kcomp + i); 

            auto tg_x_yyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 25 * kcomp + i); 

            auto tg_x_yzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 26 * kcomp + i); 

            auto tg_x_zzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 27 * kcomp + i); 

            auto tg_y_xxxxxx_0 = braBuffer.data(pidx_g_1_6_m0 + 28 * kcomp + i); 

            auto tg_y_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + 29 * kcomp + i); 

            auto tg_y_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 30 * kcomp + i); 

            auto tg_y_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 31 * kcomp + i); 

            auto tg_y_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 32 * kcomp + i); 

            auto tg_y_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 33 * kcomp + i); 

            auto tg_y_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 34 * kcomp + i); 

            auto tg_y_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 35 * kcomp + i); 

            auto tg_y_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 36 * kcomp + i); 

            auto tg_y_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 37 * kcomp + i); 

            auto tg_y_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 38 * kcomp + i); 

            auto tg_y_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 39 * kcomp + i); 

            auto tg_y_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 40 * kcomp + i); 

            auto tg_y_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 41 * kcomp + i); 

            auto tg_y_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 42 * kcomp + i); 

            auto tg_y_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 43 * kcomp + i); 

            auto tg_y_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 44 * kcomp + i); 

            auto tg_y_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 45 * kcomp + i); 

            auto tg_y_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 46 * kcomp + i); 

            auto tg_y_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 47 * kcomp + i); 

            auto tg_y_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 48 * kcomp + i); 

            auto tg_y_yyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 49 * kcomp + i); 

            auto tg_y_yyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 50 * kcomp + i); 

            auto tg_y_yyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 51 * kcomp + i); 

            auto tg_y_yyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 52 * kcomp + i); 

            auto tg_y_yyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 53 * kcomp + i); 

            auto tg_y_yzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 54 * kcomp + i); 

            auto tg_y_zzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 55 * kcomp + i); 

            auto tg_z_xxxxxx_0 = braBuffer.data(pidx_g_1_6_m0 + 56 * kcomp + i); 

            auto tg_z_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + 57 * kcomp + i); 

            auto tg_z_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 58 * kcomp + i); 

            auto tg_z_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 59 * kcomp + i); 

            auto tg_z_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 60 * kcomp + i); 

            auto tg_z_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 61 * kcomp + i); 

            auto tg_z_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 62 * kcomp + i); 

            auto tg_z_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 63 * kcomp + i); 

            auto tg_z_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 64 * kcomp + i); 

            auto tg_z_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 65 * kcomp + i); 

            auto tg_z_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 66 * kcomp + i); 

            auto tg_z_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 67 * kcomp + i); 

            auto tg_z_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 68 * kcomp + i); 

            auto tg_z_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 69 * kcomp + i); 

            auto tg_z_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 70 * kcomp + i); 

            auto tg_z_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 71 * kcomp + i); 

            auto tg_z_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 72 * kcomp + i); 

            auto tg_z_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 73 * kcomp + i); 

            auto tg_z_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 74 * kcomp + i); 

            auto tg_z_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 75 * kcomp + i); 

            auto tg_z_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 76 * kcomp + i); 

            auto tg_z_yyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 77 * kcomp + i); 

            auto tg_z_yyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 78 * kcomp + i); 

            auto tg_z_yyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 79 * kcomp + i); 

            auto tg_z_yyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 80 * kcomp + i); 

            auto tg_z_yyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 81 * kcomp + i); 

            auto tg_z_yzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 82 * kcomp + i); 

            auto tg_z_zzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 83 * kcomp + i); 

            auto tg_x_xxxxxxx_0 = braBuffer.data(pidx_g_1_7_m0 + i); 

            auto tg_x_xxxxxxy_0 = braBuffer.data(pidx_g_1_7_m0 + kcomp + i); 

            auto tg_x_xxxxxxz_0 = braBuffer.data(pidx_g_1_7_m0 + 2 * kcomp + i); 

            auto tg_x_xxxxxyy_0 = braBuffer.data(pidx_g_1_7_m0 + 3 * kcomp + i); 

            auto tg_x_xxxxxyz_0 = braBuffer.data(pidx_g_1_7_m0 + 4 * kcomp + i); 

            auto tg_x_xxxxxzz_0 = braBuffer.data(pidx_g_1_7_m0 + 5 * kcomp + i); 

            auto tg_x_xxxxyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 6 * kcomp + i); 

            auto tg_x_xxxxyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 7 * kcomp + i); 

            auto tg_x_xxxxyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 8 * kcomp + i); 

            auto tg_x_xxxxzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 9 * kcomp + i); 

            auto tg_x_xxxyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 10 * kcomp + i); 

            auto tg_x_xxxyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 11 * kcomp + i); 

            auto tg_x_xxxyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 12 * kcomp + i); 

            auto tg_x_xxxyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 13 * kcomp + i); 

            auto tg_x_xxxzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 14 * kcomp + i); 

            auto tg_x_xxyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 15 * kcomp + i); 

            auto tg_x_xxyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 16 * kcomp + i); 

            auto tg_x_xxyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 17 * kcomp + i); 

            auto tg_x_xxyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 18 * kcomp + i); 

            auto tg_x_xxyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 19 * kcomp + i); 

            auto tg_x_xxzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 20 * kcomp + i); 

            auto tg_x_xyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 21 * kcomp + i); 

            auto tg_x_xyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 22 * kcomp + i); 

            auto tg_x_xyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 23 * kcomp + i); 

            auto tg_x_xyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 24 * kcomp + i); 

            auto tg_x_xyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 25 * kcomp + i); 

            auto tg_x_xyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 26 * kcomp + i); 

            auto tg_x_xzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 27 * kcomp + i); 

            auto tg_y_xxxxxxx_0 = braBuffer.data(pidx_g_1_7_m0 + 36 * kcomp + i); 

            auto tg_y_xxxxxxy_0 = braBuffer.data(pidx_g_1_7_m0 + 37 * kcomp + i); 

            auto tg_y_xxxxxxz_0 = braBuffer.data(pidx_g_1_7_m0 + 38 * kcomp + i); 

            auto tg_y_xxxxxyy_0 = braBuffer.data(pidx_g_1_7_m0 + 39 * kcomp + i); 

            auto tg_y_xxxxxyz_0 = braBuffer.data(pidx_g_1_7_m0 + 40 * kcomp + i); 

            auto tg_y_xxxxxzz_0 = braBuffer.data(pidx_g_1_7_m0 + 41 * kcomp + i); 

            auto tg_y_xxxxyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 42 * kcomp + i); 

            auto tg_y_xxxxyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 43 * kcomp + i); 

            auto tg_y_xxxxyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 44 * kcomp + i); 

            auto tg_y_xxxxzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 45 * kcomp + i); 

            auto tg_y_xxxyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 46 * kcomp + i); 

            auto tg_y_xxxyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 47 * kcomp + i); 

            auto tg_y_xxxyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 48 * kcomp + i); 

            auto tg_y_xxxyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 49 * kcomp + i); 

            auto tg_y_xxxzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 50 * kcomp + i); 

            auto tg_y_xxyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 51 * kcomp + i); 

            auto tg_y_xxyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 52 * kcomp + i); 

            auto tg_y_xxyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 53 * kcomp + i); 

            auto tg_y_xxyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 54 * kcomp + i); 

            auto tg_y_xxyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 55 * kcomp + i); 

            auto tg_y_xxzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 56 * kcomp + i); 

            auto tg_y_xyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 57 * kcomp + i); 

            auto tg_y_xyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 58 * kcomp + i); 

            auto tg_y_xyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 59 * kcomp + i); 

            auto tg_y_xyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 60 * kcomp + i); 

            auto tg_y_xyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 61 * kcomp + i); 

            auto tg_y_xyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 62 * kcomp + i); 

            auto tg_y_xzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 63 * kcomp + i); 

            auto tg_z_xxxxxxx_0 = braBuffer.data(pidx_g_1_7_m0 + 72 * kcomp + i); 

            auto tg_z_xxxxxxy_0 = braBuffer.data(pidx_g_1_7_m0 + 73 * kcomp + i); 

            auto tg_z_xxxxxxz_0 = braBuffer.data(pidx_g_1_7_m0 + 74 * kcomp + i); 

            auto tg_z_xxxxxyy_0 = braBuffer.data(pidx_g_1_7_m0 + 75 * kcomp + i); 

            auto tg_z_xxxxxyz_0 = braBuffer.data(pidx_g_1_7_m0 + 76 * kcomp + i); 

            auto tg_z_xxxxxzz_0 = braBuffer.data(pidx_g_1_7_m0 + 77 * kcomp + i); 

            auto tg_z_xxxxyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 78 * kcomp + i); 

            auto tg_z_xxxxyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 79 * kcomp + i); 

            auto tg_z_xxxxyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 80 * kcomp + i); 

            auto tg_z_xxxxzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 81 * kcomp + i); 

            auto tg_z_xxxyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 82 * kcomp + i); 

            auto tg_z_xxxyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 83 * kcomp + i); 

            auto tg_z_xxxyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 84 * kcomp + i); 

            auto tg_z_xxxyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 85 * kcomp + i); 

            auto tg_z_xxxzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 86 * kcomp + i); 

            auto tg_z_xxyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 87 * kcomp + i); 

            auto tg_z_xxyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 88 * kcomp + i); 

            auto tg_z_xxyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 89 * kcomp + i); 

            auto tg_z_xxyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 90 * kcomp + i); 

            auto tg_z_xxyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 91 * kcomp + i); 

            auto tg_z_xxzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 92 * kcomp + i); 

            auto tg_z_xyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 93 * kcomp + i); 

            auto tg_z_xyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 94 * kcomp + i); 

            auto tg_z_xyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 95 * kcomp + i); 

            auto tg_z_xyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 96 * kcomp + i); 

            auto tg_z_xyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 97 * kcomp + i); 

            auto tg_z_xyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 98 * kcomp + i); 

            auto tg_z_xzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 99 * kcomp + i); 

            // set up pointers to integrals

            auto tg_xx_xxxxxx_0 = braBuffer.data(pidx_g_2_6_m0 + i); 

            auto tg_xx_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + kcomp + i); 

            auto tg_xx_xxxxxz_0 = braBuffer.data(pidx_g_2_6_m0 + 2 * kcomp + i); 

            auto tg_xx_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 3 * kcomp + i); 

            auto tg_xx_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 4 * kcomp + i); 

            auto tg_xx_xxxxzz_0 = braBuffer.data(pidx_g_2_6_m0 + 5 * kcomp + i); 

            auto tg_xx_xxxyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 6 * kcomp + i); 

            auto tg_xx_xxxyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 7 * kcomp + i); 

            auto tg_xx_xxxyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 8 * kcomp + i); 

            auto tg_xx_xxxzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 9 * kcomp + i); 

            auto tg_xx_xxyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 10 * kcomp + i); 

            auto tg_xx_xxyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 11 * kcomp + i); 

            auto tg_xx_xxyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 12 * kcomp + i); 

            auto tg_xx_xxyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 13 * kcomp + i); 

            auto tg_xx_xxzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 14 * kcomp + i); 

            auto tg_xx_xyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 15 * kcomp + i); 

            auto tg_xx_xyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 16 * kcomp + i); 

            auto tg_xx_xyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 17 * kcomp + i); 

            auto tg_xx_xyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 18 * kcomp + i); 

            auto tg_xx_xyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 19 * kcomp + i); 

            auto tg_xx_xzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 20 * kcomp + i); 

            auto tg_xx_yyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 21 * kcomp + i); 

            auto tg_xx_yyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 22 * kcomp + i); 

            auto tg_xx_yyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 23 * kcomp + i); 

            auto tg_xx_yyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 24 * kcomp + i); 

            auto tg_xx_yyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 25 * kcomp + i); 

            auto tg_xx_yzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 26 * kcomp + i); 

            auto tg_xx_zzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 27 * kcomp + i); 

            auto tg_xy_xxxxxx_0 = braBuffer.data(pidx_g_2_6_m0 + 28 * kcomp + i); 

            auto tg_xy_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + 29 * kcomp + i); 

            auto tg_xy_xxxxxz_0 = braBuffer.data(pidx_g_2_6_m0 + 30 * kcomp + i); 

            auto tg_xy_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 31 * kcomp + i); 

            auto tg_xy_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 32 * kcomp + i); 

            auto tg_xy_xxxxzz_0 = braBuffer.data(pidx_g_2_6_m0 + 33 * kcomp + i); 

            auto tg_xy_xxxyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 34 * kcomp + i); 

            auto tg_xy_xxxyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 35 * kcomp + i); 

            auto tg_xy_xxxyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 36 * kcomp + i); 

            auto tg_xy_xxxzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 37 * kcomp + i); 

            auto tg_xy_xxyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 38 * kcomp + i); 

            auto tg_xy_xxyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 39 * kcomp + i); 

            auto tg_xy_xxyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 40 * kcomp + i); 

            auto tg_xy_xxyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 41 * kcomp + i); 

            auto tg_xy_xxzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 42 * kcomp + i); 

            auto tg_xy_xyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 43 * kcomp + i); 

            auto tg_xy_xyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 44 * kcomp + i); 

            auto tg_xy_xyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 45 * kcomp + i); 

            auto tg_xy_xyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 46 * kcomp + i); 

            auto tg_xy_xyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 47 * kcomp + i); 

            auto tg_xy_xzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 48 * kcomp + i); 

            auto tg_xy_yyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 49 * kcomp + i); 

            auto tg_xy_yyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 50 * kcomp + i); 

            auto tg_xy_yyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 51 * kcomp + i); 

            auto tg_xy_yyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 52 * kcomp + i); 

            auto tg_xy_yyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 53 * kcomp + i); 

            auto tg_xy_yzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 54 * kcomp + i); 

            auto tg_xy_zzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 55 * kcomp + i); 

            auto tg_xz_xxxxxx_0 = braBuffer.data(pidx_g_2_6_m0 + 56 * kcomp + i); 

            auto tg_xz_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + 57 * kcomp + i); 

            auto tg_xz_xxxxxz_0 = braBuffer.data(pidx_g_2_6_m0 + 58 * kcomp + i); 

            auto tg_xz_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 59 * kcomp + i); 

            auto tg_xz_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 60 * kcomp + i); 

            auto tg_xz_xxxxzz_0 = braBuffer.data(pidx_g_2_6_m0 + 61 * kcomp + i); 

            auto tg_xz_xxxyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 62 * kcomp + i); 

            auto tg_xz_xxxyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 63 * kcomp + i); 

            auto tg_xz_xxxyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 64 * kcomp + i); 

            auto tg_xz_xxxzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 65 * kcomp + i); 

            auto tg_xz_xxyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 66 * kcomp + i); 

            auto tg_xz_xxyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 67 * kcomp + i); 

            auto tg_xz_xxyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 68 * kcomp + i); 

            auto tg_xz_xxyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 69 * kcomp + i); 

            auto tg_xz_xxzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 70 * kcomp + i); 

            auto tg_xz_xyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 71 * kcomp + i); 

            auto tg_xz_xyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 72 * kcomp + i); 

            auto tg_xz_xyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 73 * kcomp + i); 

            auto tg_xz_xyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 74 * kcomp + i); 

            auto tg_xz_xyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 75 * kcomp + i); 

            auto tg_xz_xzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 76 * kcomp + i); 

            auto tg_xz_yyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 77 * kcomp + i); 

            auto tg_xz_yyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 78 * kcomp + i); 

            auto tg_xz_yyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 79 * kcomp + i); 

            auto tg_xz_yyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 80 * kcomp + i); 

            auto tg_xz_yyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 81 * kcomp + i); 

            auto tg_xz_yzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 82 * kcomp + i); 

            auto tg_xz_zzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 83 * kcomp + i); 

            // Batch of Integrals (0,84)

            #pragma omp simd aligned(tg_x_xxxxxx_0, tg_x_xxxxxxx_0, tg_x_xxxxxxy_0, tg_x_xxxxxxz_0, \
                                     tg_x_xxxxxy_0, tg_x_xxxxxyy_0, tg_x_xxxxxyz_0, tg_x_xxxxxz_0, tg_x_xxxxxzz_0, \
                                     tg_x_xxxxyy_0, tg_x_xxxxyyy_0, tg_x_xxxxyyz_0, tg_x_xxxxyz_0, tg_x_xxxxyzz_0, \
                                     tg_x_xxxxzz_0, tg_x_xxxxzzz_0, tg_x_xxxyyy_0, tg_x_xxxyyyy_0, tg_x_xxxyyyz_0, \
                                     tg_x_xxxyyz_0, tg_x_xxxyyzz_0, tg_x_xxxyzz_0, tg_x_xxxyzzz_0, tg_x_xxxzzz_0, \
                                     tg_x_xxxzzzz_0, tg_x_xxyyyy_0, tg_x_xxyyyyy_0, tg_x_xxyyyyz_0, tg_x_xxyyyz_0, \
                                     tg_x_xxyyyzz_0, tg_x_xxyyzz_0, tg_x_xxyyzzz_0, tg_x_xxyzzz_0, tg_x_xxyzzzz_0, \
                                     tg_x_xxzzzz_0, tg_x_xxzzzzz_0, tg_x_xyyyyy_0, tg_x_xyyyyyy_0, tg_x_xyyyyyz_0, \
                                     tg_x_xyyyyz_0, tg_x_xyyyyzz_0, tg_x_xyyyzz_0, tg_x_xyyyzzz_0, tg_x_xyyzzz_0, \
                                     tg_x_xyyzzzz_0, tg_x_xyzzzz_0, tg_x_xyzzzzz_0, tg_x_xzzzzz_0, tg_x_xzzzzzz_0, \
                                     tg_x_yyyyyy_0, tg_x_yyyyyz_0, tg_x_yyyyzz_0, tg_x_yyyzzz_0, tg_x_yyzzzz_0, \
                                     tg_x_yzzzzz_0, tg_x_zzzzzz_0, tg_xx_xxxxxx_0, tg_xx_xxxxxy_0, tg_xx_xxxxxz_0, \
                                     tg_xx_xxxxyy_0, tg_xx_xxxxyz_0, tg_xx_xxxxzz_0, tg_xx_xxxyyy_0, tg_xx_xxxyyz_0, \
                                     tg_xx_xxxyzz_0, tg_xx_xxxzzz_0, tg_xx_xxyyyy_0, tg_xx_xxyyyz_0, tg_xx_xxyyzz_0, \
                                     tg_xx_xxyzzz_0, tg_xx_xxzzzz_0, tg_xx_xyyyyy_0, tg_xx_xyyyyz_0, tg_xx_xyyyzz_0, \
                                     tg_xx_xyyzzz_0, tg_xx_xyzzzz_0, tg_xx_xzzzzz_0, tg_xx_yyyyyy_0, tg_xx_yyyyyz_0, \
                                     tg_xx_yyyyzz_0, tg_xx_yyyzzz_0, tg_xx_yyzzzz_0, tg_xx_yzzzzz_0, tg_xx_zzzzzz_0, \
                                     tg_xy_xxxxxx_0, tg_xy_xxxxxy_0, tg_xy_xxxxxz_0, tg_xy_xxxxyy_0, tg_xy_xxxxyz_0, \
                                     tg_xy_xxxxzz_0, tg_xy_xxxyyy_0, tg_xy_xxxyyz_0, tg_xy_xxxyzz_0, tg_xy_xxxzzz_0, \
                                     tg_xy_xxyyyy_0, tg_xy_xxyyyz_0, tg_xy_xxyyzz_0, tg_xy_xxyzzz_0, tg_xy_xxzzzz_0, \
                                     tg_xy_xyyyyy_0, tg_xy_xyyyyz_0, tg_xy_xyyyzz_0, tg_xy_xyyzzz_0, tg_xy_xyzzzz_0, \
                                     tg_xy_xzzzzz_0, tg_xy_yyyyyy_0, tg_xy_yyyyyz_0, tg_xy_yyyyzz_0, tg_xy_yyyzzz_0, \
                                     tg_xy_yyzzzz_0, tg_xy_yzzzzz_0, tg_xy_zzzzzz_0, tg_xz_xxxxxx_0, tg_xz_xxxxxy_0, \
                                     tg_xz_xxxxxz_0, tg_xz_xxxxyy_0, tg_xz_xxxxyz_0, tg_xz_xxxxzz_0, tg_xz_xxxyyy_0, \
                                     tg_xz_xxxyyz_0, tg_xz_xxxyzz_0, tg_xz_xxxzzz_0, tg_xz_xxyyyy_0, tg_xz_xxyyyz_0, \
                                     tg_xz_xxyyzz_0, tg_xz_xxyzzz_0, tg_xz_xxzzzz_0, tg_xz_xyyyyy_0, tg_xz_xyyyyz_0, \
                                     tg_xz_xyyyzz_0, tg_xz_xyyzzz_0, tg_xz_xyzzzz_0, tg_xz_xzzzzz_0, tg_xz_yyyyyy_0, \
                                     tg_xz_yyyyyz_0, tg_xz_yyyyzz_0, tg_xz_yyyzzz_0, tg_xz_yyzzzz_0, tg_xz_yzzzzz_0, \
                                     tg_xz_zzzzzz_0, tg_y_xxxxxx_0, tg_y_xxxxxxx_0, tg_y_xxxxxxy_0, tg_y_xxxxxxz_0, \
                                     tg_y_xxxxxy_0, tg_y_xxxxxyy_0, tg_y_xxxxxyz_0, tg_y_xxxxxz_0, tg_y_xxxxxzz_0, \
                                     tg_y_xxxxyy_0, tg_y_xxxxyyy_0, tg_y_xxxxyyz_0, tg_y_xxxxyz_0, tg_y_xxxxyzz_0, \
                                     tg_y_xxxxzz_0, tg_y_xxxxzzz_0, tg_y_xxxyyy_0, tg_y_xxxyyyy_0, tg_y_xxxyyyz_0, \
                                     tg_y_xxxyyz_0, tg_y_xxxyyzz_0, tg_y_xxxyzz_0, tg_y_xxxyzzz_0, tg_y_xxxzzz_0, \
                                     tg_y_xxxzzzz_0, tg_y_xxyyyy_0, tg_y_xxyyyyy_0, tg_y_xxyyyyz_0, tg_y_xxyyyz_0, \
                                     tg_y_xxyyyzz_0, tg_y_xxyyzz_0, tg_y_xxyyzzz_0, tg_y_xxyzzz_0, tg_y_xxyzzzz_0, \
                                     tg_y_xxzzzz_0, tg_y_xxzzzzz_0, tg_y_xyyyyy_0, tg_y_xyyyyyy_0, tg_y_xyyyyyz_0, \
                                     tg_y_xyyyyz_0, tg_y_xyyyyzz_0, tg_y_xyyyzz_0, tg_y_xyyyzzz_0, tg_y_xyyzzz_0, \
                                     tg_y_xyyzzzz_0, tg_y_xyzzzz_0, tg_y_xyzzzzz_0, tg_y_xzzzzz_0, tg_y_xzzzzzz_0, \
                                     tg_y_yyyyyy_0, tg_y_yyyyyz_0, tg_y_yyyyzz_0, tg_y_yyyzzz_0, tg_y_yyzzzz_0, \
                                     tg_y_yzzzzz_0, tg_y_zzzzzz_0, tg_z_xxxxxx_0, tg_z_xxxxxxx_0, tg_z_xxxxxxy_0, \
                                     tg_z_xxxxxxz_0, tg_z_xxxxxy_0, tg_z_xxxxxyy_0, tg_z_xxxxxyz_0, tg_z_xxxxxz_0, \
                                     tg_z_xxxxxzz_0, tg_z_xxxxyy_0, tg_z_xxxxyyy_0, tg_z_xxxxyyz_0, tg_z_xxxxyz_0, \
                                     tg_z_xxxxyzz_0, tg_z_xxxxzz_0, tg_z_xxxxzzz_0, tg_z_xxxyyy_0, tg_z_xxxyyyy_0, \
                                     tg_z_xxxyyyz_0, tg_z_xxxyyz_0, tg_z_xxxyyzz_0, tg_z_xxxyzz_0, tg_z_xxxyzzz_0, \
                                     tg_z_xxxzzz_0, tg_z_xxxzzzz_0, tg_z_xxyyyy_0, tg_z_xxyyyyy_0, tg_z_xxyyyyz_0, \
                                     tg_z_xxyyyz_0, tg_z_xxyyyzz_0, tg_z_xxyyzz_0, tg_z_xxyyzzz_0, tg_z_xxyzzz_0, \
                                     tg_z_xxyzzzz_0, tg_z_xxzzzz_0, tg_z_xxzzzzz_0, tg_z_xyyyyy_0, tg_z_xyyyyyy_0, \
                                     tg_z_xyyyyyz_0, tg_z_xyyyyz_0, tg_z_xyyyyzz_0, tg_z_xyyyzz_0, tg_z_xyyyzzz_0, \
                                     tg_z_xyyzzz_0, tg_z_xyyzzzz_0, tg_z_xyzzzz_0, tg_z_xyzzzzz_0, tg_z_xzzzzz_0, \
                                     tg_z_xzzzzzz_0, tg_z_yyyyyy_0, tg_z_yyyyyz_0, tg_z_yyyyzz_0, tg_z_yyyzzz_0, \
                                     tg_z_yyzzzz_0, tg_z_yzzzzz_0, tg_z_zzzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xx_xxxxxx_0[j] = -ab_x * tg_x_xxxxxx_0[j] + tg_x_xxxxxxx_0[j];

                tg_xx_xxxxxy_0[j] = -ab_x * tg_x_xxxxxy_0[j] + tg_x_xxxxxxy_0[j];

                tg_xx_xxxxxz_0[j] = -ab_x * tg_x_xxxxxz_0[j] + tg_x_xxxxxxz_0[j];

                tg_xx_xxxxyy_0[j] = -ab_x * tg_x_xxxxyy_0[j] + tg_x_xxxxxyy_0[j];

                tg_xx_xxxxyz_0[j] = -ab_x * tg_x_xxxxyz_0[j] + tg_x_xxxxxyz_0[j];

                tg_xx_xxxxzz_0[j] = -ab_x * tg_x_xxxxzz_0[j] + tg_x_xxxxxzz_0[j];

                tg_xx_xxxyyy_0[j] = -ab_x * tg_x_xxxyyy_0[j] + tg_x_xxxxyyy_0[j];

                tg_xx_xxxyyz_0[j] = -ab_x * tg_x_xxxyyz_0[j] + tg_x_xxxxyyz_0[j];

                tg_xx_xxxyzz_0[j] = -ab_x * tg_x_xxxyzz_0[j] + tg_x_xxxxyzz_0[j];

                tg_xx_xxxzzz_0[j] = -ab_x * tg_x_xxxzzz_0[j] + tg_x_xxxxzzz_0[j];

                tg_xx_xxyyyy_0[j] = -ab_x * tg_x_xxyyyy_0[j] + tg_x_xxxyyyy_0[j];

                tg_xx_xxyyyz_0[j] = -ab_x * tg_x_xxyyyz_0[j] + tg_x_xxxyyyz_0[j];

                tg_xx_xxyyzz_0[j] = -ab_x * tg_x_xxyyzz_0[j] + tg_x_xxxyyzz_0[j];

                tg_xx_xxyzzz_0[j] = -ab_x * tg_x_xxyzzz_0[j] + tg_x_xxxyzzz_0[j];

                tg_xx_xxzzzz_0[j] = -ab_x * tg_x_xxzzzz_0[j] + tg_x_xxxzzzz_0[j];

                tg_xx_xyyyyy_0[j] = -ab_x * tg_x_xyyyyy_0[j] + tg_x_xxyyyyy_0[j];

                tg_xx_xyyyyz_0[j] = -ab_x * tg_x_xyyyyz_0[j] + tg_x_xxyyyyz_0[j];

                tg_xx_xyyyzz_0[j] = -ab_x * tg_x_xyyyzz_0[j] + tg_x_xxyyyzz_0[j];

                tg_xx_xyyzzz_0[j] = -ab_x * tg_x_xyyzzz_0[j] + tg_x_xxyyzzz_0[j];

                tg_xx_xyzzzz_0[j] = -ab_x * tg_x_xyzzzz_0[j] + tg_x_xxyzzzz_0[j];

                tg_xx_xzzzzz_0[j] = -ab_x * tg_x_xzzzzz_0[j] + tg_x_xxzzzzz_0[j];

                tg_xx_yyyyyy_0[j] = -ab_x * tg_x_yyyyyy_0[j] + tg_x_xyyyyyy_0[j];

                tg_xx_yyyyyz_0[j] = -ab_x * tg_x_yyyyyz_0[j] + tg_x_xyyyyyz_0[j];

                tg_xx_yyyyzz_0[j] = -ab_x * tg_x_yyyyzz_0[j] + tg_x_xyyyyzz_0[j];

                tg_xx_yyyzzz_0[j] = -ab_x * tg_x_yyyzzz_0[j] + tg_x_xyyyzzz_0[j];

                tg_xx_yyzzzz_0[j] = -ab_x * tg_x_yyzzzz_0[j] + tg_x_xyyzzzz_0[j];

                tg_xx_yzzzzz_0[j] = -ab_x * tg_x_yzzzzz_0[j] + tg_x_xyzzzzz_0[j];

                tg_xx_zzzzzz_0[j] = -ab_x * tg_x_zzzzzz_0[j] + tg_x_xzzzzzz_0[j];

                tg_xy_xxxxxx_0[j] = -ab_x * tg_y_xxxxxx_0[j] + tg_y_xxxxxxx_0[j];

                tg_xy_xxxxxy_0[j] = -ab_x * tg_y_xxxxxy_0[j] + tg_y_xxxxxxy_0[j];

                tg_xy_xxxxxz_0[j] = -ab_x * tg_y_xxxxxz_0[j] + tg_y_xxxxxxz_0[j];

                tg_xy_xxxxyy_0[j] = -ab_x * tg_y_xxxxyy_0[j] + tg_y_xxxxxyy_0[j];

                tg_xy_xxxxyz_0[j] = -ab_x * tg_y_xxxxyz_0[j] + tg_y_xxxxxyz_0[j];

                tg_xy_xxxxzz_0[j] = -ab_x * tg_y_xxxxzz_0[j] + tg_y_xxxxxzz_0[j];

                tg_xy_xxxyyy_0[j] = -ab_x * tg_y_xxxyyy_0[j] + tg_y_xxxxyyy_0[j];

                tg_xy_xxxyyz_0[j] = -ab_x * tg_y_xxxyyz_0[j] + tg_y_xxxxyyz_0[j];

                tg_xy_xxxyzz_0[j] = -ab_x * tg_y_xxxyzz_0[j] + tg_y_xxxxyzz_0[j];

                tg_xy_xxxzzz_0[j] = -ab_x * tg_y_xxxzzz_0[j] + tg_y_xxxxzzz_0[j];

                tg_xy_xxyyyy_0[j] = -ab_x * tg_y_xxyyyy_0[j] + tg_y_xxxyyyy_0[j];

                tg_xy_xxyyyz_0[j] = -ab_x * tg_y_xxyyyz_0[j] + tg_y_xxxyyyz_0[j];

                tg_xy_xxyyzz_0[j] = -ab_x * tg_y_xxyyzz_0[j] + tg_y_xxxyyzz_0[j];

                tg_xy_xxyzzz_0[j] = -ab_x * tg_y_xxyzzz_0[j] + tg_y_xxxyzzz_0[j];

                tg_xy_xxzzzz_0[j] = -ab_x * tg_y_xxzzzz_0[j] + tg_y_xxxzzzz_0[j];

                tg_xy_xyyyyy_0[j] = -ab_x * tg_y_xyyyyy_0[j] + tg_y_xxyyyyy_0[j];

                tg_xy_xyyyyz_0[j] = -ab_x * tg_y_xyyyyz_0[j] + tg_y_xxyyyyz_0[j];

                tg_xy_xyyyzz_0[j] = -ab_x * tg_y_xyyyzz_0[j] + tg_y_xxyyyzz_0[j];

                tg_xy_xyyzzz_0[j] = -ab_x * tg_y_xyyzzz_0[j] + tg_y_xxyyzzz_0[j];

                tg_xy_xyzzzz_0[j] = -ab_x * tg_y_xyzzzz_0[j] + tg_y_xxyzzzz_0[j];

                tg_xy_xzzzzz_0[j] = -ab_x * tg_y_xzzzzz_0[j] + tg_y_xxzzzzz_0[j];

                tg_xy_yyyyyy_0[j] = -ab_x * tg_y_yyyyyy_0[j] + tg_y_xyyyyyy_0[j];

                tg_xy_yyyyyz_0[j] = -ab_x * tg_y_yyyyyz_0[j] + tg_y_xyyyyyz_0[j];

                tg_xy_yyyyzz_0[j] = -ab_x * tg_y_yyyyzz_0[j] + tg_y_xyyyyzz_0[j];

                tg_xy_yyyzzz_0[j] = -ab_x * tg_y_yyyzzz_0[j] + tg_y_xyyyzzz_0[j];

                tg_xy_yyzzzz_0[j] = -ab_x * tg_y_yyzzzz_0[j] + tg_y_xyyzzzz_0[j];

                tg_xy_yzzzzz_0[j] = -ab_x * tg_y_yzzzzz_0[j] + tg_y_xyzzzzz_0[j];

                tg_xy_zzzzzz_0[j] = -ab_x * tg_y_zzzzzz_0[j] + tg_y_xzzzzzz_0[j];

                tg_xz_xxxxxx_0[j] = -ab_x * tg_z_xxxxxx_0[j] + tg_z_xxxxxxx_0[j];

                tg_xz_xxxxxy_0[j] = -ab_x * tg_z_xxxxxy_0[j] + tg_z_xxxxxxy_0[j];

                tg_xz_xxxxxz_0[j] = -ab_x * tg_z_xxxxxz_0[j] + tg_z_xxxxxxz_0[j];

                tg_xz_xxxxyy_0[j] = -ab_x * tg_z_xxxxyy_0[j] + tg_z_xxxxxyy_0[j];

                tg_xz_xxxxyz_0[j] = -ab_x * tg_z_xxxxyz_0[j] + tg_z_xxxxxyz_0[j];

                tg_xz_xxxxzz_0[j] = -ab_x * tg_z_xxxxzz_0[j] + tg_z_xxxxxzz_0[j];

                tg_xz_xxxyyy_0[j] = -ab_x * tg_z_xxxyyy_0[j] + tg_z_xxxxyyy_0[j];

                tg_xz_xxxyyz_0[j] = -ab_x * tg_z_xxxyyz_0[j] + tg_z_xxxxyyz_0[j];

                tg_xz_xxxyzz_0[j] = -ab_x * tg_z_xxxyzz_0[j] + tg_z_xxxxyzz_0[j];

                tg_xz_xxxzzz_0[j] = -ab_x * tg_z_xxxzzz_0[j] + tg_z_xxxxzzz_0[j];

                tg_xz_xxyyyy_0[j] = -ab_x * tg_z_xxyyyy_0[j] + tg_z_xxxyyyy_0[j];

                tg_xz_xxyyyz_0[j] = -ab_x * tg_z_xxyyyz_0[j] + tg_z_xxxyyyz_0[j];

                tg_xz_xxyyzz_0[j] = -ab_x * tg_z_xxyyzz_0[j] + tg_z_xxxyyzz_0[j];

                tg_xz_xxyzzz_0[j] = -ab_x * tg_z_xxyzzz_0[j] + tg_z_xxxyzzz_0[j];

                tg_xz_xxzzzz_0[j] = -ab_x * tg_z_xxzzzz_0[j] + tg_z_xxxzzzz_0[j];

                tg_xz_xyyyyy_0[j] = -ab_x * tg_z_xyyyyy_0[j] + tg_z_xxyyyyy_0[j];

                tg_xz_xyyyyz_0[j] = -ab_x * tg_z_xyyyyz_0[j] + tg_z_xxyyyyz_0[j];

                tg_xz_xyyyzz_0[j] = -ab_x * tg_z_xyyyzz_0[j] + tg_z_xxyyyzz_0[j];

                tg_xz_xyyzzz_0[j] = -ab_x * tg_z_xyyzzz_0[j] + tg_z_xxyyzzz_0[j];

                tg_xz_xyzzzz_0[j] = -ab_x * tg_z_xyzzzz_0[j] + tg_z_xxyzzzz_0[j];

                tg_xz_xzzzzz_0[j] = -ab_x * tg_z_xzzzzz_0[j] + tg_z_xxzzzzz_0[j];

                tg_xz_yyyyyy_0[j] = -ab_x * tg_z_yyyyyy_0[j] + tg_z_xyyyyyy_0[j];

                tg_xz_yyyyyz_0[j] = -ab_x * tg_z_yyyyyz_0[j] + tg_z_xyyyyyz_0[j];

                tg_xz_yyyyzz_0[j] = -ab_x * tg_z_yyyyzz_0[j] + tg_z_xyyyyzz_0[j];

                tg_xz_yyyzzz_0[j] = -ab_x * tg_z_yyyzzz_0[j] + tg_z_xyyyzzz_0[j];

                tg_xz_yyzzzz_0[j] = -ab_x * tg_z_yyzzzz_0[j] + tg_z_xyyzzzz_0[j];

                tg_xz_yzzzzz_0[j] = -ab_x * tg_z_yzzzzz_0[j] + tg_z_xyzzzzz_0[j];

                tg_xz_zzzzzz_0[j] = -ab_x * tg_z_zzzzzz_0[j] + tg_z_xzzzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForDIXY_84_168(      CMemBlock2D<double>& braBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& abDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetContrPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (84,168)

        // set up distances R(AB) = A - B

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_2_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_2_6_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_1_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 7, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_y_xxxxxx_0 = braBuffer.data(pidx_g_1_6_m0 + 28 * kcomp + i); 

            auto tg_y_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + 29 * kcomp + i); 

            auto tg_y_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 30 * kcomp + i); 

            auto tg_y_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 31 * kcomp + i); 

            auto tg_y_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 32 * kcomp + i); 

            auto tg_y_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 33 * kcomp + i); 

            auto tg_y_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 34 * kcomp + i); 

            auto tg_y_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 35 * kcomp + i); 

            auto tg_y_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 36 * kcomp + i); 

            auto tg_y_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 37 * kcomp + i); 

            auto tg_y_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 38 * kcomp + i); 

            auto tg_y_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 39 * kcomp + i); 

            auto tg_y_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 40 * kcomp + i); 

            auto tg_y_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 41 * kcomp + i); 

            auto tg_y_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 42 * kcomp + i); 

            auto tg_y_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 43 * kcomp + i); 

            auto tg_y_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 44 * kcomp + i); 

            auto tg_y_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 45 * kcomp + i); 

            auto tg_y_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 46 * kcomp + i); 

            auto tg_y_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 47 * kcomp + i); 

            auto tg_y_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 48 * kcomp + i); 

            auto tg_y_yyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 49 * kcomp + i); 

            auto tg_y_yyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 50 * kcomp + i); 

            auto tg_y_yyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 51 * kcomp + i); 

            auto tg_y_yyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 52 * kcomp + i); 

            auto tg_y_yyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 53 * kcomp + i); 

            auto tg_y_yzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 54 * kcomp + i); 

            auto tg_y_zzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 55 * kcomp + i); 

            auto tg_z_xxxxxx_0 = braBuffer.data(pidx_g_1_6_m0 + 56 * kcomp + i); 

            auto tg_z_xxxxxy_0 = braBuffer.data(pidx_g_1_6_m0 + 57 * kcomp + i); 

            auto tg_z_xxxxxz_0 = braBuffer.data(pidx_g_1_6_m0 + 58 * kcomp + i); 

            auto tg_z_xxxxyy_0 = braBuffer.data(pidx_g_1_6_m0 + 59 * kcomp + i); 

            auto tg_z_xxxxyz_0 = braBuffer.data(pidx_g_1_6_m0 + 60 * kcomp + i); 

            auto tg_z_xxxxzz_0 = braBuffer.data(pidx_g_1_6_m0 + 61 * kcomp + i); 

            auto tg_z_xxxyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 62 * kcomp + i); 

            auto tg_z_xxxyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 63 * kcomp + i); 

            auto tg_z_xxxyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 64 * kcomp + i); 

            auto tg_z_xxxzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 65 * kcomp + i); 

            auto tg_z_xxyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 66 * kcomp + i); 

            auto tg_z_xxyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 67 * kcomp + i); 

            auto tg_z_xxyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 68 * kcomp + i); 

            auto tg_z_xxyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 69 * kcomp + i); 

            auto tg_z_xxzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 70 * kcomp + i); 

            auto tg_z_xyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 71 * kcomp + i); 

            auto tg_z_xyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 72 * kcomp + i); 

            auto tg_z_xyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 73 * kcomp + i); 

            auto tg_z_xyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 74 * kcomp + i); 

            auto tg_z_xyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 75 * kcomp + i); 

            auto tg_z_xzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 76 * kcomp + i); 

            auto tg_z_yyyyyy_0 = braBuffer.data(pidx_g_1_6_m0 + 77 * kcomp + i); 

            auto tg_z_yyyyyz_0 = braBuffer.data(pidx_g_1_6_m0 + 78 * kcomp + i); 

            auto tg_z_yyyyzz_0 = braBuffer.data(pidx_g_1_6_m0 + 79 * kcomp + i); 

            auto tg_z_yyyzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 80 * kcomp + i); 

            auto tg_z_yyzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 81 * kcomp + i); 

            auto tg_z_yzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 82 * kcomp + i); 

            auto tg_z_zzzzzz_0 = braBuffer.data(pidx_g_1_6_m0 + 83 * kcomp + i); 

            auto tg_y_xxxxxxy_0 = braBuffer.data(pidx_g_1_7_m0 + 37 * kcomp + i); 

            auto tg_y_xxxxxyy_0 = braBuffer.data(pidx_g_1_7_m0 + 39 * kcomp + i); 

            auto tg_y_xxxxxyz_0 = braBuffer.data(pidx_g_1_7_m0 + 40 * kcomp + i); 

            auto tg_y_xxxxyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 42 * kcomp + i); 

            auto tg_y_xxxxyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 43 * kcomp + i); 

            auto tg_y_xxxxyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 44 * kcomp + i); 

            auto tg_y_xxxyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 46 * kcomp + i); 

            auto tg_y_xxxyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 47 * kcomp + i); 

            auto tg_y_xxxyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 48 * kcomp + i); 

            auto tg_y_xxxyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 49 * kcomp + i); 

            auto tg_y_xxyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 51 * kcomp + i); 

            auto tg_y_xxyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 52 * kcomp + i); 

            auto tg_y_xxyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 53 * kcomp + i); 

            auto tg_y_xxyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 54 * kcomp + i); 

            auto tg_y_xxyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 55 * kcomp + i); 

            auto tg_y_xyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 57 * kcomp + i); 

            auto tg_y_xyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 58 * kcomp + i); 

            auto tg_y_xyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 59 * kcomp + i); 

            auto tg_y_xyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 60 * kcomp + i); 

            auto tg_y_xyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 61 * kcomp + i); 

            auto tg_y_xyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 62 * kcomp + i); 

            auto tg_y_yyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 64 * kcomp + i); 

            auto tg_y_yyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 65 * kcomp + i); 

            auto tg_y_yyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 66 * kcomp + i); 

            auto tg_y_yyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 67 * kcomp + i); 

            auto tg_y_yyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 68 * kcomp + i); 

            auto tg_y_yyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 69 * kcomp + i); 

            auto tg_y_yzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 70 * kcomp + i); 

            auto tg_z_xxxxxxy_0 = braBuffer.data(pidx_g_1_7_m0 + 73 * kcomp + i); 

            auto tg_z_xxxxxxz_0 = braBuffer.data(pidx_g_1_7_m0 + 74 * kcomp + i); 

            auto tg_z_xxxxxyy_0 = braBuffer.data(pidx_g_1_7_m0 + 75 * kcomp + i); 

            auto tg_z_xxxxxyz_0 = braBuffer.data(pidx_g_1_7_m0 + 76 * kcomp + i); 

            auto tg_z_xxxxxzz_0 = braBuffer.data(pidx_g_1_7_m0 + 77 * kcomp + i); 

            auto tg_z_xxxxyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 78 * kcomp + i); 

            auto tg_z_xxxxyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 79 * kcomp + i); 

            auto tg_z_xxxxyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 80 * kcomp + i); 

            auto tg_z_xxxxzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 81 * kcomp + i); 

            auto tg_z_xxxyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 82 * kcomp + i); 

            auto tg_z_xxxyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 83 * kcomp + i); 

            auto tg_z_xxxyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 84 * kcomp + i); 

            auto tg_z_xxxyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 85 * kcomp + i); 

            auto tg_z_xxxzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 86 * kcomp + i); 

            auto tg_z_xxyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 87 * kcomp + i); 

            auto tg_z_xxyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 88 * kcomp + i); 

            auto tg_z_xxyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 89 * kcomp + i); 

            auto tg_z_xxyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 90 * kcomp + i); 

            auto tg_z_xxyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 91 * kcomp + i); 

            auto tg_z_xxzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 92 * kcomp + i); 

            auto tg_z_xyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 93 * kcomp + i); 

            auto tg_z_xyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 94 * kcomp + i); 

            auto tg_z_xyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 95 * kcomp + i); 

            auto tg_z_xyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 96 * kcomp + i); 

            auto tg_z_xyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 97 * kcomp + i); 

            auto tg_z_xyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 98 * kcomp + i); 

            auto tg_z_xzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 99 * kcomp + i); 

            auto tg_z_yyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 100 * kcomp + i); 

            auto tg_z_yyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 101 * kcomp + i); 

            auto tg_z_yyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 102 * kcomp + i); 

            auto tg_z_yyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 103 * kcomp + i); 

            auto tg_z_yyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 104 * kcomp + i); 

            auto tg_z_yyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 105 * kcomp + i); 

            auto tg_z_yzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 106 * kcomp + i); 

            auto tg_z_zzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 107 * kcomp + i); 

            // set up pointers to integrals

            auto tg_yy_xxxxxx_0 = braBuffer.data(pidx_g_2_6_m0 + 84 * kcomp + i); 

            auto tg_yy_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + 85 * kcomp + i); 

            auto tg_yy_xxxxxz_0 = braBuffer.data(pidx_g_2_6_m0 + 86 * kcomp + i); 

            auto tg_yy_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 87 * kcomp + i); 

            auto tg_yy_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 88 * kcomp + i); 

            auto tg_yy_xxxxzz_0 = braBuffer.data(pidx_g_2_6_m0 + 89 * kcomp + i); 

            auto tg_yy_xxxyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 90 * kcomp + i); 

            auto tg_yy_xxxyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 91 * kcomp + i); 

            auto tg_yy_xxxyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 92 * kcomp + i); 

            auto tg_yy_xxxzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 93 * kcomp + i); 

            auto tg_yy_xxyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 94 * kcomp + i); 

            auto tg_yy_xxyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 95 * kcomp + i); 

            auto tg_yy_xxyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 96 * kcomp + i); 

            auto tg_yy_xxyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 97 * kcomp + i); 

            auto tg_yy_xxzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 98 * kcomp + i); 

            auto tg_yy_xyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 99 * kcomp + i); 

            auto tg_yy_xyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 100 * kcomp + i); 

            auto tg_yy_xyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 101 * kcomp + i); 

            auto tg_yy_xyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 102 * kcomp + i); 

            auto tg_yy_xyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 103 * kcomp + i); 

            auto tg_yy_xzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 104 * kcomp + i); 

            auto tg_yy_yyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 105 * kcomp + i); 

            auto tg_yy_yyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 106 * kcomp + i); 

            auto tg_yy_yyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 107 * kcomp + i); 

            auto tg_yy_yyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 108 * kcomp + i); 

            auto tg_yy_yyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 109 * kcomp + i); 

            auto tg_yy_yzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 110 * kcomp + i); 

            auto tg_yy_zzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 111 * kcomp + i); 

            auto tg_yz_xxxxxx_0 = braBuffer.data(pidx_g_2_6_m0 + 112 * kcomp + i); 

            auto tg_yz_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + 113 * kcomp + i); 

            auto tg_yz_xxxxxz_0 = braBuffer.data(pidx_g_2_6_m0 + 114 * kcomp + i); 

            auto tg_yz_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 115 * kcomp + i); 

            auto tg_yz_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 116 * kcomp + i); 

            auto tg_yz_xxxxzz_0 = braBuffer.data(pidx_g_2_6_m0 + 117 * kcomp + i); 

            auto tg_yz_xxxyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 118 * kcomp + i); 

            auto tg_yz_xxxyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 119 * kcomp + i); 

            auto tg_yz_xxxyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 120 * kcomp + i); 

            auto tg_yz_xxxzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 121 * kcomp + i); 

            auto tg_yz_xxyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 122 * kcomp + i); 

            auto tg_yz_xxyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 123 * kcomp + i); 

            auto tg_yz_xxyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 124 * kcomp + i); 

            auto tg_yz_xxyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 125 * kcomp + i); 

            auto tg_yz_xxzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 126 * kcomp + i); 

            auto tg_yz_xyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 127 * kcomp + i); 

            auto tg_yz_xyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 128 * kcomp + i); 

            auto tg_yz_xyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 129 * kcomp + i); 

            auto tg_yz_xyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 130 * kcomp + i); 

            auto tg_yz_xyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 131 * kcomp + i); 

            auto tg_yz_xzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 132 * kcomp + i); 

            auto tg_yz_yyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 133 * kcomp + i); 

            auto tg_yz_yyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 134 * kcomp + i); 

            auto tg_yz_yyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 135 * kcomp + i); 

            auto tg_yz_yyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 136 * kcomp + i); 

            auto tg_yz_yyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 137 * kcomp + i); 

            auto tg_yz_yzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 138 * kcomp + i); 

            auto tg_yz_zzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 139 * kcomp + i); 

            auto tg_zz_xxxxxx_0 = braBuffer.data(pidx_g_2_6_m0 + 140 * kcomp + i); 

            auto tg_zz_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + 141 * kcomp + i); 

            auto tg_zz_xxxxxz_0 = braBuffer.data(pidx_g_2_6_m0 + 142 * kcomp + i); 

            auto tg_zz_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 143 * kcomp + i); 

            auto tg_zz_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 144 * kcomp + i); 

            auto tg_zz_xxxxzz_0 = braBuffer.data(pidx_g_2_6_m0 + 145 * kcomp + i); 

            auto tg_zz_xxxyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 146 * kcomp + i); 

            auto tg_zz_xxxyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 147 * kcomp + i); 

            auto tg_zz_xxxyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 148 * kcomp + i); 

            auto tg_zz_xxxzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 149 * kcomp + i); 

            auto tg_zz_xxyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 150 * kcomp + i); 

            auto tg_zz_xxyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 151 * kcomp + i); 

            auto tg_zz_xxyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 152 * kcomp + i); 

            auto tg_zz_xxyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 153 * kcomp + i); 

            auto tg_zz_xxzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 154 * kcomp + i); 

            auto tg_zz_xyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 155 * kcomp + i); 

            auto tg_zz_xyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 156 * kcomp + i); 

            auto tg_zz_xyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 157 * kcomp + i); 

            auto tg_zz_xyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 158 * kcomp + i); 

            auto tg_zz_xyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 159 * kcomp + i); 

            auto tg_zz_xzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 160 * kcomp + i); 

            auto tg_zz_yyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 161 * kcomp + i); 

            auto tg_zz_yyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 162 * kcomp + i); 

            auto tg_zz_yyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 163 * kcomp + i); 

            auto tg_zz_yyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 164 * kcomp + i); 

            auto tg_zz_yyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 165 * kcomp + i); 

            auto tg_zz_yzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 166 * kcomp + i); 

            auto tg_zz_zzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 167 * kcomp + i); 

            // Batch of Integrals (84,168)

            #pragma omp simd aligned(tg_y_xxxxxx_0, tg_y_xxxxxxy_0, tg_y_xxxxxy_0, tg_y_xxxxxyy_0, \
                                     tg_y_xxxxxyz_0, tg_y_xxxxxz_0, tg_y_xxxxyy_0, tg_y_xxxxyyy_0, tg_y_xxxxyyz_0, \
                                     tg_y_xxxxyz_0, tg_y_xxxxyzz_0, tg_y_xxxxzz_0, tg_y_xxxyyy_0, tg_y_xxxyyyy_0, \
                                     tg_y_xxxyyyz_0, tg_y_xxxyyz_0, tg_y_xxxyyzz_0, tg_y_xxxyzz_0, tg_y_xxxyzzz_0, \
                                     tg_y_xxxzzz_0, tg_y_xxyyyy_0, tg_y_xxyyyyy_0, tg_y_xxyyyyz_0, tg_y_xxyyyz_0, \
                                     tg_y_xxyyyzz_0, tg_y_xxyyzz_0, tg_y_xxyyzzz_0, tg_y_xxyzzz_0, tg_y_xxyzzzz_0, \
                                     tg_y_xxzzzz_0, tg_y_xyyyyy_0, tg_y_xyyyyyy_0, tg_y_xyyyyyz_0, tg_y_xyyyyz_0, \
                                     tg_y_xyyyyzz_0, tg_y_xyyyzz_0, tg_y_xyyyzzz_0, tg_y_xyyzzz_0, tg_y_xyyzzzz_0, \
                                     tg_y_xyzzzz_0, tg_y_xyzzzzz_0, tg_y_xzzzzz_0, tg_y_yyyyyy_0, tg_y_yyyyyyy_0, \
                                     tg_y_yyyyyyz_0, tg_y_yyyyyz_0, tg_y_yyyyyzz_0, tg_y_yyyyzz_0, tg_y_yyyyzzz_0, \
                                     tg_y_yyyzzz_0, tg_y_yyyzzzz_0, tg_y_yyzzzz_0, tg_y_yyzzzzz_0, tg_y_yzzzzz_0, \
                                     tg_y_yzzzzzz_0, tg_y_zzzzzz_0, tg_yy_xxxxxx_0, tg_yy_xxxxxy_0, tg_yy_xxxxxz_0, \
                                     tg_yy_xxxxyy_0, tg_yy_xxxxyz_0, tg_yy_xxxxzz_0, tg_yy_xxxyyy_0, tg_yy_xxxyyz_0, \
                                     tg_yy_xxxyzz_0, tg_yy_xxxzzz_0, tg_yy_xxyyyy_0, tg_yy_xxyyyz_0, tg_yy_xxyyzz_0, \
                                     tg_yy_xxyzzz_0, tg_yy_xxzzzz_0, tg_yy_xyyyyy_0, tg_yy_xyyyyz_0, tg_yy_xyyyzz_0, \
                                     tg_yy_xyyzzz_0, tg_yy_xyzzzz_0, tg_yy_xzzzzz_0, tg_yy_yyyyyy_0, tg_yy_yyyyyz_0, \
                                     tg_yy_yyyyzz_0, tg_yy_yyyzzz_0, tg_yy_yyzzzz_0, tg_yy_yzzzzz_0, tg_yy_zzzzzz_0, \
                                     tg_yz_xxxxxx_0, tg_yz_xxxxxy_0, tg_yz_xxxxxz_0, tg_yz_xxxxyy_0, tg_yz_xxxxyz_0, \
                                     tg_yz_xxxxzz_0, tg_yz_xxxyyy_0, tg_yz_xxxyyz_0, tg_yz_xxxyzz_0, tg_yz_xxxzzz_0, \
                                     tg_yz_xxyyyy_0, tg_yz_xxyyyz_0, tg_yz_xxyyzz_0, tg_yz_xxyzzz_0, tg_yz_xxzzzz_0, \
                                     tg_yz_xyyyyy_0, tg_yz_xyyyyz_0, tg_yz_xyyyzz_0, tg_yz_xyyzzz_0, tg_yz_xyzzzz_0, \
                                     tg_yz_xzzzzz_0, tg_yz_yyyyyy_0, tg_yz_yyyyyz_0, tg_yz_yyyyzz_0, tg_yz_yyyzzz_0, \
                                     tg_yz_yyzzzz_0, tg_yz_yzzzzz_0, tg_yz_zzzzzz_0, tg_z_xxxxxx_0, tg_z_xxxxxxy_0, \
                                     tg_z_xxxxxxz_0, tg_z_xxxxxy_0, tg_z_xxxxxyy_0, tg_z_xxxxxyz_0, tg_z_xxxxxz_0, \
                                     tg_z_xxxxxzz_0, tg_z_xxxxyy_0, tg_z_xxxxyyy_0, tg_z_xxxxyyz_0, tg_z_xxxxyz_0, \
                                     tg_z_xxxxyzz_0, tg_z_xxxxzz_0, tg_z_xxxxzzz_0, tg_z_xxxyyy_0, tg_z_xxxyyyy_0, \
                                     tg_z_xxxyyyz_0, tg_z_xxxyyz_0, tg_z_xxxyyzz_0, tg_z_xxxyzz_0, tg_z_xxxyzzz_0, \
                                     tg_z_xxxzzz_0, tg_z_xxxzzzz_0, tg_z_xxyyyy_0, tg_z_xxyyyyy_0, tg_z_xxyyyyz_0, \
                                     tg_z_xxyyyz_0, tg_z_xxyyyzz_0, tg_z_xxyyzz_0, tg_z_xxyyzzz_0, tg_z_xxyzzz_0, \
                                     tg_z_xxyzzzz_0, tg_z_xxzzzz_0, tg_z_xxzzzzz_0, tg_z_xyyyyy_0, tg_z_xyyyyyy_0, \
                                     tg_z_xyyyyyz_0, tg_z_xyyyyz_0, tg_z_xyyyyzz_0, tg_z_xyyyzz_0, tg_z_xyyyzzz_0, \
                                     tg_z_xyyzzz_0, tg_z_xyyzzzz_0, tg_z_xyzzzz_0, tg_z_xyzzzzz_0, tg_z_xzzzzz_0, \
                                     tg_z_xzzzzzz_0, tg_z_yyyyyy_0, tg_z_yyyyyyy_0, tg_z_yyyyyyz_0, tg_z_yyyyyz_0, \
                                     tg_z_yyyyyzz_0, tg_z_yyyyzz_0, tg_z_yyyyzzz_0, tg_z_yyyzzz_0, tg_z_yyyzzzz_0, \
                                     tg_z_yyzzzz_0, tg_z_yyzzzzz_0, tg_z_yzzzzz_0, tg_z_yzzzzzz_0, tg_z_zzzzzz_0, \
                                     tg_z_zzzzzzz_0, tg_zz_xxxxxx_0, tg_zz_xxxxxy_0, tg_zz_xxxxxz_0, tg_zz_xxxxyy_0, \
                                     tg_zz_xxxxyz_0, tg_zz_xxxxzz_0, tg_zz_xxxyyy_0, tg_zz_xxxyyz_0, tg_zz_xxxyzz_0, \
                                     tg_zz_xxxzzz_0, tg_zz_xxyyyy_0, tg_zz_xxyyyz_0, tg_zz_xxyyzz_0, tg_zz_xxyzzz_0, \
                                     tg_zz_xxzzzz_0, tg_zz_xyyyyy_0, tg_zz_xyyyyz_0, tg_zz_xyyyzz_0, tg_zz_xyyzzz_0, \
                                     tg_zz_xyzzzz_0, tg_zz_xzzzzz_0, tg_zz_yyyyyy_0, tg_zz_yyyyyz_0, tg_zz_yyyyzz_0, \
                                     tg_zz_yyyzzz_0, tg_zz_yyzzzz_0, tg_zz_yzzzzz_0, tg_zz_zzzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_yy_xxxxxx_0[j] = -ab_y * tg_y_xxxxxx_0[j] + tg_y_xxxxxxy_0[j];

                tg_yy_xxxxxy_0[j] = -ab_y * tg_y_xxxxxy_0[j] + tg_y_xxxxxyy_0[j];

                tg_yy_xxxxxz_0[j] = -ab_y * tg_y_xxxxxz_0[j] + tg_y_xxxxxyz_0[j];

                tg_yy_xxxxyy_0[j] = -ab_y * tg_y_xxxxyy_0[j] + tg_y_xxxxyyy_0[j];

                tg_yy_xxxxyz_0[j] = -ab_y * tg_y_xxxxyz_0[j] + tg_y_xxxxyyz_0[j];

                tg_yy_xxxxzz_0[j] = -ab_y * tg_y_xxxxzz_0[j] + tg_y_xxxxyzz_0[j];

                tg_yy_xxxyyy_0[j] = -ab_y * tg_y_xxxyyy_0[j] + tg_y_xxxyyyy_0[j];

                tg_yy_xxxyyz_0[j] = -ab_y * tg_y_xxxyyz_0[j] + tg_y_xxxyyyz_0[j];

                tg_yy_xxxyzz_0[j] = -ab_y * tg_y_xxxyzz_0[j] + tg_y_xxxyyzz_0[j];

                tg_yy_xxxzzz_0[j] = -ab_y * tg_y_xxxzzz_0[j] + tg_y_xxxyzzz_0[j];

                tg_yy_xxyyyy_0[j] = -ab_y * tg_y_xxyyyy_0[j] + tg_y_xxyyyyy_0[j];

                tg_yy_xxyyyz_0[j] = -ab_y * tg_y_xxyyyz_0[j] + tg_y_xxyyyyz_0[j];

                tg_yy_xxyyzz_0[j] = -ab_y * tg_y_xxyyzz_0[j] + tg_y_xxyyyzz_0[j];

                tg_yy_xxyzzz_0[j] = -ab_y * tg_y_xxyzzz_0[j] + tg_y_xxyyzzz_0[j];

                tg_yy_xxzzzz_0[j] = -ab_y * tg_y_xxzzzz_0[j] + tg_y_xxyzzzz_0[j];

                tg_yy_xyyyyy_0[j] = -ab_y * tg_y_xyyyyy_0[j] + tg_y_xyyyyyy_0[j];

                tg_yy_xyyyyz_0[j] = -ab_y * tg_y_xyyyyz_0[j] + tg_y_xyyyyyz_0[j];

                tg_yy_xyyyzz_0[j] = -ab_y * tg_y_xyyyzz_0[j] + tg_y_xyyyyzz_0[j];

                tg_yy_xyyzzz_0[j] = -ab_y * tg_y_xyyzzz_0[j] + tg_y_xyyyzzz_0[j];

                tg_yy_xyzzzz_0[j] = -ab_y * tg_y_xyzzzz_0[j] + tg_y_xyyzzzz_0[j];

                tg_yy_xzzzzz_0[j] = -ab_y * tg_y_xzzzzz_0[j] + tg_y_xyzzzzz_0[j];

                tg_yy_yyyyyy_0[j] = -ab_y * tg_y_yyyyyy_0[j] + tg_y_yyyyyyy_0[j];

                tg_yy_yyyyyz_0[j] = -ab_y * tg_y_yyyyyz_0[j] + tg_y_yyyyyyz_0[j];

                tg_yy_yyyyzz_0[j] = -ab_y * tg_y_yyyyzz_0[j] + tg_y_yyyyyzz_0[j];

                tg_yy_yyyzzz_0[j] = -ab_y * tg_y_yyyzzz_0[j] + tg_y_yyyyzzz_0[j];

                tg_yy_yyzzzz_0[j] = -ab_y * tg_y_yyzzzz_0[j] + tg_y_yyyzzzz_0[j];

                tg_yy_yzzzzz_0[j] = -ab_y * tg_y_yzzzzz_0[j] + tg_y_yyzzzzz_0[j];

                tg_yy_zzzzzz_0[j] = -ab_y * tg_y_zzzzzz_0[j] + tg_y_yzzzzzz_0[j];

                tg_yz_xxxxxx_0[j] = -ab_y * tg_z_xxxxxx_0[j] + tg_z_xxxxxxy_0[j];

                tg_yz_xxxxxy_0[j] = -ab_y * tg_z_xxxxxy_0[j] + tg_z_xxxxxyy_0[j];

                tg_yz_xxxxxz_0[j] = -ab_y * tg_z_xxxxxz_0[j] + tg_z_xxxxxyz_0[j];

                tg_yz_xxxxyy_0[j] = -ab_y * tg_z_xxxxyy_0[j] + tg_z_xxxxyyy_0[j];

                tg_yz_xxxxyz_0[j] = -ab_y * tg_z_xxxxyz_0[j] + tg_z_xxxxyyz_0[j];

                tg_yz_xxxxzz_0[j] = -ab_y * tg_z_xxxxzz_0[j] + tg_z_xxxxyzz_0[j];

                tg_yz_xxxyyy_0[j] = -ab_y * tg_z_xxxyyy_0[j] + tg_z_xxxyyyy_0[j];

                tg_yz_xxxyyz_0[j] = -ab_y * tg_z_xxxyyz_0[j] + tg_z_xxxyyyz_0[j];

                tg_yz_xxxyzz_0[j] = -ab_y * tg_z_xxxyzz_0[j] + tg_z_xxxyyzz_0[j];

                tg_yz_xxxzzz_0[j] = -ab_y * tg_z_xxxzzz_0[j] + tg_z_xxxyzzz_0[j];

                tg_yz_xxyyyy_0[j] = -ab_y * tg_z_xxyyyy_0[j] + tg_z_xxyyyyy_0[j];

                tg_yz_xxyyyz_0[j] = -ab_y * tg_z_xxyyyz_0[j] + tg_z_xxyyyyz_0[j];

                tg_yz_xxyyzz_0[j] = -ab_y * tg_z_xxyyzz_0[j] + tg_z_xxyyyzz_0[j];

                tg_yz_xxyzzz_0[j] = -ab_y * tg_z_xxyzzz_0[j] + tg_z_xxyyzzz_0[j];

                tg_yz_xxzzzz_0[j] = -ab_y * tg_z_xxzzzz_0[j] + tg_z_xxyzzzz_0[j];

                tg_yz_xyyyyy_0[j] = -ab_y * tg_z_xyyyyy_0[j] + tg_z_xyyyyyy_0[j];

                tg_yz_xyyyyz_0[j] = -ab_y * tg_z_xyyyyz_0[j] + tg_z_xyyyyyz_0[j];

                tg_yz_xyyyzz_0[j] = -ab_y * tg_z_xyyyzz_0[j] + tg_z_xyyyyzz_0[j];

                tg_yz_xyyzzz_0[j] = -ab_y * tg_z_xyyzzz_0[j] + tg_z_xyyyzzz_0[j];

                tg_yz_xyzzzz_0[j] = -ab_y * tg_z_xyzzzz_0[j] + tg_z_xyyzzzz_0[j];

                tg_yz_xzzzzz_0[j] = -ab_y * tg_z_xzzzzz_0[j] + tg_z_xyzzzzz_0[j];

                tg_yz_yyyyyy_0[j] = -ab_y * tg_z_yyyyyy_0[j] + tg_z_yyyyyyy_0[j];

                tg_yz_yyyyyz_0[j] = -ab_y * tg_z_yyyyyz_0[j] + tg_z_yyyyyyz_0[j];

                tg_yz_yyyyzz_0[j] = -ab_y * tg_z_yyyyzz_0[j] + tg_z_yyyyyzz_0[j];

                tg_yz_yyyzzz_0[j] = -ab_y * tg_z_yyyzzz_0[j] + tg_z_yyyyzzz_0[j];

                tg_yz_yyzzzz_0[j] = -ab_y * tg_z_yyzzzz_0[j] + tg_z_yyyzzzz_0[j];

                tg_yz_yzzzzz_0[j] = -ab_y * tg_z_yzzzzz_0[j] + tg_z_yyzzzzz_0[j];

                tg_yz_zzzzzz_0[j] = -ab_y * tg_z_zzzzzz_0[j] + tg_z_yzzzzzz_0[j];

                tg_zz_xxxxxx_0[j] = -ab_z * tg_z_xxxxxx_0[j] + tg_z_xxxxxxz_0[j];

                tg_zz_xxxxxy_0[j] = -ab_z * tg_z_xxxxxy_0[j] + tg_z_xxxxxyz_0[j];

                tg_zz_xxxxxz_0[j] = -ab_z * tg_z_xxxxxz_0[j] + tg_z_xxxxxzz_0[j];

                tg_zz_xxxxyy_0[j] = -ab_z * tg_z_xxxxyy_0[j] + tg_z_xxxxyyz_0[j];

                tg_zz_xxxxyz_0[j] = -ab_z * tg_z_xxxxyz_0[j] + tg_z_xxxxyzz_0[j];

                tg_zz_xxxxzz_0[j] = -ab_z * tg_z_xxxxzz_0[j] + tg_z_xxxxzzz_0[j];

                tg_zz_xxxyyy_0[j] = -ab_z * tg_z_xxxyyy_0[j] + tg_z_xxxyyyz_0[j];

                tg_zz_xxxyyz_0[j] = -ab_z * tg_z_xxxyyz_0[j] + tg_z_xxxyyzz_0[j];

                tg_zz_xxxyzz_0[j] = -ab_z * tg_z_xxxyzz_0[j] + tg_z_xxxyzzz_0[j];

                tg_zz_xxxzzz_0[j] = -ab_z * tg_z_xxxzzz_0[j] + tg_z_xxxzzzz_0[j];

                tg_zz_xxyyyy_0[j] = -ab_z * tg_z_xxyyyy_0[j] + tg_z_xxyyyyz_0[j];

                tg_zz_xxyyyz_0[j] = -ab_z * tg_z_xxyyyz_0[j] + tg_z_xxyyyzz_0[j];

                tg_zz_xxyyzz_0[j] = -ab_z * tg_z_xxyyzz_0[j] + tg_z_xxyyzzz_0[j];

                tg_zz_xxyzzz_0[j] = -ab_z * tg_z_xxyzzz_0[j] + tg_z_xxyzzzz_0[j];

                tg_zz_xxzzzz_0[j] = -ab_z * tg_z_xxzzzz_0[j] + tg_z_xxzzzzz_0[j];

                tg_zz_xyyyyy_0[j] = -ab_z * tg_z_xyyyyy_0[j] + tg_z_xyyyyyz_0[j];

                tg_zz_xyyyyz_0[j] = -ab_z * tg_z_xyyyyz_0[j] + tg_z_xyyyyzz_0[j];

                tg_zz_xyyyzz_0[j] = -ab_z * tg_z_xyyyzz_0[j] + tg_z_xyyyzzz_0[j];

                tg_zz_xyyzzz_0[j] = -ab_z * tg_z_xyyzzz_0[j] + tg_z_xyyzzzz_0[j];

                tg_zz_xyzzzz_0[j] = -ab_z * tg_z_xyzzzz_0[j] + tg_z_xyzzzzz_0[j];

                tg_zz_xzzzzz_0[j] = -ab_z * tg_z_xzzzzz_0[j] + tg_z_xzzzzzz_0[j];

                tg_zz_yyyyyy_0[j] = -ab_z * tg_z_yyyyyy_0[j] + tg_z_yyyyyyz_0[j];

                tg_zz_yyyyyz_0[j] = -ab_z * tg_z_yyyyyz_0[j] + tg_z_yyyyyzz_0[j];

                tg_zz_yyyyzz_0[j] = -ab_z * tg_z_yyyyzz_0[j] + tg_z_yyyyzzz_0[j];

                tg_zz_yyyzzz_0[j] = -ab_z * tg_z_yyyzzz_0[j] + tg_z_yyyzzzz_0[j];

                tg_zz_yyzzzz_0[j] = -ab_z * tg_z_yyzzzz_0[j] + tg_z_yyzzzzz_0[j];

                tg_zz_yzzzzz_0[j] = -ab_z * tg_z_yzzzzz_0[j] + tg_z_yzzzzzz_0[j];

                tg_zz_zzzzzz_0[j] = -ab_z * tg_z_zzzzzz_0[j] + tg_z_zzzzzzz_0[j];
            }
        }
    }


} // eribrrfunc namespace

