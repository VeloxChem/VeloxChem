//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionBRRRecFuncForPXYY.hpp"

#include "AngularMomentum.hpp"

namespace eribrrfunc { // eribrrfunc namespace

    void
    compElectronRepulsionForPPXY(      CMemBlock2D<double>& braBuffer,
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

        auto pidx_g_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 1, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_1_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 1, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 2, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_0_x_0 = braBuffer.data(pidx_g_0_1_m0 + i); 

            auto tg_0_y_0 = braBuffer.data(pidx_g_0_1_m0 + kcomp + i); 

            auto tg_0_z_0 = braBuffer.data(pidx_g_0_1_m0 + 2 * kcomp + i); 

            auto tg_0_xx_0 = braBuffer.data(pidx_g_0_2_m0 + i); 

            auto tg_0_xy_0 = braBuffer.data(pidx_g_0_2_m0 + kcomp + i); 

            auto tg_0_xz_0 = braBuffer.data(pidx_g_0_2_m0 + 2 * kcomp + i); 

            auto tg_0_yy_0 = braBuffer.data(pidx_g_0_2_m0 + 3 * kcomp + i); 

            auto tg_0_yz_0 = braBuffer.data(pidx_g_0_2_m0 + 4 * kcomp + i); 

            auto tg_0_zz_0 = braBuffer.data(pidx_g_0_2_m0 + 5 * kcomp + i); 

            // set up pointers to integrals

            auto tg_x_x_0 = braBuffer.data(pidx_g_1_1_m0 + i); 

            auto tg_x_y_0 = braBuffer.data(pidx_g_1_1_m0 + kcomp + i); 

            auto tg_x_z_0 = braBuffer.data(pidx_g_1_1_m0 + 2 * kcomp + i); 

            auto tg_y_x_0 = braBuffer.data(pidx_g_1_1_m0 + 3 * kcomp + i); 

            auto tg_y_y_0 = braBuffer.data(pidx_g_1_1_m0 + 4 * kcomp + i); 

            auto tg_y_z_0 = braBuffer.data(pidx_g_1_1_m0 + 5 * kcomp + i); 

            auto tg_z_x_0 = braBuffer.data(pidx_g_1_1_m0 + 6 * kcomp + i); 

            auto tg_z_y_0 = braBuffer.data(pidx_g_1_1_m0 + 7 * kcomp + i); 

            auto tg_z_z_0 = braBuffer.data(pidx_g_1_1_m0 + 8 * kcomp + i); 

            #pragma omp simd aligned(tg_0_x_0, tg_0_xx_0, tg_0_xy_0, tg_0_xz_0, tg_0_y_0, tg_0_yy_0, \
                                     tg_0_yz_0, tg_0_z_0, tg_0_zz_0, tg_x_x_0, tg_x_y_0, tg_x_z_0, tg_y_x_0, tg_y_y_0, \
                                     tg_y_z_0, tg_z_x_0, tg_z_y_0, tg_z_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_x_x_0[j] = -ab_x * tg_0_x_0[j] + tg_0_xx_0[j];

                tg_x_y_0[j] = -ab_x * tg_0_y_0[j] + tg_0_xy_0[j];

                tg_x_z_0[j] = -ab_x * tg_0_z_0[j] + tg_0_xz_0[j];

                tg_y_x_0[j] = -ab_y * tg_0_x_0[j] + tg_0_xy_0[j];

                tg_y_y_0[j] = -ab_y * tg_0_y_0[j] + tg_0_yy_0[j];

                tg_y_z_0[j] = -ab_y * tg_0_z_0[j] + tg_0_yz_0[j];

                tg_z_x_0[j] = -ab_z * tg_0_x_0[j] + tg_0_xz_0[j];

                tg_z_y_0[j] = -ab_z * tg_0_y_0[j] + tg_0_yz_0[j];

                tg_z_z_0[j] = -ab_z * tg_0_z_0[j] + tg_0_zz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForPDXY(      CMemBlock2D<double>& braBuffer,
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

        auto pidx_g_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 2, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_1_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 2, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 3, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_0_xx_0 = braBuffer.data(pidx_g_0_2_m0 + i); 

            auto tg_0_xy_0 = braBuffer.data(pidx_g_0_2_m0 + kcomp + i); 

            auto tg_0_xz_0 = braBuffer.data(pidx_g_0_2_m0 + 2 * kcomp + i); 

            auto tg_0_yy_0 = braBuffer.data(pidx_g_0_2_m0 + 3 * kcomp + i); 

            auto tg_0_yz_0 = braBuffer.data(pidx_g_0_2_m0 + 4 * kcomp + i); 

            auto tg_0_zz_0 = braBuffer.data(pidx_g_0_2_m0 + 5 * kcomp + i); 

            auto tg_0_xxx_0 = braBuffer.data(pidx_g_0_3_m0 + i); 

            auto tg_0_xxy_0 = braBuffer.data(pidx_g_0_3_m0 + kcomp + i); 

            auto tg_0_xxz_0 = braBuffer.data(pidx_g_0_3_m0 + 2 * kcomp + i); 

            auto tg_0_xyy_0 = braBuffer.data(pidx_g_0_3_m0 + 3 * kcomp + i); 

            auto tg_0_xyz_0 = braBuffer.data(pidx_g_0_3_m0 + 4 * kcomp + i); 

            auto tg_0_xzz_0 = braBuffer.data(pidx_g_0_3_m0 + 5 * kcomp + i); 

            auto tg_0_yyy_0 = braBuffer.data(pidx_g_0_3_m0 + 6 * kcomp + i); 

            auto tg_0_yyz_0 = braBuffer.data(pidx_g_0_3_m0 + 7 * kcomp + i); 

            auto tg_0_yzz_0 = braBuffer.data(pidx_g_0_3_m0 + 8 * kcomp + i); 

            auto tg_0_zzz_0 = braBuffer.data(pidx_g_0_3_m0 + 9 * kcomp + i); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(tg_0_xx_0, tg_0_xxx_0, tg_0_xxy_0, tg_0_xxz_0, tg_0_xy_0, tg_0_xyy_0, \
                                     tg_0_xyz_0, tg_0_xz_0, tg_0_xzz_0, tg_0_yy_0, tg_0_yyy_0, tg_0_yyz_0, tg_0_yz_0, \
                                     tg_0_yzz_0, tg_0_zz_0, tg_0_zzz_0, tg_x_xx_0, tg_x_xy_0, tg_x_xz_0, tg_x_yy_0, \
                                     tg_x_yz_0, tg_x_zz_0, tg_y_xx_0, tg_y_xy_0, tg_y_xz_0, tg_y_yy_0, tg_y_yz_0, \
                                     tg_y_zz_0, tg_z_xx_0, tg_z_xy_0, tg_z_xz_0, tg_z_yy_0, tg_z_yz_0, tg_z_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_x_xx_0[j] = -ab_x * tg_0_xx_0[j] + tg_0_xxx_0[j];

                tg_x_xy_0[j] = -ab_x * tg_0_xy_0[j] + tg_0_xxy_0[j];

                tg_x_xz_0[j] = -ab_x * tg_0_xz_0[j] + tg_0_xxz_0[j];

                tg_x_yy_0[j] = -ab_x * tg_0_yy_0[j] + tg_0_xyy_0[j];

                tg_x_yz_0[j] = -ab_x * tg_0_yz_0[j] + tg_0_xyz_0[j];

                tg_x_zz_0[j] = -ab_x * tg_0_zz_0[j] + tg_0_xzz_0[j];

                tg_y_xx_0[j] = -ab_y * tg_0_xx_0[j] + tg_0_xxy_0[j];

                tg_y_xy_0[j] = -ab_y * tg_0_xy_0[j] + tg_0_xyy_0[j];

                tg_y_xz_0[j] = -ab_y * tg_0_xz_0[j] + tg_0_xyz_0[j];

                tg_y_yy_0[j] = -ab_y * tg_0_yy_0[j] + tg_0_yyy_0[j];

                tg_y_yz_0[j] = -ab_y * tg_0_yz_0[j] + tg_0_yyz_0[j];

                tg_y_zz_0[j] = -ab_y * tg_0_zz_0[j] + tg_0_yzz_0[j];

                tg_z_xx_0[j] = -ab_z * tg_0_xx_0[j] + tg_0_xxz_0[j];

                tg_z_xy_0[j] = -ab_z * tg_0_xy_0[j] + tg_0_xyz_0[j];

                tg_z_xz_0[j] = -ab_z * tg_0_xz_0[j] + tg_0_xzz_0[j];

                tg_z_yy_0[j] = -ab_z * tg_0_yy_0[j] + tg_0_yyz_0[j];

                tg_z_yz_0[j] = -ab_z * tg_0_yz_0[j] + tg_0_yzz_0[j];

                tg_z_zz_0[j] = -ab_z * tg_0_zz_0[j] + tg_0_zzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForPFXY(      CMemBlock2D<double>& braBuffer,
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

        auto pidx_g_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 3, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_1_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 3, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_0_xxx_0 = braBuffer.data(pidx_g_0_3_m0 + i); 

            auto tg_0_xxy_0 = braBuffer.data(pidx_g_0_3_m0 + kcomp + i); 

            auto tg_0_xxz_0 = braBuffer.data(pidx_g_0_3_m0 + 2 * kcomp + i); 

            auto tg_0_xyy_0 = braBuffer.data(pidx_g_0_3_m0 + 3 * kcomp + i); 

            auto tg_0_xyz_0 = braBuffer.data(pidx_g_0_3_m0 + 4 * kcomp + i); 

            auto tg_0_xzz_0 = braBuffer.data(pidx_g_0_3_m0 + 5 * kcomp + i); 

            auto tg_0_yyy_0 = braBuffer.data(pidx_g_0_3_m0 + 6 * kcomp + i); 

            auto tg_0_yyz_0 = braBuffer.data(pidx_g_0_3_m0 + 7 * kcomp + i); 

            auto tg_0_yzz_0 = braBuffer.data(pidx_g_0_3_m0 + 8 * kcomp + i); 

            auto tg_0_zzz_0 = braBuffer.data(pidx_g_0_3_m0 + 9 * kcomp + i); 

            auto tg_0_xxxx_0 = braBuffer.data(pidx_g_0_4_m0 + i); 

            auto tg_0_xxxy_0 = braBuffer.data(pidx_g_0_4_m0 + kcomp + i); 

            auto tg_0_xxxz_0 = braBuffer.data(pidx_g_0_4_m0 + 2 * kcomp + i); 

            auto tg_0_xxyy_0 = braBuffer.data(pidx_g_0_4_m0 + 3 * kcomp + i); 

            auto tg_0_xxyz_0 = braBuffer.data(pidx_g_0_4_m0 + 4 * kcomp + i); 

            auto tg_0_xxzz_0 = braBuffer.data(pidx_g_0_4_m0 + 5 * kcomp + i); 

            auto tg_0_xyyy_0 = braBuffer.data(pidx_g_0_4_m0 + 6 * kcomp + i); 

            auto tg_0_xyyz_0 = braBuffer.data(pidx_g_0_4_m0 + 7 * kcomp + i); 

            auto tg_0_xyzz_0 = braBuffer.data(pidx_g_0_4_m0 + 8 * kcomp + i); 

            auto tg_0_xzzz_0 = braBuffer.data(pidx_g_0_4_m0 + 9 * kcomp + i); 

            auto tg_0_yyyy_0 = braBuffer.data(pidx_g_0_4_m0 + 10 * kcomp + i); 

            auto tg_0_yyyz_0 = braBuffer.data(pidx_g_0_4_m0 + 11 * kcomp + i); 

            auto tg_0_yyzz_0 = braBuffer.data(pidx_g_0_4_m0 + 12 * kcomp + i); 

            auto tg_0_yzzz_0 = braBuffer.data(pidx_g_0_4_m0 + 13 * kcomp + i); 

            auto tg_0_zzzz_0 = braBuffer.data(pidx_g_0_4_m0 + 14 * kcomp + i); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(tg_0_xxx_0, tg_0_xxxx_0, tg_0_xxxy_0, tg_0_xxxz_0, tg_0_xxy_0, \
                                     tg_0_xxyy_0, tg_0_xxyz_0, tg_0_xxz_0, tg_0_xxzz_0, tg_0_xyy_0, tg_0_xyyy_0, \
                                     tg_0_xyyz_0, tg_0_xyz_0, tg_0_xyzz_0, tg_0_xzz_0, tg_0_xzzz_0, tg_0_yyy_0, \
                                     tg_0_yyyy_0, tg_0_yyyz_0, tg_0_yyz_0, tg_0_yyzz_0, tg_0_yzz_0, tg_0_yzzz_0, \
                                     tg_0_zzz_0, tg_0_zzzz_0, tg_x_xxx_0, tg_x_xxy_0, tg_x_xxz_0, tg_x_xyy_0, \
                                     tg_x_xyz_0, tg_x_xzz_0, tg_x_yyy_0, tg_x_yyz_0, tg_x_yzz_0, tg_x_zzz_0, tg_y_xxx_0, \
                                     tg_y_xxy_0, tg_y_xxz_0, tg_y_xyy_0, tg_y_xyz_0, tg_y_xzz_0, tg_y_yyy_0, tg_y_yyz_0, \
                                     tg_y_yzz_0, tg_y_zzz_0, tg_z_xxx_0, tg_z_xxy_0, tg_z_xxz_0, tg_z_xyy_0, tg_z_xyz_0, \
                                     tg_z_xzz_0, tg_z_yyy_0, tg_z_yyz_0, tg_z_yzz_0, tg_z_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_x_xxx_0[j] = -ab_x * tg_0_xxx_0[j] + tg_0_xxxx_0[j];

                tg_x_xxy_0[j] = -ab_x * tg_0_xxy_0[j] + tg_0_xxxy_0[j];

                tg_x_xxz_0[j] = -ab_x * tg_0_xxz_0[j] + tg_0_xxxz_0[j];

                tg_x_xyy_0[j] = -ab_x * tg_0_xyy_0[j] + tg_0_xxyy_0[j];

                tg_x_xyz_0[j] = -ab_x * tg_0_xyz_0[j] + tg_0_xxyz_0[j];

                tg_x_xzz_0[j] = -ab_x * tg_0_xzz_0[j] + tg_0_xxzz_0[j];

                tg_x_yyy_0[j] = -ab_x * tg_0_yyy_0[j] + tg_0_xyyy_0[j];

                tg_x_yyz_0[j] = -ab_x * tg_0_yyz_0[j] + tg_0_xyyz_0[j];

                tg_x_yzz_0[j] = -ab_x * tg_0_yzz_0[j] + tg_0_xyzz_0[j];

                tg_x_zzz_0[j] = -ab_x * tg_0_zzz_0[j] + tg_0_xzzz_0[j];

                tg_y_xxx_0[j] = -ab_y * tg_0_xxx_0[j] + tg_0_xxxy_0[j];

                tg_y_xxy_0[j] = -ab_y * tg_0_xxy_0[j] + tg_0_xxyy_0[j];

                tg_y_xxz_0[j] = -ab_y * tg_0_xxz_0[j] + tg_0_xxyz_0[j];

                tg_y_xyy_0[j] = -ab_y * tg_0_xyy_0[j] + tg_0_xyyy_0[j];

                tg_y_xyz_0[j] = -ab_y * tg_0_xyz_0[j] + tg_0_xyyz_0[j];

                tg_y_xzz_0[j] = -ab_y * tg_0_xzz_0[j] + tg_0_xyzz_0[j];

                tg_y_yyy_0[j] = -ab_y * tg_0_yyy_0[j] + tg_0_yyyy_0[j];

                tg_y_yyz_0[j] = -ab_y * tg_0_yyz_0[j] + tg_0_yyyz_0[j];

                tg_y_yzz_0[j] = -ab_y * tg_0_yzz_0[j] + tg_0_yyzz_0[j];

                tg_y_zzz_0[j] = -ab_y * tg_0_zzz_0[j] + tg_0_yzzz_0[j];

                tg_z_xxx_0[j] = -ab_z * tg_0_xxx_0[j] + tg_0_xxxz_0[j];

                tg_z_xxy_0[j] = -ab_z * tg_0_xxy_0[j] + tg_0_xxyz_0[j];

                tg_z_xxz_0[j] = -ab_z * tg_0_xxz_0[j] + tg_0_xxzz_0[j];

                tg_z_xyy_0[j] = -ab_z * tg_0_xyy_0[j] + tg_0_xyyz_0[j];

                tg_z_xyz_0[j] = -ab_z * tg_0_xyz_0[j] + tg_0_xyzz_0[j];

                tg_z_xzz_0[j] = -ab_z * tg_0_xzz_0[j] + tg_0_xzzz_0[j];

                tg_z_yyy_0[j] = -ab_z * tg_0_yyy_0[j] + tg_0_yyyz_0[j];

                tg_z_yyz_0[j] = -ab_z * tg_0_yyz_0[j] + tg_0_yyzz_0[j];

                tg_z_yzz_0[j] = -ab_z * tg_0_yzz_0[j] + tg_0_yzzz_0[j];

                tg_z_zzz_0[j] = -ab_z * tg_0_zzz_0[j] + tg_0_zzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForPGXY(      CMemBlock2D<double>& braBuffer,
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

        auto pidx_g_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_1_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_0_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_0_xxxx_0 = braBuffer.data(pidx_g_0_4_m0 + i); 

            auto tg_0_xxxy_0 = braBuffer.data(pidx_g_0_4_m0 + kcomp + i); 

            auto tg_0_xxxz_0 = braBuffer.data(pidx_g_0_4_m0 + 2 * kcomp + i); 

            auto tg_0_xxyy_0 = braBuffer.data(pidx_g_0_4_m0 + 3 * kcomp + i); 

            auto tg_0_xxyz_0 = braBuffer.data(pidx_g_0_4_m0 + 4 * kcomp + i); 

            auto tg_0_xxzz_0 = braBuffer.data(pidx_g_0_4_m0 + 5 * kcomp + i); 

            auto tg_0_xyyy_0 = braBuffer.data(pidx_g_0_4_m0 + 6 * kcomp + i); 

            auto tg_0_xyyz_0 = braBuffer.data(pidx_g_0_4_m0 + 7 * kcomp + i); 

            auto tg_0_xyzz_0 = braBuffer.data(pidx_g_0_4_m0 + 8 * kcomp + i); 

            auto tg_0_xzzz_0 = braBuffer.data(pidx_g_0_4_m0 + 9 * kcomp + i); 

            auto tg_0_yyyy_0 = braBuffer.data(pidx_g_0_4_m0 + 10 * kcomp + i); 

            auto tg_0_yyyz_0 = braBuffer.data(pidx_g_0_4_m0 + 11 * kcomp + i); 

            auto tg_0_yyzz_0 = braBuffer.data(pidx_g_0_4_m0 + 12 * kcomp + i); 

            auto tg_0_yzzz_0 = braBuffer.data(pidx_g_0_4_m0 + 13 * kcomp + i); 

            auto tg_0_zzzz_0 = braBuffer.data(pidx_g_0_4_m0 + 14 * kcomp + i); 

            auto tg_0_xxxxx_0 = braBuffer.data(pidx_g_0_5_m0 + i); 

            auto tg_0_xxxxy_0 = braBuffer.data(pidx_g_0_5_m0 + kcomp + i); 

            auto tg_0_xxxxz_0 = braBuffer.data(pidx_g_0_5_m0 + 2 * kcomp + i); 

            auto tg_0_xxxyy_0 = braBuffer.data(pidx_g_0_5_m0 + 3 * kcomp + i); 

            auto tg_0_xxxyz_0 = braBuffer.data(pidx_g_0_5_m0 + 4 * kcomp + i); 

            auto tg_0_xxxzz_0 = braBuffer.data(pidx_g_0_5_m0 + 5 * kcomp + i); 

            auto tg_0_xxyyy_0 = braBuffer.data(pidx_g_0_5_m0 + 6 * kcomp + i); 

            auto tg_0_xxyyz_0 = braBuffer.data(pidx_g_0_5_m0 + 7 * kcomp + i); 

            auto tg_0_xxyzz_0 = braBuffer.data(pidx_g_0_5_m0 + 8 * kcomp + i); 

            auto tg_0_xxzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 9 * kcomp + i); 

            auto tg_0_xyyyy_0 = braBuffer.data(pidx_g_0_5_m0 + 10 * kcomp + i); 

            auto tg_0_xyyyz_0 = braBuffer.data(pidx_g_0_5_m0 + 11 * kcomp + i); 

            auto tg_0_xyyzz_0 = braBuffer.data(pidx_g_0_5_m0 + 12 * kcomp + i); 

            auto tg_0_xyzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 13 * kcomp + i); 

            auto tg_0_xzzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 14 * kcomp + i); 

            auto tg_0_yyyyy_0 = braBuffer.data(pidx_g_0_5_m0 + 15 * kcomp + i); 

            auto tg_0_yyyyz_0 = braBuffer.data(pidx_g_0_5_m0 + 16 * kcomp + i); 

            auto tg_0_yyyzz_0 = braBuffer.data(pidx_g_0_5_m0 + 17 * kcomp + i); 

            auto tg_0_yyzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 18 * kcomp + i); 

            auto tg_0_yzzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 19 * kcomp + i); 

            auto tg_0_zzzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 20 * kcomp + i); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(tg_0_xxxx_0, tg_0_xxxxx_0, tg_0_xxxxy_0, tg_0_xxxxz_0, tg_0_xxxy_0, \
                                     tg_0_xxxyy_0, tg_0_xxxyz_0, tg_0_xxxz_0, tg_0_xxxzz_0, tg_0_xxyy_0, tg_0_xxyyy_0, \
                                     tg_0_xxyyz_0, tg_0_xxyz_0, tg_0_xxyzz_0, tg_0_xxzz_0, tg_0_xxzzz_0, tg_0_xyyy_0, \
                                     tg_0_xyyyy_0, tg_0_xyyyz_0, tg_0_xyyz_0, tg_0_xyyzz_0, tg_0_xyzz_0, tg_0_xyzzz_0, \
                                     tg_0_xzzz_0, tg_0_xzzzz_0, tg_0_yyyy_0, tg_0_yyyyy_0, tg_0_yyyyz_0, tg_0_yyyz_0, \
                                     tg_0_yyyzz_0, tg_0_yyzz_0, tg_0_yyzzz_0, tg_0_yzzz_0, tg_0_yzzzz_0, tg_0_zzzz_0, \
                                     tg_0_zzzzz_0, tg_x_xxxx_0, tg_x_xxxy_0, tg_x_xxxz_0, tg_x_xxyy_0, tg_x_xxyz_0, \
                                     tg_x_xxzz_0, tg_x_xyyy_0, tg_x_xyyz_0, tg_x_xyzz_0, tg_x_xzzz_0, tg_x_yyyy_0, \
                                     tg_x_yyyz_0, tg_x_yyzz_0, tg_x_yzzz_0, tg_x_zzzz_0, tg_y_xxxx_0, tg_y_xxxy_0, \
                                     tg_y_xxxz_0, tg_y_xxyy_0, tg_y_xxyz_0, tg_y_xxzz_0, tg_y_xyyy_0, tg_y_xyyz_0, \
                                     tg_y_xyzz_0, tg_y_xzzz_0, tg_y_yyyy_0, tg_y_yyyz_0, tg_y_yyzz_0, tg_y_yzzz_0, \
                                     tg_y_zzzz_0, tg_z_xxxx_0, tg_z_xxxy_0, tg_z_xxxz_0, tg_z_xxyy_0, tg_z_xxyz_0, \
                                     tg_z_xxzz_0, tg_z_xyyy_0, tg_z_xyyz_0, tg_z_xyzz_0, tg_z_xzzz_0, tg_z_yyyy_0, \
                                     tg_z_yyyz_0, tg_z_yyzz_0, tg_z_yzzz_0, tg_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_x_xxxx_0[j] = -ab_x * tg_0_xxxx_0[j] + tg_0_xxxxx_0[j];

                tg_x_xxxy_0[j] = -ab_x * tg_0_xxxy_0[j] + tg_0_xxxxy_0[j];

                tg_x_xxxz_0[j] = -ab_x * tg_0_xxxz_0[j] + tg_0_xxxxz_0[j];

                tg_x_xxyy_0[j] = -ab_x * tg_0_xxyy_0[j] + tg_0_xxxyy_0[j];

                tg_x_xxyz_0[j] = -ab_x * tg_0_xxyz_0[j] + tg_0_xxxyz_0[j];

                tg_x_xxzz_0[j] = -ab_x * tg_0_xxzz_0[j] + tg_0_xxxzz_0[j];

                tg_x_xyyy_0[j] = -ab_x * tg_0_xyyy_0[j] + tg_0_xxyyy_0[j];

                tg_x_xyyz_0[j] = -ab_x * tg_0_xyyz_0[j] + tg_0_xxyyz_0[j];

                tg_x_xyzz_0[j] = -ab_x * tg_0_xyzz_0[j] + tg_0_xxyzz_0[j];

                tg_x_xzzz_0[j] = -ab_x * tg_0_xzzz_0[j] + tg_0_xxzzz_0[j];

                tg_x_yyyy_0[j] = -ab_x * tg_0_yyyy_0[j] + tg_0_xyyyy_0[j];

                tg_x_yyyz_0[j] = -ab_x * tg_0_yyyz_0[j] + tg_0_xyyyz_0[j];

                tg_x_yyzz_0[j] = -ab_x * tg_0_yyzz_0[j] + tg_0_xyyzz_0[j];

                tg_x_yzzz_0[j] = -ab_x * tg_0_yzzz_0[j] + tg_0_xyzzz_0[j];

                tg_x_zzzz_0[j] = -ab_x * tg_0_zzzz_0[j] + tg_0_xzzzz_0[j];

                tg_y_xxxx_0[j] = -ab_y * tg_0_xxxx_0[j] + tg_0_xxxxy_0[j];

                tg_y_xxxy_0[j] = -ab_y * tg_0_xxxy_0[j] + tg_0_xxxyy_0[j];

                tg_y_xxxz_0[j] = -ab_y * tg_0_xxxz_0[j] + tg_0_xxxyz_0[j];

                tg_y_xxyy_0[j] = -ab_y * tg_0_xxyy_0[j] + tg_0_xxyyy_0[j];

                tg_y_xxyz_0[j] = -ab_y * tg_0_xxyz_0[j] + tg_0_xxyyz_0[j];

                tg_y_xxzz_0[j] = -ab_y * tg_0_xxzz_0[j] + tg_0_xxyzz_0[j];

                tg_y_xyyy_0[j] = -ab_y * tg_0_xyyy_0[j] + tg_0_xyyyy_0[j];

                tg_y_xyyz_0[j] = -ab_y * tg_0_xyyz_0[j] + tg_0_xyyyz_0[j];

                tg_y_xyzz_0[j] = -ab_y * tg_0_xyzz_0[j] + tg_0_xyyzz_0[j];

                tg_y_xzzz_0[j] = -ab_y * tg_0_xzzz_0[j] + tg_0_xyzzz_0[j];

                tg_y_yyyy_0[j] = -ab_y * tg_0_yyyy_0[j] + tg_0_yyyyy_0[j];

                tg_y_yyyz_0[j] = -ab_y * tg_0_yyyz_0[j] + tg_0_yyyyz_0[j];

                tg_y_yyzz_0[j] = -ab_y * tg_0_yyzz_0[j] + tg_0_yyyzz_0[j];

                tg_y_yzzz_0[j] = -ab_y * tg_0_yzzz_0[j] + tg_0_yyzzz_0[j];

                tg_y_zzzz_0[j] = -ab_y * tg_0_zzzz_0[j] + tg_0_yzzzz_0[j];

                tg_z_xxxx_0[j] = -ab_z * tg_0_xxxx_0[j] + tg_0_xxxxz_0[j];

                tg_z_xxxy_0[j] = -ab_z * tg_0_xxxy_0[j] + tg_0_xxxyz_0[j];

                tg_z_xxxz_0[j] = -ab_z * tg_0_xxxz_0[j] + tg_0_xxxzz_0[j];

                tg_z_xxyy_0[j] = -ab_z * tg_0_xxyy_0[j] + tg_0_xxyyz_0[j];

                tg_z_xxyz_0[j] = -ab_z * tg_0_xxyz_0[j] + tg_0_xxyzz_0[j];

                tg_z_xxzz_0[j] = -ab_z * tg_0_xxzz_0[j] + tg_0_xxzzz_0[j];

                tg_z_xyyy_0[j] = -ab_z * tg_0_xyyy_0[j] + tg_0_xyyyz_0[j];

                tg_z_xyyz_0[j] = -ab_z * tg_0_xyyz_0[j] + tg_0_xyyzz_0[j];

                tg_z_xyzz_0[j] = -ab_z * tg_0_xyzz_0[j] + tg_0_xyzzz_0[j];

                tg_z_xzzz_0[j] = -ab_z * tg_0_xzzz_0[j] + tg_0_xzzzz_0[j];

                tg_z_yyyy_0[j] = -ab_z * tg_0_yyyy_0[j] + tg_0_yyyyz_0[j];

                tg_z_yyyz_0[j] = -ab_z * tg_0_yyyz_0[j] + tg_0_yyyzz_0[j];

                tg_z_yyzz_0[j] = -ab_z * tg_0_yyzz_0[j] + tg_0_yyzzz_0[j];

                tg_z_yzzz_0[j] = -ab_z * tg_0_yzzz_0[j] + tg_0_yzzzz_0[j];

                tg_z_zzzz_0[j] = -ab_z * tg_0_zzzz_0[j] + tg_0_zzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForPHXY(      CMemBlock2D<double>& braBuffer,
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

        auto pidx_g_1_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_1_5_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_0_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_0_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_0_xxxxx_0 = braBuffer.data(pidx_g_0_5_m0 + i); 

            auto tg_0_xxxxy_0 = braBuffer.data(pidx_g_0_5_m0 + kcomp + i); 

            auto tg_0_xxxxz_0 = braBuffer.data(pidx_g_0_5_m0 + 2 * kcomp + i); 

            auto tg_0_xxxyy_0 = braBuffer.data(pidx_g_0_5_m0 + 3 * kcomp + i); 

            auto tg_0_xxxyz_0 = braBuffer.data(pidx_g_0_5_m0 + 4 * kcomp + i); 

            auto tg_0_xxxzz_0 = braBuffer.data(pidx_g_0_5_m0 + 5 * kcomp + i); 

            auto tg_0_xxyyy_0 = braBuffer.data(pidx_g_0_5_m0 + 6 * kcomp + i); 

            auto tg_0_xxyyz_0 = braBuffer.data(pidx_g_0_5_m0 + 7 * kcomp + i); 

            auto tg_0_xxyzz_0 = braBuffer.data(pidx_g_0_5_m0 + 8 * kcomp + i); 

            auto tg_0_xxzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 9 * kcomp + i); 

            auto tg_0_xyyyy_0 = braBuffer.data(pidx_g_0_5_m0 + 10 * kcomp + i); 

            auto tg_0_xyyyz_0 = braBuffer.data(pidx_g_0_5_m0 + 11 * kcomp + i); 

            auto tg_0_xyyzz_0 = braBuffer.data(pidx_g_0_5_m0 + 12 * kcomp + i); 

            auto tg_0_xyzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 13 * kcomp + i); 

            auto tg_0_xzzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 14 * kcomp + i); 

            auto tg_0_yyyyy_0 = braBuffer.data(pidx_g_0_5_m0 + 15 * kcomp + i); 

            auto tg_0_yyyyz_0 = braBuffer.data(pidx_g_0_5_m0 + 16 * kcomp + i); 

            auto tg_0_yyyzz_0 = braBuffer.data(pidx_g_0_5_m0 + 17 * kcomp + i); 

            auto tg_0_yyzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 18 * kcomp + i); 

            auto tg_0_yzzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 19 * kcomp + i); 

            auto tg_0_zzzzz_0 = braBuffer.data(pidx_g_0_5_m0 + 20 * kcomp + i); 

            auto tg_0_xxxxxx_0 = braBuffer.data(pidx_g_0_6_m0 + i); 

            auto tg_0_xxxxxy_0 = braBuffer.data(pidx_g_0_6_m0 + kcomp + i); 

            auto tg_0_xxxxxz_0 = braBuffer.data(pidx_g_0_6_m0 + 2 * kcomp + i); 

            auto tg_0_xxxxyy_0 = braBuffer.data(pidx_g_0_6_m0 + 3 * kcomp + i); 

            auto tg_0_xxxxyz_0 = braBuffer.data(pidx_g_0_6_m0 + 4 * kcomp + i); 

            auto tg_0_xxxxzz_0 = braBuffer.data(pidx_g_0_6_m0 + 5 * kcomp + i); 

            auto tg_0_xxxyyy_0 = braBuffer.data(pidx_g_0_6_m0 + 6 * kcomp + i); 

            auto tg_0_xxxyyz_0 = braBuffer.data(pidx_g_0_6_m0 + 7 * kcomp + i); 

            auto tg_0_xxxyzz_0 = braBuffer.data(pidx_g_0_6_m0 + 8 * kcomp + i); 

            auto tg_0_xxxzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 9 * kcomp + i); 

            auto tg_0_xxyyyy_0 = braBuffer.data(pidx_g_0_6_m0 + 10 * kcomp + i); 

            auto tg_0_xxyyyz_0 = braBuffer.data(pidx_g_0_6_m0 + 11 * kcomp + i); 

            auto tg_0_xxyyzz_0 = braBuffer.data(pidx_g_0_6_m0 + 12 * kcomp + i); 

            auto tg_0_xxyzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 13 * kcomp + i); 

            auto tg_0_xxzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 14 * kcomp + i); 

            auto tg_0_xyyyyy_0 = braBuffer.data(pidx_g_0_6_m0 + 15 * kcomp + i); 

            auto tg_0_xyyyyz_0 = braBuffer.data(pidx_g_0_6_m0 + 16 * kcomp + i); 

            auto tg_0_xyyyzz_0 = braBuffer.data(pidx_g_0_6_m0 + 17 * kcomp + i); 

            auto tg_0_xyyzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 18 * kcomp + i); 

            auto tg_0_xyzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 19 * kcomp + i); 

            auto tg_0_xzzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 20 * kcomp + i); 

            auto tg_0_yyyyyy_0 = braBuffer.data(pidx_g_0_6_m0 + 21 * kcomp + i); 

            auto tg_0_yyyyyz_0 = braBuffer.data(pidx_g_0_6_m0 + 22 * kcomp + i); 

            auto tg_0_yyyyzz_0 = braBuffer.data(pidx_g_0_6_m0 + 23 * kcomp + i); 

            auto tg_0_yyyzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 24 * kcomp + i); 

            auto tg_0_yyzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 25 * kcomp + i); 

            auto tg_0_yzzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 26 * kcomp + i); 

            auto tg_0_zzzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 27 * kcomp + i); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(tg_0_xxxxx_0, tg_0_xxxxxx_0, tg_0_xxxxxy_0, tg_0_xxxxxz_0, \
                                     tg_0_xxxxy_0, tg_0_xxxxyy_0, tg_0_xxxxyz_0, tg_0_xxxxz_0, tg_0_xxxxzz_0, \
                                     tg_0_xxxyy_0, tg_0_xxxyyy_0, tg_0_xxxyyz_0, tg_0_xxxyz_0, tg_0_xxxyzz_0, \
                                     tg_0_xxxzz_0, tg_0_xxxzzz_0, tg_0_xxyyy_0, tg_0_xxyyyy_0, tg_0_xxyyyz_0, \
                                     tg_0_xxyyz_0, tg_0_xxyyzz_0, tg_0_xxyzz_0, tg_0_xxyzzz_0, tg_0_xxzzz_0, \
                                     tg_0_xxzzzz_0, tg_0_xyyyy_0, tg_0_xyyyyy_0, tg_0_xyyyyz_0, tg_0_xyyyz_0, \
                                     tg_0_xyyyzz_0, tg_0_xyyzz_0, tg_0_xyyzzz_0, tg_0_xyzzz_0, tg_0_xyzzzz_0, \
                                     tg_0_xzzzz_0, tg_0_xzzzzz_0, tg_0_yyyyy_0, tg_0_yyyyyy_0, tg_0_yyyyyz_0, \
                                     tg_0_yyyyz_0, tg_0_yyyyzz_0, tg_0_yyyzz_0, tg_0_yyyzzz_0, tg_0_yyzzz_0, \
                                     tg_0_yyzzzz_0, tg_0_yzzzz_0, tg_0_yzzzzz_0, tg_0_zzzzz_0, tg_0_zzzzzz_0, \
                                     tg_x_xxxxx_0, tg_x_xxxxy_0, tg_x_xxxxz_0, tg_x_xxxyy_0, tg_x_xxxyz_0, tg_x_xxxzz_0, \
                                     tg_x_xxyyy_0, tg_x_xxyyz_0, tg_x_xxyzz_0, tg_x_xxzzz_0, tg_x_xyyyy_0, tg_x_xyyyz_0, \
                                     tg_x_xyyzz_0, tg_x_xyzzz_0, tg_x_xzzzz_0, tg_x_yyyyy_0, tg_x_yyyyz_0, tg_x_yyyzz_0, \
                                     tg_x_yyzzz_0, tg_x_yzzzz_0, tg_x_zzzzz_0, tg_y_xxxxx_0, tg_y_xxxxy_0, tg_y_xxxxz_0, \
                                     tg_y_xxxyy_0, tg_y_xxxyz_0, tg_y_xxxzz_0, tg_y_xxyyy_0, tg_y_xxyyz_0, tg_y_xxyzz_0, \
                                     tg_y_xxzzz_0, tg_y_xyyyy_0, tg_y_xyyyz_0, tg_y_xyyzz_0, tg_y_xyzzz_0, tg_y_xzzzz_0, \
                                     tg_y_yyyyy_0, tg_y_yyyyz_0, tg_y_yyyzz_0, tg_y_yyzzz_0, tg_y_yzzzz_0, tg_y_zzzzz_0, \
                                     tg_z_xxxxx_0, tg_z_xxxxy_0, tg_z_xxxxz_0, tg_z_xxxyy_0, tg_z_xxxyz_0, tg_z_xxxzz_0, \
                                     tg_z_xxyyy_0, tg_z_xxyyz_0, tg_z_xxyzz_0, tg_z_xxzzz_0, tg_z_xyyyy_0, tg_z_xyyyz_0, \
                                     tg_z_xyyzz_0, tg_z_xyzzz_0, tg_z_xzzzz_0, tg_z_yyyyy_0, tg_z_yyyyz_0, tg_z_yyyzz_0, \
                                     tg_z_yyzzz_0, tg_z_yzzzz_0, tg_z_zzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_x_xxxxx_0[j] = -ab_x * tg_0_xxxxx_0[j] + tg_0_xxxxxx_0[j];

                tg_x_xxxxy_0[j] = -ab_x * tg_0_xxxxy_0[j] + tg_0_xxxxxy_0[j];

                tg_x_xxxxz_0[j] = -ab_x * tg_0_xxxxz_0[j] + tg_0_xxxxxz_0[j];

                tg_x_xxxyy_0[j] = -ab_x * tg_0_xxxyy_0[j] + tg_0_xxxxyy_0[j];

                tg_x_xxxyz_0[j] = -ab_x * tg_0_xxxyz_0[j] + tg_0_xxxxyz_0[j];

                tg_x_xxxzz_0[j] = -ab_x * tg_0_xxxzz_0[j] + tg_0_xxxxzz_0[j];

                tg_x_xxyyy_0[j] = -ab_x * tg_0_xxyyy_0[j] + tg_0_xxxyyy_0[j];

                tg_x_xxyyz_0[j] = -ab_x * tg_0_xxyyz_0[j] + tg_0_xxxyyz_0[j];

                tg_x_xxyzz_0[j] = -ab_x * tg_0_xxyzz_0[j] + tg_0_xxxyzz_0[j];

                tg_x_xxzzz_0[j] = -ab_x * tg_0_xxzzz_0[j] + tg_0_xxxzzz_0[j];

                tg_x_xyyyy_0[j] = -ab_x * tg_0_xyyyy_0[j] + tg_0_xxyyyy_0[j];

                tg_x_xyyyz_0[j] = -ab_x * tg_0_xyyyz_0[j] + tg_0_xxyyyz_0[j];

                tg_x_xyyzz_0[j] = -ab_x * tg_0_xyyzz_0[j] + tg_0_xxyyzz_0[j];

                tg_x_xyzzz_0[j] = -ab_x * tg_0_xyzzz_0[j] + tg_0_xxyzzz_0[j];

                tg_x_xzzzz_0[j] = -ab_x * tg_0_xzzzz_0[j] + tg_0_xxzzzz_0[j];

                tg_x_yyyyy_0[j] = -ab_x * tg_0_yyyyy_0[j] + tg_0_xyyyyy_0[j];

                tg_x_yyyyz_0[j] = -ab_x * tg_0_yyyyz_0[j] + tg_0_xyyyyz_0[j];

                tg_x_yyyzz_0[j] = -ab_x * tg_0_yyyzz_0[j] + tg_0_xyyyzz_0[j];

                tg_x_yyzzz_0[j] = -ab_x * tg_0_yyzzz_0[j] + tg_0_xyyzzz_0[j];

                tg_x_yzzzz_0[j] = -ab_x * tg_0_yzzzz_0[j] + tg_0_xyzzzz_0[j];

                tg_x_zzzzz_0[j] = -ab_x * tg_0_zzzzz_0[j] + tg_0_xzzzzz_0[j];

                tg_y_xxxxx_0[j] = -ab_y * tg_0_xxxxx_0[j] + tg_0_xxxxxy_0[j];

                tg_y_xxxxy_0[j] = -ab_y * tg_0_xxxxy_0[j] + tg_0_xxxxyy_0[j];

                tg_y_xxxxz_0[j] = -ab_y * tg_0_xxxxz_0[j] + tg_0_xxxxyz_0[j];

                tg_y_xxxyy_0[j] = -ab_y * tg_0_xxxyy_0[j] + tg_0_xxxyyy_0[j];

                tg_y_xxxyz_0[j] = -ab_y * tg_0_xxxyz_0[j] + tg_0_xxxyyz_0[j];

                tg_y_xxxzz_0[j] = -ab_y * tg_0_xxxzz_0[j] + tg_0_xxxyzz_0[j];

                tg_y_xxyyy_0[j] = -ab_y * tg_0_xxyyy_0[j] + tg_0_xxyyyy_0[j];

                tg_y_xxyyz_0[j] = -ab_y * tg_0_xxyyz_0[j] + tg_0_xxyyyz_0[j];

                tg_y_xxyzz_0[j] = -ab_y * tg_0_xxyzz_0[j] + tg_0_xxyyzz_0[j];

                tg_y_xxzzz_0[j] = -ab_y * tg_0_xxzzz_0[j] + tg_0_xxyzzz_0[j];

                tg_y_xyyyy_0[j] = -ab_y * tg_0_xyyyy_0[j] + tg_0_xyyyyy_0[j];

                tg_y_xyyyz_0[j] = -ab_y * tg_0_xyyyz_0[j] + tg_0_xyyyyz_0[j];

                tg_y_xyyzz_0[j] = -ab_y * tg_0_xyyzz_0[j] + tg_0_xyyyzz_0[j];

                tg_y_xyzzz_0[j] = -ab_y * tg_0_xyzzz_0[j] + tg_0_xyyzzz_0[j];

                tg_y_xzzzz_0[j] = -ab_y * tg_0_xzzzz_0[j] + tg_0_xyzzzz_0[j];

                tg_y_yyyyy_0[j] = -ab_y * tg_0_yyyyy_0[j] + tg_0_yyyyyy_0[j];

                tg_y_yyyyz_0[j] = -ab_y * tg_0_yyyyz_0[j] + tg_0_yyyyyz_0[j];

                tg_y_yyyzz_0[j] = -ab_y * tg_0_yyyzz_0[j] + tg_0_yyyyzz_0[j];

                tg_y_yyzzz_0[j] = -ab_y * tg_0_yyzzz_0[j] + tg_0_yyyzzz_0[j];

                tg_y_yzzzz_0[j] = -ab_y * tg_0_yzzzz_0[j] + tg_0_yyzzzz_0[j];

                tg_y_zzzzz_0[j] = -ab_y * tg_0_zzzzz_0[j] + tg_0_yzzzzz_0[j];

                tg_z_xxxxx_0[j] = -ab_z * tg_0_xxxxx_0[j] + tg_0_xxxxxz_0[j];

                tg_z_xxxxy_0[j] = -ab_z * tg_0_xxxxy_0[j] + tg_0_xxxxyz_0[j];

                tg_z_xxxxz_0[j] = -ab_z * tg_0_xxxxz_0[j] + tg_0_xxxxzz_0[j];

                tg_z_xxxyy_0[j] = -ab_z * tg_0_xxxyy_0[j] + tg_0_xxxyyz_0[j];

                tg_z_xxxyz_0[j] = -ab_z * tg_0_xxxyz_0[j] + tg_0_xxxyzz_0[j];

                tg_z_xxxzz_0[j] = -ab_z * tg_0_xxxzz_0[j] + tg_0_xxxzzz_0[j];

                tg_z_xxyyy_0[j] = -ab_z * tg_0_xxyyy_0[j] + tg_0_xxyyyz_0[j];

                tg_z_xxyyz_0[j] = -ab_z * tg_0_xxyyz_0[j] + tg_0_xxyyzz_0[j];

                tg_z_xxyzz_0[j] = -ab_z * tg_0_xxyzz_0[j] + tg_0_xxyzzz_0[j];

                tg_z_xxzzz_0[j] = -ab_z * tg_0_xxzzz_0[j] + tg_0_xxzzzz_0[j];

                tg_z_xyyyy_0[j] = -ab_z * tg_0_xyyyy_0[j] + tg_0_xyyyyz_0[j];

                tg_z_xyyyz_0[j] = -ab_z * tg_0_xyyyz_0[j] + tg_0_xyyyzz_0[j];

                tg_z_xyyzz_0[j] = -ab_z * tg_0_xyyzz_0[j] + tg_0_xyyzzz_0[j];

                tg_z_xyzzz_0[j] = -ab_z * tg_0_xyzzz_0[j] + tg_0_xyzzzz_0[j];

                tg_z_xzzzz_0[j] = -ab_z * tg_0_xzzzz_0[j] + tg_0_xzzzzz_0[j];

                tg_z_yyyyy_0[j] = -ab_z * tg_0_yyyyy_0[j] + tg_0_yyyyyz_0[j];

                tg_z_yyyyz_0[j] = -ab_z * tg_0_yyyyz_0[j] + tg_0_yyyyzz_0[j];

                tg_z_yyyzz_0[j] = -ab_z * tg_0_yyyzz_0[j] + tg_0_yyyzzz_0[j];

                tg_z_yyzzz_0[j] = -ab_z * tg_0_yyzzz_0[j] + tg_0_yyzzzz_0[j];

                tg_z_yzzzz_0[j] = -ab_z * tg_0_yzzzz_0[j] + tg_0_yzzzzz_0[j];

                tg_z_zzzzz_0[j] = -ab_z * tg_0_zzzzz_0[j] + tg_0_zzzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForPIXY(      CMemBlock2D<double>& braBuffer,
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

        auto pidx_g_1_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_1_6_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_0_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_0_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 7, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_0_xxxxxx_0 = braBuffer.data(pidx_g_0_6_m0 + i); 

            auto tg_0_xxxxxy_0 = braBuffer.data(pidx_g_0_6_m0 + kcomp + i); 

            auto tg_0_xxxxxz_0 = braBuffer.data(pidx_g_0_6_m0 + 2 * kcomp + i); 

            auto tg_0_xxxxyy_0 = braBuffer.data(pidx_g_0_6_m0 + 3 * kcomp + i); 

            auto tg_0_xxxxyz_0 = braBuffer.data(pidx_g_0_6_m0 + 4 * kcomp + i); 

            auto tg_0_xxxxzz_0 = braBuffer.data(pidx_g_0_6_m0 + 5 * kcomp + i); 

            auto tg_0_xxxyyy_0 = braBuffer.data(pidx_g_0_6_m0 + 6 * kcomp + i); 

            auto tg_0_xxxyyz_0 = braBuffer.data(pidx_g_0_6_m0 + 7 * kcomp + i); 

            auto tg_0_xxxyzz_0 = braBuffer.data(pidx_g_0_6_m0 + 8 * kcomp + i); 

            auto tg_0_xxxzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 9 * kcomp + i); 

            auto tg_0_xxyyyy_0 = braBuffer.data(pidx_g_0_6_m0 + 10 * kcomp + i); 

            auto tg_0_xxyyyz_0 = braBuffer.data(pidx_g_0_6_m0 + 11 * kcomp + i); 

            auto tg_0_xxyyzz_0 = braBuffer.data(pidx_g_0_6_m0 + 12 * kcomp + i); 

            auto tg_0_xxyzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 13 * kcomp + i); 

            auto tg_0_xxzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 14 * kcomp + i); 

            auto tg_0_xyyyyy_0 = braBuffer.data(pidx_g_0_6_m0 + 15 * kcomp + i); 

            auto tg_0_xyyyyz_0 = braBuffer.data(pidx_g_0_6_m0 + 16 * kcomp + i); 

            auto tg_0_xyyyzz_0 = braBuffer.data(pidx_g_0_6_m0 + 17 * kcomp + i); 

            auto tg_0_xyyzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 18 * kcomp + i); 

            auto tg_0_xyzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 19 * kcomp + i); 

            auto tg_0_xzzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 20 * kcomp + i); 

            auto tg_0_yyyyyy_0 = braBuffer.data(pidx_g_0_6_m0 + 21 * kcomp + i); 

            auto tg_0_yyyyyz_0 = braBuffer.data(pidx_g_0_6_m0 + 22 * kcomp + i); 

            auto tg_0_yyyyzz_0 = braBuffer.data(pidx_g_0_6_m0 + 23 * kcomp + i); 

            auto tg_0_yyyzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 24 * kcomp + i); 

            auto tg_0_yyzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 25 * kcomp + i); 

            auto tg_0_yzzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 26 * kcomp + i); 

            auto tg_0_zzzzzz_0 = braBuffer.data(pidx_g_0_6_m0 + 27 * kcomp + i); 

            auto tg_0_xxxxxxx_0 = braBuffer.data(pidx_g_0_7_m0 + i); 

            auto tg_0_xxxxxxy_0 = braBuffer.data(pidx_g_0_7_m0 + kcomp + i); 

            auto tg_0_xxxxxxz_0 = braBuffer.data(pidx_g_0_7_m0 + 2 * kcomp + i); 

            auto tg_0_xxxxxyy_0 = braBuffer.data(pidx_g_0_7_m0 + 3 * kcomp + i); 

            auto tg_0_xxxxxyz_0 = braBuffer.data(pidx_g_0_7_m0 + 4 * kcomp + i); 

            auto tg_0_xxxxxzz_0 = braBuffer.data(pidx_g_0_7_m0 + 5 * kcomp + i); 

            auto tg_0_xxxxyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 6 * kcomp + i); 

            auto tg_0_xxxxyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 7 * kcomp + i); 

            auto tg_0_xxxxyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 8 * kcomp + i); 

            auto tg_0_xxxxzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 9 * kcomp + i); 

            auto tg_0_xxxyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 10 * kcomp + i); 

            auto tg_0_xxxyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 11 * kcomp + i); 

            auto tg_0_xxxyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 12 * kcomp + i); 

            auto tg_0_xxxyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 13 * kcomp + i); 

            auto tg_0_xxxzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 14 * kcomp + i); 

            auto tg_0_xxyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 15 * kcomp + i); 

            auto tg_0_xxyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 16 * kcomp + i); 

            auto tg_0_xxyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 17 * kcomp + i); 

            auto tg_0_xxyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 18 * kcomp + i); 

            auto tg_0_xxyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 19 * kcomp + i); 

            auto tg_0_xxzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 20 * kcomp + i); 

            auto tg_0_xyyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 21 * kcomp + i); 

            auto tg_0_xyyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 22 * kcomp + i); 

            auto tg_0_xyyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 23 * kcomp + i); 

            auto tg_0_xyyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 24 * kcomp + i); 

            auto tg_0_xyyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 25 * kcomp + i); 

            auto tg_0_xyzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 26 * kcomp + i); 

            auto tg_0_xzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 27 * kcomp + i); 

            auto tg_0_yyyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 28 * kcomp + i); 

            auto tg_0_yyyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 29 * kcomp + i); 

            auto tg_0_yyyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 30 * kcomp + i); 

            auto tg_0_yyyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 31 * kcomp + i); 

            auto tg_0_yyyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 32 * kcomp + i); 

            auto tg_0_yyzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 33 * kcomp + i); 

            auto tg_0_yzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 34 * kcomp + i); 

            auto tg_0_zzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 35 * kcomp + i); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(tg_0_xxxxxx_0, tg_0_xxxxxxx_0, tg_0_xxxxxxy_0, tg_0_xxxxxxz_0, \
                                     tg_0_xxxxxy_0, tg_0_xxxxxyy_0, tg_0_xxxxxyz_0, tg_0_xxxxxz_0, tg_0_xxxxxzz_0, \
                                     tg_0_xxxxyy_0, tg_0_xxxxyyy_0, tg_0_xxxxyyz_0, tg_0_xxxxyz_0, tg_0_xxxxyzz_0, \
                                     tg_0_xxxxzz_0, tg_0_xxxxzzz_0, tg_0_xxxyyy_0, tg_0_xxxyyyy_0, tg_0_xxxyyyz_0, \
                                     tg_0_xxxyyz_0, tg_0_xxxyyzz_0, tg_0_xxxyzz_0, tg_0_xxxyzzz_0, tg_0_xxxzzz_0, \
                                     tg_0_xxxzzzz_0, tg_0_xxyyyy_0, tg_0_xxyyyyy_0, tg_0_xxyyyyz_0, tg_0_xxyyyz_0, \
                                     tg_0_xxyyyzz_0, tg_0_xxyyzz_0, tg_0_xxyyzzz_0, tg_0_xxyzzz_0, tg_0_xxyzzzz_0, \
                                     tg_0_xxzzzz_0, tg_0_xxzzzzz_0, tg_0_xyyyyy_0, tg_0_xyyyyyy_0, tg_0_xyyyyyz_0, \
                                     tg_0_xyyyyz_0, tg_0_xyyyyzz_0, tg_0_xyyyzz_0, tg_0_xyyyzzz_0, tg_0_xyyzzz_0, \
                                     tg_0_xyyzzzz_0, tg_0_xyzzzz_0, tg_0_xyzzzzz_0, tg_0_xzzzzz_0, tg_0_xzzzzzz_0, \
                                     tg_0_yyyyyy_0, tg_0_yyyyyyy_0, tg_0_yyyyyyz_0, tg_0_yyyyyz_0, tg_0_yyyyyzz_0, \
                                     tg_0_yyyyzz_0, tg_0_yyyyzzz_0, tg_0_yyyzzz_0, tg_0_yyyzzzz_0, tg_0_yyzzzz_0, \
                                     tg_0_yyzzzzz_0, tg_0_yzzzzz_0, tg_0_yzzzzzz_0, tg_0_zzzzzz_0, tg_0_zzzzzzz_0, \
                                     tg_x_xxxxxx_0, tg_x_xxxxxy_0, tg_x_xxxxxz_0, tg_x_xxxxyy_0, tg_x_xxxxyz_0, \
                                     tg_x_xxxxzz_0, tg_x_xxxyyy_0, tg_x_xxxyyz_0, tg_x_xxxyzz_0, tg_x_xxxzzz_0, \
                                     tg_x_xxyyyy_0, tg_x_xxyyyz_0, tg_x_xxyyzz_0, tg_x_xxyzzz_0, tg_x_xxzzzz_0, \
                                     tg_x_xyyyyy_0, tg_x_xyyyyz_0, tg_x_xyyyzz_0, tg_x_xyyzzz_0, tg_x_xyzzzz_0, \
                                     tg_x_xzzzzz_0, tg_x_yyyyyy_0, tg_x_yyyyyz_0, tg_x_yyyyzz_0, tg_x_yyyzzz_0, \
                                     tg_x_yyzzzz_0, tg_x_yzzzzz_0, tg_x_zzzzzz_0, tg_y_xxxxxx_0, tg_y_xxxxxy_0, \
                                     tg_y_xxxxxz_0, tg_y_xxxxyy_0, tg_y_xxxxyz_0, tg_y_xxxxzz_0, tg_y_xxxyyy_0, \
                                     tg_y_xxxyyz_0, tg_y_xxxyzz_0, tg_y_xxxzzz_0, tg_y_xxyyyy_0, tg_y_xxyyyz_0, \
                                     tg_y_xxyyzz_0, tg_y_xxyzzz_0, tg_y_xxzzzz_0, tg_y_xyyyyy_0, tg_y_xyyyyz_0, \
                                     tg_y_xyyyzz_0, tg_y_xyyzzz_0, tg_y_xyzzzz_0, tg_y_xzzzzz_0, tg_y_yyyyyy_0, \
                                     tg_y_yyyyyz_0, tg_y_yyyyzz_0, tg_y_yyyzzz_0, tg_y_yyzzzz_0, tg_y_yzzzzz_0, \
                                     tg_y_zzzzzz_0, tg_z_xxxxxx_0, tg_z_xxxxxy_0, tg_z_xxxxxz_0, tg_z_xxxxyy_0, \
                                     tg_z_xxxxyz_0, tg_z_xxxxzz_0, tg_z_xxxyyy_0, tg_z_xxxyyz_0, tg_z_xxxyzz_0, \
                                     tg_z_xxxzzz_0, tg_z_xxyyyy_0, tg_z_xxyyyz_0, tg_z_xxyyzz_0, tg_z_xxyzzz_0, \
                                     tg_z_xxzzzz_0, tg_z_xyyyyy_0, tg_z_xyyyyz_0, tg_z_xyyyzz_0, tg_z_xyyzzz_0, \
                                     tg_z_xyzzzz_0, tg_z_xzzzzz_0, tg_z_yyyyyy_0, tg_z_yyyyyz_0, tg_z_yyyyzz_0, \
                                     tg_z_yyyzzz_0, tg_z_yyzzzz_0, tg_z_yzzzzz_0, tg_z_zzzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_x_xxxxxx_0[j] = -ab_x * tg_0_xxxxxx_0[j] + tg_0_xxxxxxx_0[j];

                tg_x_xxxxxy_0[j] = -ab_x * tg_0_xxxxxy_0[j] + tg_0_xxxxxxy_0[j];

                tg_x_xxxxxz_0[j] = -ab_x * tg_0_xxxxxz_0[j] + tg_0_xxxxxxz_0[j];

                tg_x_xxxxyy_0[j] = -ab_x * tg_0_xxxxyy_0[j] + tg_0_xxxxxyy_0[j];

                tg_x_xxxxyz_0[j] = -ab_x * tg_0_xxxxyz_0[j] + tg_0_xxxxxyz_0[j];

                tg_x_xxxxzz_0[j] = -ab_x * tg_0_xxxxzz_0[j] + tg_0_xxxxxzz_0[j];

                tg_x_xxxyyy_0[j] = -ab_x * tg_0_xxxyyy_0[j] + tg_0_xxxxyyy_0[j];

                tg_x_xxxyyz_0[j] = -ab_x * tg_0_xxxyyz_0[j] + tg_0_xxxxyyz_0[j];

                tg_x_xxxyzz_0[j] = -ab_x * tg_0_xxxyzz_0[j] + tg_0_xxxxyzz_0[j];

                tg_x_xxxzzz_0[j] = -ab_x * tg_0_xxxzzz_0[j] + tg_0_xxxxzzz_0[j];

                tg_x_xxyyyy_0[j] = -ab_x * tg_0_xxyyyy_0[j] + tg_0_xxxyyyy_0[j];

                tg_x_xxyyyz_0[j] = -ab_x * tg_0_xxyyyz_0[j] + tg_0_xxxyyyz_0[j];

                tg_x_xxyyzz_0[j] = -ab_x * tg_0_xxyyzz_0[j] + tg_0_xxxyyzz_0[j];

                tg_x_xxyzzz_0[j] = -ab_x * tg_0_xxyzzz_0[j] + tg_0_xxxyzzz_0[j];

                tg_x_xxzzzz_0[j] = -ab_x * tg_0_xxzzzz_0[j] + tg_0_xxxzzzz_0[j];

                tg_x_xyyyyy_0[j] = -ab_x * tg_0_xyyyyy_0[j] + tg_0_xxyyyyy_0[j];

                tg_x_xyyyyz_0[j] = -ab_x * tg_0_xyyyyz_0[j] + tg_0_xxyyyyz_0[j];

                tg_x_xyyyzz_0[j] = -ab_x * tg_0_xyyyzz_0[j] + tg_0_xxyyyzz_0[j];

                tg_x_xyyzzz_0[j] = -ab_x * tg_0_xyyzzz_0[j] + tg_0_xxyyzzz_0[j];

                tg_x_xyzzzz_0[j] = -ab_x * tg_0_xyzzzz_0[j] + tg_0_xxyzzzz_0[j];

                tg_x_xzzzzz_0[j] = -ab_x * tg_0_xzzzzz_0[j] + tg_0_xxzzzzz_0[j];

                tg_x_yyyyyy_0[j] = -ab_x * tg_0_yyyyyy_0[j] + tg_0_xyyyyyy_0[j];

                tg_x_yyyyyz_0[j] = -ab_x * tg_0_yyyyyz_0[j] + tg_0_xyyyyyz_0[j];

                tg_x_yyyyzz_0[j] = -ab_x * tg_0_yyyyzz_0[j] + tg_0_xyyyyzz_0[j];

                tg_x_yyyzzz_0[j] = -ab_x * tg_0_yyyzzz_0[j] + tg_0_xyyyzzz_0[j];

                tg_x_yyzzzz_0[j] = -ab_x * tg_0_yyzzzz_0[j] + tg_0_xyyzzzz_0[j];

                tg_x_yzzzzz_0[j] = -ab_x * tg_0_yzzzzz_0[j] + tg_0_xyzzzzz_0[j];

                tg_x_zzzzzz_0[j] = -ab_x * tg_0_zzzzzz_0[j] + tg_0_xzzzzzz_0[j];

                tg_y_xxxxxx_0[j] = -ab_y * tg_0_xxxxxx_0[j] + tg_0_xxxxxxy_0[j];

                tg_y_xxxxxy_0[j] = -ab_y * tg_0_xxxxxy_0[j] + tg_0_xxxxxyy_0[j];

                tg_y_xxxxxz_0[j] = -ab_y * tg_0_xxxxxz_0[j] + tg_0_xxxxxyz_0[j];

                tg_y_xxxxyy_0[j] = -ab_y * tg_0_xxxxyy_0[j] + tg_0_xxxxyyy_0[j];

                tg_y_xxxxyz_0[j] = -ab_y * tg_0_xxxxyz_0[j] + tg_0_xxxxyyz_0[j];

                tg_y_xxxxzz_0[j] = -ab_y * tg_0_xxxxzz_0[j] + tg_0_xxxxyzz_0[j];

                tg_y_xxxyyy_0[j] = -ab_y * tg_0_xxxyyy_0[j] + tg_0_xxxyyyy_0[j];

                tg_y_xxxyyz_0[j] = -ab_y * tg_0_xxxyyz_0[j] + tg_0_xxxyyyz_0[j];

                tg_y_xxxyzz_0[j] = -ab_y * tg_0_xxxyzz_0[j] + tg_0_xxxyyzz_0[j];

                tg_y_xxxzzz_0[j] = -ab_y * tg_0_xxxzzz_0[j] + tg_0_xxxyzzz_0[j];

                tg_y_xxyyyy_0[j] = -ab_y * tg_0_xxyyyy_0[j] + tg_0_xxyyyyy_0[j];

                tg_y_xxyyyz_0[j] = -ab_y * tg_0_xxyyyz_0[j] + tg_0_xxyyyyz_0[j];

                tg_y_xxyyzz_0[j] = -ab_y * tg_0_xxyyzz_0[j] + tg_0_xxyyyzz_0[j];

                tg_y_xxyzzz_0[j] = -ab_y * tg_0_xxyzzz_0[j] + tg_0_xxyyzzz_0[j];

                tg_y_xxzzzz_0[j] = -ab_y * tg_0_xxzzzz_0[j] + tg_0_xxyzzzz_0[j];

                tg_y_xyyyyy_0[j] = -ab_y * tg_0_xyyyyy_0[j] + tg_0_xyyyyyy_0[j];

                tg_y_xyyyyz_0[j] = -ab_y * tg_0_xyyyyz_0[j] + tg_0_xyyyyyz_0[j];

                tg_y_xyyyzz_0[j] = -ab_y * tg_0_xyyyzz_0[j] + tg_0_xyyyyzz_0[j];

                tg_y_xyyzzz_0[j] = -ab_y * tg_0_xyyzzz_0[j] + tg_0_xyyyzzz_0[j];

                tg_y_xyzzzz_0[j] = -ab_y * tg_0_xyzzzz_0[j] + tg_0_xyyzzzz_0[j];

                tg_y_xzzzzz_0[j] = -ab_y * tg_0_xzzzzz_0[j] + tg_0_xyzzzzz_0[j];

                tg_y_yyyyyy_0[j] = -ab_y * tg_0_yyyyyy_0[j] + tg_0_yyyyyyy_0[j];

                tg_y_yyyyyz_0[j] = -ab_y * tg_0_yyyyyz_0[j] + tg_0_yyyyyyz_0[j];

                tg_y_yyyyzz_0[j] = -ab_y * tg_0_yyyyzz_0[j] + tg_0_yyyyyzz_0[j];

                tg_y_yyyzzz_0[j] = -ab_y * tg_0_yyyzzz_0[j] + tg_0_yyyyzzz_0[j];

                tg_y_yyzzzz_0[j] = -ab_y * tg_0_yyzzzz_0[j] + tg_0_yyyzzzz_0[j];

                tg_y_yzzzzz_0[j] = -ab_y * tg_0_yzzzzz_0[j] + tg_0_yyzzzzz_0[j];

                tg_y_zzzzzz_0[j] = -ab_y * tg_0_zzzzzz_0[j] + tg_0_yzzzzzz_0[j];

                tg_z_xxxxxx_0[j] = -ab_z * tg_0_xxxxxx_0[j] + tg_0_xxxxxxz_0[j];

                tg_z_xxxxxy_0[j] = -ab_z * tg_0_xxxxxy_0[j] + tg_0_xxxxxyz_0[j];

                tg_z_xxxxxz_0[j] = -ab_z * tg_0_xxxxxz_0[j] + tg_0_xxxxxzz_0[j];

                tg_z_xxxxyy_0[j] = -ab_z * tg_0_xxxxyy_0[j] + tg_0_xxxxyyz_0[j];

                tg_z_xxxxyz_0[j] = -ab_z * tg_0_xxxxyz_0[j] + tg_0_xxxxyzz_0[j];

                tg_z_xxxxzz_0[j] = -ab_z * tg_0_xxxxzz_0[j] + tg_0_xxxxzzz_0[j];

                tg_z_xxxyyy_0[j] = -ab_z * tg_0_xxxyyy_0[j] + tg_0_xxxyyyz_0[j];

                tg_z_xxxyyz_0[j] = -ab_z * tg_0_xxxyyz_0[j] + tg_0_xxxyyzz_0[j];

                tg_z_xxxyzz_0[j] = -ab_z * tg_0_xxxyzz_0[j] + tg_0_xxxyzzz_0[j];

                tg_z_xxxzzz_0[j] = -ab_z * tg_0_xxxzzz_0[j] + tg_0_xxxzzzz_0[j];

                tg_z_xxyyyy_0[j] = -ab_z * tg_0_xxyyyy_0[j] + tg_0_xxyyyyz_0[j];

                tg_z_xxyyyz_0[j] = -ab_z * tg_0_xxyyyz_0[j] + tg_0_xxyyyzz_0[j];

                tg_z_xxyyzz_0[j] = -ab_z * tg_0_xxyyzz_0[j] + tg_0_xxyyzzz_0[j];

                tg_z_xxyzzz_0[j] = -ab_z * tg_0_xxyzzz_0[j] + tg_0_xxyzzzz_0[j];

                tg_z_xxzzzz_0[j] = -ab_z * tg_0_xxzzzz_0[j] + tg_0_xxzzzzz_0[j];

                tg_z_xyyyyy_0[j] = -ab_z * tg_0_xyyyyy_0[j] + tg_0_xyyyyyz_0[j];

                tg_z_xyyyyz_0[j] = -ab_z * tg_0_xyyyyz_0[j] + tg_0_xyyyyzz_0[j];

                tg_z_xyyyzz_0[j] = -ab_z * tg_0_xyyyzz_0[j] + tg_0_xyyyzzz_0[j];

                tg_z_xyyzzz_0[j] = -ab_z * tg_0_xyyzzz_0[j] + tg_0_xyyzzzz_0[j];

                tg_z_xyzzzz_0[j] = -ab_z * tg_0_xyzzzz_0[j] + tg_0_xyzzzzz_0[j];

                tg_z_xzzzzz_0[j] = -ab_z * tg_0_xzzzzz_0[j] + tg_0_xzzzzzz_0[j];

                tg_z_yyyyyy_0[j] = -ab_z * tg_0_yyyyyy_0[j] + tg_0_yyyyyyz_0[j];

                tg_z_yyyyyz_0[j] = -ab_z * tg_0_yyyyyz_0[j] + tg_0_yyyyyzz_0[j];

                tg_z_yyyyzz_0[j] = -ab_z * tg_0_yyyyzz_0[j] + tg_0_yyyyzzz_0[j];

                tg_z_yyyzzz_0[j] = -ab_z * tg_0_yyyzzz_0[j] + tg_0_yyyzzzz_0[j];

                tg_z_yyzzzz_0[j] = -ab_z * tg_0_yyzzzz_0[j] + tg_0_yyzzzzz_0[j];

                tg_z_yzzzzz_0[j] = -ab_z * tg_0_yzzzzz_0[j] + tg_0_yzzzzzz_0[j];

                tg_z_zzzzzz_0[j] = -ab_z * tg_0_zzzzzz_0[j] + tg_0_zzzzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForPKXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        eribrrfunc::compElectronRepulsionForPKXY_0_54(braBuffer,
                                                      recursionMap,
                                                      abDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        eribrrfunc::compElectronRepulsionForPKXY_54_108(braBuffer,
                                                        recursionMap,
                                                        abDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForPKXY_0_54(      CMemBlock2D<double>& braBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& abDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,54)

        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        auto ab_y = (abDistances.data(1))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 7, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_1_7_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_0_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 7, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_0_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 8, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_0_xxxxxxx_0 = braBuffer.data(pidx_g_0_7_m0 + i); 

            auto tg_0_xxxxxxy_0 = braBuffer.data(pidx_g_0_7_m0 + kcomp + i); 

            auto tg_0_xxxxxxz_0 = braBuffer.data(pidx_g_0_7_m0 + 2 * kcomp + i); 

            auto tg_0_xxxxxyy_0 = braBuffer.data(pidx_g_0_7_m0 + 3 * kcomp + i); 

            auto tg_0_xxxxxyz_0 = braBuffer.data(pidx_g_0_7_m0 + 4 * kcomp + i); 

            auto tg_0_xxxxxzz_0 = braBuffer.data(pidx_g_0_7_m0 + 5 * kcomp + i); 

            auto tg_0_xxxxyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 6 * kcomp + i); 

            auto tg_0_xxxxyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 7 * kcomp + i); 

            auto tg_0_xxxxyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 8 * kcomp + i); 

            auto tg_0_xxxxzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 9 * kcomp + i); 

            auto tg_0_xxxyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 10 * kcomp + i); 

            auto tg_0_xxxyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 11 * kcomp + i); 

            auto tg_0_xxxyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 12 * kcomp + i); 

            auto tg_0_xxxyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 13 * kcomp + i); 

            auto tg_0_xxxzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 14 * kcomp + i); 

            auto tg_0_xxyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 15 * kcomp + i); 

            auto tg_0_xxyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 16 * kcomp + i); 

            auto tg_0_xxyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 17 * kcomp + i); 

            auto tg_0_xxyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 18 * kcomp + i); 

            auto tg_0_xxyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 19 * kcomp + i); 

            auto tg_0_xxzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 20 * kcomp + i); 

            auto tg_0_xyyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 21 * kcomp + i); 

            auto tg_0_xyyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 22 * kcomp + i); 

            auto tg_0_xyyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 23 * kcomp + i); 

            auto tg_0_xyyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 24 * kcomp + i); 

            auto tg_0_xyyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 25 * kcomp + i); 

            auto tg_0_xyzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 26 * kcomp + i); 

            auto tg_0_xzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 27 * kcomp + i); 

            auto tg_0_yyyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 28 * kcomp + i); 

            auto tg_0_yyyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 29 * kcomp + i); 

            auto tg_0_yyyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 30 * kcomp + i); 

            auto tg_0_yyyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 31 * kcomp + i); 

            auto tg_0_yyyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 32 * kcomp + i); 

            auto tg_0_yyzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 33 * kcomp + i); 

            auto tg_0_yzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 34 * kcomp + i); 

            auto tg_0_zzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 35 * kcomp + i); 

            auto tg_0_xxxxxxxx_0 = braBuffer.data(pidx_g_0_8_m0 + i); 

            auto tg_0_xxxxxxxy_0 = braBuffer.data(pidx_g_0_8_m0 + kcomp + i); 

            auto tg_0_xxxxxxxz_0 = braBuffer.data(pidx_g_0_8_m0 + 2 * kcomp + i); 

            auto tg_0_xxxxxxyy_0 = braBuffer.data(pidx_g_0_8_m0 + 3 * kcomp + i); 

            auto tg_0_xxxxxxyz_0 = braBuffer.data(pidx_g_0_8_m0 + 4 * kcomp + i); 

            auto tg_0_xxxxxxzz_0 = braBuffer.data(pidx_g_0_8_m0 + 5 * kcomp + i); 

            auto tg_0_xxxxxyyy_0 = braBuffer.data(pidx_g_0_8_m0 + 6 * kcomp + i); 

            auto tg_0_xxxxxyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 7 * kcomp + i); 

            auto tg_0_xxxxxyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 8 * kcomp + i); 

            auto tg_0_xxxxxzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 9 * kcomp + i); 

            auto tg_0_xxxxyyyy_0 = braBuffer.data(pidx_g_0_8_m0 + 10 * kcomp + i); 

            auto tg_0_xxxxyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 11 * kcomp + i); 

            auto tg_0_xxxxyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 12 * kcomp + i); 

            auto tg_0_xxxxyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 13 * kcomp + i); 

            auto tg_0_xxxxzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 14 * kcomp + i); 

            auto tg_0_xxxyyyyy_0 = braBuffer.data(pidx_g_0_8_m0 + 15 * kcomp + i); 

            auto tg_0_xxxyyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 16 * kcomp + i); 

            auto tg_0_xxxyyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 17 * kcomp + i); 

            auto tg_0_xxxyyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 18 * kcomp + i); 

            auto tg_0_xxxyzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 19 * kcomp + i); 

            auto tg_0_xxxzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 20 * kcomp + i); 

            auto tg_0_xxyyyyyy_0 = braBuffer.data(pidx_g_0_8_m0 + 21 * kcomp + i); 

            auto tg_0_xxyyyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 22 * kcomp + i); 

            auto tg_0_xxyyyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 23 * kcomp + i); 

            auto tg_0_xxyyyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 24 * kcomp + i); 

            auto tg_0_xxyyzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 25 * kcomp + i); 

            auto tg_0_xxyzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 26 * kcomp + i); 

            auto tg_0_xxzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 27 * kcomp + i); 

            auto tg_0_xyyyyyyy_0 = braBuffer.data(pidx_g_0_8_m0 + 28 * kcomp + i); 

            auto tg_0_xyyyyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 29 * kcomp + i); 

            auto tg_0_xyyyyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 30 * kcomp + i); 

            auto tg_0_xyyyyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 31 * kcomp + i); 

            auto tg_0_xyyyzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 32 * kcomp + i); 

            auto tg_0_xyyzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 33 * kcomp + i); 

            auto tg_0_xyzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 34 * kcomp + i); 

            auto tg_0_xzzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 35 * kcomp + i); 

            // set up pointers to integrals

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

            auto tg_x_yyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 28 * kcomp + i); 

            auto tg_x_yyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 29 * kcomp + i); 

            auto tg_x_yyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 30 * kcomp + i); 

            auto tg_x_yyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 31 * kcomp + i); 

            auto tg_x_yyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 32 * kcomp + i); 

            auto tg_x_yyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 33 * kcomp + i); 

            auto tg_x_yzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 34 * kcomp + i); 

            auto tg_x_zzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 35 * kcomp + i); 

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

            // Batch of Integrals (0,54)

            #pragma omp simd aligned(tg_0_xxxxxxx_0, tg_0_xxxxxxxx_0, tg_0_xxxxxxxy_0, tg_0_xxxxxxxz_0, \
                                     tg_0_xxxxxxy_0, tg_0_xxxxxxyy_0, tg_0_xxxxxxyz_0, tg_0_xxxxxxz_0, tg_0_xxxxxxzz_0, \
                                     tg_0_xxxxxyy_0, tg_0_xxxxxyyy_0, tg_0_xxxxxyyz_0, tg_0_xxxxxyz_0, tg_0_xxxxxyzz_0, \
                                     tg_0_xxxxxzz_0, tg_0_xxxxxzzz_0, tg_0_xxxxyyy_0, tg_0_xxxxyyyy_0, tg_0_xxxxyyyz_0, \
                                     tg_0_xxxxyyz_0, tg_0_xxxxyyzz_0, tg_0_xxxxyzz_0, tg_0_xxxxyzzz_0, tg_0_xxxxzzz_0, \
                                     tg_0_xxxxzzzz_0, tg_0_xxxyyyy_0, tg_0_xxxyyyyy_0, tg_0_xxxyyyyz_0, tg_0_xxxyyyz_0, \
                                     tg_0_xxxyyyzz_0, tg_0_xxxyyzz_0, tg_0_xxxyyzzz_0, tg_0_xxxyzzz_0, tg_0_xxxyzzzz_0, \
                                     tg_0_xxxzzzz_0, tg_0_xxxzzzzz_0, tg_0_xxyyyyy_0, tg_0_xxyyyyyy_0, tg_0_xxyyyyyz_0, \
                                     tg_0_xxyyyyz_0, tg_0_xxyyyyzz_0, tg_0_xxyyyzz_0, tg_0_xxyyyzzz_0, tg_0_xxyyzzz_0, \
                                     tg_0_xxyyzzzz_0, tg_0_xxyzzzz_0, tg_0_xxyzzzzz_0, tg_0_xxzzzzz_0, tg_0_xxzzzzzz_0, \
                                     tg_0_xyyyyyy_0, tg_0_xyyyyyyy_0, tg_0_xyyyyyyz_0, tg_0_xyyyyyz_0, tg_0_xyyyyyzz_0, \
                                     tg_0_xyyyyzz_0, tg_0_xyyyyzzz_0, tg_0_xyyyzzz_0, tg_0_xyyyzzzz_0, tg_0_xyyzzzz_0, \
                                     tg_0_xyyzzzzz_0, tg_0_xyzzzzz_0, tg_0_xyzzzzzz_0, tg_0_xzzzzzz_0, tg_0_xzzzzzzz_0, \
                                     tg_0_yyyyyyy_0, tg_0_yyyyyyz_0, tg_0_yyyyyzz_0, tg_0_yyyyzzz_0, tg_0_yyyzzzz_0, \
                                     tg_0_yyzzzzz_0, tg_0_yzzzzzz_0, tg_0_zzzzzzz_0, tg_x_xxxxxxx_0, tg_x_xxxxxxy_0, \
                                     tg_x_xxxxxxz_0, tg_x_xxxxxyy_0, tg_x_xxxxxyz_0, tg_x_xxxxxzz_0, tg_x_xxxxyyy_0, \
                                     tg_x_xxxxyyz_0, tg_x_xxxxyzz_0, tg_x_xxxxzzz_0, tg_x_xxxyyyy_0, tg_x_xxxyyyz_0, \
                                     tg_x_xxxyyzz_0, tg_x_xxxyzzz_0, tg_x_xxxzzzz_0, tg_x_xxyyyyy_0, tg_x_xxyyyyz_0, \
                                     tg_x_xxyyyzz_0, tg_x_xxyyzzz_0, tg_x_xxyzzzz_0, tg_x_xxzzzzz_0, tg_x_xyyyyyy_0, \
                                     tg_x_xyyyyyz_0, tg_x_xyyyyzz_0, tg_x_xyyyzzz_0, tg_x_xyyzzzz_0, tg_x_xyzzzzz_0, \
                                     tg_x_xzzzzzz_0, tg_x_yyyyyyy_0, tg_x_yyyyyyz_0, tg_x_yyyyyzz_0, tg_x_yyyyzzz_0, \
                                     tg_x_yyyzzzz_0, tg_x_yyzzzzz_0, tg_x_yzzzzzz_0, tg_x_zzzzzzz_0, tg_y_xxxxxxx_0, \
                                     tg_y_xxxxxxy_0, tg_y_xxxxxxz_0, tg_y_xxxxxyy_0, tg_y_xxxxxyz_0, tg_y_xxxxxzz_0, \
                                     tg_y_xxxxyyy_0, tg_y_xxxxyyz_0, tg_y_xxxxyzz_0, tg_y_xxxxzzz_0, tg_y_xxxyyyy_0, \
                                     tg_y_xxxyyyz_0, tg_y_xxxyyzz_0, tg_y_xxxyzzz_0, tg_y_xxxzzzz_0, tg_y_xxyyyyy_0, \
                                     tg_y_xxyyyyz_0, tg_y_xxyyyzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_x_xxxxxxx_0[j] = -ab_x * tg_0_xxxxxxx_0[j] + tg_0_xxxxxxxx_0[j];

                tg_x_xxxxxxy_0[j] = -ab_x * tg_0_xxxxxxy_0[j] + tg_0_xxxxxxxy_0[j];

                tg_x_xxxxxxz_0[j] = -ab_x * tg_0_xxxxxxz_0[j] + tg_0_xxxxxxxz_0[j];

                tg_x_xxxxxyy_0[j] = -ab_x * tg_0_xxxxxyy_0[j] + tg_0_xxxxxxyy_0[j];

                tg_x_xxxxxyz_0[j] = -ab_x * tg_0_xxxxxyz_0[j] + tg_0_xxxxxxyz_0[j];

                tg_x_xxxxxzz_0[j] = -ab_x * tg_0_xxxxxzz_0[j] + tg_0_xxxxxxzz_0[j];

                tg_x_xxxxyyy_0[j] = -ab_x * tg_0_xxxxyyy_0[j] + tg_0_xxxxxyyy_0[j];

                tg_x_xxxxyyz_0[j] = -ab_x * tg_0_xxxxyyz_0[j] + tg_0_xxxxxyyz_0[j];

                tg_x_xxxxyzz_0[j] = -ab_x * tg_0_xxxxyzz_0[j] + tg_0_xxxxxyzz_0[j];

                tg_x_xxxxzzz_0[j] = -ab_x * tg_0_xxxxzzz_0[j] + tg_0_xxxxxzzz_0[j];

                tg_x_xxxyyyy_0[j] = -ab_x * tg_0_xxxyyyy_0[j] + tg_0_xxxxyyyy_0[j];

                tg_x_xxxyyyz_0[j] = -ab_x * tg_0_xxxyyyz_0[j] + tg_0_xxxxyyyz_0[j];

                tg_x_xxxyyzz_0[j] = -ab_x * tg_0_xxxyyzz_0[j] + tg_0_xxxxyyzz_0[j];

                tg_x_xxxyzzz_0[j] = -ab_x * tg_0_xxxyzzz_0[j] + tg_0_xxxxyzzz_0[j];

                tg_x_xxxzzzz_0[j] = -ab_x * tg_0_xxxzzzz_0[j] + tg_0_xxxxzzzz_0[j];

                tg_x_xxyyyyy_0[j] = -ab_x * tg_0_xxyyyyy_0[j] + tg_0_xxxyyyyy_0[j];

                tg_x_xxyyyyz_0[j] = -ab_x * tg_0_xxyyyyz_0[j] + tg_0_xxxyyyyz_0[j];

                tg_x_xxyyyzz_0[j] = -ab_x * tg_0_xxyyyzz_0[j] + tg_0_xxxyyyzz_0[j];

                tg_x_xxyyzzz_0[j] = -ab_x * tg_0_xxyyzzz_0[j] + tg_0_xxxyyzzz_0[j];

                tg_x_xxyzzzz_0[j] = -ab_x * tg_0_xxyzzzz_0[j] + tg_0_xxxyzzzz_0[j];

                tg_x_xxzzzzz_0[j] = -ab_x * tg_0_xxzzzzz_0[j] + tg_0_xxxzzzzz_0[j];

                tg_x_xyyyyyy_0[j] = -ab_x * tg_0_xyyyyyy_0[j] + tg_0_xxyyyyyy_0[j];

                tg_x_xyyyyyz_0[j] = -ab_x * tg_0_xyyyyyz_0[j] + tg_0_xxyyyyyz_0[j];

                tg_x_xyyyyzz_0[j] = -ab_x * tg_0_xyyyyzz_0[j] + tg_0_xxyyyyzz_0[j];

                tg_x_xyyyzzz_0[j] = -ab_x * tg_0_xyyyzzz_0[j] + tg_0_xxyyyzzz_0[j];

                tg_x_xyyzzzz_0[j] = -ab_x * tg_0_xyyzzzz_0[j] + tg_0_xxyyzzzz_0[j];

                tg_x_xyzzzzz_0[j] = -ab_x * tg_0_xyzzzzz_0[j] + tg_0_xxyzzzzz_0[j];

                tg_x_xzzzzzz_0[j] = -ab_x * tg_0_xzzzzzz_0[j] + tg_0_xxzzzzzz_0[j];

                tg_x_yyyyyyy_0[j] = -ab_x * tg_0_yyyyyyy_0[j] + tg_0_xyyyyyyy_0[j];

                tg_x_yyyyyyz_0[j] = -ab_x * tg_0_yyyyyyz_0[j] + tg_0_xyyyyyyz_0[j];

                tg_x_yyyyyzz_0[j] = -ab_x * tg_0_yyyyyzz_0[j] + tg_0_xyyyyyzz_0[j];

                tg_x_yyyyzzz_0[j] = -ab_x * tg_0_yyyyzzz_0[j] + tg_0_xyyyyzzz_0[j];

                tg_x_yyyzzzz_0[j] = -ab_x * tg_0_yyyzzzz_0[j] + tg_0_xyyyzzzz_0[j];

                tg_x_yyzzzzz_0[j] = -ab_x * tg_0_yyzzzzz_0[j] + tg_0_xyyzzzzz_0[j];

                tg_x_yzzzzzz_0[j] = -ab_x * tg_0_yzzzzzz_0[j] + tg_0_xyzzzzzz_0[j];

                tg_x_zzzzzzz_0[j] = -ab_x * tg_0_zzzzzzz_0[j] + tg_0_xzzzzzzz_0[j];

                tg_y_xxxxxxx_0[j] = -ab_y * tg_0_xxxxxxx_0[j] + tg_0_xxxxxxxy_0[j];

                tg_y_xxxxxxy_0[j] = -ab_y * tg_0_xxxxxxy_0[j] + tg_0_xxxxxxyy_0[j];

                tg_y_xxxxxxz_0[j] = -ab_y * tg_0_xxxxxxz_0[j] + tg_0_xxxxxxyz_0[j];

                tg_y_xxxxxyy_0[j] = -ab_y * tg_0_xxxxxyy_0[j] + tg_0_xxxxxyyy_0[j];

                tg_y_xxxxxyz_0[j] = -ab_y * tg_0_xxxxxyz_0[j] + tg_0_xxxxxyyz_0[j];

                tg_y_xxxxxzz_0[j] = -ab_y * tg_0_xxxxxzz_0[j] + tg_0_xxxxxyzz_0[j];

                tg_y_xxxxyyy_0[j] = -ab_y * tg_0_xxxxyyy_0[j] + tg_0_xxxxyyyy_0[j];

                tg_y_xxxxyyz_0[j] = -ab_y * tg_0_xxxxyyz_0[j] + tg_0_xxxxyyyz_0[j];

                tg_y_xxxxyzz_0[j] = -ab_y * tg_0_xxxxyzz_0[j] + tg_0_xxxxyyzz_0[j];

                tg_y_xxxxzzz_0[j] = -ab_y * tg_0_xxxxzzz_0[j] + tg_0_xxxxyzzz_0[j];

                tg_y_xxxyyyy_0[j] = -ab_y * tg_0_xxxyyyy_0[j] + tg_0_xxxyyyyy_0[j];

                tg_y_xxxyyyz_0[j] = -ab_y * tg_0_xxxyyyz_0[j] + tg_0_xxxyyyyz_0[j];

                tg_y_xxxyyzz_0[j] = -ab_y * tg_0_xxxyyzz_0[j] + tg_0_xxxyyyzz_0[j];

                tg_y_xxxyzzz_0[j] = -ab_y * tg_0_xxxyzzz_0[j] + tg_0_xxxyyzzz_0[j];

                tg_y_xxxzzzz_0[j] = -ab_y * tg_0_xxxzzzz_0[j] + tg_0_xxxyzzzz_0[j];

                tg_y_xxyyyyy_0[j] = -ab_y * tg_0_xxyyyyy_0[j] + tg_0_xxyyyyyy_0[j];

                tg_y_xxyyyyz_0[j] = -ab_y * tg_0_xxyyyyz_0[j] + tg_0_xxyyyyyz_0[j];

                tg_y_xxyyyzz_0[j] = -ab_y * tg_0_xxyyyzz_0[j] + tg_0_xxyyyyzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForPKXY_54_108(      CMemBlock2D<double>& braBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& abDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetContrPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (54,108)

        // set up distances R(AB) = A - B

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {1, 7, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_1_7_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_0_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 7, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_0_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {0, 8, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

            auto tg_0_xxxxxxx_0 = braBuffer.data(pidx_g_0_7_m0 + i); 

            auto tg_0_xxxxxxy_0 = braBuffer.data(pidx_g_0_7_m0 + kcomp + i); 

            auto tg_0_xxxxxxz_0 = braBuffer.data(pidx_g_0_7_m0 + 2 * kcomp + i); 

            auto tg_0_xxxxxyy_0 = braBuffer.data(pidx_g_0_7_m0 + 3 * kcomp + i); 

            auto tg_0_xxxxxyz_0 = braBuffer.data(pidx_g_0_7_m0 + 4 * kcomp + i); 

            auto tg_0_xxxxxzz_0 = braBuffer.data(pidx_g_0_7_m0 + 5 * kcomp + i); 

            auto tg_0_xxxxyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 6 * kcomp + i); 

            auto tg_0_xxxxyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 7 * kcomp + i); 

            auto tg_0_xxxxyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 8 * kcomp + i); 

            auto tg_0_xxxxzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 9 * kcomp + i); 

            auto tg_0_xxxyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 10 * kcomp + i); 

            auto tg_0_xxxyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 11 * kcomp + i); 

            auto tg_0_xxxyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 12 * kcomp + i); 

            auto tg_0_xxxyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 13 * kcomp + i); 

            auto tg_0_xxxzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 14 * kcomp + i); 

            auto tg_0_xxyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 15 * kcomp + i); 

            auto tg_0_xxyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 16 * kcomp + i); 

            auto tg_0_xxyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 17 * kcomp + i); 

            auto tg_0_xxyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 18 * kcomp + i); 

            auto tg_0_xxyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 19 * kcomp + i); 

            auto tg_0_xxzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 20 * kcomp + i); 

            auto tg_0_xyyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 21 * kcomp + i); 

            auto tg_0_xyyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 22 * kcomp + i); 

            auto tg_0_xyyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 23 * kcomp + i); 

            auto tg_0_xyyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 24 * kcomp + i); 

            auto tg_0_xyyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 25 * kcomp + i); 

            auto tg_0_xyzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 26 * kcomp + i); 

            auto tg_0_xzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 27 * kcomp + i); 

            auto tg_0_yyyyyyy_0 = braBuffer.data(pidx_g_0_7_m0 + 28 * kcomp + i); 

            auto tg_0_yyyyyyz_0 = braBuffer.data(pidx_g_0_7_m0 + 29 * kcomp + i); 

            auto tg_0_yyyyyzz_0 = braBuffer.data(pidx_g_0_7_m0 + 30 * kcomp + i); 

            auto tg_0_yyyyzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 31 * kcomp + i); 

            auto tg_0_yyyzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 32 * kcomp + i); 

            auto tg_0_yyzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 33 * kcomp + i); 

            auto tg_0_yzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 34 * kcomp + i); 

            auto tg_0_zzzzzzz_0 = braBuffer.data(pidx_g_0_7_m0 + 35 * kcomp + i); 

            auto tg_0_xxxxxxxz_0 = braBuffer.data(pidx_g_0_8_m0 + 2 * kcomp + i); 

            auto tg_0_xxxxxxyz_0 = braBuffer.data(pidx_g_0_8_m0 + 4 * kcomp + i); 

            auto tg_0_xxxxxxzz_0 = braBuffer.data(pidx_g_0_8_m0 + 5 * kcomp + i); 

            auto tg_0_xxxxxyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 7 * kcomp + i); 

            auto tg_0_xxxxxyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 8 * kcomp + i); 

            auto tg_0_xxxxxzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 9 * kcomp + i); 

            auto tg_0_xxxxyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 11 * kcomp + i); 

            auto tg_0_xxxxyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 12 * kcomp + i); 

            auto tg_0_xxxxyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 13 * kcomp + i); 

            auto tg_0_xxxxzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 14 * kcomp + i); 

            auto tg_0_xxxyyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 16 * kcomp + i); 

            auto tg_0_xxxyyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 17 * kcomp + i); 

            auto tg_0_xxxyyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 18 * kcomp + i); 

            auto tg_0_xxxyzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 19 * kcomp + i); 

            auto tg_0_xxxzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 20 * kcomp + i); 

            auto tg_0_xxyyyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 22 * kcomp + i); 

            auto tg_0_xxyyyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 23 * kcomp + i); 

            auto tg_0_xxyyyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 24 * kcomp + i); 

            auto tg_0_xxyyzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 25 * kcomp + i); 

            auto tg_0_xxyzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 26 * kcomp + i); 

            auto tg_0_xxzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 27 * kcomp + i); 

            auto tg_0_xyyyyyyy_0 = braBuffer.data(pidx_g_0_8_m0 + 28 * kcomp + i); 

            auto tg_0_xyyyyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 29 * kcomp + i); 

            auto tg_0_xyyyyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 30 * kcomp + i); 

            auto tg_0_xyyyyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 31 * kcomp + i); 

            auto tg_0_xyyyzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 32 * kcomp + i); 

            auto tg_0_xyyzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 33 * kcomp + i); 

            auto tg_0_xyzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 34 * kcomp + i); 

            auto tg_0_xzzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 35 * kcomp + i); 

            auto tg_0_yyyyyyyy_0 = braBuffer.data(pidx_g_0_8_m0 + 36 * kcomp + i); 

            auto tg_0_yyyyyyyz_0 = braBuffer.data(pidx_g_0_8_m0 + 37 * kcomp + i); 

            auto tg_0_yyyyyyzz_0 = braBuffer.data(pidx_g_0_8_m0 + 38 * kcomp + i); 

            auto tg_0_yyyyyzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 39 * kcomp + i); 

            auto tg_0_yyyyzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 40 * kcomp + i); 

            auto tg_0_yyyzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 41 * kcomp + i); 

            auto tg_0_yyzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 42 * kcomp + i); 

            auto tg_0_yzzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 43 * kcomp + i); 

            auto tg_0_zzzzzzzz_0 = braBuffer.data(pidx_g_0_8_m0 + 44 * kcomp + i); 

            // set up pointers to integrals

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

            auto tg_y_yyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 64 * kcomp + i); 

            auto tg_y_yyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 65 * kcomp + i); 

            auto tg_y_yyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 66 * kcomp + i); 

            auto tg_y_yyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 67 * kcomp + i); 

            auto tg_y_yyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 68 * kcomp + i); 

            auto tg_y_yyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 69 * kcomp + i); 

            auto tg_y_yzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 70 * kcomp + i); 

            auto tg_y_zzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 71 * kcomp + i); 

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

            auto tg_z_yyyyyyy_0 = braBuffer.data(pidx_g_1_7_m0 + 100 * kcomp + i); 

            auto tg_z_yyyyyyz_0 = braBuffer.data(pidx_g_1_7_m0 + 101 * kcomp + i); 

            auto tg_z_yyyyyzz_0 = braBuffer.data(pidx_g_1_7_m0 + 102 * kcomp + i); 

            auto tg_z_yyyyzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 103 * kcomp + i); 

            auto tg_z_yyyzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 104 * kcomp + i); 

            auto tg_z_yyzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 105 * kcomp + i); 

            auto tg_z_yzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 106 * kcomp + i); 

            auto tg_z_zzzzzzz_0 = braBuffer.data(pidx_g_1_7_m0 + 107 * kcomp + i); 

            // Batch of Integrals (54,108)

            #pragma omp simd aligned(tg_0_xxxxxxx_0, tg_0_xxxxxxxz_0, tg_0_xxxxxxy_0, tg_0_xxxxxxyz_0, \
                                     tg_0_xxxxxxz_0, tg_0_xxxxxxzz_0, tg_0_xxxxxyy_0, tg_0_xxxxxyyz_0, tg_0_xxxxxyz_0, \
                                     tg_0_xxxxxyzz_0, tg_0_xxxxxzz_0, tg_0_xxxxxzzz_0, tg_0_xxxxyyy_0, tg_0_xxxxyyyz_0, \
                                     tg_0_xxxxyyz_0, tg_0_xxxxyyzz_0, tg_0_xxxxyzz_0, tg_0_xxxxyzzz_0, tg_0_xxxxzzz_0, \
                                     tg_0_xxxxzzzz_0, tg_0_xxxyyyy_0, tg_0_xxxyyyyz_0, tg_0_xxxyyyz_0, tg_0_xxxyyyzz_0, \
                                     tg_0_xxxyyzz_0, tg_0_xxxyyzzz_0, tg_0_xxxyzzz_0, tg_0_xxxyzzzz_0, tg_0_xxxzzzz_0, \
                                     tg_0_xxxzzzzz_0, tg_0_xxyyyyy_0, tg_0_xxyyyyyz_0, tg_0_xxyyyyz_0, tg_0_xxyyyyzz_0, \
                                     tg_0_xxyyyzz_0, tg_0_xxyyyzzz_0, tg_0_xxyyzzz_0, tg_0_xxyyzzzz_0, tg_0_xxyzzzz_0, \
                                     tg_0_xxyzzzzz_0, tg_0_xxzzzzz_0, tg_0_xxzzzzzz_0, tg_0_xyyyyyy_0, tg_0_xyyyyyyy_0, \
                                     tg_0_xyyyyyyz_0, tg_0_xyyyyyz_0, tg_0_xyyyyyzz_0, tg_0_xyyyyzz_0, tg_0_xyyyyzzz_0, \
                                     tg_0_xyyyzzz_0, tg_0_xyyyzzzz_0, tg_0_xyyzzzz_0, tg_0_xyyzzzzz_0, tg_0_xyzzzzz_0, \
                                     tg_0_xyzzzzzz_0, tg_0_xzzzzzz_0, tg_0_xzzzzzzz_0, tg_0_yyyyyyy_0, tg_0_yyyyyyyy_0, \
                                     tg_0_yyyyyyyz_0, tg_0_yyyyyyz_0, tg_0_yyyyyyzz_0, tg_0_yyyyyzz_0, tg_0_yyyyyzzz_0, \
                                     tg_0_yyyyzzz_0, tg_0_yyyyzzzz_0, tg_0_yyyzzzz_0, tg_0_yyyzzzzz_0, tg_0_yyzzzzz_0, \
                                     tg_0_yyzzzzzz_0, tg_0_yzzzzzz_0, tg_0_yzzzzzzz_0, tg_0_zzzzzzz_0, tg_0_zzzzzzzz_0, \
                                     tg_y_xxyyzzz_0, tg_y_xxyzzzz_0, tg_y_xxzzzzz_0, tg_y_xyyyyyy_0, tg_y_xyyyyyz_0, \
                                     tg_y_xyyyyzz_0, tg_y_xyyyzzz_0, tg_y_xyyzzzz_0, tg_y_xyzzzzz_0, tg_y_xzzzzzz_0, \
                                     tg_y_yyyyyyy_0, tg_y_yyyyyyz_0, tg_y_yyyyyzz_0, tg_y_yyyyzzz_0, tg_y_yyyzzzz_0, \
                                     tg_y_yyzzzzz_0, tg_y_yzzzzzz_0, tg_y_zzzzzzz_0, tg_z_xxxxxxx_0, tg_z_xxxxxxy_0, \
                                     tg_z_xxxxxxz_0, tg_z_xxxxxyy_0, tg_z_xxxxxyz_0, tg_z_xxxxxzz_0, tg_z_xxxxyyy_0, \
                                     tg_z_xxxxyyz_0, tg_z_xxxxyzz_0, tg_z_xxxxzzz_0, tg_z_xxxyyyy_0, tg_z_xxxyyyz_0, \
                                     tg_z_xxxyyzz_0, tg_z_xxxyzzz_0, tg_z_xxxzzzz_0, tg_z_xxyyyyy_0, tg_z_xxyyyyz_0, \
                                     tg_z_xxyyyzz_0, tg_z_xxyyzzz_0, tg_z_xxyzzzz_0, tg_z_xxzzzzz_0, tg_z_xyyyyyy_0, \
                                     tg_z_xyyyyyz_0, tg_z_xyyyyzz_0, tg_z_xyyyzzz_0, tg_z_xyyzzzz_0, tg_z_xyzzzzz_0, \
                                     tg_z_xzzzzzz_0, tg_z_yyyyyyy_0, tg_z_yyyyyyz_0, tg_z_yyyyyzz_0, tg_z_yyyyzzz_0, \
                                     tg_z_yyyzzzz_0, tg_z_yyzzzzz_0, tg_z_yzzzzzz_0, tg_z_zzzzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_y_xxyyzzz_0[j] = -ab_y * tg_0_xxyyzzz_0[j] + tg_0_xxyyyzzz_0[j];

                tg_y_xxyzzzz_0[j] = -ab_y * tg_0_xxyzzzz_0[j] + tg_0_xxyyzzzz_0[j];

                tg_y_xxzzzzz_0[j] = -ab_y * tg_0_xxzzzzz_0[j] + tg_0_xxyzzzzz_0[j];

                tg_y_xyyyyyy_0[j] = -ab_y * tg_0_xyyyyyy_0[j] + tg_0_xyyyyyyy_0[j];

                tg_y_xyyyyyz_0[j] = -ab_y * tg_0_xyyyyyz_0[j] + tg_0_xyyyyyyz_0[j];

                tg_y_xyyyyzz_0[j] = -ab_y * tg_0_xyyyyzz_0[j] + tg_0_xyyyyyzz_0[j];

                tg_y_xyyyzzz_0[j] = -ab_y * tg_0_xyyyzzz_0[j] + tg_0_xyyyyzzz_0[j];

                tg_y_xyyzzzz_0[j] = -ab_y * tg_0_xyyzzzz_0[j] + tg_0_xyyyzzzz_0[j];

                tg_y_xyzzzzz_0[j] = -ab_y * tg_0_xyzzzzz_0[j] + tg_0_xyyzzzzz_0[j];

                tg_y_xzzzzzz_0[j] = -ab_y * tg_0_xzzzzzz_0[j] + tg_0_xyzzzzzz_0[j];

                tg_y_yyyyyyy_0[j] = -ab_y * tg_0_yyyyyyy_0[j] + tg_0_yyyyyyyy_0[j];

                tg_y_yyyyyyz_0[j] = -ab_y * tg_0_yyyyyyz_0[j] + tg_0_yyyyyyyz_0[j];

                tg_y_yyyyyzz_0[j] = -ab_y * tg_0_yyyyyzz_0[j] + tg_0_yyyyyyzz_0[j];

                tg_y_yyyyzzz_0[j] = -ab_y * tg_0_yyyyzzz_0[j] + tg_0_yyyyyzzz_0[j];

                tg_y_yyyzzzz_0[j] = -ab_y * tg_0_yyyzzzz_0[j] + tg_0_yyyyzzzz_0[j];

                tg_y_yyzzzzz_0[j] = -ab_y * tg_0_yyzzzzz_0[j] + tg_0_yyyzzzzz_0[j];

                tg_y_yzzzzzz_0[j] = -ab_y * tg_0_yzzzzzz_0[j] + tg_0_yyzzzzzz_0[j];

                tg_y_zzzzzzz_0[j] = -ab_y * tg_0_zzzzzzz_0[j] + tg_0_yzzzzzzz_0[j];

                tg_z_xxxxxxx_0[j] = -ab_z * tg_0_xxxxxxx_0[j] + tg_0_xxxxxxxz_0[j];

                tg_z_xxxxxxy_0[j] = -ab_z * tg_0_xxxxxxy_0[j] + tg_0_xxxxxxyz_0[j];

                tg_z_xxxxxxz_0[j] = -ab_z * tg_0_xxxxxxz_0[j] + tg_0_xxxxxxzz_0[j];

                tg_z_xxxxxyy_0[j] = -ab_z * tg_0_xxxxxyy_0[j] + tg_0_xxxxxyyz_0[j];

                tg_z_xxxxxyz_0[j] = -ab_z * tg_0_xxxxxyz_0[j] + tg_0_xxxxxyzz_0[j];

                tg_z_xxxxxzz_0[j] = -ab_z * tg_0_xxxxxzz_0[j] + tg_0_xxxxxzzz_0[j];

                tg_z_xxxxyyy_0[j] = -ab_z * tg_0_xxxxyyy_0[j] + tg_0_xxxxyyyz_0[j];

                tg_z_xxxxyyz_0[j] = -ab_z * tg_0_xxxxyyz_0[j] + tg_0_xxxxyyzz_0[j];

                tg_z_xxxxyzz_0[j] = -ab_z * tg_0_xxxxyzz_0[j] + tg_0_xxxxyzzz_0[j];

                tg_z_xxxxzzz_0[j] = -ab_z * tg_0_xxxxzzz_0[j] + tg_0_xxxxzzzz_0[j];

                tg_z_xxxyyyy_0[j] = -ab_z * tg_0_xxxyyyy_0[j] + tg_0_xxxyyyyz_0[j];

                tg_z_xxxyyyz_0[j] = -ab_z * tg_0_xxxyyyz_0[j] + tg_0_xxxyyyzz_0[j];

                tg_z_xxxyyzz_0[j] = -ab_z * tg_0_xxxyyzz_0[j] + tg_0_xxxyyzzz_0[j];

                tg_z_xxxyzzz_0[j] = -ab_z * tg_0_xxxyzzz_0[j] + tg_0_xxxyzzzz_0[j];

                tg_z_xxxzzzz_0[j] = -ab_z * tg_0_xxxzzzz_0[j] + tg_0_xxxzzzzz_0[j];

                tg_z_xxyyyyy_0[j] = -ab_z * tg_0_xxyyyyy_0[j] + tg_0_xxyyyyyz_0[j];

                tg_z_xxyyyyz_0[j] = -ab_z * tg_0_xxyyyyz_0[j] + tg_0_xxyyyyzz_0[j];

                tg_z_xxyyyzz_0[j] = -ab_z * tg_0_xxyyyzz_0[j] + tg_0_xxyyyzzz_0[j];

                tg_z_xxyyzzz_0[j] = -ab_z * tg_0_xxyyzzz_0[j] + tg_0_xxyyzzzz_0[j];

                tg_z_xxyzzzz_0[j] = -ab_z * tg_0_xxyzzzz_0[j] + tg_0_xxyzzzzz_0[j];

                tg_z_xxzzzzz_0[j] = -ab_z * tg_0_xxzzzzz_0[j] + tg_0_xxzzzzzz_0[j];

                tg_z_xyyyyyy_0[j] = -ab_z * tg_0_xyyyyyy_0[j] + tg_0_xyyyyyyz_0[j];

                tg_z_xyyyyyz_0[j] = -ab_z * tg_0_xyyyyyz_0[j] + tg_0_xyyyyyzz_0[j];

                tg_z_xyyyyzz_0[j] = -ab_z * tg_0_xyyyyzz_0[j] + tg_0_xyyyyzzz_0[j];

                tg_z_xyyyzzz_0[j] = -ab_z * tg_0_xyyyzzz_0[j] + tg_0_xyyyzzzz_0[j];

                tg_z_xyyzzzz_0[j] = -ab_z * tg_0_xyyzzzz_0[j] + tg_0_xyyzzzzz_0[j];

                tg_z_xyzzzzz_0[j] = -ab_z * tg_0_xyzzzzz_0[j] + tg_0_xyzzzzzz_0[j];

                tg_z_xzzzzzz_0[j] = -ab_z * tg_0_xzzzzzz_0[j] + tg_0_xzzzzzzz_0[j];

                tg_z_yyyyyyy_0[j] = -ab_z * tg_0_yyyyyyy_0[j] + tg_0_yyyyyyyz_0[j];

                tg_z_yyyyyyz_0[j] = -ab_z * tg_0_yyyyyyz_0[j] + tg_0_yyyyyyzz_0[j];

                tg_z_yyyyyzz_0[j] = -ab_z * tg_0_yyyyyzz_0[j] + tg_0_yyyyyzzz_0[j];

                tg_z_yyyyzzz_0[j] = -ab_z * tg_0_yyyyzzz_0[j] + tg_0_yyyyzzzz_0[j];

                tg_z_yyyzzzz_0[j] = -ab_z * tg_0_yyyzzzz_0[j] + tg_0_yyyzzzzz_0[j];

                tg_z_yyzzzzz_0[j] = -ab_z * tg_0_yyzzzzz_0[j] + tg_0_yyzzzzzz_0[j];

                tg_z_yzzzzzz_0[j] = -ab_z * tg_0_yzzzzzz_0[j] + tg_0_yzzzzzzz_0[j];

                tg_z_zzzzzzz_0[j] = -ab_z * tg_0_zzzzzzz_0[j] + tg_0_zzzzzzzz_0[j];
            }
        }
    }


} // eribrrfunc namespace

