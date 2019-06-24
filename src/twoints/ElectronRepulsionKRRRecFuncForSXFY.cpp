//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionKRRRecFuncForSXFY.hpp"

#include "AngularMomentum.hpp"

namespace erikrrfunc { // erikrrfunc namespace

    void
    compElectronRepulsionForSXFF(      CMemBlock2D<double>& ketBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& cdDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to distances R(CD) = C - D

        auto cd_x = cdDistances.data(0);

        auto cd_y = cdDistances.data(1);

        auto cd_z = cdDistances.data(2);

        // set up bra side loop over intergals

        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();

        for (int32_t iang = 0; iang <= bang; iang++)
        {
            // set up index of integral

            auto pidx_g_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 3, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 3, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 4, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

                auto tg_xx_xxx_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i); 

                auto tg_xx_xxy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 1); 

                auto tg_xx_xxz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 2); 

                auto tg_xx_xyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 3); 

                auto tg_xx_xyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 4); 

                auto tg_xx_xzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 5); 

                auto tg_xx_yyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 6); 

                auto tg_xx_yyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 7); 

                auto tg_xx_yzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 8); 

                auto tg_xx_zzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 9); 

                auto tg_xy_xxx_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 10); 

                auto tg_xy_xxy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 11); 

                auto tg_xy_xxz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 12); 

                auto tg_xy_xyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 13); 

                auto tg_xy_xyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 14); 

                auto tg_xy_xzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 15); 

                auto tg_xy_yyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 16); 

                auto tg_xy_yyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 17); 

                auto tg_xy_yzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 18); 

                auto tg_xy_zzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 19); 

                auto tg_xz_xxx_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 20); 

                auto tg_xz_xxy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 21); 

                auto tg_xz_xxz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 22); 

                auto tg_xz_xyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 23); 

                auto tg_xz_xyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 24); 

                auto tg_xz_xzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 25); 

                auto tg_xz_yyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 26); 

                auto tg_xz_yyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 27); 

                auto tg_xz_yzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 28); 

                auto tg_xz_zzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 29); 

                auto tg_yy_xxx_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 30); 

                auto tg_yy_xxy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 31); 

                auto tg_yy_xxz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 32); 

                auto tg_yy_xyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 33); 

                auto tg_yy_xyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 34); 

                auto tg_yy_xzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 35); 

                auto tg_yy_yyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 36); 

                auto tg_yy_yyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 37); 

                auto tg_yy_yzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 38); 

                auto tg_yy_zzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 39); 

                auto tg_yz_xxx_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 40); 

                auto tg_yz_xxy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 41); 

                auto tg_yz_xxz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 42); 

                auto tg_yz_xyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 43); 

                auto tg_yz_xyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 44); 

                auto tg_yz_xzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 45); 

                auto tg_yz_yyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 46); 

                auto tg_yz_yyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 47); 

                auto tg_yz_yzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 48); 

                auto tg_yz_zzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 49); 

                auto tg_zz_xxx_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 50); 

                auto tg_zz_xxy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 51); 

                auto tg_zz_xxz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 52); 

                auto tg_zz_xyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 53); 

                auto tg_zz_xyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 54); 

                auto tg_zz_xzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 55); 

                auto tg_zz_yyy_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 56); 

                auto tg_zz_yyz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 57); 

                auto tg_zz_yzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 58); 

                auto tg_zz_zzz_0 = ketBuffer.data(pidx_g_2_3_m0 + 60 * i + 59); 

                auto tg_xx_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i); 

                auto tg_xx_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 1); 

                auto tg_xx_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 2); 

                auto tg_xx_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 3); 

                auto tg_xx_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 4); 

                auto tg_xx_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 5); 

                auto tg_xx_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 6); 

                auto tg_xx_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 7); 

                auto tg_xx_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 8); 

                auto tg_xx_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 9); 

                auto tg_xy_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 15); 

                auto tg_xy_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 16); 

                auto tg_xy_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 17); 

                auto tg_xy_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 18); 

                auto tg_xy_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 19); 

                auto tg_xy_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 20); 

                auto tg_xy_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 21); 

                auto tg_xy_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 22); 

                auto tg_xy_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 23); 

                auto tg_xy_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 24); 

                auto tg_xz_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 30); 

                auto tg_xz_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 31); 

                auto tg_xz_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 32); 

                auto tg_xz_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 33); 

                auto tg_xz_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 34); 

                auto tg_xz_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 35); 

                auto tg_xz_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 36); 

                auto tg_xz_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 37); 

                auto tg_xz_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 38); 

                auto tg_xz_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 39); 

                auto tg_yy_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 45); 

                auto tg_yy_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 46); 

                auto tg_yy_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 47); 

                auto tg_yy_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 48); 

                auto tg_yy_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 49); 

                auto tg_yy_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 50); 

                auto tg_yy_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 51); 

                auto tg_yy_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 52); 

                auto tg_yy_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 53); 

                auto tg_yy_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 54); 

                auto tg_yy_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 55); 

                auto tg_yy_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 56); 

                auto tg_yy_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 57); 

                auto tg_yy_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 58); 

                auto tg_yz_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 60); 

                auto tg_yz_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 61); 

                auto tg_yz_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 62); 

                auto tg_yz_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 63); 

                auto tg_yz_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 64); 

                auto tg_yz_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 65); 

                auto tg_yz_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 66); 

                auto tg_yz_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 67); 

                auto tg_yz_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 68); 

                auto tg_yz_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 69); 

                auto tg_yz_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 70); 

                auto tg_yz_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 71); 

                auto tg_yz_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 72); 

                auto tg_yz_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 73); 

                auto tg_zz_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 75); 

                auto tg_zz_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 76); 

                auto tg_zz_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 77); 

                auto tg_zz_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 78); 

                auto tg_zz_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 79); 

                auto tg_zz_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 80); 

                auto tg_zz_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 81); 

                auto tg_zz_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 82); 

                auto tg_zz_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 83); 

                auto tg_zz_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 84); 

                auto tg_zz_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 85); 

                auto tg_zz_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 86); 

                auto tg_zz_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 87); 

                auto tg_zz_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 88); 

                auto tg_zz_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 89); 

                // set up pointers to integrals

                auto tg_xxx_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i); 

                auto tg_xxx_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 1); 

                auto tg_xxx_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 2); 

                auto tg_xxx_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 3); 

                auto tg_xxx_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 4); 

                auto tg_xxx_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 5); 

                auto tg_xxx_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 6); 

                auto tg_xxx_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 7); 

                auto tg_xxx_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 8); 

                auto tg_xxx_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 9); 

                auto tg_xxy_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 10); 

                auto tg_xxy_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 11); 

                auto tg_xxy_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 12); 

                auto tg_xxy_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 13); 

                auto tg_xxy_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 14); 

                auto tg_xxy_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 15); 

                auto tg_xxy_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 16); 

                auto tg_xxy_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 17); 

                auto tg_xxy_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 18); 

                auto tg_xxy_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 19); 

                auto tg_xxz_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 20); 

                auto tg_xxz_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 21); 

                auto tg_xxz_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 22); 

                auto tg_xxz_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 23); 

                auto tg_xxz_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 24); 

                auto tg_xxz_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 25); 

                auto tg_xxz_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 26); 

                auto tg_xxz_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 27); 

                auto tg_xxz_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 28); 

                auto tg_xxz_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 29); 

                auto tg_xyy_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 30); 

                auto tg_xyy_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 31); 

                auto tg_xyy_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 32); 

                auto tg_xyy_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 33); 

                auto tg_xyy_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 34); 

                auto tg_xyy_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 35); 

                auto tg_xyy_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 36); 

                auto tg_xyy_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 37); 

                auto tg_xyy_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 38); 

                auto tg_xyy_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 39); 

                auto tg_xyz_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 40); 

                auto tg_xyz_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 41); 

                auto tg_xyz_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 42); 

                auto tg_xyz_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 43); 

                auto tg_xyz_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 44); 

                auto tg_xyz_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 45); 

                auto tg_xyz_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 46); 

                auto tg_xyz_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 47); 

                auto tg_xyz_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 48); 

                auto tg_xyz_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 49); 

                auto tg_xzz_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 50); 

                auto tg_xzz_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 51); 

                auto tg_xzz_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 52); 

                auto tg_xzz_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 53); 

                auto tg_xzz_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 54); 

                auto tg_xzz_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 55); 

                auto tg_xzz_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 56); 

                auto tg_xzz_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 57); 

                auto tg_xzz_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 58); 

                auto tg_xzz_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 59); 

                auto tg_yyy_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 60); 

                auto tg_yyy_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 61); 

                auto tg_yyy_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 62); 

                auto tg_yyy_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 63); 

                auto tg_yyy_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 64); 

                auto tg_yyy_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 65); 

                auto tg_yyy_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 66); 

                auto tg_yyy_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 67); 

                auto tg_yyy_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 68); 

                auto tg_yyy_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 69); 

                auto tg_yyz_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 70); 

                auto tg_yyz_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 71); 

                auto tg_yyz_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 72); 

                auto tg_yyz_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 73); 

                auto tg_yyz_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 74); 

                auto tg_yyz_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 75); 

                auto tg_yyz_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 76); 

                auto tg_yyz_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 77); 

                auto tg_yyz_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 78); 

                auto tg_yyz_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 79); 

                auto tg_yzz_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 80); 

                auto tg_yzz_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 81); 

                auto tg_yzz_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 82); 

                auto tg_yzz_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 83); 

                auto tg_yzz_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 84); 

                auto tg_yzz_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 85); 

                auto tg_yzz_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 86); 

                auto tg_yzz_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 87); 

                auto tg_yzz_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 88); 

                auto tg_yzz_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 89); 

                auto tg_zzz_xxx_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 90); 

                auto tg_zzz_xxy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 91); 

                auto tg_zzz_xxz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 92); 

                auto tg_zzz_xyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 93); 

                auto tg_zzz_xyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 94); 

                auto tg_zzz_xzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 95); 

                auto tg_zzz_yyy_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 96); 

                auto tg_zzz_yyz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 97); 

                auto tg_zzz_yzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 98); 

                auto tg_zzz_zzz_0 = ketBuffer.data(pidx_g_3_3_m0 + 100 * i + 99); 

                #pragma omp simd aligned(cd_x, cd_y, cd_z, tg_xx_xxx_0, tg_xx_xxxx_0, tg_xx_xxxy_0, tg_xx_xxxz_0, \
                                         tg_xx_xxy_0, tg_xx_xxyy_0, tg_xx_xxyz_0, tg_xx_xxz_0, tg_xx_xxzz_0, tg_xx_xyy_0, \
                                         tg_xx_xyyy_0, tg_xx_xyyz_0, tg_xx_xyz_0, tg_xx_xyzz_0, tg_xx_xzz_0, tg_xx_xzzz_0, \
                                         tg_xx_yyy_0, tg_xx_yyz_0, tg_xx_yzz_0, tg_xx_zzz_0, tg_xxx_xxx_0, tg_xxx_xxy_0, \
                                         tg_xxx_xxz_0, tg_xxx_xyy_0, tg_xxx_xyz_0, tg_xxx_xzz_0, tg_xxx_yyy_0, tg_xxx_yyz_0, \
                                         tg_xxx_yzz_0, tg_xxx_zzz_0, tg_xxy_xxx_0, tg_xxy_xxy_0, tg_xxy_xxz_0, tg_xxy_xyy_0, \
                                         tg_xxy_xyz_0, tg_xxy_xzz_0, tg_xxy_yyy_0, tg_xxy_yyz_0, tg_xxy_yzz_0, tg_xxy_zzz_0, \
                                         tg_xxz_xxx_0, tg_xxz_xxy_0, tg_xxz_xxz_0, tg_xxz_xyy_0, tg_xxz_xyz_0, tg_xxz_xzz_0, \
                                         tg_xxz_yyy_0, tg_xxz_yyz_0, tg_xxz_yzz_0, tg_xxz_zzz_0, tg_xy_xxx_0, tg_xy_xxxx_0, \
                                         tg_xy_xxxy_0, tg_xy_xxxz_0, tg_xy_xxy_0, tg_xy_xxyy_0, tg_xy_xxyz_0, tg_xy_xxz_0, \
                                         tg_xy_xxzz_0, tg_xy_xyy_0, tg_xy_xyyy_0, tg_xy_xyyz_0, tg_xy_xyz_0, tg_xy_xyzz_0, \
                                         tg_xy_xzz_0, tg_xy_xzzz_0, tg_xy_yyy_0, tg_xy_yyz_0, tg_xy_yzz_0, tg_xy_zzz_0, \
                                         tg_xyy_xxx_0, tg_xyy_xxy_0, tg_xyy_xxz_0, tg_xyy_xyy_0, tg_xyy_xyz_0, tg_xyy_xzz_0, \
                                         tg_xyy_yyy_0, tg_xyy_yyz_0, tg_xyy_yzz_0, tg_xyy_zzz_0, tg_xyz_xxx_0, tg_xyz_xxy_0, \
                                         tg_xyz_xxz_0, tg_xyz_xyy_0, tg_xyz_xyz_0, tg_xyz_xzz_0, tg_xyz_yyy_0, tg_xyz_yyz_0, \
                                         tg_xyz_yzz_0, tg_xyz_zzz_0, tg_xz_xxx_0, tg_xz_xxxx_0, tg_xz_xxxy_0, tg_xz_xxxz_0, \
                                         tg_xz_xxy_0, tg_xz_xxyy_0, tg_xz_xxyz_0, tg_xz_xxz_0, tg_xz_xxzz_0, tg_xz_xyy_0, \
                                         tg_xz_xyyy_0, tg_xz_xyyz_0, tg_xz_xyz_0, tg_xz_xyzz_0, tg_xz_xzz_0, tg_xz_xzzz_0, \
                                         tg_xz_yyy_0, tg_xz_yyz_0, tg_xz_yzz_0, tg_xz_zzz_0, tg_xzz_xxx_0, tg_xzz_xxy_0, \
                                         tg_xzz_xxz_0, tg_xzz_xyy_0, tg_xzz_xyz_0, tg_xzz_xzz_0, tg_xzz_yyy_0, tg_xzz_yyz_0, \
                                         tg_xzz_yzz_0, tg_xzz_zzz_0, tg_yy_xxx_0, tg_yy_xxxx_0, tg_yy_xxxy_0, tg_yy_xxxz_0, \
                                         tg_yy_xxy_0, tg_yy_xxyy_0, tg_yy_xxyz_0, tg_yy_xxz_0, tg_yy_xxzz_0, tg_yy_xyy_0, \
                                         tg_yy_xyyy_0, tg_yy_xyyz_0, tg_yy_xyz_0, tg_yy_xyzz_0, tg_yy_xzz_0, tg_yy_xzzz_0, \
                                         tg_yy_yyy_0, tg_yy_yyyy_0, tg_yy_yyyz_0, tg_yy_yyz_0, tg_yy_yyzz_0, tg_yy_yzz_0, \
                                         tg_yy_yzzz_0, tg_yy_zzz_0, tg_yyy_xxx_0, tg_yyy_xxy_0, tg_yyy_xxz_0, tg_yyy_xyy_0, \
                                         tg_yyy_xyz_0, tg_yyy_xzz_0, tg_yyy_yyy_0, tg_yyy_yyz_0, tg_yyy_yzz_0, tg_yyy_zzz_0, \
                                         tg_yyz_xxx_0, tg_yyz_xxy_0, tg_yyz_xxz_0, tg_yyz_xyy_0, tg_yyz_xyz_0, tg_yyz_xzz_0, \
                                         tg_yyz_yyy_0, tg_yyz_yyz_0, tg_yyz_yzz_0, tg_yyz_zzz_0, tg_yz_xxx_0, tg_yz_xxxx_0, \
                                         tg_yz_xxxy_0, tg_yz_xxxz_0, tg_yz_xxy_0, tg_yz_xxyy_0, tg_yz_xxyz_0, tg_yz_xxz_0, \
                                         tg_yz_xxzz_0, tg_yz_xyy_0, tg_yz_xyyy_0, tg_yz_xyyz_0, tg_yz_xyz_0, tg_yz_xyzz_0, \
                                         tg_yz_xzz_0, tg_yz_xzzz_0, tg_yz_yyy_0, tg_yz_yyyy_0, tg_yz_yyyz_0, tg_yz_yyz_0, \
                                         tg_yz_yyzz_0, tg_yz_yzz_0, tg_yz_yzzz_0, tg_yz_zzz_0, tg_yzz_xxx_0, tg_yzz_xxy_0, \
                                         tg_yzz_xxz_0, tg_yzz_xyy_0, tg_yzz_xyz_0, tg_yzz_xzz_0, tg_yzz_yyy_0, tg_yzz_yyz_0, \
                                         tg_yzz_yzz_0, tg_yzz_zzz_0, tg_zz_xxx_0, tg_zz_xxxx_0, tg_zz_xxxy_0, tg_zz_xxxz_0, \
                                         tg_zz_xxy_0, tg_zz_xxyy_0, tg_zz_xxyz_0, tg_zz_xxz_0, tg_zz_xxzz_0, tg_zz_xyy_0, \
                                         tg_zz_xyyy_0, tg_zz_xyyz_0, tg_zz_xyz_0, tg_zz_xyzz_0, tg_zz_xzz_0, tg_zz_xzzz_0, \
                                         tg_zz_yyy_0, tg_zz_yyyy_0, tg_zz_yyyz_0, tg_zz_yyz_0, tg_zz_yyzz_0, tg_zz_yzz_0, \
                                         tg_zz_yzzz_0, tg_zz_zzz_0, tg_zz_zzzz_0, tg_zzz_xxx_0, tg_zzz_xxy_0, tg_zzz_xxz_0, \
                                         tg_zzz_xyy_0, tg_zzz_xyz_0, tg_zzz_xzz_0, tg_zzz_yyy_0, tg_zzz_yyz_0, tg_zzz_yzz_0, \
                                         tg_zzz_zzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nKetContrPairs; j++)
                {
                    tg_xxx_xxx_0[j] = -cd_x[j] * tg_xx_xxx_0[j] + tg_xx_xxxx_0[j];

                    tg_xxx_xxy_0[j] = -cd_x[j] * tg_xx_xxy_0[j] + tg_xx_xxxy_0[j];

                    tg_xxx_xxz_0[j] = -cd_x[j] * tg_xx_xxz_0[j] + tg_xx_xxxz_0[j];

                    tg_xxx_xyy_0[j] = -cd_x[j] * tg_xx_xyy_0[j] + tg_xx_xxyy_0[j];

                    tg_xxx_xyz_0[j] = -cd_x[j] * tg_xx_xyz_0[j] + tg_xx_xxyz_0[j];

                    tg_xxx_xzz_0[j] = -cd_x[j] * tg_xx_xzz_0[j] + tg_xx_xxzz_0[j];

                    tg_xxx_yyy_0[j] = -cd_x[j] * tg_xx_yyy_0[j] + tg_xx_xyyy_0[j];

                    tg_xxx_yyz_0[j] = -cd_x[j] * tg_xx_yyz_0[j] + tg_xx_xyyz_0[j];

                    tg_xxx_yzz_0[j] = -cd_x[j] * tg_xx_yzz_0[j] + tg_xx_xyzz_0[j];

                    tg_xxx_zzz_0[j] = -cd_x[j] * tg_xx_zzz_0[j] + tg_xx_xzzz_0[j];

                    tg_xxy_xxx_0[j] = -cd_x[j] * tg_xy_xxx_0[j] + tg_xy_xxxx_0[j];

                    tg_xxy_xxy_0[j] = -cd_x[j] * tg_xy_xxy_0[j] + tg_xy_xxxy_0[j];

                    tg_xxy_xxz_0[j] = -cd_x[j] * tg_xy_xxz_0[j] + tg_xy_xxxz_0[j];

                    tg_xxy_xyy_0[j] = -cd_x[j] * tg_xy_xyy_0[j] + tg_xy_xxyy_0[j];

                    tg_xxy_xyz_0[j] = -cd_x[j] * tg_xy_xyz_0[j] + tg_xy_xxyz_0[j];

                    tg_xxy_xzz_0[j] = -cd_x[j] * tg_xy_xzz_0[j] + tg_xy_xxzz_0[j];

                    tg_xxy_yyy_0[j] = -cd_x[j] * tg_xy_yyy_0[j] + tg_xy_xyyy_0[j];

                    tg_xxy_yyz_0[j] = -cd_x[j] * tg_xy_yyz_0[j] + tg_xy_xyyz_0[j];

                    tg_xxy_yzz_0[j] = -cd_x[j] * tg_xy_yzz_0[j] + tg_xy_xyzz_0[j];

                    tg_xxy_zzz_0[j] = -cd_x[j] * tg_xy_zzz_0[j] + tg_xy_xzzz_0[j];

                    tg_xxz_xxx_0[j] = -cd_x[j] * tg_xz_xxx_0[j] + tg_xz_xxxx_0[j];

                    tg_xxz_xxy_0[j] = -cd_x[j] * tg_xz_xxy_0[j] + tg_xz_xxxy_0[j];

                    tg_xxz_xxz_0[j] = -cd_x[j] * tg_xz_xxz_0[j] + tg_xz_xxxz_0[j];

                    tg_xxz_xyy_0[j] = -cd_x[j] * tg_xz_xyy_0[j] + tg_xz_xxyy_0[j];

                    tg_xxz_xyz_0[j] = -cd_x[j] * tg_xz_xyz_0[j] + tg_xz_xxyz_0[j];

                    tg_xxz_xzz_0[j] = -cd_x[j] * tg_xz_xzz_0[j] + tg_xz_xxzz_0[j];

                    tg_xxz_yyy_0[j] = -cd_x[j] * tg_xz_yyy_0[j] + tg_xz_xyyy_0[j];

                    tg_xxz_yyz_0[j] = -cd_x[j] * tg_xz_yyz_0[j] + tg_xz_xyyz_0[j];

                    tg_xxz_yzz_0[j] = -cd_x[j] * tg_xz_yzz_0[j] + tg_xz_xyzz_0[j];

                    tg_xxz_zzz_0[j] = -cd_x[j] * tg_xz_zzz_0[j] + tg_xz_xzzz_0[j];

                    tg_xyy_xxx_0[j] = -cd_x[j] * tg_yy_xxx_0[j] + tg_yy_xxxx_0[j];

                    tg_xyy_xxy_0[j] = -cd_x[j] * tg_yy_xxy_0[j] + tg_yy_xxxy_0[j];

                    tg_xyy_xxz_0[j] = -cd_x[j] * tg_yy_xxz_0[j] + tg_yy_xxxz_0[j];

                    tg_xyy_xyy_0[j] = -cd_x[j] * tg_yy_xyy_0[j] + tg_yy_xxyy_0[j];

                    tg_xyy_xyz_0[j] = -cd_x[j] * tg_yy_xyz_0[j] + tg_yy_xxyz_0[j];

                    tg_xyy_xzz_0[j] = -cd_x[j] * tg_yy_xzz_0[j] + tg_yy_xxzz_0[j];

                    tg_xyy_yyy_0[j] = -cd_x[j] * tg_yy_yyy_0[j] + tg_yy_xyyy_0[j];

                    tg_xyy_yyz_0[j] = -cd_x[j] * tg_yy_yyz_0[j] + tg_yy_xyyz_0[j];

                    tg_xyy_yzz_0[j] = -cd_x[j] * tg_yy_yzz_0[j] + tg_yy_xyzz_0[j];

                    tg_xyy_zzz_0[j] = -cd_x[j] * tg_yy_zzz_0[j] + tg_yy_xzzz_0[j];

                    tg_xyz_xxx_0[j] = -cd_x[j] * tg_yz_xxx_0[j] + tg_yz_xxxx_0[j];

                    tg_xyz_xxy_0[j] = -cd_x[j] * tg_yz_xxy_0[j] + tg_yz_xxxy_0[j];

                    tg_xyz_xxz_0[j] = -cd_x[j] * tg_yz_xxz_0[j] + tg_yz_xxxz_0[j];

                    tg_xyz_xyy_0[j] = -cd_x[j] * tg_yz_xyy_0[j] + tg_yz_xxyy_0[j];

                    tg_xyz_xyz_0[j] = -cd_x[j] * tg_yz_xyz_0[j] + tg_yz_xxyz_0[j];

                    tg_xyz_xzz_0[j] = -cd_x[j] * tg_yz_xzz_0[j] + tg_yz_xxzz_0[j];

                    tg_xyz_yyy_0[j] = -cd_x[j] * tg_yz_yyy_0[j] + tg_yz_xyyy_0[j];

                    tg_xyz_yyz_0[j] = -cd_x[j] * tg_yz_yyz_0[j] + tg_yz_xyyz_0[j];

                    tg_xyz_yzz_0[j] = -cd_x[j] * tg_yz_yzz_0[j] + tg_yz_xyzz_0[j];

                    tg_xyz_zzz_0[j] = -cd_x[j] * tg_yz_zzz_0[j] + tg_yz_xzzz_0[j];

                    tg_xzz_xxx_0[j] = -cd_x[j] * tg_zz_xxx_0[j] + tg_zz_xxxx_0[j];

                    tg_xzz_xxy_0[j] = -cd_x[j] * tg_zz_xxy_0[j] + tg_zz_xxxy_0[j];

                    tg_xzz_xxz_0[j] = -cd_x[j] * tg_zz_xxz_0[j] + tg_zz_xxxz_0[j];

                    tg_xzz_xyy_0[j] = -cd_x[j] * tg_zz_xyy_0[j] + tg_zz_xxyy_0[j];

                    tg_xzz_xyz_0[j] = -cd_x[j] * tg_zz_xyz_0[j] + tg_zz_xxyz_0[j];

                    tg_xzz_xzz_0[j] = -cd_x[j] * tg_zz_xzz_0[j] + tg_zz_xxzz_0[j];

                    tg_xzz_yyy_0[j] = -cd_x[j] * tg_zz_yyy_0[j] + tg_zz_xyyy_0[j];

                    tg_xzz_yyz_0[j] = -cd_x[j] * tg_zz_yyz_0[j] + tg_zz_xyyz_0[j];

                    tg_xzz_yzz_0[j] = -cd_x[j] * tg_zz_yzz_0[j] + tg_zz_xyzz_0[j];

                    tg_xzz_zzz_0[j] = -cd_x[j] * tg_zz_zzz_0[j] + tg_zz_xzzz_0[j];

                    tg_yyy_xxx_0[j] = -cd_y[j] * tg_yy_xxx_0[j] + tg_yy_xxxy_0[j];

                    tg_yyy_xxy_0[j] = -cd_y[j] * tg_yy_xxy_0[j] + tg_yy_xxyy_0[j];

                    tg_yyy_xxz_0[j] = -cd_y[j] * tg_yy_xxz_0[j] + tg_yy_xxyz_0[j];

                    tg_yyy_xyy_0[j] = -cd_y[j] * tg_yy_xyy_0[j] + tg_yy_xyyy_0[j];

                    tg_yyy_xyz_0[j] = -cd_y[j] * tg_yy_xyz_0[j] + tg_yy_xyyz_0[j];

                    tg_yyy_xzz_0[j] = -cd_y[j] * tg_yy_xzz_0[j] + tg_yy_xyzz_0[j];

                    tg_yyy_yyy_0[j] = -cd_y[j] * tg_yy_yyy_0[j] + tg_yy_yyyy_0[j];

                    tg_yyy_yyz_0[j] = -cd_y[j] * tg_yy_yyz_0[j] + tg_yy_yyyz_0[j];

                    tg_yyy_yzz_0[j] = -cd_y[j] * tg_yy_yzz_0[j] + tg_yy_yyzz_0[j];

                    tg_yyy_zzz_0[j] = -cd_y[j] * tg_yy_zzz_0[j] + tg_yy_yzzz_0[j];

                    tg_yyz_xxx_0[j] = -cd_y[j] * tg_yz_xxx_0[j] + tg_yz_xxxy_0[j];

                    tg_yyz_xxy_0[j] = -cd_y[j] * tg_yz_xxy_0[j] + tg_yz_xxyy_0[j];

                    tg_yyz_xxz_0[j] = -cd_y[j] * tg_yz_xxz_0[j] + tg_yz_xxyz_0[j];

                    tg_yyz_xyy_0[j] = -cd_y[j] * tg_yz_xyy_0[j] + tg_yz_xyyy_0[j];

                    tg_yyz_xyz_0[j] = -cd_y[j] * tg_yz_xyz_0[j] + tg_yz_xyyz_0[j];

                    tg_yyz_xzz_0[j] = -cd_y[j] * tg_yz_xzz_0[j] + tg_yz_xyzz_0[j];

                    tg_yyz_yyy_0[j] = -cd_y[j] * tg_yz_yyy_0[j] + tg_yz_yyyy_0[j];

                    tg_yyz_yyz_0[j] = -cd_y[j] * tg_yz_yyz_0[j] + tg_yz_yyyz_0[j];

                    tg_yyz_yzz_0[j] = -cd_y[j] * tg_yz_yzz_0[j] + tg_yz_yyzz_0[j];

                    tg_yyz_zzz_0[j] = -cd_y[j] * tg_yz_zzz_0[j] + tg_yz_yzzz_0[j];

                    tg_yzz_xxx_0[j] = -cd_y[j] * tg_zz_xxx_0[j] + tg_zz_xxxy_0[j];

                    tg_yzz_xxy_0[j] = -cd_y[j] * tg_zz_xxy_0[j] + tg_zz_xxyy_0[j];

                    tg_yzz_xxz_0[j] = -cd_y[j] * tg_zz_xxz_0[j] + tg_zz_xxyz_0[j];

                    tg_yzz_xyy_0[j] = -cd_y[j] * tg_zz_xyy_0[j] + tg_zz_xyyy_0[j];

                    tg_yzz_xyz_0[j] = -cd_y[j] * tg_zz_xyz_0[j] + tg_zz_xyyz_0[j];

                    tg_yzz_xzz_0[j] = -cd_y[j] * tg_zz_xzz_0[j] + tg_zz_xyzz_0[j];

                    tg_yzz_yyy_0[j] = -cd_y[j] * tg_zz_yyy_0[j] + tg_zz_yyyy_0[j];

                    tg_yzz_yyz_0[j] = -cd_y[j] * tg_zz_yyz_0[j] + tg_zz_yyyz_0[j];

                    tg_yzz_yzz_0[j] = -cd_y[j] * tg_zz_yzz_0[j] + tg_zz_yyzz_0[j];

                    tg_yzz_zzz_0[j] = -cd_y[j] * tg_zz_zzz_0[j] + tg_zz_yzzz_0[j];

                    tg_zzz_xxx_0[j] = -cd_z[j] * tg_zz_xxx_0[j] + tg_zz_xxxz_0[j];

                    tg_zzz_xxy_0[j] = -cd_z[j] * tg_zz_xxy_0[j] + tg_zz_xxyz_0[j];

                    tg_zzz_xxz_0[j] = -cd_z[j] * tg_zz_xxz_0[j] + tg_zz_xxzz_0[j];

                    tg_zzz_xyy_0[j] = -cd_z[j] * tg_zz_xyy_0[j] + tg_zz_xyyz_0[j];

                    tg_zzz_xyz_0[j] = -cd_z[j] * tg_zz_xyz_0[j] + tg_zz_xyzz_0[j];

                    tg_zzz_xzz_0[j] = -cd_z[j] * tg_zz_xzz_0[j] + tg_zz_xzzz_0[j];

                    tg_zzz_yyy_0[j] = -cd_z[j] * tg_zz_yyy_0[j] + tg_zz_yyyz_0[j];

                    tg_zzz_yyz_0[j] = -cd_z[j] * tg_zz_yyz_0[j] + tg_zz_yyzz_0[j];

                    tg_zzz_yzz_0[j] = -cd_z[j] * tg_zz_yzz_0[j] + tg_zz_yzzz_0[j];

                    tg_zzz_zzz_0[j] = -cd_z[j] * tg_zz_zzz_0[j] + tg_zz_zzzz_0[j];
                }
            }
        }
    }

    void
    compElectronRepulsionForSXFG(      CMemBlock2D<double>& ketBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& cdDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        erikrrfunc::compElectronRepulsionForSXFG_0_75(ketBuffer,
                                                      recursionMap,
                                                      cdDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        erikrrfunc::compElectronRepulsionForSXFG_75_150(ketBuffer,
                                                        recursionMap,
                                                        cdDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForSXFG_0_75(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,75)

        // set up pointers to distances R(CD) = C - D

        auto cd_x = cdDistances.data(0);

        // set up bra side loop over intergals

        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();

        for (int32_t iang = 0; iang <= bang; iang++)
        {
            // set up index of integral

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 4, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 4, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 5, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

                auto tg_xx_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i); 

                auto tg_xx_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 1); 

                auto tg_xx_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 2); 

                auto tg_xx_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 3); 

                auto tg_xx_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 4); 

                auto tg_xx_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 5); 

                auto tg_xx_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 6); 

                auto tg_xx_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 7); 

                auto tg_xx_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 8); 

                auto tg_xx_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 9); 

                auto tg_xx_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 10); 

                auto tg_xx_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 11); 

                auto tg_xx_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 12); 

                auto tg_xx_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 13); 

                auto tg_xx_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 14); 

                auto tg_xy_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 15); 

                auto tg_xy_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 16); 

                auto tg_xy_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 17); 

                auto tg_xy_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 18); 

                auto tg_xy_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 19); 

                auto tg_xy_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 20); 

                auto tg_xy_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 21); 

                auto tg_xy_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 22); 

                auto tg_xy_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 23); 

                auto tg_xy_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 24); 

                auto tg_xy_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 25); 

                auto tg_xy_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 26); 

                auto tg_xy_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 27); 

                auto tg_xy_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 28); 

                auto tg_xy_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 29); 

                auto tg_xz_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 30); 

                auto tg_xz_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 31); 

                auto tg_xz_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 32); 

                auto tg_xz_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 33); 

                auto tg_xz_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 34); 

                auto tg_xz_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 35); 

                auto tg_xz_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 36); 

                auto tg_xz_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 37); 

                auto tg_xz_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 38); 

                auto tg_xz_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 39); 

                auto tg_xz_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 40); 

                auto tg_xz_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 41); 

                auto tg_xz_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 42); 

                auto tg_xz_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 43); 

                auto tg_xz_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 44); 

                auto tg_yy_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 45); 

                auto tg_yy_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 46); 

                auto tg_yy_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 47); 

                auto tg_yy_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 48); 

                auto tg_yy_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 49); 

                auto tg_yy_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 50); 

                auto tg_yy_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 51); 

                auto tg_yy_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 52); 

                auto tg_yy_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 53); 

                auto tg_yy_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 54); 

                auto tg_yy_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 55); 

                auto tg_yy_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 56); 

                auto tg_yy_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 57); 

                auto tg_yy_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 58); 

                auto tg_yy_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 59); 

                auto tg_yz_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 60); 

                auto tg_yz_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 61); 

                auto tg_yz_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 62); 

                auto tg_yz_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 63); 

                auto tg_yz_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 64); 

                auto tg_yz_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 65); 

                auto tg_yz_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 66); 

                auto tg_yz_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 67); 

                auto tg_yz_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 68); 

                auto tg_yz_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 69); 

                auto tg_yz_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 70); 

                auto tg_yz_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 71); 

                auto tg_yz_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 72); 

                auto tg_yz_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 73); 

                auto tg_yz_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 74); 

                auto tg_xx_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i); 

                auto tg_xx_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 1); 

                auto tg_xx_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 2); 

                auto tg_xx_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 3); 

                auto tg_xx_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 4); 

                auto tg_xx_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 5); 

                auto tg_xx_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 6); 

                auto tg_xx_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 7); 

                auto tg_xx_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 8); 

                auto tg_xx_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 9); 

                auto tg_xx_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 10); 

                auto tg_xx_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 11); 

                auto tg_xx_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 12); 

                auto tg_xx_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 13); 

                auto tg_xx_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 14); 

                auto tg_xy_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 21); 

                auto tg_xy_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 22); 

                auto tg_xy_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 23); 

                auto tg_xy_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 24); 

                auto tg_xy_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 25); 

                auto tg_xy_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 26); 

                auto tg_xy_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 27); 

                auto tg_xy_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 28); 

                auto tg_xy_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 29); 

                auto tg_xy_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 30); 

                auto tg_xy_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 31); 

                auto tg_xy_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 32); 

                auto tg_xy_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 33); 

                auto tg_xy_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 34); 

                auto tg_xy_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 35); 

                auto tg_xz_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 42); 

                auto tg_xz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 43); 

                auto tg_xz_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 44); 

                auto tg_xz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 45); 

                auto tg_xz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 46); 

                auto tg_xz_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 47); 

                auto tg_xz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 48); 

                auto tg_xz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 49); 

                auto tg_xz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 50); 

                auto tg_xz_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 51); 

                auto tg_xz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 52); 

                auto tg_xz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 53); 

                auto tg_xz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 54); 

                auto tg_xz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 55); 

                auto tg_xz_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 56); 

                auto tg_yy_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 63); 

                auto tg_yy_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 64); 

                auto tg_yy_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 65); 

                auto tg_yy_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 66); 

                auto tg_yy_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 67); 

                auto tg_yy_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 68); 

                auto tg_yy_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 69); 

                auto tg_yy_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 70); 

                auto tg_yy_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 71); 

                auto tg_yy_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 72); 

                auto tg_yy_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 73); 

                auto tg_yy_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 74); 

                auto tg_yy_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 75); 

                auto tg_yy_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 76); 

                auto tg_yy_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 77); 

                auto tg_yz_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 84); 

                auto tg_yz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 85); 

                auto tg_yz_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 86); 

                auto tg_yz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 87); 

                auto tg_yz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 88); 

                auto tg_yz_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 89); 

                auto tg_yz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 90); 

                auto tg_yz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 91); 

                auto tg_yz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 92); 

                auto tg_yz_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 93); 

                auto tg_yz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 94); 

                auto tg_yz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 95); 

                auto tg_yz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 96); 

                auto tg_yz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 97); 

                auto tg_yz_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 98); 

                // set up pointers to integrals

                auto tg_xxx_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i); 

                auto tg_xxx_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 1); 

                auto tg_xxx_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 2); 

                auto tg_xxx_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 3); 

                auto tg_xxx_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 4); 

                auto tg_xxx_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 5); 

                auto tg_xxx_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 6); 

                auto tg_xxx_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 7); 

                auto tg_xxx_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 8); 

                auto tg_xxx_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 9); 

                auto tg_xxx_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 10); 

                auto tg_xxx_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 11); 

                auto tg_xxx_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 12); 

                auto tg_xxx_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 13); 

                auto tg_xxx_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 14); 

                auto tg_xxy_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 15); 

                auto tg_xxy_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 16); 

                auto tg_xxy_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 17); 

                auto tg_xxy_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 18); 

                auto tg_xxy_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 19); 

                auto tg_xxy_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 20); 

                auto tg_xxy_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 21); 

                auto tg_xxy_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 22); 

                auto tg_xxy_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 23); 

                auto tg_xxy_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 24); 

                auto tg_xxy_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 25); 

                auto tg_xxy_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 26); 

                auto tg_xxy_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 27); 

                auto tg_xxy_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 28); 

                auto tg_xxy_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 29); 

                auto tg_xxz_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 30); 

                auto tg_xxz_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 31); 

                auto tg_xxz_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 32); 

                auto tg_xxz_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 33); 

                auto tg_xxz_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 34); 

                auto tg_xxz_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 35); 

                auto tg_xxz_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 36); 

                auto tg_xxz_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 37); 

                auto tg_xxz_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 38); 

                auto tg_xxz_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 39); 

                auto tg_xxz_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 40); 

                auto tg_xxz_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 41); 

                auto tg_xxz_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 42); 

                auto tg_xxz_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 43); 

                auto tg_xxz_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 44); 

                auto tg_xyy_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 45); 

                auto tg_xyy_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 46); 

                auto tg_xyy_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 47); 

                auto tg_xyy_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 48); 

                auto tg_xyy_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 49); 

                auto tg_xyy_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 50); 

                auto tg_xyy_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 51); 

                auto tg_xyy_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 52); 

                auto tg_xyy_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 53); 

                auto tg_xyy_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 54); 

                auto tg_xyy_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 55); 

                auto tg_xyy_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 56); 

                auto tg_xyy_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 57); 

                auto tg_xyy_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 58); 

                auto tg_xyy_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 59); 

                auto tg_xyz_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 60); 

                auto tg_xyz_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 61); 

                auto tg_xyz_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 62); 

                auto tg_xyz_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 63); 

                auto tg_xyz_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 64); 

                auto tg_xyz_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 65); 

                auto tg_xyz_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 66); 

                auto tg_xyz_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 67); 

                auto tg_xyz_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 68); 

                auto tg_xyz_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 69); 

                auto tg_xyz_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 70); 

                auto tg_xyz_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 71); 

                auto tg_xyz_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 72); 

                auto tg_xyz_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 73); 

                auto tg_xyz_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 74); 

                // Batch of Integrals (0,75)

                #pragma omp simd aligned(cd_x, tg_xx_xxxx_0, tg_xx_xxxxx_0, tg_xx_xxxxy_0, tg_xx_xxxxz_0, \
                                         tg_xx_xxxy_0, tg_xx_xxxyy_0, tg_xx_xxxyz_0, tg_xx_xxxz_0, tg_xx_xxxzz_0, \
                                         tg_xx_xxyy_0, tg_xx_xxyyy_0, tg_xx_xxyyz_0, tg_xx_xxyz_0, tg_xx_xxyzz_0, \
                                         tg_xx_xxzz_0, tg_xx_xxzzz_0, tg_xx_xyyy_0, tg_xx_xyyyy_0, tg_xx_xyyyz_0, \
                                         tg_xx_xyyz_0, tg_xx_xyyzz_0, tg_xx_xyzz_0, tg_xx_xyzzz_0, tg_xx_xzzz_0, \
                                         tg_xx_xzzzz_0, tg_xx_yyyy_0, tg_xx_yyyz_0, tg_xx_yyzz_0, tg_xx_yzzz_0, tg_xx_zzzz_0, \
                                         tg_xxx_xxxx_0, tg_xxx_xxxy_0, tg_xxx_xxxz_0, tg_xxx_xxyy_0, tg_xxx_xxyz_0, \
                                         tg_xxx_xxzz_0, tg_xxx_xyyy_0, tg_xxx_xyyz_0, tg_xxx_xyzz_0, tg_xxx_xzzz_0, \
                                         tg_xxx_yyyy_0, tg_xxx_yyyz_0, tg_xxx_yyzz_0, tg_xxx_yzzz_0, tg_xxx_zzzz_0, \
                                         tg_xxy_xxxx_0, tg_xxy_xxxy_0, tg_xxy_xxxz_0, tg_xxy_xxyy_0, tg_xxy_xxyz_0, \
                                         tg_xxy_xxzz_0, tg_xxy_xyyy_0, tg_xxy_xyyz_0, tg_xxy_xyzz_0, tg_xxy_xzzz_0, \
                                         tg_xxy_yyyy_0, tg_xxy_yyyz_0, tg_xxy_yyzz_0, tg_xxy_yzzz_0, tg_xxy_zzzz_0, \
                                         tg_xxz_xxxx_0, tg_xxz_xxxy_0, tg_xxz_xxxz_0, tg_xxz_xxyy_0, tg_xxz_xxyz_0, \
                                         tg_xxz_xxzz_0, tg_xxz_xyyy_0, tg_xxz_xyyz_0, tg_xxz_xyzz_0, tg_xxz_xzzz_0, \
                                         tg_xxz_yyyy_0, tg_xxz_yyyz_0, tg_xxz_yyzz_0, tg_xxz_yzzz_0, tg_xxz_zzzz_0, \
                                         tg_xy_xxxx_0, tg_xy_xxxxx_0, tg_xy_xxxxy_0, tg_xy_xxxxz_0, tg_xy_xxxy_0, \
                                         tg_xy_xxxyy_0, tg_xy_xxxyz_0, tg_xy_xxxz_0, tg_xy_xxxzz_0, tg_xy_xxyy_0, \
                                         tg_xy_xxyyy_0, tg_xy_xxyyz_0, tg_xy_xxyz_0, tg_xy_xxyzz_0, tg_xy_xxzz_0, \
                                         tg_xy_xxzzz_0, tg_xy_xyyy_0, tg_xy_xyyyy_0, tg_xy_xyyyz_0, tg_xy_xyyz_0, \
                                         tg_xy_xyyzz_0, tg_xy_xyzz_0, tg_xy_xyzzz_0, tg_xy_xzzz_0, tg_xy_xzzzz_0, \
                                         tg_xy_yyyy_0, tg_xy_yyyz_0, tg_xy_yyzz_0, tg_xy_yzzz_0, tg_xy_zzzz_0, \
                                         tg_xyy_xxxx_0, tg_xyy_xxxy_0, tg_xyy_xxxz_0, tg_xyy_xxyy_0, tg_xyy_xxyz_0, \
                                         tg_xyy_xxzz_0, tg_xyy_xyyy_0, tg_xyy_xyyz_0, tg_xyy_xyzz_0, tg_xyy_xzzz_0, \
                                         tg_xyy_yyyy_0, tg_xyy_yyyz_0, tg_xyy_yyzz_0, tg_xyy_yzzz_0, tg_xyy_zzzz_0, \
                                         tg_xyz_xxxx_0, tg_xyz_xxxy_0, tg_xyz_xxxz_0, tg_xyz_xxyy_0, tg_xyz_xxyz_0, \
                                         tg_xyz_xxzz_0, tg_xyz_xyyy_0, tg_xyz_xyyz_0, tg_xyz_xyzz_0, tg_xyz_xzzz_0, \
                                         tg_xyz_yyyy_0, tg_xyz_yyyz_0, tg_xyz_yyzz_0, tg_xyz_yzzz_0, tg_xyz_zzzz_0, \
                                         tg_xz_xxxx_0, tg_xz_xxxxx_0, tg_xz_xxxxy_0, tg_xz_xxxxz_0, tg_xz_xxxy_0, \
                                         tg_xz_xxxyy_0, tg_xz_xxxyz_0, tg_xz_xxxz_0, tg_xz_xxxzz_0, tg_xz_xxyy_0, \
                                         tg_xz_xxyyy_0, tg_xz_xxyyz_0, tg_xz_xxyz_0, tg_xz_xxyzz_0, tg_xz_xxzz_0, \
                                         tg_xz_xxzzz_0, tg_xz_xyyy_0, tg_xz_xyyyy_0, tg_xz_xyyyz_0, tg_xz_xyyz_0, \
                                         tg_xz_xyyzz_0, tg_xz_xyzz_0, tg_xz_xyzzz_0, tg_xz_xzzz_0, tg_xz_xzzzz_0, \
                                         tg_xz_yyyy_0, tg_xz_yyyz_0, tg_xz_yyzz_0, tg_xz_yzzz_0, tg_xz_zzzz_0, tg_yy_xxxx_0, \
                                         tg_yy_xxxxx_0, tg_yy_xxxxy_0, tg_yy_xxxxz_0, tg_yy_xxxy_0, tg_yy_xxxyy_0, \
                                         tg_yy_xxxyz_0, tg_yy_xxxz_0, tg_yy_xxxzz_0, tg_yy_xxyy_0, tg_yy_xxyyy_0, \
                                         tg_yy_xxyyz_0, tg_yy_xxyz_0, tg_yy_xxyzz_0, tg_yy_xxzz_0, tg_yy_xxzzz_0, \
                                         tg_yy_xyyy_0, tg_yy_xyyyy_0, tg_yy_xyyyz_0, tg_yy_xyyz_0, tg_yy_xyyzz_0, \
                                         tg_yy_xyzz_0, tg_yy_xyzzz_0, tg_yy_xzzz_0, tg_yy_xzzzz_0, tg_yy_yyyy_0, \
                                         tg_yy_yyyz_0, tg_yy_yyzz_0, tg_yy_yzzz_0, tg_yy_zzzz_0, tg_yz_xxxx_0, \
                                         tg_yz_xxxxx_0, tg_yz_xxxxy_0, tg_yz_xxxxz_0, tg_yz_xxxy_0, tg_yz_xxxyy_0, \
                                         tg_yz_xxxyz_0, tg_yz_xxxz_0, tg_yz_xxxzz_0, tg_yz_xxyy_0, tg_yz_xxyyy_0, \
                                         tg_yz_xxyyz_0, tg_yz_xxyz_0, tg_yz_xxyzz_0, tg_yz_xxzz_0, tg_yz_xxzzz_0, \
                                         tg_yz_xyyy_0, tg_yz_xyyyy_0, tg_yz_xyyyz_0, tg_yz_xyyz_0, tg_yz_xyyzz_0, \
                                         tg_yz_xyzz_0, tg_yz_xyzzz_0, tg_yz_xzzz_0, tg_yz_xzzzz_0, tg_yz_yyyy_0, \
                                         tg_yz_yyyz_0, tg_yz_yyzz_0, tg_yz_yzzz_0, tg_yz_zzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nKetContrPairs; j++)
                {
                    tg_xxx_xxxx_0[j] = -cd_x[j] * tg_xx_xxxx_0[j] + tg_xx_xxxxx_0[j];

                    tg_xxx_xxxy_0[j] = -cd_x[j] * tg_xx_xxxy_0[j] + tg_xx_xxxxy_0[j];

                    tg_xxx_xxxz_0[j] = -cd_x[j] * tg_xx_xxxz_0[j] + tg_xx_xxxxz_0[j];

                    tg_xxx_xxyy_0[j] = -cd_x[j] * tg_xx_xxyy_0[j] + tg_xx_xxxyy_0[j];

                    tg_xxx_xxyz_0[j] = -cd_x[j] * tg_xx_xxyz_0[j] + tg_xx_xxxyz_0[j];

                    tg_xxx_xxzz_0[j] = -cd_x[j] * tg_xx_xxzz_0[j] + tg_xx_xxxzz_0[j];

                    tg_xxx_xyyy_0[j] = -cd_x[j] * tg_xx_xyyy_0[j] + tg_xx_xxyyy_0[j];

                    tg_xxx_xyyz_0[j] = -cd_x[j] * tg_xx_xyyz_0[j] + tg_xx_xxyyz_0[j];

                    tg_xxx_xyzz_0[j] = -cd_x[j] * tg_xx_xyzz_0[j] + tg_xx_xxyzz_0[j];

                    tg_xxx_xzzz_0[j] = -cd_x[j] * tg_xx_xzzz_0[j] + tg_xx_xxzzz_0[j];

                    tg_xxx_yyyy_0[j] = -cd_x[j] * tg_xx_yyyy_0[j] + tg_xx_xyyyy_0[j];

                    tg_xxx_yyyz_0[j] = -cd_x[j] * tg_xx_yyyz_0[j] + tg_xx_xyyyz_0[j];

                    tg_xxx_yyzz_0[j] = -cd_x[j] * tg_xx_yyzz_0[j] + tg_xx_xyyzz_0[j];

                    tg_xxx_yzzz_0[j] = -cd_x[j] * tg_xx_yzzz_0[j] + tg_xx_xyzzz_0[j];

                    tg_xxx_zzzz_0[j] = -cd_x[j] * tg_xx_zzzz_0[j] + tg_xx_xzzzz_0[j];

                    tg_xxy_xxxx_0[j] = -cd_x[j] * tg_xy_xxxx_0[j] + tg_xy_xxxxx_0[j];

                    tg_xxy_xxxy_0[j] = -cd_x[j] * tg_xy_xxxy_0[j] + tg_xy_xxxxy_0[j];

                    tg_xxy_xxxz_0[j] = -cd_x[j] * tg_xy_xxxz_0[j] + tg_xy_xxxxz_0[j];

                    tg_xxy_xxyy_0[j] = -cd_x[j] * tg_xy_xxyy_0[j] + tg_xy_xxxyy_0[j];

                    tg_xxy_xxyz_0[j] = -cd_x[j] * tg_xy_xxyz_0[j] + tg_xy_xxxyz_0[j];

                    tg_xxy_xxzz_0[j] = -cd_x[j] * tg_xy_xxzz_0[j] + tg_xy_xxxzz_0[j];

                    tg_xxy_xyyy_0[j] = -cd_x[j] * tg_xy_xyyy_0[j] + tg_xy_xxyyy_0[j];

                    tg_xxy_xyyz_0[j] = -cd_x[j] * tg_xy_xyyz_0[j] + tg_xy_xxyyz_0[j];

                    tg_xxy_xyzz_0[j] = -cd_x[j] * tg_xy_xyzz_0[j] + tg_xy_xxyzz_0[j];

                    tg_xxy_xzzz_0[j] = -cd_x[j] * tg_xy_xzzz_0[j] + tg_xy_xxzzz_0[j];

                    tg_xxy_yyyy_0[j] = -cd_x[j] * tg_xy_yyyy_0[j] + tg_xy_xyyyy_0[j];

                    tg_xxy_yyyz_0[j] = -cd_x[j] * tg_xy_yyyz_0[j] + tg_xy_xyyyz_0[j];

                    tg_xxy_yyzz_0[j] = -cd_x[j] * tg_xy_yyzz_0[j] + tg_xy_xyyzz_0[j];

                    tg_xxy_yzzz_0[j] = -cd_x[j] * tg_xy_yzzz_0[j] + tg_xy_xyzzz_0[j];

                    tg_xxy_zzzz_0[j] = -cd_x[j] * tg_xy_zzzz_0[j] + tg_xy_xzzzz_0[j];

                    tg_xxz_xxxx_0[j] = -cd_x[j] * tg_xz_xxxx_0[j] + tg_xz_xxxxx_0[j];

                    tg_xxz_xxxy_0[j] = -cd_x[j] * tg_xz_xxxy_0[j] + tg_xz_xxxxy_0[j];

                    tg_xxz_xxxz_0[j] = -cd_x[j] * tg_xz_xxxz_0[j] + tg_xz_xxxxz_0[j];

                    tg_xxz_xxyy_0[j] = -cd_x[j] * tg_xz_xxyy_0[j] + tg_xz_xxxyy_0[j];

                    tg_xxz_xxyz_0[j] = -cd_x[j] * tg_xz_xxyz_0[j] + tg_xz_xxxyz_0[j];

                    tg_xxz_xxzz_0[j] = -cd_x[j] * tg_xz_xxzz_0[j] + tg_xz_xxxzz_0[j];

                    tg_xxz_xyyy_0[j] = -cd_x[j] * tg_xz_xyyy_0[j] + tg_xz_xxyyy_0[j];

                    tg_xxz_xyyz_0[j] = -cd_x[j] * tg_xz_xyyz_0[j] + tg_xz_xxyyz_0[j];

                    tg_xxz_xyzz_0[j] = -cd_x[j] * tg_xz_xyzz_0[j] + tg_xz_xxyzz_0[j];

                    tg_xxz_xzzz_0[j] = -cd_x[j] * tg_xz_xzzz_0[j] + tg_xz_xxzzz_0[j];

                    tg_xxz_yyyy_0[j] = -cd_x[j] * tg_xz_yyyy_0[j] + tg_xz_xyyyy_0[j];

                    tg_xxz_yyyz_0[j] = -cd_x[j] * tg_xz_yyyz_0[j] + tg_xz_xyyyz_0[j];

                    tg_xxz_yyzz_0[j] = -cd_x[j] * tg_xz_yyzz_0[j] + tg_xz_xyyzz_0[j];

                    tg_xxz_yzzz_0[j] = -cd_x[j] * tg_xz_yzzz_0[j] + tg_xz_xyzzz_0[j];

                    tg_xxz_zzzz_0[j] = -cd_x[j] * tg_xz_zzzz_0[j] + tg_xz_xzzzz_0[j];

                    tg_xyy_xxxx_0[j] = -cd_x[j] * tg_yy_xxxx_0[j] + tg_yy_xxxxx_0[j];

                    tg_xyy_xxxy_0[j] = -cd_x[j] * tg_yy_xxxy_0[j] + tg_yy_xxxxy_0[j];

                    tg_xyy_xxxz_0[j] = -cd_x[j] * tg_yy_xxxz_0[j] + tg_yy_xxxxz_0[j];

                    tg_xyy_xxyy_0[j] = -cd_x[j] * tg_yy_xxyy_0[j] + tg_yy_xxxyy_0[j];

                    tg_xyy_xxyz_0[j] = -cd_x[j] * tg_yy_xxyz_0[j] + tg_yy_xxxyz_0[j];

                    tg_xyy_xxzz_0[j] = -cd_x[j] * tg_yy_xxzz_0[j] + tg_yy_xxxzz_0[j];

                    tg_xyy_xyyy_0[j] = -cd_x[j] * tg_yy_xyyy_0[j] + tg_yy_xxyyy_0[j];

                    tg_xyy_xyyz_0[j] = -cd_x[j] * tg_yy_xyyz_0[j] + tg_yy_xxyyz_0[j];

                    tg_xyy_xyzz_0[j] = -cd_x[j] * tg_yy_xyzz_0[j] + tg_yy_xxyzz_0[j];

                    tg_xyy_xzzz_0[j] = -cd_x[j] * tg_yy_xzzz_0[j] + tg_yy_xxzzz_0[j];

                    tg_xyy_yyyy_0[j] = -cd_x[j] * tg_yy_yyyy_0[j] + tg_yy_xyyyy_0[j];

                    tg_xyy_yyyz_0[j] = -cd_x[j] * tg_yy_yyyz_0[j] + tg_yy_xyyyz_0[j];

                    tg_xyy_yyzz_0[j] = -cd_x[j] * tg_yy_yyzz_0[j] + tg_yy_xyyzz_0[j];

                    tg_xyy_yzzz_0[j] = -cd_x[j] * tg_yy_yzzz_0[j] + tg_yy_xyzzz_0[j];

                    tg_xyy_zzzz_0[j] = -cd_x[j] * tg_yy_zzzz_0[j] + tg_yy_xzzzz_0[j];

                    tg_xyz_xxxx_0[j] = -cd_x[j] * tg_yz_xxxx_0[j] + tg_yz_xxxxx_0[j];

                    tg_xyz_xxxy_0[j] = -cd_x[j] * tg_yz_xxxy_0[j] + tg_yz_xxxxy_0[j];

                    tg_xyz_xxxz_0[j] = -cd_x[j] * tg_yz_xxxz_0[j] + tg_yz_xxxxz_0[j];

                    tg_xyz_xxyy_0[j] = -cd_x[j] * tg_yz_xxyy_0[j] + tg_yz_xxxyy_0[j];

                    tg_xyz_xxyz_0[j] = -cd_x[j] * tg_yz_xxyz_0[j] + tg_yz_xxxyz_0[j];

                    tg_xyz_xxzz_0[j] = -cd_x[j] * tg_yz_xxzz_0[j] + tg_yz_xxxzz_0[j];

                    tg_xyz_xyyy_0[j] = -cd_x[j] * tg_yz_xyyy_0[j] + tg_yz_xxyyy_0[j];

                    tg_xyz_xyyz_0[j] = -cd_x[j] * tg_yz_xyyz_0[j] + tg_yz_xxyyz_0[j];

                    tg_xyz_xyzz_0[j] = -cd_x[j] * tg_yz_xyzz_0[j] + tg_yz_xxyzz_0[j];

                    tg_xyz_xzzz_0[j] = -cd_x[j] * tg_yz_xzzz_0[j] + tg_yz_xxzzz_0[j];

                    tg_xyz_yyyy_0[j] = -cd_x[j] * tg_yz_yyyy_0[j] + tg_yz_xyyyy_0[j];

                    tg_xyz_yyyz_0[j] = -cd_x[j] * tg_yz_yyyz_0[j] + tg_yz_xyyyz_0[j];

                    tg_xyz_yyzz_0[j] = -cd_x[j] * tg_yz_yyzz_0[j] + tg_yz_xyyzz_0[j];

                    tg_xyz_yzzz_0[j] = -cd_x[j] * tg_yz_yzzz_0[j] + tg_yz_xyzzz_0[j];

                    tg_xyz_zzzz_0[j] = -cd_x[j] * tg_yz_zzzz_0[j] + tg_yz_xzzzz_0[j];
                }
            }
        }
    }

    void
    compElectronRepulsionForSXFG_75_150(      CMemBlock2D<double>& ketBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& cdDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetContrPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (75,150)

        // set up pointers to distances R(CD) = C - D

        auto cd_x = cdDistances.data(0);

        auto cd_y = cdDistances.data(1);

        auto cd_z = cdDistances.data(2);

        // set up bra side loop over intergals

        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();

        for (int32_t iang = 0; iang <= bang; iang++)
        {
            // set up index of integral

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 4, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 4, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 5, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

                auto tg_yy_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 45); 

                auto tg_yy_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 46); 

                auto tg_yy_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 47); 

                auto tg_yy_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 48); 

                auto tg_yy_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 49); 

                auto tg_yy_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 50); 

                auto tg_yy_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 51); 

                auto tg_yy_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 52); 

                auto tg_yy_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 53); 

                auto tg_yy_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 54); 

                auto tg_yy_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 55); 

                auto tg_yy_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 56); 

                auto tg_yy_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 57); 

                auto tg_yy_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 58); 

                auto tg_yy_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 59); 

                auto tg_yz_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 60); 

                auto tg_yz_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 61); 

                auto tg_yz_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 62); 

                auto tg_yz_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 63); 

                auto tg_yz_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 64); 

                auto tg_yz_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 65); 

                auto tg_yz_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 66); 

                auto tg_yz_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 67); 

                auto tg_yz_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 68); 

                auto tg_yz_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 69); 

                auto tg_yz_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 70); 

                auto tg_yz_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 71); 

                auto tg_yz_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 72); 

                auto tg_yz_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 73); 

                auto tg_yz_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 74); 

                auto tg_zz_xxxx_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 75); 

                auto tg_zz_xxxy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 76); 

                auto tg_zz_xxxz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 77); 

                auto tg_zz_xxyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 78); 

                auto tg_zz_xxyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 79); 

                auto tg_zz_xxzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 80); 

                auto tg_zz_xyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 81); 

                auto tg_zz_xyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 82); 

                auto tg_zz_xyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 83); 

                auto tg_zz_xzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 84); 

                auto tg_zz_yyyy_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 85); 

                auto tg_zz_yyyz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 86); 

                auto tg_zz_yyzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 87); 

                auto tg_zz_yzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 88); 

                auto tg_zz_zzzz_0 = ketBuffer.data(pidx_g_2_4_m0 + 90 * i + 89); 

                auto tg_yy_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 64); 

                auto tg_yy_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 66); 

                auto tg_yy_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 67); 

                auto tg_yy_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 69); 

                auto tg_yy_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 70); 

                auto tg_yy_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 71); 

                auto tg_yy_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 73); 

                auto tg_yy_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 74); 

                auto tg_yy_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 75); 

                auto tg_yy_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 76); 

                auto tg_yy_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 78); 

                auto tg_yy_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 79); 

                auto tg_yy_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 80); 

                auto tg_yy_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 81); 

                auto tg_yy_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 82); 

                auto tg_yz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 85); 

                auto tg_yz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 87); 

                auto tg_yz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 88); 

                auto tg_yz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 90); 

                auto tg_yz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 91); 

                auto tg_yz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 92); 

                auto tg_yz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 94); 

                auto tg_yz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 95); 

                auto tg_yz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 96); 

                auto tg_yz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 97); 

                auto tg_yz_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 99); 

                auto tg_yz_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 100); 

                auto tg_yz_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 101); 

                auto tg_yz_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 102); 

                auto tg_yz_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 103); 

                auto tg_zz_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 105); 

                auto tg_zz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 106); 

                auto tg_zz_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 107); 

                auto tg_zz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 108); 

                auto tg_zz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 109); 

                auto tg_zz_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 110); 

                auto tg_zz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 111); 

                auto tg_zz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 112); 

                auto tg_zz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 113); 

                auto tg_zz_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 114); 

                auto tg_zz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 115); 

                auto tg_zz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 116); 

                auto tg_zz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 117); 

                auto tg_zz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 118); 

                auto tg_zz_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 119); 

                auto tg_zz_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 120); 

                auto tg_zz_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 121); 

                auto tg_zz_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 122); 

                auto tg_zz_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 123); 

                auto tg_zz_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 124); 

                auto tg_zz_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 125); 

                // set up pointers to integrals

                auto tg_xzz_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 75); 

                auto tg_xzz_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 76); 

                auto tg_xzz_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 77); 

                auto tg_xzz_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 78); 

                auto tg_xzz_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 79); 

                auto tg_xzz_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 80); 

                auto tg_xzz_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 81); 

                auto tg_xzz_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 82); 

                auto tg_xzz_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 83); 

                auto tg_xzz_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 84); 

                auto tg_xzz_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 85); 

                auto tg_xzz_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 86); 

                auto tg_xzz_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 87); 

                auto tg_xzz_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 88); 

                auto tg_xzz_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 89); 

                auto tg_yyy_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 90); 

                auto tg_yyy_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 91); 

                auto tg_yyy_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 92); 

                auto tg_yyy_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 93); 

                auto tg_yyy_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 94); 

                auto tg_yyy_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 95); 

                auto tg_yyy_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 96); 

                auto tg_yyy_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 97); 

                auto tg_yyy_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 98); 

                auto tg_yyy_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 99); 

                auto tg_yyy_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 100); 

                auto tg_yyy_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 101); 

                auto tg_yyy_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 102); 

                auto tg_yyy_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 103); 

                auto tg_yyy_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 104); 

                auto tg_yyz_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 105); 

                auto tg_yyz_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 106); 

                auto tg_yyz_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 107); 

                auto tg_yyz_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 108); 

                auto tg_yyz_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 109); 

                auto tg_yyz_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 110); 

                auto tg_yyz_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 111); 

                auto tg_yyz_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 112); 

                auto tg_yyz_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 113); 

                auto tg_yyz_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 114); 

                auto tg_yyz_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 115); 

                auto tg_yyz_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 116); 

                auto tg_yyz_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 117); 

                auto tg_yyz_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 118); 

                auto tg_yyz_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 119); 

                auto tg_yzz_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 120); 

                auto tg_yzz_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 121); 

                auto tg_yzz_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 122); 

                auto tg_yzz_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 123); 

                auto tg_yzz_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 124); 

                auto tg_yzz_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 125); 

                auto tg_yzz_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 126); 

                auto tg_yzz_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 127); 

                auto tg_yzz_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 128); 

                auto tg_yzz_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 129); 

                auto tg_yzz_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 130); 

                auto tg_yzz_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 131); 

                auto tg_yzz_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 132); 

                auto tg_yzz_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 133); 

                auto tg_yzz_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 134); 

                auto tg_zzz_xxxx_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 135); 

                auto tg_zzz_xxxy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 136); 

                auto tg_zzz_xxxz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 137); 

                auto tg_zzz_xxyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 138); 

                auto tg_zzz_xxyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 139); 

                auto tg_zzz_xxzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 140); 

                auto tg_zzz_xyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 141); 

                auto tg_zzz_xyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 142); 

                auto tg_zzz_xyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 143); 

                auto tg_zzz_xzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 144); 

                auto tg_zzz_yyyy_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 145); 

                auto tg_zzz_yyyz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 146); 

                auto tg_zzz_yyzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 147); 

                auto tg_zzz_yzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 148); 

                auto tg_zzz_zzzz_0 = ketBuffer.data(pidx_g_3_4_m0 + 150 * i + 149); 

                // Batch of Integrals (75,150)

                #pragma omp simd aligned(cd_x, cd_y, cd_z, tg_xzz_xxxx_0, tg_xzz_xxxy_0, tg_xzz_xxxz_0, \
                                         tg_xzz_xxyy_0, tg_xzz_xxyz_0, tg_xzz_xxzz_0, tg_xzz_xyyy_0, tg_xzz_xyyz_0, \
                                         tg_xzz_xyzz_0, tg_xzz_xzzz_0, tg_xzz_yyyy_0, tg_xzz_yyyz_0, tg_xzz_yyzz_0, \
                                         tg_xzz_yzzz_0, tg_xzz_zzzz_0, tg_yy_xxxx_0, tg_yy_xxxxy_0, tg_yy_xxxy_0, \
                                         tg_yy_xxxyy_0, tg_yy_xxxyz_0, tg_yy_xxxz_0, tg_yy_xxyy_0, tg_yy_xxyyy_0, \
                                         tg_yy_xxyyz_0, tg_yy_xxyz_0, tg_yy_xxyzz_0, tg_yy_xxzz_0, tg_yy_xyyy_0, \
                                         tg_yy_xyyyy_0, tg_yy_xyyyz_0, tg_yy_xyyz_0, tg_yy_xyyzz_0, tg_yy_xyzz_0, \
                                         tg_yy_xyzzz_0, tg_yy_xzzz_0, tg_yy_yyyy_0, tg_yy_yyyyy_0, tg_yy_yyyyz_0, \
                                         tg_yy_yyyz_0, tg_yy_yyyzz_0, tg_yy_yyzz_0, tg_yy_yyzzz_0, tg_yy_yzzz_0, \
                                         tg_yy_yzzzz_0, tg_yy_zzzz_0, tg_yyy_xxxx_0, tg_yyy_xxxy_0, tg_yyy_xxxz_0, \
                                         tg_yyy_xxyy_0, tg_yyy_xxyz_0, tg_yyy_xxzz_0, tg_yyy_xyyy_0, tg_yyy_xyyz_0, \
                                         tg_yyy_xyzz_0, tg_yyy_xzzz_0, tg_yyy_yyyy_0, tg_yyy_yyyz_0, tg_yyy_yyzz_0, \
                                         tg_yyy_yzzz_0, tg_yyy_zzzz_0, tg_yyz_xxxx_0, tg_yyz_xxxy_0, tg_yyz_xxxz_0, \
                                         tg_yyz_xxyy_0, tg_yyz_xxyz_0, tg_yyz_xxzz_0, tg_yyz_xyyy_0, tg_yyz_xyyz_0, \
                                         tg_yyz_xyzz_0, tg_yyz_xzzz_0, tg_yyz_yyyy_0, tg_yyz_yyyz_0, tg_yyz_yyzz_0, \
                                         tg_yyz_yzzz_0, tg_yyz_zzzz_0, tg_yz_xxxx_0, tg_yz_xxxxy_0, tg_yz_xxxy_0, \
                                         tg_yz_xxxyy_0, tg_yz_xxxyz_0, tg_yz_xxxz_0, tg_yz_xxyy_0, tg_yz_xxyyy_0, \
                                         tg_yz_xxyyz_0, tg_yz_xxyz_0, tg_yz_xxyzz_0, tg_yz_xxzz_0, tg_yz_xyyy_0, \
                                         tg_yz_xyyyy_0, tg_yz_xyyyz_0, tg_yz_xyyz_0, tg_yz_xyyzz_0, tg_yz_xyzz_0, \
                                         tg_yz_xyzzz_0, tg_yz_xzzz_0, tg_yz_yyyy_0, tg_yz_yyyyy_0, tg_yz_yyyyz_0, \
                                         tg_yz_yyyz_0, tg_yz_yyyzz_0, tg_yz_yyzz_0, tg_yz_yyzzz_0, tg_yz_yzzz_0, \
                                         tg_yz_yzzzz_0, tg_yz_zzzz_0, tg_yzz_xxxx_0, tg_yzz_xxxy_0, tg_yzz_xxxz_0, \
                                         tg_yzz_xxyy_0, tg_yzz_xxyz_0, tg_yzz_xxzz_0, tg_yzz_xyyy_0, tg_yzz_xyyz_0, \
                                         tg_yzz_xyzz_0, tg_yzz_xzzz_0, tg_yzz_yyyy_0, tg_yzz_yyyz_0, tg_yzz_yyzz_0, \
                                         tg_yzz_yzzz_0, tg_yzz_zzzz_0, tg_zz_xxxx_0, tg_zz_xxxxx_0, tg_zz_xxxxy_0, \
                                         tg_zz_xxxxz_0, tg_zz_xxxy_0, tg_zz_xxxyy_0, tg_zz_xxxyz_0, tg_zz_xxxz_0, \
                                         tg_zz_xxxzz_0, tg_zz_xxyy_0, tg_zz_xxyyy_0, tg_zz_xxyyz_0, tg_zz_xxyz_0, \
                                         tg_zz_xxyzz_0, tg_zz_xxzz_0, tg_zz_xxzzz_0, tg_zz_xyyy_0, tg_zz_xyyyy_0, \
                                         tg_zz_xyyyz_0, tg_zz_xyyz_0, tg_zz_xyyzz_0, tg_zz_xyzz_0, tg_zz_xyzzz_0, \
                                         tg_zz_xzzz_0, tg_zz_xzzzz_0, tg_zz_yyyy_0, tg_zz_yyyyy_0, tg_zz_yyyyz_0, \
                                         tg_zz_yyyz_0, tg_zz_yyyzz_0, tg_zz_yyzz_0, tg_zz_yyzzz_0, tg_zz_yzzz_0, \
                                         tg_zz_yzzzz_0, tg_zz_zzzz_0, tg_zz_zzzzz_0, tg_zzz_xxxx_0, tg_zzz_xxxy_0, \
                                         tg_zzz_xxxz_0, tg_zzz_xxyy_0, tg_zzz_xxyz_0, tg_zzz_xxzz_0, tg_zzz_xyyy_0, \
                                         tg_zzz_xyyz_0, tg_zzz_xyzz_0, tg_zzz_xzzz_0, tg_zzz_yyyy_0, tg_zzz_yyyz_0, \
                                         tg_zzz_yyzz_0, tg_zzz_yzzz_0, tg_zzz_zzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nKetContrPairs; j++)
                {
                    tg_xzz_xxxx_0[j] = -cd_x[j] * tg_zz_xxxx_0[j] + tg_zz_xxxxx_0[j];

                    tg_xzz_xxxy_0[j] = -cd_x[j] * tg_zz_xxxy_0[j] + tg_zz_xxxxy_0[j];

                    tg_xzz_xxxz_0[j] = -cd_x[j] * tg_zz_xxxz_0[j] + tg_zz_xxxxz_0[j];

                    tg_xzz_xxyy_0[j] = -cd_x[j] * tg_zz_xxyy_0[j] + tg_zz_xxxyy_0[j];

                    tg_xzz_xxyz_0[j] = -cd_x[j] * tg_zz_xxyz_0[j] + tg_zz_xxxyz_0[j];

                    tg_xzz_xxzz_0[j] = -cd_x[j] * tg_zz_xxzz_0[j] + tg_zz_xxxzz_0[j];

                    tg_xzz_xyyy_0[j] = -cd_x[j] * tg_zz_xyyy_0[j] + tg_zz_xxyyy_0[j];

                    tg_xzz_xyyz_0[j] = -cd_x[j] * tg_zz_xyyz_0[j] + tg_zz_xxyyz_0[j];

                    tg_xzz_xyzz_0[j] = -cd_x[j] * tg_zz_xyzz_0[j] + tg_zz_xxyzz_0[j];

                    tg_xzz_xzzz_0[j] = -cd_x[j] * tg_zz_xzzz_0[j] + tg_zz_xxzzz_0[j];

                    tg_xzz_yyyy_0[j] = -cd_x[j] * tg_zz_yyyy_0[j] + tg_zz_xyyyy_0[j];

                    tg_xzz_yyyz_0[j] = -cd_x[j] * tg_zz_yyyz_0[j] + tg_zz_xyyyz_0[j];

                    tg_xzz_yyzz_0[j] = -cd_x[j] * tg_zz_yyzz_0[j] + tg_zz_xyyzz_0[j];

                    tg_xzz_yzzz_0[j] = -cd_x[j] * tg_zz_yzzz_0[j] + tg_zz_xyzzz_0[j];

                    tg_xzz_zzzz_0[j] = -cd_x[j] * tg_zz_zzzz_0[j] + tg_zz_xzzzz_0[j];

                    tg_yyy_xxxx_0[j] = -cd_y[j] * tg_yy_xxxx_0[j] + tg_yy_xxxxy_0[j];

                    tg_yyy_xxxy_0[j] = -cd_y[j] * tg_yy_xxxy_0[j] + tg_yy_xxxyy_0[j];

                    tg_yyy_xxxz_0[j] = -cd_y[j] * tg_yy_xxxz_0[j] + tg_yy_xxxyz_0[j];

                    tg_yyy_xxyy_0[j] = -cd_y[j] * tg_yy_xxyy_0[j] + tg_yy_xxyyy_0[j];

                    tg_yyy_xxyz_0[j] = -cd_y[j] * tg_yy_xxyz_0[j] + tg_yy_xxyyz_0[j];

                    tg_yyy_xxzz_0[j] = -cd_y[j] * tg_yy_xxzz_0[j] + tg_yy_xxyzz_0[j];

                    tg_yyy_xyyy_0[j] = -cd_y[j] * tg_yy_xyyy_0[j] + tg_yy_xyyyy_0[j];

                    tg_yyy_xyyz_0[j] = -cd_y[j] * tg_yy_xyyz_0[j] + tg_yy_xyyyz_0[j];

                    tg_yyy_xyzz_0[j] = -cd_y[j] * tg_yy_xyzz_0[j] + tg_yy_xyyzz_0[j];

                    tg_yyy_xzzz_0[j] = -cd_y[j] * tg_yy_xzzz_0[j] + tg_yy_xyzzz_0[j];

                    tg_yyy_yyyy_0[j] = -cd_y[j] * tg_yy_yyyy_0[j] + tg_yy_yyyyy_0[j];

                    tg_yyy_yyyz_0[j] = -cd_y[j] * tg_yy_yyyz_0[j] + tg_yy_yyyyz_0[j];

                    tg_yyy_yyzz_0[j] = -cd_y[j] * tg_yy_yyzz_0[j] + tg_yy_yyyzz_0[j];

                    tg_yyy_yzzz_0[j] = -cd_y[j] * tg_yy_yzzz_0[j] + tg_yy_yyzzz_0[j];

                    tg_yyy_zzzz_0[j] = -cd_y[j] * tg_yy_zzzz_0[j] + tg_yy_yzzzz_0[j];

                    tg_yyz_xxxx_0[j] = -cd_y[j] * tg_yz_xxxx_0[j] + tg_yz_xxxxy_0[j];

                    tg_yyz_xxxy_0[j] = -cd_y[j] * tg_yz_xxxy_0[j] + tg_yz_xxxyy_0[j];

                    tg_yyz_xxxz_0[j] = -cd_y[j] * tg_yz_xxxz_0[j] + tg_yz_xxxyz_0[j];

                    tg_yyz_xxyy_0[j] = -cd_y[j] * tg_yz_xxyy_0[j] + tg_yz_xxyyy_0[j];

                    tg_yyz_xxyz_0[j] = -cd_y[j] * tg_yz_xxyz_0[j] + tg_yz_xxyyz_0[j];

                    tg_yyz_xxzz_0[j] = -cd_y[j] * tg_yz_xxzz_0[j] + tg_yz_xxyzz_0[j];

                    tg_yyz_xyyy_0[j] = -cd_y[j] * tg_yz_xyyy_0[j] + tg_yz_xyyyy_0[j];

                    tg_yyz_xyyz_0[j] = -cd_y[j] * tg_yz_xyyz_0[j] + tg_yz_xyyyz_0[j];

                    tg_yyz_xyzz_0[j] = -cd_y[j] * tg_yz_xyzz_0[j] + tg_yz_xyyzz_0[j];

                    tg_yyz_xzzz_0[j] = -cd_y[j] * tg_yz_xzzz_0[j] + tg_yz_xyzzz_0[j];

                    tg_yyz_yyyy_0[j] = -cd_y[j] * tg_yz_yyyy_0[j] + tg_yz_yyyyy_0[j];

                    tg_yyz_yyyz_0[j] = -cd_y[j] * tg_yz_yyyz_0[j] + tg_yz_yyyyz_0[j];

                    tg_yyz_yyzz_0[j] = -cd_y[j] * tg_yz_yyzz_0[j] + tg_yz_yyyzz_0[j];

                    tg_yyz_yzzz_0[j] = -cd_y[j] * tg_yz_yzzz_0[j] + tg_yz_yyzzz_0[j];

                    tg_yyz_zzzz_0[j] = -cd_y[j] * tg_yz_zzzz_0[j] + tg_yz_yzzzz_0[j];

                    tg_yzz_xxxx_0[j] = -cd_y[j] * tg_zz_xxxx_0[j] + tg_zz_xxxxy_0[j];

                    tg_yzz_xxxy_0[j] = -cd_y[j] * tg_zz_xxxy_0[j] + tg_zz_xxxyy_0[j];

                    tg_yzz_xxxz_0[j] = -cd_y[j] * tg_zz_xxxz_0[j] + tg_zz_xxxyz_0[j];

                    tg_yzz_xxyy_0[j] = -cd_y[j] * tg_zz_xxyy_0[j] + tg_zz_xxyyy_0[j];

                    tg_yzz_xxyz_0[j] = -cd_y[j] * tg_zz_xxyz_0[j] + tg_zz_xxyyz_0[j];

                    tg_yzz_xxzz_0[j] = -cd_y[j] * tg_zz_xxzz_0[j] + tg_zz_xxyzz_0[j];

                    tg_yzz_xyyy_0[j] = -cd_y[j] * tg_zz_xyyy_0[j] + tg_zz_xyyyy_0[j];

                    tg_yzz_xyyz_0[j] = -cd_y[j] * tg_zz_xyyz_0[j] + tg_zz_xyyyz_0[j];

                    tg_yzz_xyzz_0[j] = -cd_y[j] * tg_zz_xyzz_0[j] + tg_zz_xyyzz_0[j];

                    tg_yzz_xzzz_0[j] = -cd_y[j] * tg_zz_xzzz_0[j] + tg_zz_xyzzz_0[j];

                    tg_yzz_yyyy_0[j] = -cd_y[j] * tg_zz_yyyy_0[j] + tg_zz_yyyyy_0[j];

                    tg_yzz_yyyz_0[j] = -cd_y[j] * tg_zz_yyyz_0[j] + tg_zz_yyyyz_0[j];

                    tg_yzz_yyzz_0[j] = -cd_y[j] * tg_zz_yyzz_0[j] + tg_zz_yyyzz_0[j];

                    tg_yzz_yzzz_0[j] = -cd_y[j] * tg_zz_yzzz_0[j] + tg_zz_yyzzz_0[j];

                    tg_yzz_zzzz_0[j] = -cd_y[j] * tg_zz_zzzz_0[j] + tg_zz_yzzzz_0[j];

                    tg_zzz_xxxx_0[j] = -cd_z[j] * tg_zz_xxxx_0[j] + tg_zz_xxxxz_0[j];

                    tg_zzz_xxxy_0[j] = -cd_z[j] * tg_zz_xxxy_0[j] + tg_zz_xxxyz_0[j];

                    tg_zzz_xxxz_0[j] = -cd_z[j] * tg_zz_xxxz_0[j] + tg_zz_xxxzz_0[j];

                    tg_zzz_xxyy_0[j] = -cd_z[j] * tg_zz_xxyy_0[j] + tg_zz_xxyyz_0[j];

                    tg_zzz_xxyz_0[j] = -cd_z[j] * tg_zz_xxyz_0[j] + tg_zz_xxyzz_0[j];

                    tg_zzz_xxzz_0[j] = -cd_z[j] * tg_zz_xxzz_0[j] + tg_zz_xxzzz_0[j];

                    tg_zzz_xyyy_0[j] = -cd_z[j] * tg_zz_xyyy_0[j] + tg_zz_xyyyz_0[j];

                    tg_zzz_xyyz_0[j] = -cd_z[j] * tg_zz_xyyz_0[j] + tg_zz_xyyzz_0[j];

                    tg_zzz_xyzz_0[j] = -cd_z[j] * tg_zz_xyzz_0[j] + tg_zz_xyzzz_0[j];

                    tg_zzz_xzzz_0[j] = -cd_z[j] * tg_zz_xzzz_0[j] + tg_zz_xzzzz_0[j];

                    tg_zzz_yyyy_0[j] = -cd_z[j] * tg_zz_yyyy_0[j] + tg_zz_yyyyz_0[j];

                    tg_zzz_yyyz_0[j] = -cd_z[j] * tg_zz_yyyz_0[j] + tg_zz_yyyzz_0[j];

                    tg_zzz_yyzz_0[j] = -cd_z[j] * tg_zz_yyzz_0[j] + tg_zz_yyzzz_0[j];

                    tg_zzz_yzzz_0[j] = -cd_z[j] * tg_zz_yzzz_0[j] + tg_zz_yzzzz_0[j];

                    tg_zzz_zzzz_0[j] = -cd_z[j] * tg_zz_zzzz_0[j] + tg_zz_zzzzz_0[j];
                }
            }
        }
    }

    void
    compElectronRepulsionForSXFH(      CMemBlock2D<double>& ketBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& cdDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        erikrrfunc::compElectronRepulsionForSXFH_0_70(ketBuffer,
                                                      recursionMap,
                                                      cdDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        erikrrfunc::compElectronRepulsionForSXFH_70_140(ketBuffer,
                                                        recursionMap,
                                                        cdDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 

        erikrrfunc::compElectronRepulsionForSXFH_140_210(ketBuffer,
                                                         recursionMap,
                                                         cdDistances,
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetContrPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSXFH_0_70(      CMemBlock2D<double>& ketBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& cdDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,70)

        // set up pointers to distances R(CD) = C - D

        auto cd_x = cdDistances.data(0);

        // set up bra side loop over intergals

        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();

        for (int32_t iang = 0; iang <= bang; iang++)
        {
            // set up index of integral

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 5, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 5, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_2_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 6, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

                auto tg_xx_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i); 

                auto tg_xx_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 1); 

                auto tg_xx_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 2); 

                auto tg_xx_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 3); 

                auto tg_xx_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 4); 

                auto tg_xx_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 5); 

                auto tg_xx_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 6); 

                auto tg_xx_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 7); 

                auto tg_xx_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 8); 

                auto tg_xx_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 9); 

                auto tg_xx_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 10); 

                auto tg_xx_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 11); 

                auto tg_xx_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 12); 

                auto tg_xx_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 13); 

                auto tg_xx_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 14); 

                auto tg_xx_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 15); 

                auto tg_xx_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 16); 

                auto tg_xx_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 17); 

                auto tg_xx_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 18); 

                auto tg_xx_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 19); 

                auto tg_xx_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 20); 

                auto tg_xy_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 21); 

                auto tg_xy_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 22); 

                auto tg_xy_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 23); 

                auto tg_xy_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 24); 

                auto tg_xy_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 25); 

                auto tg_xy_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 26); 

                auto tg_xy_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 27); 

                auto tg_xy_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 28); 

                auto tg_xy_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 29); 

                auto tg_xy_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 30); 

                auto tg_xy_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 31); 

                auto tg_xy_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 32); 

                auto tg_xy_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 33); 

                auto tg_xy_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 34); 

                auto tg_xy_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 35); 

                auto tg_xy_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 36); 

                auto tg_xy_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 37); 

                auto tg_xy_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 38); 

                auto tg_xy_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 39); 

                auto tg_xy_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 40); 

                auto tg_xy_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 41); 

                auto tg_xz_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 42); 

                auto tg_xz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 43); 

                auto tg_xz_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 44); 

                auto tg_xz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 45); 

                auto tg_xz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 46); 

                auto tg_xz_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 47); 

                auto tg_xz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 48); 

                auto tg_xz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 49); 

                auto tg_xz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 50); 

                auto tg_xz_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 51); 

                auto tg_xz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 52); 

                auto tg_xz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 53); 

                auto tg_xz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 54); 

                auto tg_xz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 55); 

                auto tg_xz_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 56); 

                auto tg_xz_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 57); 

                auto tg_xz_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 58); 

                auto tg_xz_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 59); 

                auto tg_xz_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 60); 

                auto tg_xz_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 61); 

                auto tg_xz_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 62); 

                auto tg_yy_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 63); 

                auto tg_yy_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 64); 

                auto tg_yy_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 65); 

                auto tg_yy_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 66); 

                auto tg_yy_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 67); 

                auto tg_yy_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 68); 

                auto tg_yy_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 69); 

                auto tg_xx_xxxxxx_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i); 

                auto tg_xx_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 1); 

                auto tg_xx_xxxxxz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 2); 

                auto tg_xx_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 3); 

                auto tg_xx_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 4); 

                auto tg_xx_xxxxzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 5); 

                auto tg_xx_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 6); 

                auto tg_xx_xxxyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 7); 

                auto tg_xx_xxxyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 8); 

                auto tg_xx_xxxzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 9); 

                auto tg_xx_xxyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 10); 

                auto tg_xx_xxyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 11); 

                auto tg_xx_xxyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 12); 

                auto tg_xx_xxyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 13); 

                auto tg_xx_xxzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 14); 

                auto tg_xx_xyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 15); 

                auto tg_xx_xyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 16); 

                auto tg_xx_xyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 17); 

                auto tg_xx_xyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 18); 

                auto tg_xx_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 19); 

                auto tg_xx_xzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 20); 

                auto tg_xy_xxxxxx_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 28); 

                auto tg_xy_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 29); 

                auto tg_xy_xxxxxz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 30); 

                auto tg_xy_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 31); 

                auto tg_xy_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 32); 

                auto tg_xy_xxxxzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 33); 

                auto tg_xy_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 34); 

                auto tg_xy_xxxyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 35); 

                auto tg_xy_xxxyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 36); 

                auto tg_xy_xxxzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 37); 

                auto tg_xy_xxyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 38); 

                auto tg_xy_xxyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 39); 

                auto tg_xy_xxyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 40); 

                auto tg_xy_xxyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 41); 

                auto tg_xy_xxzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 42); 

                auto tg_xy_xyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 43); 

                auto tg_xy_xyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 44); 

                auto tg_xy_xyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 45); 

                auto tg_xy_xyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 46); 

                auto tg_xy_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 47); 

                auto tg_xy_xzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 48); 

                auto tg_xz_xxxxxx_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 56); 

                auto tg_xz_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 57); 

                auto tg_xz_xxxxxz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 58); 

                auto tg_xz_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 59); 

                auto tg_xz_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 60); 

                auto tg_xz_xxxxzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 61); 

                auto tg_xz_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 62); 

                auto tg_xz_xxxyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 63); 

                auto tg_xz_xxxyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 64); 

                auto tg_xz_xxxzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 65); 

                auto tg_xz_xxyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 66); 

                auto tg_xz_xxyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 67); 

                auto tg_xz_xxyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 68); 

                auto tg_xz_xxyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 69); 

                auto tg_xz_xxzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 70); 

                auto tg_xz_xyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 71); 

                auto tg_xz_xyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 72); 

                auto tg_xz_xyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 73); 

                auto tg_xz_xyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 74); 

                auto tg_xz_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 75); 

                auto tg_xz_xzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 76); 

                auto tg_yy_xxxxxx_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 84); 

                auto tg_yy_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 85); 

                auto tg_yy_xxxxxz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 86); 

                auto tg_yy_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 87); 

                auto tg_yy_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 88); 

                auto tg_yy_xxxxzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 89); 

                auto tg_yy_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 90); 

                // set up pointers to integrals

                auto tg_xxx_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i); 

                auto tg_xxx_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 1); 

                auto tg_xxx_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 2); 

                auto tg_xxx_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 3); 

                auto tg_xxx_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 4); 

                auto tg_xxx_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 5); 

                auto tg_xxx_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 6); 

                auto tg_xxx_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 7); 

                auto tg_xxx_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 8); 

                auto tg_xxx_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 9); 

                auto tg_xxx_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 10); 

                auto tg_xxx_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 11); 

                auto tg_xxx_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 12); 

                auto tg_xxx_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 13); 

                auto tg_xxx_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 14); 

                auto tg_xxx_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 15); 

                auto tg_xxx_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 16); 

                auto tg_xxx_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 17); 

                auto tg_xxx_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 18); 

                auto tg_xxx_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 19); 

                auto tg_xxx_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 20); 

                auto tg_xxy_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 21); 

                auto tg_xxy_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 22); 

                auto tg_xxy_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 23); 

                auto tg_xxy_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 24); 

                auto tg_xxy_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 25); 

                auto tg_xxy_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 26); 

                auto tg_xxy_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 27); 

                auto tg_xxy_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 28); 

                auto tg_xxy_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 29); 

                auto tg_xxy_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 30); 

                auto tg_xxy_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 31); 

                auto tg_xxy_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 32); 

                auto tg_xxy_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 33); 

                auto tg_xxy_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 34); 

                auto tg_xxy_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 35); 

                auto tg_xxy_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 36); 

                auto tg_xxy_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 37); 

                auto tg_xxy_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 38); 

                auto tg_xxy_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 39); 

                auto tg_xxy_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 40); 

                auto tg_xxy_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 41); 

                auto tg_xxz_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 42); 

                auto tg_xxz_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 43); 

                auto tg_xxz_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 44); 

                auto tg_xxz_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 45); 

                auto tg_xxz_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 46); 

                auto tg_xxz_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 47); 

                auto tg_xxz_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 48); 

                auto tg_xxz_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 49); 

                auto tg_xxz_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 50); 

                auto tg_xxz_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 51); 

                auto tg_xxz_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 52); 

                auto tg_xxz_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 53); 

                auto tg_xxz_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 54); 

                auto tg_xxz_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 55); 

                auto tg_xxz_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 56); 

                auto tg_xxz_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 57); 

                auto tg_xxz_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 58); 

                auto tg_xxz_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 59); 

                auto tg_xxz_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 60); 

                auto tg_xxz_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 61); 

                auto tg_xxz_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 62); 

                auto tg_xyy_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 63); 

                auto tg_xyy_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 64); 

                auto tg_xyy_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 65); 

                auto tg_xyy_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 66); 

                auto tg_xyy_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 67); 

                auto tg_xyy_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 68); 

                auto tg_xyy_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 69); 

                // Batch of Integrals (0,70)

                #pragma omp simd aligned(cd_x, tg_xx_xxxxx_0, tg_xx_xxxxxx_0, tg_xx_xxxxxy_0, tg_xx_xxxxxz_0, \
                                         tg_xx_xxxxy_0, tg_xx_xxxxyy_0, tg_xx_xxxxyz_0, tg_xx_xxxxz_0, tg_xx_xxxxzz_0, \
                                         tg_xx_xxxyy_0, tg_xx_xxxyyy_0, tg_xx_xxxyyz_0, tg_xx_xxxyz_0, tg_xx_xxxyzz_0, \
                                         tg_xx_xxxzz_0, tg_xx_xxxzzz_0, tg_xx_xxyyy_0, tg_xx_xxyyyy_0, tg_xx_xxyyyz_0, \
                                         tg_xx_xxyyz_0, tg_xx_xxyyzz_0, tg_xx_xxyzz_0, tg_xx_xxyzzz_0, tg_xx_xxzzz_0, \
                                         tg_xx_xxzzzz_0, tg_xx_xyyyy_0, tg_xx_xyyyyy_0, tg_xx_xyyyyz_0, tg_xx_xyyyz_0, \
                                         tg_xx_xyyyzz_0, tg_xx_xyyzz_0, tg_xx_xyyzzz_0, tg_xx_xyzzz_0, tg_xx_xyzzzz_0, \
                                         tg_xx_xzzzz_0, tg_xx_xzzzzz_0, tg_xx_yyyyy_0, tg_xx_yyyyz_0, tg_xx_yyyzz_0, \
                                         tg_xx_yyzzz_0, tg_xx_yzzzz_0, tg_xx_zzzzz_0, tg_xxx_xxxxx_0, tg_xxx_xxxxy_0, \
                                         tg_xxx_xxxxz_0, tg_xxx_xxxyy_0, tg_xxx_xxxyz_0, tg_xxx_xxxzz_0, tg_xxx_xxyyy_0, \
                                         tg_xxx_xxyyz_0, tg_xxx_xxyzz_0, tg_xxx_xxzzz_0, tg_xxx_xyyyy_0, tg_xxx_xyyyz_0, \
                                         tg_xxx_xyyzz_0, tg_xxx_xyzzz_0, tg_xxx_xzzzz_0, tg_xxx_yyyyy_0, tg_xxx_yyyyz_0, \
                                         tg_xxx_yyyzz_0, tg_xxx_yyzzz_0, tg_xxx_yzzzz_0, tg_xxx_zzzzz_0, tg_xxy_xxxxx_0, \
                                         tg_xxy_xxxxy_0, tg_xxy_xxxxz_0, tg_xxy_xxxyy_0, tg_xxy_xxxyz_0, tg_xxy_xxxzz_0, \
                                         tg_xxy_xxyyy_0, tg_xxy_xxyyz_0, tg_xxy_xxyzz_0, tg_xxy_xxzzz_0, tg_xxy_xyyyy_0, \
                                         tg_xxy_xyyyz_0, tg_xxy_xyyzz_0, tg_xxy_xyzzz_0, tg_xxy_xzzzz_0, tg_xxy_yyyyy_0, \
                                         tg_xxy_yyyyz_0, tg_xxy_yyyzz_0, tg_xxy_yyzzz_0, tg_xxy_yzzzz_0, tg_xxy_zzzzz_0, \
                                         tg_xxz_xxxxx_0, tg_xxz_xxxxy_0, tg_xxz_xxxxz_0, tg_xxz_xxxyy_0, tg_xxz_xxxyz_0, \
                                         tg_xxz_xxxzz_0, tg_xxz_xxyyy_0, tg_xxz_xxyyz_0, tg_xxz_xxyzz_0, tg_xxz_xxzzz_0, \
                                         tg_xxz_xyyyy_0, tg_xxz_xyyyz_0, tg_xxz_xyyzz_0, tg_xxz_xyzzz_0, tg_xxz_xzzzz_0, \
                                         tg_xxz_yyyyy_0, tg_xxz_yyyyz_0, tg_xxz_yyyzz_0, tg_xxz_yyzzz_0, tg_xxz_yzzzz_0, \
                                         tg_xxz_zzzzz_0, tg_xy_xxxxx_0, tg_xy_xxxxxx_0, tg_xy_xxxxxy_0, tg_xy_xxxxxz_0, \
                                         tg_xy_xxxxy_0, tg_xy_xxxxyy_0, tg_xy_xxxxyz_0, tg_xy_xxxxz_0, tg_xy_xxxxzz_0, \
                                         tg_xy_xxxyy_0, tg_xy_xxxyyy_0, tg_xy_xxxyyz_0, tg_xy_xxxyz_0, tg_xy_xxxyzz_0, \
                                         tg_xy_xxxzz_0, tg_xy_xxxzzz_0, tg_xy_xxyyy_0, tg_xy_xxyyyy_0, tg_xy_xxyyyz_0, \
                                         tg_xy_xxyyz_0, tg_xy_xxyyzz_0, tg_xy_xxyzz_0, tg_xy_xxyzzz_0, tg_xy_xxzzz_0, \
                                         tg_xy_xxzzzz_0, tg_xy_xyyyy_0, tg_xy_xyyyyy_0, tg_xy_xyyyyz_0, tg_xy_xyyyz_0, \
                                         tg_xy_xyyyzz_0, tg_xy_xyyzz_0, tg_xy_xyyzzz_0, tg_xy_xyzzz_0, tg_xy_xyzzzz_0, \
                                         tg_xy_xzzzz_0, tg_xy_xzzzzz_0, tg_xy_yyyyy_0, tg_xy_yyyyz_0, tg_xy_yyyzz_0, \
                                         tg_xy_yyzzz_0, tg_xy_yzzzz_0, tg_xy_zzzzz_0, tg_xyy_xxxxx_0, tg_xyy_xxxxy_0, \
                                         tg_xyy_xxxxz_0, tg_xyy_xxxyy_0, tg_xyy_xxxyz_0, tg_xyy_xxxzz_0, tg_xyy_xxyyy_0, \
                                         tg_xz_xxxxx_0, tg_xz_xxxxxx_0, tg_xz_xxxxxy_0, tg_xz_xxxxxz_0, tg_xz_xxxxy_0, \
                                         tg_xz_xxxxyy_0, tg_xz_xxxxyz_0, tg_xz_xxxxz_0, tg_xz_xxxxzz_0, tg_xz_xxxyy_0, \
                                         tg_xz_xxxyyy_0, tg_xz_xxxyyz_0, tg_xz_xxxyz_0, tg_xz_xxxyzz_0, tg_xz_xxxzz_0, \
                                         tg_xz_xxxzzz_0, tg_xz_xxyyy_0, tg_xz_xxyyyy_0, tg_xz_xxyyyz_0, tg_xz_xxyyz_0, \
                                         tg_xz_xxyyzz_0, tg_xz_xxyzz_0, tg_xz_xxyzzz_0, tg_xz_xxzzz_0, tg_xz_xxzzzz_0, \
                                         tg_xz_xyyyy_0, tg_xz_xyyyyy_0, tg_xz_xyyyyz_0, tg_xz_xyyyz_0, tg_xz_xyyyzz_0, \
                                         tg_xz_xyyzz_0, tg_xz_xyyzzz_0, tg_xz_xyzzz_0, tg_xz_xyzzzz_0, tg_xz_xzzzz_0, \
                                         tg_xz_xzzzzz_0, tg_xz_yyyyy_0, tg_xz_yyyyz_0, tg_xz_yyyzz_0, tg_xz_yyzzz_0, \
                                         tg_xz_yzzzz_0, tg_xz_zzzzz_0, tg_yy_xxxxx_0, tg_yy_xxxxxx_0, tg_yy_xxxxxy_0, \
                                         tg_yy_xxxxxz_0, tg_yy_xxxxy_0, tg_yy_xxxxyy_0, tg_yy_xxxxyz_0, tg_yy_xxxxz_0, \
                                         tg_yy_xxxxzz_0, tg_yy_xxxyy_0, tg_yy_xxxyyy_0, tg_yy_xxxyz_0, tg_yy_xxxzz_0, \
                                         tg_yy_xxyyy_0: VLX_ALIGN)
                for (int32_t j = 0; j < nKetContrPairs; j++)
                {
                    tg_xxx_xxxxx_0[j] = -cd_x[j] * tg_xx_xxxxx_0[j] + tg_xx_xxxxxx_0[j];

                    tg_xxx_xxxxy_0[j] = -cd_x[j] * tg_xx_xxxxy_0[j] + tg_xx_xxxxxy_0[j];

                    tg_xxx_xxxxz_0[j] = -cd_x[j] * tg_xx_xxxxz_0[j] + tg_xx_xxxxxz_0[j];

                    tg_xxx_xxxyy_0[j] = -cd_x[j] * tg_xx_xxxyy_0[j] + tg_xx_xxxxyy_0[j];

                    tg_xxx_xxxyz_0[j] = -cd_x[j] * tg_xx_xxxyz_0[j] + tg_xx_xxxxyz_0[j];

                    tg_xxx_xxxzz_0[j] = -cd_x[j] * tg_xx_xxxzz_0[j] + tg_xx_xxxxzz_0[j];

                    tg_xxx_xxyyy_0[j] = -cd_x[j] * tg_xx_xxyyy_0[j] + tg_xx_xxxyyy_0[j];

                    tg_xxx_xxyyz_0[j] = -cd_x[j] * tg_xx_xxyyz_0[j] + tg_xx_xxxyyz_0[j];

                    tg_xxx_xxyzz_0[j] = -cd_x[j] * tg_xx_xxyzz_0[j] + tg_xx_xxxyzz_0[j];

                    tg_xxx_xxzzz_0[j] = -cd_x[j] * tg_xx_xxzzz_0[j] + tg_xx_xxxzzz_0[j];

                    tg_xxx_xyyyy_0[j] = -cd_x[j] * tg_xx_xyyyy_0[j] + tg_xx_xxyyyy_0[j];

                    tg_xxx_xyyyz_0[j] = -cd_x[j] * tg_xx_xyyyz_0[j] + tg_xx_xxyyyz_0[j];

                    tg_xxx_xyyzz_0[j] = -cd_x[j] * tg_xx_xyyzz_0[j] + tg_xx_xxyyzz_0[j];

                    tg_xxx_xyzzz_0[j] = -cd_x[j] * tg_xx_xyzzz_0[j] + tg_xx_xxyzzz_0[j];

                    tg_xxx_xzzzz_0[j] = -cd_x[j] * tg_xx_xzzzz_0[j] + tg_xx_xxzzzz_0[j];

                    tg_xxx_yyyyy_0[j] = -cd_x[j] * tg_xx_yyyyy_0[j] + tg_xx_xyyyyy_0[j];

                    tg_xxx_yyyyz_0[j] = -cd_x[j] * tg_xx_yyyyz_0[j] + tg_xx_xyyyyz_0[j];

                    tg_xxx_yyyzz_0[j] = -cd_x[j] * tg_xx_yyyzz_0[j] + tg_xx_xyyyzz_0[j];

                    tg_xxx_yyzzz_0[j] = -cd_x[j] * tg_xx_yyzzz_0[j] + tg_xx_xyyzzz_0[j];

                    tg_xxx_yzzzz_0[j] = -cd_x[j] * tg_xx_yzzzz_0[j] + tg_xx_xyzzzz_0[j];

                    tg_xxx_zzzzz_0[j] = -cd_x[j] * tg_xx_zzzzz_0[j] + tg_xx_xzzzzz_0[j];

                    tg_xxy_xxxxx_0[j] = -cd_x[j] * tg_xy_xxxxx_0[j] + tg_xy_xxxxxx_0[j];

                    tg_xxy_xxxxy_0[j] = -cd_x[j] * tg_xy_xxxxy_0[j] + tg_xy_xxxxxy_0[j];

                    tg_xxy_xxxxz_0[j] = -cd_x[j] * tg_xy_xxxxz_0[j] + tg_xy_xxxxxz_0[j];

                    tg_xxy_xxxyy_0[j] = -cd_x[j] * tg_xy_xxxyy_0[j] + tg_xy_xxxxyy_0[j];

                    tg_xxy_xxxyz_0[j] = -cd_x[j] * tg_xy_xxxyz_0[j] + tg_xy_xxxxyz_0[j];

                    tg_xxy_xxxzz_0[j] = -cd_x[j] * tg_xy_xxxzz_0[j] + tg_xy_xxxxzz_0[j];

                    tg_xxy_xxyyy_0[j] = -cd_x[j] * tg_xy_xxyyy_0[j] + tg_xy_xxxyyy_0[j];

                    tg_xxy_xxyyz_0[j] = -cd_x[j] * tg_xy_xxyyz_0[j] + tg_xy_xxxyyz_0[j];

                    tg_xxy_xxyzz_0[j] = -cd_x[j] * tg_xy_xxyzz_0[j] + tg_xy_xxxyzz_0[j];

                    tg_xxy_xxzzz_0[j] = -cd_x[j] * tg_xy_xxzzz_0[j] + tg_xy_xxxzzz_0[j];

                    tg_xxy_xyyyy_0[j] = -cd_x[j] * tg_xy_xyyyy_0[j] + tg_xy_xxyyyy_0[j];

                    tg_xxy_xyyyz_0[j] = -cd_x[j] * tg_xy_xyyyz_0[j] + tg_xy_xxyyyz_0[j];

                    tg_xxy_xyyzz_0[j] = -cd_x[j] * tg_xy_xyyzz_0[j] + tg_xy_xxyyzz_0[j];

                    tg_xxy_xyzzz_0[j] = -cd_x[j] * tg_xy_xyzzz_0[j] + tg_xy_xxyzzz_0[j];

                    tg_xxy_xzzzz_0[j] = -cd_x[j] * tg_xy_xzzzz_0[j] + tg_xy_xxzzzz_0[j];

                    tg_xxy_yyyyy_0[j] = -cd_x[j] * tg_xy_yyyyy_0[j] + tg_xy_xyyyyy_0[j];

                    tg_xxy_yyyyz_0[j] = -cd_x[j] * tg_xy_yyyyz_0[j] + tg_xy_xyyyyz_0[j];

                    tg_xxy_yyyzz_0[j] = -cd_x[j] * tg_xy_yyyzz_0[j] + tg_xy_xyyyzz_0[j];

                    tg_xxy_yyzzz_0[j] = -cd_x[j] * tg_xy_yyzzz_0[j] + tg_xy_xyyzzz_0[j];

                    tg_xxy_yzzzz_0[j] = -cd_x[j] * tg_xy_yzzzz_0[j] + tg_xy_xyzzzz_0[j];

                    tg_xxy_zzzzz_0[j] = -cd_x[j] * tg_xy_zzzzz_0[j] + tg_xy_xzzzzz_0[j];

                    tg_xxz_xxxxx_0[j] = -cd_x[j] * tg_xz_xxxxx_0[j] + tg_xz_xxxxxx_0[j];

                    tg_xxz_xxxxy_0[j] = -cd_x[j] * tg_xz_xxxxy_0[j] + tg_xz_xxxxxy_0[j];

                    tg_xxz_xxxxz_0[j] = -cd_x[j] * tg_xz_xxxxz_0[j] + tg_xz_xxxxxz_0[j];

                    tg_xxz_xxxyy_0[j] = -cd_x[j] * tg_xz_xxxyy_0[j] + tg_xz_xxxxyy_0[j];

                    tg_xxz_xxxyz_0[j] = -cd_x[j] * tg_xz_xxxyz_0[j] + tg_xz_xxxxyz_0[j];

                    tg_xxz_xxxzz_0[j] = -cd_x[j] * tg_xz_xxxzz_0[j] + tg_xz_xxxxzz_0[j];

                    tg_xxz_xxyyy_0[j] = -cd_x[j] * tg_xz_xxyyy_0[j] + tg_xz_xxxyyy_0[j];

                    tg_xxz_xxyyz_0[j] = -cd_x[j] * tg_xz_xxyyz_0[j] + tg_xz_xxxyyz_0[j];

                    tg_xxz_xxyzz_0[j] = -cd_x[j] * tg_xz_xxyzz_0[j] + tg_xz_xxxyzz_0[j];

                    tg_xxz_xxzzz_0[j] = -cd_x[j] * tg_xz_xxzzz_0[j] + tg_xz_xxxzzz_0[j];

                    tg_xxz_xyyyy_0[j] = -cd_x[j] * tg_xz_xyyyy_0[j] + tg_xz_xxyyyy_0[j];

                    tg_xxz_xyyyz_0[j] = -cd_x[j] * tg_xz_xyyyz_0[j] + tg_xz_xxyyyz_0[j];

                    tg_xxz_xyyzz_0[j] = -cd_x[j] * tg_xz_xyyzz_0[j] + tg_xz_xxyyzz_0[j];

                    tg_xxz_xyzzz_0[j] = -cd_x[j] * tg_xz_xyzzz_0[j] + tg_xz_xxyzzz_0[j];

                    tg_xxz_xzzzz_0[j] = -cd_x[j] * tg_xz_xzzzz_0[j] + tg_xz_xxzzzz_0[j];

                    tg_xxz_yyyyy_0[j] = -cd_x[j] * tg_xz_yyyyy_0[j] + tg_xz_xyyyyy_0[j];

                    tg_xxz_yyyyz_0[j] = -cd_x[j] * tg_xz_yyyyz_0[j] + tg_xz_xyyyyz_0[j];

                    tg_xxz_yyyzz_0[j] = -cd_x[j] * tg_xz_yyyzz_0[j] + tg_xz_xyyyzz_0[j];

                    tg_xxz_yyzzz_0[j] = -cd_x[j] * tg_xz_yyzzz_0[j] + tg_xz_xyyzzz_0[j];

                    tg_xxz_yzzzz_0[j] = -cd_x[j] * tg_xz_yzzzz_0[j] + tg_xz_xyzzzz_0[j];

                    tg_xxz_zzzzz_0[j] = -cd_x[j] * tg_xz_zzzzz_0[j] + tg_xz_xzzzzz_0[j];

                    tg_xyy_xxxxx_0[j] = -cd_x[j] * tg_yy_xxxxx_0[j] + tg_yy_xxxxxx_0[j];

                    tg_xyy_xxxxy_0[j] = -cd_x[j] * tg_yy_xxxxy_0[j] + tg_yy_xxxxxy_0[j];

                    tg_xyy_xxxxz_0[j] = -cd_x[j] * tg_yy_xxxxz_0[j] + tg_yy_xxxxxz_0[j];

                    tg_xyy_xxxyy_0[j] = -cd_x[j] * tg_yy_xxxyy_0[j] + tg_yy_xxxxyy_0[j];

                    tg_xyy_xxxyz_0[j] = -cd_x[j] * tg_yy_xxxyz_0[j] + tg_yy_xxxxyz_0[j];

                    tg_xyy_xxxzz_0[j] = -cd_x[j] * tg_yy_xxxzz_0[j] + tg_yy_xxxxzz_0[j];

                    tg_xyy_xxyyy_0[j] = -cd_x[j] * tg_yy_xxyyy_0[j] + tg_yy_xxxyyy_0[j];
                }
            }
        }
    }

    void
    compElectronRepulsionForSXFH_70_140(      CMemBlock2D<double>& ketBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& cdDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetContrPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (70,140)

        // set up pointers to distances R(CD) = C - D

        auto cd_x = cdDistances.data(0);

        auto cd_y = cdDistances.data(1);

        // set up bra side loop over intergals

        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();

        for (int32_t iang = 0; iang <= bang; iang++)
        {
            // set up index of integral

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 5, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 5, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_2_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 6, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

                auto tg_yy_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 63); 

                auto tg_yy_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 64); 

                auto tg_yy_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 65); 

                auto tg_yy_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 66); 

                auto tg_yy_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 67); 

                auto tg_yy_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 68); 

                auto tg_yy_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 69); 

                auto tg_yy_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 70); 

                auto tg_yy_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 71); 

                auto tg_yy_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 72); 

                auto tg_yy_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 73); 

                auto tg_yy_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 74); 

                auto tg_yy_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 75); 

                auto tg_yy_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 76); 

                auto tg_yy_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 77); 

                auto tg_yy_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 78); 

                auto tg_yy_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 79); 

                auto tg_yy_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 80); 

                auto tg_yy_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 81); 

                auto tg_yy_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 82); 

                auto tg_yy_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 83); 

                auto tg_yz_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 84); 

                auto tg_yz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 85); 

                auto tg_yz_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 86); 

                auto tg_yz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 87); 

                auto tg_yz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 88); 

                auto tg_yz_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 89); 

                auto tg_yz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 90); 

                auto tg_yz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 91); 

                auto tg_yz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 92); 

                auto tg_yz_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 93); 

                auto tg_yz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 94); 

                auto tg_yz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 95); 

                auto tg_yz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 96); 

                auto tg_yz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 97); 

                auto tg_yz_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 98); 

                auto tg_yz_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 99); 

                auto tg_yz_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 100); 

                auto tg_yz_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 101); 

                auto tg_yz_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 102); 

                auto tg_yz_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 103); 

                auto tg_yz_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 104); 

                auto tg_zz_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 105); 

                auto tg_zz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 106); 

                auto tg_zz_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 107); 

                auto tg_zz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 108); 

                auto tg_zz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 109); 

                auto tg_zz_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 110); 

                auto tg_zz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 111); 

                auto tg_zz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 112); 

                auto tg_zz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 113); 

                auto tg_zz_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 114); 

                auto tg_zz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 115); 

                auto tg_zz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 116); 

                auto tg_zz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 117); 

                auto tg_zz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 118); 

                auto tg_zz_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 119); 

                auto tg_zz_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 120); 

                auto tg_zz_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 121); 

                auto tg_zz_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 122); 

                auto tg_zz_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 123); 

                auto tg_zz_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 124); 

                auto tg_zz_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 125); 

                auto tg_yy_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 85); 

                auto tg_yy_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 87); 

                auto tg_yy_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 88); 

                auto tg_yy_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 90); 

                auto tg_yy_xxxyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 91); 

                auto tg_yy_xxxyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 92); 

                auto tg_yy_xxxzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 93); 

                auto tg_yy_xxyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 94); 

                auto tg_yy_xxyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 95); 

                auto tg_yy_xxyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 96); 

                auto tg_yy_xxyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 97); 

                auto tg_yy_xxzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 98); 

                auto tg_yy_xyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 99); 

                auto tg_yy_xyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 100); 

                auto tg_yy_xyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 101); 

                auto tg_yy_xyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 102); 

                auto tg_yy_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 103); 

                auto tg_yy_xzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 104); 

                auto tg_yz_xxxxxx_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 112); 

                auto tg_yz_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 113); 

                auto tg_yz_xxxxxz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 114); 

                auto tg_yz_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 115); 

                auto tg_yz_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 116); 

                auto tg_yz_xxxxzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 117); 

                auto tg_yz_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 118); 

                auto tg_yz_xxxyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 119); 

                auto tg_yz_xxxyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 120); 

                auto tg_yz_xxxzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 121); 

                auto tg_yz_xxyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 122); 

                auto tg_yz_xxyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 123); 

                auto tg_yz_xxyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 124); 

                auto tg_yz_xxyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 125); 

                auto tg_yz_xxzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 126); 

                auto tg_yz_xyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 127); 

                auto tg_yz_xyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 128); 

                auto tg_yz_xyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 129); 

                auto tg_yz_xyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 130); 

                auto tg_yz_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 131); 

                auto tg_yz_xzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 132); 

                auto tg_zz_xxxxxx_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 140); 

                auto tg_zz_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 141); 

                auto tg_zz_xxxxxz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 142); 

                auto tg_zz_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 143); 

                auto tg_zz_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 144); 

                auto tg_zz_xxxxzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 145); 

                auto tg_zz_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 146); 

                auto tg_zz_xxxyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 147); 

                auto tg_zz_xxxyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 148); 

                auto tg_zz_xxxzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 149); 

                auto tg_zz_xxyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 150); 

                auto tg_zz_xxyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 151); 

                auto tg_zz_xxyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 152); 

                auto tg_zz_xxyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 153); 

                auto tg_zz_xxzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 154); 

                auto tg_zz_xyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 155); 

                auto tg_zz_xyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 156); 

                auto tg_zz_xyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 157); 

                auto tg_zz_xyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 158); 

                auto tg_zz_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 159); 

                auto tg_zz_xzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 160); 

                // set up pointers to integrals

                auto tg_xyy_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 70); 

                auto tg_xyy_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 71); 

                auto tg_xyy_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 72); 

                auto tg_xyy_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 73); 

                auto tg_xyy_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 74); 

                auto tg_xyy_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 75); 

                auto tg_xyy_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 76); 

                auto tg_xyy_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 77); 

                auto tg_xyy_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 78); 

                auto tg_xyy_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 79); 

                auto tg_xyy_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 80); 

                auto tg_xyy_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 81); 

                auto tg_xyy_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 82); 

                auto tg_xyy_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 83); 

                auto tg_xyz_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 84); 

                auto tg_xyz_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 85); 

                auto tg_xyz_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 86); 

                auto tg_xyz_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 87); 

                auto tg_xyz_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 88); 

                auto tg_xyz_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 89); 

                auto tg_xyz_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 90); 

                auto tg_xyz_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 91); 

                auto tg_xyz_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 92); 

                auto tg_xyz_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 93); 

                auto tg_xyz_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 94); 

                auto tg_xyz_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 95); 

                auto tg_xyz_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 96); 

                auto tg_xyz_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 97); 

                auto tg_xyz_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 98); 

                auto tg_xyz_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 99); 

                auto tg_xyz_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 100); 

                auto tg_xyz_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 101); 

                auto tg_xyz_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 102); 

                auto tg_xyz_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 103); 

                auto tg_xyz_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 104); 

                auto tg_xzz_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 105); 

                auto tg_xzz_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 106); 

                auto tg_xzz_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 107); 

                auto tg_xzz_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 108); 

                auto tg_xzz_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 109); 

                auto tg_xzz_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 110); 

                auto tg_xzz_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 111); 

                auto tg_xzz_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 112); 

                auto tg_xzz_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 113); 

                auto tg_xzz_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 114); 

                auto tg_xzz_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 115); 

                auto tg_xzz_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 116); 

                auto tg_xzz_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 117); 

                auto tg_xzz_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 118); 

                auto tg_xzz_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 119); 

                auto tg_xzz_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 120); 

                auto tg_xzz_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 121); 

                auto tg_xzz_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 122); 

                auto tg_xzz_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 123); 

                auto tg_xzz_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 124); 

                auto tg_xzz_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 125); 

                auto tg_yyy_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 126); 

                auto tg_yyy_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 127); 

                auto tg_yyy_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 128); 

                auto tg_yyy_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 129); 

                auto tg_yyy_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 130); 

                auto tg_yyy_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 131); 

                auto tg_yyy_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 132); 

                auto tg_yyy_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 133); 

                auto tg_yyy_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 134); 

                auto tg_yyy_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 135); 

                auto tg_yyy_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 136); 

                auto tg_yyy_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 137); 

                auto tg_yyy_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 138); 

                auto tg_yyy_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 139); 

                // Batch of Integrals (70,140)

                #pragma omp simd aligned(cd_x, cd_y, tg_xyy_xxyyz_0, tg_xyy_xxyzz_0, tg_xyy_xxzzz_0, \
                                         tg_xyy_xyyyy_0, tg_xyy_xyyyz_0, tg_xyy_xyyzz_0, tg_xyy_xyzzz_0, tg_xyy_xzzzz_0, \
                                         tg_xyy_yyyyy_0, tg_xyy_yyyyz_0, tg_xyy_yyyzz_0, tg_xyy_yyzzz_0, tg_xyy_yzzzz_0, \
                                         tg_xyy_zzzzz_0, tg_xyz_xxxxx_0, tg_xyz_xxxxy_0, tg_xyz_xxxxz_0, tg_xyz_xxxyy_0, \
                                         tg_xyz_xxxyz_0, tg_xyz_xxxzz_0, tg_xyz_xxyyy_0, tg_xyz_xxyyz_0, tg_xyz_xxyzz_0, \
                                         tg_xyz_xxzzz_0, tg_xyz_xyyyy_0, tg_xyz_xyyyz_0, tg_xyz_xyyzz_0, tg_xyz_xyzzz_0, \
                                         tg_xyz_xzzzz_0, tg_xyz_yyyyy_0, tg_xyz_yyyyz_0, tg_xyz_yyyzz_0, tg_xyz_yyzzz_0, \
                                         tg_xyz_yzzzz_0, tg_xyz_zzzzz_0, tg_xzz_xxxxx_0, tg_xzz_xxxxy_0, tg_xzz_xxxxz_0, \
                                         tg_xzz_xxxyy_0, tg_xzz_xxxyz_0, tg_xzz_xxxzz_0, tg_xzz_xxyyy_0, tg_xzz_xxyyz_0, \
                                         tg_xzz_xxyzz_0, tg_xzz_xxzzz_0, tg_xzz_xyyyy_0, tg_xzz_xyyyz_0, tg_xzz_xyyzz_0, \
                                         tg_xzz_xyzzz_0, tg_xzz_xzzzz_0, tg_xzz_yyyyy_0, tg_xzz_yyyyz_0, tg_xzz_yyyzz_0, \
                                         tg_xzz_yyzzz_0, tg_xzz_yzzzz_0, tg_xzz_zzzzz_0, tg_yy_xxxxx_0, tg_yy_xxxxxy_0, \
                                         tg_yy_xxxxy_0, tg_yy_xxxxyy_0, tg_yy_xxxxyz_0, tg_yy_xxxxz_0, tg_yy_xxxyy_0, \
                                         tg_yy_xxxyyy_0, tg_yy_xxxyyz_0, tg_yy_xxxyz_0, tg_yy_xxxyzz_0, tg_yy_xxxzz_0, \
                                         tg_yy_xxxzzz_0, tg_yy_xxyyy_0, tg_yy_xxyyyy_0, tg_yy_xxyyyz_0, tg_yy_xxyyz_0, \
                                         tg_yy_xxyyzz_0, tg_yy_xxyzz_0, tg_yy_xxyzzz_0, tg_yy_xxzzz_0, tg_yy_xxzzzz_0, \
                                         tg_yy_xyyyy_0, tg_yy_xyyyyy_0, tg_yy_xyyyyz_0, tg_yy_xyyyz_0, tg_yy_xyyyzz_0, \
                                         tg_yy_xyyzz_0, tg_yy_xyyzzz_0, tg_yy_xyzzz_0, tg_yy_xyzzzz_0, tg_yy_xzzzz_0, \
                                         tg_yy_xzzzzz_0, tg_yy_yyyyy_0, tg_yy_yyyyz_0, tg_yy_yyyzz_0, tg_yy_yyzzz_0, \
                                         tg_yy_yzzzz_0, tg_yy_zzzzz_0, tg_yyy_xxxxx_0, tg_yyy_xxxxy_0, tg_yyy_xxxxz_0, \
                                         tg_yyy_xxxyy_0, tg_yyy_xxxyz_0, tg_yyy_xxxzz_0, tg_yyy_xxyyy_0, tg_yyy_xxyyz_0, \
                                         tg_yyy_xxyzz_0, tg_yyy_xxzzz_0, tg_yyy_xyyyy_0, tg_yyy_xyyyz_0, tg_yyy_xyyzz_0, \
                                         tg_yyy_xyzzz_0, tg_yz_xxxxx_0, tg_yz_xxxxxx_0, tg_yz_xxxxxy_0, tg_yz_xxxxxz_0, \
                                         tg_yz_xxxxy_0, tg_yz_xxxxyy_0, tg_yz_xxxxyz_0, tg_yz_xxxxz_0, tg_yz_xxxxzz_0, \
                                         tg_yz_xxxyy_0, tg_yz_xxxyyy_0, tg_yz_xxxyyz_0, tg_yz_xxxyz_0, tg_yz_xxxyzz_0, \
                                         tg_yz_xxxzz_0, tg_yz_xxxzzz_0, tg_yz_xxyyy_0, tg_yz_xxyyyy_0, tg_yz_xxyyyz_0, \
                                         tg_yz_xxyyz_0, tg_yz_xxyyzz_0, tg_yz_xxyzz_0, tg_yz_xxyzzz_0, tg_yz_xxzzz_0, \
                                         tg_yz_xxzzzz_0, tg_yz_xyyyy_0, tg_yz_xyyyyy_0, tg_yz_xyyyyz_0, tg_yz_xyyyz_0, \
                                         tg_yz_xyyyzz_0, tg_yz_xyyzz_0, tg_yz_xyyzzz_0, tg_yz_xyzzz_0, tg_yz_xyzzzz_0, \
                                         tg_yz_xzzzz_0, tg_yz_xzzzzz_0, tg_yz_yyyyy_0, tg_yz_yyyyz_0, tg_yz_yyyzz_0, \
                                         tg_yz_yyzzz_0, tg_yz_yzzzz_0, tg_yz_zzzzz_0, tg_zz_xxxxx_0, tg_zz_xxxxxx_0, \
                                         tg_zz_xxxxxy_0, tg_zz_xxxxxz_0, tg_zz_xxxxy_0, tg_zz_xxxxyy_0, tg_zz_xxxxyz_0, \
                                         tg_zz_xxxxz_0, tg_zz_xxxxzz_0, tg_zz_xxxyy_0, tg_zz_xxxyyy_0, tg_zz_xxxyyz_0, \
                                         tg_zz_xxxyz_0, tg_zz_xxxyzz_0, tg_zz_xxxzz_0, tg_zz_xxxzzz_0, tg_zz_xxyyy_0, \
                                         tg_zz_xxyyyy_0, tg_zz_xxyyyz_0, tg_zz_xxyyz_0, tg_zz_xxyyzz_0, tg_zz_xxyzz_0, \
                                         tg_zz_xxyzzz_0, tg_zz_xxzzz_0, tg_zz_xxzzzz_0, tg_zz_xyyyy_0, tg_zz_xyyyyy_0, \
                                         tg_zz_xyyyyz_0, tg_zz_xyyyz_0, tg_zz_xyyyzz_0, tg_zz_xyyzz_0, tg_zz_xyyzzz_0, \
                                         tg_zz_xyzzz_0, tg_zz_xyzzzz_0, tg_zz_xzzzz_0, tg_zz_xzzzzz_0, tg_zz_yyyyy_0, \
                                         tg_zz_yyyyz_0, tg_zz_yyyzz_0, tg_zz_yyzzz_0, tg_zz_yzzzz_0, tg_zz_zzzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nKetContrPairs; j++)
                {
                    tg_xyy_xxyyz_0[j] = -cd_x[j] * tg_yy_xxyyz_0[j] + tg_yy_xxxyyz_0[j];

                    tg_xyy_xxyzz_0[j] = -cd_x[j] * tg_yy_xxyzz_0[j] + tg_yy_xxxyzz_0[j];

                    tg_xyy_xxzzz_0[j] = -cd_x[j] * tg_yy_xxzzz_0[j] + tg_yy_xxxzzz_0[j];

                    tg_xyy_xyyyy_0[j] = -cd_x[j] * tg_yy_xyyyy_0[j] + tg_yy_xxyyyy_0[j];

                    tg_xyy_xyyyz_0[j] = -cd_x[j] * tg_yy_xyyyz_0[j] + tg_yy_xxyyyz_0[j];

                    tg_xyy_xyyzz_0[j] = -cd_x[j] * tg_yy_xyyzz_0[j] + tg_yy_xxyyzz_0[j];

                    tg_xyy_xyzzz_0[j] = -cd_x[j] * tg_yy_xyzzz_0[j] + tg_yy_xxyzzz_0[j];

                    tg_xyy_xzzzz_0[j] = -cd_x[j] * tg_yy_xzzzz_0[j] + tg_yy_xxzzzz_0[j];

                    tg_xyy_yyyyy_0[j] = -cd_x[j] * tg_yy_yyyyy_0[j] + tg_yy_xyyyyy_0[j];

                    tg_xyy_yyyyz_0[j] = -cd_x[j] * tg_yy_yyyyz_0[j] + tg_yy_xyyyyz_0[j];

                    tg_xyy_yyyzz_0[j] = -cd_x[j] * tg_yy_yyyzz_0[j] + tg_yy_xyyyzz_0[j];

                    tg_xyy_yyzzz_0[j] = -cd_x[j] * tg_yy_yyzzz_0[j] + tg_yy_xyyzzz_0[j];

                    tg_xyy_yzzzz_0[j] = -cd_x[j] * tg_yy_yzzzz_0[j] + tg_yy_xyzzzz_0[j];

                    tg_xyy_zzzzz_0[j] = -cd_x[j] * tg_yy_zzzzz_0[j] + tg_yy_xzzzzz_0[j];

                    tg_xyz_xxxxx_0[j] = -cd_x[j] * tg_yz_xxxxx_0[j] + tg_yz_xxxxxx_0[j];

                    tg_xyz_xxxxy_0[j] = -cd_x[j] * tg_yz_xxxxy_0[j] + tg_yz_xxxxxy_0[j];

                    tg_xyz_xxxxz_0[j] = -cd_x[j] * tg_yz_xxxxz_0[j] + tg_yz_xxxxxz_0[j];

                    tg_xyz_xxxyy_0[j] = -cd_x[j] * tg_yz_xxxyy_0[j] + tg_yz_xxxxyy_0[j];

                    tg_xyz_xxxyz_0[j] = -cd_x[j] * tg_yz_xxxyz_0[j] + tg_yz_xxxxyz_0[j];

                    tg_xyz_xxxzz_0[j] = -cd_x[j] * tg_yz_xxxzz_0[j] + tg_yz_xxxxzz_0[j];

                    tg_xyz_xxyyy_0[j] = -cd_x[j] * tg_yz_xxyyy_0[j] + tg_yz_xxxyyy_0[j];

                    tg_xyz_xxyyz_0[j] = -cd_x[j] * tg_yz_xxyyz_0[j] + tg_yz_xxxyyz_0[j];

                    tg_xyz_xxyzz_0[j] = -cd_x[j] * tg_yz_xxyzz_0[j] + tg_yz_xxxyzz_0[j];

                    tg_xyz_xxzzz_0[j] = -cd_x[j] * tg_yz_xxzzz_0[j] + tg_yz_xxxzzz_0[j];

                    tg_xyz_xyyyy_0[j] = -cd_x[j] * tg_yz_xyyyy_0[j] + tg_yz_xxyyyy_0[j];

                    tg_xyz_xyyyz_0[j] = -cd_x[j] * tg_yz_xyyyz_0[j] + tg_yz_xxyyyz_0[j];

                    tg_xyz_xyyzz_0[j] = -cd_x[j] * tg_yz_xyyzz_0[j] + tg_yz_xxyyzz_0[j];

                    tg_xyz_xyzzz_0[j] = -cd_x[j] * tg_yz_xyzzz_0[j] + tg_yz_xxyzzz_0[j];

                    tg_xyz_xzzzz_0[j] = -cd_x[j] * tg_yz_xzzzz_0[j] + tg_yz_xxzzzz_0[j];

                    tg_xyz_yyyyy_0[j] = -cd_x[j] * tg_yz_yyyyy_0[j] + tg_yz_xyyyyy_0[j];

                    tg_xyz_yyyyz_0[j] = -cd_x[j] * tg_yz_yyyyz_0[j] + tg_yz_xyyyyz_0[j];

                    tg_xyz_yyyzz_0[j] = -cd_x[j] * tg_yz_yyyzz_0[j] + tg_yz_xyyyzz_0[j];

                    tg_xyz_yyzzz_0[j] = -cd_x[j] * tg_yz_yyzzz_0[j] + tg_yz_xyyzzz_0[j];

                    tg_xyz_yzzzz_0[j] = -cd_x[j] * tg_yz_yzzzz_0[j] + tg_yz_xyzzzz_0[j];

                    tg_xyz_zzzzz_0[j] = -cd_x[j] * tg_yz_zzzzz_0[j] + tg_yz_xzzzzz_0[j];

                    tg_xzz_xxxxx_0[j] = -cd_x[j] * tg_zz_xxxxx_0[j] + tg_zz_xxxxxx_0[j];

                    tg_xzz_xxxxy_0[j] = -cd_x[j] * tg_zz_xxxxy_0[j] + tg_zz_xxxxxy_0[j];

                    tg_xzz_xxxxz_0[j] = -cd_x[j] * tg_zz_xxxxz_0[j] + tg_zz_xxxxxz_0[j];

                    tg_xzz_xxxyy_0[j] = -cd_x[j] * tg_zz_xxxyy_0[j] + tg_zz_xxxxyy_0[j];

                    tg_xzz_xxxyz_0[j] = -cd_x[j] * tg_zz_xxxyz_0[j] + tg_zz_xxxxyz_0[j];

                    tg_xzz_xxxzz_0[j] = -cd_x[j] * tg_zz_xxxzz_0[j] + tg_zz_xxxxzz_0[j];

                    tg_xzz_xxyyy_0[j] = -cd_x[j] * tg_zz_xxyyy_0[j] + tg_zz_xxxyyy_0[j];

                    tg_xzz_xxyyz_0[j] = -cd_x[j] * tg_zz_xxyyz_0[j] + tg_zz_xxxyyz_0[j];

                    tg_xzz_xxyzz_0[j] = -cd_x[j] * tg_zz_xxyzz_0[j] + tg_zz_xxxyzz_0[j];

                    tg_xzz_xxzzz_0[j] = -cd_x[j] * tg_zz_xxzzz_0[j] + tg_zz_xxxzzz_0[j];

                    tg_xzz_xyyyy_0[j] = -cd_x[j] * tg_zz_xyyyy_0[j] + tg_zz_xxyyyy_0[j];

                    tg_xzz_xyyyz_0[j] = -cd_x[j] * tg_zz_xyyyz_0[j] + tg_zz_xxyyyz_0[j];

                    tg_xzz_xyyzz_0[j] = -cd_x[j] * tg_zz_xyyzz_0[j] + tg_zz_xxyyzz_0[j];

                    tg_xzz_xyzzz_0[j] = -cd_x[j] * tg_zz_xyzzz_0[j] + tg_zz_xxyzzz_0[j];

                    tg_xzz_xzzzz_0[j] = -cd_x[j] * tg_zz_xzzzz_0[j] + tg_zz_xxzzzz_0[j];

                    tg_xzz_yyyyy_0[j] = -cd_x[j] * tg_zz_yyyyy_0[j] + tg_zz_xyyyyy_0[j];

                    tg_xzz_yyyyz_0[j] = -cd_x[j] * tg_zz_yyyyz_0[j] + tg_zz_xyyyyz_0[j];

                    tg_xzz_yyyzz_0[j] = -cd_x[j] * tg_zz_yyyzz_0[j] + tg_zz_xyyyzz_0[j];

                    tg_xzz_yyzzz_0[j] = -cd_x[j] * tg_zz_yyzzz_0[j] + tg_zz_xyyzzz_0[j];

                    tg_xzz_yzzzz_0[j] = -cd_x[j] * tg_zz_yzzzz_0[j] + tg_zz_xyzzzz_0[j];

                    tg_xzz_zzzzz_0[j] = -cd_x[j] * tg_zz_zzzzz_0[j] + tg_zz_xzzzzz_0[j];

                    tg_yyy_xxxxx_0[j] = -cd_y[j] * tg_yy_xxxxx_0[j] + tg_yy_xxxxxy_0[j];

                    tg_yyy_xxxxy_0[j] = -cd_y[j] * tg_yy_xxxxy_0[j] + tg_yy_xxxxyy_0[j];

                    tg_yyy_xxxxz_0[j] = -cd_y[j] * tg_yy_xxxxz_0[j] + tg_yy_xxxxyz_0[j];

                    tg_yyy_xxxyy_0[j] = -cd_y[j] * tg_yy_xxxyy_0[j] + tg_yy_xxxyyy_0[j];

                    tg_yyy_xxxyz_0[j] = -cd_y[j] * tg_yy_xxxyz_0[j] + tg_yy_xxxyyz_0[j];

                    tg_yyy_xxxzz_0[j] = -cd_y[j] * tg_yy_xxxzz_0[j] + tg_yy_xxxyzz_0[j];

                    tg_yyy_xxyyy_0[j] = -cd_y[j] * tg_yy_xxyyy_0[j] + tg_yy_xxyyyy_0[j];

                    tg_yyy_xxyyz_0[j] = -cd_y[j] * tg_yy_xxyyz_0[j] + tg_yy_xxyyyz_0[j];

                    tg_yyy_xxyzz_0[j] = -cd_y[j] * tg_yy_xxyzz_0[j] + tg_yy_xxyyzz_0[j];

                    tg_yyy_xxzzz_0[j] = -cd_y[j] * tg_yy_xxzzz_0[j] + tg_yy_xxyzzz_0[j];

                    tg_yyy_xyyyy_0[j] = -cd_y[j] * tg_yy_xyyyy_0[j] + tg_yy_xyyyyy_0[j];

                    tg_yyy_xyyyz_0[j] = -cd_y[j] * tg_yy_xyyyz_0[j] + tg_yy_xyyyyz_0[j];

                    tg_yyy_xyyzz_0[j] = -cd_y[j] * tg_yy_xyyzz_0[j] + tg_yy_xyyyzz_0[j];

                    tg_yyy_xyzzz_0[j] = -cd_y[j] * tg_yy_xyzzz_0[j] + tg_yy_xyyzzz_0[j];
                }
            }
        }
    }

    void
    compElectronRepulsionForSXFH_140_210(      CMemBlock2D<double>& ketBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& cdDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetContrPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (140,210)

        // set up pointers to distances R(CD) = C - D

        auto cd_y = cdDistances.data(1);

        auto cd_z = cdDistances.data(2);

        // set up bra side loop over intergals

        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();

        for (int32_t iang = 0; iang <= bang; iang++)
        {
            // set up index of integral

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 5, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 5, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_2_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {2, 6, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

                auto tg_yy_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 77); 

                auto tg_yy_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 78); 

                auto tg_yy_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 79); 

                auto tg_yy_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 80); 

                auto tg_yy_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 81); 

                auto tg_yy_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 82); 

                auto tg_yy_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 83); 

                auto tg_yz_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 84); 

                auto tg_yz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 85); 

                auto tg_yz_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 86); 

                auto tg_yz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 87); 

                auto tg_yz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 88); 

                auto tg_yz_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 89); 

                auto tg_yz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 90); 

                auto tg_yz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 91); 

                auto tg_yz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 92); 

                auto tg_yz_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 93); 

                auto tg_yz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 94); 

                auto tg_yz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 95); 

                auto tg_yz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 96); 

                auto tg_yz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 97); 

                auto tg_yz_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 98); 

                auto tg_yz_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 99); 

                auto tg_yz_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 100); 

                auto tg_yz_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 101); 

                auto tg_yz_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 102); 

                auto tg_yz_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 103); 

                auto tg_yz_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 104); 

                auto tg_zz_xxxxx_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 105); 

                auto tg_zz_xxxxy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 106); 

                auto tg_zz_xxxxz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 107); 

                auto tg_zz_xxxyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 108); 

                auto tg_zz_xxxyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 109); 

                auto tg_zz_xxxzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 110); 

                auto tg_zz_xxyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 111); 

                auto tg_zz_xxyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 112); 

                auto tg_zz_xxyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 113); 

                auto tg_zz_xxzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 114); 

                auto tg_zz_xyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 115); 

                auto tg_zz_xyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 116); 

                auto tg_zz_xyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 117); 

                auto tg_zz_xyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 118); 

                auto tg_zz_xzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 119); 

                auto tg_zz_yyyyy_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 120); 

                auto tg_zz_yyyyz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 121); 

                auto tg_zz_yyyzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 122); 

                auto tg_zz_yyzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 123); 

                auto tg_zz_yzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 124); 

                auto tg_zz_zzzzz_0 = ketBuffer.data(pidx_g_2_5_m0 + 126 * i + 125); 

                auto tg_yy_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 103); 

                auto tg_yy_yyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 105); 

                auto tg_yy_yyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 106); 

                auto tg_yy_yyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 107); 

                auto tg_yy_yyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 108); 

                auto tg_yy_yyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 109); 

                auto tg_yy_yzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 110); 

                auto tg_yz_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 113); 

                auto tg_yz_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 115); 

                auto tg_yz_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 116); 

                auto tg_yz_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 118); 

                auto tg_yz_xxxyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 119); 

                auto tg_yz_xxxyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 120); 

                auto tg_yz_xxyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 122); 

                auto tg_yz_xxyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 123); 

                auto tg_yz_xxyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 124); 

                auto tg_yz_xxyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 125); 

                auto tg_yz_xyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 127); 

                auto tg_yz_xyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 128); 

                auto tg_yz_xyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 129); 

                auto tg_yz_xyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 130); 

                auto tg_yz_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 131); 

                auto tg_yz_yyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 133); 

                auto tg_yz_yyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 134); 

                auto tg_yz_yyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 135); 

                auto tg_yz_yyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 136); 

                auto tg_yz_yyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 137); 

                auto tg_yz_yzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 138); 

                auto tg_zz_xxxxxy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 141); 

                auto tg_zz_xxxxxz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 142); 

                auto tg_zz_xxxxyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 143); 

                auto tg_zz_xxxxyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 144); 

                auto tg_zz_xxxxzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 145); 

                auto tg_zz_xxxyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 146); 

                auto tg_zz_xxxyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 147); 

                auto tg_zz_xxxyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 148); 

                auto tg_zz_xxxzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 149); 

                auto tg_zz_xxyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 150); 

                auto tg_zz_xxyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 151); 

                auto tg_zz_xxyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 152); 

                auto tg_zz_xxyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 153); 

                auto tg_zz_xxzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 154); 

                auto tg_zz_xyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 155); 

                auto tg_zz_xyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 156); 

                auto tg_zz_xyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 157); 

                auto tg_zz_xyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 158); 

                auto tg_zz_xyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 159); 

                auto tg_zz_xzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 160); 

                auto tg_zz_yyyyyy_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 161); 

                auto tg_zz_yyyyyz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 162); 

                auto tg_zz_yyyyzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 163); 

                auto tg_zz_yyyzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 164); 

                auto tg_zz_yyzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 165); 

                auto tg_zz_yzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 166); 

                auto tg_zz_zzzzzz_0 = ketBuffer.data(pidx_g_2_6_m0 + 168 * i + 167); 

                // set up pointers to integrals

                auto tg_yyy_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 140); 

                auto tg_yyy_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 141); 

                auto tg_yyy_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 142); 

                auto tg_yyy_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 143); 

                auto tg_yyy_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 144); 

                auto tg_yyy_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 145); 

                auto tg_yyy_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 146); 

                auto tg_yyz_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 147); 

                auto tg_yyz_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 148); 

                auto tg_yyz_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 149); 

                auto tg_yyz_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 150); 

                auto tg_yyz_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 151); 

                auto tg_yyz_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 152); 

                auto tg_yyz_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 153); 

                auto tg_yyz_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 154); 

                auto tg_yyz_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 155); 

                auto tg_yyz_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 156); 

                auto tg_yyz_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 157); 

                auto tg_yyz_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 158); 

                auto tg_yyz_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 159); 

                auto tg_yyz_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 160); 

                auto tg_yyz_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 161); 

                auto tg_yyz_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 162); 

                auto tg_yyz_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 163); 

                auto tg_yyz_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 164); 

                auto tg_yyz_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 165); 

                auto tg_yyz_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 166); 

                auto tg_yyz_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 167); 

                auto tg_yzz_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 168); 

                auto tg_yzz_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 169); 

                auto tg_yzz_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 170); 

                auto tg_yzz_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 171); 

                auto tg_yzz_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 172); 

                auto tg_yzz_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 173); 

                auto tg_yzz_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 174); 

                auto tg_yzz_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 175); 

                auto tg_yzz_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 176); 

                auto tg_yzz_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 177); 

                auto tg_yzz_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 178); 

                auto tg_yzz_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 179); 

                auto tg_yzz_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 180); 

                auto tg_yzz_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 181); 

                auto tg_yzz_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 182); 

                auto tg_yzz_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 183); 

                auto tg_yzz_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 184); 

                auto tg_yzz_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 185); 

                auto tg_yzz_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 186); 

                auto tg_yzz_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 187); 

                auto tg_yzz_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 188); 

                auto tg_zzz_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 189); 

                auto tg_zzz_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 190); 

                auto tg_zzz_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 191); 

                auto tg_zzz_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 192); 

                auto tg_zzz_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 193); 

                auto tg_zzz_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 194); 

                auto tg_zzz_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 195); 

                auto tg_zzz_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 196); 

                auto tg_zzz_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 197); 

                auto tg_zzz_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 198); 

                auto tg_zzz_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 199); 

                auto tg_zzz_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 200); 

                auto tg_zzz_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 201); 

                auto tg_zzz_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 202); 

                auto tg_zzz_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 203); 

                auto tg_zzz_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 204); 

                auto tg_zzz_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 205); 

                auto tg_zzz_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 206); 

                auto tg_zzz_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 207); 

                auto tg_zzz_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 208); 

                auto tg_zzz_zzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 209); 

                // Batch of Integrals (140,210)

                #pragma omp simd aligned(cd_y, cd_z, tg_yy_xyzzzz_0, tg_yy_xzzzz_0, tg_yy_yyyyy_0, \
                                         tg_yy_yyyyyy_0, tg_yy_yyyyyz_0, tg_yy_yyyyz_0, tg_yy_yyyyzz_0, tg_yy_yyyzz_0, \
                                         tg_yy_yyyzzz_0, tg_yy_yyzzz_0, tg_yy_yyzzzz_0, tg_yy_yzzzz_0, tg_yy_yzzzzz_0, \
                                         tg_yy_zzzzz_0, tg_yyy_xzzzz_0, tg_yyy_yyyyy_0, tg_yyy_yyyyz_0, tg_yyy_yyyzz_0, \
                                         tg_yyy_yyzzz_0, tg_yyy_yzzzz_0, tg_yyy_zzzzz_0, tg_yyz_xxxxx_0, tg_yyz_xxxxy_0, \
                                         tg_yyz_xxxxz_0, tg_yyz_xxxyy_0, tg_yyz_xxxyz_0, tg_yyz_xxxzz_0, tg_yyz_xxyyy_0, \
                                         tg_yyz_xxyyz_0, tg_yyz_xxyzz_0, tg_yyz_xxzzz_0, tg_yyz_xyyyy_0, tg_yyz_xyyyz_0, \
                                         tg_yyz_xyyzz_0, tg_yyz_xyzzz_0, tg_yyz_xzzzz_0, tg_yyz_yyyyy_0, tg_yyz_yyyyz_0, \
                                         tg_yyz_yyyzz_0, tg_yyz_yyzzz_0, tg_yyz_yzzzz_0, tg_yyz_zzzzz_0, tg_yz_xxxxx_0, \
                                         tg_yz_xxxxxy_0, tg_yz_xxxxy_0, tg_yz_xxxxyy_0, tg_yz_xxxxyz_0, tg_yz_xxxxz_0, \
                                         tg_yz_xxxyy_0, tg_yz_xxxyyy_0, tg_yz_xxxyyz_0, tg_yz_xxxyz_0, tg_yz_xxxyzz_0, \
                                         tg_yz_xxxzz_0, tg_yz_xxyyy_0, tg_yz_xxyyyy_0, tg_yz_xxyyyz_0, tg_yz_xxyyz_0, \
                                         tg_yz_xxyyzz_0, tg_yz_xxyzz_0, tg_yz_xxyzzz_0, tg_yz_xxzzz_0, tg_yz_xyyyy_0, \
                                         tg_yz_xyyyyy_0, tg_yz_xyyyyz_0, tg_yz_xyyyz_0, tg_yz_xyyyzz_0, tg_yz_xyyzz_0, \
                                         tg_yz_xyyzzz_0, tg_yz_xyzzz_0, tg_yz_xyzzzz_0, tg_yz_xzzzz_0, tg_yz_yyyyy_0, \
                                         tg_yz_yyyyyy_0, tg_yz_yyyyyz_0, tg_yz_yyyyz_0, tg_yz_yyyyzz_0, tg_yz_yyyzz_0, \
                                         tg_yz_yyyzzz_0, tg_yz_yyzzz_0, tg_yz_yyzzzz_0, tg_yz_yzzzz_0, tg_yz_yzzzzz_0, \
                                         tg_yz_zzzzz_0, tg_yzz_xxxxx_0, tg_yzz_xxxxy_0, tg_yzz_xxxxz_0, tg_yzz_xxxyy_0, \
                                         tg_yzz_xxxyz_0, tg_yzz_xxxzz_0, tg_yzz_xxyyy_0, tg_yzz_xxyyz_0, tg_yzz_xxyzz_0, \
                                         tg_yzz_xxzzz_0, tg_yzz_xyyyy_0, tg_yzz_xyyyz_0, tg_yzz_xyyzz_0, tg_yzz_xyzzz_0, \
                                         tg_yzz_xzzzz_0, tg_yzz_yyyyy_0, tg_yzz_yyyyz_0, tg_yzz_yyyzz_0, tg_yzz_yyzzz_0, \
                                         tg_yzz_yzzzz_0, tg_yzz_zzzzz_0, tg_zz_xxxxx_0, tg_zz_xxxxxy_0, tg_zz_xxxxxz_0, \
                                         tg_zz_xxxxy_0, tg_zz_xxxxyy_0, tg_zz_xxxxyz_0, tg_zz_xxxxz_0, tg_zz_xxxxzz_0, \
                                         tg_zz_xxxyy_0, tg_zz_xxxyyy_0, tg_zz_xxxyyz_0, tg_zz_xxxyz_0, tg_zz_xxxyzz_0, \
                                         tg_zz_xxxzz_0, tg_zz_xxxzzz_0, tg_zz_xxyyy_0, tg_zz_xxyyyy_0, tg_zz_xxyyyz_0, \
                                         tg_zz_xxyyz_0, tg_zz_xxyyzz_0, tg_zz_xxyzz_0, tg_zz_xxyzzz_0, tg_zz_xxzzz_0, \
                                         tg_zz_xxzzzz_0, tg_zz_xyyyy_0, tg_zz_xyyyyy_0, tg_zz_xyyyyz_0, tg_zz_xyyyz_0, \
                                         tg_zz_xyyyzz_0, tg_zz_xyyzz_0, tg_zz_xyyzzz_0, tg_zz_xyzzz_0, tg_zz_xyzzzz_0, \
                                         tg_zz_xzzzz_0, tg_zz_xzzzzz_0, tg_zz_yyyyy_0, tg_zz_yyyyyy_0, tg_zz_yyyyyz_0, \
                                         tg_zz_yyyyz_0, tg_zz_yyyyzz_0, tg_zz_yyyzz_0, tg_zz_yyyzzz_0, tg_zz_yyzzz_0, \
                                         tg_zz_yyzzzz_0, tg_zz_yzzzz_0, tg_zz_yzzzzz_0, tg_zz_zzzzz_0, tg_zz_zzzzzz_0, \
                                         tg_zzz_xxxxx_0, tg_zzz_xxxxy_0, tg_zzz_xxxxz_0, tg_zzz_xxxyy_0, tg_zzz_xxxyz_0, \
                                         tg_zzz_xxxzz_0, tg_zzz_xxyyy_0, tg_zzz_xxyyz_0, tg_zzz_xxyzz_0, tg_zzz_xxzzz_0, \
                                         tg_zzz_xyyyy_0, tg_zzz_xyyyz_0, tg_zzz_xyyzz_0, tg_zzz_xyzzz_0, tg_zzz_xzzzz_0, \
                                         tg_zzz_yyyyy_0, tg_zzz_yyyyz_0, tg_zzz_yyyzz_0, tg_zzz_yyzzz_0, tg_zzz_yzzzz_0, \
                                         tg_zzz_zzzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nKetContrPairs; j++)
                {
                    tg_yyy_xzzzz_0[j] = -cd_y[j] * tg_yy_xzzzz_0[j] + tg_yy_xyzzzz_0[j];

                    tg_yyy_yyyyy_0[j] = -cd_y[j] * tg_yy_yyyyy_0[j] + tg_yy_yyyyyy_0[j];

                    tg_yyy_yyyyz_0[j] = -cd_y[j] * tg_yy_yyyyz_0[j] + tg_yy_yyyyyz_0[j];

                    tg_yyy_yyyzz_0[j] = -cd_y[j] * tg_yy_yyyzz_0[j] + tg_yy_yyyyzz_0[j];

                    tg_yyy_yyzzz_0[j] = -cd_y[j] * tg_yy_yyzzz_0[j] + tg_yy_yyyzzz_0[j];

                    tg_yyy_yzzzz_0[j] = -cd_y[j] * tg_yy_yzzzz_0[j] + tg_yy_yyzzzz_0[j];

                    tg_yyy_zzzzz_0[j] = -cd_y[j] * tg_yy_zzzzz_0[j] + tg_yy_yzzzzz_0[j];

                    tg_yyz_xxxxx_0[j] = -cd_y[j] * tg_yz_xxxxx_0[j] + tg_yz_xxxxxy_0[j];

                    tg_yyz_xxxxy_0[j] = -cd_y[j] * tg_yz_xxxxy_0[j] + tg_yz_xxxxyy_0[j];

                    tg_yyz_xxxxz_0[j] = -cd_y[j] * tg_yz_xxxxz_0[j] + tg_yz_xxxxyz_0[j];

                    tg_yyz_xxxyy_0[j] = -cd_y[j] * tg_yz_xxxyy_0[j] + tg_yz_xxxyyy_0[j];

                    tg_yyz_xxxyz_0[j] = -cd_y[j] * tg_yz_xxxyz_0[j] + tg_yz_xxxyyz_0[j];

                    tg_yyz_xxxzz_0[j] = -cd_y[j] * tg_yz_xxxzz_0[j] + tg_yz_xxxyzz_0[j];

                    tg_yyz_xxyyy_0[j] = -cd_y[j] * tg_yz_xxyyy_0[j] + tg_yz_xxyyyy_0[j];

                    tg_yyz_xxyyz_0[j] = -cd_y[j] * tg_yz_xxyyz_0[j] + tg_yz_xxyyyz_0[j];

                    tg_yyz_xxyzz_0[j] = -cd_y[j] * tg_yz_xxyzz_0[j] + tg_yz_xxyyzz_0[j];

                    tg_yyz_xxzzz_0[j] = -cd_y[j] * tg_yz_xxzzz_0[j] + tg_yz_xxyzzz_0[j];

                    tg_yyz_xyyyy_0[j] = -cd_y[j] * tg_yz_xyyyy_0[j] + tg_yz_xyyyyy_0[j];

                    tg_yyz_xyyyz_0[j] = -cd_y[j] * tg_yz_xyyyz_0[j] + tg_yz_xyyyyz_0[j];

                    tg_yyz_xyyzz_0[j] = -cd_y[j] * tg_yz_xyyzz_0[j] + tg_yz_xyyyzz_0[j];

                    tg_yyz_xyzzz_0[j] = -cd_y[j] * tg_yz_xyzzz_0[j] + tg_yz_xyyzzz_0[j];

                    tg_yyz_xzzzz_0[j] = -cd_y[j] * tg_yz_xzzzz_0[j] + tg_yz_xyzzzz_0[j];

                    tg_yyz_yyyyy_0[j] = -cd_y[j] * tg_yz_yyyyy_0[j] + tg_yz_yyyyyy_0[j];

                    tg_yyz_yyyyz_0[j] = -cd_y[j] * tg_yz_yyyyz_0[j] + tg_yz_yyyyyz_0[j];

                    tg_yyz_yyyzz_0[j] = -cd_y[j] * tg_yz_yyyzz_0[j] + tg_yz_yyyyzz_0[j];

                    tg_yyz_yyzzz_0[j] = -cd_y[j] * tg_yz_yyzzz_0[j] + tg_yz_yyyzzz_0[j];

                    tg_yyz_yzzzz_0[j] = -cd_y[j] * tg_yz_yzzzz_0[j] + tg_yz_yyzzzz_0[j];

                    tg_yyz_zzzzz_0[j] = -cd_y[j] * tg_yz_zzzzz_0[j] + tg_yz_yzzzzz_0[j];

                    tg_yzz_xxxxx_0[j] = -cd_y[j] * tg_zz_xxxxx_0[j] + tg_zz_xxxxxy_0[j];

                    tg_yzz_xxxxy_0[j] = -cd_y[j] * tg_zz_xxxxy_0[j] + tg_zz_xxxxyy_0[j];

                    tg_yzz_xxxxz_0[j] = -cd_y[j] * tg_zz_xxxxz_0[j] + tg_zz_xxxxyz_0[j];

                    tg_yzz_xxxyy_0[j] = -cd_y[j] * tg_zz_xxxyy_0[j] + tg_zz_xxxyyy_0[j];

                    tg_yzz_xxxyz_0[j] = -cd_y[j] * tg_zz_xxxyz_0[j] + tg_zz_xxxyyz_0[j];

                    tg_yzz_xxxzz_0[j] = -cd_y[j] * tg_zz_xxxzz_0[j] + tg_zz_xxxyzz_0[j];

                    tg_yzz_xxyyy_0[j] = -cd_y[j] * tg_zz_xxyyy_0[j] + tg_zz_xxyyyy_0[j];

                    tg_yzz_xxyyz_0[j] = -cd_y[j] * tg_zz_xxyyz_0[j] + tg_zz_xxyyyz_0[j];

                    tg_yzz_xxyzz_0[j] = -cd_y[j] * tg_zz_xxyzz_0[j] + tg_zz_xxyyzz_0[j];

                    tg_yzz_xxzzz_0[j] = -cd_y[j] * tg_zz_xxzzz_0[j] + tg_zz_xxyzzz_0[j];

                    tg_yzz_xyyyy_0[j] = -cd_y[j] * tg_zz_xyyyy_0[j] + tg_zz_xyyyyy_0[j];

                    tg_yzz_xyyyz_0[j] = -cd_y[j] * tg_zz_xyyyz_0[j] + tg_zz_xyyyyz_0[j];

                    tg_yzz_xyyzz_0[j] = -cd_y[j] * tg_zz_xyyzz_0[j] + tg_zz_xyyyzz_0[j];

                    tg_yzz_xyzzz_0[j] = -cd_y[j] * tg_zz_xyzzz_0[j] + tg_zz_xyyzzz_0[j];

                    tg_yzz_xzzzz_0[j] = -cd_y[j] * tg_zz_xzzzz_0[j] + tg_zz_xyzzzz_0[j];

                    tg_yzz_yyyyy_0[j] = -cd_y[j] * tg_zz_yyyyy_0[j] + tg_zz_yyyyyy_0[j];

                    tg_yzz_yyyyz_0[j] = -cd_y[j] * tg_zz_yyyyz_0[j] + tg_zz_yyyyyz_0[j];

                    tg_yzz_yyyzz_0[j] = -cd_y[j] * tg_zz_yyyzz_0[j] + tg_zz_yyyyzz_0[j];

                    tg_yzz_yyzzz_0[j] = -cd_y[j] * tg_zz_yyzzz_0[j] + tg_zz_yyyzzz_0[j];

                    tg_yzz_yzzzz_0[j] = -cd_y[j] * tg_zz_yzzzz_0[j] + tg_zz_yyzzzz_0[j];

                    tg_yzz_zzzzz_0[j] = -cd_y[j] * tg_zz_zzzzz_0[j] + tg_zz_yzzzzz_0[j];

                    tg_zzz_xxxxx_0[j] = -cd_z[j] * tg_zz_xxxxx_0[j] + tg_zz_xxxxxz_0[j];

                    tg_zzz_xxxxy_0[j] = -cd_z[j] * tg_zz_xxxxy_0[j] + tg_zz_xxxxyz_0[j];

                    tg_zzz_xxxxz_0[j] = -cd_z[j] * tg_zz_xxxxz_0[j] + tg_zz_xxxxzz_0[j];

                    tg_zzz_xxxyy_0[j] = -cd_z[j] * tg_zz_xxxyy_0[j] + tg_zz_xxxyyz_0[j];

                    tg_zzz_xxxyz_0[j] = -cd_z[j] * tg_zz_xxxyz_0[j] + tg_zz_xxxyzz_0[j];

                    tg_zzz_xxxzz_0[j] = -cd_z[j] * tg_zz_xxxzz_0[j] + tg_zz_xxxzzz_0[j];

                    tg_zzz_xxyyy_0[j] = -cd_z[j] * tg_zz_xxyyy_0[j] + tg_zz_xxyyyz_0[j];

                    tg_zzz_xxyyz_0[j] = -cd_z[j] * tg_zz_xxyyz_0[j] + tg_zz_xxyyzz_0[j];

                    tg_zzz_xxyzz_0[j] = -cd_z[j] * tg_zz_xxyzz_0[j] + tg_zz_xxyzzz_0[j];

                    tg_zzz_xxzzz_0[j] = -cd_z[j] * tg_zz_xxzzz_0[j] + tg_zz_xxzzzz_0[j];

                    tg_zzz_xyyyy_0[j] = -cd_z[j] * tg_zz_xyyyy_0[j] + tg_zz_xyyyyz_0[j];

                    tg_zzz_xyyyz_0[j] = -cd_z[j] * tg_zz_xyyyz_0[j] + tg_zz_xyyyzz_0[j];

                    tg_zzz_xyyzz_0[j] = -cd_z[j] * tg_zz_xyyzz_0[j] + tg_zz_xyyzzz_0[j];

                    tg_zzz_xyzzz_0[j] = -cd_z[j] * tg_zz_xyzzz_0[j] + tg_zz_xyzzzz_0[j];

                    tg_zzz_xzzzz_0[j] = -cd_z[j] * tg_zz_xzzzz_0[j] + tg_zz_xzzzzz_0[j];

                    tg_zzz_yyyyy_0[j] = -cd_z[j] * tg_zz_yyyyy_0[j] + tg_zz_yyyyyz_0[j];

                    tg_zzz_yyyyz_0[j] = -cd_z[j] * tg_zz_yyyyz_0[j] + tg_zz_yyyyzz_0[j];

                    tg_zzz_yyyzz_0[j] = -cd_z[j] * tg_zz_yyyzz_0[j] + tg_zz_yyyzzz_0[j];

                    tg_zzz_yyzzz_0[j] = -cd_z[j] * tg_zz_yyzzz_0[j] + tg_zz_yyzzzz_0[j];

                    tg_zzz_yzzzz_0[j] = -cd_z[j] * tg_zz_yzzzz_0[j] + tg_zz_yzzzzz_0[j];

                    tg_zzz_zzzzz_0[j] = -cd_z[j] * tg_zz_zzzzz_0[j] + tg_zz_zzzzzz_0[j];
                }
            }
        }
    }


} // erikrrfunc namespace

