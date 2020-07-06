//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionBRRRecFuncForFXYY.hpp"

#include "AngularMomentum.hpp"

namespace eribrrfunc { // eribrrfunc namespace

    void
    compElectronRepulsionForFFXY(      CMemBlock2D<double>& braBuffer,
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

        auto pidx_g_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 3, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_3_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 3, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            // set up pointers to integrals

            auto tg_xxx_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + i); 

            auto tg_xxx_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + kcomp + i); 

            auto tg_xxx_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 2 * kcomp + i); 

            auto tg_xxx_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 3 * kcomp + i); 

            auto tg_xxx_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 4 * kcomp + i); 

            auto tg_xxx_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 5 * kcomp + i); 

            auto tg_xxx_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 6 * kcomp + i); 

            auto tg_xxx_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 7 * kcomp + i); 

            auto tg_xxx_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 8 * kcomp + i); 

            auto tg_xxx_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 9 * kcomp + i); 

            auto tg_xxy_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 10 * kcomp + i); 

            auto tg_xxy_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 11 * kcomp + i); 

            auto tg_xxy_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 12 * kcomp + i); 

            auto tg_xxy_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 13 * kcomp + i); 

            auto tg_xxy_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 14 * kcomp + i); 

            auto tg_xxy_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 15 * kcomp + i); 

            auto tg_xxy_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 16 * kcomp + i); 

            auto tg_xxy_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 17 * kcomp + i); 

            auto tg_xxy_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 18 * kcomp + i); 

            auto tg_xxy_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 19 * kcomp + i); 

            auto tg_xxz_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 20 * kcomp + i); 

            auto tg_xxz_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 21 * kcomp + i); 

            auto tg_xxz_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 22 * kcomp + i); 

            auto tg_xxz_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 23 * kcomp + i); 

            auto tg_xxz_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 24 * kcomp + i); 

            auto tg_xxz_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 25 * kcomp + i); 

            auto tg_xxz_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 26 * kcomp + i); 

            auto tg_xxz_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 27 * kcomp + i); 

            auto tg_xxz_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 28 * kcomp + i); 

            auto tg_xxz_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 29 * kcomp + i); 

            auto tg_xyy_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 30 * kcomp + i); 

            auto tg_xyy_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 31 * kcomp + i); 

            auto tg_xyy_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 32 * kcomp + i); 

            auto tg_xyy_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 33 * kcomp + i); 

            auto tg_xyy_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 34 * kcomp + i); 

            auto tg_xyy_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 35 * kcomp + i); 

            auto tg_xyy_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 36 * kcomp + i); 

            auto tg_xyy_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 37 * kcomp + i); 

            auto tg_xyy_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 38 * kcomp + i); 

            auto tg_xyy_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 39 * kcomp + i); 

            auto tg_xyz_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 40 * kcomp + i); 

            auto tg_xyz_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 41 * kcomp + i); 

            auto tg_xyz_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 42 * kcomp + i); 

            auto tg_xyz_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 43 * kcomp + i); 

            auto tg_xyz_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 44 * kcomp + i); 

            auto tg_xyz_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 45 * kcomp + i); 

            auto tg_xyz_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 46 * kcomp + i); 

            auto tg_xyz_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 47 * kcomp + i); 

            auto tg_xyz_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 48 * kcomp + i); 

            auto tg_xyz_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 49 * kcomp + i); 

            auto tg_xzz_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 50 * kcomp + i); 

            auto tg_xzz_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 51 * kcomp + i); 

            auto tg_xzz_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 52 * kcomp + i); 

            auto tg_xzz_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 53 * kcomp + i); 

            auto tg_xzz_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 54 * kcomp + i); 

            auto tg_xzz_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 55 * kcomp + i); 

            auto tg_xzz_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 56 * kcomp + i); 

            auto tg_xzz_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 57 * kcomp + i); 

            auto tg_xzz_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 58 * kcomp + i); 

            auto tg_xzz_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 59 * kcomp + i); 

            auto tg_yyy_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 60 * kcomp + i); 

            auto tg_yyy_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 61 * kcomp + i); 

            auto tg_yyy_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 62 * kcomp + i); 

            auto tg_yyy_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 63 * kcomp + i); 

            auto tg_yyy_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 64 * kcomp + i); 

            auto tg_yyy_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 65 * kcomp + i); 

            auto tg_yyy_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 66 * kcomp + i); 

            auto tg_yyy_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 67 * kcomp + i); 

            auto tg_yyy_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 68 * kcomp + i); 

            auto tg_yyy_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 69 * kcomp + i); 

            auto tg_yyz_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 70 * kcomp + i); 

            auto tg_yyz_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 71 * kcomp + i); 

            auto tg_yyz_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 72 * kcomp + i); 

            auto tg_yyz_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 73 * kcomp + i); 

            auto tg_yyz_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 74 * kcomp + i); 

            auto tg_yyz_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 75 * kcomp + i); 

            auto tg_yyz_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 76 * kcomp + i); 

            auto tg_yyz_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 77 * kcomp + i); 

            auto tg_yyz_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 78 * kcomp + i); 

            auto tg_yyz_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 79 * kcomp + i); 

            auto tg_yzz_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 80 * kcomp + i); 

            auto tg_yzz_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 81 * kcomp + i); 

            auto tg_yzz_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 82 * kcomp + i); 

            auto tg_yzz_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 83 * kcomp + i); 

            auto tg_yzz_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 84 * kcomp + i); 

            auto tg_yzz_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 85 * kcomp + i); 

            auto tg_yzz_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 86 * kcomp + i); 

            auto tg_yzz_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 87 * kcomp + i); 

            auto tg_yzz_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 88 * kcomp + i); 

            auto tg_yzz_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 89 * kcomp + i); 

            auto tg_zzz_xxx_0 = braBuffer.data(pidx_g_3_3_m0 + 90 * kcomp + i); 

            auto tg_zzz_xxy_0 = braBuffer.data(pidx_g_3_3_m0 + 91 * kcomp + i); 

            auto tg_zzz_xxz_0 = braBuffer.data(pidx_g_3_3_m0 + 92 * kcomp + i); 

            auto tg_zzz_xyy_0 = braBuffer.data(pidx_g_3_3_m0 + 93 * kcomp + i); 

            auto tg_zzz_xyz_0 = braBuffer.data(pidx_g_3_3_m0 + 94 * kcomp + i); 

            auto tg_zzz_xzz_0 = braBuffer.data(pidx_g_3_3_m0 + 95 * kcomp + i); 

            auto tg_zzz_yyy_0 = braBuffer.data(pidx_g_3_3_m0 + 96 * kcomp + i); 

            auto tg_zzz_yyz_0 = braBuffer.data(pidx_g_3_3_m0 + 97 * kcomp + i); 

            auto tg_zzz_yzz_0 = braBuffer.data(pidx_g_3_3_m0 + 98 * kcomp + i); 

            auto tg_zzz_zzz_0 = braBuffer.data(pidx_g_3_3_m0 + 99 * kcomp + i); 

            #pragma omp simd aligned(tg_xx_xxx_0, tg_xx_xxxx_0, tg_xx_xxxy_0, tg_xx_xxxz_0, tg_xx_xxy_0, \
                                     tg_xx_xxyy_0, tg_xx_xxyz_0, tg_xx_xxz_0, tg_xx_xxzz_0, tg_xx_xyy_0, tg_xx_xyyy_0, \
                                     tg_xx_xyyz_0, tg_xx_xyz_0, tg_xx_xyzz_0, tg_xx_xzz_0, tg_xx_xzzz_0, tg_xx_yyy_0, \
                                     tg_xx_yyz_0, tg_xx_yzz_0, tg_xx_zzz_0, tg_xxx_xxx_0, tg_xxx_xxy_0, tg_xxx_xxz_0, \
                                     tg_xxx_xyy_0, tg_xxx_xyz_0, tg_xxx_xzz_0, tg_xxx_yyy_0, tg_xxx_yyz_0, tg_xxx_yzz_0, \
                                     tg_xxx_zzz_0, tg_xxy_xxx_0, tg_xxy_xxy_0, tg_xxy_xxz_0, tg_xxy_xyy_0, tg_xxy_xyz_0, \
                                     tg_xxy_xzz_0, tg_xxy_yyy_0, tg_xxy_yyz_0, tg_xxy_yzz_0, tg_xxy_zzz_0, tg_xxz_xxx_0, \
                                     tg_xxz_xxy_0, tg_xxz_xxz_0, tg_xxz_xyy_0, tg_xxz_xyz_0, tg_xxz_xzz_0, tg_xxz_yyy_0, \
                                     tg_xxz_yyz_0, tg_xxz_yzz_0, tg_xxz_zzz_0, tg_xy_xxx_0, tg_xy_xxxx_0, tg_xy_xxxy_0, \
                                     tg_xy_xxxz_0, tg_xy_xxy_0, tg_xy_xxyy_0, tg_xy_xxyz_0, tg_xy_xxz_0, tg_xy_xxzz_0, \
                                     tg_xy_xyy_0, tg_xy_xyyy_0, tg_xy_xyyz_0, tg_xy_xyz_0, tg_xy_xyzz_0, tg_xy_xzz_0, \
                                     tg_xy_xzzz_0, tg_xy_yyy_0, tg_xy_yyz_0, tg_xy_yzz_0, tg_xy_zzz_0, tg_xyy_xxx_0, \
                                     tg_xyy_xxy_0, tg_xyy_xxz_0, tg_xyy_xyy_0, tg_xyy_xyz_0, tg_xyy_xzz_0, tg_xyy_yyy_0, \
                                     tg_xyy_yyz_0, tg_xyy_yzz_0, tg_xyy_zzz_0, tg_xyz_xxx_0, tg_xyz_xxy_0, tg_xyz_xxz_0, \
                                     tg_xyz_xyy_0, tg_xyz_xyz_0, tg_xyz_xzz_0, tg_xyz_yyy_0, tg_xyz_yyz_0, tg_xyz_yzz_0, \
                                     tg_xyz_zzz_0, tg_xz_xxx_0, tg_xz_xxxx_0, tg_xz_xxxy_0, tg_xz_xxxz_0, tg_xz_xxy_0, \
                                     tg_xz_xxyy_0, tg_xz_xxyz_0, tg_xz_xxz_0, tg_xz_xxzz_0, tg_xz_xyy_0, tg_xz_xyyy_0, \
                                     tg_xz_xyyz_0, tg_xz_xyz_0, tg_xz_xyzz_0, tg_xz_xzz_0, tg_xz_xzzz_0, tg_xz_yyy_0, \
                                     tg_xz_yyz_0, tg_xz_yzz_0, tg_xz_zzz_0, tg_xzz_xxx_0, tg_xzz_xxy_0, tg_xzz_xxz_0, \
                                     tg_xzz_xyy_0, tg_xzz_xyz_0, tg_xzz_xzz_0, tg_xzz_yyy_0, tg_xzz_yyz_0, tg_xzz_yzz_0, \
                                     tg_xzz_zzz_0, tg_yy_xxx_0, tg_yy_xxxx_0, tg_yy_xxxy_0, tg_yy_xxxz_0, tg_yy_xxy_0, \
                                     tg_yy_xxyy_0, tg_yy_xxyz_0, tg_yy_xxz_0, tg_yy_xxzz_0, tg_yy_xyy_0, tg_yy_xyyy_0, \
                                     tg_yy_xyyz_0, tg_yy_xyz_0, tg_yy_xyzz_0, tg_yy_xzz_0, tg_yy_xzzz_0, tg_yy_yyy_0, \
                                     tg_yy_yyyy_0, tg_yy_yyyz_0, tg_yy_yyz_0, tg_yy_yyzz_0, tg_yy_yzz_0, tg_yy_yzzz_0, \
                                     tg_yy_zzz_0, tg_yyy_xxx_0, tg_yyy_xxy_0, tg_yyy_xxz_0, tg_yyy_xyy_0, tg_yyy_xyz_0, \
                                     tg_yyy_xzz_0, tg_yyy_yyy_0, tg_yyy_yyz_0, tg_yyy_yzz_0, tg_yyy_zzz_0, tg_yyz_xxx_0, \
                                     tg_yyz_xxy_0, tg_yyz_xxz_0, tg_yyz_xyy_0, tg_yyz_xyz_0, tg_yyz_xzz_0, tg_yyz_yyy_0, \
                                     tg_yyz_yyz_0, tg_yyz_yzz_0, tg_yyz_zzz_0, tg_yz_xxx_0, tg_yz_xxxx_0, tg_yz_xxxy_0, \
                                     tg_yz_xxxz_0, tg_yz_xxy_0, tg_yz_xxyy_0, tg_yz_xxyz_0, tg_yz_xxz_0, tg_yz_xxzz_0, \
                                     tg_yz_xyy_0, tg_yz_xyyy_0, tg_yz_xyyz_0, tg_yz_xyz_0, tg_yz_xyzz_0, tg_yz_xzz_0, \
                                     tg_yz_xzzz_0, tg_yz_yyy_0, tg_yz_yyyy_0, tg_yz_yyyz_0, tg_yz_yyz_0, tg_yz_yyzz_0, \
                                     tg_yz_yzz_0, tg_yz_yzzz_0, tg_yz_zzz_0, tg_yzz_xxx_0, tg_yzz_xxy_0, tg_yzz_xxz_0, \
                                     tg_yzz_xyy_0, tg_yzz_xyz_0, tg_yzz_xzz_0, tg_yzz_yyy_0, tg_yzz_yyz_0, tg_yzz_yzz_0, \
                                     tg_yzz_zzz_0, tg_zz_xxx_0, tg_zz_xxxx_0, tg_zz_xxxy_0, tg_zz_xxxz_0, tg_zz_xxy_0, \
                                     tg_zz_xxyy_0, tg_zz_xxyz_0, tg_zz_xxz_0, tg_zz_xxzz_0, tg_zz_xyy_0, tg_zz_xyyy_0, \
                                     tg_zz_xyyz_0, tg_zz_xyz_0, tg_zz_xyzz_0, tg_zz_xzz_0, tg_zz_xzzz_0, tg_zz_yyy_0, \
                                     tg_zz_yyyy_0, tg_zz_yyyz_0, tg_zz_yyz_0, tg_zz_yyzz_0, tg_zz_yzz_0, tg_zz_yzzz_0, \
                                     tg_zz_zzz_0, tg_zz_zzzz_0, tg_zzz_xxx_0, tg_zzz_xxy_0, tg_zzz_xxz_0, tg_zzz_xyy_0, \
                                     tg_zzz_xyz_0, tg_zzz_xzz_0, tg_zzz_yyy_0, tg_zzz_yyz_0, tg_zzz_yzz_0, tg_zzz_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xxx_xxx_0[j] = -ab_x * tg_xx_xxx_0[j] + tg_xx_xxxx_0[j];

                tg_xxx_xxy_0[j] = -ab_x * tg_xx_xxy_0[j] + tg_xx_xxxy_0[j];

                tg_xxx_xxz_0[j] = -ab_x * tg_xx_xxz_0[j] + tg_xx_xxxz_0[j];

                tg_xxx_xyy_0[j] = -ab_x * tg_xx_xyy_0[j] + tg_xx_xxyy_0[j];

                tg_xxx_xyz_0[j] = -ab_x * tg_xx_xyz_0[j] + tg_xx_xxyz_0[j];

                tg_xxx_xzz_0[j] = -ab_x * tg_xx_xzz_0[j] + tg_xx_xxzz_0[j];

                tg_xxx_yyy_0[j] = -ab_x * tg_xx_yyy_0[j] + tg_xx_xyyy_0[j];

                tg_xxx_yyz_0[j] = -ab_x * tg_xx_yyz_0[j] + tg_xx_xyyz_0[j];

                tg_xxx_yzz_0[j] = -ab_x * tg_xx_yzz_0[j] + tg_xx_xyzz_0[j];

                tg_xxx_zzz_0[j] = -ab_x * tg_xx_zzz_0[j] + tg_xx_xzzz_0[j];

                tg_xxy_xxx_0[j] = -ab_x * tg_xy_xxx_0[j] + tg_xy_xxxx_0[j];

                tg_xxy_xxy_0[j] = -ab_x * tg_xy_xxy_0[j] + tg_xy_xxxy_0[j];

                tg_xxy_xxz_0[j] = -ab_x * tg_xy_xxz_0[j] + tg_xy_xxxz_0[j];

                tg_xxy_xyy_0[j] = -ab_x * tg_xy_xyy_0[j] + tg_xy_xxyy_0[j];

                tg_xxy_xyz_0[j] = -ab_x * tg_xy_xyz_0[j] + tg_xy_xxyz_0[j];

                tg_xxy_xzz_0[j] = -ab_x * tg_xy_xzz_0[j] + tg_xy_xxzz_0[j];

                tg_xxy_yyy_0[j] = -ab_x * tg_xy_yyy_0[j] + tg_xy_xyyy_0[j];

                tg_xxy_yyz_0[j] = -ab_x * tg_xy_yyz_0[j] + tg_xy_xyyz_0[j];

                tg_xxy_yzz_0[j] = -ab_x * tg_xy_yzz_0[j] + tg_xy_xyzz_0[j];

                tg_xxy_zzz_0[j] = -ab_x * tg_xy_zzz_0[j] + tg_xy_xzzz_0[j];

                tg_xxz_xxx_0[j] = -ab_x * tg_xz_xxx_0[j] + tg_xz_xxxx_0[j];

                tg_xxz_xxy_0[j] = -ab_x * tg_xz_xxy_0[j] + tg_xz_xxxy_0[j];

                tg_xxz_xxz_0[j] = -ab_x * tg_xz_xxz_0[j] + tg_xz_xxxz_0[j];

                tg_xxz_xyy_0[j] = -ab_x * tg_xz_xyy_0[j] + tg_xz_xxyy_0[j];

                tg_xxz_xyz_0[j] = -ab_x * tg_xz_xyz_0[j] + tg_xz_xxyz_0[j];

                tg_xxz_xzz_0[j] = -ab_x * tg_xz_xzz_0[j] + tg_xz_xxzz_0[j];

                tg_xxz_yyy_0[j] = -ab_x * tg_xz_yyy_0[j] + tg_xz_xyyy_0[j];

                tg_xxz_yyz_0[j] = -ab_x * tg_xz_yyz_0[j] + tg_xz_xyyz_0[j];

                tg_xxz_yzz_0[j] = -ab_x * tg_xz_yzz_0[j] + tg_xz_xyzz_0[j];

                tg_xxz_zzz_0[j] = -ab_x * tg_xz_zzz_0[j] + tg_xz_xzzz_0[j];

                tg_xyy_xxx_0[j] = -ab_x * tg_yy_xxx_0[j] + tg_yy_xxxx_0[j];

                tg_xyy_xxy_0[j] = -ab_x * tg_yy_xxy_0[j] + tg_yy_xxxy_0[j];

                tg_xyy_xxz_0[j] = -ab_x * tg_yy_xxz_0[j] + tg_yy_xxxz_0[j];

                tg_xyy_xyy_0[j] = -ab_x * tg_yy_xyy_0[j] + tg_yy_xxyy_0[j];

                tg_xyy_xyz_0[j] = -ab_x * tg_yy_xyz_0[j] + tg_yy_xxyz_0[j];

                tg_xyy_xzz_0[j] = -ab_x * tg_yy_xzz_0[j] + tg_yy_xxzz_0[j];

                tg_xyy_yyy_0[j] = -ab_x * tg_yy_yyy_0[j] + tg_yy_xyyy_0[j];

                tg_xyy_yyz_0[j] = -ab_x * tg_yy_yyz_0[j] + tg_yy_xyyz_0[j];

                tg_xyy_yzz_0[j] = -ab_x * tg_yy_yzz_0[j] + tg_yy_xyzz_0[j];

                tg_xyy_zzz_0[j] = -ab_x * tg_yy_zzz_0[j] + tg_yy_xzzz_0[j];

                tg_xyz_xxx_0[j] = -ab_x * tg_yz_xxx_0[j] + tg_yz_xxxx_0[j];

                tg_xyz_xxy_0[j] = -ab_x * tg_yz_xxy_0[j] + tg_yz_xxxy_0[j];

                tg_xyz_xxz_0[j] = -ab_x * tg_yz_xxz_0[j] + tg_yz_xxxz_0[j];

                tg_xyz_xyy_0[j] = -ab_x * tg_yz_xyy_0[j] + tg_yz_xxyy_0[j];

                tg_xyz_xyz_0[j] = -ab_x * tg_yz_xyz_0[j] + tg_yz_xxyz_0[j];

                tg_xyz_xzz_0[j] = -ab_x * tg_yz_xzz_0[j] + tg_yz_xxzz_0[j];

                tg_xyz_yyy_0[j] = -ab_x * tg_yz_yyy_0[j] + tg_yz_xyyy_0[j];

                tg_xyz_yyz_0[j] = -ab_x * tg_yz_yyz_0[j] + tg_yz_xyyz_0[j];

                tg_xyz_yzz_0[j] = -ab_x * tg_yz_yzz_0[j] + tg_yz_xyzz_0[j];

                tg_xyz_zzz_0[j] = -ab_x * tg_yz_zzz_0[j] + tg_yz_xzzz_0[j];

                tg_xzz_xxx_0[j] = -ab_x * tg_zz_xxx_0[j] + tg_zz_xxxx_0[j];

                tg_xzz_xxy_0[j] = -ab_x * tg_zz_xxy_0[j] + tg_zz_xxxy_0[j];

                tg_xzz_xxz_0[j] = -ab_x * tg_zz_xxz_0[j] + tg_zz_xxxz_0[j];

                tg_xzz_xyy_0[j] = -ab_x * tg_zz_xyy_0[j] + tg_zz_xxyy_0[j];

                tg_xzz_xyz_0[j] = -ab_x * tg_zz_xyz_0[j] + tg_zz_xxyz_0[j];

                tg_xzz_xzz_0[j] = -ab_x * tg_zz_xzz_0[j] + tg_zz_xxzz_0[j];

                tg_xzz_yyy_0[j] = -ab_x * tg_zz_yyy_0[j] + tg_zz_xyyy_0[j];

                tg_xzz_yyz_0[j] = -ab_x * tg_zz_yyz_0[j] + tg_zz_xyyz_0[j];

                tg_xzz_yzz_0[j] = -ab_x * tg_zz_yzz_0[j] + tg_zz_xyzz_0[j];

                tg_xzz_zzz_0[j] = -ab_x * tg_zz_zzz_0[j] + tg_zz_xzzz_0[j];

                tg_yyy_xxx_0[j] = -ab_y * tg_yy_xxx_0[j] + tg_yy_xxxy_0[j];

                tg_yyy_xxy_0[j] = -ab_y * tg_yy_xxy_0[j] + tg_yy_xxyy_0[j];

                tg_yyy_xxz_0[j] = -ab_y * tg_yy_xxz_0[j] + tg_yy_xxyz_0[j];

                tg_yyy_xyy_0[j] = -ab_y * tg_yy_xyy_0[j] + tg_yy_xyyy_0[j];

                tg_yyy_xyz_0[j] = -ab_y * tg_yy_xyz_0[j] + tg_yy_xyyz_0[j];

                tg_yyy_xzz_0[j] = -ab_y * tg_yy_xzz_0[j] + tg_yy_xyzz_0[j];

                tg_yyy_yyy_0[j] = -ab_y * tg_yy_yyy_0[j] + tg_yy_yyyy_0[j];

                tg_yyy_yyz_0[j] = -ab_y * tg_yy_yyz_0[j] + tg_yy_yyyz_0[j];

                tg_yyy_yzz_0[j] = -ab_y * tg_yy_yzz_0[j] + tg_yy_yyzz_0[j];

                tg_yyy_zzz_0[j] = -ab_y * tg_yy_zzz_0[j] + tg_yy_yzzz_0[j];

                tg_yyz_xxx_0[j] = -ab_y * tg_yz_xxx_0[j] + tg_yz_xxxy_0[j];

                tg_yyz_xxy_0[j] = -ab_y * tg_yz_xxy_0[j] + tg_yz_xxyy_0[j];

                tg_yyz_xxz_0[j] = -ab_y * tg_yz_xxz_0[j] + tg_yz_xxyz_0[j];

                tg_yyz_xyy_0[j] = -ab_y * tg_yz_xyy_0[j] + tg_yz_xyyy_0[j];

                tg_yyz_xyz_0[j] = -ab_y * tg_yz_xyz_0[j] + tg_yz_xyyz_0[j];

                tg_yyz_xzz_0[j] = -ab_y * tg_yz_xzz_0[j] + tg_yz_xyzz_0[j];

                tg_yyz_yyy_0[j] = -ab_y * tg_yz_yyy_0[j] + tg_yz_yyyy_0[j];

                tg_yyz_yyz_0[j] = -ab_y * tg_yz_yyz_0[j] + tg_yz_yyyz_0[j];

                tg_yyz_yzz_0[j] = -ab_y * tg_yz_yzz_0[j] + tg_yz_yyzz_0[j];

                tg_yyz_zzz_0[j] = -ab_y * tg_yz_zzz_0[j] + tg_yz_yzzz_0[j];

                tg_yzz_xxx_0[j] = -ab_y * tg_zz_xxx_0[j] + tg_zz_xxxy_0[j];

                tg_yzz_xxy_0[j] = -ab_y * tg_zz_xxy_0[j] + tg_zz_xxyy_0[j];

                tg_yzz_xxz_0[j] = -ab_y * tg_zz_xxz_0[j] + tg_zz_xxyz_0[j];

                tg_yzz_xyy_0[j] = -ab_y * tg_zz_xyy_0[j] + tg_zz_xyyy_0[j];

                tg_yzz_xyz_0[j] = -ab_y * tg_zz_xyz_0[j] + tg_zz_xyyz_0[j];

                tg_yzz_xzz_0[j] = -ab_y * tg_zz_xzz_0[j] + tg_zz_xyzz_0[j];

                tg_yzz_yyy_0[j] = -ab_y * tg_zz_yyy_0[j] + tg_zz_yyyy_0[j];

                tg_yzz_yyz_0[j] = -ab_y * tg_zz_yyz_0[j] + tg_zz_yyyz_0[j];

                tg_yzz_yzz_0[j] = -ab_y * tg_zz_yzz_0[j] + tg_zz_yyzz_0[j];

                tg_yzz_zzz_0[j] = -ab_y * tg_zz_zzz_0[j] + tg_zz_yzzz_0[j];

                tg_zzz_xxx_0[j] = -ab_z * tg_zz_xxx_0[j] + tg_zz_xxxz_0[j];

                tg_zzz_xxy_0[j] = -ab_z * tg_zz_xxy_0[j] + tg_zz_xxyz_0[j];

                tg_zzz_xxz_0[j] = -ab_z * tg_zz_xxz_0[j] + tg_zz_xxzz_0[j];

                tg_zzz_xyy_0[j] = -ab_z * tg_zz_xyy_0[j] + tg_zz_xyyz_0[j];

                tg_zzz_xyz_0[j] = -ab_z * tg_zz_xyz_0[j] + tg_zz_xyzz_0[j];

                tg_zzz_xzz_0[j] = -ab_z * tg_zz_xzz_0[j] + tg_zz_xzzz_0[j];

                tg_zzz_yyy_0[j] = -ab_z * tg_zz_yyy_0[j] + tg_zz_yyyz_0[j];

                tg_zzz_yyz_0[j] = -ab_z * tg_zz_yyz_0[j] + tg_zz_yyzz_0[j];

                tg_zzz_yzz_0[j] = -ab_z * tg_zz_yzz_0[j] + tg_zz_yzzz_0[j];

                tg_zzz_zzz_0[j] = -ab_z * tg_zz_zzz_0[j] + tg_zz_zzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForFGXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        eribrrfunc::compElectronRepulsionForFGXY_0_75(braBuffer,
                                                      recursionMap,
                                                      abDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        eribrrfunc::compElectronRepulsionForFGXY_75_150(braBuffer,
                                                        recursionMap,
                                                        abDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForFGXY_0_75(      CMemBlock2D<double>& braBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& abDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,75)

        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_3_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            // set up pointers to integrals

            auto tg_xxx_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + i); 

            auto tg_xxx_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + kcomp + i); 

            auto tg_xxx_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 2 * kcomp + i); 

            auto tg_xxx_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 3 * kcomp + i); 

            auto tg_xxx_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 4 * kcomp + i); 

            auto tg_xxx_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 5 * kcomp + i); 

            auto tg_xxx_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 6 * kcomp + i); 

            auto tg_xxx_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 7 * kcomp + i); 

            auto tg_xxx_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 8 * kcomp + i); 

            auto tg_xxx_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 9 * kcomp + i); 

            auto tg_xxx_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 10 * kcomp + i); 

            auto tg_xxx_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 11 * kcomp + i); 

            auto tg_xxx_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 12 * kcomp + i); 

            auto tg_xxx_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 13 * kcomp + i); 

            auto tg_xxx_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 14 * kcomp + i); 

            auto tg_xxy_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 15 * kcomp + i); 

            auto tg_xxy_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 16 * kcomp + i); 

            auto tg_xxy_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 17 * kcomp + i); 

            auto tg_xxy_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 18 * kcomp + i); 

            auto tg_xxy_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 19 * kcomp + i); 

            auto tg_xxy_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 20 * kcomp + i); 

            auto tg_xxy_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 21 * kcomp + i); 

            auto tg_xxy_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 22 * kcomp + i); 

            auto tg_xxy_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 23 * kcomp + i); 

            auto tg_xxy_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 24 * kcomp + i); 

            auto tg_xxy_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 25 * kcomp + i); 

            auto tg_xxy_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 26 * kcomp + i); 

            auto tg_xxy_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 27 * kcomp + i); 

            auto tg_xxy_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 28 * kcomp + i); 

            auto tg_xxy_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 29 * kcomp + i); 

            auto tg_xxz_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 30 * kcomp + i); 

            auto tg_xxz_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 31 * kcomp + i); 

            auto tg_xxz_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 32 * kcomp + i); 

            auto tg_xxz_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 33 * kcomp + i); 

            auto tg_xxz_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 34 * kcomp + i); 

            auto tg_xxz_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 35 * kcomp + i); 

            auto tg_xxz_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 36 * kcomp + i); 

            auto tg_xxz_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 37 * kcomp + i); 

            auto tg_xxz_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 38 * kcomp + i); 

            auto tg_xxz_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 39 * kcomp + i); 

            auto tg_xxz_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 40 * kcomp + i); 

            auto tg_xxz_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 41 * kcomp + i); 

            auto tg_xxz_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 42 * kcomp + i); 

            auto tg_xxz_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 43 * kcomp + i); 

            auto tg_xxz_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 44 * kcomp + i); 

            auto tg_xyy_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 45 * kcomp + i); 

            auto tg_xyy_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 46 * kcomp + i); 

            auto tg_xyy_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 47 * kcomp + i); 

            auto tg_xyy_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 48 * kcomp + i); 

            auto tg_xyy_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 49 * kcomp + i); 

            auto tg_xyy_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 50 * kcomp + i); 

            auto tg_xyy_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 51 * kcomp + i); 

            auto tg_xyy_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 52 * kcomp + i); 

            auto tg_xyy_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 53 * kcomp + i); 

            auto tg_xyy_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 54 * kcomp + i); 

            auto tg_xyy_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 55 * kcomp + i); 

            auto tg_xyy_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 56 * kcomp + i); 

            auto tg_xyy_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 57 * kcomp + i); 

            auto tg_xyy_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 58 * kcomp + i); 

            auto tg_xyy_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 59 * kcomp + i); 

            auto tg_xyz_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 60 * kcomp + i); 

            auto tg_xyz_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 61 * kcomp + i); 

            auto tg_xyz_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 62 * kcomp + i); 

            auto tg_xyz_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 63 * kcomp + i); 

            auto tg_xyz_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 64 * kcomp + i); 

            auto tg_xyz_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 65 * kcomp + i); 

            auto tg_xyz_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 66 * kcomp + i); 

            auto tg_xyz_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 67 * kcomp + i); 

            auto tg_xyz_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 68 * kcomp + i); 

            auto tg_xyz_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 69 * kcomp + i); 

            auto tg_xyz_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 70 * kcomp + i); 

            auto tg_xyz_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 71 * kcomp + i); 

            auto tg_xyz_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 72 * kcomp + i); 

            auto tg_xyz_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 73 * kcomp + i); 

            auto tg_xyz_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 74 * kcomp + i); 

            // Batch of Integrals (0,75)

            #pragma omp simd aligned(tg_xx_xxxx_0, tg_xx_xxxxx_0, tg_xx_xxxxy_0, tg_xx_xxxxz_0, \
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
                tg_xxx_xxxx_0[j] = -ab_x * tg_xx_xxxx_0[j] + tg_xx_xxxxx_0[j];

                tg_xxx_xxxy_0[j] = -ab_x * tg_xx_xxxy_0[j] + tg_xx_xxxxy_0[j];

                tg_xxx_xxxz_0[j] = -ab_x * tg_xx_xxxz_0[j] + tg_xx_xxxxz_0[j];

                tg_xxx_xxyy_0[j] = -ab_x * tg_xx_xxyy_0[j] + tg_xx_xxxyy_0[j];

                tg_xxx_xxyz_0[j] = -ab_x * tg_xx_xxyz_0[j] + tg_xx_xxxyz_0[j];

                tg_xxx_xxzz_0[j] = -ab_x * tg_xx_xxzz_0[j] + tg_xx_xxxzz_0[j];

                tg_xxx_xyyy_0[j] = -ab_x * tg_xx_xyyy_0[j] + tg_xx_xxyyy_0[j];

                tg_xxx_xyyz_0[j] = -ab_x * tg_xx_xyyz_0[j] + tg_xx_xxyyz_0[j];

                tg_xxx_xyzz_0[j] = -ab_x * tg_xx_xyzz_0[j] + tg_xx_xxyzz_0[j];

                tg_xxx_xzzz_0[j] = -ab_x * tg_xx_xzzz_0[j] + tg_xx_xxzzz_0[j];

                tg_xxx_yyyy_0[j] = -ab_x * tg_xx_yyyy_0[j] + tg_xx_xyyyy_0[j];

                tg_xxx_yyyz_0[j] = -ab_x * tg_xx_yyyz_0[j] + tg_xx_xyyyz_0[j];

                tg_xxx_yyzz_0[j] = -ab_x * tg_xx_yyzz_0[j] + tg_xx_xyyzz_0[j];

                tg_xxx_yzzz_0[j] = -ab_x * tg_xx_yzzz_0[j] + tg_xx_xyzzz_0[j];

                tg_xxx_zzzz_0[j] = -ab_x * tg_xx_zzzz_0[j] + tg_xx_xzzzz_0[j];

                tg_xxy_xxxx_0[j] = -ab_x * tg_xy_xxxx_0[j] + tg_xy_xxxxx_0[j];

                tg_xxy_xxxy_0[j] = -ab_x * tg_xy_xxxy_0[j] + tg_xy_xxxxy_0[j];

                tg_xxy_xxxz_0[j] = -ab_x * tg_xy_xxxz_0[j] + tg_xy_xxxxz_0[j];

                tg_xxy_xxyy_0[j] = -ab_x * tg_xy_xxyy_0[j] + tg_xy_xxxyy_0[j];

                tg_xxy_xxyz_0[j] = -ab_x * tg_xy_xxyz_0[j] + tg_xy_xxxyz_0[j];

                tg_xxy_xxzz_0[j] = -ab_x * tg_xy_xxzz_0[j] + tg_xy_xxxzz_0[j];

                tg_xxy_xyyy_0[j] = -ab_x * tg_xy_xyyy_0[j] + tg_xy_xxyyy_0[j];

                tg_xxy_xyyz_0[j] = -ab_x * tg_xy_xyyz_0[j] + tg_xy_xxyyz_0[j];

                tg_xxy_xyzz_0[j] = -ab_x * tg_xy_xyzz_0[j] + tg_xy_xxyzz_0[j];

                tg_xxy_xzzz_0[j] = -ab_x * tg_xy_xzzz_0[j] + tg_xy_xxzzz_0[j];

                tg_xxy_yyyy_0[j] = -ab_x * tg_xy_yyyy_0[j] + tg_xy_xyyyy_0[j];

                tg_xxy_yyyz_0[j] = -ab_x * tg_xy_yyyz_0[j] + tg_xy_xyyyz_0[j];

                tg_xxy_yyzz_0[j] = -ab_x * tg_xy_yyzz_0[j] + tg_xy_xyyzz_0[j];

                tg_xxy_yzzz_0[j] = -ab_x * tg_xy_yzzz_0[j] + tg_xy_xyzzz_0[j];

                tg_xxy_zzzz_0[j] = -ab_x * tg_xy_zzzz_0[j] + tg_xy_xzzzz_0[j];

                tg_xxz_xxxx_0[j] = -ab_x * tg_xz_xxxx_0[j] + tg_xz_xxxxx_0[j];

                tg_xxz_xxxy_0[j] = -ab_x * tg_xz_xxxy_0[j] + tg_xz_xxxxy_0[j];

                tg_xxz_xxxz_0[j] = -ab_x * tg_xz_xxxz_0[j] + tg_xz_xxxxz_0[j];

                tg_xxz_xxyy_0[j] = -ab_x * tg_xz_xxyy_0[j] + tg_xz_xxxyy_0[j];

                tg_xxz_xxyz_0[j] = -ab_x * tg_xz_xxyz_0[j] + tg_xz_xxxyz_0[j];

                tg_xxz_xxzz_0[j] = -ab_x * tg_xz_xxzz_0[j] + tg_xz_xxxzz_0[j];

                tg_xxz_xyyy_0[j] = -ab_x * tg_xz_xyyy_0[j] + tg_xz_xxyyy_0[j];

                tg_xxz_xyyz_0[j] = -ab_x * tg_xz_xyyz_0[j] + tg_xz_xxyyz_0[j];

                tg_xxz_xyzz_0[j] = -ab_x * tg_xz_xyzz_0[j] + tg_xz_xxyzz_0[j];

                tg_xxz_xzzz_0[j] = -ab_x * tg_xz_xzzz_0[j] + tg_xz_xxzzz_0[j];

                tg_xxz_yyyy_0[j] = -ab_x * tg_xz_yyyy_0[j] + tg_xz_xyyyy_0[j];

                tg_xxz_yyyz_0[j] = -ab_x * tg_xz_yyyz_0[j] + tg_xz_xyyyz_0[j];

                tg_xxz_yyzz_0[j] = -ab_x * tg_xz_yyzz_0[j] + tg_xz_xyyzz_0[j];

                tg_xxz_yzzz_0[j] = -ab_x * tg_xz_yzzz_0[j] + tg_xz_xyzzz_0[j];

                tg_xxz_zzzz_0[j] = -ab_x * tg_xz_zzzz_0[j] + tg_xz_xzzzz_0[j];

                tg_xyy_xxxx_0[j] = -ab_x * tg_yy_xxxx_0[j] + tg_yy_xxxxx_0[j];

                tg_xyy_xxxy_0[j] = -ab_x * tg_yy_xxxy_0[j] + tg_yy_xxxxy_0[j];

                tg_xyy_xxxz_0[j] = -ab_x * tg_yy_xxxz_0[j] + tg_yy_xxxxz_0[j];

                tg_xyy_xxyy_0[j] = -ab_x * tg_yy_xxyy_0[j] + tg_yy_xxxyy_0[j];

                tg_xyy_xxyz_0[j] = -ab_x * tg_yy_xxyz_0[j] + tg_yy_xxxyz_0[j];

                tg_xyy_xxzz_0[j] = -ab_x * tg_yy_xxzz_0[j] + tg_yy_xxxzz_0[j];

                tg_xyy_xyyy_0[j] = -ab_x * tg_yy_xyyy_0[j] + tg_yy_xxyyy_0[j];

                tg_xyy_xyyz_0[j] = -ab_x * tg_yy_xyyz_0[j] + tg_yy_xxyyz_0[j];

                tg_xyy_xyzz_0[j] = -ab_x * tg_yy_xyzz_0[j] + tg_yy_xxyzz_0[j];

                tg_xyy_xzzz_0[j] = -ab_x * tg_yy_xzzz_0[j] + tg_yy_xxzzz_0[j];

                tg_xyy_yyyy_0[j] = -ab_x * tg_yy_yyyy_0[j] + tg_yy_xyyyy_0[j];

                tg_xyy_yyyz_0[j] = -ab_x * tg_yy_yyyz_0[j] + tg_yy_xyyyz_0[j];

                tg_xyy_yyzz_0[j] = -ab_x * tg_yy_yyzz_0[j] + tg_yy_xyyzz_0[j];

                tg_xyy_yzzz_0[j] = -ab_x * tg_yy_yzzz_0[j] + tg_yy_xyzzz_0[j];

                tg_xyy_zzzz_0[j] = -ab_x * tg_yy_zzzz_0[j] + tg_yy_xzzzz_0[j];

                tg_xyz_xxxx_0[j] = -ab_x * tg_yz_xxxx_0[j] + tg_yz_xxxxx_0[j];

                tg_xyz_xxxy_0[j] = -ab_x * tg_yz_xxxy_0[j] + tg_yz_xxxxy_0[j];

                tg_xyz_xxxz_0[j] = -ab_x * tg_yz_xxxz_0[j] + tg_yz_xxxxz_0[j];

                tg_xyz_xxyy_0[j] = -ab_x * tg_yz_xxyy_0[j] + tg_yz_xxxyy_0[j];

                tg_xyz_xxyz_0[j] = -ab_x * tg_yz_xxyz_0[j] + tg_yz_xxxyz_0[j];

                tg_xyz_xxzz_0[j] = -ab_x * tg_yz_xxzz_0[j] + tg_yz_xxxzz_0[j];

                tg_xyz_xyyy_0[j] = -ab_x * tg_yz_xyyy_0[j] + tg_yz_xxyyy_0[j];

                tg_xyz_xyyz_0[j] = -ab_x * tg_yz_xyyz_0[j] + tg_yz_xxyyz_0[j];

                tg_xyz_xyzz_0[j] = -ab_x * tg_yz_xyzz_0[j] + tg_yz_xxyzz_0[j];

                tg_xyz_xzzz_0[j] = -ab_x * tg_yz_xzzz_0[j] + tg_yz_xxzzz_0[j];

                tg_xyz_yyyy_0[j] = -ab_x * tg_yz_yyyy_0[j] + tg_yz_xyyyy_0[j];

                tg_xyz_yyyz_0[j] = -ab_x * tg_yz_yyyz_0[j] + tg_yz_xyyyz_0[j];

                tg_xyz_yyzz_0[j] = -ab_x * tg_yz_yyzz_0[j] + tg_yz_xyyzz_0[j];

                tg_xyz_yzzz_0[j] = -ab_x * tg_yz_yzzz_0[j] + tg_yz_xyzzz_0[j];

                tg_xyz_zzzz_0[j] = -ab_x * tg_yz_zzzz_0[j] + tg_yz_xzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForFGXY_75_150(      CMemBlock2D<double>& braBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& abDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetContrPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (75,150)

        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_3_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            auto tg_yy_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + 64 * kcomp + i); 

            auto tg_yy_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 66 * kcomp + i); 

            auto tg_yy_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 67 * kcomp + i); 

            auto tg_yy_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 69 * kcomp + i); 

            auto tg_yy_xxyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 70 * kcomp + i); 

            auto tg_yy_xxyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 71 * kcomp + i); 

            auto tg_yy_xyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 73 * kcomp + i); 

            auto tg_yy_xyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 74 * kcomp + i); 

            auto tg_yy_xyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 75 * kcomp + i); 

            auto tg_yy_xyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 76 * kcomp + i); 

            auto tg_yy_yyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 78 * kcomp + i); 

            auto tg_yy_yyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 79 * kcomp + i); 

            auto tg_yy_yyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 80 * kcomp + i); 

            auto tg_yy_yyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 81 * kcomp + i); 

            auto tg_yy_yzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 82 * kcomp + i); 

            auto tg_yz_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + 85 * kcomp + i); 

            auto tg_yz_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 87 * kcomp + i); 

            auto tg_yz_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 88 * kcomp + i); 

            auto tg_yz_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 90 * kcomp + i); 

            auto tg_yz_xxyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 91 * kcomp + i); 

            auto tg_yz_xxyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 92 * kcomp + i); 

            auto tg_yz_xyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 94 * kcomp + i); 

            auto tg_yz_xyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 95 * kcomp + i); 

            auto tg_yz_xyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 96 * kcomp + i); 

            auto tg_yz_xyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 97 * kcomp + i); 

            auto tg_yz_yyyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 99 * kcomp + i); 

            auto tg_yz_yyyyz_0 = braBuffer.data(pidx_g_2_5_m0 + 100 * kcomp + i); 

            auto tg_yz_yyyzz_0 = braBuffer.data(pidx_g_2_5_m0 + 101 * kcomp + i); 

            auto tg_yz_yyzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 102 * kcomp + i); 

            auto tg_yz_yzzzz_0 = braBuffer.data(pidx_g_2_5_m0 + 103 * kcomp + i); 

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

            // set up pointers to integrals

            auto tg_xzz_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 75 * kcomp + i); 

            auto tg_xzz_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 76 * kcomp + i); 

            auto tg_xzz_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 77 * kcomp + i); 

            auto tg_xzz_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 78 * kcomp + i); 

            auto tg_xzz_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 79 * kcomp + i); 

            auto tg_xzz_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 80 * kcomp + i); 

            auto tg_xzz_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 81 * kcomp + i); 

            auto tg_xzz_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 82 * kcomp + i); 

            auto tg_xzz_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 83 * kcomp + i); 

            auto tg_xzz_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 84 * kcomp + i); 

            auto tg_xzz_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 85 * kcomp + i); 

            auto tg_xzz_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 86 * kcomp + i); 

            auto tg_xzz_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 87 * kcomp + i); 

            auto tg_xzz_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 88 * kcomp + i); 

            auto tg_xzz_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 89 * kcomp + i); 

            auto tg_yyy_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 90 * kcomp + i); 

            auto tg_yyy_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 91 * kcomp + i); 

            auto tg_yyy_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 92 * kcomp + i); 

            auto tg_yyy_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 93 * kcomp + i); 

            auto tg_yyy_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 94 * kcomp + i); 

            auto tg_yyy_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 95 * kcomp + i); 

            auto tg_yyy_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 96 * kcomp + i); 

            auto tg_yyy_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 97 * kcomp + i); 

            auto tg_yyy_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 98 * kcomp + i); 

            auto tg_yyy_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 99 * kcomp + i); 

            auto tg_yyy_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 100 * kcomp + i); 

            auto tg_yyy_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 101 * kcomp + i); 

            auto tg_yyy_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 102 * kcomp + i); 

            auto tg_yyy_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 103 * kcomp + i); 

            auto tg_yyy_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 104 * kcomp + i); 

            auto tg_yyz_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 105 * kcomp + i); 

            auto tg_yyz_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 106 * kcomp + i); 

            auto tg_yyz_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 107 * kcomp + i); 

            auto tg_yyz_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 108 * kcomp + i); 

            auto tg_yyz_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 109 * kcomp + i); 

            auto tg_yyz_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 110 * kcomp + i); 

            auto tg_yyz_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 111 * kcomp + i); 

            auto tg_yyz_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 112 * kcomp + i); 

            auto tg_yyz_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 113 * kcomp + i); 

            auto tg_yyz_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 114 * kcomp + i); 

            auto tg_yyz_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 115 * kcomp + i); 

            auto tg_yyz_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 116 * kcomp + i); 

            auto tg_yyz_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 117 * kcomp + i); 

            auto tg_yyz_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 118 * kcomp + i); 

            auto tg_yyz_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 119 * kcomp + i); 

            auto tg_yzz_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 120 * kcomp + i); 

            auto tg_yzz_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 121 * kcomp + i); 

            auto tg_yzz_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 122 * kcomp + i); 

            auto tg_yzz_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 123 * kcomp + i); 

            auto tg_yzz_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 124 * kcomp + i); 

            auto tg_yzz_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 125 * kcomp + i); 

            auto tg_yzz_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 126 * kcomp + i); 

            auto tg_yzz_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 127 * kcomp + i); 

            auto tg_yzz_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 128 * kcomp + i); 

            auto tg_yzz_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 129 * kcomp + i); 

            auto tg_yzz_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 130 * kcomp + i); 

            auto tg_yzz_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 131 * kcomp + i); 

            auto tg_yzz_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 132 * kcomp + i); 

            auto tg_yzz_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 133 * kcomp + i); 

            auto tg_yzz_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 134 * kcomp + i); 

            auto tg_zzz_xxxx_0 = braBuffer.data(pidx_g_3_4_m0 + 135 * kcomp + i); 

            auto tg_zzz_xxxy_0 = braBuffer.data(pidx_g_3_4_m0 + 136 * kcomp + i); 

            auto tg_zzz_xxxz_0 = braBuffer.data(pidx_g_3_4_m0 + 137 * kcomp + i); 

            auto tg_zzz_xxyy_0 = braBuffer.data(pidx_g_3_4_m0 + 138 * kcomp + i); 

            auto tg_zzz_xxyz_0 = braBuffer.data(pidx_g_3_4_m0 + 139 * kcomp + i); 

            auto tg_zzz_xxzz_0 = braBuffer.data(pidx_g_3_4_m0 + 140 * kcomp + i); 

            auto tg_zzz_xyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 141 * kcomp + i); 

            auto tg_zzz_xyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 142 * kcomp + i); 

            auto tg_zzz_xyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 143 * kcomp + i); 

            auto tg_zzz_xzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 144 * kcomp + i); 

            auto tg_zzz_yyyy_0 = braBuffer.data(pidx_g_3_4_m0 + 145 * kcomp + i); 

            auto tg_zzz_yyyz_0 = braBuffer.data(pidx_g_3_4_m0 + 146 * kcomp + i); 

            auto tg_zzz_yyzz_0 = braBuffer.data(pidx_g_3_4_m0 + 147 * kcomp + i); 

            auto tg_zzz_yzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 148 * kcomp + i); 

            auto tg_zzz_zzzz_0 = braBuffer.data(pidx_g_3_4_m0 + 149 * kcomp + i); 

            // Batch of Integrals (75,150)

            #pragma omp simd aligned(tg_xzz_xxxx_0, tg_xzz_xxxy_0, tg_xzz_xxxz_0, tg_xzz_xxyy_0, \
                                     tg_xzz_xxyz_0, tg_xzz_xxzz_0, tg_xzz_xyyy_0, tg_xzz_xyyz_0, tg_xzz_xyzz_0, \
                                     tg_xzz_xzzz_0, tg_xzz_yyyy_0, tg_xzz_yyyz_0, tg_xzz_yyzz_0, tg_xzz_yzzz_0, \
                                     tg_xzz_zzzz_0, tg_yy_xxxx_0, tg_yy_xxxxy_0, tg_yy_xxxy_0, tg_yy_xxxyy_0, \
                                     tg_yy_xxxyz_0, tg_yy_xxxz_0, tg_yy_xxyy_0, tg_yy_xxyyy_0, tg_yy_xxyyz_0, \
                                     tg_yy_xxyz_0, tg_yy_xxyzz_0, tg_yy_xxzz_0, tg_yy_xyyy_0, tg_yy_xyyyy_0, \
                                     tg_yy_xyyyz_0, tg_yy_xyyz_0, tg_yy_xyyzz_0, tg_yy_xyzz_0, tg_yy_xyzzz_0, \
                                     tg_yy_xzzz_0, tg_yy_yyyy_0, tg_yy_yyyyy_0, tg_yy_yyyyz_0, tg_yy_yyyz_0, \
                                     tg_yy_yyyzz_0, tg_yy_yyzz_0, tg_yy_yyzzz_0, tg_yy_yzzz_0, tg_yy_yzzzz_0, \
                                     tg_yy_zzzz_0, tg_yyy_xxxx_0, tg_yyy_xxxy_0, tg_yyy_xxxz_0, tg_yyy_xxyy_0, \
                                     tg_yyy_xxyz_0, tg_yyy_xxzz_0, tg_yyy_xyyy_0, tg_yyy_xyyz_0, tg_yyy_xyzz_0, \
                                     tg_yyy_xzzz_0, tg_yyy_yyyy_0, tg_yyy_yyyz_0, tg_yyy_yyzz_0, tg_yyy_yzzz_0, \
                                     tg_yyy_zzzz_0, tg_yyz_xxxx_0, tg_yyz_xxxy_0, tg_yyz_xxxz_0, tg_yyz_xxyy_0, \
                                     tg_yyz_xxyz_0, tg_yyz_xxzz_0, tg_yyz_xyyy_0, tg_yyz_xyyz_0, tg_yyz_xyzz_0, \
                                     tg_yyz_xzzz_0, tg_yyz_yyyy_0, tg_yyz_yyyz_0, tg_yyz_yyzz_0, tg_yyz_yzzz_0, \
                                     tg_yyz_zzzz_0, tg_yz_xxxx_0, tg_yz_xxxxy_0, tg_yz_xxxy_0, tg_yz_xxxyy_0, \
                                     tg_yz_xxxyz_0, tg_yz_xxxz_0, tg_yz_xxyy_0, tg_yz_xxyyy_0, tg_yz_xxyyz_0, \
                                     tg_yz_xxyz_0, tg_yz_xxyzz_0, tg_yz_xxzz_0, tg_yz_xyyy_0, tg_yz_xyyyy_0, \
                                     tg_yz_xyyyz_0, tg_yz_xyyz_0, tg_yz_xyyzz_0, tg_yz_xyzz_0, tg_yz_xyzzz_0, \
                                     tg_yz_xzzz_0, tg_yz_yyyy_0, tg_yz_yyyyy_0, tg_yz_yyyyz_0, tg_yz_yyyz_0, \
                                     tg_yz_yyyzz_0, tg_yz_yyzz_0, tg_yz_yyzzz_0, tg_yz_yzzz_0, tg_yz_yzzzz_0, \
                                     tg_yz_zzzz_0, tg_yzz_xxxx_0, tg_yzz_xxxy_0, tg_yzz_xxxz_0, tg_yzz_xxyy_0, \
                                     tg_yzz_xxyz_0, tg_yzz_xxzz_0, tg_yzz_xyyy_0, tg_yzz_xyyz_0, tg_yzz_xyzz_0, \
                                     tg_yzz_xzzz_0, tg_yzz_yyyy_0, tg_yzz_yyyz_0, tg_yzz_yyzz_0, tg_yzz_yzzz_0, \
                                     tg_yzz_zzzz_0, tg_zz_xxxx_0, tg_zz_xxxxx_0, tg_zz_xxxxy_0, tg_zz_xxxxz_0, \
                                     tg_zz_xxxy_0, tg_zz_xxxyy_0, tg_zz_xxxyz_0, tg_zz_xxxz_0, tg_zz_xxxzz_0, \
                                     tg_zz_xxyy_0, tg_zz_xxyyy_0, tg_zz_xxyyz_0, tg_zz_xxyz_0, tg_zz_xxyzz_0, \
                                     tg_zz_xxzz_0, tg_zz_xxzzz_0, tg_zz_xyyy_0, tg_zz_xyyyy_0, tg_zz_xyyyz_0, \
                                     tg_zz_xyyz_0, tg_zz_xyyzz_0, tg_zz_xyzz_0, tg_zz_xyzzz_0, tg_zz_xzzz_0, \
                                     tg_zz_xzzzz_0, tg_zz_yyyy_0, tg_zz_yyyyy_0, tg_zz_yyyyz_0, tg_zz_yyyz_0, \
                                     tg_zz_yyyzz_0, tg_zz_yyzz_0, tg_zz_yyzzz_0, tg_zz_yzzz_0, tg_zz_yzzzz_0, \
                                     tg_zz_zzzz_0, tg_zz_zzzzz_0, tg_zzz_xxxx_0, tg_zzz_xxxy_0, tg_zzz_xxxz_0, \
                                     tg_zzz_xxyy_0, tg_zzz_xxyz_0, tg_zzz_xxzz_0, tg_zzz_xyyy_0, tg_zzz_xyyz_0, \
                                     tg_zzz_xyzz_0, tg_zzz_xzzz_0, tg_zzz_yyyy_0, tg_zzz_yyyz_0, tg_zzz_yyzz_0, \
                                     tg_zzz_yzzz_0, tg_zzz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xzz_xxxx_0[j] = -ab_x * tg_zz_xxxx_0[j] + tg_zz_xxxxx_0[j];

                tg_xzz_xxxy_0[j] = -ab_x * tg_zz_xxxy_0[j] + tg_zz_xxxxy_0[j];

                tg_xzz_xxxz_0[j] = -ab_x * tg_zz_xxxz_0[j] + tg_zz_xxxxz_0[j];

                tg_xzz_xxyy_0[j] = -ab_x * tg_zz_xxyy_0[j] + tg_zz_xxxyy_0[j];

                tg_xzz_xxyz_0[j] = -ab_x * tg_zz_xxyz_0[j] + tg_zz_xxxyz_0[j];

                tg_xzz_xxzz_0[j] = -ab_x * tg_zz_xxzz_0[j] + tg_zz_xxxzz_0[j];

                tg_xzz_xyyy_0[j] = -ab_x * tg_zz_xyyy_0[j] + tg_zz_xxyyy_0[j];

                tg_xzz_xyyz_0[j] = -ab_x * tg_zz_xyyz_0[j] + tg_zz_xxyyz_0[j];

                tg_xzz_xyzz_0[j] = -ab_x * tg_zz_xyzz_0[j] + tg_zz_xxyzz_0[j];

                tg_xzz_xzzz_0[j] = -ab_x * tg_zz_xzzz_0[j] + tg_zz_xxzzz_0[j];

                tg_xzz_yyyy_0[j] = -ab_x * tg_zz_yyyy_0[j] + tg_zz_xyyyy_0[j];

                tg_xzz_yyyz_0[j] = -ab_x * tg_zz_yyyz_0[j] + tg_zz_xyyyz_0[j];

                tg_xzz_yyzz_0[j] = -ab_x * tg_zz_yyzz_0[j] + tg_zz_xyyzz_0[j];

                tg_xzz_yzzz_0[j] = -ab_x * tg_zz_yzzz_0[j] + tg_zz_xyzzz_0[j];

                tg_xzz_zzzz_0[j] = -ab_x * tg_zz_zzzz_0[j] + tg_zz_xzzzz_0[j];

                tg_yyy_xxxx_0[j] = -ab_y * tg_yy_xxxx_0[j] + tg_yy_xxxxy_0[j];

                tg_yyy_xxxy_0[j] = -ab_y * tg_yy_xxxy_0[j] + tg_yy_xxxyy_0[j];

                tg_yyy_xxxz_0[j] = -ab_y * tg_yy_xxxz_0[j] + tg_yy_xxxyz_0[j];

                tg_yyy_xxyy_0[j] = -ab_y * tg_yy_xxyy_0[j] + tg_yy_xxyyy_0[j];

                tg_yyy_xxyz_0[j] = -ab_y * tg_yy_xxyz_0[j] + tg_yy_xxyyz_0[j];

                tg_yyy_xxzz_0[j] = -ab_y * tg_yy_xxzz_0[j] + tg_yy_xxyzz_0[j];

                tg_yyy_xyyy_0[j] = -ab_y * tg_yy_xyyy_0[j] + tg_yy_xyyyy_0[j];

                tg_yyy_xyyz_0[j] = -ab_y * tg_yy_xyyz_0[j] + tg_yy_xyyyz_0[j];

                tg_yyy_xyzz_0[j] = -ab_y * tg_yy_xyzz_0[j] + tg_yy_xyyzz_0[j];

                tg_yyy_xzzz_0[j] = -ab_y * tg_yy_xzzz_0[j] + tg_yy_xyzzz_0[j];

                tg_yyy_yyyy_0[j] = -ab_y * tg_yy_yyyy_0[j] + tg_yy_yyyyy_0[j];

                tg_yyy_yyyz_0[j] = -ab_y * tg_yy_yyyz_0[j] + tg_yy_yyyyz_0[j];

                tg_yyy_yyzz_0[j] = -ab_y * tg_yy_yyzz_0[j] + tg_yy_yyyzz_0[j];

                tg_yyy_yzzz_0[j] = -ab_y * tg_yy_yzzz_0[j] + tg_yy_yyzzz_0[j];

                tg_yyy_zzzz_0[j] = -ab_y * tg_yy_zzzz_0[j] + tg_yy_yzzzz_0[j];

                tg_yyz_xxxx_0[j] = -ab_y * tg_yz_xxxx_0[j] + tg_yz_xxxxy_0[j];

                tg_yyz_xxxy_0[j] = -ab_y * tg_yz_xxxy_0[j] + tg_yz_xxxyy_0[j];

                tg_yyz_xxxz_0[j] = -ab_y * tg_yz_xxxz_0[j] + tg_yz_xxxyz_0[j];

                tg_yyz_xxyy_0[j] = -ab_y * tg_yz_xxyy_0[j] + tg_yz_xxyyy_0[j];

                tg_yyz_xxyz_0[j] = -ab_y * tg_yz_xxyz_0[j] + tg_yz_xxyyz_0[j];

                tg_yyz_xxzz_0[j] = -ab_y * tg_yz_xxzz_0[j] + tg_yz_xxyzz_0[j];

                tg_yyz_xyyy_0[j] = -ab_y * tg_yz_xyyy_0[j] + tg_yz_xyyyy_0[j];

                tg_yyz_xyyz_0[j] = -ab_y * tg_yz_xyyz_0[j] + tg_yz_xyyyz_0[j];

                tg_yyz_xyzz_0[j] = -ab_y * tg_yz_xyzz_0[j] + tg_yz_xyyzz_0[j];

                tg_yyz_xzzz_0[j] = -ab_y * tg_yz_xzzz_0[j] + tg_yz_xyzzz_0[j];

                tg_yyz_yyyy_0[j] = -ab_y * tg_yz_yyyy_0[j] + tg_yz_yyyyy_0[j];

                tg_yyz_yyyz_0[j] = -ab_y * tg_yz_yyyz_0[j] + tg_yz_yyyyz_0[j];

                tg_yyz_yyzz_0[j] = -ab_y * tg_yz_yyzz_0[j] + tg_yz_yyyzz_0[j];

                tg_yyz_yzzz_0[j] = -ab_y * tg_yz_yzzz_0[j] + tg_yz_yyzzz_0[j];

                tg_yyz_zzzz_0[j] = -ab_y * tg_yz_zzzz_0[j] + tg_yz_yzzzz_0[j];

                tg_yzz_xxxx_0[j] = -ab_y * tg_zz_xxxx_0[j] + tg_zz_xxxxy_0[j];

                tg_yzz_xxxy_0[j] = -ab_y * tg_zz_xxxy_0[j] + tg_zz_xxxyy_0[j];

                tg_yzz_xxxz_0[j] = -ab_y * tg_zz_xxxz_0[j] + tg_zz_xxxyz_0[j];

                tg_yzz_xxyy_0[j] = -ab_y * tg_zz_xxyy_0[j] + tg_zz_xxyyy_0[j];

                tg_yzz_xxyz_0[j] = -ab_y * tg_zz_xxyz_0[j] + tg_zz_xxyyz_0[j];

                tg_yzz_xxzz_0[j] = -ab_y * tg_zz_xxzz_0[j] + tg_zz_xxyzz_0[j];

                tg_yzz_xyyy_0[j] = -ab_y * tg_zz_xyyy_0[j] + tg_zz_xyyyy_0[j];

                tg_yzz_xyyz_0[j] = -ab_y * tg_zz_xyyz_0[j] + tg_zz_xyyyz_0[j];

                tg_yzz_xyzz_0[j] = -ab_y * tg_zz_xyzz_0[j] + tg_zz_xyyzz_0[j];

                tg_yzz_xzzz_0[j] = -ab_y * tg_zz_xzzz_0[j] + tg_zz_xyzzz_0[j];

                tg_yzz_yyyy_0[j] = -ab_y * tg_zz_yyyy_0[j] + tg_zz_yyyyy_0[j];

                tg_yzz_yyyz_0[j] = -ab_y * tg_zz_yyyz_0[j] + tg_zz_yyyyz_0[j];

                tg_yzz_yyzz_0[j] = -ab_y * tg_zz_yyzz_0[j] + tg_zz_yyyzz_0[j];

                tg_yzz_yzzz_0[j] = -ab_y * tg_zz_yzzz_0[j] + tg_zz_yyzzz_0[j];

                tg_yzz_zzzz_0[j] = -ab_y * tg_zz_zzzz_0[j] + tg_zz_yzzzz_0[j];

                tg_zzz_xxxx_0[j] = -ab_z * tg_zz_xxxx_0[j] + tg_zz_xxxxz_0[j];

                tg_zzz_xxxy_0[j] = -ab_z * tg_zz_xxxy_0[j] + tg_zz_xxxyz_0[j];

                tg_zzz_xxxz_0[j] = -ab_z * tg_zz_xxxz_0[j] + tg_zz_xxxzz_0[j];

                tg_zzz_xxyy_0[j] = -ab_z * tg_zz_xxyy_0[j] + tg_zz_xxyyz_0[j];

                tg_zzz_xxyz_0[j] = -ab_z * tg_zz_xxyz_0[j] + tg_zz_xxyzz_0[j];

                tg_zzz_xxzz_0[j] = -ab_z * tg_zz_xxzz_0[j] + tg_zz_xxzzz_0[j];

                tg_zzz_xyyy_0[j] = -ab_z * tg_zz_xyyy_0[j] + tg_zz_xyyyz_0[j];

                tg_zzz_xyyz_0[j] = -ab_z * tg_zz_xyyz_0[j] + tg_zz_xyyzz_0[j];

                tg_zzz_xyzz_0[j] = -ab_z * tg_zz_xyzz_0[j] + tg_zz_xyzzz_0[j];

                tg_zzz_xzzz_0[j] = -ab_z * tg_zz_xzzz_0[j] + tg_zz_xzzzz_0[j];

                tg_zzz_yyyy_0[j] = -ab_z * tg_zz_yyyy_0[j] + tg_zz_yyyyz_0[j];

                tg_zzz_yyyz_0[j] = -ab_z * tg_zz_yyyz_0[j] + tg_zz_yyyzz_0[j];

                tg_zzz_yyzz_0[j] = -ab_z * tg_zz_yyzz_0[j] + tg_zz_yyzzz_0[j];

                tg_zzz_yzzz_0[j] = -ab_z * tg_zz_yzzz_0[j] + tg_zz_yzzzz_0[j];

                tg_zzz_zzzz_0[j] = -ab_z * tg_zz_zzzz_0[j] + tg_zz_zzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForFHXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        eribrrfunc::compElectronRepulsionForFHXY_0_70(braBuffer,
                                                      recursionMap,
                                                      abDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        eribrrfunc::compElectronRepulsionForFHXY_70_140(braBuffer,
                                                        recursionMap,
                                                        abDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 

        eribrrfunc::compElectronRepulsionForFHXY_140_210(braBuffer,
                                                         recursionMap,
                                                         abDistances,
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetContrPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForFHXY_0_70(      CMemBlock2D<double>& braBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& abDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,70)

        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_3_5_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_2_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            auto tg_yy_xxxxx_0 = braBuffer.data(pidx_g_2_5_m0 + 63 * kcomp + i); 

            auto tg_yy_xxxxy_0 = braBuffer.data(pidx_g_2_5_m0 + 64 * kcomp + i); 

            auto tg_yy_xxxxz_0 = braBuffer.data(pidx_g_2_5_m0 + 65 * kcomp + i); 

            auto tg_yy_xxxyy_0 = braBuffer.data(pidx_g_2_5_m0 + 66 * kcomp + i); 

            auto tg_yy_xxxyz_0 = braBuffer.data(pidx_g_2_5_m0 + 67 * kcomp + i); 

            auto tg_yy_xxxzz_0 = braBuffer.data(pidx_g_2_5_m0 + 68 * kcomp + i); 

            auto tg_yy_xxyyy_0 = braBuffer.data(pidx_g_2_5_m0 + 69 * kcomp + i); 

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

            auto tg_yy_xxxxxx_0 = braBuffer.data(pidx_g_2_6_m0 + 84 * kcomp + i); 

            auto tg_yy_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + 85 * kcomp + i); 

            auto tg_yy_xxxxxz_0 = braBuffer.data(pidx_g_2_6_m0 + 86 * kcomp + i); 

            auto tg_yy_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 87 * kcomp + i); 

            auto tg_yy_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 88 * kcomp + i); 

            auto tg_yy_xxxxzz_0 = braBuffer.data(pidx_g_2_6_m0 + 89 * kcomp + i); 

            auto tg_yy_xxxyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 90 * kcomp + i); 

            // set up pointers to integrals

            auto tg_xxx_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + i); 

            auto tg_xxx_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + kcomp + i); 

            auto tg_xxx_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 2 * kcomp + i); 

            auto tg_xxx_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 3 * kcomp + i); 

            auto tg_xxx_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 4 * kcomp + i); 

            auto tg_xxx_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 5 * kcomp + i); 

            auto tg_xxx_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 6 * kcomp + i); 

            auto tg_xxx_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 7 * kcomp + i); 

            auto tg_xxx_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 8 * kcomp + i); 

            auto tg_xxx_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 9 * kcomp + i); 

            auto tg_xxx_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 10 * kcomp + i); 

            auto tg_xxx_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 11 * kcomp + i); 

            auto tg_xxx_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 12 * kcomp + i); 

            auto tg_xxx_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 13 * kcomp + i); 

            auto tg_xxx_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 14 * kcomp + i); 

            auto tg_xxx_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 15 * kcomp + i); 

            auto tg_xxx_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 16 * kcomp + i); 

            auto tg_xxx_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 17 * kcomp + i); 

            auto tg_xxx_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 18 * kcomp + i); 

            auto tg_xxx_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 19 * kcomp + i); 

            auto tg_xxx_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 20 * kcomp + i); 

            auto tg_xxy_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 21 * kcomp + i); 

            auto tg_xxy_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 22 * kcomp + i); 

            auto tg_xxy_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 23 * kcomp + i); 

            auto tg_xxy_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 24 * kcomp + i); 

            auto tg_xxy_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 25 * kcomp + i); 

            auto tg_xxy_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 26 * kcomp + i); 

            auto tg_xxy_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 27 * kcomp + i); 

            auto tg_xxy_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 28 * kcomp + i); 

            auto tg_xxy_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 29 * kcomp + i); 

            auto tg_xxy_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 30 * kcomp + i); 

            auto tg_xxy_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 31 * kcomp + i); 

            auto tg_xxy_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 32 * kcomp + i); 

            auto tg_xxy_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 33 * kcomp + i); 

            auto tg_xxy_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 34 * kcomp + i); 

            auto tg_xxy_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 35 * kcomp + i); 

            auto tg_xxy_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 36 * kcomp + i); 

            auto tg_xxy_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 37 * kcomp + i); 

            auto tg_xxy_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 38 * kcomp + i); 

            auto tg_xxy_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 39 * kcomp + i); 

            auto tg_xxy_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 40 * kcomp + i); 

            auto tg_xxy_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 41 * kcomp + i); 

            auto tg_xxz_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 42 * kcomp + i); 

            auto tg_xxz_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 43 * kcomp + i); 

            auto tg_xxz_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 44 * kcomp + i); 

            auto tg_xxz_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 45 * kcomp + i); 

            auto tg_xxz_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 46 * kcomp + i); 

            auto tg_xxz_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 47 * kcomp + i); 

            auto tg_xxz_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 48 * kcomp + i); 

            auto tg_xxz_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 49 * kcomp + i); 

            auto tg_xxz_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 50 * kcomp + i); 

            auto tg_xxz_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 51 * kcomp + i); 

            auto tg_xxz_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 52 * kcomp + i); 

            auto tg_xxz_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 53 * kcomp + i); 

            auto tg_xxz_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 54 * kcomp + i); 

            auto tg_xxz_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 55 * kcomp + i); 

            auto tg_xxz_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 56 * kcomp + i); 

            auto tg_xxz_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 57 * kcomp + i); 

            auto tg_xxz_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 58 * kcomp + i); 

            auto tg_xxz_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 59 * kcomp + i); 

            auto tg_xxz_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 60 * kcomp + i); 

            auto tg_xxz_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 61 * kcomp + i); 

            auto tg_xxz_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 62 * kcomp + i); 

            auto tg_xyy_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 63 * kcomp + i); 

            auto tg_xyy_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 64 * kcomp + i); 

            auto tg_xyy_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 65 * kcomp + i); 

            auto tg_xyy_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 66 * kcomp + i); 

            auto tg_xyy_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 67 * kcomp + i); 

            auto tg_xyy_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 68 * kcomp + i); 

            auto tg_xyy_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 69 * kcomp + i); 

            // Batch of Integrals (0,70)

            #pragma omp simd aligned(tg_xx_xxxxx_0, tg_xx_xxxxxx_0, tg_xx_xxxxxy_0, tg_xx_xxxxxz_0, \
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
                tg_xxx_xxxxx_0[j] = -ab_x * tg_xx_xxxxx_0[j] + tg_xx_xxxxxx_0[j];

                tg_xxx_xxxxy_0[j] = -ab_x * tg_xx_xxxxy_0[j] + tg_xx_xxxxxy_0[j];

                tg_xxx_xxxxz_0[j] = -ab_x * tg_xx_xxxxz_0[j] + tg_xx_xxxxxz_0[j];

                tg_xxx_xxxyy_0[j] = -ab_x * tg_xx_xxxyy_0[j] + tg_xx_xxxxyy_0[j];

                tg_xxx_xxxyz_0[j] = -ab_x * tg_xx_xxxyz_0[j] + tg_xx_xxxxyz_0[j];

                tg_xxx_xxxzz_0[j] = -ab_x * tg_xx_xxxzz_0[j] + tg_xx_xxxxzz_0[j];

                tg_xxx_xxyyy_0[j] = -ab_x * tg_xx_xxyyy_0[j] + tg_xx_xxxyyy_0[j];

                tg_xxx_xxyyz_0[j] = -ab_x * tg_xx_xxyyz_0[j] + tg_xx_xxxyyz_0[j];

                tg_xxx_xxyzz_0[j] = -ab_x * tg_xx_xxyzz_0[j] + tg_xx_xxxyzz_0[j];

                tg_xxx_xxzzz_0[j] = -ab_x * tg_xx_xxzzz_0[j] + tg_xx_xxxzzz_0[j];

                tg_xxx_xyyyy_0[j] = -ab_x * tg_xx_xyyyy_0[j] + tg_xx_xxyyyy_0[j];

                tg_xxx_xyyyz_0[j] = -ab_x * tg_xx_xyyyz_0[j] + tg_xx_xxyyyz_0[j];

                tg_xxx_xyyzz_0[j] = -ab_x * tg_xx_xyyzz_0[j] + tg_xx_xxyyzz_0[j];

                tg_xxx_xyzzz_0[j] = -ab_x * tg_xx_xyzzz_0[j] + tg_xx_xxyzzz_0[j];

                tg_xxx_xzzzz_0[j] = -ab_x * tg_xx_xzzzz_0[j] + tg_xx_xxzzzz_0[j];

                tg_xxx_yyyyy_0[j] = -ab_x * tg_xx_yyyyy_0[j] + tg_xx_xyyyyy_0[j];

                tg_xxx_yyyyz_0[j] = -ab_x * tg_xx_yyyyz_0[j] + tg_xx_xyyyyz_0[j];

                tg_xxx_yyyzz_0[j] = -ab_x * tg_xx_yyyzz_0[j] + tg_xx_xyyyzz_0[j];

                tg_xxx_yyzzz_0[j] = -ab_x * tg_xx_yyzzz_0[j] + tg_xx_xyyzzz_0[j];

                tg_xxx_yzzzz_0[j] = -ab_x * tg_xx_yzzzz_0[j] + tg_xx_xyzzzz_0[j];

                tg_xxx_zzzzz_0[j] = -ab_x * tg_xx_zzzzz_0[j] + tg_xx_xzzzzz_0[j];

                tg_xxy_xxxxx_0[j] = -ab_x * tg_xy_xxxxx_0[j] + tg_xy_xxxxxx_0[j];

                tg_xxy_xxxxy_0[j] = -ab_x * tg_xy_xxxxy_0[j] + tg_xy_xxxxxy_0[j];

                tg_xxy_xxxxz_0[j] = -ab_x * tg_xy_xxxxz_0[j] + tg_xy_xxxxxz_0[j];

                tg_xxy_xxxyy_0[j] = -ab_x * tg_xy_xxxyy_0[j] + tg_xy_xxxxyy_0[j];

                tg_xxy_xxxyz_0[j] = -ab_x * tg_xy_xxxyz_0[j] + tg_xy_xxxxyz_0[j];

                tg_xxy_xxxzz_0[j] = -ab_x * tg_xy_xxxzz_0[j] + tg_xy_xxxxzz_0[j];

                tg_xxy_xxyyy_0[j] = -ab_x * tg_xy_xxyyy_0[j] + tg_xy_xxxyyy_0[j];

                tg_xxy_xxyyz_0[j] = -ab_x * tg_xy_xxyyz_0[j] + tg_xy_xxxyyz_0[j];

                tg_xxy_xxyzz_0[j] = -ab_x * tg_xy_xxyzz_0[j] + tg_xy_xxxyzz_0[j];

                tg_xxy_xxzzz_0[j] = -ab_x * tg_xy_xxzzz_0[j] + tg_xy_xxxzzz_0[j];

                tg_xxy_xyyyy_0[j] = -ab_x * tg_xy_xyyyy_0[j] + tg_xy_xxyyyy_0[j];

                tg_xxy_xyyyz_0[j] = -ab_x * tg_xy_xyyyz_0[j] + tg_xy_xxyyyz_0[j];

                tg_xxy_xyyzz_0[j] = -ab_x * tg_xy_xyyzz_0[j] + tg_xy_xxyyzz_0[j];

                tg_xxy_xyzzz_0[j] = -ab_x * tg_xy_xyzzz_0[j] + tg_xy_xxyzzz_0[j];

                tg_xxy_xzzzz_0[j] = -ab_x * tg_xy_xzzzz_0[j] + tg_xy_xxzzzz_0[j];

                tg_xxy_yyyyy_0[j] = -ab_x * tg_xy_yyyyy_0[j] + tg_xy_xyyyyy_0[j];

                tg_xxy_yyyyz_0[j] = -ab_x * tg_xy_yyyyz_0[j] + tg_xy_xyyyyz_0[j];

                tg_xxy_yyyzz_0[j] = -ab_x * tg_xy_yyyzz_0[j] + tg_xy_xyyyzz_0[j];

                tg_xxy_yyzzz_0[j] = -ab_x * tg_xy_yyzzz_0[j] + tg_xy_xyyzzz_0[j];

                tg_xxy_yzzzz_0[j] = -ab_x * tg_xy_yzzzz_0[j] + tg_xy_xyzzzz_0[j];

                tg_xxy_zzzzz_0[j] = -ab_x * tg_xy_zzzzz_0[j] + tg_xy_xzzzzz_0[j];

                tg_xxz_xxxxx_0[j] = -ab_x * tg_xz_xxxxx_0[j] + tg_xz_xxxxxx_0[j];

                tg_xxz_xxxxy_0[j] = -ab_x * tg_xz_xxxxy_0[j] + tg_xz_xxxxxy_0[j];

                tg_xxz_xxxxz_0[j] = -ab_x * tg_xz_xxxxz_0[j] + tg_xz_xxxxxz_0[j];

                tg_xxz_xxxyy_0[j] = -ab_x * tg_xz_xxxyy_0[j] + tg_xz_xxxxyy_0[j];

                tg_xxz_xxxyz_0[j] = -ab_x * tg_xz_xxxyz_0[j] + tg_xz_xxxxyz_0[j];

                tg_xxz_xxxzz_0[j] = -ab_x * tg_xz_xxxzz_0[j] + tg_xz_xxxxzz_0[j];

                tg_xxz_xxyyy_0[j] = -ab_x * tg_xz_xxyyy_0[j] + tg_xz_xxxyyy_0[j];

                tg_xxz_xxyyz_0[j] = -ab_x * tg_xz_xxyyz_0[j] + tg_xz_xxxyyz_0[j];

                tg_xxz_xxyzz_0[j] = -ab_x * tg_xz_xxyzz_0[j] + tg_xz_xxxyzz_0[j];

                tg_xxz_xxzzz_0[j] = -ab_x * tg_xz_xxzzz_0[j] + tg_xz_xxxzzz_0[j];

                tg_xxz_xyyyy_0[j] = -ab_x * tg_xz_xyyyy_0[j] + tg_xz_xxyyyy_0[j];

                tg_xxz_xyyyz_0[j] = -ab_x * tg_xz_xyyyz_0[j] + tg_xz_xxyyyz_0[j];

                tg_xxz_xyyzz_0[j] = -ab_x * tg_xz_xyyzz_0[j] + tg_xz_xxyyzz_0[j];

                tg_xxz_xyzzz_0[j] = -ab_x * tg_xz_xyzzz_0[j] + tg_xz_xxyzzz_0[j];

                tg_xxz_xzzzz_0[j] = -ab_x * tg_xz_xzzzz_0[j] + tg_xz_xxzzzz_0[j];

                tg_xxz_yyyyy_0[j] = -ab_x * tg_xz_yyyyy_0[j] + tg_xz_xyyyyy_0[j];

                tg_xxz_yyyyz_0[j] = -ab_x * tg_xz_yyyyz_0[j] + tg_xz_xyyyyz_0[j];

                tg_xxz_yyyzz_0[j] = -ab_x * tg_xz_yyyzz_0[j] + tg_xz_xyyyzz_0[j];

                tg_xxz_yyzzz_0[j] = -ab_x * tg_xz_yyzzz_0[j] + tg_xz_xyyzzz_0[j];

                tg_xxz_yzzzz_0[j] = -ab_x * tg_xz_yzzzz_0[j] + tg_xz_xyzzzz_0[j];

                tg_xxz_zzzzz_0[j] = -ab_x * tg_xz_zzzzz_0[j] + tg_xz_xzzzzz_0[j];

                tg_xyy_xxxxx_0[j] = -ab_x * tg_yy_xxxxx_0[j] + tg_yy_xxxxxx_0[j];

                tg_xyy_xxxxy_0[j] = -ab_x * tg_yy_xxxxy_0[j] + tg_yy_xxxxxy_0[j];

                tg_xyy_xxxxz_0[j] = -ab_x * tg_yy_xxxxz_0[j] + tg_yy_xxxxxz_0[j];

                tg_xyy_xxxyy_0[j] = -ab_x * tg_yy_xxxyy_0[j] + tg_yy_xxxxyy_0[j];

                tg_xyy_xxxyz_0[j] = -ab_x * tg_yy_xxxyz_0[j] + tg_yy_xxxxyz_0[j];

                tg_xyy_xxxzz_0[j] = -ab_x * tg_yy_xxxzz_0[j] + tg_yy_xxxxzz_0[j];

                tg_xyy_xxyyy_0[j] = -ab_x * tg_yy_xxyyy_0[j] + tg_yy_xxxyyy_0[j];
            }
        }
    }

    void
    compElectronRepulsionForFHXY_70_140(      CMemBlock2D<double>& braBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& abDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetContrPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (70,140)

        // set up distances R(AB) = A - B

        auto ab_x = (abDistances.data(0))[iContrPair];

        auto ab_y = (abDistances.data(1))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_3_5_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_2_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            auto tg_yy_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + 85 * kcomp + i); 

            auto tg_yy_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 87 * kcomp + i); 

            auto tg_yy_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 88 * kcomp + i); 

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

            // set up pointers to integrals

            auto tg_xyy_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 70 * kcomp + i); 

            auto tg_xyy_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 71 * kcomp + i); 

            auto tg_xyy_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 72 * kcomp + i); 

            auto tg_xyy_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 73 * kcomp + i); 

            auto tg_xyy_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 74 * kcomp + i); 

            auto tg_xyy_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 75 * kcomp + i); 

            auto tg_xyy_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 76 * kcomp + i); 

            auto tg_xyy_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 77 * kcomp + i); 

            auto tg_xyy_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 78 * kcomp + i); 

            auto tg_xyy_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 79 * kcomp + i); 

            auto tg_xyy_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 80 * kcomp + i); 

            auto tg_xyy_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 81 * kcomp + i); 

            auto tg_xyy_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 82 * kcomp + i); 

            auto tg_xyy_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 83 * kcomp + i); 

            auto tg_xyz_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 84 * kcomp + i); 

            auto tg_xyz_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 85 * kcomp + i); 

            auto tg_xyz_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 86 * kcomp + i); 

            auto tg_xyz_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 87 * kcomp + i); 

            auto tg_xyz_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 88 * kcomp + i); 

            auto tg_xyz_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 89 * kcomp + i); 

            auto tg_xyz_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 90 * kcomp + i); 

            auto tg_xyz_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 91 * kcomp + i); 

            auto tg_xyz_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 92 * kcomp + i); 

            auto tg_xyz_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 93 * kcomp + i); 

            auto tg_xyz_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 94 * kcomp + i); 

            auto tg_xyz_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 95 * kcomp + i); 

            auto tg_xyz_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 96 * kcomp + i); 

            auto tg_xyz_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 97 * kcomp + i); 

            auto tg_xyz_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 98 * kcomp + i); 

            auto tg_xyz_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 99 * kcomp + i); 

            auto tg_xyz_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 100 * kcomp + i); 

            auto tg_xyz_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 101 * kcomp + i); 

            auto tg_xyz_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 102 * kcomp + i); 

            auto tg_xyz_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 103 * kcomp + i); 

            auto tg_xyz_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 104 * kcomp + i); 

            auto tg_xzz_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 105 * kcomp + i); 

            auto tg_xzz_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 106 * kcomp + i); 

            auto tg_xzz_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 107 * kcomp + i); 

            auto tg_xzz_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 108 * kcomp + i); 

            auto tg_xzz_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 109 * kcomp + i); 

            auto tg_xzz_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 110 * kcomp + i); 

            auto tg_xzz_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 111 * kcomp + i); 

            auto tg_xzz_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 112 * kcomp + i); 

            auto tg_xzz_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 113 * kcomp + i); 

            auto tg_xzz_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 114 * kcomp + i); 

            auto tg_xzz_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 115 * kcomp + i); 

            auto tg_xzz_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 116 * kcomp + i); 

            auto tg_xzz_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 117 * kcomp + i); 

            auto tg_xzz_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 118 * kcomp + i); 

            auto tg_xzz_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 119 * kcomp + i); 

            auto tg_xzz_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 120 * kcomp + i); 

            auto tg_xzz_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 121 * kcomp + i); 

            auto tg_xzz_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 122 * kcomp + i); 

            auto tg_xzz_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 123 * kcomp + i); 

            auto tg_xzz_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 124 * kcomp + i); 

            auto tg_xzz_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 125 * kcomp + i); 

            auto tg_yyy_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 126 * kcomp + i); 

            auto tg_yyy_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 127 * kcomp + i); 

            auto tg_yyy_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 128 * kcomp + i); 

            auto tg_yyy_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 129 * kcomp + i); 

            auto tg_yyy_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 130 * kcomp + i); 

            auto tg_yyy_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 131 * kcomp + i); 

            auto tg_yyy_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 132 * kcomp + i); 

            auto tg_yyy_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 133 * kcomp + i); 

            auto tg_yyy_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 134 * kcomp + i); 

            auto tg_yyy_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 135 * kcomp + i); 

            auto tg_yyy_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 136 * kcomp + i); 

            auto tg_yyy_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 137 * kcomp + i); 

            auto tg_yyy_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 138 * kcomp + i); 

            auto tg_yyy_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 139 * kcomp + i); 

            // Batch of Integrals (70,140)

            #pragma omp simd aligned(tg_xyy_xxyyz_0, tg_xyy_xxyzz_0, tg_xyy_xxzzz_0, tg_xyy_xyyyy_0, \
                                     tg_xyy_xyyyz_0, tg_xyy_xyyzz_0, tg_xyy_xyzzz_0, tg_xyy_xzzzz_0, tg_xyy_yyyyy_0, \
                                     tg_xyy_yyyyz_0, tg_xyy_yyyzz_0, tg_xyy_yyzzz_0, tg_xyy_yzzzz_0, tg_xyy_zzzzz_0, \
                                     tg_xyz_xxxxx_0, tg_xyz_xxxxy_0, tg_xyz_xxxxz_0, tg_xyz_xxxyy_0, tg_xyz_xxxyz_0, \
                                     tg_xyz_xxxzz_0, tg_xyz_xxyyy_0, tg_xyz_xxyyz_0, tg_xyz_xxyzz_0, tg_xyz_xxzzz_0, \
                                     tg_xyz_xyyyy_0, tg_xyz_xyyyz_0, tg_xyz_xyyzz_0, tg_xyz_xyzzz_0, tg_xyz_xzzzz_0, \
                                     tg_xyz_yyyyy_0, tg_xyz_yyyyz_0, tg_xyz_yyyzz_0, tg_xyz_yyzzz_0, tg_xyz_yzzzz_0, \
                                     tg_xyz_zzzzz_0, tg_xzz_xxxxx_0, tg_xzz_xxxxy_0, tg_xzz_xxxxz_0, tg_xzz_xxxyy_0, \
                                     tg_xzz_xxxyz_0, tg_xzz_xxxzz_0, tg_xzz_xxyyy_0, tg_xzz_xxyyz_0, tg_xzz_xxyzz_0, \
                                     tg_xzz_xxzzz_0, tg_xzz_xyyyy_0, tg_xzz_xyyyz_0, tg_xzz_xyyzz_0, tg_xzz_xyzzz_0, \
                                     tg_xzz_xzzzz_0, tg_xzz_yyyyy_0, tg_xzz_yyyyz_0, tg_xzz_yyyzz_0, tg_xzz_yyzzz_0, \
                                     tg_xzz_yzzzz_0, tg_xzz_zzzzz_0, tg_yy_xxxxx_0, tg_yy_xxxxxy_0, tg_yy_xxxxy_0, \
                                     tg_yy_xxxxyy_0, tg_yy_xxxxyz_0, tg_yy_xxxxz_0, tg_yy_xxxyy_0, tg_yy_xxxyyy_0, \
                                     tg_yy_xxxyyz_0, tg_yy_xxxyz_0, tg_yy_xxxyzz_0, tg_yy_xxxzz_0, tg_yy_xxxzzz_0, \
                                     tg_yy_xxyyy_0, tg_yy_xxyyyy_0, tg_yy_xxyyyz_0, tg_yy_xxyyz_0, tg_yy_xxyyzz_0, \
                                     tg_yy_xxyzz_0, tg_yy_xxyzzz_0, tg_yy_xxzzz_0, tg_yy_xxzzzz_0, tg_yy_xyyyy_0, \
                                     tg_yy_xyyyyy_0, tg_yy_xyyyyz_0, tg_yy_xyyyz_0, tg_yy_xyyyzz_0, tg_yy_xyyzz_0, \
                                     tg_yy_xyyzzz_0, tg_yy_xyzzz_0, tg_yy_xyzzzz_0, tg_yy_xzzzz_0, tg_yy_xzzzzz_0, \
                                     tg_yy_yyyyy_0, tg_yy_yyyyz_0, tg_yy_yyyzz_0, tg_yy_yyzzz_0, tg_yy_yzzzz_0, \
                                     tg_yy_zzzzz_0, tg_yyy_xxxxx_0, tg_yyy_xxxxy_0, tg_yyy_xxxxz_0, tg_yyy_xxxyy_0, \
                                     tg_yyy_xxxyz_0, tg_yyy_xxxzz_0, tg_yyy_xxyyy_0, tg_yyy_xxyyz_0, tg_yyy_xxyzz_0, \
                                     tg_yyy_xxzzz_0, tg_yyy_xyyyy_0, tg_yyy_xyyyz_0, tg_yyy_xyyzz_0, tg_yyy_xyzzz_0, \
                                     tg_yz_xxxxx_0, tg_yz_xxxxxx_0, tg_yz_xxxxxy_0, tg_yz_xxxxxz_0, tg_yz_xxxxy_0, \
                                     tg_yz_xxxxyy_0, tg_yz_xxxxyz_0, tg_yz_xxxxz_0, tg_yz_xxxxzz_0, tg_yz_xxxyy_0, \
                                     tg_yz_xxxyyy_0, tg_yz_xxxyyz_0, tg_yz_xxxyz_0, tg_yz_xxxyzz_0, tg_yz_xxxzz_0, \
                                     tg_yz_xxxzzz_0, tg_yz_xxyyy_0, tg_yz_xxyyyy_0, tg_yz_xxyyyz_0, tg_yz_xxyyz_0, \
                                     tg_yz_xxyyzz_0, tg_yz_xxyzz_0, tg_yz_xxyzzz_0, tg_yz_xxzzz_0, tg_yz_xxzzzz_0, \
                                     tg_yz_xyyyy_0, tg_yz_xyyyyy_0, tg_yz_xyyyyz_0, tg_yz_xyyyz_0, tg_yz_xyyyzz_0, \
                                     tg_yz_xyyzz_0, tg_yz_xyyzzz_0, tg_yz_xyzzz_0, tg_yz_xyzzzz_0, tg_yz_xzzzz_0, \
                                     tg_yz_xzzzzz_0, tg_yz_yyyyy_0, tg_yz_yyyyz_0, tg_yz_yyyzz_0, tg_yz_yyzzz_0, \
                                     tg_yz_yzzzz_0, tg_yz_zzzzz_0, tg_zz_xxxxx_0, tg_zz_xxxxxx_0, tg_zz_xxxxxy_0, \
                                     tg_zz_xxxxxz_0, tg_zz_xxxxy_0, tg_zz_xxxxyy_0, tg_zz_xxxxyz_0, tg_zz_xxxxz_0, \
                                     tg_zz_xxxxzz_0, tg_zz_xxxyy_0, tg_zz_xxxyyy_0, tg_zz_xxxyyz_0, tg_zz_xxxyz_0, \
                                     tg_zz_xxxyzz_0, tg_zz_xxxzz_0, tg_zz_xxxzzz_0, tg_zz_xxyyy_0, tg_zz_xxyyyy_0, \
                                     tg_zz_xxyyyz_0, tg_zz_xxyyz_0, tg_zz_xxyyzz_0, tg_zz_xxyzz_0, tg_zz_xxyzzz_0, \
                                     tg_zz_xxzzz_0, tg_zz_xxzzzz_0, tg_zz_xyyyy_0, tg_zz_xyyyyy_0, tg_zz_xyyyyz_0, \
                                     tg_zz_xyyyz_0, tg_zz_xyyyzz_0, tg_zz_xyyzz_0, tg_zz_xyyzzz_0, tg_zz_xyzzz_0, \
                                     tg_zz_xyzzzz_0, tg_zz_xzzzz_0, tg_zz_xzzzzz_0, tg_zz_yyyyy_0, tg_zz_yyyyz_0, \
                                     tg_zz_yyyzz_0, tg_zz_yyzzz_0, tg_zz_yzzzz_0, tg_zz_zzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xyy_xxyyz_0[j] = -ab_x * tg_yy_xxyyz_0[j] + tg_yy_xxxyyz_0[j];

                tg_xyy_xxyzz_0[j] = -ab_x * tg_yy_xxyzz_0[j] + tg_yy_xxxyzz_0[j];

                tg_xyy_xxzzz_0[j] = -ab_x * tg_yy_xxzzz_0[j] + tg_yy_xxxzzz_0[j];

                tg_xyy_xyyyy_0[j] = -ab_x * tg_yy_xyyyy_0[j] + tg_yy_xxyyyy_0[j];

                tg_xyy_xyyyz_0[j] = -ab_x * tg_yy_xyyyz_0[j] + tg_yy_xxyyyz_0[j];

                tg_xyy_xyyzz_0[j] = -ab_x * tg_yy_xyyzz_0[j] + tg_yy_xxyyzz_0[j];

                tg_xyy_xyzzz_0[j] = -ab_x * tg_yy_xyzzz_0[j] + tg_yy_xxyzzz_0[j];

                tg_xyy_xzzzz_0[j] = -ab_x * tg_yy_xzzzz_0[j] + tg_yy_xxzzzz_0[j];

                tg_xyy_yyyyy_0[j] = -ab_x * tg_yy_yyyyy_0[j] + tg_yy_xyyyyy_0[j];

                tg_xyy_yyyyz_0[j] = -ab_x * tg_yy_yyyyz_0[j] + tg_yy_xyyyyz_0[j];

                tg_xyy_yyyzz_0[j] = -ab_x * tg_yy_yyyzz_0[j] + tg_yy_xyyyzz_0[j];

                tg_xyy_yyzzz_0[j] = -ab_x * tg_yy_yyzzz_0[j] + tg_yy_xyyzzz_0[j];

                tg_xyy_yzzzz_0[j] = -ab_x * tg_yy_yzzzz_0[j] + tg_yy_xyzzzz_0[j];

                tg_xyy_zzzzz_0[j] = -ab_x * tg_yy_zzzzz_0[j] + tg_yy_xzzzzz_0[j];

                tg_xyz_xxxxx_0[j] = -ab_x * tg_yz_xxxxx_0[j] + tg_yz_xxxxxx_0[j];

                tg_xyz_xxxxy_0[j] = -ab_x * tg_yz_xxxxy_0[j] + tg_yz_xxxxxy_0[j];

                tg_xyz_xxxxz_0[j] = -ab_x * tg_yz_xxxxz_0[j] + tg_yz_xxxxxz_0[j];

                tg_xyz_xxxyy_0[j] = -ab_x * tg_yz_xxxyy_0[j] + tg_yz_xxxxyy_0[j];

                tg_xyz_xxxyz_0[j] = -ab_x * tg_yz_xxxyz_0[j] + tg_yz_xxxxyz_0[j];

                tg_xyz_xxxzz_0[j] = -ab_x * tg_yz_xxxzz_0[j] + tg_yz_xxxxzz_0[j];

                tg_xyz_xxyyy_0[j] = -ab_x * tg_yz_xxyyy_0[j] + tg_yz_xxxyyy_0[j];

                tg_xyz_xxyyz_0[j] = -ab_x * tg_yz_xxyyz_0[j] + tg_yz_xxxyyz_0[j];

                tg_xyz_xxyzz_0[j] = -ab_x * tg_yz_xxyzz_0[j] + tg_yz_xxxyzz_0[j];

                tg_xyz_xxzzz_0[j] = -ab_x * tg_yz_xxzzz_0[j] + tg_yz_xxxzzz_0[j];

                tg_xyz_xyyyy_0[j] = -ab_x * tg_yz_xyyyy_0[j] + tg_yz_xxyyyy_0[j];

                tg_xyz_xyyyz_0[j] = -ab_x * tg_yz_xyyyz_0[j] + tg_yz_xxyyyz_0[j];

                tg_xyz_xyyzz_0[j] = -ab_x * tg_yz_xyyzz_0[j] + tg_yz_xxyyzz_0[j];

                tg_xyz_xyzzz_0[j] = -ab_x * tg_yz_xyzzz_0[j] + tg_yz_xxyzzz_0[j];

                tg_xyz_xzzzz_0[j] = -ab_x * tg_yz_xzzzz_0[j] + tg_yz_xxzzzz_0[j];

                tg_xyz_yyyyy_0[j] = -ab_x * tg_yz_yyyyy_0[j] + tg_yz_xyyyyy_0[j];

                tg_xyz_yyyyz_0[j] = -ab_x * tg_yz_yyyyz_0[j] + tg_yz_xyyyyz_0[j];

                tg_xyz_yyyzz_0[j] = -ab_x * tg_yz_yyyzz_0[j] + tg_yz_xyyyzz_0[j];

                tg_xyz_yyzzz_0[j] = -ab_x * tg_yz_yyzzz_0[j] + tg_yz_xyyzzz_0[j];

                tg_xyz_yzzzz_0[j] = -ab_x * tg_yz_yzzzz_0[j] + tg_yz_xyzzzz_0[j];

                tg_xyz_zzzzz_0[j] = -ab_x * tg_yz_zzzzz_0[j] + tg_yz_xzzzzz_0[j];

                tg_xzz_xxxxx_0[j] = -ab_x * tg_zz_xxxxx_0[j] + tg_zz_xxxxxx_0[j];

                tg_xzz_xxxxy_0[j] = -ab_x * tg_zz_xxxxy_0[j] + tg_zz_xxxxxy_0[j];

                tg_xzz_xxxxz_0[j] = -ab_x * tg_zz_xxxxz_0[j] + tg_zz_xxxxxz_0[j];

                tg_xzz_xxxyy_0[j] = -ab_x * tg_zz_xxxyy_0[j] + tg_zz_xxxxyy_0[j];

                tg_xzz_xxxyz_0[j] = -ab_x * tg_zz_xxxyz_0[j] + tg_zz_xxxxyz_0[j];

                tg_xzz_xxxzz_0[j] = -ab_x * tg_zz_xxxzz_0[j] + tg_zz_xxxxzz_0[j];

                tg_xzz_xxyyy_0[j] = -ab_x * tg_zz_xxyyy_0[j] + tg_zz_xxxyyy_0[j];

                tg_xzz_xxyyz_0[j] = -ab_x * tg_zz_xxyyz_0[j] + tg_zz_xxxyyz_0[j];

                tg_xzz_xxyzz_0[j] = -ab_x * tg_zz_xxyzz_0[j] + tg_zz_xxxyzz_0[j];

                tg_xzz_xxzzz_0[j] = -ab_x * tg_zz_xxzzz_0[j] + tg_zz_xxxzzz_0[j];

                tg_xzz_xyyyy_0[j] = -ab_x * tg_zz_xyyyy_0[j] + tg_zz_xxyyyy_0[j];

                tg_xzz_xyyyz_0[j] = -ab_x * tg_zz_xyyyz_0[j] + tg_zz_xxyyyz_0[j];

                tg_xzz_xyyzz_0[j] = -ab_x * tg_zz_xyyzz_0[j] + tg_zz_xxyyzz_0[j];

                tg_xzz_xyzzz_0[j] = -ab_x * tg_zz_xyzzz_0[j] + tg_zz_xxyzzz_0[j];

                tg_xzz_xzzzz_0[j] = -ab_x * tg_zz_xzzzz_0[j] + tg_zz_xxzzzz_0[j];

                tg_xzz_yyyyy_0[j] = -ab_x * tg_zz_yyyyy_0[j] + tg_zz_xyyyyy_0[j];

                tg_xzz_yyyyz_0[j] = -ab_x * tg_zz_yyyyz_0[j] + tg_zz_xyyyyz_0[j];

                tg_xzz_yyyzz_0[j] = -ab_x * tg_zz_yyyzz_0[j] + tg_zz_xyyyzz_0[j];

                tg_xzz_yyzzz_0[j] = -ab_x * tg_zz_yyzzz_0[j] + tg_zz_xyyzzz_0[j];

                tg_xzz_yzzzz_0[j] = -ab_x * tg_zz_yzzzz_0[j] + tg_zz_xyzzzz_0[j];

                tg_xzz_zzzzz_0[j] = -ab_x * tg_zz_zzzzz_0[j] + tg_zz_xzzzzz_0[j];

                tg_yyy_xxxxx_0[j] = -ab_y * tg_yy_xxxxx_0[j] + tg_yy_xxxxxy_0[j];

                tg_yyy_xxxxy_0[j] = -ab_y * tg_yy_xxxxy_0[j] + tg_yy_xxxxyy_0[j];

                tg_yyy_xxxxz_0[j] = -ab_y * tg_yy_xxxxz_0[j] + tg_yy_xxxxyz_0[j];

                tg_yyy_xxxyy_0[j] = -ab_y * tg_yy_xxxyy_0[j] + tg_yy_xxxyyy_0[j];

                tg_yyy_xxxyz_0[j] = -ab_y * tg_yy_xxxyz_0[j] + tg_yy_xxxyyz_0[j];

                tg_yyy_xxxzz_0[j] = -ab_y * tg_yy_xxxzz_0[j] + tg_yy_xxxyzz_0[j];

                tg_yyy_xxyyy_0[j] = -ab_y * tg_yy_xxyyy_0[j] + tg_yy_xxyyyy_0[j];

                tg_yyy_xxyyz_0[j] = -ab_y * tg_yy_xxyyz_0[j] + tg_yy_xxyyyz_0[j];

                tg_yyy_xxyzz_0[j] = -ab_y * tg_yy_xxyzz_0[j] + tg_yy_xxyyzz_0[j];

                tg_yyy_xxzzz_0[j] = -ab_y * tg_yy_xxzzz_0[j] + tg_yy_xxyzzz_0[j];

                tg_yyy_xyyyy_0[j] = -ab_y * tg_yy_xyyyy_0[j] + tg_yy_xyyyyy_0[j];

                tg_yyy_xyyyz_0[j] = -ab_y * tg_yy_xyyyz_0[j] + tg_yy_xyyyyz_0[j];

                tg_yyy_xyyzz_0[j] = -ab_y * tg_yy_xyyzz_0[j] + tg_yy_xyyyzz_0[j];

                tg_yyy_xyzzz_0[j] = -ab_y * tg_yy_xyzzz_0[j] + tg_yy_xyyzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForFHXY_140_210(      CMemBlock2D<double>& braBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& abDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetContrPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (140,210)

        // set up distances R(AB) = A - B

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_3_5_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_2_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {2, 6, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            auto tg_yy_xyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 103 * kcomp + i); 

            auto tg_yy_yyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 105 * kcomp + i); 

            auto tg_yy_yyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 106 * kcomp + i); 

            auto tg_yy_yyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 107 * kcomp + i); 

            auto tg_yy_yyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 108 * kcomp + i); 

            auto tg_yy_yyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 109 * kcomp + i); 

            auto tg_yy_yzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 110 * kcomp + i); 

            auto tg_yz_xxxxxy_0 = braBuffer.data(pidx_g_2_6_m0 + 113 * kcomp + i); 

            auto tg_yz_xxxxyy_0 = braBuffer.data(pidx_g_2_6_m0 + 115 * kcomp + i); 

            auto tg_yz_xxxxyz_0 = braBuffer.data(pidx_g_2_6_m0 + 116 * kcomp + i); 

            auto tg_yz_xxxyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 118 * kcomp + i); 

            auto tg_yz_xxxyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 119 * kcomp + i); 

            auto tg_yz_xxxyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 120 * kcomp + i); 

            auto tg_yz_xxyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 122 * kcomp + i); 

            auto tg_yz_xxyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 123 * kcomp + i); 

            auto tg_yz_xxyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 124 * kcomp + i); 

            auto tg_yz_xxyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 125 * kcomp + i); 

            auto tg_yz_xyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 127 * kcomp + i); 

            auto tg_yz_xyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 128 * kcomp + i); 

            auto tg_yz_xyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 129 * kcomp + i); 

            auto tg_yz_xyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 130 * kcomp + i); 

            auto tg_yz_xyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 131 * kcomp + i); 

            auto tg_yz_yyyyyy_0 = braBuffer.data(pidx_g_2_6_m0 + 133 * kcomp + i); 

            auto tg_yz_yyyyyz_0 = braBuffer.data(pidx_g_2_6_m0 + 134 * kcomp + i); 

            auto tg_yz_yyyyzz_0 = braBuffer.data(pidx_g_2_6_m0 + 135 * kcomp + i); 

            auto tg_yz_yyyzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 136 * kcomp + i); 

            auto tg_yz_yyzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 137 * kcomp + i); 

            auto tg_yz_yzzzzz_0 = braBuffer.data(pidx_g_2_6_m0 + 138 * kcomp + i); 

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

            // set up pointers to integrals

            auto tg_yyy_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 140 * kcomp + i); 

            auto tg_yyy_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 141 * kcomp + i); 

            auto tg_yyy_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 142 * kcomp + i); 

            auto tg_yyy_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 143 * kcomp + i); 

            auto tg_yyy_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 144 * kcomp + i); 

            auto tg_yyy_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 145 * kcomp + i); 

            auto tg_yyy_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 146 * kcomp + i); 

            auto tg_yyz_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 147 * kcomp + i); 

            auto tg_yyz_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 148 * kcomp + i); 

            auto tg_yyz_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 149 * kcomp + i); 

            auto tg_yyz_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 150 * kcomp + i); 

            auto tg_yyz_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 151 * kcomp + i); 

            auto tg_yyz_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 152 * kcomp + i); 

            auto tg_yyz_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 153 * kcomp + i); 

            auto tg_yyz_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 154 * kcomp + i); 

            auto tg_yyz_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 155 * kcomp + i); 

            auto tg_yyz_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 156 * kcomp + i); 

            auto tg_yyz_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 157 * kcomp + i); 

            auto tg_yyz_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 158 * kcomp + i); 

            auto tg_yyz_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 159 * kcomp + i); 

            auto tg_yyz_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 160 * kcomp + i); 

            auto tg_yyz_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 161 * kcomp + i); 

            auto tg_yyz_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 162 * kcomp + i); 

            auto tg_yyz_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 163 * kcomp + i); 

            auto tg_yyz_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 164 * kcomp + i); 

            auto tg_yyz_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 165 * kcomp + i); 

            auto tg_yyz_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 166 * kcomp + i); 

            auto tg_yyz_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 167 * kcomp + i); 

            auto tg_yzz_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 168 * kcomp + i); 

            auto tg_yzz_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 169 * kcomp + i); 

            auto tg_yzz_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 170 * kcomp + i); 

            auto tg_yzz_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 171 * kcomp + i); 

            auto tg_yzz_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 172 * kcomp + i); 

            auto tg_yzz_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 173 * kcomp + i); 

            auto tg_yzz_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 174 * kcomp + i); 

            auto tg_yzz_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 175 * kcomp + i); 

            auto tg_yzz_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 176 * kcomp + i); 

            auto tg_yzz_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 177 * kcomp + i); 

            auto tg_yzz_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 178 * kcomp + i); 

            auto tg_yzz_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 179 * kcomp + i); 

            auto tg_yzz_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 180 * kcomp + i); 

            auto tg_yzz_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 181 * kcomp + i); 

            auto tg_yzz_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 182 * kcomp + i); 

            auto tg_yzz_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 183 * kcomp + i); 

            auto tg_yzz_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 184 * kcomp + i); 

            auto tg_yzz_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 185 * kcomp + i); 

            auto tg_yzz_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 186 * kcomp + i); 

            auto tg_yzz_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 187 * kcomp + i); 

            auto tg_yzz_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 188 * kcomp + i); 

            auto tg_zzz_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 189 * kcomp + i); 

            auto tg_zzz_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 190 * kcomp + i); 

            auto tg_zzz_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 191 * kcomp + i); 

            auto tg_zzz_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 192 * kcomp + i); 

            auto tg_zzz_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 193 * kcomp + i); 

            auto tg_zzz_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 194 * kcomp + i); 

            auto tg_zzz_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 195 * kcomp + i); 

            auto tg_zzz_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 196 * kcomp + i); 

            auto tg_zzz_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 197 * kcomp + i); 

            auto tg_zzz_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 198 * kcomp + i); 

            auto tg_zzz_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 199 * kcomp + i); 

            auto tg_zzz_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 200 * kcomp + i); 

            auto tg_zzz_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 201 * kcomp + i); 

            auto tg_zzz_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 202 * kcomp + i); 

            auto tg_zzz_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 203 * kcomp + i); 

            auto tg_zzz_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 204 * kcomp + i); 

            auto tg_zzz_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 205 * kcomp + i); 

            auto tg_zzz_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 206 * kcomp + i); 

            auto tg_zzz_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 207 * kcomp + i); 

            auto tg_zzz_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 208 * kcomp + i); 

            auto tg_zzz_zzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 209 * kcomp + i); 

            // Batch of Integrals (140,210)

            #pragma omp simd aligned(tg_yy_xyzzzz_0, tg_yy_xzzzz_0, tg_yy_yyyyy_0, tg_yy_yyyyyy_0, \
                                     tg_yy_yyyyyz_0, tg_yy_yyyyz_0, tg_yy_yyyyzz_0, tg_yy_yyyzz_0, tg_yy_yyyzzz_0, \
                                     tg_yy_yyzzz_0, tg_yy_yyzzzz_0, tg_yy_yzzzz_0, tg_yy_yzzzzz_0, tg_yy_zzzzz_0, \
                                     tg_yyy_xzzzz_0, tg_yyy_yyyyy_0, tg_yyy_yyyyz_0, tg_yyy_yyyzz_0, tg_yyy_yyzzz_0, \
                                     tg_yyy_yzzzz_0, tg_yyy_zzzzz_0, tg_yyz_xxxxx_0, tg_yyz_xxxxy_0, tg_yyz_xxxxz_0, \
                                     tg_yyz_xxxyy_0, tg_yyz_xxxyz_0, tg_yyz_xxxzz_0, tg_yyz_xxyyy_0, tg_yyz_xxyyz_0, \
                                     tg_yyz_xxyzz_0, tg_yyz_xxzzz_0, tg_yyz_xyyyy_0, tg_yyz_xyyyz_0, tg_yyz_xyyzz_0, \
                                     tg_yyz_xyzzz_0, tg_yyz_xzzzz_0, tg_yyz_yyyyy_0, tg_yyz_yyyyz_0, tg_yyz_yyyzz_0, \
                                     tg_yyz_yyzzz_0, tg_yyz_yzzzz_0, tg_yyz_zzzzz_0, tg_yz_xxxxx_0, tg_yz_xxxxxy_0, \
                                     tg_yz_xxxxy_0, tg_yz_xxxxyy_0, tg_yz_xxxxyz_0, tg_yz_xxxxz_0, tg_yz_xxxyy_0, \
                                     tg_yz_xxxyyy_0, tg_yz_xxxyyz_0, tg_yz_xxxyz_0, tg_yz_xxxyzz_0, tg_yz_xxxzz_0, \
                                     tg_yz_xxyyy_0, tg_yz_xxyyyy_0, tg_yz_xxyyyz_0, tg_yz_xxyyz_0, tg_yz_xxyyzz_0, \
                                     tg_yz_xxyzz_0, tg_yz_xxyzzz_0, tg_yz_xxzzz_0, tg_yz_xyyyy_0, tg_yz_xyyyyy_0, \
                                     tg_yz_xyyyyz_0, tg_yz_xyyyz_0, tg_yz_xyyyzz_0, tg_yz_xyyzz_0, tg_yz_xyyzzz_0, \
                                     tg_yz_xyzzz_0, tg_yz_xyzzzz_0, tg_yz_xzzzz_0, tg_yz_yyyyy_0, tg_yz_yyyyyy_0, \
                                     tg_yz_yyyyyz_0, tg_yz_yyyyz_0, tg_yz_yyyyzz_0, tg_yz_yyyzz_0, tg_yz_yyyzzz_0, \
                                     tg_yz_yyzzz_0, tg_yz_yyzzzz_0, tg_yz_yzzzz_0, tg_yz_yzzzzz_0, tg_yz_zzzzz_0, \
                                     tg_yzz_xxxxx_0, tg_yzz_xxxxy_0, tg_yzz_xxxxz_0, tg_yzz_xxxyy_0, tg_yzz_xxxyz_0, \
                                     tg_yzz_xxxzz_0, tg_yzz_xxyyy_0, tg_yzz_xxyyz_0, tg_yzz_xxyzz_0, tg_yzz_xxzzz_0, \
                                     tg_yzz_xyyyy_0, tg_yzz_xyyyz_0, tg_yzz_xyyzz_0, tg_yzz_xyzzz_0, tg_yzz_xzzzz_0, \
                                     tg_yzz_yyyyy_0, tg_yzz_yyyyz_0, tg_yzz_yyyzz_0, tg_yzz_yyzzz_0, tg_yzz_yzzzz_0, \
                                     tg_yzz_zzzzz_0, tg_zz_xxxxx_0, tg_zz_xxxxxy_0, tg_zz_xxxxxz_0, tg_zz_xxxxy_0, \
                                     tg_zz_xxxxyy_0, tg_zz_xxxxyz_0, tg_zz_xxxxz_0, tg_zz_xxxxzz_0, tg_zz_xxxyy_0, \
                                     tg_zz_xxxyyy_0, tg_zz_xxxyyz_0, tg_zz_xxxyz_0, tg_zz_xxxyzz_0, tg_zz_xxxzz_0, \
                                     tg_zz_xxxzzz_0, tg_zz_xxyyy_0, tg_zz_xxyyyy_0, tg_zz_xxyyyz_0, tg_zz_xxyyz_0, \
                                     tg_zz_xxyyzz_0, tg_zz_xxyzz_0, tg_zz_xxyzzz_0, tg_zz_xxzzz_0, tg_zz_xxzzzz_0, \
                                     tg_zz_xyyyy_0, tg_zz_xyyyyy_0, tg_zz_xyyyyz_0, tg_zz_xyyyz_0, tg_zz_xyyyzz_0, \
                                     tg_zz_xyyzz_0, tg_zz_xyyzzz_0, tg_zz_xyzzz_0, tg_zz_xyzzzz_0, tg_zz_xzzzz_0, \
                                     tg_zz_xzzzzz_0, tg_zz_yyyyy_0, tg_zz_yyyyyy_0, tg_zz_yyyyyz_0, tg_zz_yyyyz_0, \
                                     tg_zz_yyyyzz_0, tg_zz_yyyzz_0, tg_zz_yyyzzz_0, tg_zz_yyzzz_0, tg_zz_yyzzzz_0, \
                                     tg_zz_yzzzz_0, tg_zz_yzzzzz_0, tg_zz_zzzzz_0, tg_zz_zzzzzz_0, tg_zzz_xxxxx_0, \
                                     tg_zzz_xxxxy_0, tg_zzz_xxxxz_0, tg_zzz_xxxyy_0, tg_zzz_xxxyz_0, tg_zzz_xxxzz_0, \
                                     tg_zzz_xxyyy_0, tg_zzz_xxyyz_0, tg_zzz_xxyzz_0, tg_zzz_xxzzz_0, tg_zzz_xyyyy_0, \
                                     tg_zzz_xyyyz_0, tg_zzz_xyyzz_0, tg_zzz_xyzzz_0, tg_zzz_xzzzz_0, tg_zzz_yyyyy_0, \
                                     tg_zzz_yyyyz_0, tg_zzz_yyyzz_0, tg_zzz_yyzzz_0, tg_zzz_yzzzz_0, tg_zzz_zzzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_yyy_xzzzz_0[j] = -ab_y * tg_yy_xzzzz_0[j] + tg_yy_xyzzzz_0[j];

                tg_yyy_yyyyy_0[j] = -ab_y * tg_yy_yyyyy_0[j] + tg_yy_yyyyyy_0[j];

                tg_yyy_yyyyz_0[j] = -ab_y * tg_yy_yyyyz_0[j] + tg_yy_yyyyyz_0[j];

                tg_yyy_yyyzz_0[j] = -ab_y * tg_yy_yyyzz_0[j] + tg_yy_yyyyzz_0[j];

                tg_yyy_yyzzz_0[j] = -ab_y * tg_yy_yyzzz_0[j] + tg_yy_yyyzzz_0[j];

                tg_yyy_yzzzz_0[j] = -ab_y * tg_yy_yzzzz_0[j] + tg_yy_yyzzzz_0[j];

                tg_yyy_zzzzz_0[j] = -ab_y * tg_yy_zzzzz_0[j] + tg_yy_yzzzzz_0[j];

                tg_yyz_xxxxx_0[j] = -ab_y * tg_yz_xxxxx_0[j] + tg_yz_xxxxxy_0[j];

                tg_yyz_xxxxy_0[j] = -ab_y * tg_yz_xxxxy_0[j] + tg_yz_xxxxyy_0[j];

                tg_yyz_xxxxz_0[j] = -ab_y * tg_yz_xxxxz_0[j] + tg_yz_xxxxyz_0[j];

                tg_yyz_xxxyy_0[j] = -ab_y * tg_yz_xxxyy_0[j] + tg_yz_xxxyyy_0[j];

                tg_yyz_xxxyz_0[j] = -ab_y * tg_yz_xxxyz_0[j] + tg_yz_xxxyyz_0[j];

                tg_yyz_xxxzz_0[j] = -ab_y * tg_yz_xxxzz_0[j] + tg_yz_xxxyzz_0[j];

                tg_yyz_xxyyy_0[j] = -ab_y * tg_yz_xxyyy_0[j] + tg_yz_xxyyyy_0[j];

                tg_yyz_xxyyz_0[j] = -ab_y * tg_yz_xxyyz_0[j] + tg_yz_xxyyyz_0[j];

                tg_yyz_xxyzz_0[j] = -ab_y * tg_yz_xxyzz_0[j] + tg_yz_xxyyzz_0[j];

                tg_yyz_xxzzz_0[j] = -ab_y * tg_yz_xxzzz_0[j] + tg_yz_xxyzzz_0[j];

                tg_yyz_xyyyy_0[j] = -ab_y * tg_yz_xyyyy_0[j] + tg_yz_xyyyyy_0[j];

                tg_yyz_xyyyz_0[j] = -ab_y * tg_yz_xyyyz_0[j] + tg_yz_xyyyyz_0[j];

                tg_yyz_xyyzz_0[j] = -ab_y * tg_yz_xyyzz_0[j] + tg_yz_xyyyzz_0[j];

                tg_yyz_xyzzz_0[j] = -ab_y * tg_yz_xyzzz_0[j] + tg_yz_xyyzzz_0[j];

                tg_yyz_xzzzz_0[j] = -ab_y * tg_yz_xzzzz_0[j] + tg_yz_xyzzzz_0[j];

                tg_yyz_yyyyy_0[j] = -ab_y * tg_yz_yyyyy_0[j] + tg_yz_yyyyyy_0[j];

                tg_yyz_yyyyz_0[j] = -ab_y * tg_yz_yyyyz_0[j] + tg_yz_yyyyyz_0[j];

                tg_yyz_yyyzz_0[j] = -ab_y * tg_yz_yyyzz_0[j] + tg_yz_yyyyzz_0[j];

                tg_yyz_yyzzz_0[j] = -ab_y * tg_yz_yyzzz_0[j] + tg_yz_yyyzzz_0[j];

                tg_yyz_yzzzz_0[j] = -ab_y * tg_yz_yzzzz_0[j] + tg_yz_yyzzzz_0[j];

                tg_yyz_zzzzz_0[j] = -ab_y * tg_yz_zzzzz_0[j] + tg_yz_yzzzzz_0[j];

                tg_yzz_xxxxx_0[j] = -ab_y * tg_zz_xxxxx_0[j] + tg_zz_xxxxxy_0[j];

                tg_yzz_xxxxy_0[j] = -ab_y * tg_zz_xxxxy_0[j] + tg_zz_xxxxyy_0[j];

                tg_yzz_xxxxz_0[j] = -ab_y * tg_zz_xxxxz_0[j] + tg_zz_xxxxyz_0[j];

                tg_yzz_xxxyy_0[j] = -ab_y * tg_zz_xxxyy_0[j] + tg_zz_xxxyyy_0[j];

                tg_yzz_xxxyz_0[j] = -ab_y * tg_zz_xxxyz_0[j] + tg_zz_xxxyyz_0[j];

                tg_yzz_xxxzz_0[j] = -ab_y * tg_zz_xxxzz_0[j] + tg_zz_xxxyzz_0[j];

                tg_yzz_xxyyy_0[j] = -ab_y * tg_zz_xxyyy_0[j] + tg_zz_xxyyyy_0[j];

                tg_yzz_xxyyz_0[j] = -ab_y * tg_zz_xxyyz_0[j] + tg_zz_xxyyyz_0[j];

                tg_yzz_xxyzz_0[j] = -ab_y * tg_zz_xxyzz_0[j] + tg_zz_xxyyzz_0[j];

                tg_yzz_xxzzz_0[j] = -ab_y * tg_zz_xxzzz_0[j] + tg_zz_xxyzzz_0[j];

                tg_yzz_xyyyy_0[j] = -ab_y * tg_zz_xyyyy_0[j] + tg_zz_xyyyyy_0[j];

                tg_yzz_xyyyz_0[j] = -ab_y * tg_zz_xyyyz_0[j] + tg_zz_xyyyyz_0[j];

                tg_yzz_xyyzz_0[j] = -ab_y * tg_zz_xyyzz_0[j] + tg_zz_xyyyzz_0[j];

                tg_yzz_xyzzz_0[j] = -ab_y * tg_zz_xyzzz_0[j] + tg_zz_xyyzzz_0[j];

                tg_yzz_xzzzz_0[j] = -ab_y * tg_zz_xzzzz_0[j] + tg_zz_xyzzzz_0[j];

                tg_yzz_yyyyy_0[j] = -ab_y * tg_zz_yyyyy_0[j] + tg_zz_yyyyyy_0[j];

                tg_yzz_yyyyz_0[j] = -ab_y * tg_zz_yyyyz_0[j] + tg_zz_yyyyyz_0[j];

                tg_yzz_yyyzz_0[j] = -ab_y * tg_zz_yyyzz_0[j] + tg_zz_yyyyzz_0[j];

                tg_yzz_yyzzz_0[j] = -ab_y * tg_zz_yyzzz_0[j] + tg_zz_yyyzzz_0[j];

                tg_yzz_yzzzz_0[j] = -ab_y * tg_zz_yzzzz_0[j] + tg_zz_yyzzzz_0[j];

                tg_yzz_zzzzz_0[j] = -ab_y * tg_zz_zzzzz_0[j] + tg_zz_yzzzzz_0[j];

                tg_zzz_xxxxx_0[j] = -ab_z * tg_zz_xxxxx_0[j] + tg_zz_xxxxxz_0[j];

                tg_zzz_xxxxy_0[j] = -ab_z * tg_zz_xxxxy_0[j] + tg_zz_xxxxyz_0[j];

                tg_zzz_xxxxz_0[j] = -ab_z * tg_zz_xxxxz_0[j] + tg_zz_xxxxzz_0[j];

                tg_zzz_xxxyy_0[j] = -ab_z * tg_zz_xxxyy_0[j] + tg_zz_xxxyyz_0[j];

                tg_zzz_xxxyz_0[j] = -ab_z * tg_zz_xxxyz_0[j] + tg_zz_xxxyzz_0[j];

                tg_zzz_xxxzz_0[j] = -ab_z * tg_zz_xxxzz_0[j] + tg_zz_xxxzzz_0[j];

                tg_zzz_xxyyy_0[j] = -ab_z * tg_zz_xxyyy_0[j] + tg_zz_xxyyyz_0[j];

                tg_zzz_xxyyz_0[j] = -ab_z * tg_zz_xxyyz_0[j] + tg_zz_xxyyzz_0[j];

                tg_zzz_xxyzz_0[j] = -ab_z * tg_zz_xxyzz_0[j] + tg_zz_xxyzzz_0[j];

                tg_zzz_xxzzz_0[j] = -ab_z * tg_zz_xxzzz_0[j] + tg_zz_xxzzzz_0[j];

                tg_zzz_xyyyy_0[j] = -ab_z * tg_zz_xyyyy_0[j] + tg_zz_xyyyyz_0[j];

                tg_zzz_xyyyz_0[j] = -ab_z * tg_zz_xyyyz_0[j] + tg_zz_xyyyzz_0[j];

                tg_zzz_xyyzz_0[j] = -ab_z * tg_zz_xyyzz_0[j] + tg_zz_xyyzzz_0[j];

                tg_zzz_xyzzz_0[j] = -ab_z * tg_zz_xyzzz_0[j] + tg_zz_xyzzzz_0[j];

                tg_zzz_xzzzz_0[j] = -ab_z * tg_zz_xzzzz_0[j] + tg_zz_xzzzzz_0[j];

                tg_zzz_yyyyy_0[j] = -ab_z * tg_zz_yyyyy_0[j] + tg_zz_yyyyyz_0[j];

                tg_zzz_yyyyz_0[j] = -ab_z * tg_zz_yyyyz_0[j] + tg_zz_yyyyzz_0[j];

                tg_zzz_yyyzz_0[j] = -ab_z * tg_zz_yyyzz_0[j] + tg_zz_yyyzzz_0[j];

                tg_zzz_yyzzz_0[j] = -ab_z * tg_zz_yyzzz_0[j] + tg_zz_yyzzzz_0[j];

                tg_zzz_yzzzz_0[j] = -ab_z * tg_zz_yzzzz_0[j] + tg_zz_yzzzzz_0[j];

                tg_zzz_zzzzz_0[j] = -ab_z * tg_zz_zzzzz_0[j] + tg_zz_zzzzzz_0[j];
            }
        }
    }


} // eribrrfunc namespace

