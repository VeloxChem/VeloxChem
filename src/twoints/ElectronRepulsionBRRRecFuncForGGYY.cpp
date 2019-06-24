//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionBRRRecFuncForGGYY.hpp"

#include "AngularMomentum.hpp"

namespace eribrrfunc { // eribrrfunc namespace

    void
    compElectronRepulsionForGGXY(      CMemBlock2D<double>& braBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& abDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        eribrrfunc::compElectronRepulsionForGGXY_0_75(braBuffer,
                                                      recursionMap,
                                                      abDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        eribrrfunc::compElectronRepulsionForGGXY_75_150(braBuffer,
                                                        recursionMap,
                                                        abDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 

        eribrrfunc::compElectronRepulsionForGGXY_150_225(braBuffer,
                                                         recursionMap,
                                                         abDistances,
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetContrPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForGGXY_0_75(      CMemBlock2D<double>& braBuffer,
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

        auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {4, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_4_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            auto tg_xyy_xxxxx_0 = braBuffer.data(pidx_g_3_5_m0 + 63 * kcomp + i); 

            auto tg_xyy_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 64 * kcomp + i); 

            auto tg_xyy_xxxxz_0 = braBuffer.data(pidx_g_3_5_m0 + 65 * kcomp + i); 

            auto tg_xyy_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 66 * kcomp + i); 

            auto tg_xyy_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 67 * kcomp + i); 

            auto tg_xyy_xxxzz_0 = braBuffer.data(pidx_g_3_5_m0 + 68 * kcomp + i); 

            auto tg_xyy_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 69 * kcomp + i); 

            auto tg_xyy_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 70 * kcomp + i); 

            auto tg_xyy_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 71 * kcomp + i); 

            auto tg_xyy_xxzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 72 * kcomp + i); 

            auto tg_xyy_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 73 * kcomp + i); 

            auto tg_xyy_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 74 * kcomp + i); 

            auto tg_xyy_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 75 * kcomp + i); 

            auto tg_xyy_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 76 * kcomp + i); 

            auto tg_xyy_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 77 * kcomp + i); 

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

            // set up pointers to integrals

            auto tg_xxxx_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + i); 

            auto tg_xxxx_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + kcomp + i); 

            auto tg_xxxx_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 2 * kcomp + i); 

            auto tg_xxxx_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 3 * kcomp + i); 

            auto tg_xxxx_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 4 * kcomp + i); 

            auto tg_xxxx_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 5 * kcomp + i); 

            auto tg_xxxx_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 6 * kcomp + i); 

            auto tg_xxxx_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 7 * kcomp + i); 

            auto tg_xxxx_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 8 * kcomp + i); 

            auto tg_xxxx_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 9 * kcomp + i); 

            auto tg_xxxx_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 10 * kcomp + i); 

            auto tg_xxxx_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 11 * kcomp + i); 

            auto tg_xxxx_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 12 * kcomp + i); 

            auto tg_xxxx_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 13 * kcomp + i); 

            auto tg_xxxx_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 14 * kcomp + i); 

            auto tg_xxxy_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 15 * kcomp + i); 

            auto tg_xxxy_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 16 * kcomp + i); 

            auto tg_xxxy_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 17 * kcomp + i); 

            auto tg_xxxy_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 18 * kcomp + i); 

            auto tg_xxxy_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 19 * kcomp + i); 

            auto tg_xxxy_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 20 * kcomp + i); 

            auto tg_xxxy_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 21 * kcomp + i); 

            auto tg_xxxy_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 22 * kcomp + i); 

            auto tg_xxxy_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 23 * kcomp + i); 

            auto tg_xxxy_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 24 * kcomp + i); 

            auto tg_xxxy_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 25 * kcomp + i); 

            auto tg_xxxy_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 26 * kcomp + i); 

            auto tg_xxxy_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 27 * kcomp + i); 

            auto tg_xxxy_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 28 * kcomp + i); 

            auto tg_xxxy_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 29 * kcomp + i); 

            auto tg_xxxz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 30 * kcomp + i); 

            auto tg_xxxz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 31 * kcomp + i); 

            auto tg_xxxz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 32 * kcomp + i); 

            auto tg_xxxz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 33 * kcomp + i); 

            auto tg_xxxz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 34 * kcomp + i); 

            auto tg_xxxz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 35 * kcomp + i); 

            auto tg_xxxz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 36 * kcomp + i); 

            auto tg_xxxz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 37 * kcomp + i); 

            auto tg_xxxz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 38 * kcomp + i); 

            auto tg_xxxz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 39 * kcomp + i); 

            auto tg_xxxz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 40 * kcomp + i); 

            auto tg_xxxz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 41 * kcomp + i); 

            auto tg_xxxz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 42 * kcomp + i); 

            auto tg_xxxz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 43 * kcomp + i); 

            auto tg_xxxz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 44 * kcomp + i); 

            auto tg_xxyy_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 45 * kcomp + i); 

            auto tg_xxyy_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 46 * kcomp + i); 

            auto tg_xxyy_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 47 * kcomp + i); 

            auto tg_xxyy_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 48 * kcomp + i); 

            auto tg_xxyy_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 49 * kcomp + i); 

            auto tg_xxyy_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 50 * kcomp + i); 

            auto tg_xxyy_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 51 * kcomp + i); 

            auto tg_xxyy_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 52 * kcomp + i); 

            auto tg_xxyy_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 53 * kcomp + i); 

            auto tg_xxyy_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 54 * kcomp + i); 

            auto tg_xxyy_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 55 * kcomp + i); 

            auto tg_xxyy_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 56 * kcomp + i); 

            auto tg_xxyy_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 57 * kcomp + i); 

            auto tg_xxyy_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 58 * kcomp + i); 

            auto tg_xxyy_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 59 * kcomp + i); 

            auto tg_xxyz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 60 * kcomp + i); 

            auto tg_xxyz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 61 * kcomp + i); 

            auto tg_xxyz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 62 * kcomp + i); 

            auto tg_xxyz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 63 * kcomp + i); 

            auto tg_xxyz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 64 * kcomp + i); 

            auto tg_xxyz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 65 * kcomp + i); 

            auto tg_xxyz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 66 * kcomp + i); 

            auto tg_xxyz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 67 * kcomp + i); 

            auto tg_xxyz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 68 * kcomp + i); 

            auto tg_xxyz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 69 * kcomp + i); 

            auto tg_xxyz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 70 * kcomp + i); 

            auto tg_xxyz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 71 * kcomp + i); 

            auto tg_xxyz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 72 * kcomp + i); 

            auto tg_xxyz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 73 * kcomp + i); 

            auto tg_xxyz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 74 * kcomp + i); 

            // Batch of Integrals (0,75)

            #pragma omp simd aligned(tg_xxx_xxxx_0, tg_xxx_xxxxx_0, tg_xxx_xxxxy_0, tg_xxx_xxxxz_0, \
                                     tg_xxx_xxxy_0, tg_xxx_xxxyy_0, tg_xxx_xxxyz_0, tg_xxx_xxxz_0, tg_xxx_xxxzz_0, \
                                     tg_xxx_xxyy_0, tg_xxx_xxyyy_0, tg_xxx_xxyyz_0, tg_xxx_xxyz_0, tg_xxx_xxyzz_0, \
                                     tg_xxx_xxzz_0, tg_xxx_xxzzz_0, tg_xxx_xyyy_0, tg_xxx_xyyyy_0, tg_xxx_xyyyz_0, \
                                     tg_xxx_xyyz_0, tg_xxx_xyyzz_0, tg_xxx_xyzz_0, tg_xxx_xyzzz_0, tg_xxx_xzzz_0, \
                                     tg_xxx_xzzzz_0, tg_xxx_yyyy_0, tg_xxx_yyyz_0, tg_xxx_yyzz_0, tg_xxx_yzzz_0, \
                                     tg_xxx_zzzz_0, tg_xxxx_xxxx_0, tg_xxxx_xxxy_0, tg_xxxx_xxxz_0, tg_xxxx_xxyy_0, \
                                     tg_xxxx_xxyz_0, tg_xxxx_xxzz_0, tg_xxxx_xyyy_0, tg_xxxx_xyyz_0, tg_xxxx_xyzz_0, \
                                     tg_xxxx_xzzz_0, tg_xxxx_yyyy_0, tg_xxxx_yyyz_0, tg_xxxx_yyzz_0, tg_xxxx_yzzz_0, \
                                     tg_xxxx_zzzz_0, tg_xxxy_xxxx_0, tg_xxxy_xxxy_0, tg_xxxy_xxxz_0, tg_xxxy_xxyy_0, \
                                     tg_xxxy_xxyz_0, tg_xxxy_xxzz_0, tg_xxxy_xyyy_0, tg_xxxy_xyyz_0, tg_xxxy_xyzz_0, \
                                     tg_xxxy_xzzz_0, tg_xxxy_yyyy_0, tg_xxxy_yyyz_0, tg_xxxy_yyzz_0, tg_xxxy_yzzz_0, \
                                     tg_xxxy_zzzz_0, tg_xxxz_xxxx_0, tg_xxxz_xxxy_0, tg_xxxz_xxxz_0, tg_xxxz_xxyy_0, \
                                     tg_xxxz_xxyz_0, tg_xxxz_xxzz_0, tg_xxxz_xyyy_0, tg_xxxz_xyyz_0, tg_xxxz_xyzz_0, \
                                     tg_xxxz_xzzz_0, tg_xxxz_yyyy_0, tg_xxxz_yyyz_0, tg_xxxz_yyzz_0, tg_xxxz_yzzz_0, \
                                     tg_xxxz_zzzz_0, tg_xxy_xxxx_0, tg_xxy_xxxxx_0, tg_xxy_xxxxy_0, tg_xxy_xxxxz_0, \
                                     tg_xxy_xxxy_0, tg_xxy_xxxyy_0, tg_xxy_xxxyz_0, tg_xxy_xxxz_0, tg_xxy_xxxzz_0, \
                                     tg_xxy_xxyy_0, tg_xxy_xxyyy_0, tg_xxy_xxyyz_0, tg_xxy_xxyz_0, tg_xxy_xxyzz_0, \
                                     tg_xxy_xxzz_0, tg_xxy_xxzzz_0, tg_xxy_xyyy_0, tg_xxy_xyyyy_0, tg_xxy_xyyyz_0, \
                                     tg_xxy_xyyz_0, tg_xxy_xyyzz_0, tg_xxy_xyzz_0, tg_xxy_xyzzz_0, tg_xxy_xzzz_0, \
                                     tg_xxy_xzzzz_0, tg_xxy_yyyy_0, tg_xxy_yyyz_0, tg_xxy_yyzz_0, tg_xxy_yzzz_0, \
                                     tg_xxy_zzzz_0, tg_xxyy_xxxx_0, tg_xxyy_xxxy_0, tg_xxyy_xxxz_0, tg_xxyy_xxyy_0, \
                                     tg_xxyy_xxyz_0, tg_xxyy_xxzz_0, tg_xxyy_xyyy_0, tg_xxyy_xyyz_0, tg_xxyy_xyzz_0, \
                                     tg_xxyy_xzzz_0, tg_xxyy_yyyy_0, tg_xxyy_yyyz_0, tg_xxyy_yyzz_0, tg_xxyy_yzzz_0, \
                                     tg_xxyy_zzzz_0, tg_xxyz_xxxx_0, tg_xxyz_xxxy_0, tg_xxyz_xxxz_0, tg_xxyz_xxyy_0, \
                                     tg_xxyz_xxyz_0, tg_xxyz_xxzz_0, tg_xxyz_xyyy_0, tg_xxyz_xyyz_0, tg_xxyz_xyzz_0, \
                                     tg_xxyz_xzzz_0, tg_xxyz_yyyy_0, tg_xxyz_yyyz_0, tg_xxyz_yyzz_0, tg_xxyz_yzzz_0, \
                                     tg_xxyz_zzzz_0, tg_xxz_xxxx_0, tg_xxz_xxxxx_0, tg_xxz_xxxxy_0, tg_xxz_xxxxz_0, \
                                     tg_xxz_xxxy_0, tg_xxz_xxxyy_0, tg_xxz_xxxyz_0, tg_xxz_xxxz_0, tg_xxz_xxxzz_0, \
                                     tg_xxz_xxyy_0, tg_xxz_xxyyy_0, tg_xxz_xxyyz_0, tg_xxz_xxyz_0, tg_xxz_xxyzz_0, \
                                     tg_xxz_xxzz_0, tg_xxz_xxzzz_0, tg_xxz_xyyy_0, tg_xxz_xyyyy_0, tg_xxz_xyyyz_0, \
                                     tg_xxz_xyyz_0, tg_xxz_xyyzz_0, tg_xxz_xyzz_0, tg_xxz_xyzzz_0, tg_xxz_xzzz_0, \
                                     tg_xxz_xzzzz_0, tg_xxz_yyyy_0, tg_xxz_yyyz_0, tg_xxz_yyzz_0, tg_xxz_yzzz_0, \
                                     tg_xxz_zzzz_0, tg_xyy_xxxx_0, tg_xyy_xxxxx_0, tg_xyy_xxxxy_0, tg_xyy_xxxxz_0, \
                                     tg_xyy_xxxy_0, tg_xyy_xxxyy_0, tg_xyy_xxxyz_0, tg_xyy_xxxz_0, tg_xyy_xxxzz_0, \
                                     tg_xyy_xxyy_0, tg_xyy_xxyyy_0, tg_xyy_xxyyz_0, tg_xyy_xxyz_0, tg_xyy_xxyzz_0, \
                                     tg_xyy_xxzz_0, tg_xyy_xxzzz_0, tg_xyy_xyyy_0, tg_xyy_xyyyy_0, tg_xyy_xyyyz_0, \
                                     tg_xyy_xyyz_0, tg_xyy_xyyzz_0, tg_xyy_xyzz_0, tg_xyy_xyzzz_0, tg_xyy_xzzz_0, \
                                     tg_xyy_xzzzz_0, tg_xyy_yyyy_0, tg_xyy_yyyz_0, tg_xyy_yyzz_0, tg_xyy_yzzz_0, \
                                     tg_xyy_zzzz_0, tg_xyz_xxxx_0, tg_xyz_xxxxx_0, tg_xyz_xxxxy_0, tg_xyz_xxxxz_0, \
                                     tg_xyz_xxxy_0, tg_xyz_xxxyy_0, tg_xyz_xxxyz_0, tg_xyz_xxxz_0, tg_xyz_xxxzz_0, \
                                     tg_xyz_xxyy_0, tg_xyz_xxyyy_0, tg_xyz_xxyyz_0, tg_xyz_xxyz_0, tg_xyz_xxyzz_0, \
                                     tg_xyz_xxzz_0, tg_xyz_xxzzz_0, tg_xyz_xyyy_0, tg_xyz_xyyyy_0, tg_xyz_xyyyz_0, \
                                     tg_xyz_xyyz_0, tg_xyz_xyyzz_0, tg_xyz_xyzz_0, tg_xyz_xyzzz_0, tg_xyz_xzzz_0, \
                                     tg_xyz_xzzzz_0, tg_xyz_yyyy_0, tg_xyz_yyyz_0, tg_xyz_yyzz_0, tg_xyz_yzzz_0, \
                                     tg_xyz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xxxx_xxxx_0[j] = -ab_x * tg_xxx_xxxx_0[j] + tg_xxx_xxxxx_0[j];

                tg_xxxx_xxxy_0[j] = -ab_x * tg_xxx_xxxy_0[j] + tg_xxx_xxxxy_0[j];

                tg_xxxx_xxxz_0[j] = -ab_x * tg_xxx_xxxz_0[j] + tg_xxx_xxxxz_0[j];

                tg_xxxx_xxyy_0[j] = -ab_x * tg_xxx_xxyy_0[j] + tg_xxx_xxxyy_0[j];

                tg_xxxx_xxyz_0[j] = -ab_x * tg_xxx_xxyz_0[j] + tg_xxx_xxxyz_0[j];

                tg_xxxx_xxzz_0[j] = -ab_x * tg_xxx_xxzz_0[j] + tg_xxx_xxxzz_0[j];

                tg_xxxx_xyyy_0[j] = -ab_x * tg_xxx_xyyy_0[j] + tg_xxx_xxyyy_0[j];

                tg_xxxx_xyyz_0[j] = -ab_x * tg_xxx_xyyz_0[j] + tg_xxx_xxyyz_0[j];

                tg_xxxx_xyzz_0[j] = -ab_x * tg_xxx_xyzz_0[j] + tg_xxx_xxyzz_0[j];

                tg_xxxx_xzzz_0[j] = -ab_x * tg_xxx_xzzz_0[j] + tg_xxx_xxzzz_0[j];

                tg_xxxx_yyyy_0[j] = -ab_x * tg_xxx_yyyy_0[j] + tg_xxx_xyyyy_0[j];

                tg_xxxx_yyyz_0[j] = -ab_x * tg_xxx_yyyz_0[j] + tg_xxx_xyyyz_0[j];

                tg_xxxx_yyzz_0[j] = -ab_x * tg_xxx_yyzz_0[j] + tg_xxx_xyyzz_0[j];

                tg_xxxx_yzzz_0[j] = -ab_x * tg_xxx_yzzz_0[j] + tg_xxx_xyzzz_0[j];

                tg_xxxx_zzzz_0[j] = -ab_x * tg_xxx_zzzz_0[j] + tg_xxx_xzzzz_0[j];

                tg_xxxy_xxxx_0[j] = -ab_x * tg_xxy_xxxx_0[j] + tg_xxy_xxxxx_0[j];

                tg_xxxy_xxxy_0[j] = -ab_x * tg_xxy_xxxy_0[j] + tg_xxy_xxxxy_0[j];

                tg_xxxy_xxxz_0[j] = -ab_x * tg_xxy_xxxz_0[j] + tg_xxy_xxxxz_0[j];

                tg_xxxy_xxyy_0[j] = -ab_x * tg_xxy_xxyy_0[j] + tg_xxy_xxxyy_0[j];

                tg_xxxy_xxyz_0[j] = -ab_x * tg_xxy_xxyz_0[j] + tg_xxy_xxxyz_0[j];

                tg_xxxy_xxzz_0[j] = -ab_x * tg_xxy_xxzz_0[j] + tg_xxy_xxxzz_0[j];

                tg_xxxy_xyyy_0[j] = -ab_x * tg_xxy_xyyy_0[j] + tg_xxy_xxyyy_0[j];

                tg_xxxy_xyyz_0[j] = -ab_x * tg_xxy_xyyz_0[j] + tg_xxy_xxyyz_0[j];

                tg_xxxy_xyzz_0[j] = -ab_x * tg_xxy_xyzz_0[j] + tg_xxy_xxyzz_0[j];

                tg_xxxy_xzzz_0[j] = -ab_x * tg_xxy_xzzz_0[j] + tg_xxy_xxzzz_0[j];

                tg_xxxy_yyyy_0[j] = -ab_x * tg_xxy_yyyy_0[j] + tg_xxy_xyyyy_0[j];

                tg_xxxy_yyyz_0[j] = -ab_x * tg_xxy_yyyz_0[j] + tg_xxy_xyyyz_0[j];

                tg_xxxy_yyzz_0[j] = -ab_x * tg_xxy_yyzz_0[j] + tg_xxy_xyyzz_0[j];

                tg_xxxy_yzzz_0[j] = -ab_x * tg_xxy_yzzz_0[j] + tg_xxy_xyzzz_0[j];

                tg_xxxy_zzzz_0[j] = -ab_x * tg_xxy_zzzz_0[j] + tg_xxy_xzzzz_0[j];

                tg_xxxz_xxxx_0[j] = -ab_x * tg_xxz_xxxx_0[j] + tg_xxz_xxxxx_0[j];

                tg_xxxz_xxxy_0[j] = -ab_x * tg_xxz_xxxy_0[j] + tg_xxz_xxxxy_0[j];

                tg_xxxz_xxxz_0[j] = -ab_x * tg_xxz_xxxz_0[j] + tg_xxz_xxxxz_0[j];

                tg_xxxz_xxyy_0[j] = -ab_x * tg_xxz_xxyy_0[j] + tg_xxz_xxxyy_0[j];

                tg_xxxz_xxyz_0[j] = -ab_x * tg_xxz_xxyz_0[j] + tg_xxz_xxxyz_0[j];

                tg_xxxz_xxzz_0[j] = -ab_x * tg_xxz_xxzz_0[j] + tg_xxz_xxxzz_0[j];

                tg_xxxz_xyyy_0[j] = -ab_x * tg_xxz_xyyy_0[j] + tg_xxz_xxyyy_0[j];

                tg_xxxz_xyyz_0[j] = -ab_x * tg_xxz_xyyz_0[j] + tg_xxz_xxyyz_0[j];

                tg_xxxz_xyzz_0[j] = -ab_x * tg_xxz_xyzz_0[j] + tg_xxz_xxyzz_0[j];

                tg_xxxz_xzzz_0[j] = -ab_x * tg_xxz_xzzz_0[j] + tg_xxz_xxzzz_0[j];

                tg_xxxz_yyyy_0[j] = -ab_x * tg_xxz_yyyy_0[j] + tg_xxz_xyyyy_0[j];

                tg_xxxz_yyyz_0[j] = -ab_x * tg_xxz_yyyz_0[j] + tg_xxz_xyyyz_0[j];

                tg_xxxz_yyzz_0[j] = -ab_x * tg_xxz_yyzz_0[j] + tg_xxz_xyyzz_0[j];

                tg_xxxz_yzzz_0[j] = -ab_x * tg_xxz_yzzz_0[j] + tg_xxz_xyzzz_0[j];

                tg_xxxz_zzzz_0[j] = -ab_x * tg_xxz_zzzz_0[j] + tg_xxz_xzzzz_0[j];

                tg_xxyy_xxxx_0[j] = -ab_x * tg_xyy_xxxx_0[j] + tg_xyy_xxxxx_0[j];

                tg_xxyy_xxxy_0[j] = -ab_x * tg_xyy_xxxy_0[j] + tg_xyy_xxxxy_0[j];

                tg_xxyy_xxxz_0[j] = -ab_x * tg_xyy_xxxz_0[j] + tg_xyy_xxxxz_0[j];

                tg_xxyy_xxyy_0[j] = -ab_x * tg_xyy_xxyy_0[j] + tg_xyy_xxxyy_0[j];

                tg_xxyy_xxyz_0[j] = -ab_x * tg_xyy_xxyz_0[j] + tg_xyy_xxxyz_0[j];

                tg_xxyy_xxzz_0[j] = -ab_x * tg_xyy_xxzz_0[j] + tg_xyy_xxxzz_0[j];

                tg_xxyy_xyyy_0[j] = -ab_x * tg_xyy_xyyy_0[j] + tg_xyy_xxyyy_0[j];

                tg_xxyy_xyyz_0[j] = -ab_x * tg_xyy_xyyz_0[j] + tg_xyy_xxyyz_0[j];

                tg_xxyy_xyzz_0[j] = -ab_x * tg_xyy_xyzz_0[j] + tg_xyy_xxyzz_0[j];

                tg_xxyy_xzzz_0[j] = -ab_x * tg_xyy_xzzz_0[j] + tg_xyy_xxzzz_0[j];

                tg_xxyy_yyyy_0[j] = -ab_x * tg_xyy_yyyy_0[j] + tg_xyy_xyyyy_0[j];

                tg_xxyy_yyyz_0[j] = -ab_x * tg_xyy_yyyz_0[j] + tg_xyy_xyyyz_0[j];

                tg_xxyy_yyzz_0[j] = -ab_x * tg_xyy_yyzz_0[j] + tg_xyy_xyyzz_0[j];

                tg_xxyy_yzzz_0[j] = -ab_x * tg_xyy_yzzz_0[j] + tg_xyy_xyzzz_0[j];

                tg_xxyy_zzzz_0[j] = -ab_x * tg_xyy_zzzz_0[j] + tg_xyy_xzzzz_0[j];

                tg_xxyz_xxxx_0[j] = -ab_x * tg_xyz_xxxx_0[j] + tg_xyz_xxxxx_0[j];

                tg_xxyz_xxxy_0[j] = -ab_x * tg_xyz_xxxy_0[j] + tg_xyz_xxxxy_0[j];

                tg_xxyz_xxxz_0[j] = -ab_x * tg_xyz_xxxz_0[j] + tg_xyz_xxxxz_0[j];

                tg_xxyz_xxyy_0[j] = -ab_x * tg_xyz_xxyy_0[j] + tg_xyz_xxxyy_0[j];

                tg_xxyz_xxyz_0[j] = -ab_x * tg_xyz_xxyz_0[j] + tg_xyz_xxxyz_0[j];

                tg_xxyz_xxzz_0[j] = -ab_x * tg_xyz_xxzz_0[j] + tg_xyz_xxxzz_0[j];

                tg_xxyz_xyyy_0[j] = -ab_x * tg_xyz_xyyy_0[j] + tg_xyz_xxyyy_0[j];

                tg_xxyz_xyyz_0[j] = -ab_x * tg_xyz_xyyz_0[j] + tg_xyz_xxyyz_0[j];

                tg_xxyz_xyzz_0[j] = -ab_x * tg_xyz_xyzz_0[j] + tg_xyz_xxyzz_0[j];

                tg_xxyz_xzzz_0[j] = -ab_x * tg_xyz_xzzz_0[j] + tg_xyz_xxzzz_0[j];

                tg_xxyz_yyyy_0[j] = -ab_x * tg_xyz_yyyy_0[j] + tg_xyz_xyyyy_0[j];

                tg_xxyz_yyyz_0[j] = -ab_x * tg_xyz_yyyz_0[j] + tg_xyz_xyyyz_0[j];

                tg_xxyz_yyzz_0[j] = -ab_x * tg_xyz_yyzz_0[j] + tg_xyz_xyyzz_0[j];

                tg_xxyz_yzzz_0[j] = -ab_x * tg_xyz_yzzz_0[j] + tg_xyz_xyzzz_0[j];

                tg_xxyz_zzzz_0[j] = -ab_x * tg_xyz_zzzz_0[j] + tg_xyz_xzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForGGXY_75_150(      CMemBlock2D<double>& braBuffer,
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

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {4, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_4_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            auto tg_yyy_xzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 140 * kcomp + i); 

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

            // set up pointers to integrals

            auto tg_xxzz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 75 * kcomp + i); 

            auto tg_xxzz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 76 * kcomp + i); 

            auto tg_xxzz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 77 * kcomp + i); 

            auto tg_xxzz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 78 * kcomp + i); 

            auto tg_xxzz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 79 * kcomp + i); 

            auto tg_xxzz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 80 * kcomp + i); 

            auto tg_xxzz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 81 * kcomp + i); 

            auto tg_xxzz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 82 * kcomp + i); 

            auto tg_xxzz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 83 * kcomp + i); 

            auto tg_xxzz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 84 * kcomp + i); 

            auto tg_xxzz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 85 * kcomp + i); 

            auto tg_xxzz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 86 * kcomp + i); 

            auto tg_xxzz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 87 * kcomp + i); 

            auto tg_xxzz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 88 * kcomp + i); 

            auto tg_xxzz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 89 * kcomp + i); 

            auto tg_xyyy_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 90 * kcomp + i); 

            auto tg_xyyy_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 91 * kcomp + i); 

            auto tg_xyyy_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 92 * kcomp + i); 

            auto tg_xyyy_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 93 * kcomp + i); 

            auto tg_xyyy_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 94 * kcomp + i); 

            auto tg_xyyy_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 95 * kcomp + i); 

            auto tg_xyyy_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 96 * kcomp + i); 

            auto tg_xyyy_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 97 * kcomp + i); 

            auto tg_xyyy_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 98 * kcomp + i); 

            auto tg_xyyy_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 99 * kcomp + i); 

            auto tg_xyyy_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 100 * kcomp + i); 

            auto tg_xyyy_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 101 * kcomp + i); 

            auto tg_xyyy_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 102 * kcomp + i); 

            auto tg_xyyy_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 103 * kcomp + i); 

            auto tg_xyyy_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 104 * kcomp + i); 

            auto tg_xyyz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 105 * kcomp + i); 

            auto tg_xyyz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 106 * kcomp + i); 

            auto tg_xyyz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 107 * kcomp + i); 

            auto tg_xyyz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 108 * kcomp + i); 

            auto tg_xyyz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 109 * kcomp + i); 

            auto tg_xyyz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 110 * kcomp + i); 

            auto tg_xyyz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 111 * kcomp + i); 

            auto tg_xyyz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 112 * kcomp + i); 

            auto tg_xyyz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 113 * kcomp + i); 

            auto tg_xyyz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 114 * kcomp + i); 

            auto tg_xyyz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 115 * kcomp + i); 

            auto tg_xyyz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 116 * kcomp + i); 

            auto tg_xyyz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 117 * kcomp + i); 

            auto tg_xyyz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 118 * kcomp + i); 

            auto tg_xyyz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 119 * kcomp + i); 

            auto tg_xyzz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 120 * kcomp + i); 

            auto tg_xyzz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 121 * kcomp + i); 

            auto tg_xyzz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 122 * kcomp + i); 

            auto tg_xyzz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 123 * kcomp + i); 

            auto tg_xyzz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 124 * kcomp + i); 

            auto tg_xyzz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 125 * kcomp + i); 

            auto tg_xyzz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 126 * kcomp + i); 

            auto tg_xyzz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 127 * kcomp + i); 

            auto tg_xyzz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 128 * kcomp + i); 

            auto tg_xyzz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 129 * kcomp + i); 

            auto tg_xyzz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 130 * kcomp + i); 

            auto tg_xyzz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 131 * kcomp + i); 

            auto tg_xyzz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 132 * kcomp + i); 

            auto tg_xyzz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 133 * kcomp + i); 

            auto tg_xyzz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 134 * kcomp + i); 

            auto tg_xzzz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 135 * kcomp + i); 

            auto tg_xzzz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 136 * kcomp + i); 

            auto tg_xzzz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 137 * kcomp + i); 

            auto tg_xzzz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 138 * kcomp + i); 

            auto tg_xzzz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 139 * kcomp + i); 

            auto tg_xzzz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 140 * kcomp + i); 

            auto tg_xzzz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 141 * kcomp + i); 

            auto tg_xzzz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 142 * kcomp + i); 

            auto tg_xzzz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 143 * kcomp + i); 

            auto tg_xzzz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 144 * kcomp + i); 

            auto tg_xzzz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 145 * kcomp + i); 

            auto tg_xzzz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 146 * kcomp + i); 

            auto tg_xzzz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 147 * kcomp + i); 

            auto tg_xzzz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 148 * kcomp + i); 

            auto tg_xzzz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 149 * kcomp + i); 

            // Batch of Integrals (75,150)

            #pragma omp simd aligned(tg_xxzz_xxxx_0, tg_xxzz_xxxy_0, tg_xxzz_xxxz_0, tg_xxzz_xxyy_0, \
                                     tg_xxzz_xxyz_0, tg_xxzz_xxzz_0, tg_xxzz_xyyy_0, tg_xxzz_xyyz_0, tg_xxzz_xyzz_0, \
                                     tg_xxzz_xzzz_0, tg_xxzz_yyyy_0, tg_xxzz_yyyz_0, tg_xxzz_yyzz_0, tg_xxzz_yzzz_0, \
                                     tg_xxzz_zzzz_0, tg_xyyy_xxxx_0, tg_xyyy_xxxy_0, tg_xyyy_xxxz_0, tg_xyyy_xxyy_0, \
                                     tg_xyyy_xxyz_0, tg_xyyy_xxzz_0, tg_xyyy_xyyy_0, tg_xyyy_xyyz_0, tg_xyyy_xyzz_0, \
                                     tg_xyyy_xzzz_0, tg_xyyy_yyyy_0, tg_xyyy_yyyz_0, tg_xyyy_yyzz_0, tg_xyyy_yzzz_0, \
                                     tg_xyyy_zzzz_0, tg_xyyz_xxxx_0, tg_xyyz_xxxy_0, tg_xyyz_xxxz_0, tg_xyyz_xxyy_0, \
                                     tg_xyyz_xxyz_0, tg_xyyz_xxzz_0, tg_xyyz_xyyy_0, tg_xyyz_xyyz_0, tg_xyyz_xyzz_0, \
                                     tg_xyyz_xzzz_0, tg_xyyz_yyyy_0, tg_xyyz_yyyz_0, tg_xyyz_yyzz_0, tg_xyyz_yzzz_0, \
                                     tg_xyyz_zzzz_0, tg_xyzz_xxxx_0, tg_xyzz_xxxy_0, tg_xyzz_xxxz_0, tg_xyzz_xxyy_0, \
                                     tg_xyzz_xxyz_0, tg_xyzz_xxzz_0, tg_xyzz_xyyy_0, tg_xyzz_xyyz_0, tg_xyzz_xyzz_0, \
                                     tg_xyzz_xzzz_0, tg_xyzz_yyyy_0, tg_xyzz_yyyz_0, tg_xyzz_yyzz_0, tg_xyzz_yzzz_0, \
                                     tg_xyzz_zzzz_0, tg_xzz_xxxx_0, tg_xzz_xxxxx_0, tg_xzz_xxxxy_0, tg_xzz_xxxxz_0, \
                                     tg_xzz_xxxy_0, tg_xzz_xxxyy_0, tg_xzz_xxxyz_0, tg_xzz_xxxz_0, tg_xzz_xxxzz_0, \
                                     tg_xzz_xxyy_0, tg_xzz_xxyyy_0, tg_xzz_xxyyz_0, tg_xzz_xxyz_0, tg_xzz_xxyzz_0, \
                                     tg_xzz_xxzz_0, tg_xzz_xxzzz_0, tg_xzz_xyyy_0, tg_xzz_xyyyy_0, tg_xzz_xyyyz_0, \
                                     tg_xzz_xyyz_0, tg_xzz_xyyzz_0, tg_xzz_xyzz_0, tg_xzz_xyzzz_0, tg_xzz_xzzz_0, \
                                     tg_xzz_xzzzz_0, tg_xzz_yyyy_0, tg_xzz_yyyz_0, tg_xzz_yyzz_0, tg_xzz_yzzz_0, \
                                     tg_xzz_zzzz_0, tg_xzzz_xxxx_0, tg_xzzz_xxxy_0, tg_xzzz_xxxz_0, tg_xzzz_xxyy_0, \
                                     tg_xzzz_xxyz_0, tg_xzzz_xxzz_0, tg_xzzz_xyyy_0, tg_xzzz_xyyz_0, tg_xzzz_xyzz_0, \
                                     tg_xzzz_xzzz_0, tg_xzzz_yyyy_0, tg_xzzz_yyyz_0, tg_xzzz_yyzz_0, tg_xzzz_yzzz_0, \
                                     tg_xzzz_zzzz_0, tg_yyy_xxxx_0, tg_yyy_xxxxx_0, tg_yyy_xxxxy_0, tg_yyy_xxxxz_0, \
                                     tg_yyy_xxxy_0, tg_yyy_xxxyy_0, tg_yyy_xxxyz_0, tg_yyy_xxxz_0, tg_yyy_xxxzz_0, \
                                     tg_yyy_xxyy_0, tg_yyy_xxyyy_0, tg_yyy_xxyyz_0, tg_yyy_xxyz_0, tg_yyy_xxyzz_0, \
                                     tg_yyy_xxzz_0, tg_yyy_xxzzz_0, tg_yyy_xyyy_0, tg_yyy_xyyyy_0, tg_yyy_xyyyz_0, \
                                     tg_yyy_xyyz_0, tg_yyy_xyyzz_0, tg_yyy_xyzz_0, tg_yyy_xyzzz_0, tg_yyy_xzzz_0, \
                                     tg_yyy_xzzzz_0, tg_yyy_yyyy_0, tg_yyy_yyyz_0, tg_yyy_yyzz_0, tg_yyy_yzzz_0, \
                                     tg_yyy_zzzz_0, tg_yyz_xxxx_0, tg_yyz_xxxxx_0, tg_yyz_xxxxy_0, tg_yyz_xxxxz_0, \
                                     tg_yyz_xxxy_0, tg_yyz_xxxyy_0, tg_yyz_xxxyz_0, tg_yyz_xxxz_0, tg_yyz_xxxzz_0, \
                                     tg_yyz_xxyy_0, tg_yyz_xxyyy_0, tg_yyz_xxyyz_0, tg_yyz_xxyz_0, tg_yyz_xxyzz_0, \
                                     tg_yyz_xxzz_0, tg_yyz_xxzzz_0, tg_yyz_xyyy_0, tg_yyz_xyyyy_0, tg_yyz_xyyyz_0, \
                                     tg_yyz_xyyz_0, tg_yyz_xyyzz_0, tg_yyz_xyzz_0, tg_yyz_xyzzz_0, tg_yyz_xzzz_0, \
                                     tg_yyz_xzzzz_0, tg_yyz_yyyy_0, tg_yyz_yyyz_0, tg_yyz_yyzz_0, tg_yyz_yzzz_0, \
                                     tg_yyz_zzzz_0, tg_yzz_xxxx_0, tg_yzz_xxxxx_0, tg_yzz_xxxxy_0, tg_yzz_xxxxz_0, \
                                     tg_yzz_xxxy_0, tg_yzz_xxxyy_0, tg_yzz_xxxyz_0, tg_yzz_xxxz_0, tg_yzz_xxxzz_0, \
                                     tg_yzz_xxyy_0, tg_yzz_xxyyy_0, tg_yzz_xxyyz_0, tg_yzz_xxyz_0, tg_yzz_xxyzz_0, \
                                     tg_yzz_xxzz_0, tg_yzz_xxzzz_0, tg_yzz_xyyy_0, tg_yzz_xyyyy_0, tg_yzz_xyyyz_0, \
                                     tg_yzz_xyyz_0, tg_yzz_xyyzz_0, tg_yzz_xyzz_0, tg_yzz_xyzzz_0, tg_yzz_xzzz_0, \
                                     tg_yzz_xzzzz_0, tg_yzz_yyyy_0, tg_yzz_yyyz_0, tg_yzz_yyzz_0, tg_yzz_yzzz_0, \
                                     tg_yzz_zzzz_0, tg_zzz_xxxx_0, tg_zzz_xxxxx_0, tg_zzz_xxxxy_0, tg_zzz_xxxxz_0, \
                                     tg_zzz_xxxy_0, tg_zzz_xxxyy_0, tg_zzz_xxxyz_0, tg_zzz_xxxz_0, tg_zzz_xxxzz_0, \
                                     tg_zzz_xxyy_0, tg_zzz_xxyyy_0, tg_zzz_xxyyz_0, tg_zzz_xxyz_0, tg_zzz_xxyzz_0, \
                                     tg_zzz_xxzz_0, tg_zzz_xxzzz_0, tg_zzz_xyyy_0, tg_zzz_xyyyy_0, tg_zzz_xyyyz_0, \
                                     tg_zzz_xyyz_0, tg_zzz_xyyzz_0, tg_zzz_xyzz_0, tg_zzz_xyzzz_0, tg_zzz_xzzz_0, \
                                     tg_zzz_xzzzz_0, tg_zzz_yyyy_0, tg_zzz_yyyz_0, tg_zzz_yyzz_0, tg_zzz_yzzz_0, \
                                     tg_zzz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_xxzz_xxxx_0[j] = -ab_x * tg_xzz_xxxx_0[j] + tg_xzz_xxxxx_0[j];

                tg_xxzz_xxxy_0[j] = -ab_x * tg_xzz_xxxy_0[j] + tg_xzz_xxxxy_0[j];

                tg_xxzz_xxxz_0[j] = -ab_x * tg_xzz_xxxz_0[j] + tg_xzz_xxxxz_0[j];

                tg_xxzz_xxyy_0[j] = -ab_x * tg_xzz_xxyy_0[j] + tg_xzz_xxxyy_0[j];

                tg_xxzz_xxyz_0[j] = -ab_x * tg_xzz_xxyz_0[j] + tg_xzz_xxxyz_0[j];

                tg_xxzz_xxzz_0[j] = -ab_x * tg_xzz_xxzz_0[j] + tg_xzz_xxxzz_0[j];

                tg_xxzz_xyyy_0[j] = -ab_x * tg_xzz_xyyy_0[j] + tg_xzz_xxyyy_0[j];

                tg_xxzz_xyyz_0[j] = -ab_x * tg_xzz_xyyz_0[j] + tg_xzz_xxyyz_0[j];

                tg_xxzz_xyzz_0[j] = -ab_x * tg_xzz_xyzz_0[j] + tg_xzz_xxyzz_0[j];

                tg_xxzz_xzzz_0[j] = -ab_x * tg_xzz_xzzz_0[j] + tg_xzz_xxzzz_0[j];

                tg_xxzz_yyyy_0[j] = -ab_x * tg_xzz_yyyy_0[j] + tg_xzz_xyyyy_0[j];

                tg_xxzz_yyyz_0[j] = -ab_x * tg_xzz_yyyz_0[j] + tg_xzz_xyyyz_0[j];

                tg_xxzz_yyzz_0[j] = -ab_x * tg_xzz_yyzz_0[j] + tg_xzz_xyyzz_0[j];

                tg_xxzz_yzzz_0[j] = -ab_x * tg_xzz_yzzz_0[j] + tg_xzz_xyzzz_0[j];

                tg_xxzz_zzzz_0[j] = -ab_x * tg_xzz_zzzz_0[j] + tg_xzz_xzzzz_0[j];

                tg_xyyy_xxxx_0[j] = -ab_x * tg_yyy_xxxx_0[j] + tg_yyy_xxxxx_0[j];

                tg_xyyy_xxxy_0[j] = -ab_x * tg_yyy_xxxy_0[j] + tg_yyy_xxxxy_0[j];

                tg_xyyy_xxxz_0[j] = -ab_x * tg_yyy_xxxz_0[j] + tg_yyy_xxxxz_0[j];

                tg_xyyy_xxyy_0[j] = -ab_x * tg_yyy_xxyy_0[j] + tg_yyy_xxxyy_0[j];

                tg_xyyy_xxyz_0[j] = -ab_x * tg_yyy_xxyz_0[j] + tg_yyy_xxxyz_0[j];

                tg_xyyy_xxzz_0[j] = -ab_x * tg_yyy_xxzz_0[j] + tg_yyy_xxxzz_0[j];

                tg_xyyy_xyyy_0[j] = -ab_x * tg_yyy_xyyy_0[j] + tg_yyy_xxyyy_0[j];

                tg_xyyy_xyyz_0[j] = -ab_x * tg_yyy_xyyz_0[j] + tg_yyy_xxyyz_0[j];

                tg_xyyy_xyzz_0[j] = -ab_x * tg_yyy_xyzz_0[j] + tg_yyy_xxyzz_0[j];

                tg_xyyy_xzzz_0[j] = -ab_x * tg_yyy_xzzz_0[j] + tg_yyy_xxzzz_0[j];

                tg_xyyy_yyyy_0[j] = -ab_x * tg_yyy_yyyy_0[j] + tg_yyy_xyyyy_0[j];

                tg_xyyy_yyyz_0[j] = -ab_x * tg_yyy_yyyz_0[j] + tg_yyy_xyyyz_0[j];

                tg_xyyy_yyzz_0[j] = -ab_x * tg_yyy_yyzz_0[j] + tg_yyy_xyyzz_0[j];

                tg_xyyy_yzzz_0[j] = -ab_x * tg_yyy_yzzz_0[j] + tg_yyy_xyzzz_0[j];

                tg_xyyy_zzzz_0[j] = -ab_x * tg_yyy_zzzz_0[j] + tg_yyy_xzzzz_0[j];

                tg_xyyz_xxxx_0[j] = -ab_x * tg_yyz_xxxx_0[j] + tg_yyz_xxxxx_0[j];

                tg_xyyz_xxxy_0[j] = -ab_x * tg_yyz_xxxy_0[j] + tg_yyz_xxxxy_0[j];

                tg_xyyz_xxxz_0[j] = -ab_x * tg_yyz_xxxz_0[j] + tg_yyz_xxxxz_0[j];

                tg_xyyz_xxyy_0[j] = -ab_x * tg_yyz_xxyy_0[j] + tg_yyz_xxxyy_0[j];

                tg_xyyz_xxyz_0[j] = -ab_x * tg_yyz_xxyz_0[j] + tg_yyz_xxxyz_0[j];

                tg_xyyz_xxzz_0[j] = -ab_x * tg_yyz_xxzz_0[j] + tg_yyz_xxxzz_0[j];

                tg_xyyz_xyyy_0[j] = -ab_x * tg_yyz_xyyy_0[j] + tg_yyz_xxyyy_0[j];

                tg_xyyz_xyyz_0[j] = -ab_x * tg_yyz_xyyz_0[j] + tg_yyz_xxyyz_0[j];

                tg_xyyz_xyzz_0[j] = -ab_x * tg_yyz_xyzz_0[j] + tg_yyz_xxyzz_0[j];

                tg_xyyz_xzzz_0[j] = -ab_x * tg_yyz_xzzz_0[j] + tg_yyz_xxzzz_0[j];

                tg_xyyz_yyyy_0[j] = -ab_x * tg_yyz_yyyy_0[j] + tg_yyz_xyyyy_0[j];

                tg_xyyz_yyyz_0[j] = -ab_x * tg_yyz_yyyz_0[j] + tg_yyz_xyyyz_0[j];

                tg_xyyz_yyzz_0[j] = -ab_x * tg_yyz_yyzz_0[j] + tg_yyz_xyyzz_0[j];

                tg_xyyz_yzzz_0[j] = -ab_x * tg_yyz_yzzz_0[j] + tg_yyz_xyzzz_0[j];

                tg_xyyz_zzzz_0[j] = -ab_x * tg_yyz_zzzz_0[j] + tg_yyz_xzzzz_0[j];

                tg_xyzz_xxxx_0[j] = -ab_x * tg_yzz_xxxx_0[j] + tg_yzz_xxxxx_0[j];

                tg_xyzz_xxxy_0[j] = -ab_x * tg_yzz_xxxy_0[j] + tg_yzz_xxxxy_0[j];

                tg_xyzz_xxxz_0[j] = -ab_x * tg_yzz_xxxz_0[j] + tg_yzz_xxxxz_0[j];

                tg_xyzz_xxyy_0[j] = -ab_x * tg_yzz_xxyy_0[j] + tg_yzz_xxxyy_0[j];

                tg_xyzz_xxyz_0[j] = -ab_x * tg_yzz_xxyz_0[j] + tg_yzz_xxxyz_0[j];

                tg_xyzz_xxzz_0[j] = -ab_x * tg_yzz_xxzz_0[j] + tg_yzz_xxxzz_0[j];

                tg_xyzz_xyyy_0[j] = -ab_x * tg_yzz_xyyy_0[j] + tg_yzz_xxyyy_0[j];

                tg_xyzz_xyyz_0[j] = -ab_x * tg_yzz_xyyz_0[j] + tg_yzz_xxyyz_0[j];

                tg_xyzz_xyzz_0[j] = -ab_x * tg_yzz_xyzz_0[j] + tg_yzz_xxyzz_0[j];

                tg_xyzz_xzzz_0[j] = -ab_x * tg_yzz_xzzz_0[j] + tg_yzz_xxzzz_0[j];

                tg_xyzz_yyyy_0[j] = -ab_x * tg_yzz_yyyy_0[j] + tg_yzz_xyyyy_0[j];

                tg_xyzz_yyyz_0[j] = -ab_x * tg_yzz_yyyz_0[j] + tg_yzz_xyyyz_0[j];

                tg_xyzz_yyzz_0[j] = -ab_x * tg_yzz_yyzz_0[j] + tg_yzz_xyyzz_0[j];

                tg_xyzz_yzzz_0[j] = -ab_x * tg_yzz_yzzz_0[j] + tg_yzz_xyzzz_0[j];

                tg_xyzz_zzzz_0[j] = -ab_x * tg_yzz_zzzz_0[j] + tg_yzz_xzzzz_0[j];

                tg_xzzz_xxxx_0[j] = -ab_x * tg_zzz_xxxx_0[j] + tg_zzz_xxxxx_0[j];

                tg_xzzz_xxxy_0[j] = -ab_x * tg_zzz_xxxy_0[j] + tg_zzz_xxxxy_0[j];

                tg_xzzz_xxxz_0[j] = -ab_x * tg_zzz_xxxz_0[j] + tg_zzz_xxxxz_0[j];

                tg_xzzz_xxyy_0[j] = -ab_x * tg_zzz_xxyy_0[j] + tg_zzz_xxxyy_0[j];

                tg_xzzz_xxyz_0[j] = -ab_x * tg_zzz_xxyz_0[j] + tg_zzz_xxxyz_0[j];

                tg_xzzz_xxzz_0[j] = -ab_x * tg_zzz_xxzz_0[j] + tg_zzz_xxxzz_0[j];

                tg_xzzz_xyyy_0[j] = -ab_x * tg_zzz_xyyy_0[j] + tg_zzz_xxyyy_0[j];

                tg_xzzz_xyyz_0[j] = -ab_x * tg_zzz_xyyz_0[j] + tg_zzz_xxyyz_0[j];

                tg_xzzz_xyzz_0[j] = -ab_x * tg_zzz_xyzz_0[j] + tg_zzz_xxyzz_0[j];

                tg_xzzz_xzzz_0[j] = -ab_x * tg_zzz_xzzz_0[j] + tg_zzz_xxzzz_0[j];

                tg_xzzz_yyyy_0[j] = -ab_x * tg_zzz_yyyy_0[j] + tg_zzz_xyyyy_0[j];

                tg_xzzz_yyyz_0[j] = -ab_x * tg_zzz_yyyz_0[j] + tg_zzz_xyyyz_0[j];

                tg_xzzz_yyzz_0[j] = -ab_x * tg_zzz_yyzz_0[j] + tg_zzz_xyyzz_0[j];

                tg_xzzz_yzzz_0[j] = -ab_x * tg_zzz_yzzz_0[j] + tg_zzz_xyzzz_0[j];

                tg_xzzz_zzzz_0[j] = -ab_x * tg_zzz_zzzz_0[j] + tg_zzz_xzzzz_0[j];
            }
        }
    }

    void
    compElectronRepulsionForGGXY_150_225(      CMemBlock2D<double>& braBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& abDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetContrPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (150,225)

        // set up distances R(AB) = A - B

        auto ab_y = (abDistances.data(1))[iContrPair];

        auto ab_z = (abDistances.data(2))[iContrPair];

        // set up ket side loop over intergals

        auto angc = ketGtoPairsBlock.getBraAngularMomentum();

        auto angd = ketGtoPairsBlock.getKetAngularMomentum();

        // set up index of integral

        auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {4, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        // check if integral is needed in recursion expansion

        if (pidx_g_4_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 4, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                         {3, 5, -1, -1}, {angc, angd, -1, -1}, 
                                                         2, 2, 0));

        auto kcomp = angmom::to_SphericalComponents(angc, angd);

        for (int32_t i = 0; i < kcomp; i++)
        {
            // set up pointers to auxilary integrals

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

            auto tg_yyy_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 127 * kcomp + i); 

            auto tg_yyy_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 129 * kcomp + i); 

            auto tg_yyy_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 130 * kcomp + i); 

            auto tg_yyy_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 132 * kcomp + i); 

            auto tg_yyy_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 133 * kcomp + i); 

            auto tg_yyy_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 134 * kcomp + i); 

            auto tg_yyy_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 136 * kcomp + i); 

            auto tg_yyy_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 137 * kcomp + i); 

            auto tg_yyy_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 138 * kcomp + i); 

            auto tg_yyy_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 139 * kcomp + i); 

            auto tg_yyy_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 141 * kcomp + i); 

            auto tg_yyy_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 142 * kcomp + i); 

            auto tg_yyy_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 143 * kcomp + i); 

            auto tg_yyy_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 144 * kcomp + i); 

            auto tg_yyy_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 145 * kcomp + i); 

            auto tg_yyz_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 148 * kcomp + i); 

            auto tg_yyz_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 150 * kcomp + i); 

            auto tg_yyz_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 151 * kcomp + i); 

            auto tg_yyz_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 153 * kcomp + i); 

            auto tg_yyz_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 154 * kcomp + i); 

            auto tg_yyz_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 155 * kcomp + i); 

            auto tg_yyz_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 157 * kcomp + i); 

            auto tg_yyz_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 158 * kcomp + i); 

            auto tg_yyz_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 159 * kcomp + i); 

            auto tg_yyz_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 160 * kcomp + i); 

            auto tg_yyz_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 162 * kcomp + i); 

            auto tg_yyz_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 163 * kcomp + i); 

            auto tg_yyz_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 164 * kcomp + i); 

            auto tg_yyz_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 165 * kcomp + i); 

            auto tg_yyz_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 166 * kcomp + i); 

            auto tg_yzz_xxxxy_0 = braBuffer.data(pidx_g_3_5_m0 + 169 * kcomp + i); 

            auto tg_yzz_xxxyy_0 = braBuffer.data(pidx_g_3_5_m0 + 171 * kcomp + i); 

            auto tg_yzz_xxxyz_0 = braBuffer.data(pidx_g_3_5_m0 + 172 * kcomp + i); 

            auto tg_yzz_xxyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 174 * kcomp + i); 

            auto tg_yzz_xxyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 175 * kcomp + i); 

            auto tg_yzz_xxyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 176 * kcomp + i); 

            auto tg_yzz_xyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 178 * kcomp + i); 

            auto tg_yzz_xyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 179 * kcomp + i); 

            auto tg_yzz_xyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 180 * kcomp + i); 

            auto tg_yzz_xyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 181 * kcomp + i); 

            auto tg_yzz_yyyyy_0 = braBuffer.data(pidx_g_3_5_m0 + 183 * kcomp + i); 

            auto tg_yzz_yyyyz_0 = braBuffer.data(pidx_g_3_5_m0 + 184 * kcomp + i); 

            auto tg_yzz_yyyzz_0 = braBuffer.data(pidx_g_3_5_m0 + 185 * kcomp + i); 

            auto tg_yzz_yyzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 186 * kcomp + i); 

            auto tg_yzz_yzzzz_0 = braBuffer.data(pidx_g_3_5_m0 + 187 * kcomp + i); 

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

            // set up pointers to integrals

            auto tg_yyyy_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 150 * kcomp + i); 

            auto tg_yyyy_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 151 * kcomp + i); 

            auto tg_yyyy_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 152 * kcomp + i); 

            auto tg_yyyy_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 153 * kcomp + i); 

            auto tg_yyyy_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 154 * kcomp + i); 

            auto tg_yyyy_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 155 * kcomp + i); 

            auto tg_yyyy_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 156 * kcomp + i); 

            auto tg_yyyy_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 157 * kcomp + i); 

            auto tg_yyyy_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 158 * kcomp + i); 

            auto tg_yyyy_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 159 * kcomp + i); 

            auto tg_yyyy_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 160 * kcomp + i); 

            auto tg_yyyy_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 161 * kcomp + i); 

            auto tg_yyyy_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 162 * kcomp + i); 

            auto tg_yyyy_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 163 * kcomp + i); 

            auto tg_yyyy_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 164 * kcomp + i); 

            auto tg_yyyz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 165 * kcomp + i); 

            auto tg_yyyz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 166 * kcomp + i); 

            auto tg_yyyz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 167 * kcomp + i); 

            auto tg_yyyz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 168 * kcomp + i); 

            auto tg_yyyz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 169 * kcomp + i); 

            auto tg_yyyz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 170 * kcomp + i); 

            auto tg_yyyz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 171 * kcomp + i); 

            auto tg_yyyz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 172 * kcomp + i); 

            auto tg_yyyz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 173 * kcomp + i); 

            auto tg_yyyz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 174 * kcomp + i); 

            auto tg_yyyz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 175 * kcomp + i); 

            auto tg_yyyz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 176 * kcomp + i); 

            auto tg_yyyz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 177 * kcomp + i); 

            auto tg_yyyz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 178 * kcomp + i); 

            auto tg_yyyz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 179 * kcomp + i); 

            auto tg_yyzz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 180 * kcomp + i); 

            auto tg_yyzz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 181 * kcomp + i); 

            auto tg_yyzz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 182 * kcomp + i); 

            auto tg_yyzz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 183 * kcomp + i); 

            auto tg_yyzz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 184 * kcomp + i); 

            auto tg_yyzz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 185 * kcomp + i); 

            auto tg_yyzz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 186 * kcomp + i); 

            auto tg_yyzz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 187 * kcomp + i); 

            auto tg_yyzz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 188 * kcomp + i); 

            auto tg_yyzz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 189 * kcomp + i); 

            auto tg_yyzz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 190 * kcomp + i); 

            auto tg_yyzz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 191 * kcomp + i); 

            auto tg_yyzz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 192 * kcomp + i); 

            auto tg_yyzz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 193 * kcomp + i); 

            auto tg_yyzz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 194 * kcomp + i); 

            auto tg_yzzz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 195 * kcomp + i); 

            auto tg_yzzz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 196 * kcomp + i); 

            auto tg_yzzz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 197 * kcomp + i); 

            auto tg_yzzz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 198 * kcomp + i); 

            auto tg_yzzz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 199 * kcomp + i); 

            auto tg_yzzz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 200 * kcomp + i); 

            auto tg_yzzz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 201 * kcomp + i); 

            auto tg_yzzz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 202 * kcomp + i); 

            auto tg_yzzz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 203 * kcomp + i); 

            auto tg_yzzz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 204 * kcomp + i); 

            auto tg_yzzz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 205 * kcomp + i); 

            auto tg_yzzz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 206 * kcomp + i); 

            auto tg_yzzz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 207 * kcomp + i); 

            auto tg_yzzz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 208 * kcomp + i); 

            auto tg_yzzz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 209 * kcomp + i); 

            auto tg_zzzz_xxxx_0 = braBuffer.data(pidx_g_4_4_m0 + 210 * kcomp + i); 

            auto tg_zzzz_xxxy_0 = braBuffer.data(pidx_g_4_4_m0 + 211 * kcomp + i); 

            auto tg_zzzz_xxxz_0 = braBuffer.data(pidx_g_4_4_m0 + 212 * kcomp + i); 

            auto tg_zzzz_xxyy_0 = braBuffer.data(pidx_g_4_4_m0 + 213 * kcomp + i); 

            auto tg_zzzz_xxyz_0 = braBuffer.data(pidx_g_4_4_m0 + 214 * kcomp + i); 

            auto tg_zzzz_xxzz_0 = braBuffer.data(pidx_g_4_4_m0 + 215 * kcomp + i); 

            auto tg_zzzz_xyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 216 * kcomp + i); 

            auto tg_zzzz_xyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 217 * kcomp + i); 

            auto tg_zzzz_xyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 218 * kcomp + i); 

            auto tg_zzzz_xzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 219 * kcomp + i); 

            auto tg_zzzz_yyyy_0 = braBuffer.data(pidx_g_4_4_m0 + 220 * kcomp + i); 

            auto tg_zzzz_yyyz_0 = braBuffer.data(pidx_g_4_4_m0 + 221 * kcomp + i); 

            auto tg_zzzz_yyzz_0 = braBuffer.data(pidx_g_4_4_m0 + 222 * kcomp + i); 

            auto tg_zzzz_yzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 223 * kcomp + i); 

            auto tg_zzzz_zzzz_0 = braBuffer.data(pidx_g_4_4_m0 + 224 * kcomp + i); 

            // Batch of Integrals (150,225)

            #pragma omp simd aligned(tg_yyy_xxxx_0, tg_yyy_xxxxy_0, tg_yyy_xxxy_0, tg_yyy_xxxyy_0, \
                                     tg_yyy_xxxyz_0, tg_yyy_xxxz_0, tg_yyy_xxyy_0, tg_yyy_xxyyy_0, tg_yyy_xxyyz_0, \
                                     tg_yyy_xxyz_0, tg_yyy_xxyzz_0, tg_yyy_xxzz_0, tg_yyy_xyyy_0, tg_yyy_xyyyy_0, \
                                     tg_yyy_xyyyz_0, tg_yyy_xyyz_0, tg_yyy_xyyzz_0, tg_yyy_xyzz_0, tg_yyy_xyzzz_0, \
                                     tg_yyy_xzzz_0, tg_yyy_yyyy_0, tg_yyy_yyyyy_0, tg_yyy_yyyyz_0, tg_yyy_yyyz_0, \
                                     tg_yyy_yyyzz_0, tg_yyy_yyzz_0, tg_yyy_yyzzz_0, tg_yyy_yzzz_0, tg_yyy_yzzzz_0, \
                                     tg_yyy_zzzz_0, tg_yyyy_xxxx_0, tg_yyyy_xxxy_0, tg_yyyy_xxxz_0, tg_yyyy_xxyy_0, \
                                     tg_yyyy_xxyz_0, tg_yyyy_xxzz_0, tg_yyyy_xyyy_0, tg_yyyy_xyyz_0, tg_yyyy_xyzz_0, \
                                     tg_yyyy_xzzz_0, tg_yyyy_yyyy_0, tg_yyyy_yyyz_0, tg_yyyy_yyzz_0, tg_yyyy_yzzz_0, \
                                     tg_yyyy_zzzz_0, tg_yyyz_xxxx_0, tg_yyyz_xxxy_0, tg_yyyz_xxxz_0, tg_yyyz_xxyy_0, \
                                     tg_yyyz_xxyz_0, tg_yyyz_xxzz_0, tg_yyyz_xyyy_0, tg_yyyz_xyyz_0, tg_yyyz_xyzz_0, \
                                     tg_yyyz_xzzz_0, tg_yyyz_yyyy_0, tg_yyyz_yyyz_0, tg_yyyz_yyzz_0, tg_yyyz_yzzz_0, \
                                     tg_yyyz_zzzz_0, tg_yyz_xxxx_0, tg_yyz_xxxxy_0, tg_yyz_xxxy_0, tg_yyz_xxxyy_0, \
                                     tg_yyz_xxxyz_0, tg_yyz_xxxz_0, tg_yyz_xxyy_0, tg_yyz_xxyyy_0, tg_yyz_xxyyz_0, \
                                     tg_yyz_xxyz_0, tg_yyz_xxyzz_0, tg_yyz_xxzz_0, tg_yyz_xyyy_0, tg_yyz_xyyyy_0, \
                                     tg_yyz_xyyyz_0, tg_yyz_xyyz_0, tg_yyz_xyyzz_0, tg_yyz_xyzz_0, tg_yyz_xyzzz_0, \
                                     tg_yyz_xzzz_0, tg_yyz_yyyy_0, tg_yyz_yyyyy_0, tg_yyz_yyyyz_0, tg_yyz_yyyz_0, \
                                     tg_yyz_yyyzz_0, tg_yyz_yyzz_0, tg_yyz_yyzzz_0, tg_yyz_yzzz_0, tg_yyz_yzzzz_0, \
                                     tg_yyz_zzzz_0, tg_yyzz_xxxx_0, tg_yyzz_xxxy_0, tg_yyzz_xxxz_0, tg_yyzz_xxyy_0, \
                                     tg_yyzz_xxyz_0, tg_yyzz_xxzz_0, tg_yyzz_xyyy_0, tg_yyzz_xyyz_0, tg_yyzz_xyzz_0, \
                                     tg_yyzz_xzzz_0, tg_yyzz_yyyy_0, tg_yyzz_yyyz_0, tg_yyzz_yyzz_0, tg_yyzz_yzzz_0, \
                                     tg_yyzz_zzzz_0, tg_yzz_xxxx_0, tg_yzz_xxxxy_0, tg_yzz_xxxy_0, tg_yzz_xxxyy_0, \
                                     tg_yzz_xxxyz_0, tg_yzz_xxxz_0, tg_yzz_xxyy_0, tg_yzz_xxyyy_0, tg_yzz_xxyyz_0, \
                                     tg_yzz_xxyz_0, tg_yzz_xxyzz_0, tg_yzz_xxzz_0, tg_yzz_xyyy_0, tg_yzz_xyyyy_0, \
                                     tg_yzz_xyyyz_0, tg_yzz_xyyz_0, tg_yzz_xyyzz_0, tg_yzz_xyzz_0, tg_yzz_xyzzz_0, \
                                     tg_yzz_xzzz_0, tg_yzz_yyyy_0, tg_yzz_yyyyy_0, tg_yzz_yyyyz_0, tg_yzz_yyyz_0, \
                                     tg_yzz_yyyzz_0, tg_yzz_yyzz_0, tg_yzz_yyzzz_0, tg_yzz_yzzz_0, tg_yzz_yzzzz_0, \
                                     tg_yzz_zzzz_0, tg_yzzz_xxxx_0, tg_yzzz_xxxy_0, tg_yzzz_xxxz_0, tg_yzzz_xxyy_0, \
                                     tg_yzzz_xxyz_0, tg_yzzz_xxzz_0, tg_yzzz_xyyy_0, tg_yzzz_xyyz_0, tg_yzzz_xyzz_0, \
                                     tg_yzzz_xzzz_0, tg_yzzz_yyyy_0, tg_yzzz_yyyz_0, tg_yzzz_yyzz_0, tg_yzzz_yzzz_0, \
                                     tg_yzzz_zzzz_0, tg_zzz_xxxx_0, tg_zzz_xxxxy_0, tg_zzz_xxxxz_0, tg_zzz_xxxy_0, \
                                     tg_zzz_xxxyy_0, tg_zzz_xxxyz_0, tg_zzz_xxxz_0, tg_zzz_xxxzz_0, tg_zzz_xxyy_0, \
                                     tg_zzz_xxyyy_0, tg_zzz_xxyyz_0, tg_zzz_xxyz_0, tg_zzz_xxyzz_0, tg_zzz_xxzz_0, \
                                     tg_zzz_xxzzz_0, tg_zzz_xyyy_0, tg_zzz_xyyyy_0, tg_zzz_xyyyz_0, tg_zzz_xyyz_0, \
                                     tg_zzz_xyyzz_0, tg_zzz_xyzz_0, tg_zzz_xyzzz_0, tg_zzz_xzzz_0, tg_zzz_xzzzz_0, \
                                     tg_zzz_yyyy_0, tg_zzz_yyyyy_0, tg_zzz_yyyyz_0, tg_zzz_yyyz_0, tg_zzz_yyyzz_0, \
                                     tg_zzz_yyzz_0, tg_zzz_yyzzz_0, tg_zzz_yzzz_0, tg_zzz_yzzzz_0, tg_zzz_zzzz_0, \
                                     tg_zzz_zzzzz_0, tg_zzzz_xxxx_0, tg_zzzz_xxxy_0, tg_zzzz_xxxz_0, tg_zzzz_xxyy_0, \
                                     tg_zzzz_xxyz_0, tg_zzzz_xxzz_0, tg_zzzz_xyyy_0, tg_zzzz_xyyz_0, tg_zzzz_xyzz_0, \
                                     tg_zzzz_xzzz_0, tg_zzzz_yyyy_0, tg_zzzz_yyyz_0, tg_zzzz_yyzz_0, tg_zzzz_yzzz_0, \
                                     tg_zzzz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nKetContrPairs; j++)
            {
                tg_yyyy_xxxx_0[j] = -ab_y * tg_yyy_xxxx_0[j] + tg_yyy_xxxxy_0[j];

                tg_yyyy_xxxy_0[j] = -ab_y * tg_yyy_xxxy_0[j] + tg_yyy_xxxyy_0[j];

                tg_yyyy_xxxz_0[j] = -ab_y * tg_yyy_xxxz_0[j] + tg_yyy_xxxyz_0[j];

                tg_yyyy_xxyy_0[j] = -ab_y * tg_yyy_xxyy_0[j] + tg_yyy_xxyyy_0[j];

                tg_yyyy_xxyz_0[j] = -ab_y * tg_yyy_xxyz_0[j] + tg_yyy_xxyyz_0[j];

                tg_yyyy_xxzz_0[j] = -ab_y * tg_yyy_xxzz_0[j] + tg_yyy_xxyzz_0[j];

                tg_yyyy_xyyy_0[j] = -ab_y * tg_yyy_xyyy_0[j] + tg_yyy_xyyyy_0[j];

                tg_yyyy_xyyz_0[j] = -ab_y * tg_yyy_xyyz_0[j] + tg_yyy_xyyyz_0[j];

                tg_yyyy_xyzz_0[j] = -ab_y * tg_yyy_xyzz_0[j] + tg_yyy_xyyzz_0[j];

                tg_yyyy_xzzz_0[j] = -ab_y * tg_yyy_xzzz_0[j] + tg_yyy_xyzzz_0[j];

                tg_yyyy_yyyy_0[j] = -ab_y * tg_yyy_yyyy_0[j] + tg_yyy_yyyyy_0[j];

                tg_yyyy_yyyz_0[j] = -ab_y * tg_yyy_yyyz_0[j] + tg_yyy_yyyyz_0[j];

                tg_yyyy_yyzz_0[j] = -ab_y * tg_yyy_yyzz_0[j] + tg_yyy_yyyzz_0[j];

                tg_yyyy_yzzz_0[j] = -ab_y * tg_yyy_yzzz_0[j] + tg_yyy_yyzzz_0[j];

                tg_yyyy_zzzz_0[j] = -ab_y * tg_yyy_zzzz_0[j] + tg_yyy_yzzzz_0[j];

                tg_yyyz_xxxx_0[j] = -ab_y * tg_yyz_xxxx_0[j] + tg_yyz_xxxxy_0[j];

                tg_yyyz_xxxy_0[j] = -ab_y * tg_yyz_xxxy_0[j] + tg_yyz_xxxyy_0[j];

                tg_yyyz_xxxz_0[j] = -ab_y * tg_yyz_xxxz_0[j] + tg_yyz_xxxyz_0[j];

                tg_yyyz_xxyy_0[j] = -ab_y * tg_yyz_xxyy_0[j] + tg_yyz_xxyyy_0[j];

                tg_yyyz_xxyz_0[j] = -ab_y * tg_yyz_xxyz_0[j] + tg_yyz_xxyyz_0[j];

                tg_yyyz_xxzz_0[j] = -ab_y * tg_yyz_xxzz_0[j] + tg_yyz_xxyzz_0[j];

                tg_yyyz_xyyy_0[j] = -ab_y * tg_yyz_xyyy_0[j] + tg_yyz_xyyyy_0[j];

                tg_yyyz_xyyz_0[j] = -ab_y * tg_yyz_xyyz_0[j] + tg_yyz_xyyyz_0[j];

                tg_yyyz_xyzz_0[j] = -ab_y * tg_yyz_xyzz_0[j] + tg_yyz_xyyzz_0[j];

                tg_yyyz_xzzz_0[j] = -ab_y * tg_yyz_xzzz_0[j] + tg_yyz_xyzzz_0[j];

                tg_yyyz_yyyy_0[j] = -ab_y * tg_yyz_yyyy_0[j] + tg_yyz_yyyyy_0[j];

                tg_yyyz_yyyz_0[j] = -ab_y * tg_yyz_yyyz_0[j] + tg_yyz_yyyyz_0[j];

                tg_yyyz_yyzz_0[j] = -ab_y * tg_yyz_yyzz_0[j] + tg_yyz_yyyzz_0[j];

                tg_yyyz_yzzz_0[j] = -ab_y * tg_yyz_yzzz_0[j] + tg_yyz_yyzzz_0[j];

                tg_yyyz_zzzz_0[j] = -ab_y * tg_yyz_zzzz_0[j] + tg_yyz_yzzzz_0[j];

                tg_yyzz_xxxx_0[j] = -ab_y * tg_yzz_xxxx_0[j] + tg_yzz_xxxxy_0[j];

                tg_yyzz_xxxy_0[j] = -ab_y * tg_yzz_xxxy_0[j] + tg_yzz_xxxyy_0[j];

                tg_yyzz_xxxz_0[j] = -ab_y * tg_yzz_xxxz_0[j] + tg_yzz_xxxyz_0[j];

                tg_yyzz_xxyy_0[j] = -ab_y * tg_yzz_xxyy_0[j] + tg_yzz_xxyyy_0[j];

                tg_yyzz_xxyz_0[j] = -ab_y * tg_yzz_xxyz_0[j] + tg_yzz_xxyyz_0[j];

                tg_yyzz_xxzz_0[j] = -ab_y * tg_yzz_xxzz_0[j] + tg_yzz_xxyzz_0[j];

                tg_yyzz_xyyy_0[j] = -ab_y * tg_yzz_xyyy_0[j] + tg_yzz_xyyyy_0[j];

                tg_yyzz_xyyz_0[j] = -ab_y * tg_yzz_xyyz_0[j] + tg_yzz_xyyyz_0[j];

                tg_yyzz_xyzz_0[j] = -ab_y * tg_yzz_xyzz_0[j] + tg_yzz_xyyzz_0[j];

                tg_yyzz_xzzz_0[j] = -ab_y * tg_yzz_xzzz_0[j] + tg_yzz_xyzzz_0[j];

                tg_yyzz_yyyy_0[j] = -ab_y * tg_yzz_yyyy_0[j] + tg_yzz_yyyyy_0[j];

                tg_yyzz_yyyz_0[j] = -ab_y * tg_yzz_yyyz_0[j] + tg_yzz_yyyyz_0[j];

                tg_yyzz_yyzz_0[j] = -ab_y * tg_yzz_yyzz_0[j] + tg_yzz_yyyzz_0[j];

                tg_yyzz_yzzz_0[j] = -ab_y * tg_yzz_yzzz_0[j] + tg_yzz_yyzzz_0[j];

                tg_yyzz_zzzz_0[j] = -ab_y * tg_yzz_zzzz_0[j] + tg_yzz_yzzzz_0[j];

                tg_yzzz_xxxx_0[j] = -ab_y * tg_zzz_xxxx_0[j] + tg_zzz_xxxxy_0[j];

                tg_yzzz_xxxy_0[j] = -ab_y * tg_zzz_xxxy_0[j] + tg_zzz_xxxyy_0[j];

                tg_yzzz_xxxz_0[j] = -ab_y * tg_zzz_xxxz_0[j] + tg_zzz_xxxyz_0[j];

                tg_yzzz_xxyy_0[j] = -ab_y * tg_zzz_xxyy_0[j] + tg_zzz_xxyyy_0[j];

                tg_yzzz_xxyz_0[j] = -ab_y * tg_zzz_xxyz_0[j] + tg_zzz_xxyyz_0[j];

                tg_yzzz_xxzz_0[j] = -ab_y * tg_zzz_xxzz_0[j] + tg_zzz_xxyzz_0[j];

                tg_yzzz_xyyy_0[j] = -ab_y * tg_zzz_xyyy_0[j] + tg_zzz_xyyyy_0[j];

                tg_yzzz_xyyz_0[j] = -ab_y * tg_zzz_xyyz_0[j] + tg_zzz_xyyyz_0[j];

                tg_yzzz_xyzz_0[j] = -ab_y * tg_zzz_xyzz_0[j] + tg_zzz_xyyzz_0[j];

                tg_yzzz_xzzz_0[j] = -ab_y * tg_zzz_xzzz_0[j] + tg_zzz_xyzzz_0[j];

                tg_yzzz_yyyy_0[j] = -ab_y * tg_zzz_yyyy_0[j] + tg_zzz_yyyyy_0[j];

                tg_yzzz_yyyz_0[j] = -ab_y * tg_zzz_yyyz_0[j] + tg_zzz_yyyyz_0[j];

                tg_yzzz_yyzz_0[j] = -ab_y * tg_zzz_yyzz_0[j] + tg_zzz_yyyzz_0[j];

                tg_yzzz_yzzz_0[j] = -ab_y * tg_zzz_yzzz_0[j] + tg_zzz_yyzzz_0[j];

                tg_yzzz_zzzz_0[j] = -ab_y * tg_zzz_zzzz_0[j] + tg_zzz_yzzzz_0[j];

                tg_zzzz_xxxx_0[j] = -ab_z * tg_zzz_xxxx_0[j] + tg_zzz_xxxxz_0[j];

                tg_zzzz_xxxy_0[j] = -ab_z * tg_zzz_xxxy_0[j] + tg_zzz_xxxyz_0[j];

                tg_zzzz_xxxz_0[j] = -ab_z * tg_zzz_xxxz_0[j] + tg_zzz_xxxzz_0[j];

                tg_zzzz_xxyy_0[j] = -ab_z * tg_zzz_xxyy_0[j] + tg_zzz_xxyyz_0[j];

                tg_zzzz_xxyz_0[j] = -ab_z * tg_zzz_xxyz_0[j] + tg_zzz_xxyzz_0[j];

                tg_zzzz_xxzz_0[j] = -ab_z * tg_zzz_xxzz_0[j] + tg_zzz_xxzzz_0[j];

                tg_zzzz_xyyy_0[j] = -ab_z * tg_zzz_xyyy_0[j] + tg_zzz_xyyyz_0[j];

                tg_zzzz_xyyz_0[j] = -ab_z * tg_zzz_xyyz_0[j] + tg_zzz_xyyzz_0[j];

                tg_zzzz_xyzz_0[j] = -ab_z * tg_zzz_xyzz_0[j] + tg_zzz_xyzzz_0[j];

                tg_zzzz_xzzz_0[j] = -ab_z * tg_zzz_xzzz_0[j] + tg_zzz_xzzzz_0[j];

                tg_zzzz_yyyy_0[j] = -ab_z * tg_zzz_yyyy_0[j] + tg_zzz_yyyyz_0[j];

                tg_zzzz_yyyz_0[j] = -ab_z * tg_zzz_yyyz_0[j] + tg_zzz_yyyzz_0[j];

                tg_zzzz_yyzz_0[j] = -ab_z * tg_zzz_yyzz_0[j] + tg_zzz_yyzzz_0[j];

                tg_zzzz_yzzz_0[j] = -ab_z * tg_zzz_yzzz_0[j] + tg_zzz_yzzzz_0[j];

                tg_zzzz_zzzz_0[j] = -ab_z * tg_zzz_zzzz_0[j] + tg_zzz_zzzzz_0[j];
            }
        }
    }


} // eribrrfunc namespace

