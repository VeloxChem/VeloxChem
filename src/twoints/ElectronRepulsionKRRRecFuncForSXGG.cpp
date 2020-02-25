//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionKRRRecFuncForSXGG.hpp"

#include "AngularMomentum.hpp"

namespace erikrrfunc { // erikrrfunc namespace

    void
    compElectronRepulsionForSXGG(      CMemBlock2D<double>& ketBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& cdDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetContrPairs,
                                 const int32_t              iContrPair)
    {
        erikrrfunc::compElectronRepulsionForSXGG_0_75(ketBuffer,
                                                      recursionMap,
                                                      cdDistances,
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetContrPairs,
                                                      iContrPair); 

        erikrrfunc::compElectronRepulsionForSXGG_75_150(ketBuffer,
                                                        recursionMap,
                                                        cdDistances,
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetContrPairs,
                                                        iContrPair); 

        erikrrfunc::compElectronRepulsionForSXGG_150_225(ketBuffer,
                                                         recursionMap,
                                                         cdDistances,
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetContrPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSXGG_0_75(      CMemBlock2D<double>& ketBuffer,
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

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {4, 4, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 4, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 5, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

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

                auto tg_xyy_xxxxx_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 63); 

                auto tg_xyy_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 64); 

                auto tg_xyy_xxxxz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 65); 

                auto tg_xyy_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 66); 

                auto tg_xyy_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 67); 

                auto tg_xyy_xxxzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 68); 

                auto tg_xyy_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 69); 

                auto tg_xyy_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 70); 

                auto tg_xyy_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 71); 

                auto tg_xyy_xxzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 72); 

                auto tg_xyy_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 73); 

                auto tg_xyy_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 74); 

                auto tg_xyy_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 75); 

                auto tg_xyy_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 76); 

                auto tg_xyy_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 77); 

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

                // set up pointers to integrals

                auto tg_xxxx_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i); 

                auto tg_xxxx_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 1); 

                auto tg_xxxx_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 2); 

                auto tg_xxxx_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 3); 

                auto tg_xxxx_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 4); 

                auto tg_xxxx_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 5); 

                auto tg_xxxx_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 6); 

                auto tg_xxxx_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 7); 

                auto tg_xxxx_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 8); 

                auto tg_xxxx_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 9); 

                auto tg_xxxx_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 10); 

                auto tg_xxxx_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 11); 

                auto tg_xxxx_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 12); 

                auto tg_xxxx_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 13); 

                auto tg_xxxx_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 14); 

                auto tg_xxxy_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 15); 

                auto tg_xxxy_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 16); 

                auto tg_xxxy_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 17); 

                auto tg_xxxy_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 18); 

                auto tg_xxxy_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 19); 

                auto tg_xxxy_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 20); 

                auto tg_xxxy_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 21); 

                auto tg_xxxy_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 22); 

                auto tg_xxxy_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 23); 

                auto tg_xxxy_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 24); 

                auto tg_xxxy_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 25); 

                auto tg_xxxy_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 26); 

                auto tg_xxxy_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 27); 

                auto tg_xxxy_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 28); 

                auto tg_xxxy_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 29); 

                auto tg_xxxz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 30); 

                auto tg_xxxz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 31); 

                auto tg_xxxz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 32); 

                auto tg_xxxz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 33); 

                auto tg_xxxz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 34); 

                auto tg_xxxz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 35); 

                auto tg_xxxz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 36); 

                auto tg_xxxz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 37); 

                auto tg_xxxz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 38); 

                auto tg_xxxz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 39); 

                auto tg_xxxz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 40); 

                auto tg_xxxz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 41); 

                auto tg_xxxz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 42); 

                auto tg_xxxz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 43); 

                auto tg_xxxz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 44); 

                auto tg_xxyy_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 45); 

                auto tg_xxyy_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 46); 

                auto tg_xxyy_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 47); 

                auto tg_xxyy_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 48); 

                auto tg_xxyy_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 49); 

                auto tg_xxyy_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 50); 

                auto tg_xxyy_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 51); 

                auto tg_xxyy_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 52); 

                auto tg_xxyy_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 53); 

                auto tg_xxyy_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 54); 

                auto tg_xxyy_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 55); 

                auto tg_xxyy_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 56); 

                auto tg_xxyy_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 57); 

                auto tg_xxyy_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 58); 

                auto tg_xxyy_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 59); 

                auto tg_xxyz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 60); 

                auto tg_xxyz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 61); 

                auto tg_xxyz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 62); 

                auto tg_xxyz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 63); 

                auto tg_xxyz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 64); 

                auto tg_xxyz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 65); 

                auto tg_xxyz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 66); 

                auto tg_xxyz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 67); 

                auto tg_xxyz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 68); 

                auto tg_xxyz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 69); 

                auto tg_xxyz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 70); 

                auto tg_xxyz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 71); 

                auto tg_xxyz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 72); 

                auto tg_xxyz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 73); 

                auto tg_xxyz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 74); 

                // Batch of Integrals (0,75)

                #pragma omp simd aligned(cd_x, tg_xxx_xxxx_0, tg_xxx_xxxxx_0, tg_xxx_xxxxy_0, tg_xxx_xxxxz_0, \
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
                    double fr = -cd_x[j]; 

                    tg_xxxx_xxxx_0[j] = fr * tg_xxx_xxxx_0[j] + tg_xxx_xxxxx_0[j];

                    tg_xxxx_xxxy_0[j] = fr * tg_xxx_xxxy_0[j] + tg_xxx_xxxxy_0[j];

                    tg_xxxx_xxxz_0[j] = fr * tg_xxx_xxxz_0[j] + tg_xxx_xxxxz_0[j];

                    tg_xxxx_xxyy_0[j] = fr * tg_xxx_xxyy_0[j] + tg_xxx_xxxyy_0[j];

                    tg_xxxx_xxyz_0[j] = fr * tg_xxx_xxyz_0[j] + tg_xxx_xxxyz_0[j];

                    tg_xxxx_xxzz_0[j] = fr * tg_xxx_xxzz_0[j] + tg_xxx_xxxzz_0[j];

                    tg_xxxx_xyyy_0[j] = fr * tg_xxx_xyyy_0[j] + tg_xxx_xxyyy_0[j];

                    tg_xxxx_xyyz_0[j] = fr * tg_xxx_xyyz_0[j] + tg_xxx_xxyyz_0[j];

                    tg_xxxx_xyzz_0[j] = fr * tg_xxx_xyzz_0[j] + tg_xxx_xxyzz_0[j];

                    tg_xxxx_xzzz_0[j] = fr * tg_xxx_xzzz_0[j] + tg_xxx_xxzzz_0[j];

                    tg_xxxx_yyyy_0[j] = fr * tg_xxx_yyyy_0[j] + tg_xxx_xyyyy_0[j];

                    tg_xxxx_yyyz_0[j] = fr * tg_xxx_yyyz_0[j] + tg_xxx_xyyyz_0[j];

                    tg_xxxx_yyzz_0[j] = fr * tg_xxx_yyzz_0[j] + tg_xxx_xyyzz_0[j];

                    tg_xxxx_yzzz_0[j] = fr * tg_xxx_yzzz_0[j] + tg_xxx_xyzzz_0[j];

                    tg_xxxx_zzzz_0[j] = fr * tg_xxx_zzzz_0[j] + tg_xxx_xzzzz_0[j];

                    tg_xxxy_xxxx_0[j] = fr * tg_xxy_xxxx_0[j] + tg_xxy_xxxxx_0[j];

                    tg_xxxy_xxxy_0[j] = fr * tg_xxy_xxxy_0[j] + tg_xxy_xxxxy_0[j];

                    tg_xxxy_xxxz_0[j] = fr * tg_xxy_xxxz_0[j] + tg_xxy_xxxxz_0[j];

                    tg_xxxy_xxyy_0[j] = fr * tg_xxy_xxyy_0[j] + tg_xxy_xxxyy_0[j];

                    tg_xxxy_xxyz_0[j] = fr * tg_xxy_xxyz_0[j] + tg_xxy_xxxyz_0[j];

                    tg_xxxy_xxzz_0[j] = fr * tg_xxy_xxzz_0[j] + tg_xxy_xxxzz_0[j];

                    tg_xxxy_xyyy_0[j] = fr * tg_xxy_xyyy_0[j] + tg_xxy_xxyyy_0[j];

                    tg_xxxy_xyyz_0[j] = fr * tg_xxy_xyyz_0[j] + tg_xxy_xxyyz_0[j];

                    tg_xxxy_xyzz_0[j] = fr * tg_xxy_xyzz_0[j] + tg_xxy_xxyzz_0[j];

                    tg_xxxy_xzzz_0[j] = fr * tg_xxy_xzzz_0[j] + tg_xxy_xxzzz_0[j];

                    tg_xxxy_yyyy_0[j] = fr * tg_xxy_yyyy_0[j] + tg_xxy_xyyyy_0[j];

                    tg_xxxy_yyyz_0[j] = fr * tg_xxy_yyyz_0[j] + tg_xxy_xyyyz_0[j];

                    tg_xxxy_yyzz_0[j] = fr * tg_xxy_yyzz_0[j] + tg_xxy_xyyzz_0[j];

                    tg_xxxy_yzzz_0[j] = fr * tg_xxy_yzzz_0[j] + tg_xxy_xyzzz_0[j];

                    tg_xxxy_zzzz_0[j] = fr * tg_xxy_zzzz_0[j] + tg_xxy_xzzzz_0[j];

                    tg_xxxz_xxxx_0[j] = fr * tg_xxz_xxxx_0[j] + tg_xxz_xxxxx_0[j];

                    tg_xxxz_xxxy_0[j] = fr * tg_xxz_xxxy_0[j] + tg_xxz_xxxxy_0[j];

                    tg_xxxz_xxxz_0[j] = fr * tg_xxz_xxxz_0[j] + tg_xxz_xxxxz_0[j];

                    tg_xxxz_xxyy_0[j] = fr * tg_xxz_xxyy_0[j] + tg_xxz_xxxyy_0[j];

                    tg_xxxz_xxyz_0[j] = fr * tg_xxz_xxyz_0[j] + tg_xxz_xxxyz_0[j];

                    tg_xxxz_xxzz_0[j] = fr * tg_xxz_xxzz_0[j] + tg_xxz_xxxzz_0[j];

                    tg_xxxz_xyyy_0[j] = fr * tg_xxz_xyyy_0[j] + tg_xxz_xxyyy_0[j];

                    tg_xxxz_xyyz_0[j] = fr * tg_xxz_xyyz_0[j] + tg_xxz_xxyyz_0[j];

                    tg_xxxz_xyzz_0[j] = fr * tg_xxz_xyzz_0[j] + tg_xxz_xxyzz_0[j];

                    tg_xxxz_xzzz_0[j] = fr * tg_xxz_xzzz_0[j] + tg_xxz_xxzzz_0[j];

                    tg_xxxz_yyyy_0[j] = fr * tg_xxz_yyyy_0[j] + tg_xxz_xyyyy_0[j];

                    tg_xxxz_yyyz_0[j] = fr * tg_xxz_yyyz_0[j] + tg_xxz_xyyyz_0[j];

                    tg_xxxz_yyzz_0[j] = fr * tg_xxz_yyzz_0[j] + tg_xxz_xyyzz_0[j];

                    tg_xxxz_yzzz_0[j] = fr * tg_xxz_yzzz_0[j] + tg_xxz_xyzzz_0[j];

                    tg_xxxz_zzzz_0[j] = fr * tg_xxz_zzzz_0[j] + tg_xxz_xzzzz_0[j];

                    tg_xxyy_xxxx_0[j] = fr * tg_xyy_xxxx_0[j] + tg_xyy_xxxxx_0[j];

                    tg_xxyy_xxxy_0[j] = fr * tg_xyy_xxxy_0[j] + tg_xyy_xxxxy_0[j];

                    tg_xxyy_xxxz_0[j] = fr * tg_xyy_xxxz_0[j] + tg_xyy_xxxxz_0[j];

                    tg_xxyy_xxyy_0[j] = fr * tg_xyy_xxyy_0[j] + tg_xyy_xxxyy_0[j];

                    tg_xxyy_xxyz_0[j] = fr * tg_xyy_xxyz_0[j] + tg_xyy_xxxyz_0[j];

                    tg_xxyy_xxzz_0[j] = fr * tg_xyy_xxzz_0[j] + tg_xyy_xxxzz_0[j];

                    tg_xxyy_xyyy_0[j] = fr * tg_xyy_xyyy_0[j] + tg_xyy_xxyyy_0[j];

                    tg_xxyy_xyyz_0[j] = fr * tg_xyy_xyyz_0[j] + tg_xyy_xxyyz_0[j];

                    tg_xxyy_xyzz_0[j] = fr * tg_xyy_xyzz_0[j] + tg_xyy_xxyzz_0[j];

                    tg_xxyy_xzzz_0[j] = fr * tg_xyy_xzzz_0[j] + tg_xyy_xxzzz_0[j];

                    tg_xxyy_yyyy_0[j] = fr * tg_xyy_yyyy_0[j] + tg_xyy_xyyyy_0[j];

                    tg_xxyy_yyyz_0[j] = fr * tg_xyy_yyyz_0[j] + tg_xyy_xyyyz_0[j];

                    tg_xxyy_yyzz_0[j] = fr * tg_xyy_yyzz_0[j] + tg_xyy_xyyzz_0[j];

                    tg_xxyy_yzzz_0[j] = fr * tg_xyy_yzzz_0[j] + tg_xyy_xyzzz_0[j];

                    tg_xxyy_zzzz_0[j] = fr * tg_xyy_zzzz_0[j] + tg_xyy_xzzzz_0[j];

                    tg_xxyz_xxxx_0[j] = fr * tg_xyz_xxxx_0[j] + tg_xyz_xxxxx_0[j];

                    tg_xxyz_xxxy_0[j] = fr * tg_xyz_xxxy_0[j] + tg_xyz_xxxxy_0[j];

                    tg_xxyz_xxxz_0[j] = fr * tg_xyz_xxxz_0[j] + tg_xyz_xxxxz_0[j];

                    tg_xxyz_xxyy_0[j] = fr * tg_xyz_xxyy_0[j] + tg_xyz_xxxyy_0[j];

                    tg_xxyz_xxyz_0[j] = fr * tg_xyz_xxyz_0[j] + tg_xyz_xxxyz_0[j];

                    tg_xxyz_xxzz_0[j] = fr * tg_xyz_xxzz_0[j] + tg_xyz_xxxzz_0[j];

                    tg_xxyz_xyyy_0[j] = fr * tg_xyz_xyyy_0[j] + tg_xyz_xxyyy_0[j];

                    tg_xxyz_xyyz_0[j] = fr * tg_xyz_xyyz_0[j] + tg_xyz_xxyyz_0[j];

                    tg_xxyz_xyzz_0[j] = fr * tg_xyz_xyzz_0[j] + tg_xyz_xxyzz_0[j];

                    tg_xxyz_xzzz_0[j] = fr * tg_xyz_xzzz_0[j] + tg_xyz_xxzzz_0[j];

                    tg_xxyz_yyyy_0[j] = fr * tg_xyz_yyyy_0[j] + tg_xyz_xyyyy_0[j];

                    tg_xxyz_yyyz_0[j] = fr * tg_xyz_yyyz_0[j] + tg_xyz_xyyyz_0[j];

                    tg_xxyz_yyzz_0[j] = fr * tg_xyz_yyzz_0[j] + tg_xyz_xyyzz_0[j];

                    tg_xxyz_yzzz_0[j] = fr * tg_xyz_yzzz_0[j] + tg_xyz_xyzzz_0[j];

                    tg_xxyz_zzzz_0[j] = fr * tg_xyz_zzzz_0[j] + tg_xyz_xzzzz_0[j];
                }
            }
        }
    }

    void
    compElectronRepulsionForSXGG_75_150(      CMemBlock2D<double>& ketBuffer,
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

        // set up bra side loop over intergals

        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();

        for (int32_t iang = 0; iang <= bang; iang++)
        {
            // set up index of integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {4, 4, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 4, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 5, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

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

                auto tg_yyy_xzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 140); 

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

                // set up pointers to integrals

                auto tg_xxzz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 75); 

                auto tg_xxzz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 76); 

                auto tg_xxzz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 77); 

                auto tg_xxzz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 78); 

                auto tg_xxzz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 79); 

                auto tg_xxzz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 80); 

                auto tg_xxzz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 81); 

                auto tg_xxzz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 82); 

                auto tg_xxzz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 83); 

                auto tg_xxzz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 84); 

                auto tg_xxzz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 85); 

                auto tg_xxzz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 86); 

                auto tg_xxzz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 87); 

                auto tg_xxzz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 88); 

                auto tg_xxzz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 89); 

                auto tg_xyyy_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 90); 

                auto tg_xyyy_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 91); 

                auto tg_xyyy_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 92); 

                auto tg_xyyy_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 93); 

                auto tg_xyyy_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 94); 

                auto tg_xyyy_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 95); 

                auto tg_xyyy_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 96); 

                auto tg_xyyy_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 97); 

                auto tg_xyyy_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 98); 

                auto tg_xyyy_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 99); 

                auto tg_xyyy_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 100); 

                auto tg_xyyy_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 101); 

                auto tg_xyyy_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 102); 

                auto tg_xyyy_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 103); 

                auto tg_xyyy_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 104); 

                auto tg_xyyz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 105); 

                auto tg_xyyz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 106); 

                auto tg_xyyz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 107); 

                auto tg_xyyz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 108); 

                auto tg_xyyz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 109); 

                auto tg_xyyz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 110); 

                auto tg_xyyz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 111); 

                auto tg_xyyz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 112); 

                auto tg_xyyz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 113); 

                auto tg_xyyz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 114); 

                auto tg_xyyz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 115); 

                auto tg_xyyz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 116); 

                auto tg_xyyz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 117); 

                auto tg_xyyz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 118); 

                auto tg_xyyz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 119); 

                auto tg_xyzz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 120); 

                auto tg_xyzz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 121); 

                auto tg_xyzz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 122); 

                auto tg_xyzz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 123); 

                auto tg_xyzz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 124); 

                auto tg_xyzz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 125); 

                auto tg_xyzz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 126); 

                auto tg_xyzz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 127); 

                auto tg_xyzz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 128); 

                auto tg_xyzz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 129); 

                auto tg_xyzz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 130); 

                auto tg_xyzz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 131); 

                auto tg_xyzz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 132); 

                auto tg_xyzz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 133); 

                auto tg_xyzz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 134); 

                auto tg_xzzz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 135); 

                auto tg_xzzz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 136); 

                auto tg_xzzz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 137); 

                auto tg_xzzz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 138); 

                auto tg_xzzz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 139); 

                auto tg_xzzz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 140); 

                auto tg_xzzz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 141); 

                auto tg_xzzz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 142); 

                auto tg_xzzz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 143); 

                auto tg_xzzz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 144); 

                auto tg_xzzz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 145); 

                auto tg_xzzz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 146); 

                auto tg_xzzz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 147); 

                auto tg_xzzz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 148); 

                auto tg_xzzz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 149); 

                // Batch of Integrals (75,150)

                #pragma omp simd aligned(cd_x, tg_xxzz_xxxx_0, tg_xxzz_xxxy_0, tg_xxzz_xxxz_0, tg_xxzz_xxyy_0, \
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
                    double fr = -cd_x[j]; 

                    tg_xxzz_xxxx_0[j] = fr * tg_xzz_xxxx_0[j] + tg_xzz_xxxxx_0[j];

                    tg_xxzz_xxxy_0[j] = fr * tg_xzz_xxxy_0[j] + tg_xzz_xxxxy_0[j];

                    tg_xxzz_xxxz_0[j] = fr * tg_xzz_xxxz_0[j] + tg_xzz_xxxxz_0[j];

                    tg_xxzz_xxyy_0[j] = fr * tg_xzz_xxyy_0[j] + tg_xzz_xxxyy_0[j];

                    tg_xxzz_xxyz_0[j] = fr * tg_xzz_xxyz_0[j] + tg_xzz_xxxyz_0[j];

                    tg_xxzz_xxzz_0[j] = fr * tg_xzz_xxzz_0[j] + tg_xzz_xxxzz_0[j];

                    tg_xxzz_xyyy_0[j] = fr * tg_xzz_xyyy_0[j] + tg_xzz_xxyyy_0[j];

                    tg_xxzz_xyyz_0[j] = fr * tg_xzz_xyyz_0[j] + tg_xzz_xxyyz_0[j];

                    tg_xxzz_xyzz_0[j] = fr * tg_xzz_xyzz_0[j] + tg_xzz_xxyzz_0[j];

                    tg_xxzz_xzzz_0[j] = fr * tg_xzz_xzzz_0[j] + tg_xzz_xxzzz_0[j];

                    tg_xxzz_yyyy_0[j] = fr * tg_xzz_yyyy_0[j] + tg_xzz_xyyyy_0[j];

                    tg_xxzz_yyyz_0[j] = fr * tg_xzz_yyyz_0[j] + tg_xzz_xyyyz_0[j];

                    tg_xxzz_yyzz_0[j] = fr * tg_xzz_yyzz_0[j] + tg_xzz_xyyzz_0[j];

                    tg_xxzz_yzzz_0[j] = fr * tg_xzz_yzzz_0[j] + tg_xzz_xyzzz_0[j];

                    tg_xxzz_zzzz_0[j] = fr * tg_xzz_zzzz_0[j] + tg_xzz_xzzzz_0[j];

                    tg_xyyy_xxxx_0[j] = fr * tg_yyy_xxxx_0[j] + tg_yyy_xxxxx_0[j];

                    tg_xyyy_xxxy_0[j] = fr * tg_yyy_xxxy_0[j] + tg_yyy_xxxxy_0[j];

                    tg_xyyy_xxxz_0[j] = fr * tg_yyy_xxxz_0[j] + tg_yyy_xxxxz_0[j];

                    tg_xyyy_xxyy_0[j] = fr * tg_yyy_xxyy_0[j] + tg_yyy_xxxyy_0[j];

                    tg_xyyy_xxyz_0[j] = fr * tg_yyy_xxyz_0[j] + tg_yyy_xxxyz_0[j];

                    tg_xyyy_xxzz_0[j] = fr * tg_yyy_xxzz_0[j] + tg_yyy_xxxzz_0[j];

                    tg_xyyy_xyyy_0[j] = fr * tg_yyy_xyyy_0[j] + tg_yyy_xxyyy_0[j];

                    tg_xyyy_xyyz_0[j] = fr * tg_yyy_xyyz_0[j] + tg_yyy_xxyyz_0[j];

                    tg_xyyy_xyzz_0[j] = fr * tg_yyy_xyzz_0[j] + tg_yyy_xxyzz_0[j];

                    tg_xyyy_xzzz_0[j] = fr * tg_yyy_xzzz_0[j] + tg_yyy_xxzzz_0[j];

                    tg_xyyy_yyyy_0[j] = fr * tg_yyy_yyyy_0[j] + tg_yyy_xyyyy_0[j];

                    tg_xyyy_yyyz_0[j] = fr * tg_yyy_yyyz_0[j] + tg_yyy_xyyyz_0[j];

                    tg_xyyy_yyzz_0[j] = fr * tg_yyy_yyzz_0[j] + tg_yyy_xyyzz_0[j];

                    tg_xyyy_yzzz_0[j] = fr * tg_yyy_yzzz_0[j] + tg_yyy_xyzzz_0[j];

                    tg_xyyy_zzzz_0[j] = fr * tg_yyy_zzzz_0[j] + tg_yyy_xzzzz_0[j];

                    tg_xyyz_xxxx_0[j] = fr * tg_yyz_xxxx_0[j] + tg_yyz_xxxxx_0[j];

                    tg_xyyz_xxxy_0[j] = fr * tg_yyz_xxxy_0[j] + tg_yyz_xxxxy_0[j];

                    tg_xyyz_xxxz_0[j] = fr * tg_yyz_xxxz_0[j] + tg_yyz_xxxxz_0[j];

                    tg_xyyz_xxyy_0[j] = fr * tg_yyz_xxyy_0[j] + tg_yyz_xxxyy_0[j];

                    tg_xyyz_xxyz_0[j] = fr * tg_yyz_xxyz_0[j] + tg_yyz_xxxyz_0[j];

                    tg_xyyz_xxzz_0[j] = fr * tg_yyz_xxzz_0[j] + tg_yyz_xxxzz_0[j];

                    tg_xyyz_xyyy_0[j] = fr * tg_yyz_xyyy_0[j] + tg_yyz_xxyyy_0[j];

                    tg_xyyz_xyyz_0[j] = fr * tg_yyz_xyyz_0[j] + tg_yyz_xxyyz_0[j];

                    tg_xyyz_xyzz_0[j] = fr * tg_yyz_xyzz_0[j] + tg_yyz_xxyzz_0[j];

                    tg_xyyz_xzzz_0[j] = fr * tg_yyz_xzzz_0[j] + tg_yyz_xxzzz_0[j];

                    tg_xyyz_yyyy_0[j] = fr * tg_yyz_yyyy_0[j] + tg_yyz_xyyyy_0[j];

                    tg_xyyz_yyyz_0[j] = fr * tg_yyz_yyyz_0[j] + tg_yyz_xyyyz_0[j];

                    tg_xyyz_yyzz_0[j] = fr * tg_yyz_yyzz_0[j] + tg_yyz_xyyzz_0[j];

                    tg_xyyz_yzzz_0[j] = fr * tg_yyz_yzzz_0[j] + tg_yyz_xyzzz_0[j];

                    tg_xyyz_zzzz_0[j] = fr * tg_yyz_zzzz_0[j] + tg_yyz_xzzzz_0[j];

                    tg_xyzz_xxxx_0[j] = fr * tg_yzz_xxxx_0[j] + tg_yzz_xxxxx_0[j];

                    tg_xyzz_xxxy_0[j] = fr * tg_yzz_xxxy_0[j] + tg_yzz_xxxxy_0[j];

                    tg_xyzz_xxxz_0[j] = fr * tg_yzz_xxxz_0[j] + tg_yzz_xxxxz_0[j];

                    tg_xyzz_xxyy_0[j] = fr * tg_yzz_xxyy_0[j] + tg_yzz_xxxyy_0[j];

                    tg_xyzz_xxyz_0[j] = fr * tg_yzz_xxyz_0[j] + tg_yzz_xxxyz_0[j];

                    tg_xyzz_xxzz_0[j] = fr * tg_yzz_xxzz_0[j] + tg_yzz_xxxzz_0[j];

                    tg_xyzz_xyyy_0[j] = fr * tg_yzz_xyyy_0[j] + tg_yzz_xxyyy_0[j];

                    tg_xyzz_xyyz_0[j] = fr * tg_yzz_xyyz_0[j] + tg_yzz_xxyyz_0[j];

                    tg_xyzz_xyzz_0[j] = fr * tg_yzz_xyzz_0[j] + tg_yzz_xxyzz_0[j];

                    tg_xyzz_xzzz_0[j] = fr * tg_yzz_xzzz_0[j] + tg_yzz_xxzzz_0[j];

                    tg_xyzz_yyyy_0[j] = fr * tg_yzz_yyyy_0[j] + tg_yzz_xyyyy_0[j];

                    tg_xyzz_yyyz_0[j] = fr * tg_yzz_yyyz_0[j] + tg_yzz_xyyyz_0[j];

                    tg_xyzz_yyzz_0[j] = fr * tg_yzz_yyzz_0[j] + tg_yzz_xyyzz_0[j];

                    tg_xyzz_yzzz_0[j] = fr * tg_yzz_yzzz_0[j] + tg_yzz_xyzzz_0[j];

                    tg_xyzz_zzzz_0[j] = fr * tg_yzz_zzzz_0[j] + tg_yzz_xzzzz_0[j];

                    tg_xzzz_xxxx_0[j] = fr * tg_zzz_xxxx_0[j] + tg_zzz_xxxxx_0[j];

                    tg_xzzz_xxxy_0[j] = fr * tg_zzz_xxxy_0[j] + tg_zzz_xxxxy_0[j];

                    tg_xzzz_xxxz_0[j] = fr * tg_zzz_xxxz_0[j] + tg_zzz_xxxxz_0[j];

                    tg_xzzz_xxyy_0[j] = fr * tg_zzz_xxyy_0[j] + tg_zzz_xxxyy_0[j];

                    tg_xzzz_xxyz_0[j] = fr * tg_zzz_xxyz_0[j] + tg_zzz_xxxyz_0[j];

                    tg_xzzz_xxzz_0[j] = fr * tg_zzz_xxzz_0[j] + tg_zzz_xxxzz_0[j];

                    tg_xzzz_xyyy_0[j] = fr * tg_zzz_xyyy_0[j] + tg_zzz_xxyyy_0[j];

                    tg_xzzz_xyyz_0[j] = fr * tg_zzz_xyyz_0[j] + tg_zzz_xxyyz_0[j];

                    tg_xzzz_xyzz_0[j] = fr * tg_zzz_xyzz_0[j] + tg_zzz_xxyzz_0[j];

                    tg_xzzz_xzzz_0[j] = fr * tg_zzz_xzzz_0[j] + tg_zzz_xxzzz_0[j];

                    tg_xzzz_yyyy_0[j] = fr * tg_zzz_yyyy_0[j] + tg_zzz_xyyyy_0[j];

                    tg_xzzz_yyyz_0[j] = fr * tg_zzz_yyyz_0[j] + tg_zzz_xyyyz_0[j];

                    tg_xzzz_yyzz_0[j] = fr * tg_zzz_yyzz_0[j] + tg_zzz_xyyzz_0[j];

                    tg_xzzz_yzzz_0[j] = fr * tg_zzz_yzzz_0[j] + tg_zzz_xyzzz_0[j];

                    tg_xzzz_zzzz_0[j] = fr * tg_zzz_zzzz_0[j] + tg_zzz_xzzzz_0[j];
                }
            }
        }
    }

    void
    compElectronRepulsionForSXGG_150_225(      CMemBlock2D<double>& ketBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& cdDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetContrPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (150,225)

        // set up pointers to distances R(CD) = C - D

        auto cd_y = cdDistances.data(1);

        auto cd_z = cdDistances.data(2);

        // set up bra side loop over intergals

        auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();

        for (int32_t iang = 0; iang <= bang; iang++)
        {
            // set up index of integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {4, 4, -1, -1}, 
                                                             1, 2, 0));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 4, -1, -1}, 
                                                             1, 2, 0));

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {iang, -1, -1, -1}, {3, 5, -1, -1}, 
                                                             1, 2, 0));

            auto bcomp = angmom::to_CartesianComponents(iang);

            for (int32_t i = 0; i < bcomp; i++)
            {
                // set up pointers to auxilary integrals

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

                auto tg_yyy_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 127); 

                auto tg_yyy_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 129); 

                auto tg_yyy_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 130); 

                auto tg_yyy_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 132); 

                auto tg_yyy_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 133); 

                auto tg_yyy_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 134); 

                auto tg_yyy_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 136); 

                auto tg_yyy_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 137); 

                auto tg_yyy_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 138); 

                auto tg_yyy_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 139); 

                auto tg_yyy_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 141); 

                auto tg_yyy_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 142); 

                auto tg_yyy_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 143); 

                auto tg_yyy_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 144); 

                auto tg_yyy_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 145); 

                auto tg_yyz_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 148); 

                auto tg_yyz_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 150); 

                auto tg_yyz_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 151); 

                auto tg_yyz_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 153); 

                auto tg_yyz_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 154); 

                auto tg_yyz_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 155); 

                auto tg_yyz_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 157); 

                auto tg_yyz_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 158); 

                auto tg_yyz_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 159); 

                auto tg_yyz_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 160); 

                auto tg_yyz_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 162); 

                auto tg_yyz_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 163); 

                auto tg_yyz_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 164); 

                auto tg_yyz_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 165); 

                auto tg_yyz_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 166); 

                auto tg_yzz_xxxxy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 169); 

                auto tg_yzz_xxxyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 171); 

                auto tg_yzz_xxxyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 172); 

                auto tg_yzz_xxyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 174); 

                auto tg_yzz_xxyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 175); 

                auto tg_yzz_xxyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 176); 

                auto tg_yzz_xyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 178); 

                auto tg_yzz_xyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 179); 

                auto tg_yzz_xyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 180); 

                auto tg_yzz_xyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 181); 

                auto tg_yzz_yyyyy_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 183); 

                auto tg_yzz_yyyyz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 184); 

                auto tg_yzz_yyyzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 185); 

                auto tg_yzz_yyzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 186); 

                auto tg_yzz_yzzzz_0 = ketBuffer.data(pidx_g_3_5_m0 + 210 * i + 187); 

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

                // set up pointers to integrals

                auto tg_yyyy_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 150); 

                auto tg_yyyy_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 151); 

                auto tg_yyyy_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 152); 

                auto tg_yyyy_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 153); 

                auto tg_yyyy_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 154); 

                auto tg_yyyy_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 155); 

                auto tg_yyyy_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 156); 

                auto tg_yyyy_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 157); 

                auto tg_yyyy_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 158); 

                auto tg_yyyy_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 159); 

                auto tg_yyyy_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 160); 

                auto tg_yyyy_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 161); 

                auto tg_yyyy_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 162); 

                auto tg_yyyy_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 163); 

                auto tg_yyyy_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 164); 

                auto tg_yyyz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 165); 

                auto tg_yyyz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 166); 

                auto tg_yyyz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 167); 

                auto tg_yyyz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 168); 

                auto tg_yyyz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 169); 

                auto tg_yyyz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 170); 

                auto tg_yyyz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 171); 

                auto tg_yyyz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 172); 

                auto tg_yyyz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 173); 

                auto tg_yyyz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 174); 

                auto tg_yyyz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 175); 

                auto tg_yyyz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 176); 

                auto tg_yyyz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 177); 

                auto tg_yyyz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 178); 

                auto tg_yyyz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 179); 

                auto tg_yyzz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 180); 

                auto tg_yyzz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 181); 

                auto tg_yyzz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 182); 

                auto tg_yyzz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 183); 

                auto tg_yyzz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 184); 

                auto tg_yyzz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 185); 

                auto tg_yyzz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 186); 

                auto tg_yyzz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 187); 

                auto tg_yyzz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 188); 

                auto tg_yyzz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 189); 

                auto tg_yyzz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 190); 

                auto tg_yyzz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 191); 

                auto tg_yyzz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 192); 

                auto tg_yyzz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 193); 

                auto tg_yyzz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 194); 

                auto tg_yzzz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 195); 

                auto tg_yzzz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 196); 

                auto tg_yzzz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 197); 

                auto tg_yzzz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 198); 

                auto tg_yzzz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 199); 

                auto tg_yzzz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 200); 

                auto tg_yzzz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 201); 

                auto tg_yzzz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 202); 

                auto tg_yzzz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 203); 

                auto tg_yzzz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 204); 

                auto tg_yzzz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 205); 

                auto tg_yzzz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 206); 

                auto tg_yzzz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 207); 

                auto tg_yzzz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 208); 

                auto tg_yzzz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 209); 

                auto tg_zzzz_xxxx_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 210); 

                auto tg_zzzz_xxxy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 211); 

                auto tg_zzzz_xxxz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 212); 

                auto tg_zzzz_xxyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 213); 

                auto tg_zzzz_xxyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 214); 

                auto tg_zzzz_xxzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 215); 

                auto tg_zzzz_xyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 216); 

                auto tg_zzzz_xyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 217); 

                auto tg_zzzz_xyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 218); 

                auto tg_zzzz_xzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 219); 

                auto tg_zzzz_yyyy_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 220); 

                auto tg_zzzz_yyyz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 221); 

                auto tg_zzzz_yyzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 222); 

                auto tg_zzzz_yzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 223); 

                auto tg_zzzz_zzzz_0 = ketBuffer.data(pidx_g_4_4_m0 + 225 * i + 224); 

                // Batch of Integrals (150,225)

                #pragma omp simd aligned(cd_y, cd_z, tg_yyy_xxxx_0, tg_yyy_xxxxy_0, tg_yyy_xxxy_0, \
                                         tg_yyy_xxxyy_0, tg_yyy_xxxyz_0, tg_yyy_xxxz_0, tg_yyy_xxyy_0, tg_yyy_xxyyy_0, \
                                         tg_yyy_xxyyz_0, tg_yyy_xxyz_0, tg_yyy_xxyzz_0, tg_yyy_xxzz_0, tg_yyy_xyyy_0, \
                                         tg_yyy_xyyyy_0, tg_yyy_xyyyz_0, tg_yyy_xyyz_0, tg_yyy_xyyzz_0, tg_yyy_xyzz_0, \
                                         tg_yyy_xyzzz_0, tg_yyy_xzzz_0, tg_yyy_yyyy_0, tg_yyy_yyyyy_0, tg_yyy_yyyyz_0, \
                                         tg_yyy_yyyz_0, tg_yyy_yyyzz_0, tg_yyy_yyzz_0, tg_yyy_yyzzz_0, tg_yyy_yzzz_0, \
                                         tg_yyy_yzzzz_0, tg_yyy_zzzz_0, tg_yyyy_xxxx_0, tg_yyyy_xxxy_0, tg_yyyy_xxxz_0, \
                                         tg_yyyy_xxyy_0, tg_yyyy_xxyz_0, tg_yyyy_xxzz_0, tg_yyyy_xyyy_0, tg_yyyy_xyyz_0, \
                                         tg_yyyy_xyzz_0, tg_yyyy_xzzz_0, tg_yyyy_yyyy_0, tg_yyyy_yyyz_0, tg_yyyy_yyzz_0, \
                                         tg_yyyy_yzzz_0, tg_yyyy_zzzz_0, tg_yyyz_xxxx_0, tg_yyyz_xxxy_0, tg_yyyz_xxxz_0, \
                                         tg_yyyz_xxyy_0, tg_yyyz_xxyz_0, tg_yyyz_xxzz_0, tg_yyyz_xyyy_0, tg_yyyz_xyyz_0, \
                                         tg_yyyz_xyzz_0, tg_yyyz_xzzz_0, tg_yyyz_yyyy_0, tg_yyyz_yyyz_0, tg_yyyz_yyzz_0, \
                                         tg_yyyz_yzzz_0, tg_yyyz_zzzz_0, tg_yyz_xxxx_0, tg_yyz_xxxxy_0, tg_yyz_xxxy_0, \
                                         tg_yyz_xxxyy_0, tg_yyz_xxxyz_0, tg_yyz_xxxz_0, tg_yyz_xxyy_0, tg_yyz_xxyyy_0, \
                                         tg_yyz_xxyyz_0, tg_yyz_xxyz_0, tg_yyz_xxyzz_0, tg_yyz_xxzz_0, tg_yyz_xyyy_0, \
                                         tg_yyz_xyyyy_0, tg_yyz_xyyyz_0, tg_yyz_xyyz_0, tg_yyz_xyyzz_0, tg_yyz_xyzz_0, \
                                         tg_yyz_xyzzz_0, tg_yyz_xzzz_0, tg_yyz_yyyy_0, tg_yyz_yyyyy_0, tg_yyz_yyyyz_0, \
                                         tg_yyz_yyyz_0, tg_yyz_yyyzz_0, tg_yyz_yyzz_0, tg_yyz_yyzzz_0, tg_yyz_yzzz_0, \
                                         tg_yyz_yzzzz_0, tg_yyz_zzzz_0, tg_yyzz_xxxx_0, tg_yyzz_xxxy_0, tg_yyzz_xxxz_0, \
                                         tg_yyzz_xxyy_0, tg_yyzz_xxyz_0, tg_yyzz_xxzz_0, tg_yyzz_xyyy_0, tg_yyzz_xyyz_0, \
                                         tg_yyzz_xyzz_0, tg_yyzz_xzzz_0, tg_yyzz_yyyy_0, tg_yyzz_yyyz_0, tg_yyzz_yyzz_0, \
                                         tg_yyzz_yzzz_0, tg_yyzz_zzzz_0, tg_yzz_xxxx_0, tg_yzz_xxxxy_0, tg_yzz_xxxy_0, \
                                         tg_yzz_xxxyy_0, tg_yzz_xxxyz_0, tg_yzz_xxxz_0, tg_yzz_xxyy_0, tg_yzz_xxyyy_0, \
                                         tg_yzz_xxyyz_0, tg_yzz_xxyz_0, tg_yzz_xxyzz_0, tg_yzz_xxzz_0, tg_yzz_xyyy_0, \
                                         tg_yzz_xyyyy_0, tg_yzz_xyyyz_0, tg_yzz_xyyz_0, tg_yzz_xyyzz_0, tg_yzz_xyzz_0, \
                                         tg_yzz_xyzzz_0, tg_yzz_xzzz_0, tg_yzz_yyyy_0, tg_yzz_yyyyy_0, tg_yzz_yyyyz_0, \
                                         tg_yzz_yyyz_0, tg_yzz_yyyzz_0, tg_yzz_yyzz_0, tg_yzz_yyzzz_0, tg_yzz_yzzz_0, \
                                         tg_yzz_yzzzz_0, tg_yzz_zzzz_0, tg_yzzz_xxxx_0, tg_yzzz_xxxy_0, tg_yzzz_xxxz_0, \
                                         tg_yzzz_xxyy_0, tg_yzzz_xxyz_0, tg_yzzz_xxzz_0, tg_yzzz_xyyy_0, tg_yzzz_xyyz_0, \
                                         tg_yzzz_xyzz_0, tg_yzzz_xzzz_0, tg_yzzz_yyyy_0, tg_yzzz_yyyz_0, tg_yzzz_yyzz_0, \
                                         tg_yzzz_yzzz_0, tg_yzzz_zzzz_0, tg_zzz_xxxx_0, tg_zzz_xxxxy_0, tg_zzz_xxxxz_0, \
                                         tg_zzz_xxxy_0, tg_zzz_xxxyy_0, tg_zzz_xxxyz_0, tg_zzz_xxxz_0, tg_zzz_xxxzz_0, \
                                         tg_zzz_xxyy_0, tg_zzz_xxyyy_0, tg_zzz_xxyyz_0, tg_zzz_xxyz_0, tg_zzz_xxyzz_0, \
                                         tg_zzz_xxzz_0, tg_zzz_xxzzz_0, tg_zzz_xyyy_0, tg_zzz_xyyyy_0, tg_zzz_xyyyz_0, \
                                         tg_zzz_xyyz_0, tg_zzz_xyyzz_0, tg_zzz_xyzz_0, tg_zzz_xyzzz_0, tg_zzz_xzzz_0, \
                                         tg_zzz_xzzzz_0, tg_zzz_yyyy_0, tg_zzz_yyyyy_0, tg_zzz_yyyyz_0, tg_zzz_yyyz_0, \
                                         tg_zzz_yyyzz_0, tg_zzz_yyzz_0, tg_zzz_yyzzz_0, tg_zzz_yzzz_0, tg_zzz_yzzzz_0, \
                                         tg_zzz_zzzz_0, tg_zzz_zzzzz_0, tg_zzzz_xxxx_0, tg_zzzz_xxxy_0, tg_zzzz_xxxz_0, \
                                         tg_zzzz_xxyy_0, tg_zzzz_xxyz_0, tg_zzzz_xxzz_0, tg_zzzz_xyyy_0, tg_zzzz_xyyz_0, \
                                         tg_zzzz_xyzz_0, tg_zzzz_xzzz_0, tg_zzzz_yyyy_0, tg_zzzz_yyyz_0, tg_zzzz_yyzz_0, \
                                         tg_zzzz_yzzz_0, tg_zzzz_zzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nKetContrPairs; j++)
                {
                    double fr = -cd_y[j]; 

                    tg_yyyy_xxxx_0[j] = fr * tg_yyy_xxxx_0[j] + tg_yyy_xxxxy_0[j];

                    tg_yyyy_xxxy_0[j] = fr * tg_yyy_xxxy_0[j] + tg_yyy_xxxyy_0[j];

                    tg_yyyy_xxxz_0[j] = fr * tg_yyy_xxxz_0[j] + tg_yyy_xxxyz_0[j];

                    tg_yyyy_xxyy_0[j] = fr * tg_yyy_xxyy_0[j] + tg_yyy_xxyyy_0[j];

                    tg_yyyy_xxyz_0[j] = fr * tg_yyy_xxyz_0[j] + tg_yyy_xxyyz_0[j];

                    tg_yyyy_xxzz_0[j] = fr * tg_yyy_xxzz_0[j] + tg_yyy_xxyzz_0[j];

                    tg_yyyy_xyyy_0[j] = fr * tg_yyy_xyyy_0[j] + tg_yyy_xyyyy_0[j];

                    tg_yyyy_xyyz_0[j] = fr * tg_yyy_xyyz_0[j] + tg_yyy_xyyyz_0[j];

                    tg_yyyy_xyzz_0[j] = fr * tg_yyy_xyzz_0[j] + tg_yyy_xyyzz_0[j];

                    tg_yyyy_xzzz_0[j] = fr * tg_yyy_xzzz_0[j] + tg_yyy_xyzzz_0[j];

                    tg_yyyy_yyyy_0[j] = fr * tg_yyy_yyyy_0[j] + tg_yyy_yyyyy_0[j];

                    tg_yyyy_yyyz_0[j] = fr * tg_yyy_yyyz_0[j] + tg_yyy_yyyyz_0[j];

                    tg_yyyy_yyzz_0[j] = fr * tg_yyy_yyzz_0[j] + tg_yyy_yyyzz_0[j];

                    tg_yyyy_yzzz_0[j] = fr * tg_yyy_yzzz_0[j] + tg_yyy_yyzzz_0[j];

                    tg_yyyy_zzzz_0[j] = fr * tg_yyy_zzzz_0[j] + tg_yyy_yzzzz_0[j];

                    tg_yyyz_xxxx_0[j] = fr * tg_yyz_xxxx_0[j] + tg_yyz_xxxxy_0[j];

                    tg_yyyz_xxxy_0[j] = fr * tg_yyz_xxxy_0[j] + tg_yyz_xxxyy_0[j];

                    tg_yyyz_xxxz_0[j] = fr * tg_yyz_xxxz_0[j] + tg_yyz_xxxyz_0[j];

                    tg_yyyz_xxyy_0[j] = fr * tg_yyz_xxyy_0[j] + tg_yyz_xxyyy_0[j];

                    tg_yyyz_xxyz_0[j] = fr * tg_yyz_xxyz_0[j] + tg_yyz_xxyyz_0[j];

                    tg_yyyz_xxzz_0[j] = fr * tg_yyz_xxzz_0[j] + tg_yyz_xxyzz_0[j];

                    tg_yyyz_xyyy_0[j] = fr * tg_yyz_xyyy_0[j] + tg_yyz_xyyyy_0[j];

                    tg_yyyz_xyyz_0[j] = fr * tg_yyz_xyyz_0[j] + tg_yyz_xyyyz_0[j];

                    tg_yyyz_xyzz_0[j] = fr * tg_yyz_xyzz_0[j] + tg_yyz_xyyzz_0[j];

                    tg_yyyz_xzzz_0[j] = fr * tg_yyz_xzzz_0[j] + tg_yyz_xyzzz_0[j];

                    tg_yyyz_yyyy_0[j] = fr * tg_yyz_yyyy_0[j] + tg_yyz_yyyyy_0[j];

                    tg_yyyz_yyyz_0[j] = fr * tg_yyz_yyyz_0[j] + tg_yyz_yyyyz_0[j];

                    tg_yyyz_yyzz_0[j] = fr * tg_yyz_yyzz_0[j] + tg_yyz_yyyzz_0[j];

                    tg_yyyz_yzzz_0[j] = fr * tg_yyz_yzzz_0[j] + tg_yyz_yyzzz_0[j];

                    tg_yyyz_zzzz_0[j] = fr * tg_yyz_zzzz_0[j] + tg_yyz_yzzzz_0[j];

                    tg_yyzz_xxxx_0[j] = fr * tg_yzz_xxxx_0[j] + tg_yzz_xxxxy_0[j];

                    tg_yyzz_xxxy_0[j] = fr * tg_yzz_xxxy_0[j] + tg_yzz_xxxyy_0[j];

                    tg_yyzz_xxxz_0[j] = fr * tg_yzz_xxxz_0[j] + tg_yzz_xxxyz_0[j];

                    tg_yyzz_xxyy_0[j] = fr * tg_yzz_xxyy_0[j] + tg_yzz_xxyyy_0[j];

                    tg_yyzz_xxyz_0[j] = fr * tg_yzz_xxyz_0[j] + tg_yzz_xxyyz_0[j];

                    tg_yyzz_xxzz_0[j] = fr * tg_yzz_xxzz_0[j] + tg_yzz_xxyzz_0[j];

                    tg_yyzz_xyyy_0[j] = fr * tg_yzz_xyyy_0[j] + tg_yzz_xyyyy_0[j];

                    tg_yyzz_xyyz_0[j] = fr * tg_yzz_xyyz_0[j] + tg_yzz_xyyyz_0[j];

                    tg_yyzz_xyzz_0[j] = fr * tg_yzz_xyzz_0[j] + tg_yzz_xyyzz_0[j];

                    tg_yyzz_xzzz_0[j] = fr * tg_yzz_xzzz_0[j] + tg_yzz_xyzzz_0[j];

                    tg_yyzz_yyyy_0[j] = fr * tg_yzz_yyyy_0[j] + tg_yzz_yyyyy_0[j];

                    tg_yyzz_yyyz_0[j] = fr * tg_yzz_yyyz_0[j] + tg_yzz_yyyyz_0[j];

                    tg_yyzz_yyzz_0[j] = fr * tg_yzz_yyzz_0[j] + tg_yzz_yyyzz_0[j];

                    tg_yyzz_yzzz_0[j] = fr * tg_yzz_yzzz_0[j] + tg_yzz_yyzzz_0[j];

                    tg_yyzz_zzzz_0[j] = fr * tg_yzz_zzzz_0[j] + tg_yzz_yzzzz_0[j];

                    tg_yzzz_xxxx_0[j] = fr * tg_zzz_xxxx_0[j] + tg_zzz_xxxxy_0[j];

                    tg_yzzz_xxxy_0[j] = fr * tg_zzz_xxxy_0[j] + tg_zzz_xxxyy_0[j];

                    tg_yzzz_xxxz_0[j] = fr * tg_zzz_xxxz_0[j] + tg_zzz_xxxyz_0[j];

                    tg_yzzz_xxyy_0[j] = fr * tg_zzz_xxyy_0[j] + tg_zzz_xxyyy_0[j];

                    tg_yzzz_xxyz_0[j] = fr * tg_zzz_xxyz_0[j] + tg_zzz_xxyyz_0[j];

                    tg_yzzz_xxzz_0[j] = fr * tg_zzz_xxzz_0[j] + tg_zzz_xxyzz_0[j];

                    tg_yzzz_xyyy_0[j] = fr * tg_zzz_xyyy_0[j] + tg_zzz_xyyyy_0[j];

                    tg_yzzz_xyyz_0[j] = fr * tg_zzz_xyyz_0[j] + tg_zzz_xyyyz_0[j];

                    tg_yzzz_xyzz_0[j] = fr * tg_zzz_xyzz_0[j] + tg_zzz_xyyzz_0[j];

                    tg_yzzz_xzzz_0[j] = fr * tg_zzz_xzzz_0[j] + tg_zzz_xyzzz_0[j];

                    tg_yzzz_yyyy_0[j] = fr * tg_zzz_yyyy_0[j] + tg_zzz_yyyyy_0[j];

                    tg_yzzz_yyyz_0[j] = fr * tg_zzz_yyyz_0[j] + tg_zzz_yyyyz_0[j];

                    tg_yzzz_yyzz_0[j] = fr * tg_zzz_yyzz_0[j] + tg_zzz_yyyzz_0[j];

                    tg_yzzz_yzzz_0[j] = fr * tg_zzz_yzzz_0[j] + tg_zzz_yyzzz_0[j];

                    tg_yzzz_zzzz_0[j] = fr * tg_zzz_zzzz_0[j] + tg_zzz_yzzzz_0[j];

                    fr = -cd_z[j]; 

                    tg_zzzz_xxxx_0[j] = fr * tg_zzz_xxxx_0[j] + tg_zzz_xxxxz_0[j];

                    tg_zzzz_xxxy_0[j] = fr * tg_zzz_xxxy_0[j] + tg_zzz_xxxyz_0[j];

                    tg_zzzz_xxxz_0[j] = fr * tg_zzz_xxxz_0[j] + tg_zzz_xxxzz_0[j];

                    tg_zzzz_xxyy_0[j] = fr * tg_zzz_xxyy_0[j] + tg_zzz_xxyyz_0[j];

                    tg_zzzz_xxyz_0[j] = fr * tg_zzz_xxyz_0[j] + tg_zzz_xxyzz_0[j];

                    tg_zzzz_xxzz_0[j] = fr * tg_zzz_xxzz_0[j] + tg_zzz_xxzzz_0[j];

                    tg_zzzz_xyyy_0[j] = fr * tg_zzz_xyyy_0[j] + tg_zzz_xyyyz_0[j];

                    tg_zzzz_xyyz_0[j] = fr * tg_zzz_xyyz_0[j] + tg_zzz_xyyzz_0[j];

                    tg_zzzz_xyzz_0[j] = fr * tg_zzz_xyzz_0[j] + tg_zzz_xyzzz_0[j];

                    tg_zzzz_xzzz_0[j] = fr * tg_zzz_xzzz_0[j] + tg_zzz_xzzzz_0[j];

                    tg_zzzz_yyyy_0[j] = fr * tg_zzz_yyyy_0[j] + tg_zzz_yyyyz_0[j];

                    tg_zzzz_yyyz_0[j] = fr * tg_zzz_yyyz_0[j] + tg_zzz_yyyzz_0[j];

                    tg_zzzz_yyzz_0[j] = fr * tg_zzz_yyzz_0[j] + tg_zzz_yyzzz_0[j];

                    tg_zzzz_yzzz_0[j] = fr * tg_zzz_yzzz_0[j] + tg_zzz_yzzzz_0[j];

                    tg_zzzz_zzzz_0[j] = fr * tg_zzz_zzzz_0[j] + tg_zzz_zzzzz_0[j];
                }
            }
        }
    }


} // erikrrfunc namespace

