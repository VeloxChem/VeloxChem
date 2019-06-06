//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricDipoleRecFuncForGF.hpp"

namespace ediprecfunc { // ediprecfunc namespace

    void
    compElectricDipoleForGF(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForGF_0_50(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForGF_50_100(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        ediprecfunc::compElectricDipoleForGF_100_150(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForGF_150_200(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForGF_200_250(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForGF_250_300(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForGF_300_350(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForGF_350_400(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForGF_400_450(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 
    }

    void
    compElectricDipoleForGF_0_50(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (0,50)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_xxx_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx); 

            auto tdy_xxx_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx); 

            auto tdz_xxx_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx); 

            auto tdx_xxx_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 1); 

            auto tdy_xxx_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 1); 

            auto tdz_xxx_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 1); 

            auto tdx_xxx_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 2); 

            auto tdy_xxx_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 2); 

            auto tdz_xxx_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 2); 

            auto tdx_xxx_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 3); 

            auto tdy_xxx_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 3); 

            auto tdz_xxx_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 3); 

            auto tdx_xxx_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 4); 

            auto tdy_xxx_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 4); 

            auto tdz_xxx_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 4); 

            auto tdx_xxx_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 5); 

            auto tdy_xxx_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 5); 

            auto tdz_xxx_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 5); 

            auto tdx_xxx_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 6); 

            auto tdy_xxx_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 6); 

            auto tdz_xxx_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 6); 

            auto tdx_xxx_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 7); 

            auto tdy_xxx_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 7); 

            auto tdz_xxx_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 7); 

            auto tdx_xxx_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 8); 

            auto tdy_xxx_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 8); 

            auto tdz_xxx_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 8); 

            auto tdx_xxx_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 9); 

            auto tdy_xxx_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 9); 

            auto tdz_xxx_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 9); 

            auto tdx_xxy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 10); 

            auto tdy_xxy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 10); 

            auto tdz_xxy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 10); 

            auto tdx_xxy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 11); 

            auto tdy_xxy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 11); 

            auto tdz_xxy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 11); 

            auto tdx_xxy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 12); 

            auto tdy_xxy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 12); 

            auto tdz_xxy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 12); 

            auto tdx_xxy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 13); 

            auto tdy_xxy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 13); 

            auto tdz_xxy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 13); 

            auto tdx_xxy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 14); 

            auto tdy_xxy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 14); 

            auto tdz_xxy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 14); 

            auto tdx_xxy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 15); 

            auto tdy_xxy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 15); 

            auto tdz_xxy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 15); 

            auto tdx_xxy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 16); 

            auto tdy_xxy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 16); 

            auto tdx_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx); 

            auto tdy_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx); 

            auto tdz_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx); 

            auto tdx_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 1); 

            auto tdy_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 1); 

            auto tdz_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 1); 

            auto tdx_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 2); 

            auto tdy_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 2); 

            auto tdz_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 2); 

            auto tdx_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 3); 

            auto tdy_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 3); 

            auto tdz_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 3); 

            auto tdx_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 4); 

            auto tdy_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 4); 

            auto tdz_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 4); 

            auto tdx_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 5); 

            auto tdy_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 5); 

            auto tdz_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 5); 

            auto tdx_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 6); 

            auto tdy_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 6); 

            auto tdz_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 6); 

            auto tdx_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 7); 

            auto tdy_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 7); 

            auto tdz_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 7); 

            auto tdx_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 8); 

            auto tdy_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 8); 

            auto tdz_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 8); 

            auto tdx_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 9); 

            auto tdy_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 9); 

            auto tdz_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 9); 

            auto tdx_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 10); 

            auto tdy_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 10); 

            auto tdz_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 10); 

            auto tdx_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 11); 

            auto tdy_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 11); 

            auto tdz_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 11); 

            auto tdx_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 12); 

            auto tdy_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 12); 

            auto tdz_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 12); 

            auto tdx_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 13); 

            auto tdy_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 13); 

            auto tdz_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 13); 

            auto tdx_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 14); 

            auto tdy_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 14); 

            auto tdz_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 14); 

            auto tdx_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 15); 

            auto tdy_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 15); 

            auto tdz_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 15); 

            auto tdx_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 16); 

            auto tdy_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 16); 

            auto tdx_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx); 

            auto tdy_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx); 

            auto tdz_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx); 

            auto tdx_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 1); 

            auto tdy_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 1); 

            auto tdz_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 1); 

            auto tdx_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 2); 

            auto tdy_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 2); 

            auto tdz_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 2); 

            auto tdx_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 3); 

            auto tdy_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 3); 

            auto tdz_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 3); 

            auto tdx_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 4); 

            auto tdy_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 4); 

            auto tdz_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 4); 

            auto tdx_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 5); 

            auto tdy_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 5); 

            auto tdz_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 5); 

            auto tdx_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 6); 

            auto tdy_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 6); 

            auto tdz_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 6); 

            auto tdx_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 7); 

            auto tdy_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 7); 

            auto tdz_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 7); 

            auto tdx_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 8); 

            auto tdy_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 8); 

            auto tdz_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 8); 

            auto tdx_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 9); 

            auto tdy_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 9); 

            auto tdz_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 9); 

            auto tdx_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 10); 

            auto tdy_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 10); 

            auto tdz_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 10); 

            auto tdx_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 11); 

            auto tdy_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 11); 

            auto tdz_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 11); 

            auto ts_xxx_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx); 

            auto ts_xxx_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 1); 

            auto ts_xxx_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 2); 

            auto ts_xxx_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 3); 

            auto ts_xxx_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 4); 

            auto ts_xxx_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 5); 

            auto ts_xxx_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 6); 

            auto ts_xxx_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 7); 

            auto ts_xxx_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 8); 

            auto ts_xxx_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 9); 

            auto ts_xxy_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 10); 

            auto ts_xxy_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 11); 

            auto ts_xxy_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 12); 

            auto ts_xxy_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 13); 

            auto ts_xxy_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 14); 

            auto ts_xxy_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 15); 

            auto ts_xxy_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 16); 

            // set up pointers to integrals

            auto tdx_xxxx_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx); 

            auto tdy_xxxx_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx); 

            auto tdz_xxxx_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx); 

            auto tdx_xxxx_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 1); 

            auto tdy_xxxx_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 1); 

            auto tdz_xxxx_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 1); 

            auto tdx_xxxx_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 2); 

            auto tdy_xxxx_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 2); 

            auto tdz_xxxx_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 2); 

            auto tdx_xxxx_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 3); 

            auto tdy_xxxx_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 3); 

            auto tdz_xxxx_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 3); 

            auto tdx_xxxx_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 4); 

            auto tdy_xxxx_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 4); 

            auto tdz_xxxx_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 4); 

            auto tdx_xxxx_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 5); 

            auto tdy_xxxx_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 5); 

            auto tdz_xxxx_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 5); 

            auto tdx_xxxx_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 6); 

            auto tdy_xxxx_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 6); 

            auto tdz_xxxx_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 6); 

            auto tdx_xxxx_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 7); 

            auto tdy_xxxx_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 7); 

            auto tdz_xxxx_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 7); 

            auto tdx_xxxx_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 8); 

            auto tdy_xxxx_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 8); 

            auto tdz_xxxx_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 8); 

            auto tdx_xxxx_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 9); 

            auto tdy_xxxx_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 9); 

            auto tdz_xxxx_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 9); 

            auto tdx_xxxy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 10); 

            auto tdy_xxxy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 10); 

            auto tdz_xxxy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 10); 

            auto tdx_xxxy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 11); 

            auto tdy_xxxy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 11); 

            auto tdz_xxxy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 11); 

            auto tdx_xxxy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 12); 

            auto tdy_xxxy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 12); 

            auto tdz_xxxy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 12); 

            auto tdx_xxxy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 13); 

            auto tdy_xxxy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 13); 

            auto tdz_xxxy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 13); 

            auto tdx_xxxy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 14); 

            auto tdy_xxxy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 14); 

            auto tdz_xxxy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 14); 

            auto tdx_xxxy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 15); 

            auto tdy_xxxy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 15); 

            auto tdz_xxxy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 15); 

            auto tdx_xxxy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 16); 

            auto tdy_xxxy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 16); 

            // Batch of Integrals (0,50)

            #pragma omp simd aligned(fx, pa_x, tdx_xx_xxx_0, tdx_xx_xxy_0, tdx_xx_xxz_0, tdx_xx_xyy_0, \
                                     tdx_xx_xyz_0, tdx_xx_xzz_0, tdx_xx_yyy_0, tdx_xx_yyz_0, tdx_xx_yzz_0, tdx_xx_zzz_0, \
                                     tdx_xxx_xx_0, tdx_xxx_xxx_0, tdx_xxx_xxy_0, tdx_xxx_xxz_0, tdx_xxx_xy_0, \
                                     tdx_xxx_xyy_0, tdx_xxx_xyz_0, tdx_xxx_xz_0, tdx_xxx_xzz_0, tdx_xxx_yy_0, \
                                     tdx_xxx_yyy_0, tdx_xxx_yyz_0, tdx_xxx_yz_0, tdx_xxx_yzz_0, tdx_xxx_zz_0, \
                                     tdx_xxx_zzz_0, tdx_xxxx_xxx_0, tdx_xxxx_xxy_0, tdx_xxxx_xxz_0, tdx_xxxx_xyy_0, \
                                     tdx_xxxx_xyz_0, tdx_xxxx_xzz_0, tdx_xxxx_yyy_0, tdx_xxxx_yyz_0, tdx_xxxx_yzz_0, \
                                     tdx_xxxx_zzz_0, tdx_xxxy_xxx_0, tdx_xxxy_xxy_0, tdx_xxxy_xxz_0, tdx_xxxy_xyy_0, \
                                     tdx_xxxy_xyz_0, tdx_xxxy_xzz_0, tdx_xxxy_yyy_0, tdx_xxy_xx_0, tdx_xxy_xxx_0, \
                                     tdx_xxy_xxy_0, tdx_xxy_xxz_0, tdx_xxy_xy_0, tdx_xxy_xyy_0, tdx_xxy_xyz_0, \
                                     tdx_xxy_xz_0, tdx_xxy_xzz_0, tdx_xxy_yy_0, tdx_xxy_yyy_0, tdx_xxy_yz_0, \
                                     tdx_xxy_zz_0, tdx_xy_xxx_0, tdx_xy_xxy_0, tdx_xy_xxz_0, tdx_xy_xyy_0, tdx_xy_xyz_0, \
                                     tdx_xy_xzz_0, tdx_xy_yyy_0, tdy_xx_xxx_0, tdy_xx_xxy_0, tdy_xx_xxz_0, tdy_xx_xyy_0, \
                                     tdy_xx_xyz_0, tdy_xx_xzz_0, tdy_xx_yyy_0, tdy_xx_yyz_0, tdy_xx_yzz_0, tdy_xx_zzz_0, \
                                     tdy_xxx_xx_0, tdy_xxx_xxx_0, tdy_xxx_xxy_0, tdy_xxx_xxz_0, tdy_xxx_xy_0, \
                                     tdy_xxx_xyy_0, tdy_xxx_xyz_0, tdy_xxx_xz_0, tdy_xxx_xzz_0, tdy_xxx_yy_0, \
                                     tdy_xxx_yyy_0, tdy_xxx_yyz_0, tdy_xxx_yz_0, tdy_xxx_yzz_0, tdy_xxx_zz_0, \
                                     tdy_xxx_zzz_0, tdy_xxxx_xxx_0, tdy_xxxx_xxy_0, tdy_xxxx_xxz_0, tdy_xxxx_xyy_0, \
                                     tdy_xxxx_xyz_0, tdy_xxxx_xzz_0, tdy_xxxx_yyy_0, tdy_xxxx_yyz_0, tdy_xxxx_yzz_0, \
                                     tdy_xxxx_zzz_0, tdy_xxxy_xxx_0, tdy_xxxy_xxy_0, tdy_xxxy_xxz_0, tdy_xxxy_xyy_0, \
                                     tdy_xxxy_xyz_0, tdy_xxxy_xzz_0, tdy_xxxy_yyy_0, tdy_xxy_xx_0, tdy_xxy_xxx_0, \
                                     tdy_xxy_xxy_0, tdy_xxy_xxz_0, tdy_xxy_xy_0, tdy_xxy_xyy_0, tdy_xxy_xyz_0, \
                                     tdy_xxy_xz_0, tdy_xxy_xzz_0, tdy_xxy_yy_0, tdy_xxy_yyy_0, tdy_xxy_yz_0, \
                                     tdy_xxy_zz_0, tdy_xy_xxx_0, tdy_xy_xxy_0, tdy_xy_xxz_0, tdy_xy_xyy_0, tdy_xy_xyz_0, \
                                     tdy_xy_xzz_0, tdy_xy_yyy_0, tdz_xx_xxx_0, tdz_xx_xxy_0, tdz_xx_xxz_0, tdz_xx_xyy_0, \
                                     tdz_xx_xyz_0, tdz_xx_xzz_0, tdz_xx_yyy_0, tdz_xx_yyz_0, tdz_xx_yzz_0, tdz_xx_zzz_0, \
                                     tdz_xxx_xx_0, tdz_xxx_xxx_0, tdz_xxx_xxy_0, tdz_xxx_xxz_0, tdz_xxx_xy_0, \
                                     tdz_xxx_xyy_0, tdz_xxx_xyz_0, tdz_xxx_xz_0, tdz_xxx_xzz_0, tdz_xxx_yy_0, \
                                     tdz_xxx_yyy_0, tdz_xxx_yyz_0, tdz_xxx_yz_0, tdz_xxx_yzz_0, tdz_xxx_zz_0, \
                                     tdz_xxx_zzz_0, tdz_xxxx_xxx_0, tdz_xxxx_xxy_0, tdz_xxxx_xxz_0, tdz_xxxx_xyy_0, \
                                     tdz_xxxx_xyz_0, tdz_xxxx_xzz_0, tdz_xxxx_yyy_0, tdz_xxxx_yyz_0, tdz_xxxx_yzz_0, \
                                     tdz_xxxx_zzz_0, tdz_xxxy_xxx_0, tdz_xxxy_xxy_0, tdz_xxxy_xxz_0, tdz_xxxy_xyy_0, \
                                     tdz_xxxy_xyz_0, tdz_xxxy_xzz_0, tdz_xxy_xx_0, tdz_xxy_xxx_0, tdz_xxy_xxy_0, \
                                     tdz_xxy_xxz_0, tdz_xxy_xy_0, tdz_xxy_xyy_0, tdz_xxy_xyz_0, tdz_xxy_xz_0, \
                                     tdz_xxy_xzz_0, tdz_xxy_yy_0, tdz_xxy_yz_0, tdz_xxy_zz_0, tdz_xy_xxx_0, tdz_xy_xxy_0, \
                                     tdz_xy_xxz_0, tdz_xy_xyy_0, tdz_xy_xyz_0, tdz_xy_xzz_0, ts_xxx_xxx_0, ts_xxx_xxy_0, \
                                     ts_xxx_xxz_0, ts_xxx_xyy_0, ts_xxx_xyz_0, ts_xxx_xzz_0, ts_xxx_yyy_0, ts_xxx_yyz_0, \
                                     ts_xxx_yzz_0, ts_xxx_zzz_0, ts_xxy_xxx_0, ts_xxy_xxy_0, ts_xxy_xxz_0, ts_xxy_xyy_0, \
                                     ts_xxy_xyz_0, ts_xxy_xzz_0, ts_xxy_yyy_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxxx_xxx_0[j] = pa_x[j] * tdx_xxx_xxx_0[j] + 1.5 * fl1_fx * tdx_xx_xxx_0[j] + 1.5 * fl1_fx * tdx_xxx_xx_0[j] + 0.5 * fl1_fx * ts_xxx_xxx_0[j];

                tdy_xxxx_xxx_0[j] = pa_x[j] * tdy_xxx_xxx_0[j] + 1.5 * fl1_fx * tdy_xx_xxx_0[j] + 1.5 * fl1_fx * tdy_xxx_xx_0[j];

                tdz_xxxx_xxx_0[j] = pa_x[j] * tdz_xxx_xxx_0[j] + 1.5 * fl1_fx * tdz_xx_xxx_0[j] + 1.5 * fl1_fx * tdz_xxx_xx_0[j];

                tdx_xxxx_xxy_0[j] = pa_x[j] * tdx_xxx_xxy_0[j] + 1.5 * fl1_fx * tdx_xx_xxy_0[j] + fl1_fx * tdx_xxx_xy_0[j] + 0.5 * fl1_fx * ts_xxx_xxy_0[j];

                tdy_xxxx_xxy_0[j] = pa_x[j] * tdy_xxx_xxy_0[j] + 1.5 * fl1_fx * tdy_xx_xxy_0[j] + fl1_fx * tdy_xxx_xy_0[j];

                tdz_xxxx_xxy_0[j] = pa_x[j] * tdz_xxx_xxy_0[j] + 1.5 * fl1_fx * tdz_xx_xxy_0[j] + fl1_fx * tdz_xxx_xy_0[j];

                tdx_xxxx_xxz_0[j] = pa_x[j] * tdx_xxx_xxz_0[j] + 1.5 * fl1_fx * tdx_xx_xxz_0[j] + fl1_fx * tdx_xxx_xz_0[j] + 0.5 * fl1_fx * ts_xxx_xxz_0[j];

                tdy_xxxx_xxz_0[j] = pa_x[j] * tdy_xxx_xxz_0[j] + 1.5 * fl1_fx * tdy_xx_xxz_0[j] + fl1_fx * tdy_xxx_xz_0[j];

                tdz_xxxx_xxz_0[j] = pa_x[j] * tdz_xxx_xxz_0[j] + 1.5 * fl1_fx * tdz_xx_xxz_0[j] + fl1_fx * tdz_xxx_xz_0[j];

                tdx_xxxx_xyy_0[j] = pa_x[j] * tdx_xxx_xyy_0[j] + 1.5 * fl1_fx * tdx_xx_xyy_0[j] + 0.5 * fl1_fx * tdx_xxx_yy_0[j] + 0.5 * fl1_fx * ts_xxx_xyy_0[j];

                tdy_xxxx_xyy_0[j] = pa_x[j] * tdy_xxx_xyy_0[j] + 1.5 * fl1_fx * tdy_xx_xyy_0[j] + 0.5 * fl1_fx * tdy_xxx_yy_0[j];

                tdz_xxxx_xyy_0[j] = pa_x[j] * tdz_xxx_xyy_0[j] + 1.5 * fl1_fx * tdz_xx_xyy_0[j] + 0.5 * fl1_fx * tdz_xxx_yy_0[j];

                tdx_xxxx_xyz_0[j] = pa_x[j] * tdx_xxx_xyz_0[j] + 1.5 * fl1_fx * tdx_xx_xyz_0[j] + 0.5 * fl1_fx * tdx_xxx_yz_0[j] + 0.5 * fl1_fx * ts_xxx_xyz_0[j];

                tdy_xxxx_xyz_0[j] = pa_x[j] * tdy_xxx_xyz_0[j] + 1.5 * fl1_fx * tdy_xx_xyz_0[j] + 0.5 * fl1_fx * tdy_xxx_yz_0[j];

                tdz_xxxx_xyz_0[j] = pa_x[j] * tdz_xxx_xyz_0[j] + 1.5 * fl1_fx * tdz_xx_xyz_0[j] + 0.5 * fl1_fx * tdz_xxx_yz_0[j];

                tdx_xxxx_xzz_0[j] = pa_x[j] * tdx_xxx_xzz_0[j] + 1.5 * fl1_fx * tdx_xx_xzz_0[j] + 0.5 * fl1_fx * tdx_xxx_zz_0[j] + 0.5 * fl1_fx * ts_xxx_xzz_0[j];

                tdy_xxxx_xzz_0[j] = pa_x[j] * tdy_xxx_xzz_0[j] + 1.5 * fl1_fx * tdy_xx_xzz_0[j] + 0.5 * fl1_fx * tdy_xxx_zz_0[j];

                tdz_xxxx_xzz_0[j] = pa_x[j] * tdz_xxx_xzz_0[j] + 1.5 * fl1_fx * tdz_xx_xzz_0[j] + 0.5 * fl1_fx * tdz_xxx_zz_0[j];

                tdx_xxxx_yyy_0[j] = pa_x[j] * tdx_xxx_yyy_0[j] + 1.5 * fl1_fx * tdx_xx_yyy_0[j] + 0.5 * fl1_fx * ts_xxx_yyy_0[j];

                tdy_xxxx_yyy_0[j] = pa_x[j] * tdy_xxx_yyy_0[j] + 1.5 * fl1_fx * tdy_xx_yyy_0[j];

                tdz_xxxx_yyy_0[j] = pa_x[j] * tdz_xxx_yyy_0[j] + 1.5 * fl1_fx * tdz_xx_yyy_0[j];

                tdx_xxxx_yyz_0[j] = pa_x[j] * tdx_xxx_yyz_0[j] + 1.5 * fl1_fx * tdx_xx_yyz_0[j] + 0.5 * fl1_fx * ts_xxx_yyz_0[j];

                tdy_xxxx_yyz_0[j] = pa_x[j] * tdy_xxx_yyz_0[j] + 1.5 * fl1_fx * tdy_xx_yyz_0[j];

                tdz_xxxx_yyz_0[j] = pa_x[j] * tdz_xxx_yyz_0[j] + 1.5 * fl1_fx * tdz_xx_yyz_0[j];

                tdx_xxxx_yzz_0[j] = pa_x[j] * tdx_xxx_yzz_0[j] + 1.5 * fl1_fx * tdx_xx_yzz_0[j] + 0.5 * fl1_fx * ts_xxx_yzz_0[j];

                tdy_xxxx_yzz_0[j] = pa_x[j] * tdy_xxx_yzz_0[j] + 1.5 * fl1_fx * tdy_xx_yzz_0[j];

                tdz_xxxx_yzz_0[j] = pa_x[j] * tdz_xxx_yzz_0[j] + 1.5 * fl1_fx * tdz_xx_yzz_0[j];

                tdx_xxxx_zzz_0[j] = pa_x[j] * tdx_xxx_zzz_0[j] + 1.5 * fl1_fx * tdx_xx_zzz_0[j] + 0.5 * fl1_fx * ts_xxx_zzz_0[j];

                tdy_xxxx_zzz_0[j] = pa_x[j] * tdy_xxx_zzz_0[j] + 1.5 * fl1_fx * tdy_xx_zzz_0[j];

                tdz_xxxx_zzz_0[j] = pa_x[j] * tdz_xxx_zzz_0[j] + 1.5 * fl1_fx * tdz_xx_zzz_0[j];

                tdx_xxxy_xxx_0[j] = pa_x[j] * tdx_xxy_xxx_0[j] + fl1_fx * tdx_xy_xxx_0[j] + 1.5 * fl1_fx * tdx_xxy_xx_0[j] + 0.5 * fl1_fx * ts_xxy_xxx_0[j];

                tdy_xxxy_xxx_0[j] = pa_x[j] * tdy_xxy_xxx_0[j] + fl1_fx * tdy_xy_xxx_0[j] + 1.5 * fl1_fx * tdy_xxy_xx_0[j];

                tdz_xxxy_xxx_0[j] = pa_x[j] * tdz_xxy_xxx_0[j] + fl1_fx * tdz_xy_xxx_0[j] + 1.5 * fl1_fx * tdz_xxy_xx_0[j];

                tdx_xxxy_xxy_0[j] = pa_x[j] * tdx_xxy_xxy_0[j] + fl1_fx * tdx_xy_xxy_0[j] + fl1_fx * tdx_xxy_xy_0[j] + 0.5 * fl1_fx * ts_xxy_xxy_0[j];

                tdy_xxxy_xxy_0[j] = pa_x[j] * tdy_xxy_xxy_0[j] + fl1_fx * tdy_xy_xxy_0[j] + fl1_fx * tdy_xxy_xy_0[j];

                tdz_xxxy_xxy_0[j] = pa_x[j] * tdz_xxy_xxy_0[j] + fl1_fx * tdz_xy_xxy_0[j] + fl1_fx * tdz_xxy_xy_0[j];

                tdx_xxxy_xxz_0[j] = pa_x[j] * tdx_xxy_xxz_0[j] + fl1_fx * tdx_xy_xxz_0[j] + fl1_fx * tdx_xxy_xz_0[j] + 0.5 * fl1_fx * ts_xxy_xxz_0[j];

                tdy_xxxy_xxz_0[j] = pa_x[j] * tdy_xxy_xxz_0[j] + fl1_fx * tdy_xy_xxz_0[j] + fl1_fx * tdy_xxy_xz_0[j];

                tdz_xxxy_xxz_0[j] = pa_x[j] * tdz_xxy_xxz_0[j] + fl1_fx * tdz_xy_xxz_0[j] + fl1_fx * tdz_xxy_xz_0[j];

                tdx_xxxy_xyy_0[j] = pa_x[j] * tdx_xxy_xyy_0[j] + fl1_fx * tdx_xy_xyy_0[j] + 0.5 * fl1_fx * tdx_xxy_yy_0[j] + 0.5 * fl1_fx * ts_xxy_xyy_0[j];

                tdy_xxxy_xyy_0[j] = pa_x[j] * tdy_xxy_xyy_0[j] + fl1_fx * tdy_xy_xyy_0[j] + 0.5 * fl1_fx * tdy_xxy_yy_0[j];

                tdz_xxxy_xyy_0[j] = pa_x[j] * tdz_xxy_xyy_0[j] + fl1_fx * tdz_xy_xyy_0[j] + 0.5 * fl1_fx * tdz_xxy_yy_0[j];

                tdx_xxxy_xyz_0[j] = pa_x[j] * tdx_xxy_xyz_0[j] + fl1_fx * tdx_xy_xyz_0[j] + 0.5 * fl1_fx * tdx_xxy_yz_0[j] + 0.5 * fl1_fx * ts_xxy_xyz_0[j];

                tdy_xxxy_xyz_0[j] = pa_x[j] * tdy_xxy_xyz_0[j] + fl1_fx * tdy_xy_xyz_0[j] + 0.5 * fl1_fx * tdy_xxy_yz_0[j];

                tdz_xxxy_xyz_0[j] = pa_x[j] * tdz_xxy_xyz_0[j] + fl1_fx * tdz_xy_xyz_0[j] + 0.5 * fl1_fx * tdz_xxy_yz_0[j];

                tdx_xxxy_xzz_0[j] = pa_x[j] * tdx_xxy_xzz_0[j] + fl1_fx * tdx_xy_xzz_0[j] + 0.5 * fl1_fx * tdx_xxy_zz_0[j] + 0.5 * fl1_fx * ts_xxy_xzz_0[j];

                tdy_xxxy_xzz_0[j] = pa_x[j] * tdy_xxy_xzz_0[j] + fl1_fx * tdy_xy_xzz_0[j] + 0.5 * fl1_fx * tdy_xxy_zz_0[j];

                tdz_xxxy_xzz_0[j] = pa_x[j] * tdz_xxy_xzz_0[j] + fl1_fx * tdz_xy_xzz_0[j] + 0.5 * fl1_fx * tdz_xxy_zz_0[j];

                tdx_xxxy_yyy_0[j] = pa_x[j] * tdx_xxy_yyy_0[j] + fl1_fx * tdx_xy_yyy_0[j] + 0.5 * fl1_fx * ts_xxy_yyy_0[j];

                tdy_xxxy_yyy_0[j] = pa_x[j] * tdy_xxy_yyy_0[j] + fl1_fx * tdy_xy_yyy_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGF_50_100(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (50,100)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdz_xxy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 16); 

            auto tdx_xxy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 17); 

            auto tdy_xxy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 17); 

            auto tdz_xxy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 17); 

            auto tdx_xxy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 18); 

            auto tdy_xxy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 18); 

            auto tdz_xxy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 18); 

            auto tdx_xxy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 19); 

            auto tdy_xxy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 19); 

            auto tdz_xxy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 19); 

            auto tdx_xxz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 20); 

            auto tdy_xxz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 20); 

            auto tdz_xxz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 20); 

            auto tdx_xxz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 21); 

            auto tdy_xxz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 21); 

            auto tdz_xxz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 21); 

            auto tdx_xxz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 22); 

            auto tdy_xxz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 22); 

            auto tdz_xxz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 22); 

            auto tdx_xxz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 23); 

            auto tdy_xxz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 23); 

            auto tdz_xxz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 23); 

            auto tdx_xxz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 24); 

            auto tdy_xxz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 24); 

            auto tdz_xxz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 24); 

            auto tdx_xxz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 25); 

            auto tdy_xxz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 25); 

            auto tdz_xxz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 25); 

            auto tdx_xxz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 26); 

            auto tdy_xxz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 26); 

            auto tdz_xxz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 26); 

            auto tdx_xxz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 27); 

            auto tdy_xxz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 27); 

            auto tdz_xxz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 27); 

            auto tdx_xxz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 28); 

            auto tdy_xxz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 28); 

            auto tdz_xxz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 28); 

            auto tdx_xxz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 29); 

            auto tdy_xxz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 29); 

            auto tdz_xxz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 29); 

            auto tdx_xyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 30); 

            auto tdy_xyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 30); 

            auto tdz_xyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 30); 

            auto tdx_xyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 31); 

            auto tdy_xyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 31); 

            auto tdz_xyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 31); 

            auto tdx_xyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 32); 

            auto tdy_xyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 32); 

            auto tdz_xyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 32); 

            auto tdx_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 33); 

            auto tdz_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 16); 

            auto tdx_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 17); 

            auto tdy_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 17); 

            auto tdz_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 17); 

            auto tdx_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 18); 

            auto tdy_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 18); 

            auto tdz_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 18); 

            auto tdx_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 19); 

            auto tdy_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 19); 

            auto tdz_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 19); 

            auto tdx_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 20); 

            auto tdy_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 20); 

            auto tdz_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 20); 

            auto tdx_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 21); 

            auto tdy_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 21); 

            auto tdz_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 21); 

            auto tdx_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 22); 

            auto tdy_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 22); 

            auto tdz_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 22); 

            auto tdx_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 23); 

            auto tdy_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 23); 

            auto tdz_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 23); 

            auto tdx_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 24); 

            auto tdy_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 24); 

            auto tdz_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 24); 

            auto tdx_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 25); 

            auto tdy_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 25); 

            auto tdz_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 25); 

            auto tdx_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 26); 

            auto tdy_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 26); 

            auto tdz_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 26); 

            auto tdx_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 27); 

            auto tdy_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 27); 

            auto tdz_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 27); 

            auto tdx_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 28); 

            auto tdy_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 28); 

            auto tdz_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 28); 

            auto tdx_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 29); 

            auto tdy_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 29); 

            auto tdz_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 29); 

            auto tdx_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 30); 

            auto tdy_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 30); 

            auto tdz_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 30); 

            auto tdx_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 31); 

            auto tdy_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 31); 

            auto tdz_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 31); 

            auto tdx_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 32); 

            auto tdy_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 32); 

            auto tdz_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 32); 

            auto tdx_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 33); 

            auto tdx_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 12); 

            auto tdy_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 12); 

            auto tdz_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 12); 

            auto tdx_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 13); 

            auto tdy_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 13); 

            auto tdz_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 13); 

            auto tdx_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 14); 

            auto tdy_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 14); 

            auto tdz_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 14); 

            auto tdx_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 15); 

            auto tdy_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 15); 

            auto tdz_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 15); 

            auto tdx_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 16); 

            auto tdy_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 16); 

            auto tdz_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 16); 

            auto tdx_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 17); 

            auto tdy_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 17); 

            auto tdz_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 17); 

            auto tdx_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 18); 

            auto tdy_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 18); 

            auto tdz_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 18); 

            auto tdx_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 19); 

            auto tdy_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 19); 

            auto tdz_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 19); 

            auto tdx_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 20); 

            auto tdy_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 20); 

            auto tdz_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 20); 

            auto tdx_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 21); 

            auto ts_xxy_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 17); 

            auto ts_xxy_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 18); 

            auto ts_xxy_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 19); 

            auto ts_xxz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 20); 

            auto ts_xxz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 21); 

            auto ts_xxz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 22); 

            auto ts_xxz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 23); 

            auto ts_xxz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 24); 

            auto ts_xxz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 25); 

            auto ts_xxz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 26); 

            auto ts_xxz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 27); 

            auto ts_xxz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 28); 

            auto ts_xxz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 29); 

            auto ts_xyy_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 30); 

            auto ts_xyy_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 31); 

            auto ts_xyy_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 32); 

            auto ts_xyy_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 33); 

            // set up pointers to integrals

            auto tdz_xxxy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 16); 

            auto tdx_xxxy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 17); 

            auto tdy_xxxy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 17); 

            auto tdz_xxxy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 17); 

            auto tdx_xxxy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 18); 

            auto tdy_xxxy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 18); 

            auto tdz_xxxy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 18); 

            auto tdx_xxxy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 19); 

            auto tdy_xxxy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 19); 

            auto tdz_xxxy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 19); 

            auto tdx_xxxz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 20); 

            auto tdy_xxxz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 20); 

            auto tdz_xxxz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 20); 

            auto tdx_xxxz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 21); 

            auto tdy_xxxz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 21); 

            auto tdz_xxxz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 21); 

            auto tdx_xxxz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 22); 

            auto tdy_xxxz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 22); 

            auto tdz_xxxz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 22); 

            auto tdx_xxxz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 23); 

            auto tdy_xxxz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 23); 

            auto tdz_xxxz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 23); 

            auto tdx_xxxz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 24); 

            auto tdy_xxxz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 24); 

            auto tdz_xxxz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 24); 

            auto tdx_xxxz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 25); 

            auto tdy_xxxz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 25); 

            auto tdz_xxxz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 25); 

            auto tdx_xxxz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 26); 

            auto tdy_xxxz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 26); 

            auto tdz_xxxz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 26); 

            auto tdx_xxxz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 27); 

            auto tdy_xxxz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 27); 

            auto tdz_xxxz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 27); 

            auto tdx_xxxz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 28); 

            auto tdy_xxxz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 28); 

            auto tdz_xxxz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 28); 

            auto tdx_xxxz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 29); 

            auto tdy_xxxz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 29); 

            auto tdz_xxxz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 29); 

            auto tdx_xxyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 30); 

            auto tdy_xxyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 30); 

            auto tdz_xxyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 30); 

            auto tdx_xxyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 31); 

            auto tdy_xxyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 31); 

            auto tdz_xxyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 31); 

            auto tdx_xxyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 32); 

            auto tdy_xxyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 32); 

            auto tdz_xxyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 32); 

            auto tdx_xxyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 33); 

            // Batch of Integrals (50,100)

            #pragma omp simd aligned(fx, pa_x, tdx_xxxy_yyz_0, tdx_xxxy_yzz_0, tdx_xxxy_zzz_0, \
                                     tdx_xxxz_xxx_0, tdx_xxxz_xxy_0, tdx_xxxz_xxz_0, tdx_xxxz_xyy_0, tdx_xxxz_xyz_0, \
                                     tdx_xxxz_xzz_0, tdx_xxxz_yyy_0, tdx_xxxz_yyz_0, tdx_xxxz_yzz_0, tdx_xxxz_zzz_0, \
                                     tdx_xxy_yyz_0, tdx_xxy_yzz_0, tdx_xxy_zzz_0, tdx_xxyy_xxx_0, tdx_xxyy_xxy_0, \
                                     tdx_xxyy_xxz_0, tdx_xxyy_xyy_0, tdx_xxz_xx_0, tdx_xxz_xxx_0, tdx_xxz_xxy_0, \
                                     tdx_xxz_xxz_0, tdx_xxz_xy_0, tdx_xxz_xyy_0, tdx_xxz_xyz_0, tdx_xxz_xz_0, \
                                     tdx_xxz_xzz_0, tdx_xxz_yy_0, tdx_xxz_yyy_0, tdx_xxz_yyz_0, tdx_xxz_yz_0, \
                                     tdx_xxz_yzz_0, tdx_xxz_zz_0, tdx_xxz_zzz_0, tdx_xy_yyz_0, tdx_xy_yzz_0, \
                                     tdx_xy_zzz_0, tdx_xyy_xx_0, tdx_xyy_xxx_0, tdx_xyy_xxy_0, tdx_xyy_xxz_0, \
                                     tdx_xyy_xy_0, tdx_xyy_xyy_0, tdx_xyy_xz_0, tdx_xyy_yy_0, tdx_xz_xxx_0, \
                                     tdx_xz_xxy_0, tdx_xz_xxz_0, tdx_xz_xyy_0, tdx_xz_xyz_0, tdx_xz_xzz_0, tdx_xz_yyy_0, \
                                     tdx_xz_yyz_0, tdx_xz_yzz_0, tdx_xz_zzz_0, tdx_yy_xxx_0, tdx_yy_xxy_0, tdx_yy_xxz_0, \
                                     tdx_yy_xyy_0, tdy_xxxy_yyz_0, tdy_xxxy_yzz_0, tdy_xxxy_zzz_0, tdy_xxxz_xxx_0, \
                                     tdy_xxxz_xxy_0, tdy_xxxz_xxz_0, tdy_xxxz_xyy_0, tdy_xxxz_xyz_0, tdy_xxxz_xzz_0, \
                                     tdy_xxxz_yyy_0, tdy_xxxz_yyz_0, tdy_xxxz_yzz_0, tdy_xxxz_zzz_0, tdy_xxy_yyz_0, \
                                     tdy_xxy_yzz_0, tdy_xxy_zzz_0, tdy_xxyy_xxx_0, tdy_xxyy_xxy_0, tdy_xxyy_xxz_0, \
                                     tdy_xxz_xx_0, tdy_xxz_xxx_0, tdy_xxz_xxy_0, tdy_xxz_xxz_0, tdy_xxz_xy_0, \
                                     tdy_xxz_xyy_0, tdy_xxz_xyz_0, tdy_xxz_xz_0, tdy_xxz_xzz_0, tdy_xxz_yy_0, \
                                     tdy_xxz_yyy_0, tdy_xxz_yyz_0, tdy_xxz_yz_0, tdy_xxz_yzz_0, tdy_xxz_zz_0, \
                                     tdy_xxz_zzz_0, tdy_xy_yyz_0, tdy_xy_yzz_0, tdy_xy_zzz_0, tdy_xyy_xx_0, \
                                     tdy_xyy_xxx_0, tdy_xyy_xxy_0, tdy_xyy_xxz_0, tdy_xyy_xy_0, tdy_xyy_xz_0, \
                                     tdy_xz_xxx_0, tdy_xz_xxy_0, tdy_xz_xxz_0, tdy_xz_xyy_0, tdy_xz_xyz_0, tdy_xz_xzz_0, \
                                     tdy_xz_yyy_0, tdy_xz_yyz_0, tdy_xz_yzz_0, tdy_xz_zzz_0, tdy_yy_xxx_0, tdy_yy_xxy_0, \
                                     tdy_yy_xxz_0, tdz_xxxy_yyy_0, tdz_xxxy_yyz_0, tdz_xxxy_yzz_0, tdz_xxxy_zzz_0, \
                                     tdz_xxxz_xxx_0, tdz_xxxz_xxy_0, tdz_xxxz_xxz_0, tdz_xxxz_xyy_0, tdz_xxxz_xyz_0, \
                                     tdz_xxxz_xzz_0, tdz_xxxz_yyy_0, tdz_xxxz_yyz_0, tdz_xxxz_yzz_0, tdz_xxxz_zzz_0, \
                                     tdz_xxy_yyy_0, tdz_xxy_yyz_0, tdz_xxy_yzz_0, tdz_xxy_zzz_0, tdz_xxyy_xxx_0, \
                                     tdz_xxyy_xxy_0, tdz_xxyy_xxz_0, tdz_xxz_xx_0, tdz_xxz_xxx_0, tdz_xxz_xxy_0, \
                                     tdz_xxz_xxz_0, tdz_xxz_xy_0, tdz_xxz_xyy_0, tdz_xxz_xyz_0, tdz_xxz_xz_0, \
                                     tdz_xxz_xzz_0, tdz_xxz_yy_0, tdz_xxz_yyy_0, tdz_xxz_yyz_0, tdz_xxz_yz_0, \
                                     tdz_xxz_yzz_0, tdz_xxz_zz_0, tdz_xxz_zzz_0, tdz_xy_yyy_0, tdz_xy_yyz_0, \
                                     tdz_xy_yzz_0, tdz_xy_zzz_0, tdz_xyy_xx_0, tdz_xyy_xxx_0, tdz_xyy_xxy_0, \
                                     tdz_xyy_xxz_0, tdz_xyy_xy_0, tdz_xyy_xz_0, tdz_xz_xxx_0, tdz_xz_xxy_0, tdz_xz_xxz_0, \
                                     tdz_xz_xyy_0, tdz_xz_xyz_0, tdz_xz_xzz_0, tdz_xz_yyy_0, tdz_xz_yyz_0, tdz_xz_yzz_0, \
                                     tdz_xz_zzz_0, tdz_yy_xxx_0, tdz_yy_xxy_0, tdz_yy_xxz_0, ts_xxy_yyz_0, ts_xxy_yzz_0, \
                                     ts_xxy_zzz_0, ts_xxz_xxx_0, ts_xxz_xxy_0, ts_xxz_xxz_0, ts_xxz_xyy_0, ts_xxz_xyz_0, \
                                     ts_xxz_xzz_0, ts_xxz_yyy_0, ts_xxz_yyz_0, ts_xxz_yzz_0, ts_xxz_zzz_0, ts_xyy_xxx_0, \
                                     ts_xyy_xxy_0, ts_xyy_xxz_0, ts_xyy_xyy_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdz_xxxy_yyy_0[j] = pa_x[j] * tdz_xxy_yyy_0[j] + fl1_fx * tdz_xy_yyy_0[j];

                tdx_xxxy_yyz_0[j] = pa_x[j] * tdx_xxy_yyz_0[j] + fl1_fx * tdx_xy_yyz_0[j] + 0.5 * fl1_fx * ts_xxy_yyz_0[j];

                tdy_xxxy_yyz_0[j] = pa_x[j] * tdy_xxy_yyz_0[j] + fl1_fx * tdy_xy_yyz_0[j];

                tdz_xxxy_yyz_0[j] = pa_x[j] * tdz_xxy_yyz_0[j] + fl1_fx * tdz_xy_yyz_0[j];

                tdx_xxxy_yzz_0[j] = pa_x[j] * tdx_xxy_yzz_0[j] + fl1_fx * tdx_xy_yzz_0[j] + 0.5 * fl1_fx * ts_xxy_yzz_0[j];

                tdy_xxxy_yzz_0[j] = pa_x[j] * tdy_xxy_yzz_0[j] + fl1_fx * tdy_xy_yzz_0[j];

                tdz_xxxy_yzz_0[j] = pa_x[j] * tdz_xxy_yzz_0[j] + fl1_fx * tdz_xy_yzz_0[j];

                tdx_xxxy_zzz_0[j] = pa_x[j] * tdx_xxy_zzz_0[j] + fl1_fx * tdx_xy_zzz_0[j] + 0.5 * fl1_fx * ts_xxy_zzz_0[j];

                tdy_xxxy_zzz_0[j] = pa_x[j] * tdy_xxy_zzz_0[j] + fl1_fx * tdy_xy_zzz_0[j];

                tdz_xxxy_zzz_0[j] = pa_x[j] * tdz_xxy_zzz_0[j] + fl1_fx * tdz_xy_zzz_0[j];

                tdx_xxxz_xxx_0[j] = pa_x[j] * tdx_xxz_xxx_0[j] + fl1_fx * tdx_xz_xxx_0[j] + 1.5 * fl1_fx * tdx_xxz_xx_0[j] + 0.5 * fl1_fx * ts_xxz_xxx_0[j];

                tdy_xxxz_xxx_0[j] = pa_x[j] * tdy_xxz_xxx_0[j] + fl1_fx * tdy_xz_xxx_0[j] + 1.5 * fl1_fx * tdy_xxz_xx_0[j];

                tdz_xxxz_xxx_0[j] = pa_x[j] * tdz_xxz_xxx_0[j] + fl1_fx * tdz_xz_xxx_0[j] + 1.5 * fl1_fx * tdz_xxz_xx_0[j];

                tdx_xxxz_xxy_0[j] = pa_x[j] * tdx_xxz_xxy_0[j] + fl1_fx * tdx_xz_xxy_0[j] + fl1_fx * tdx_xxz_xy_0[j] + 0.5 * fl1_fx * ts_xxz_xxy_0[j];

                tdy_xxxz_xxy_0[j] = pa_x[j] * tdy_xxz_xxy_0[j] + fl1_fx * tdy_xz_xxy_0[j] + fl1_fx * tdy_xxz_xy_0[j];

                tdz_xxxz_xxy_0[j] = pa_x[j] * tdz_xxz_xxy_0[j] + fl1_fx * tdz_xz_xxy_0[j] + fl1_fx * tdz_xxz_xy_0[j];

                tdx_xxxz_xxz_0[j] = pa_x[j] * tdx_xxz_xxz_0[j] + fl1_fx * tdx_xz_xxz_0[j] + fl1_fx * tdx_xxz_xz_0[j] + 0.5 * fl1_fx * ts_xxz_xxz_0[j];

                tdy_xxxz_xxz_0[j] = pa_x[j] * tdy_xxz_xxz_0[j] + fl1_fx * tdy_xz_xxz_0[j] + fl1_fx * tdy_xxz_xz_0[j];

                tdz_xxxz_xxz_0[j] = pa_x[j] * tdz_xxz_xxz_0[j] + fl1_fx * tdz_xz_xxz_0[j] + fl1_fx * tdz_xxz_xz_0[j];

                tdx_xxxz_xyy_0[j] = pa_x[j] * tdx_xxz_xyy_0[j] + fl1_fx * tdx_xz_xyy_0[j] + 0.5 * fl1_fx * tdx_xxz_yy_0[j] + 0.5 * fl1_fx * ts_xxz_xyy_0[j];

                tdy_xxxz_xyy_0[j] = pa_x[j] * tdy_xxz_xyy_0[j] + fl1_fx * tdy_xz_xyy_0[j] + 0.5 * fl1_fx * tdy_xxz_yy_0[j];

                tdz_xxxz_xyy_0[j] = pa_x[j] * tdz_xxz_xyy_0[j] + fl1_fx * tdz_xz_xyy_0[j] + 0.5 * fl1_fx * tdz_xxz_yy_0[j];

                tdx_xxxz_xyz_0[j] = pa_x[j] * tdx_xxz_xyz_0[j] + fl1_fx * tdx_xz_xyz_0[j] + 0.5 * fl1_fx * tdx_xxz_yz_0[j] + 0.5 * fl1_fx * ts_xxz_xyz_0[j];

                tdy_xxxz_xyz_0[j] = pa_x[j] * tdy_xxz_xyz_0[j] + fl1_fx * tdy_xz_xyz_0[j] + 0.5 * fl1_fx * tdy_xxz_yz_0[j];

                tdz_xxxz_xyz_0[j] = pa_x[j] * tdz_xxz_xyz_0[j] + fl1_fx * tdz_xz_xyz_0[j] + 0.5 * fl1_fx * tdz_xxz_yz_0[j];

                tdx_xxxz_xzz_0[j] = pa_x[j] * tdx_xxz_xzz_0[j] + fl1_fx * tdx_xz_xzz_0[j] + 0.5 * fl1_fx * tdx_xxz_zz_0[j] + 0.5 * fl1_fx * ts_xxz_xzz_0[j];

                tdy_xxxz_xzz_0[j] = pa_x[j] * tdy_xxz_xzz_0[j] + fl1_fx * tdy_xz_xzz_0[j] + 0.5 * fl1_fx * tdy_xxz_zz_0[j];

                tdz_xxxz_xzz_0[j] = pa_x[j] * tdz_xxz_xzz_0[j] + fl1_fx * tdz_xz_xzz_0[j] + 0.5 * fl1_fx * tdz_xxz_zz_0[j];

                tdx_xxxz_yyy_0[j] = pa_x[j] * tdx_xxz_yyy_0[j] + fl1_fx * tdx_xz_yyy_0[j] + 0.5 * fl1_fx * ts_xxz_yyy_0[j];

                tdy_xxxz_yyy_0[j] = pa_x[j] * tdy_xxz_yyy_0[j] + fl1_fx * tdy_xz_yyy_0[j];

                tdz_xxxz_yyy_0[j] = pa_x[j] * tdz_xxz_yyy_0[j] + fl1_fx * tdz_xz_yyy_0[j];

                tdx_xxxz_yyz_0[j] = pa_x[j] * tdx_xxz_yyz_0[j] + fl1_fx * tdx_xz_yyz_0[j] + 0.5 * fl1_fx * ts_xxz_yyz_0[j];

                tdy_xxxz_yyz_0[j] = pa_x[j] * tdy_xxz_yyz_0[j] + fl1_fx * tdy_xz_yyz_0[j];

                tdz_xxxz_yyz_0[j] = pa_x[j] * tdz_xxz_yyz_0[j] + fl1_fx * tdz_xz_yyz_0[j];

                tdx_xxxz_yzz_0[j] = pa_x[j] * tdx_xxz_yzz_0[j] + fl1_fx * tdx_xz_yzz_0[j] + 0.5 * fl1_fx * ts_xxz_yzz_0[j];

                tdy_xxxz_yzz_0[j] = pa_x[j] * tdy_xxz_yzz_0[j] + fl1_fx * tdy_xz_yzz_0[j];

                tdz_xxxz_yzz_0[j] = pa_x[j] * tdz_xxz_yzz_0[j] + fl1_fx * tdz_xz_yzz_0[j];

                tdx_xxxz_zzz_0[j] = pa_x[j] * tdx_xxz_zzz_0[j] + fl1_fx * tdx_xz_zzz_0[j] + 0.5 * fl1_fx * ts_xxz_zzz_0[j];

                tdy_xxxz_zzz_0[j] = pa_x[j] * tdy_xxz_zzz_0[j] + fl1_fx * tdy_xz_zzz_0[j];

                tdz_xxxz_zzz_0[j] = pa_x[j] * tdz_xxz_zzz_0[j] + fl1_fx * tdz_xz_zzz_0[j];

                tdx_xxyy_xxx_0[j] = pa_x[j] * tdx_xyy_xxx_0[j] + 0.5 * fl1_fx * tdx_yy_xxx_0[j] + 1.5 * fl1_fx * tdx_xyy_xx_0[j] + 0.5 * fl1_fx * ts_xyy_xxx_0[j];

                tdy_xxyy_xxx_0[j] = pa_x[j] * tdy_xyy_xxx_0[j] + 0.5 * fl1_fx * tdy_yy_xxx_0[j] + 1.5 * fl1_fx * tdy_xyy_xx_0[j];

                tdz_xxyy_xxx_0[j] = pa_x[j] * tdz_xyy_xxx_0[j] + 0.5 * fl1_fx * tdz_yy_xxx_0[j] + 1.5 * fl1_fx * tdz_xyy_xx_0[j];

                tdx_xxyy_xxy_0[j] = pa_x[j] * tdx_xyy_xxy_0[j] + 0.5 * fl1_fx * tdx_yy_xxy_0[j] + fl1_fx * tdx_xyy_xy_0[j] + 0.5 * fl1_fx * ts_xyy_xxy_0[j];

                tdy_xxyy_xxy_0[j] = pa_x[j] * tdy_xyy_xxy_0[j] + 0.5 * fl1_fx * tdy_yy_xxy_0[j] + fl1_fx * tdy_xyy_xy_0[j];

                tdz_xxyy_xxy_0[j] = pa_x[j] * tdz_xyy_xxy_0[j] + 0.5 * fl1_fx * tdz_yy_xxy_0[j] + fl1_fx * tdz_xyy_xy_0[j];

                tdx_xxyy_xxz_0[j] = pa_x[j] * tdx_xyy_xxz_0[j] + 0.5 * fl1_fx * tdx_yy_xxz_0[j] + fl1_fx * tdx_xyy_xz_0[j] + 0.5 * fl1_fx * ts_xyy_xxz_0[j];

                tdy_xxyy_xxz_0[j] = pa_x[j] * tdy_xyy_xxz_0[j] + 0.5 * fl1_fx * tdy_yy_xxz_0[j] + fl1_fx * tdy_xyy_xz_0[j];

                tdz_xxyy_xxz_0[j] = pa_x[j] * tdz_xyy_xxz_0[j] + 0.5 * fl1_fx * tdz_yy_xxz_0[j] + fl1_fx * tdz_xyy_xz_0[j];

                tdx_xxyy_xyy_0[j] = pa_x[j] * tdx_xyy_xyy_0[j] + 0.5 * fl1_fx * tdx_yy_xyy_0[j] + 0.5 * fl1_fx * tdx_xyy_yy_0[j] + 0.5 * fl1_fx * ts_xyy_xyy_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGF_100_150(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (100,150)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdy_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 33); 

            auto tdz_xyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 33); 

            auto tdx_xyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 34); 

            auto tdy_xyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 34); 

            auto tdz_xyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 34); 

            auto tdx_xyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 35); 

            auto tdy_xyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 35); 

            auto tdz_xyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 35); 

            auto tdx_xyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 36); 

            auto tdy_xyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 36); 

            auto tdz_xyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 36); 

            auto tdx_xyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 37); 

            auto tdy_xyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 37); 

            auto tdz_xyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 37); 

            auto tdx_xyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 38); 

            auto tdy_xyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 38); 

            auto tdz_xyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 38); 

            auto tdx_xyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 39); 

            auto tdy_xyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 39); 

            auto tdz_xyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 39); 

            auto tdx_xyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 40); 

            auto tdy_xyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 40); 

            auto tdz_xyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 40); 

            auto tdx_xyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 41); 

            auto tdy_xyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 41); 

            auto tdz_xyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 41); 

            auto tdx_xyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 42); 

            auto tdy_xyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 42); 

            auto tdz_xyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 42); 

            auto tdx_xyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 43); 

            auto tdy_xyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 43); 

            auto tdz_xyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 43); 

            auto tdx_xyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 44); 

            auto tdy_xyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 44); 

            auto tdz_xyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 44); 

            auto tdx_xyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 45); 

            auto tdy_xyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 45); 

            auto tdz_xyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 45); 

            auto tdx_xyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 46); 

            auto tdy_xyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 46); 

            auto tdz_xyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 46); 

            auto tdx_xyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 47); 

            auto tdy_xyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 47); 

            auto tdz_xyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 47); 

            auto tdx_xyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 48); 

            auto tdy_xyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 48); 

            auto tdz_xyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 48); 

            auto tdx_xyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 49); 

            auto tdy_xyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 49); 

            auto tdz_xyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 49); 

            auto tdy_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 33); 

            auto tdz_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 33); 

            auto tdx_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 34); 

            auto tdy_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 34); 

            auto tdz_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 34); 

            auto tdx_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 35); 

            auto tdy_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 35); 

            auto tdz_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 35); 

            auto tdx_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 36); 

            auto tdy_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 37); 

            auto tdy_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 38); 

            auto tdy_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 39); 

            auto tdy_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 40); 

            auto tdy_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 41); 

            auto tdy_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdx_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 42); 

            auto tdy_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 43); 

            auto tdy_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 44); 

            auto tdy_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 44); 

            auto tdx_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 45); 

            auto tdy_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 45); 

            auto tdz_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 46); 

            auto tdy_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdz_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 47); 

            auto tdy_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 48); 

            auto tdy_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 49); 

            auto tdy_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdy_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 21); 

            auto tdz_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 21); 

            auto tdx_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 22); 

            auto tdy_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 22); 

            auto tdz_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 22); 

            auto tdx_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 23); 

            auto tdy_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 23); 

            auto tdz_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 23); 

            auto tdx_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 24); 

            auto tdy_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 24); 

            auto tdz_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 24); 

            auto tdx_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 25); 

            auto tdy_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 25); 

            auto tdz_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 25); 

            auto tdx_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 26); 

            auto tdy_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 26); 

            auto tdz_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 26); 

            auto tdx_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 27); 

            auto tdy_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 27); 

            auto tdz_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 27); 

            auto tdx_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 28); 

            auto tdy_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 28); 

            auto tdz_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 28); 

            auto tdx_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 29); 

            auto tdy_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 29); 

            auto tdz_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 29); 

            auto ts_xyy_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 34); 

            auto ts_xyy_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 35); 

            auto ts_xyy_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 36); 

            auto ts_xyy_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 37); 

            auto ts_xyy_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 38); 

            auto ts_xyy_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 39); 

            auto ts_xyz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 40); 

            auto ts_xyz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 41); 

            auto ts_xyz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 42); 

            auto ts_xyz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 43); 

            auto ts_xyz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 44); 

            auto ts_xyz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 45); 

            auto ts_xyz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 46); 

            auto ts_xyz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 47); 

            auto ts_xyz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 48); 

            auto ts_xyz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 49); 

            // set up pointers to integrals

            auto tdy_xxyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 33); 

            auto tdz_xxyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 33); 

            auto tdx_xxyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 34); 

            auto tdy_xxyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 34); 

            auto tdz_xxyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 34); 

            auto tdx_xxyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 35); 

            auto tdy_xxyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 35); 

            auto tdz_xxyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 35); 

            auto tdx_xxyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 36); 

            auto tdy_xxyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 36); 

            auto tdz_xxyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 36); 

            auto tdx_xxyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 37); 

            auto tdy_xxyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 37); 

            auto tdz_xxyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 37); 

            auto tdx_xxyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 38); 

            auto tdy_xxyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 38); 

            auto tdz_xxyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 38); 

            auto tdx_xxyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 39); 

            auto tdy_xxyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 39); 

            auto tdz_xxyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 39); 

            auto tdx_xxyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 40); 

            auto tdy_xxyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 40); 

            auto tdz_xxyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 40); 

            auto tdx_xxyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 41); 

            auto tdy_xxyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 41); 

            auto tdz_xxyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 41); 

            auto tdx_xxyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 42); 

            auto tdy_xxyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 42); 

            auto tdz_xxyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 42); 

            auto tdx_xxyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 43); 

            auto tdy_xxyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 43); 

            auto tdz_xxyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 43); 

            auto tdx_xxyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 44); 

            auto tdy_xxyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 44); 

            auto tdz_xxyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 44); 

            auto tdx_xxyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 45); 

            auto tdy_xxyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 45); 

            auto tdz_xxyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 45); 

            auto tdx_xxyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 46); 

            auto tdy_xxyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 46); 

            auto tdz_xxyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 46); 

            auto tdx_xxyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 47); 

            auto tdy_xxyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 47); 

            auto tdz_xxyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 47); 

            auto tdx_xxyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 48); 

            auto tdy_xxyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 48); 

            auto tdz_xxyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 48); 

            auto tdx_xxyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 49); 

            auto tdy_xxyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 49); 

            auto tdz_xxyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 49); 

            // Batch of Integrals (100,150)

            #pragma omp simd aligned(fx, pa_x, tdx_xxyy_xyz_0, tdx_xxyy_xzz_0, tdx_xxyy_yyy_0, \
                                     tdx_xxyy_yyz_0, tdx_xxyy_yzz_0, tdx_xxyy_zzz_0, tdx_xxyz_xxx_0, tdx_xxyz_xxy_0, \
                                     tdx_xxyz_xxz_0, tdx_xxyz_xyy_0, tdx_xxyz_xyz_0, tdx_xxyz_xzz_0, tdx_xxyz_yyy_0, \
                                     tdx_xxyz_yyz_0, tdx_xxyz_yzz_0, tdx_xxyz_zzz_0, tdx_xyy_xyz_0, tdx_xyy_xzz_0, \
                                     tdx_xyy_yyy_0, tdx_xyy_yyz_0, tdx_xyy_yz_0, tdx_xyy_yzz_0, tdx_xyy_zz_0, \
                                     tdx_xyy_zzz_0, tdx_xyz_xx_0, tdx_xyz_xxx_0, tdx_xyz_xxy_0, tdx_xyz_xxz_0, \
                                     tdx_xyz_xy_0, tdx_xyz_xyy_0, tdx_xyz_xyz_0, tdx_xyz_xz_0, tdx_xyz_xzz_0, \
                                     tdx_xyz_yy_0, tdx_xyz_yyy_0, tdx_xyz_yyz_0, tdx_xyz_yz_0, tdx_xyz_yzz_0, \
                                     tdx_xyz_zz_0, tdx_xyz_zzz_0, tdx_yy_xyz_0, tdx_yy_xzz_0, tdx_yy_yyy_0, \
                                     tdx_yy_yyz_0, tdx_yy_yzz_0, tdx_yy_zzz_0, tdx_yz_xxx_0, tdx_yz_xxy_0, tdx_yz_xxz_0, \
                                     tdx_yz_xyy_0, tdx_yz_xyz_0, tdx_yz_xzz_0, tdx_yz_yyy_0, tdx_yz_yyz_0, tdx_yz_yzz_0, \
                                     tdx_yz_zzz_0, tdy_xxyy_xyy_0, tdy_xxyy_xyz_0, tdy_xxyy_xzz_0, tdy_xxyy_yyy_0, \
                                     tdy_xxyy_yyz_0, tdy_xxyy_yzz_0, tdy_xxyy_zzz_0, tdy_xxyz_xxx_0, tdy_xxyz_xxy_0, \
                                     tdy_xxyz_xxz_0, tdy_xxyz_xyy_0, tdy_xxyz_xyz_0, tdy_xxyz_xzz_0, tdy_xxyz_yyy_0, \
                                     tdy_xxyz_yyz_0, tdy_xxyz_yzz_0, tdy_xxyz_zzz_0, tdy_xyy_xyy_0, tdy_xyy_xyz_0, \
                                     tdy_xyy_xzz_0, tdy_xyy_yy_0, tdy_xyy_yyy_0, tdy_xyy_yyz_0, tdy_xyy_yz_0, \
                                     tdy_xyy_yzz_0, tdy_xyy_zz_0, tdy_xyy_zzz_0, tdy_xyz_xx_0, tdy_xyz_xxx_0, \
                                     tdy_xyz_xxy_0, tdy_xyz_xxz_0, tdy_xyz_xy_0, tdy_xyz_xyy_0, tdy_xyz_xyz_0, \
                                     tdy_xyz_xz_0, tdy_xyz_xzz_0, tdy_xyz_yy_0, tdy_xyz_yyy_0, tdy_xyz_yyz_0, \
                                     tdy_xyz_yz_0, tdy_xyz_yzz_0, tdy_xyz_zz_0, tdy_xyz_zzz_0, tdy_yy_xyy_0, \
                                     tdy_yy_xyz_0, tdy_yy_xzz_0, tdy_yy_yyy_0, tdy_yy_yyz_0, tdy_yy_yzz_0, tdy_yy_zzz_0, \
                                     tdy_yz_xxx_0, tdy_yz_xxy_0, tdy_yz_xxz_0, tdy_yz_xyy_0, tdy_yz_xyz_0, tdy_yz_xzz_0, \
                                     tdy_yz_yyy_0, tdy_yz_yyz_0, tdy_yz_yzz_0, tdy_yz_zzz_0, tdz_xxyy_xyy_0, \
                                     tdz_xxyy_xyz_0, tdz_xxyy_xzz_0, tdz_xxyy_yyy_0, tdz_xxyy_yyz_0, tdz_xxyy_yzz_0, \
                                     tdz_xxyy_zzz_0, tdz_xxyz_xxx_0, tdz_xxyz_xxy_0, tdz_xxyz_xxz_0, tdz_xxyz_xyy_0, \
                                     tdz_xxyz_xyz_0, tdz_xxyz_xzz_0, tdz_xxyz_yyy_0, tdz_xxyz_yyz_0, tdz_xxyz_yzz_0, \
                                     tdz_xxyz_zzz_0, tdz_xyy_xyy_0, tdz_xyy_xyz_0, tdz_xyy_xzz_0, tdz_xyy_yy_0, \
                                     tdz_xyy_yyy_0, tdz_xyy_yyz_0, tdz_xyy_yz_0, tdz_xyy_yzz_0, tdz_xyy_zz_0, \
                                     tdz_xyy_zzz_0, tdz_xyz_xx_0, tdz_xyz_xxx_0, tdz_xyz_xxy_0, tdz_xyz_xxz_0, \
                                     tdz_xyz_xy_0, tdz_xyz_xyy_0, tdz_xyz_xyz_0, tdz_xyz_xz_0, tdz_xyz_xzz_0, \
                                     tdz_xyz_yy_0, tdz_xyz_yyy_0, tdz_xyz_yyz_0, tdz_xyz_yz_0, tdz_xyz_yzz_0, \
                                     tdz_xyz_zz_0, tdz_xyz_zzz_0, tdz_yy_xyy_0, tdz_yy_xyz_0, tdz_yy_xzz_0, \
                                     tdz_yy_yyy_0, tdz_yy_yyz_0, tdz_yy_yzz_0, tdz_yy_zzz_0, tdz_yz_xxx_0, tdz_yz_xxy_0, \
                                     tdz_yz_xxz_0, tdz_yz_xyy_0, tdz_yz_xyz_0, tdz_yz_xzz_0, tdz_yz_yyy_0, tdz_yz_yyz_0, \
                                     tdz_yz_yzz_0, tdz_yz_zzz_0, ts_xyy_xyz_0, ts_xyy_xzz_0, ts_xyy_yyy_0, ts_xyy_yyz_0, \
                                     ts_xyy_yzz_0, ts_xyy_zzz_0, ts_xyz_xxx_0, ts_xyz_xxy_0, ts_xyz_xxz_0, ts_xyz_xyy_0, \
                                     ts_xyz_xyz_0, ts_xyz_xzz_0, ts_xyz_yyy_0, ts_xyz_yyz_0, ts_xyz_yzz_0, ts_xyz_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdy_xxyy_xyy_0[j] = pa_x[j] * tdy_xyy_xyy_0[j] + 0.5 * fl1_fx * tdy_yy_xyy_0[j] + 0.5 * fl1_fx * tdy_xyy_yy_0[j];

                tdz_xxyy_xyy_0[j] = pa_x[j] * tdz_xyy_xyy_0[j] + 0.5 * fl1_fx * tdz_yy_xyy_0[j] + 0.5 * fl1_fx * tdz_xyy_yy_0[j];

                tdx_xxyy_xyz_0[j] = pa_x[j] * tdx_xyy_xyz_0[j] + 0.5 * fl1_fx * tdx_yy_xyz_0[j] + 0.5 * fl1_fx * tdx_xyy_yz_0[j] + 0.5 * fl1_fx * ts_xyy_xyz_0[j];

                tdy_xxyy_xyz_0[j] = pa_x[j] * tdy_xyy_xyz_0[j] + 0.5 * fl1_fx * tdy_yy_xyz_0[j] + 0.5 * fl1_fx * tdy_xyy_yz_0[j];

                tdz_xxyy_xyz_0[j] = pa_x[j] * tdz_xyy_xyz_0[j] + 0.5 * fl1_fx * tdz_yy_xyz_0[j] + 0.5 * fl1_fx * tdz_xyy_yz_0[j];

                tdx_xxyy_xzz_0[j] = pa_x[j] * tdx_xyy_xzz_0[j] + 0.5 * fl1_fx * tdx_yy_xzz_0[j] + 0.5 * fl1_fx * tdx_xyy_zz_0[j] + 0.5 * fl1_fx * ts_xyy_xzz_0[j];

                tdy_xxyy_xzz_0[j] = pa_x[j] * tdy_xyy_xzz_0[j] + 0.5 * fl1_fx * tdy_yy_xzz_0[j] + 0.5 * fl1_fx * tdy_xyy_zz_0[j];

                tdz_xxyy_xzz_0[j] = pa_x[j] * tdz_xyy_xzz_0[j] + 0.5 * fl1_fx * tdz_yy_xzz_0[j] + 0.5 * fl1_fx * tdz_xyy_zz_0[j];

                tdx_xxyy_yyy_0[j] = pa_x[j] * tdx_xyy_yyy_0[j] + 0.5 * fl1_fx * tdx_yy_yyy_0[j] + 0.5 * fl1_fx * ts_xyy_yyy_0[j];

                tdy_xxyy_yyy_0[j] = pa_x[j] * tdy_xyy_yyy_0[j] + 0.5 * fl1_fx * tdy_yy_yyy_0[j];

                tdz_xxyy_yyy_0[j] = pa_x[j] * tdz_xyy_yyy_0[j] + 0.5 * fl1_fx * tdz_yy_yyy_0[j];

                tdx_xxyy_yyz_0[j] = pa_x[j] * tdx_xyy_yyz_0[j] + 0.5 * fl1_fx * tdx_yy_yyz_0[j] + 0.5 * fl1_fx * ts_xyy_yyz_0[j];

                tdy_xxyy_yyz_0[j] = pa_x[j] * tdy_xyy_yyz_0[j] + 0.5 * fl1_fx * tdy_yy_yyz_0[j];

                tdz_xxyy_yyz_0[j] = pa_x[j] * tdz_xyy_yyz_0[j] + 0.5 * fl1_fx * tdz_yy_yyz_0[j];

                tdx_xxyy_yzz_0[j] = pa_x[j] * tdx_xyy_yzz_0[j] + 0.5 * fl1_fx * tdx_yy_yzz_0[j] + 0.5 * fl1_fx * ts_xyy_yzz_0[j];

                tdy_xxyy_yzz_0[j] = pa_x[j] * tdy_xyy_yzz_0[j] + 0.5 * fl1_fx * tdy_yy_yzz_0[j];

                tdz_xxyy_yzz_0[j] = pa_x[j] * tdz_xyy_yzz_0[j] + 0.5 * fl1_fx * tdz_yy_yzz_0[j];

                tdx_xxyy_zzz_0[j] = pa_x[j] * tdx_xyy_zzz_0[j] + 0.5 * fl1_fx * tdx_yy_zzz_0[j] + 0.5 * fl1_fx * ts_xyy_zzz_0[j];

                tdy_xxyy_zzz_0[j] = pa_x[j] * tdy_xyy_zzz_0[j] + 0.5 * fl1_fx * tdy_yy_zzz_0[j];

                tdz_xxyy_zzz_0[j] = pa_x[j] * tdz_xyy_zzz_0[j] + 0.5 * fl1_fx * tdz_yy_zzz_0[j];

                tdx_xxyz_xxx_0[j] = pa_x[j] * tdx_xyz_xxx_0[j] + 0.5 * fl1_fx * tdx_yz_xxx_0[j] + 1.5 * fl1_fx * tdx_xyz_xx_0[j] + 0.5 * fl1_fx * ts_xyz_xxx_0[j];

                tdy_xxyz_xxx_0[j] = pa_x[j] * tdy_xyz_xxx_0[j] + 0.5 * fl1_fx * tdy_yz_xxx_0[j] + 1.5 * fl1_fx * tdy_xyz_xx_0[j];

                tdz_xxyz_xxx_0[j] = pa_x[j] * tdz_xyz_xxx_0[j] + 0.5 * fl1_fx * tdz_yz_xxx_0[j] + 1.5 * fl1_fx * tdz_xyz_xx_0[j];

                tdx_xxyz_xxy_0[j] = pa_x[j] * tdx_xyz_xxy_0[j] + 0.5 * fl1_fx * tdx_yz_xxy_0[j] + fl1_fx * tdx_xyz_xy_0[j] + 0.5 * fl1_fx * ts_xyz_xxy_0[j];

                tdy_xxyz_xxy_0[j] = pa_x[j] * tdy_xyz_xxy_0[j] + 0.5 * fl1_fx * tdy_yz_xxy_0[j] + fl1_fx * tdy_xyz_xy_0[j];

                tdz_xxyz_xxy_0[j] = pa_x[j] * tdz_xyz_xxy_0[j] + 0.5 * fl1_fx * tdz_yz_xxy_0[j] + fl1_fx * tdz_xyz_xy_0[j];

                tdx_xxyz_xxz_0[j] = pa_x[j] * tdx_xyz_xxz_0[j] + 0.5 * fl1_fx * tdx_yz_xxz_0[j] + fl1_fx * tdx_xyz_xz_0[j] + 0.5 * fl1_fx * ts_xyz_xxz_0[j];

                tdy_xxyz_xxz_0[j] = pa_x[j] * tdy_xyz_xxz_0[j] + 0.5 * fl1_fx * tdy_yz_xxz_0[j] + fl1_fx * tdy_xyz_xz_0[j];

                tdz_xxyz_xxz_0[j] = pa_x[j] * tdz_xyz_xxz_0[j] + 0.5 * fl1_fx * tdz_yz_xxz_0[j] + fl1_fx * tdz_xyz_xz_0[j];

                tdx_xxyz_xyy_0[j] = pa_x[j] * tdx_xyz_xyy_0[j] + 0.5 * fl1_fx * tdx_yz_xyy_0[j] + 0.5 * fl1_fx * tdx_xyz_yy_0[j] + 0.5 * fl1_fx * ts_xyz_xyy_0[j];

                tdy_xxyz_xyy_0[j] = pa_x[j] * tdy_xyz_xyy_0[j] + 0.5 * fl1_fx * tdy_yz_xyy_0[j] + 0.5 * fl1_fx * tdy_xyz_yy_0[j];

                tdz_xxyz_xyy_0[j] = pa_x[j] * tdz_xyz_xyy_0[j] + 0.5 * fl1_fx * tdz_yz_xyy_0[j] + 0.5 * fl1_fx * tdz_xyz_yy_0[j];

                tdx_xxyz_xyz_0[j] = pa_x[j] * tdx_xyz_xyz_0[j] + 0.5 * fl1_fx * tdx_yz_xyz_0[j] + 0.5 * fl1_fx * tdx_xyz_yz_0[j] + 0.5 * fl1_fx * ts_xyz_xyz_0[j];

                tdy_xxyz_xyz_0[j] = pa_x[j] * tdy_xyz_xyz_0[j] + 0.5 * fl1_fx * tdy_yz_xyz_0[j] + 0.5 * fl1_fx * tdy_xyz_yz_0[j];

                tdz_xxyz_xyz_0[j] = pa_x[j] * tdz_xyz_xyz_0[j] + 0.5 * fl1_fx * tdz_yz_xyz_0[j] + 0.5 * fl1_fx * tdz_xyz_yz_0[j];

                tdx_xxyz_xzz_0[j] = pa_x[j] * tdx_xyz_xzz_0[j] + 0.5 * fl1_fx * tdx_yz_xzz_0[j] + 0.5 * fl1_fx * tdx_xyz_zz_0[j] + 0.5 * fl1_fx * ts_xyz_xzz_0[j];

                tdy_xxyz_xzz_0[j] = pa_x[j] * tdy_xyz_xzz_0[j] + 0.5 * fl1_fx * tdy_yz_xzz_0[j] + 0.5 * fl1_fx * tdy_xyz_zz_0[j];

                tdz_xxyz_xzz_0[j] = pa_x[j] * tdz_xyz_xzz_0[j] + 0.5 * fl1_fx * tdz_yz_xzz_0[j] + 0.5 * fl1_fx * tdz_xyz_zz_0[j];

                tdx_xxyz_yyy_0[j] = pa_x[j] * tdx_xyz_yyy_0[j] + 0.5 * fl1_fx * tdx_yz_yyy_0[j] + 0.5 * fl1_fx * ts_xyz_yyy_0[j];

                tdy_xxyz_yyy_0[j] = pa_x[j] * tdy_xyz_yyy_0[j] + 0.5 * fl1_fx * tdy_yz_yyy_0[j];

                tdz_xxyz_yyy_0[j] = pa_x[j] * tdz_xyz_yyy_0[j] + 0.5 * fl1_fx * tdz_yz_yyy_0[j];

                tdx_xxyz_yyz_0[j] = pa_x[j] * tdx_xyz_yyz_0[j] + 0.5 * fl1_fx * tdx_yz_yyz_0[j] + 0.5 * fl1_fx * ts_xyz_yyz_0[j];

                tdy_xxyz_yyz_0[j] = pa_x[j] * tdy_xyz_yyz_0[j] + 0.5 * fl1_fx * tdy_yz_yyz_0[j];

                tdz_xxyz_yyz_0[j] = pa_x[j] * tdz_xyz_yyz_0[j] + 0.5 * fl1_fx * tdz_yz_yyz_0[j];

                tdx_xxyz_yzz_0[j] = pa_x[j] * tdx_xyz_yzz_0[j] + 0.5 * fl1_fx * tdx_yz_yzz_0[j] + 0.5 * fl1_fx * ts_xyz_yzz_0[j];

                tdy_xxyz_yzz_0[j] = pa_x[j] * tdy_xyz_yzz_0[j] + 0.5 * fl1_fx * tdy_yz_yzz_0[j];

                tdz_xxyz_yzz_0[j] = pa_x[j] * tdz_xyz_yzz_0[j] + 0.5 * fl1_fx * tdz_yz_yzz_0[j];

                tdx_xxyz_zzz_0[j] = pa_x[j] * tdx_xyz_zzz_0[j] + 0.5 * fl1_fx * tdx_yz_zzz_0[j] + 0.5 * fl1_fx * ts_xyz_zzz_0[j];

                tdy_xxyz_zzz_0[j] = pa_x[j] * tdy_xyz_zzz_0[j] + 0.5 * fl1_fx * tdy_yz_zzz_0[j];

                tdz_xxyz_zzz_0[j] = pa_x[j] * tdz_xyz_zzz_0[j] + 0.5 * fl1_fx * tdz_yz_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGF_150_200(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (150,200)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_xzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 50); 

            auto tdy_xzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 50); 

            auto tdz_xzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 50); 

            auto tdx_xzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 51); 

            auto tdy_xzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 51); 

            auto tdz_xzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 51); 

            auto tdx_xzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 52); 

            auto tdy_xzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 52); 

            auto tdz_xzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 52); 

            auto tdx_xzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 53); 

            auto tdy_xzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 53); 

            auto tdz_xzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 53); 

            auto tdx_xzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 54); 

            auto tdy_xzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 54); 

            auto tdz_xzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 54); 

            auto tdx_xzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 55); 

            auto tdy_xzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 55); 

            auto tdz_xzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 55); 

            auto tdx_xzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 56); 

            auto tdy_xzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 56); 

            auto tdz_xzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 56); 

            auto tdx_xzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 57); 

            auto tdy_xzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 57); 

            auto tdz_xzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 57); 

            auto tdx_xzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 58); 

            auto tdy_xzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 58); 

            auto tdz_xzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 58); 

            auto tdx_xzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 59); 

            auto tdy_xzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 59); 

            auto tdz_xzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 59); 

            auto tdx_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 60); 

            auto tdy_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 60); 

            auto tdz_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 60); 

            auto tdx_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 61); 

            auto tdy_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 61); 

            auto tdz_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 61); 

            auto tdx_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 62); 

            auto tdy_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 62); 

            auto tdz_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 62); 

            auto tdx_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 63); 

            auto tdy_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 63); 

            auto tdz_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 63); 

            auto tdx_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 64); 

            auto tdy_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 64); 

            auto tdz_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 64); 

            auto tdx_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 65); 

            auto tdy_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 65); 

            auto tdz_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 65); 

            auto tdx_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 66); 

            auto tdy_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 66); 

            auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50); 

            auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51); 

            auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52); 

            auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 53); 

            auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53); 

            auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53); 

            auto tdx_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 54); 

            auto tdy_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdx_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 55); 

            auto tdy_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 55); 

            auto tdz_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 55); 

            auto tdx_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 56); 

            auto tdy_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 56); 

            auto tdz_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 56); 

            auto tdx_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 57); 

            auto tdy_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 57); 

            auto tdz_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 57); 

            auto tdx_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 58); 

            auto tdy_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 58); 

            auto tdz_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 58); 

            auto tdx_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 59); 

            auto tdy_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 59); 

            auto tdz_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 59); 

            auto tdx_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 30); 

            auto tdy_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 30); 

            auto tdz_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 30); 

            auto tdx_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 31); 

            auto tdy_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 31); 

            auto tdz_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 31); 

            auto tdx_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 32); 

            auto tdy_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 32); 

            auto tdz_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 32); 

            auto tdx_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 33); 

            auto tdy_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 33); 

            auto tdz_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 33); 

            auto tdx_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 34); 

            auto tdy_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 34); 

            auto tdz_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 34); 

            auto tdx_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 35); 

            auto tdy_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 35); 

            auto tdz_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 35); 

            auto tdx_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 36); 

            auto tdy_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 37); 

            auto tdy_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 38); 

            auto tdy_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 39); 

            auto tdy_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 40); 

            auto tdy_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 41); 

            auto tdy_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto ts_xzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 50); 

            auto ts_xzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 51); 

            auto ts_xzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 52); 

            auto ts_xzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 53); 

            auto ts_xzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 54); 

            auto ts_xzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 55); 

            auto ts_xzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 56); 

            auto ts_xzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 57); 

            auto ts_xzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 58); 

            auto ts_xzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 59); 

            auto ts_yyy_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 60); 

            auto ts_yyy_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 61); 

            auto ts_yyy_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 62); 

            auto ts_yyy_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 63); 

            auto ts_yyy_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 64); 

            auto ts_yyy_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 65); 

            auto ts_yyy_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 66); 

            // set up pointers to integrals

            auto tdx_xxzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 50); 

            auto tdy_xxzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 50); 

            auto tdz_xxzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 50); 

            auto tdx_xxzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 51); 

            auto tdy_xxzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 51); 

            auto tdz_xxzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 51); 

            auto tdx_xxzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 52); 

            auto tdy_xxzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 52); 

            auto tdz_xxzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 52); 

            auto tdx_xxzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 53); 

            auto tdy_xxzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 53); 

            auto tdz_xxzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 53); 

            auto tdx_xxzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 54); 

            auto tdy_xxzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 54); 

            auto tdz_xxzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 54); 

            auto tdx_xxzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 55); 

            auto tdy_xxzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 55); 

            auto tdz_xxzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 55); 

            auto tdx_xxzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 56); 

            auto tdy_xxzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 56); 

            auto tdz_xxzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 56); 

            auto tdx_xxzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 57); 

            auto tdy_xxzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 57); 

            auto tdz_xxzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 57); 

            auto tdx_xxzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 58); 

            auto tdy_xxzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 58); 

            auto tdz_xxzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 58); 

            auto tdx_xxzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 59); 

            auto tdy_xxzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 59); 

            auto tdz_xxzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 59); 

            auto tdx_xyyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 60); 

            auto tdy_xyyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 60); 

            auto tdz_xyyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 60); 

            auto tdx_xyyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 61); 

            auto tdy_xyyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 61); 

            auto tdz_xyyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 61); 

            auto tdx_xyyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 62); 

            auto tdy_xyyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 62); 

            auto tdz_xyyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 62); 

            auto tdx_xyyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 63); 

            auto tdy_xyyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 63); 

            auto tdz_xyyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 63); 

            auto tdx_xyyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 64); 

            auto tdy_xyyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 64); 

            auto tdz_xyyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 64); 

            auto tdx_xyyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 65); 

            auto tdy_xyyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 65); 

            auto tdz_xyyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 65); 

            auto tdx_xyyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 66); 

            auto tdy_xyyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 66); 

            // Batch of Integrals (150,200)

            #pragma omp simd aligned(fx, pa_x, tdx_xxzz_xxx_0, tdx_xxzz_xxy_0, tdx_xxzz_xxz_0, \
                                     tdx_xxzz_xyy_0, tdx_xxzz_xyz_0, tdx_xxzz_xzz_0, tdx_xxzz_yyy_0, tdx_xxzz_yyz_0, \
                                     tdx_xxzz_yzz_0, tdx_xxzz_zzz_0, tdx_xyyy_xxx_0, tdx_xyyy_xxy_0, tdx_xyyy_xxz_0, \
                                     tdx_xyyy_xyy_0, tdx_xyyy_xyz_0, tdx_xyyy_xzz_0, tdx_xyyy_yyy_0, tdx_xzz_xx_0, \
                                     tdx_xzz_xxx_0, tdx_xzz_xxy_0, tdx_xzz_xxz_0, tdx_xzz_xy_0, tdx_xzz_xyy_0, \
                                     tdx_xzz_xyz_0, tdx_xzz_xz_0, tdx_xzz_xzz_0, tdx_xzz_yy_0, tdx_xzz_yyy_0, \
                                     tdx_xzz_yyz_0, tdx_xzz_yz_0, tdx_xzz_yzz_0, tdx_xzz_zz_0, tdx_xzz_zzz_0, \
                                     tdx_yyy_xx_0, tdx_yyy_xxx_0, tdx_yyy_xxy_0, tdx_yyy_xxz_0, tdx_yyy_xy_0, \
                                     tdx_yyy_xyy_0, tdx_yyy_xyz_0, tdx_yyy_xz_0, tdx_yyy_xzz_0, tdx_yyy_yy_0, \
                                     tdx_yyy_yyy_0, tdx_yyy_yz_0, tdx_yyy_zz_0, tdx_zz_xxx_0, tdx_zz_xxy_0, tdx_zz_xxz_0, \
                                     tdx_zz_xyy_0, tdx_zz_xyz_0, tdx_zz_xzz_0, tdx_zz_yyy_0, tdx_zz_yyz_0, tdx_zz_yzz_0, \
                                     tdx_zz_zzz_0, tdy_xxzz_xxx_0, tdy_xxzz_xxy_0, tdy_xxzz_xxz_0, tdy_xxzz_xyy_0, \
                                     tdy_xxzz_xyz_0, tdy_xxzz_xzz_0, tdy_xxzz_yyy_0, tdy_xxzz_yyz_0, tdy_xxzz_yzz_0, \
                                     tdy_xxzz_zzz_0, tdy_xyyy_xxx_0, tdy_xyyy_xxy_0, tdy_xyyy_xxz_0, tdy_xyyy_xyy_0, \
                                     tdy_xyyy_xyz_0, tdy_xyyy_xzz_0, tdy_xyyy_yyy_0, tdy_xzz_xx_0, tdy_xzz_xxx_0, \
                                     tdy_xzz_xxy_0, tdy_xzz_xxz_0, tdy_xzz_xy_0, tdy_xzz_xyy_0, tdy_xzz_xyz_0, \
                                     tdy_xzz_xz_0, tdy_xzz_xzz_0, tdy_xzz_yy_0, tdy_xzz_yyy_0, tdy_xzz_yyz_0, \
                                     tdy_xzz_yz_0, tdy_xzz_yzz_0, tdy_xzz_zz_0, tdy_xzz_zzz_0, tdy_yyy_xx_0, \
                                     tdy_yyy_xxx_0, tdy_yyy_xxy_0, tdy_yyy_xxz_0, tdy_yyy_xy_0, tdy_yyy_xyy_0, \
                                     tdy_yyy_xyz_0, tdy_yyy_xz_0, tdy_yyy_xzz_0, tdy_yyy_yy_0, tdy_yyy_yyy_0, \
                                     tdy_yyy_yz_0, tdy_yyy_zz_0, tdy_zz_xxx_0, tdy_zz_xxy_0, tdy_zz_xxz_0, tdy_zz_xyy_0, \
                                     tdy_zz_xyz_0, tdy_zz_xzz_0, tdy_zz_yyy_0, tdy_zz_yyz_0, tdy_zz_yzz_0, tdy_zz_zzz_0, \
                                     tdz_xxzz_xxx_0, tdz_xxzz_xxy_0, tdz_xxzz_xxz_0, tdz_xxzz_xyy_0, tdz_xxzz_xyz_0, \
                                     tdz_xxzz_xzz_0, tdz_xxzz_yyy_0, tdz_xxzz_yyz_0, tdz_xxzz_yzz_0, tdz_xxzz_zzz_0, \
                                     tdz_xyyy_xxx_0, tdz_xyyy_xxy_0, tdz_xyyy_xxz_0, tdz_xyyy_xyy_0, tdz_xyyy_xyz_0, \
                                     tdz_xyyy_xzz_0, tdz_xzz_xx_0, tdz_xzz_xxx_0, tdz_xzz_xxy_0, tdz_xzz_xxz_0, \
                                     tdz_xzz_xy_0, tdz_xzz_xyy_0, tdz_xzz_xyz_0, tdz_xzz_xz_0, tdz_xzz_xzz_0, \
                                     tdz_xzz_yy_0, tdz_xzz_yyy_0, tdz_xzz_yyz_0, tdz_xzz_yz_0, tdz_xzz_yzz_0, \
                                     tdz_xzz_zz_0, tdz_xzz_zzz_0, tdz_yyy_xx_0, tdz_yyy_xxx_0, tdz_yyy_xxy_0, \
                                     tdz_yyy_xxz_0, tdz_yyy_xy_0, tdz_yyy_xyy_0, tdz_yyy_xyz_0, tdz_yyy_xz_0, \
                                     tdz_yyy_xzz_0, tdz_yyy_yy_0, tdz_yyy_yz_0, tdz_yyy_zz_0, tdz_zz_xxx_0, tdz_zz_xxy_0, \
                                     tdz_zz_xxz_0, tdz_zz_xyy_0, tdz_zz_xyz_0, tdz_zz_xzz_0, tdz_zz_yyy_0, tdz_zz_yyz_0, \
                                     tdz_zz_yzz_0, tdz_zz_zzz_0, ts_xzz_xxx_0, ts_xzz_xxy_0, ts_xzz_xxz_0, ts_xzz_xyy_0, \
                                     ts_xzz_xyz_0, ts_xzz_xzz_0, ts_xzz_yyy_0, ts_xzz_yyz_0, ts_xzz_yzz_0, ts_xzz_zzz_0, \
                                     ts_yyy_xxx_0, ts_yyy_xxy_0, ts_yyy_xxz_0, ts_yyy_xyy_0, ts_yyy_xyz_0, ts_yyy_xzz_0, \
                                     ts_yyy_yyy_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxzz_xxx_0[j] = pa_x[j] * tdx_xzz_xxx_0[j] + 0.5 * fl1_fx * tdx_zz_xxx_0[j] + 1.5 * fl1_fx * tdx_xzz_xx_0[j] + 0.5 * fl1_fx * ts_xzz_xxx_0[j];

                tdy_xxzz_xxx_0[j] = pa_x[j] * tdy_xzz_xxx_0[j] + 0.5 * fl1_fx * tdy_zz_xxx_0[j] + 1.5 * fl1_fx * tdy_xzz_xx_0[j];

                tdz_xxzz_xxx_0[j] = pa_x[j] * tdz_xzz_xxx_0[j] + 0.5 * fl1_fx * tdz_zz_xxx_0[j] + 1.5 * fl1_fx * tdz_xzz_xx_0[j];

                tdx_xxzz_xxy_0[j] = pa_x[j] * tdx_xzz_xxy_0[j] + 0.5 * fl1_fx * tdx_zz_xxy_0[j] + fl1_fx * tdx_xzz_xy_0[j] + 0.5 * fl1_fx * ts_xzz_xxy_0[j];

                tdy_xxzz_xxy_0[j] = pa_x[j] * tdy_xzz_xxy_0[j] + 0.5 * fl1_fx * tdy_zz_xxy_0[j] + fl1_fx * tdy_xzz_xy_0[j];

                tdz_xxzz_xxy_0[j] = pa_x[j] * tdz_xzz_xxy_0[j] + 0.5 * fl1_fx * tdz_zz_xxy_0[j] + fl1_fx * tdz_xzz_xy_0[j];

                tdx_xxzz_xxz_0[j] = pa_x[j] * tdx_xzz_xxz_0[j] + 0.5 * fl1_fx * tdx_zz_xxz_0[j] + fl1_fx * tdx_xzz_xz_0[j] + 0.5 * fl1_fx * ts_xzz_xxz_0[j];

                tdy_xxzz_xxz_0[j] = pa_x[j] * tdy_xzz_xxz_0[j] + 0.5 * fl1_fx * tdy_zz_xxz_0[j] + fl1_fx * tdy_xzz_xz_0[j];

                tdz_xxzz_xxz_0[j] = pa_x[j] * tdz_xzz_xxz_0[j] + 0.5 * fl1_fx * tdz_zz_xxz_0[j] + fl1_fx * tdz_xzz_xz_0[j];

                tdx_xxzz_xyy_0[j] = pa_x[j] * tdx_xzz_xyy_0[j] + 0.5 * fl1_fx * tdx_zz_xyy_0[j] + 0.5 * fl1_fx * tdx_xzz_yy_0[j] + 0.5 * fl1_fx * ts_xzz_xyy_0[j];

                tdy_xxzz_xyy_0[j] = pa_x[j] * tdy_xzz_xyy_0[j] + 0.5 * fl1_fx * tdy_zz_xyy_0[j] + 0.5 * fl1_fx * tdy_xzz_yy_0[j];

                tdz_xxzz_xyy_0[j] = pa_x[j] * tdz_xzz_xyy_0[j] + 0.5 * fl1_fx * tdz_zz_xyy_0[j] + 0.5 * fl1_fx * tdz_xzz_yy_0[j];

                tdx_xxzz_xyz_0[j] = pa_x[j] * tdx_xzz_xyz_0[j] + 0.5 * fl1_fx * tdx_zz_xyz_0[j] + 0.5 * fl1_fx * tdx_xzz_yz_0[j] + 0.5 * fl1_fx * ts_xzz_xyz_0[j];

                tdy_xxzz_xyz_0[j] = pa_x[j] * tdy_xzz_xyz_0[j] + 0.5 * fl1_fx * tdy_zz_xyz_0[j] + 0.5 * fl1_fx * tdy_xzz_yz_0[j];

                tdz_xxzz_xyz_0[j] = pa_x[j] * tdz_xzz_xyz_0[j] + 0.5 * fl1_fx * tdz_zz_xyz_0[j] + 0.5 * fl1_fx * tdz_xzz_yz_0[j];

                tdx_xxzz_xzz_0[j] = pa_x[j] * tdx_xzz_xzz_0[j] + 0.5 * fl1_fx * tdx_zz_xzz_0[j] + 0.5 * fl1_fx * tdx_xzz_zz_0[j] + 0.5 * fl1_fx * ts_xzz_xzz_0[j];

                tdy_xxzz_xzz_0[j] = pa_x[j] * tdy_xzz_xzz_0[j] + 0.5 * fl1_fx * tdy_zz_xzz_0[j] + 0.5 * fl1_fx * tdy_xzz_zz_0[j];

                tdz_xxzz_xzz_0[j] = pa_x[j] * tdz_xzz_xzz_0[j] + 0.5 * fl1_fx * tdz_zz_xzz_0[j] + 0.5 * fl1_fx * tdz_xzz_zz_0[j];

                tdx_xxzz_yyy_0[j] = pa_x[j] * tdx_xzz_yyy_0[j] + 0.5 * fl1_fx * tdx_zz_yyy_0[j] + 0.5 * fl1_fx * ts_xzz_yyy_0[j];

                tdy_xxzz_yyy_0[j] = pa_x[j] * tdy_xzz_yyy_0[j] + 0.5 * fl1_fx * tdy_zz_yyy_0[j];

                tdz_xxzz_yyy_0[j] = pa_x[j] * tdz_xzz_yyy_0[j] + 0.5 * fl1_fx * tdz_zz_yyy_0[j];

                tdx_xxzz_yyz_0[j] = pa_x[j] * tdx_xzz_yyz_0[j] + 0.5 * fl1_fx * tdx_zz_yyz_0[j] + 0.5 * fl1_fx * ts_xzz_yyz_0[j];

                tdy_xxzz_yyz_0[j] = pa_x[j] * tdy_xzz_yyz_0[j] + 0.5 * fl1_fx * tdy_zz_yyz_0[j];

                tdz_xxzz_yyz_0[j] = pa_x[j] * tdz_xzz_yyz_0[j] + 0.5 * fl1_fx * tdz_zz_yyz_0[j];

                tdx_xxzz_yzz_0[j] = pa_x[j] * tdx_xzz_yzz_0[j] + 0.5 * fl1_fx * tdx_zz_yzz_0[j] + 0.5 * fl1_fx * ts_xzz_yzz_0[j];

                tdy_xxzz_yzz_0[j] = pa_x[j] * tdy_xzz_yzz_0[j] + 0.5 * fl1_fx * tdy_zz_yzz_0[j];

                tdz_xxzz_yzz_0[j] = pa_x[j] * tdz_xzz_yzz_0[j] + 0.5 * fl1_fx * tdz_zz_yzz_0[j];

                tdx_xxzz_zzz_0[j] = pa_x[j] * tdx_xzz_zzz_0[j] + 0.5 * fl1_fx * tdx_zz_zzz_0[j] + 0.5 * fl1_fx * ts_xzz_zzz_0[j];

                tdy_xxzz_zzz_0[j] = pa_x[j] * tdy_xzz_zzz_0[j] + 0.5 * fl1_fx * tdy_zz_zzz_0[j];

                tdz_xxzz_zzz_0[j] = pa_x[j] * tdz_xzz_zzz_0[j] + 0.5 * fl1_fx * tdz_zz_zzz_0[j];

                tdx_xyyy_xxx_0[j] = pa_x[j] * tdx_yyy_xxx_0[j] + 1.5 * fl1_fx * tdx_yyy_xx_0[j] + 0.5 * fl1_fx * ts_yyy_xxx_0[j];

                tdy_xyyy_xxx_0[j] = pa_x[j] * tdy_yyy_xxx_0[j] + 1.5 * fl1_fx * tdy_yyy_xx_0[j];

                tdz_xyyy_xxx_0[j] = pa_x[j] * tdz_yyy_xxx_0[j] + 1.5 * fl1_fx * tdz_yyy_xx_0[j];

                tdx_xyyy_xxy_0[j] = pa_x[j] * tdx_yyy_xxy_0[j] + fl1_fx * tdx_yyy_xy_0[j] + 0.5 * fl1_fx * ts_yyy_xxy_0[j];

                tdy_xyyy_xxy_0[j] = pa_x[j] * tdy_yyy_xxy_0[j] + fl1_fx * tdy_yyy_xy_0[j];

                tdz_xyyy_xxy_0[j] = pa_x[j] * tdz_yyy_xxy_0[j] + fl1_fx * tdz_yyy_xy_0[j];

                tdx_xyyy_xxz_0[j] = pa_x[j] * tdx_yyy_xxz_0[j] + fl1_fx * tdx_yyy_xz_0[j] + 0.5 * fl1_fx * ts_yyy_xxz_0[j];

                tdy_xyyy_xxz_0[j] = pa_x[j] * tdy_yyy_xxz_0[j] + fl1_fx * tdy_yyy_xz_0[j];

                tdz_xyyy_xxz_0[j] = pa_x[j] * tdz_yyy_xxz_0[j] + fl1_fx * tdz_yyy_xz_0[j];

                tdx_xyyy_xyy_0[j] = pa_x[j] * tdx_yyy_xyy_0[j] + 0.5 * fl1_fx * tdx_yyy_yy_0[j] + 0.5 * fl1_fx * ts_yyy_xyy_0[j];

                tdy_xyyy_xyy_0[j] = pa_x[j] * tdy_yyy_xyy_0[j] + 0.5 * fl1_fx * tdy_yyy_yy_0[j];

                tdz_xyyy_xyy_0[j] = pa_x[j] * tdz_yyy_xyy_0[j] + 0.5 * fl1_fx * tdz_yyy_yy_0[j];

                tdx_xyyy_xyz_0[j] = pa_x[j] * tdx_yyy_xyz_0[j] + 0.5 * fl1_fx * tdx_yyy_yz_0[j] + 0.5 * fl1_fx * ts_yyy_xyz_0[j];

                tdy_xyyy_xyz_0[j] = pa_x[j] * tdy_yyy_xyz_0[j] + 0.5 * fl1_fx * tdy_yyy_yz_0[j];

                tdz_xyyy_xyz_0[j] = pa_x[j] * tdz_yyy_xyz_0[j] + 0.5 * fl1_fx * tdz_yyy_yz_0[j];

                tdx_xyyy_xzz_0[j] = pa_x[j] * tdx_yyy_xzz_0[j] + 0.5 * fl1_fx * tdx_yyy_zz_0[j] + 0.5 * fl1_fx * ts_yyy_xzz_0[j];

                tdy_xyyy_xzz_0[j] = pa_x[j] * tdy_yyy_xzz_0[j] + 0.5 * fl1_fx * tdy_yyy_zz_0[j];

                tdz_xyyy_xzz_0[j] = pa_x[j] * tdz_yyy_xzz_0[j] + 0.5 * fl1_fx * tdz_yyy_zz_0[j];

                tdx_xyyy_yyy_0[j] = pa_x[j] * tdx_yyy_yyy_0[j] + 0.5 * fl1_fx * ts_yyy_yyy_0[j];

                tdy_xyyy_yyy_0[j] = pa_x[j] * tdy_yyy_yyy_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGF_200_250(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (200,250)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdz_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 66); 

            auto tdx_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 67); 

            auto tdy_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 67); 

            auto tdz_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 67); 

            auto tdx_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 68); 

            auto tdy_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 68); 

            auto tdz_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 68); 

            auto tdx_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 69); 

            auto tdy_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 69); 

            auto tdz_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 69); 

            auto tdx_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 70); 

            auto tdy_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 70); 

            auto tdz_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 70); 

            auto tdx_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 71); 

            auto tdy_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 71); 

            auto tdz_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 71); 

            auto tdx_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 72); 

            auto tdy_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 72); 

            auto tdz_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 72); 

            auto tdx_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 73); 

            auto tdy_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 73); 

            auto tdz_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 73); 

            auto tdx_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 74); 

            auto tdy_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 74); 

            auto tdz_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 74); 

            auto tdx_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 75); 

            auto tdy_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 75); 

            auto tdz_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 75); 

            auto tdx_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 76); 

            auto tdy_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 76); 

            auto tdz_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 76); 

            auto tdx_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 77); 

            auto tdy_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 77); 

            auto tdz_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 77); 

            auto tdx_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 78); 

            auto tdy_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 78); 

            auto tdz_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 78); 

            auto tdx_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 79); 

            auto tdy_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 79); 

            auto tdz_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 79); 

            auto tdx_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 80); 

            auto tdy_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 80); 

            auto tdz_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 80); 

            auto tdx_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 81); 

            auto tdy_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 81); 

            auto tdz_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 81); 

            auto tdx_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 82); 

            auto tdy_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 82); 

            auto tdz_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 82); 

            auto tdx_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 83); 

            auto tdx_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 42); 

            auto tdy_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 43); 

            auto tdy_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 44); 

            auto tdy_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            auto tdx_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 45); 

            auto tdy_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tdz_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 46); 

            auto tdy_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdz_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 47); 

            auto tdy_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 48); 

            auto tdy_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 49); 

            auto tdy_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdx_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 50); 

            auto tdy_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 51); 

            auto ts_yyy_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 67); 

            auto ts_yyy_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 68); 

            auto ts_yyy_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 69); 

            auto ts_yyz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 70); 

            auto ts_yyz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 71); 

            auto ts_yyz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 72); 

            auto ts_yyz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 73); 

            auto ts_yyz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 74); 

            auto ts_yyz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 75); 

            auto ts_yyz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 76); 

            auto ts_yyz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 77); 

            auto ts_yyz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 78); 

            auto ts_yyz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 79); 

            auto ts_yzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 80); 

            auto ts_yzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 81); 

            auto ts_yzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 82); 

            auto ts_yzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 83); 

            // set up pointers to integrals

            auto tdz_xyyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 66); 

            auto tdx_xyyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 67); 

            auto tdy_xyyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 67); 

            auto tdz_xyyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 67); 

            auto tdx_xyyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 68); 

            auto tdy_xyyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 68); 

            auto tdz_xyyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 68); 

            auto tdx_xyyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 69); 

            auto tdy_xyyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 69); 

            auto tdz_xyyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 69); 

            auto tdx_xyyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 70); 

            auto tdy_xyyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 70); 

            auto tdz_xyyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 70); 

            auto tdx_xyyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 71); 

            auto tdy_xyyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 71); 

            auto tdz_xyyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 71); 

            auto tdx_xyyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 72); 

            auto tdy_xyyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 72); 

            auto tdz_xyyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 72); 

            auto tdx_xyyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 73); 

            auto tdy_xyyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 73); 

            auto tdz_xyyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 73); 

            auto tdx_xyyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 74); 

            auto tdy_xyyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 74); 

            auto tdz_xyyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 74); 

            auto tdx_xyyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 75); 

            auto tdy_xyyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 75); 

            auto tdz_xyyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 75); 

            auto tdx_xyyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 76); 

            auto tdy_xyyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 76); 

            auto tdz_xyyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 76); 

            auto tdx_xyyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 77); 

            auto tdy_xyyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 77); 

            auto tdz_xyyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 77); 

            auto tdx_xyyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 78); 

            auto tdy_xyyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 78); 

            auto tdz_xyyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 78); 

            auto tdx_xyyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 79); 

            auto tdy_xyyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 79); 

            auto tdz_xyyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 79); 

            auto tdx_xyzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 80); 

            auto tdy_xyzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 80); 

            auto tdz_xyzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 80); 

            auto tdx_xyzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 81); 

            auto tdy_xyzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 81); 

            auto tdz_xyzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 81); 

            auto tdx_xyzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 82); 

            auto tdy_xyzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 82); 

            auto tdz_xyzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 82); 

            auto tdx_xyzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 83); 

            // Batch of Integrals (200,250)

            #pragma omp simd aligned(fx, pa_x, tdx_xyyy_yyz_0, tdx_xyyy_yzz_0, tdx_xyyy_zzz_0, \
                                     tdx_xyyz_xxx_0, tdx_xyyz_xxy_0, tdx_xyyz_xxz_0, tdx_xyyz_xyy_0, tdx_xyyz_xyz_0, \
                                     tdx_xyyz_xzz_0, tdx_xyyz_yyy_0, tdx_xyyz_yyz_0, tdx_xyyz_yzz_0, tdx_xyyz_zzz_0, \
                                     tdx_xyzz_xxx_0, tdx_xyzz_xxy_0, tdx_xyzz_xxz_0, tdx_xyzz_xyy_0, tdx_yyy_yyz_0, \
                                     tdx_yyy_yzz_0, tdx_yyy_zzz_0, tdx_yyz_xx_0, tdx_yyz_xxx_0, tdx_yyz_xxy_0, \
                                     tdx_yyz_xxz_0, tdx_yyz_xy_0, tdx_yyz_xyy_0, tdx_yyz_xyz_0, tdx_yyz_xz_0, \
                                     tdx_yyz_xzz_0, tdx_yyz_yy_0, tdx_yyz_yyy_0, tdx_yyz_yyz_0, tdx_yyz_yz_0, \
                                     tdx_yyz_yzz_0, tdx_yyz_zz_0, tdx_yyz_zzz_0, tdx_yzz_xx_0, tdx_yzz_xxx_0, \
                                     tdx_yzz_xxy_0, tdx_yzz_xxz_0, tdx_yzz_xy_0, tdx_yzz_xyy_0, tdx_yzz_xz_0, \
                                     tdx_yzz_yy_0, tdy_xyyy_yyz_0, tdy_xyyy_yzz_0, tdy_xyyy_zzz_0, tdy_xyyz_xxx_0, \
                                     tdy_xyyz_xxy_0, tdy_xyyz_xxz_0, tdy_xyyz_xyy_0, tdy_xyyz_xyz_0, tdy_xyyz_xzz_0, \
                                     tdy_xyyz_yyy_0, tdy_xyyz_yyz_0, tdy_xyyz_yzz_0, tdy_xyyz_zzz_0, tdy_xyzz_xxx_0, \
                                     tdy_xyzz_xxy_0, tdy_xyzz_xxz_0, tdy_yyy_yyz_0, tdy_yyy_yzz_0, tdy_yyy_zzz_0, \
                                     tdy_yyz_xx_0, tdy_yyz_xxx_0, tdy_yyz_xxy_0, tdy_yyz_xxz_0, tdy_yyz_xy_0, \
                                     tdy_yyz_xyy_0, tdy_yyz_xyz_0, tdy_yyz_xz_0, tdy_yyz_xzz_0, tdy_yyz_yy_0, \
                                     tdy_yyz_yyy_0, tdy_yyz_yyz_0, tdy_yyz_yz_0, tdy_yyz_yzz_0, tdy_yyz_zz_0, \
                                     tdy_yyz_zzz_0, tdy_yzz_xx_0, tdy_yzz_xxx_0, tdy_yzz_xxy_0, tdy_yzz_xxz_0, \
                                     tdy_yzz_xy_0, tdy_yzz_xz_0, tdz_xyyy_yyy_0, tdz_xyyy_yyz_0, tdz_xyyy_yzz_0, \
                                     tdz_xyyy_zzz_0, tdz_xyyz_xxx_0, tdz_xyyz_xxy_0, tdz_xyyz_xxz_0, tdz_xyyz_xyy_0, \
                                     tdz_xyyz_xyz_0, tdz_xyyz_xzz_0, tdz_xyyz_yyy_0, tdz_xyyz_yyz_0, tdz_xyyz_yzz_0, \
                                     tdz_xyyz_zzz_0, tdz_xyzz_xxx_0, tdz_xyzz_xxy_0, tdz_xyzz_xxz_0, tdz_yyy_yyy_0, \
                                     tdz_yyy_yyz_0, tdz_yyy_yzz_0, tdz_yyy_zzz_0, tdz_yyz_xx_0, tdz_yyz_xxx_0, \
                                     tdz_yyz_xxy_0, tdz_yyz_xxz_0, tdz_yyz_xy_0, tdz_yyz_xyy_0, tdz_yyz_xyz_0, \
                                     tdz_yyz_xz_0, tdz_yyz_xzz_0, tdz_yyz_yy_0, tdz_yyz_yyy_0, tdz_yyz_yyz_0, \
                                     tdz_yyz_yz_0, tdz_yyz_yzz_0, tdz_yyz_zz_0, tdz_yyz_zzz_0, tdz_yzz_xx_0, \
                                     tdz_yzz_xxx_0, tdz_yzz_xxy_0, tdz_yzz_xxz_0, tdz_yzz_xy_0, tdz_yzz_xz_0, \
                                     ts_yyy_yyz_0, ts_yyy_yzz_0, ts_yyy_zzz_0, ts_yyz_xxx_0, ts_yyz_xxy_0, ts_yyz_xxz_0, \
                                     ts_yyz_xyy_0, ts_yyz_xyz_0, ts_yyz_xzz_0, ts_yyz_yyy_0, ts_yyz_yyz_0, ts_yyz_yzz_0, \
                                     ts_yyz_zzz_0, ts_yzz_xxx_0, ts_yzz_xxy_0, ts_yzz_xxz_0, ts_yzz_xyy_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdz_xyyy_yyy_0[j] = pa_x[j] * tdz_yyy_yyy_0[j];

                tdx_xyyy_yyz_0[j] = pa_x[j] * tdx_yyy_yyz_0[j] + 0.5 * fl1_fx * ts_yyy_yyz_0[j];

                tdy_xyyy_yyz_0[j] = pa_x[j] * tdy_yyy_yyz_0[j];

                tdz_xyyy_yyz_0[j] = pa_x[j] * tdz_yyy_yyz_0[j];

                tdx_xyyy_yzz_0[j] = pa_x[j] * tdx_yyy_yzz_0[j] + 0.5 * fl1_fx * ts_yyy_yzz_0[j];

                tdy_xyyy_yzz_0[j] = pa_x[j] * tdy_yyy_yzz_0[j];

                tdz_xyyy_yzz_0[j] = pa_x[j] * tdz_yyy_yzz_0[j];

                tdx_xyyy_zzz_0[j] = pa_x[j] * tdx_yyy_zzz_0[j] + 0.5 * fl1_fx * ts_yyy_zzz_0[j];

                tdy_xyyy_zzz_0[j] = pa_x[j] * tdy_yyy_zzz_0[j];

                tdz_xyyy_zzz_0[j] = pa_x[j] * tdz_yyy_zzz_0[j];

                tdx_xyyz_xxx_0[j] = pa_x[j] * tdx_yyz_xxx_0[j] + 1.5 * fl1_fx * tdx_yyz_xx_0[j] + 0.5 * fl1_fx * ts_yyz_xxx_0[j];

                tdy_xyyz_xxx_0[j] = pa_x[j] * tdy_yyz_xxx_0[j] + 1.5 * fl1_fx * tdy_yyz_xx_0[j];

                tdz_xyyz_xxx_0[j] = pa_x[j] * tdz_yyz_xxx_0[j] + 1.5 * fl1_fx * tdz_yyz_xx_0[j];

                tdx_xyyz_xxy_0[j] = pa_x[j] * tdx_yyz_xxy_0[j] + fl1_fx * tdx_yyz_xy_0[j] + 0.5 * fl1_fx * ts_yyz_xxy_0[j];

                tdy_xyyz_xxy_0[j] = pa_x[j] * tdy_yyz_xxy_0[j] + fl1_fx * tdy_yyz_xy_0[j];

                tdz_xyyz_xxy_0[j] = pa_x[j] * tdz_yyz_xxy_0[j] + fl1_fx * tdz_yyz_xy_0[j];

                tdx_xyyz_xxz_0[j] = pa_x[j] * tdx_yyz_xxz_0[j] + fl1_fx * tdx_yyz_xz_0[j] + 0.5 * fl1_fx * ts_yyz_xxz_0[j];

                tdy_xyyz_xxz_0[j] = pa_x[j] * tdy_yyz_xxz_0[j] + fl1_fx * tdy_yyz_xz_0[j];

                tdz_xyyz_xxz_0[j] = pa_x[j] * tdz_yyz_xxz_0[j] + fl1_fx * tdz_yyz_xz_0[j];

                tdx_xyyz_xyy_0[j] = pa_x[j] * tdx_yyz_xyy_0[j] + 0.5 * fl1_fx * tdx_yyz_yy_0[j] + 0.5 * fl1_fx * ts_yyz_xyy_0[j];

                tdy_xyyz_xyy_0[j] = pa_x[j] * tdy_yyz_xyy_0[j] + 0.5 * fl1_fx * tdy_yyz_yy_0[j];

                tdz_xyyz_xyy_0[j] = pa_x[j] * tdz_yyz_xyy_0[j] + 0.5 * fl1_fx * tdz_yyz_yy_0[j];

                tdx_xyyz_xyz_0[j] = pa_x[j] * tdx_yyz_xyz_0[j] + 0.5 * fl1_fx * tdx_yyz_yz_0[j] + 0.5 * fl1_fx * ts_yyz_xyz_0[j];

                tdy_xyyz_xyz_0[j] = pa_x[j] * tdy_yyz_xyz_0[j] + 0.5 * fl1_fx * tdy_yyz_yz_0[j];

                tdz_xyyz_xyz_0[j] = pa_x[j] * tdz_yyz_xyz_0[j] + 0.5 * fl1_fx * tdz_yyz_yz_0[j];

                tdx_xyyz_xzz_0[j] = pa_x[j] * tdx_yyz_xzz_0[j] + 0.5 * fl1_fx * tdx_yyz_zz_0[j] + 0.5 * fl1_fx * ts_yyz_xzz_0[j];

                tdy_xyyz_xzz_0[j] = pa_x[j] * tdy_yyz_xzz_0[j] + 0.5 * fl1_fx * tdy_yyz_zz_0[j];

                tdz_xyyz_xzz_0[j] = pa_x[j] * tdz_yyz_xzz_0[j] + 0.5 * fl1_fx * tdz_yyz_zz_0[j];

                tdx_xyyz_yyy_0[j] = pa_x[j] * tdx_yyz_yyy_0[j] + 0.5 * fl1_fx * ts_yyz_yyy_0[j];

                tdy_xyyz_yyy_0[j] = pa_x[j] * tdy_yyz_yyy_0[j];

                tdz_xyyz_yyy_0[j] = pa_x[j] * tdz_yyz_yyy_0[j];

                tdx_xyyz_yyz_0[j] = pa_x[j] * tdx_yyz_yyz_0[j] + 0.5 * fl1_fx * ts_yyz_yyz_0[j];

                tdy_xyyz_yyz_0[j] = pa_x[j] * tdy_yyz_yyz_0[j];

                tdz_xyyz_yyz_0[j] = pa_x[j] * tdz_yyz_yyz_0[j];

                tdx_xyyz_yzz_0[j] = pa_x[j] * tdx_yyz_yzz_0[j] + 0.5 * fl1_fx * ts_yyz_yzz_0[j];

                tdy_xyyz_yzz_0[j] = pa_x[j] * tdy_yyz_yzz_0[j];

                tdz_xyyz_yzz_0[j] = pa_x[j] * tdz_yyz_yzz_0[j];

                tdx_xyyz_zzz_0[j] = pa_x[j] * tdx_yyz_zzz_0[j] + 0.5 * fl1_fx * ts_yyz_zzz_0[j];

                tdy_xyyz_zzz_0[j] = pa_x[j] * tdy_yyz_zzz_0[j];

                tdz_xyyz_zzz_0[j] = pa_x[j] * tdz_yyz_zzz_0[j];

                tdx_xyzz_xxx_0[j] = pa_x[j] * tdx_yzz_xxx_0[j] + 1.5 * fl1_fx * tdx_yzz_xx_0[j] + 0.5 * fl1_fx * ts_yzz_xxx_0[j];

                tdy_xyzz_xxx_0[j] = pa_x[j] * tdy_yzz_xxx_0[j] + 1.5 * fl1_fx * tdy_yzz_xx_0[j];

                tdz_xyzz_xxx_0[j] = pa_x[j] * tdz_yzz_xxx_0[j] + 1.5 * fl1_fx * tdz_yzz_xx_0[j];

                tdx_xyzz_xxy_0[j] = pa_x[j] * tdx_yzz_xxy_0[j] + fl1_fx * tdx_yzz_xy_0[j] + 0.5 * fl1_fx * ts_yzz_xxy_0[j];

                tdy_xyzz_xxy_0[j] = pa_x[j] * tdy_yzz_xxy_0[j] + fl1_fx * tdy_yzz_xy_0[j];

                tdz_xyzz_xxy_0[j] = pa_x[j] * tdz_yzz_xxy_0[j] + fl1_fx * tdz_yzz_xy_0[j];

                tdx_xyzz_xxz_0[j] = pa_x[j] * tdx_yzz_xxz_0[j] + fl1_fx * tdx_yzz_xz_0[j] + 0.5 * fl1_fx * ts_yzz_xxz_0[j];

                tdy_xyzz_xxz_0[j] = pa_x[j] * tdy_yzz_xxz_0[j] + fl1_fx * tdy_yzz_xz_0[j];

                tdz_xyzz_xxz_0[j] = pa_x[j] * tdz_yzz_xxz_0[j] + fl1_fx * tdz_yzz_xz_0[j];

                tdx_xyzz_xyy_0[j] = pa_x[j] * tdx_yzz_xyy_0[j] + 0.5 * fl1_fx * tdx_yzz_yy_0[j] + 0.5 * fl1_fx * ts_yzz_xyy_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGF_250_300(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (250,300)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdy_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 83); 

            auto tdz_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 83); 

            auto tdx_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 84); 

            auto tdy_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 84); 

            auto tdz_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 84); 

            auto tdx_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 85); 

            auto tdy_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 85); 

            auto tdz_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 85); 

            auto tdx_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 86); 

            auto tdy_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 86); 

            auto tdz_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 86); 

            auto tdx_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 87); 

            auto tdy_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 87); 

            auto tdz_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 87); 

            auto tdx_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 88); 

            auto tdy_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 88); 

            auto tdz_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 88); 

            auto tdx_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 89); 

            auto tdy_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 89); 

            auto tdz_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 89); 

            auto tdx_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 90); 

            auto tdy_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 90); 

            auto tdz_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 90); 

            auto tdx_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 91); 

            auto tdy_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 91); 

            auto tdz_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 91); 

            auto tdx_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 92); 

            auto tdy_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 92); 

            auto tdz_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 92); 

            auto tdx_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 93); 

            auto tdy_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 93); 

            auto tdz_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 93); 

            auto tdx_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 94); 

            auto tdy_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 94); 

            auto tdz_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 94); 

            auto tdx_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 95); 

            auto tdy_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 95); 

            auto tdz_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 95); 

            auto tdx_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 96); 

            auto tdy_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 96); 

            auto tdz_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 96); 

            auto tdx_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 97); 

            auto tdy_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 97); 

            auto tdz_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 97); 

            auto tdx_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 98); 

            auto tdy_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 98); 

            auto tdz_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 98); 

            auto tdx_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 99); 

            auto tdy_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 99); 

            auto tdz_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 99); 

            auto tdy_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 52); 

            auto tdy_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 53); 

            auto tdy_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 53); 

            auto tdz_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tdx_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 54); 

            auto tdy_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdx_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 55); 

            auto tdy_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tdz_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tdx_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 56); 

            auto tdy_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tdz_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tdx_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 57); 

            auto tdy_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tdz_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tdx_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 58); 

            auto tdy_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tdz_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tdx_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 59); 

            auto tdy_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tdz_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            auto ts_yzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 84); 

            auto ts_yzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 85); 

            auto ts_yzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 86); 

            auto ts_yzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 87); 

            auto ts_yzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 88); 

            auto ts_yzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 89); 

            auto ts_zzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 90); 

            auto ts_zzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 91); 

            auto ts_zzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 92); 

            auto ts_zzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 93); 

            auto ts_zzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 94); 

            auto ts_zzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 95); 

            auto ts_zzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 96); 

            auto ts_zzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 97); 

            auto ts_zzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 98); 

            auto ts_zzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 99); 

            // set up pointers to integrals

            auto tdy_xyzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 83); 

            auto tdz_xyzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 83); 

            auto tdx_xyzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 84); 

            auto tdy_xyzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 84); 

            auto tdz_xyzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 84); 

            auto tdx_xyzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 85); 

            auto tdy_xyzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 85); 

            auto tdz_xyzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 85); 

            auto tdx_xyzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 86); 

            auto tdy_xyzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 86); 

            auto tdz_xyzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 86); 

            auto tdx_xyzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 87); 

            auto tdy_xyzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 87); 

            auto tdz_xyzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 87); 

            auto tdx_xyzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 88); 

            auto tdy_xyzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 88); 

            auto tdz_xyzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 88); 

            auto tdx_xyzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 89); 

            auto tdy_xyzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 89); 

            auto tdz_xyzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 89); 

            auto tdx_xzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 90); 

            auto tdy_xzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 90); 

            auto tdz_xzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 90); 

            auto tdx_xzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 91); 

            auto tdy_xzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 91); 

            auto tdz_xzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 91); 

            auto tdx_xzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 92); 

            auto tdy_xzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 92); 

            auto tdz_xzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 92); 

            auto tdx_xzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 93); 

            auto tdy_xzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 93); 

            auto tdz_xzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 93); 

            auto tdx_xzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 94); 

            auto tdy_xzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 94); 

            auto tdz_xzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 94); 

            auto tdx_xzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 95); 

            auto tdy_xzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 95); 

            auto tdz_xzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 95); 

            auto tdx_xzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 96); 

            auto tdy_xzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 96); 

            auto tdz_xzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 96); 

            auto tdx_xzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 97); 

            auto tdy_xzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 97); 

            auto tdz_xzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 97); 

            auto tdx_xzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 98); 

            auto tdy_xzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 98); 

            auto tdz_xzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 98); 

            auto tdx_xzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 99); 

            auto tdy_xzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 99); 

            auto tdz_xzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 99); 

            // Batch of Integrals (250,300)

            #pragma omp simd aligned(fx, pa_x, tdx_xyzz_xyz_0, tdx_xyzz_xzz_0, tdx_xyzz_yyy_0, \
                                     tdx_xyzz_yyz_0, tdx_xyzz_yzz_0, tdx_xyzz_zzz_0, tdx_xzzz_xxx_0, tdx_xzzz_xxy_0, \
                                     tdx_xzzz_xxz_0, tdx_xzzz_xyy_0, tdx_xzzz_xyz_0, tdx_xzzz_xzz_0, tdx_xzzz_yyy_0, \
                                     tdx_xzzz_yyz_0, tdx_xzzz_yzz_0, tdx_xzzz_zzz_0, tdx_yzz_xyz_0, tdx_yzz_xzz_0, \
                                     tdx_yzz_yyy_0, tdx_yzz_yyz_0, tdx_yzz_yz_0, tdx_yzz_yzz_0, tdx_yzz_zz_0, \
                                     tdx_yzz_zzz_0, tdx_zzz_xx_0, tdx_zzz_xxx_0, tdx_zzz_xxy_0, tdx_zzz_xxz_0, \
                                     tdx_zzz_xy_0, tdx_zzz_xyy_0, tdx_zzz_xyz_0, tdx_zzz_xz_0, tdx_zzz_xzz_0, \
                                     tdx_zzz_yy_0, tdx_zzz_yyy_0, tdx_zzz_yyz_0, tdx_zzz_yz_0, tdx_zzz_yzz_0, \
                                     tdx_zzz_zz_0, tdx_zzz_zzz_0, tdy_xyzz_xyy_0, tdy_xyzz_xyz_0, tdy_xyzz_xzz_0, \
                                     tdy_xyzz_yyy_0, tdy_xyzz_yyz_0, tdy_xyzz_yzz_0, tdy_xyzz_zzz_0, tdy_xzzz_xxx_0, \
                                     tdy_xzzz_xxy_0, tdy_xzzz_xxz_0, tdy_xzzz_xyy_0, tdy_xzzz_xyz_0, tdy_xzzz_xzz_0, \
                                     tdy_xzzz_yyy_0, tdy_xzzz_yyz_0, tdy_xzzz_yzz_0, tdy_xzzz_zzz_0, tdy_yzz_xyy_0, \
                                     tdy_yzz_xyz_0, tdy_yzz_xzz_0, tdy_yzz_yy_0, tdy_yzz_yyy_0, tdy_yzz_yyz_0, \
                                     tdy_yzz_yz_0, tdy_yzz_yzz_0, tdy_yzz_zz_0, tdy_yzz_zzz_0, tdy_zzz_xx_0, \
                                     tdy_zzz_xxx_0, tdy_zzz_xxy_0, tdy_zzz_xxz_0, tdy_zzz_xy_0, tdy_zzz_xyy_0, \
                                     tdy_zzz_xyz_0, tdy_zzz_xz_0, tdy_zzz_xzz_0, tdy_zzz_yy_0, tdy_zzz_yyy_0, \
                                     tdy_zzz_yyz_0, tdy_zzz_yz_0, tdy_zzz_yzz_0, tdy_zzz_zz_0, tdy_zzz_zzz_0, \
                                     tdz_xyzz_xyy_0, tdz_xyzz_xyz_0, tdz_xyzz_xzz_0, tdz_xyzz_yyy_0, tdz_xyzz_yyz_0, \
                                     tdz_xyzz_yzz_0, tdz_xyzz_zzz_0, tdz_xzzz_xxx_0, tdz_xzzz_xxy_0, tdz_xzzz_xxz_0, \
                                     tdz_xzzz_xyy_0, tdz_xzzz_xyz_0, tdz_xzzz_xzz_0, tdz_xzzz_yyy_0, tdz_xzzz_yyz_0, \
                                     tdz_xzzz_yzz_0, tdz_xzzz_zzz_0, tdz_yzz_xyy_0, tdz_yzz_xyz_0, tdz_yzz_xzz_0, \
                                     tdz_yzz_yy_0, tdz_yzz_yyy_0, tdz_yzz_yyz_0, tdz_yzz_yz_0, tdz_yzz_yzz_0, \
                                     tdz_yzz_zz_0, tdz_yzz_zzz_0, tdz_zzz_xx_0, tdz_zzz_xxx_0, tdz_zzz_xxy_0, \
                                     tdz_zzz_xxz_0, tdz_zzz_xy_0, tdz_zzz_xyy_0, tdz_zzz_xyz_0, tdz_zzz_xz_0, \
                                     tdz_zzz_xzz_0, tdz_zzz_yy_0, tdz_zzz_yyy_0, tdz_zzz_yyz_0, tdz_zzz_yz_0, \
                                     tdz_zzz_yzz_0, tdz_zzz_zz_0, tdz_zzz_zzz_0, ts_yzz_xyz_0, ts_yzz_xzz_0, \
                                     ts_yzz_yyy_0, ts_yzz_yyz_0, ts_yzz_yzz_0, ts_yzz_zzz_0, ts_zzz_xxx_0, ts_zzz_xxy_0, \
                                     ts_zzz_xxz_0, ts_zzz_xyy_0, ts_zzz_xyz_0, ts_zzz_xzz_0, ts_zzz_yyy_0, ts_zzz_yyz_0, \
                                     ts_zzz_yzz_0, ts_zzz_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdy_xyzz_xyy_0[j] = pa_x[j] * tdy_yzz_xyy_0[j] + 0.5 * fl1_fx * tdy_yzz_yy_0[j];

                tdz_xyzz_xyy_0[j] = pa_x[j] * tdz_yzz_xyy_0[j] + 0.5 * fl1_fx * tdz_yzz_yy_0[j];

                tdx_xyzz_xyz_0[j] = pa_x[j] * tdx_yzz_xyz_0[j] + 0.5 * fl1_fx * tdx_yzz_yz_0[j] + 0.5 * fl1_fx * ts_yzz_xyz_0[j];

                tdy_xyzz_xyz_0[j] = pa_x[j] * tdy_yzz_xyz_0[j] + 0.5 * fl1_fx * tdy_yzz_yz_0[j];

                tdz_xyzz_xyz_0[j] = pa_x[j] * tdz_yzz_xyz_0[j] + 0.5 * fl1_fx * tdz_yzz_yz_0[j];

                tdx_xyzz_xzz_0[j] = pa_x[j] * tdx_yzz_xzz_0[j] + 0.5 * fl1_fx * tdx_yzz_zz_0[j] + 0.5 * fl1_fx * ts_yzz_xzz_0[j];

                tdy_xyzz_xzz_0[j] = pa_x[j] * tdy_yzz_xzz_0[j] + 0.5 * fl1_fx * tdy_yzz_zz_0[j];

                tdz_xyzz_xzz_0[j] = pa_x[j] * tdz_yzz_xzz_0[j] + 0.5 * fl1_fx * tdz_yzz_zz_0[j];

                tdx_xyzz_yyy_0[j] = pa_x[j] * tdx_yzz_yyy_0[j] + 0.5 * fl1_fx * ts_yzz_yyy_0[j];

                tdy_xyzz_yyy_0[j] = pa_x[j] * tdy_yzz_yyy_0[j];

                tdz_xyzz_yyy_0[j] = pa_x[j] * tdz_yzz_yyy_0[j];

                tdx_xyzz_yyz_0[j] = pa_x[j] * tdx_yzz_yyz_0[j] + 0.5 * fl1_fx * ts_yzz_yyz_0[j];

                tdy_xyzz_yyz_0[j] = pa_x[j] * tdy_yzz_yyz_0[j];

                tdz_xyzz_yyz_0[j] = pa_x[j] * tdz_yzz_yyz_0[j];

                tdx_xyzz_yzz_0[j] = pa_x[j] * tdx_yzz_yzz_0[j] + 0.5 * fl1_fx * ts_yzz_yzz_0[j];

                tdy_xyzz_yzz_0[j] = pa_x[j] * tdy_yzz_yzz_0[j];

                tdz_xyzz_yzz_0[j] = pa_x[j] * tdz_yzz_yzz_0[j];

                tdx_xyzz_zzz_0[j] = pa_x[j] * tdx_yzz_zzz_0[j] + 0.5 * fl1_fx * ts_yzz_zzz_0[j];

                tdy_xyzz_zzz_0[j] = pa_x[j] * tdy_yzz_zzz_0[j];

                tdz_xyzz_zzz_0[j] = pa_x[j] * tdz_yzz_zzz_0[j];

                tdx_xzzz_xxx_0[j] = pa_x[j] * tdx_zzz_xxx_0[j] + 1.5 * fl1_fx * tdx_zzz_xx_0[j] + 0.5 * fl1_fx * ts_zzz_xxx_0[j];

                tdy_xzzz_xxx_0[j] = pa_x[j] * tdy_zzz_xxx_0[j] + 1.5 * fl1_fx * tdy_zzz_xx_0[j];

                tdz_xzzz_xxx_0[j] = pa_x[j] * tdz_zzz_xxx_0[j] + 1.5 * fl1_fx * tdz_zzz_xx_0[j];

                tdx_xzzz_xxy_0[j] = pa_x[j] * tdx_zzz_xxy_0[j] + fl1_fx * tdx_zzz_xy_0[j] + 0.5 * fl1_fx * ts_zzz_xxy_0[j];

                tdy_xzzz_xxy_0[j] = pa_x[j] * tdy_zzz_xxy_0[j] + fl1_fx * tdy_zzz_xy_0[j];

                tdz_xzzz_xxy_0[j] = pa_x[j] * tdz_zzz_xxy_0[j] + fl1_fx * tdz_zzz_xy_0[j];

                tdx_xzzz_xxz_0[j] = pa_x[j] * tdx_zzz_xxz_0[j] + fl1_fx * tdx_zzz_xz_0[j] + 0.5 * fl1_fx * ts_zzz_xxz_0[j];

                tdy_xzzz_xxz_0[j] = pa_x[j] * tdy_zzz_xxz_0[j] + fl1_fx * tdy_zzz_xz_0[j];

                tdz_xzzz_xxz_0[j] = pa_x[j] * tdz_zzz_xxz_0[j] + fl1_fx * tdz_zzz_xz_0[j];

                tdx_xzzz_xyy_0[j] = pa_x[j] * tdx_zzz_xyy_0[j] + 0.5 * fl1_fx * tdx_zzz_yy_0[j] + 0.5 * fl1_fx * ts_zzz_xyy_0[j];

                tdy_xzzz_xyy_0[j] = pa_x[j] * tdy_zzz_xyy_0[j] + 0.5 * fl1_fx * tdy_zzz_yy_0[j];

                tdz_xzzz_xyy_0[j] = pa_x[j] * tdz_zzz_xyy_0[j] + 0.5 * fl1_fx * tdz_zzz_yy_0[j];

                tdx_xzzz_xyz_0[j] = pa_x[j] * tdx_zzz_xyz_0[j] + 0.5 * fl1_fx * tdx_zzz_yz_0[j] + 0.5 * fl1_fx * ts_zzz_xyz_0[j];

                tdy_xzzz_xyz_0[j] = pa_x[j] * tdy_zzz_xyz_0[j] + 0.5 * fl1_fx * tdy_zzz_yz_0[j];

                tdz_xzzz_xyz_0[j] = pa_x[j] * tdz_zzz_xyz_0[j] + 0.5 * fl1_fx * tdz_zzz_yz_0[j];

                tdx_xzzz_xzz_0[j] = pa_x[j] * tdx_zzz_xzz_0[j] + 0.5 * fl1_fx * tdx_zzz_zz_0[j] + 0.5 * fl1_fx * ts_zzz_xzz_0[j];

                tdy_xzzz_xzz_0[j] = pa_x[j] * tdy_zzz_xzz_0[j] + 0.5 * fl1_fx * tdy_zzz_zz_0[j];

                tdz_xzzz_xzz_0[j] = pa_x[j] * tdz_zzz_xzz_0[j] + 0.5 * fl1_fx * tdz_zzz_zz_0[j];

                tdx_xzzz_yyy_0[j] = pa_x[j] * tdx_zzz_yyy_0[j] + 0.5 * fl1_fx * ts_zzz_yyy_0[j];

                tdy_xzzz_yyy_0[j] = pa_x[j] * tdy_zzz_yyy_0[j];

                tdz_xzzz_yyy_0[j] = pa_x[j] * tdz_zzz_yyy_0[j];

                tdx_xzzz_yyz_0[j] = pa_x[j] * tdx_zzz_yyz_0[j] + 0.5 * fl1_fx * ts_zzz_yyz_0[j];

                tdy_xzzz_yyz_0[j] = pa_x[j] * tdy_zzz_yyz_0[j];

                tdz_xzzz_yyz_0[j] = pa_x[j] * tdz_zzz_yyz_0[j];

                tdx_xzzz_yzz_0[j] = pa_x[j] * tdx_zzz_yzz_0[j] + 0.5 * fl1_fx * ts_zzz_yzz_0[j];

                tdy_xzzz_yzz_0[j] = pa_x[j] * tdy_zzz_yzz_0[j];

                tdz_xzzz_yzz_0[j] = pa_x[j] * tdz_zzz_yzz_0[j];

                tdx_xzzz_zzz_0[j] = pa_x[j] * tdx_zzz_zzz_0[j] + 0.5 * fl1_fx * ts_zzz_zzz_0[j];

                tdy_xzzz_zzz_0[j] = pa_x[j] * tdy_zzz_zzz_0[j];

                tdz_xzzz_zzz_0[j] = pa_x[j] * tdz_zzz_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGF_300_350(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (300,350)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tdx_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 60); 

            auto tdy_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 60); 

            auto tdz_yyy_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 60); 

            auto tdx_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 61); 

            auto tdy_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 61); 

            auto tdz_yyy_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 61); 

            auto tdx_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 62); 

            auto tdy_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 62); 

            auto tdz_yyy_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 62); 

            auto tdx_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 63); 

            auto tdy_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 63); 

            auto tdz_yyy_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 63); 

            auto tdx_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 64); 

            auto tdy_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 64); 

            auto tdz_yyy_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 64); 

            auto tdx_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 65); 

            auto tdy_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 65); 

            auto tdz_yyy_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 65); 

            auto tdx_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 66); 

            auto tdy_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 66); 

            auto tdz_yyy_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 66); 

            auto tdx_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 67); 

            auto tdy_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 67); 

            auto tdz_yyy_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 67); 

            auto tdx_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 68); 

            auto tdy_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 68); 

            auto tdz_yyy_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 68); 

            auto tdx_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 69); 

            auto tdy_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 69); 

            auto tdz_yyy_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 69); 

            auto tdx_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 70); 

            auto tdy_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 70); 

            auto tdz_yyz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 70); 

            auto tdx_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 71); 

            auto tdy_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 71); 

            auto tdz_yyz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 71); 

            auto tdx_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 72); 

            auto tdy_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 72); 

            auto tdz_yyz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 72); 

            auto tdx_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 73); 

            auto tdy_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 73); 

            auto tdz_yyz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 73); 

            auto tdx_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 74); 

            auto tdy_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 74); 

            auto tdz_yyz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 74); 

            auto tdx_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 75); 

            auto tdy_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 75); 

            auto tdz_yyz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 75); 

            auto tdx_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 76); 

            auto tdy_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 76); 

            auto tdx_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 30); 

            auto tdy_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 30); 

            auto tdz_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 30); 

            auto tdx_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 31); 

            auto tdy_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 31); 

            auto tdz_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 31); 

            auto tdx_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 32); 

            auto tdy_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 32); 

            auto tdz_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 32); 

            auto tdx_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 33); 

            auto tdy_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 33); 

            auto tdz_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 33); 

            auto tdx_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 34); 

            auto tdy_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 34); 

            auto tdz_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 34); 

            auto tdx_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 35); 

            auto tdy_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 35); 

            auto tdz_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 35); 

            auto tdx_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 36); 

            auto tdy_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 37); 

            auto tdy_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 38); 

            auto tdy_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 39); 

            auto tdy_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 40); 

            auto tdy_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 41); 

            auto tdy_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdx_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 42); 

            auto tdy_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 43); 

            auto tdy_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 44); 

            auto tdy_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 44); 

            auto tdx_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 45); 

            auto tdy_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 45); 

            auto tdz_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 46); 

            auto tdy_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdx_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 36); 

            auto tdy_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 37); 

            auto tdy_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 38); 

            auto tdy_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 39); 

            auto tdy_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 40); 

            auto tdy_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 41); 

            auto tdy_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdx_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 42); 

            auto tdy_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 43); 

            auto tdy_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 44); 

            auto tdy_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            auto tdx_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 45); 

            auto tdy_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto ts_yyy_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 60); 

            auto ts_yyy_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 61); 

            auto ts_yyy_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 62); 

            auto ts_yyy_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 63); 

            auto ts_yyy_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 64); 

            auto ts_yyy_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 65); 

            auto ts_yyy_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 66); 

            auto ts_yyy_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 67); 

            auto ts_yyy_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 68); 

            auto ts_yyy_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 69); 

            auto ts_yyz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 70); 

            auto ts_yyz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 71); 

            auto ts_yyz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 72); 

            auto ts_yyz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 73); 

            auto ts_yyz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 74); 

            auto ts_yyz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 75); 

            auto ts_yyz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 76); 

            // set up pointers to integrals

            auto tdx_yyyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 100); 

            auto tdy_yyyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 100); 

            auto tdz_yyyy_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 100); 

            auto tdx_yyyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 101); 

            auto tdy_yyyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 101); 

            auto tdz_yyyy_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 101); 

            auto tdx_yyyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 102); 

            auto tdy_yyyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 102); 

            auto tdz_yyyy_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 102); 

            auto tdx_yyyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 103); 

            auto tdy_yyyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 103); 

            auto tdz_yyyy_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 103); 

            auto tdx_yyyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 104); 

            auto tdy_yyyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 104); 

            auto tdz_yyyy_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 104); 

            auto tdx_yyyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 105); 

            auto tdy_yyyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 105); 

            auto tdz_yyyy_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 105); 

            auto tdx_yyyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 106); 

            auto tdy_yyyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 106); 

            auto tdz_yyyy_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 106); 

            auto tdx_yyyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 107); 

            auto tdy_yyyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 107); 

            auto tdz_yyyy_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 107); 

            auto tdx_yyyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 108); 

            auto tdy_yyyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 108); 

            auto tdz_yyyy_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 108); 

            auto tdx_yyyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 109); 

            auto tdy_yyyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 109); 

            auto tdz_yyyy_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 109); 

            auto tdx_yyyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 110); 

            auto tdy_yyyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 110); 

            auto tdz_yyyz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 110); 

            auto tdx_yyyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 111); 

            auto tdy_yyyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 111); 

            auto tdz_yyyz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 111); 

            auto tdx_yyyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 112); 

            auto tdy_yyyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 112); 

            auto tdz_yyyz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 112); 

            auto tdx_yyyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 113); 

            auto tdy_yyyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 113); 

            auto tdz_yyyz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 113); 

            auto tdx_yyyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 114); 

            auto tdy_yyyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 114); 

            auto tdz_yyyz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 114); 

            auto tdx_yyyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 115); 

            auto tdy_yyyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 115); 

            auto tdz_yyyz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 115); 

            auto tdx_yyyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 116); 

            auto tdy_yyyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 116); 

            // Batch of Integrals (300,350)

            #pragma omp simd aligned(fx, pa_y, tdx_yy_xxx_0, tdx_yy_xxy_0, tdx_yy_xxz_0, tdx_yy_xyy_0, \
                                     tdx_yy_xyz_0, tdx_yy_xzz_0, tdx_yy_yyy_0, tdx_yy_yyz_0, tdx_yy_yzz_0, tdx_yy_zzz_0, \
                                     tdx_yyy_xx_0, tdx_yyy_xxx_0, tdx_yyy_xxy_0, tdx_yyy_xxz_0, tdx_yyy_xy_0, \
                                     tdx_yyy_xyy_0, tdx_yyy_xyz_0, tdx_yyy_xz_0, tdx_yyy_xzz_0, tdx_yyy_yy_0, \
                                     tdx_yyy_yyy_0, tdx_yyy_yyz_0, tdx_yyy_yz_0, tdx_yyy_yzz_0, tdx_yyy_zz_0, \
                                     tdx_yyy_zzz_0, tdx_yyyy_xxx_0, tdx_yyyy_xxy_0, tdx_yyyy_xxz_0, tdx_yyyy_xyy_0, \
                                     tdx_yyyy_xyz_0, tdx_yyyy_xzz_0, tdx_yyyy_yyy_0, tdx_yyyy_yyz_0, tdx_yyyy_yzz_0, \
                                     tdx_yyyy_zzz_0, tdx_yyyz_xxx_0, tdx_yyyz_xxy_0, tdx_yyyz_xxz_0, tdx_yyyz_xyy_0, \
                                     tdx_yyyz_xyz_0, tdx_yyyz_xzz_0, tdx_yyyz_yyy_0, tdx_yyz_xx_0, tdx_yyz_xxx_0, \
                                     tdx_yyz_xxy_0, tdx_yyz_xxz_0, tdx_yyz_xy_0, tdx_yyz_xyy_0, tdx_yyz_xyz_0, \
                                     tdx_yyz_xz_0, tdx_yyz_xzz_0, tdx_yyz_yy_0, tdx_yyz_yyy_0, tdx_yz_xxx_0, \
                                     tdx_yz_xxy_0, tdx_yz_xxz_0, tdx_yz_xyy_0, tdx_yz_xyz_0, tdx_yz_xzz_0, tdx_yz_yyy_0, \
                                     tdy_yy_xxx_0, tdy_yy_xxy_0, tdy_yy_xxz_0, tdy_yy_xyy_0, tdy_yy_xyz_0, tdy_yy_xzz_0, \
                                     tdy_yy_yyy_0, tdy_yy_yyz_0, tdy_yy_yzz_0, tdy_yy_zzz_0, tdy_yyy_xx_0, \
                                     tdy_yyy_xxx_0, tdy_yyy_xxy_0, tdy_yyy_xxz_0, tdy_yyy_xy_0, tdy_yyy_xyy_0, \
                                     tdy_yyy_xyz_0, tdy_yyy_xz_0, tdy_yyy_xzz_0, tdy_yyy_yy_0, tdy_yyy_yyy_0, \
                                     tdy_yyy_yyz_0, tdy_yyy_yz_0, tdy_yyy_yzz_0, tdy_yyy_zz_0, tdy_yyy_zzz_0, \
                                     tdy_yyyy_xxx_0, tdy_yyyy_xxy_0, tdy_yyyy_xxz_0, tdy_yyyy_xyy_0, tdy_yyyy_xyz_0, \
                                     tdy_yyyy_xzz_0, tdy_yyyy_yyy_0, tdy_yyyy_yyz_0, tdy_yyyy_yzz_0, tdy_yyyy_zzz_0, \
                                     tdy_yyyz_xxx_0, tdy_yyyz_xxy_0, tdy_yyyz_xxz_0, tdy_yyyz_xyy_0, tdy_yyyz_xyz_0, \
                                     tdy_yyyz_xzz_0, tdy_yyyz_yyy_0, tdy_yyz_xx_0, tdy_yyz_xxx_0, tdy_yyz_xxy_0, \
                                     tdy_yyz_xxz_0, tdy_yyz_xy_0, tdy_yyz_xyy_0, tdy_yyz_xyz_0, tdy_yyz_xz_0, \
                                     tdy_yyz_xzz_0, tdy_yyz_yy_0, tdy_yyz_yyy_0, tdy_yz_xxx_0, tdy_yz_xxy_0, \
                                     tdy_yz_xxz_0, tdy_yz_xyy_0, tdy_yz_xyz_0, tdy_yz_xzz_0, tdy_yz_yyy_0, tdz_yy_xxx_0, \
                                     tdz_yy_xxy_0, tdz_yy_xxz_0, tdz_yy_xyy_0, tdz_yy_xyz_0, tdz_yy_xzz_0, tdz_yy_yyy_0, \
                                     tdz_yy_yyz_0, tdz_yy_yzz_0, tdz_yy_zzz_0, tdz_yyy_xx_0, tdz_yyy_xxx_0, \
                                     tdz_yyy_xxy_0, tdz_yyy_xxz_0, tdz_yyy_xy_0, tdz_yyy_xyy_0, tdz_yyy_xyz_0, \
                                     tdz_yyy_xz_0, tdz_yyy_xzz_0, tdz_yyy_yy_0, tdz_yyy_yyy_0, tdz_yyy_yyz_0, \
                                     tdz_yyy_yz_0, tdz_yyy_yzz_0, tdz_yyy_zz_0, tdz_yyy_zzz_0, tdz_yyyy_xxx_0, \
                                     tdz_yyyy_xxy_0, tdz_yyyy_xxz_0, tdz_yyyy_xyy_0, tdz_yyyy_xyz_0, tdz_yyyy_xzz_0, \
                                     tdz_yyyy_yyy_0, tdz_yyyy_yyz_0, tdz_yyyy_yzz_0, tdz_yyyy_zzz_0, tdz_yyyz_xxx_0, \
                                     tdz_yyyz_xxy_0, tdz_yyyz_xxz_0, tdz_yyyz_xyy_0, tdz_yyyz_xyz_0, tdz_yyyz_xzz_0, \
                                     tdz_yyz_xx_0, tdz_yyz_xxx_0, tdz_yyz_xxy_0, tdz_yyz_xxz_0, tdz_yyz_xy_0, \
                                     tdz_yyz_xyy_0, tdz_yyz_xyz_0, tdz_yyz_xz_0, tdz_yyz_xzz_0, tdz_yz_xxx_0, \
                                     tdz_yz_xxy_0, tdz_yz_xxz_0, tdz_yz_xyy_0, tdz_yz_xyz_0, tdz_yz_xzz_0, ts_yyy_xxx_0, \
                                     ts_yyy_xxy_0, ts_yyy_xxz_0, ts_yyy_xyy_0, ts_yyy_xyz_0, ts_yyy_xzz_0, ts_yyy_yyy_0, \
                                     ts_yyy_yyz_0, ts_yyy_yzz_0, ts_yyy_zzz_0, ts_yyz_xxx_0, ts_yyz_xxy_0, ts_yyz_xxz_0, \
                                     ts_yyz_xyy_0, ts_yyz_xyz_0, ts_yyz_xzz_0, ts_yyz_yyy_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yyyy_xxx_0[j] = pa_y[j] * tdx_yyy_xxx_0[j] + 1.5 * fl1_fx * tdx_yy_xxx_0[j];

                tdy_yyyy_xxx_0[j] = pa_y[j] * tdy_yyy_xxx_0[j] + 1.5 * fl1_fx * tdy_yy_xxx_0[j] + 0.5 * fl1_fx * ts_yyy_xxx_0[j];

                tdz_yyyy_xxx_0[j] = pa_y[j] * tdz_yyy_xxx_0[j] + 1.5 * fl1_fx * tdz_yy_xxx_0[j];

                tdx_yyyy_xxy_0[j] = pa_y[j] * tdx_yyy_xxy_0[j] + 1.5 * fl1_fx * tdx_yy_xxy_0[j] + 0.5 * fl1_fx * tdx_yyy_xx_0[j];

                tdy_yyyy_xxy_0[j] = pa_y[j] * tdy_yyy_xxy_0[j] + 1.5 * fl1_fx * tdy_yy_xxy_0[j] + 0.5 * fl1_fx * tdy_yyy_xx_0[j] + 0.5 * fl1_fx * ts_yyy_xxy_0[j];

                tdz_yyyy_xxy_0[j] = pa_y[j] * tdz_yyy_xxy_0[j] + 1.5 * fl1_fx * tdz_yy_xxy_0[j] + 0.5 * fl1_fx * tdz_yyy_xx_0[j];

                tdx_yyyy_xxz_0[j] = pa_y[j] * tdx_yyy_xxz_0[j] + 1.5 * fl1_fx * tdx_yy_xxz_0[j];

                tdy_yyyy_xxz_0[j] = pa_y[j] * tdy_yyy_xxz_0[j] + 1.5 * fl1_fx * tdy_yy_xxz_0[j] + 0.5 * fl1_fx * ts_yyy_xxz_0[j];

                tdz_yyyy_xxz_0[j] = pa_y[j] * tdz_yyy_xxz_0[j] + 1.5 * fl1_fx * tdz_yy_xxz_0[j];

                tdx_yyyy_xyy_0[j] = pa_y[j] * tdx_yyy_xyy_0[j] + 1.5 * fl1_fx * tdx_yy_xyy_0[j] + fl1_fx * tdx_yyy_xy_0[j];

                tdy_yyyy_xyy_0[j] = pa_y[j] * tdy_yyy_xyy_0[j] + 1.5 * fl1_fx * tdy_yy_xyy_0[j] + fl1_fx * tdy_yyy_xy_0[j] + 0.5 * fl1_fx * ts_yyy_xyy_0[j];

                tdz_yyyy_xyy_0[j] = pa_y[j] * tdz_yyy_xyy_0[j] + 1.5 * fl1_fx * tdz_yy_xyy_0[j] + fl1_fx * tdz_yyy_xy_0[j];

                tdx_yyyy_xyz_0[j] = pa_y[j] * tdx_yyy_xyz_0[j] + 1.5 * fl1_fx * tdx_yy_xyz_0[j] + 0.5 * fl1_fx * tdx_yyy_xz_0[j];

                tdy_yyyy_xyz_0[j] = pa_y[j] * tdy_yyy_xyz_0[j] + 1.5 * fl1_fx * tdy_yy_xyz_0[j] + 0.5 * fl1_fx * tdy_yyy_xz_0[j] + 0.5 * fl1_fx * ts_yyy_xyz_0[j];

                tdz_yyyy_xyz_0[j] = pa_y[j] * tdz_yyy_xyz_0[j] + 1.5 * fl1_fx * tdz_yy_xyz_0[j] + 0.5 * fl1_fx * tdz_yyy_xz_0[j];

                tdx_yyyy_xzz_0[j] = pa_y[j] * tdx_yyy_xzz_0[j] + 1.5 * fl1_fx * tdx_yy_xzz_0[j];

                tdy_yyyy_xzz_0[j] = pa_y[j] * tdy_yyy_xzz_0[j] + 1.5 * fl1_fx * tdy_yy_xzz_0[j] + 0.5 * fl1_fx * ts_yyy_xzz_0[j];

                tdz_yyyy_xzz_0[j] = pa_y[j] * tdz_yyy_xzz_0[j] + 1.5 * fl1_fx * tdz_yy_xzz_0[j];

                tdx_yyyy_yyy_0[j] = pa_y[j] * tdx_yyy_yyy_0[j] + 1.5 * fl1_fx * tdx_yy_yyy_0[j] + 1.5 * fl1_fx * tdx_yyy_yy_0[j];

                tdy_yyyy_yyy_0[j] = pa_y[j] * tdy_yyy_yyy_0[j] + 1.5 * fl1_fx * tdy_yy_yyy_0[j] + 1.5 * fl1_fx * tdy_yyy_yy_0[j] + 0.5 * fl1_fx * ts_yyy_yyy_0[j];

                tdz_yyyy_yyy_0[j] = pa_y[j] * tdz_yyy_yyy_0[j] + 1.5 * fl1_fx * tdz_yy_yyy_0[j] + 1.5 * fl1_fx * tdz_yyy_yy_0[j];

                tdx_yyyy_yyz_0[j] = pa_y[j] * tdx_yyy_yyz_0[j] + 1.5 * fl1_fx * tdx_yy_yyz_0[j] + fl1_fx * tdx_yyy_yz_0[j];

                tdy_yyyy_yyz_0[j] = pa_y[j] * tdy_yyy_yyz_0[j] + 1.5 * fl1_fx * tdy_yy_yyz_0[j] + fl1_fx * tdy_yyy_yz_0[j] + 0.5 * fl1_fx * ts_yyy_yyz_0[j];

                tdz_yyyy_yyz_0[j] = pa_y[j] * tdz_yyy_yyz_0[j] + 1.5 * fl1_fx * tdz_yy_yyz_0[j] + fl1_fx * tdz_yyy_yz_0[j];

                tdx_yyyy_yzz_0[j] = pa_y[j] * tdx_yyy_yzz_0[j] + 1.5 * fl1_fx * tdx_yy_yzz_0[j] + 0.5 * fl1_fx * tdx_yyy_zz_0[j];

                tdy_yyyy_yzz_0[j] = pa_y[j] * tdy_yyy_yzz_0[j] + 1.5 * fl1_fx * tdy_yy_yzz_0[j] + 0.5 * fl1_fx * tdy_yyy_zz_0[j] + 0.5 * fl1_fx * ts_yyy_yzz_0[j];

                tdz_yyyy_yzz_0[j] = pa_y[j] * tdz_yyy_yzz_0[j] + 1.5 * fl1_fx * tdz_yy_yzz_0[j] + 0.5 * fl1_fx * tdz_yyy_zz_0[j];

                tdx_yyyy_zzz_0[j] = pa_y[j] * tdx_yyy_zzz_0[j] + 1.5 * fl1_fx * tdx_yy_zzz_0[j];

                tdy_yyyy_zzz_0[j] = pa_y[j] * tdy_yyy_zzz_0[j] + 1.5 * fl1_fx * tdy_yy_zzz_0[j] + 0.5 * fl1_fx * ts_yyy_zzz_0[j];

                tdz_yyyy_zzz_0[j] = pa_y[j] * tdz_yyy_zzz_0[j] + 1.5 * fl1_fx * tdz_yy_zzz_0[j];

                tdx_yyyz_xxx_0[j] = pa_y[j] * tdx_yyz_xxx_0[j] + fl1_fx * tdx_yz_xxx_0[j];

                tdy_yyyz_xxx_0[j] = pa_y[j] * tdy_yyz_xxx_0[j] + fl1_fx * tdy_yz_xxx_0[j] + 0.5 * fl1_fx * ts_yyz_xxx_0[j];

                tdz_yyyz_xxx_0[j] = pa_y[j] * tdz_yyz_xxx_0[j] + fl1_fx * tdz_yz_xxx_0[j];

                tdx_yyyz_xxy_0[j] = pa_y[j] * tdx_yyz_xxy_0[j] + fl1_fx * tdx_yz_xxy_0[j] + 0.5 * fl1_fx * tdx_yyz_xx_0[j];

                tdy_yyyz_xxy_0[j] = pa_y[j] * tdy_yyz_xxy_0[j] + fl1_fx * tdy_yz_xxy_0[j] + 0.5 * fl1_fx * tdy_yyz_xx_0[j] + 0.5 * fl1_fx * ts_yyz_xxy_0[j];

                tdz_yyyz_xxy_0[j] = pa_y[j] * tdz_yyz_xxy_0[j] + fl1_fx * tdz_yz_xxy_0[j] + 0.5 * fl1_fx * tdz_yyz_xx_0[j];

                tdx_yyyz_xxz_0[j] = pa_y[j] * tdx_yyz_xxz_0[j] + fl1_fx * tdx_yz_xxz_0[j];

                tdy_yyyz_xxz_0[j] = pa_y[j] * tdy_yyz_xxz_0[j] + fl1_fx * tdy_yz_xxz_0[j] + 0.5 * fl1_fx * ts_yyz_xxz_0[j];

                tdz_yyyz_xxz_0[j] = pa_y[j] * tdz_yyz_xxz_0[j] + fl1_fx * tdz_yz_xxz_0[j];

                tdx_yyyz_xyy_0[j] = pa_y[j] * tdx_yyz_xyy_0[j] + fl1_fx * tdx_yz_xyy_0[j] + fl1_fx * tdx_yyz_xy_0[j];

                tdy_yyyz_xyy_0[j] = pa_y[j] * tdy_yyz_xyy_0[j] + fl1_fx * tdy_yz_xyy_0[j] + fl1_fx * tdy_yyz_xy_0[j] + 0.5 * fl1_fx * ts_yyz_xyy_0[j];

                tdz_yyyz_xyy_0[j] = pa_y[j] * tdz_yyz_xyy_0[j] + fl1_fx * tdz_yz_xyy_0[j] + fl1_fx * tdz_yyz_xy_0[j];

                tdx_yyyz_xyz_0[j] = pa_y[j] * tdx_yyz_xyz_0[j] + fl1_fx * tdx_yz_xyz_0[j] + 0.5 * fl1_fx * tdx_yyz_xz_0[j];

                tdy_yyyz_xyz_0[j] = pa_y[j] * tdy_yyz_xyz_0[j] + fl1_fx * tdy_yz_xyz_0[j] + 0.5 * fl1_fx * tdy_yyz_xz_0[j] + 0.5 * fl1_fx * ts_yyz_xyz_0[j];

                tdz_yyyz_xyz_0[j] = pa_y[j] * tdz_yyz_xyz_0[j] + fl1_fx * tdz_yz_xyz_0[j] + 0.5 * fl1_fx * tdz_yyz_xz_0[j];

                tdx_yyyz_xzz_0[j] = pa_y[j] * tdx_yyz_xzz_0[j] + fl1_fx * tdx_yz_xzz_0[j];

                tdy_yyyz_xzz_0[j] = pa_y[j] * tdy_yyz_xzz_0[j] + fl1_fx * tdy_yz_xzz_0[j] + 0.5 * fl1_fx * ts_yyz_xzz_0[j];

                tdz_yyyz_xzz_0[j] = pa_y[j] * tdz_yyz_xzz_0[j] + fl1_fx * tdz_yz_xzz_0[j];

                tdx_yyyz_yyy_0[j] = pa_y[j] * tdx_yyz_yyy_0[j] + fl1_fx * tdx_yz_yyy_0[j] + 1.5 * fl1_fx * tdx_yyz_yy_0[j];

                tdy_yyyz_yyy_0[j] = pa_y[j] * tdy_yyz_yyy_0[j] + fl1_fx * tdy_yz_yyy_0[j] + 1.5 * fl1_fx * tdy_yyz_yy_0[j] + 0.5 * fl1_fx * ts_yyz_yyy_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGF_350_400(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (350,400)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tdz_yyz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 76); 

            auto tdx_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 77); 

            auto tdy_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 77); 

            auto tdz_yyz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 77); 

            auto tdx_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 78); 

            auto tdy_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 78); 

            auto tdz_yyz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 78); 

            auto tdx_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 79); 

            auto tdy_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 79); 

            auto tdz_yyz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 79); 

            auto tdx_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 80); 

            auto tdy_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 80); 

            auto tdz_yzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 80); 

            auto tdx_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 81); 

            auto tdy_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 81); 

            auto tdz_yzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 81); 

            auto tdx_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 82); 

            auto tdy_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 82); 

            auto tdz_yzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 82); 

            auto tdx_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 83); 

            auto tdy_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 83); 

            auto tdz_yzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 83); 

            auto tdx_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 84); 

            auto tdy_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 84); 

            auto tdz_yzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 84); 

            auto tdx_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 85); 

            auto tdy_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 85); 

            auto tdz_yzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 85); 

            auto tdx_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 86); 

            auto tdy_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 86); 

            auto tdz_yzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 86); 

            auto tdx_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 87); 

            auto tdy_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 87); 

            auto tdz_yzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 87); 

            auto tdx_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 88); 

            auto tdy_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 88); 

            auto tdz_yzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 88); 

            auto tdx_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 89); 

            auto tdy_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 89); 

            auto tdz_yzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 89); 

            auto tdx_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 90); 

            auto tdy_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 90); 

            auto tdz_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 90); 

            auto tdx_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 91); 

            auto tdy_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 91); 

            auto tdz_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 91); 

            auto tdx_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 92); 

            auto tdy_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 92); 

            auto tdz_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 92); 

            auto tdx_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 93); 

            auto tdz_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 47); 

            auto tdy_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 48); 

            auto tdy_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 49); 

            auto tdy_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50); 

            auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51); 

            auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52); 

            auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 53); 

            auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53); 

            auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53); 

            auto tdx_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 54); 

            auto tdy_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdx_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 55); 

            auto tdy_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 55); 

            auto tdz_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 55); 

            auto tdx_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 56); 

            auto tdy_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 56); 

            auto tdz_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 56); 

            auto tdx_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 57); 

            auto tdy_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 57); 

            auto tdz_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 57); 

            auto tdx_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 58); 

            auto tdy_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 58); 

            auto tdz_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 58); 

            auto tdx_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 59); 

            auto tdy_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 59); 

            auto tdz_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 59); 

            auto tdz_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 46); 

            auto tdy_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdz_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 47); 

            auto tdy_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 48); 

            auto tdy_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 49); 

            auto tdy_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdx_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 50); 

            auto tdy_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 51); 

            auto tdy_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 52); 

            auto tdy_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 53); 

            auto tdy_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 53); 

            auto tdz_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tdx_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 54); 

            auto tdy_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdx_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 55); 

            auto ts_yyz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 77); 

            auto ts_yyz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 78); 

            auto ts_yyz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 79); 

            auto ts_yzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 80); 

            auto ts_yzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 81); 

            auto ts_yzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 82); 

            auto ts_yzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 83); 

            auto ts_yzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 84); 

            auto ts_yzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 85); 

            auto ts_yzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 86); 

            auto ts_yzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 87); 

            auto ts_yzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 88); 

            auto ts_yzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 89); 

            auto ts_zzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 90); 

            auto ts_zzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 91); 

            auto ts_zzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 92); 

            // set up pointers to integrals

            auto tdz_yyyz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 116); 

            auto tdx_yyyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 117); 

            auto tdy_yyyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 117); 

            auto tdz_yyyz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 117); 

            auto tdx_yyyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 118); 

            auto tdy_yyyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 118); 

            auto tdz_yyyz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 118); 

            auto tdx_yyyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 119); 

            auto tdy_yyyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 119); 

            auto tdz_yyyz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 119); 

            auto tdx_yyzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 120); 

            auto tdy_yyzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 120); 

            auto tdz_yyzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 120); 

            auto tdx_yyzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 121); 

            auto tdy_yyzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 121); 

            auto tdz_yyzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 121); 

            auto tdx_yyzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 122); 

            auto tdy_yyzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 122); 

            auto tdz_yyzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 122); 

            auto tdx_yyzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 123); 

            auto tdy_yyzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 123); 

            auto tdz_yyzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 123); 

            auto tdx_yyzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 124); 

            auto tdy_yyzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 124); 

            auto tdz_yyzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 124); 

            auto tdx_yyzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 125); 

            auto tdy_yyzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 125); 

            auto tdz_yyzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 125); 

            auto tdx_yyzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 126); 

            auto tdy_yyzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 126); 

            auto tdz_yyzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 126); 

            auto tdx_yyzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 127); 

            auto tdy_yyzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 127); 

            auto tdz_yyzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 127); 

            auto tdx_yyzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 128); 

            auto tdy_yyzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 128); 

            auto tdz_yyzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 128); 

            auto tdx_yyzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 129); 

            auto tdy_yyzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 129); 

            auto tdz_yyzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 129); 

            auto tdx_yzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 130); 

            auto tdy_yzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 130); 

            auto tdz_yzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 130); 

            auto tdx_yzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 131); 

            auto tdy_yzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 131); 

            auto tdz_yzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 131); 

            auto tdx_yzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 132); 

            auto tdy_yzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 132); 

            auto tdz_yzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 132); 

            auto tdx_yzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 133); 

            // Batch of Integrals (350,400)

            #pragma omp simd aligned(fx, pa_y, tdx_yyyz_yyz_0, tdx_yyyz_yzz_0, tdx_yyyz_zzz_0, \
                                     tdx_yyz_yyz_0, tdx_yyz_yz_0, tdx_yyz_yzz_0, tdx_yyz_zz_0, tdx_yyz_zzz_0, \
                                     tdx_yyzz_xxx_0, tdx_yyzz_xxy_0, tdx_yyzz_xxz_0, tdx_yyzz_xyy_0, tdx_yyzz_xyz_0, \
                                     tdx_yyzz_xzz_0, tdx_yyzz_yyy_0, tdx_yyzz_yyz_0, tdx_yyzz_yzz_0, tdx_yyzz_zzz_0, \
                                     tdx_yz_yyz_0, tdx_yz_yzz_0, tdx_yz_zzz_0, tdx_yzz_xx_0, tdx_yzz_xxx_0, \
                                     tdx_yzz_xxy_0, tdx_yzz_xxz_0, tdx_yzz_xy_0, tdx_yzz_xyy_0, tdx_yzz_xyz_0, \
                                     tdx_yzz_xz_0, tdx_yzz_xzz_0, tdx_yzz_yy_0, tdx_yzz_yyy_0, tdx_yzz_yyz_0, \
                                     tdx_yzz_yz_0, tdx_yzz_yzz_0, tdx_yzz_zz_0, tdx_yzz_zzz_0, tdx_yzzz_xxx_0, \
                                     tdx_yzzz_xxy_0, tdx_yzzz_xxz_0, tdx_yzzz_xyy_0, tdx_zz_xxx_0, tdx_zz_xxy_0, \
                                     tdx_zz_xxz_0, tdx_zz_xyy_0, tdx_zz_xyz_0, tdx_zz_xzz_0, tdx_zz_yyy_0, tdx_zz_yyz_0, \
                                     tdx_zz_yzz_0, tdx_zz_zzz_0, tdx_zzz_xx_0, tdx_zzz_xxx_0, tdx_zzz_xxy_0, \
                                     tdx_zzz_xxz_0, tdx_zzz_xy_0, tdx_zzz_xyy_0, tdy_yyyz_yyz_0, tdy_yyyz_yzz_0, \
                                     tdy_yyyz_zzz_0, tdy_yyz_yyz_0, tdy_yyz_yz_0, tdy_yyz_yzz_0, tdy_yyz_zz_0, \
                                     tdy_yyz_zzz_0, tdy_yyzz_xxx_0, tdy_yyzz_xxy_0, tdy_yyzz_xxz_0, tdy_yyzz_xyy_0, \
                                     tdy_yyzz_xyz_0, tdy_yyzz_xzz_0, tdy_yyzz_yyy_0, tdy_yyzz_yyz_0, tdy_yyzz_yzz_0, \
                                     tdy_yyzz_zzz_0, tdy_yz_yyz_0, tdy_yz_yzz_0, tdy_yz_zzz_0, tdy_yzz_xx_0, \
                                     tdy_yzz_xxx_0, tdy_yzz_xxy_0, tdy_yzz_xxz_0, tdy_yzz_xy_0, tdy_yzz_xyy_0, \
                                     tdy_yzz_xyz_0, tdy_yzz_xz_0, tdy_yzz_xzz_0, tdy_yzz_yy_0, tdy_yzz_yyy_0, \
                                     tdy_yzz_yyz_0, tdy_yzz_yz_0, tdy_yzz_yzz_0, tdy_yzz_zz_0, tdy_yzz_zzz_0, \
                                     tdy_yzzz_xxx_0, tdy_yzzz_xxy_0, tdy_yzzz_xxz_0, tdy_zz_xxx_0, tdy_zz_xxy_0, \
                                     tdy_zz_xxz_0, tdy_zz_xyy_0, tdy_zz_xyz_0, tdy_zz_xzz_0, tdy_zz_yyy_0, tdy_zz_yyz_0, \
                                     tdy_zz_yzz_0, tdy_zz_zzz_0, tdy_zzz_xx_0, tdy_zzz_xxx_0, tdy_zzz_xxy_0, \
                                     tdy_zzz_xxz_0, tdz_yyyz_yyy_0, tdz_yyyz_yyz_0, tdz_yyyz_yzz_0, tdz_yyyz_zzz_0, \
                                     tdz_yyz_yy_0, tdz_yyz_yyy_0, tdz_yyz_yyz_0, tdz_yyz_yz_0, tdz_yyz_yzz_0, \
                                     tdz_yyz_zz_0, tdz_yyz_zzz_0, tdz_yyzz_xxx_0, tdz_yyzz_xxy_0, tdz_yyzz_xxz_0, \
                                     tdz_yyzz_xyy_0, tdz_yyzz_xyz_0, tdz_yyzz_xzz_0, tdz_yyzz_yyy_0, tdz_yyzz_yyz_0, \
                                     tdz_yyzz_yzz_0, tdz_yyzz_zzz_0, tdz_yz_yyy_0, tdz_yz_yyz_0, tdz_yz_yzz_0, \
                                     tdz_yz_zzz_0, tdz_yzz_xx_0, tdz_yzz_xxx_0, tdz_yzz_xxy_0, tdz_yzz_xxz_0, \
                                     tdz_yzz_xy_0, tdz_yzz_xyy_0, tdz_yzz_xyz_0, tdz_yzz_xz_0, tdz_yzz_xzz_0, \
                                     tdz_yzz_yy_0, tdz_yzz_yyy_0, tdz_yzz_yyz_0, tdz_yzz_yz_0, tdz_yzz_yzz_0, \
                                     tdz_yzz_zz_0, tdz_yzz_zzz_0, tdz_yzzz_xxx_0, tdz_yzzz_xxy_0, tdz_yzzz_xxz_0, \
                                     tdz_zz_xxx_0, tdz_zz_xxy_0, tdz_zz_xxz_0, tdz_zz_xyy_0, tdz_zz_xyz_0, tdz_zz_xzz_0, \
                                     tdz_zz_yyy_0, tdz_zz_yyz_0, tdz_zz_yzz_0, tdz_zz_zzz_0, tdz_zzz_xx_0, \
                                     tdz_zzz_xxx_0, tdz_zzz_xxy_0, tdz_zzz_xxz_0, ts_yyz_yyz_0, ts_yyz_yzz_0, \
                                     ts_yyz_zzz_0, ts_yzz_xxx_0, ts_yzz_xxy_0, ts_yzz_xxz_0, ts_yzz_xyy_0, ts_yzz_xyz_0, \
                                     ts_yzz_xzz_0, ts_yzz_yyy_0, ts_yzz_yyz_0, ts_yzz_yzz_0, ts_yzz_zzz_0, ts_zzz_xxx_0, \
                                     ts_zzz_xxy_0, ts_zzz_xxz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdz_yyyz_yyy_0[j] = pa_y[j] * tdz_yyz_yyy_0[j] + fl1_fx * tdz_yz_yyy_0[j] + 1.5 * fl1_fx * tdz_yyz_yy_0[j];

                tdx_yyyz_yyz_0[j] = pa_y[j] * tdx_yyz_yyz_0[j] + fl1_fx * tdx_yz_yyz_0[j] + fl1_fx * tdx_yyz_yz_0[j];

                tdy_yyyz_yyz_0[j] = pa_y[j] * tdy_yyz_yyz_0[j] + fl1_fx * tdy_yz_yyz_0[j] + fl1_fx * tdy_yyz_yz_0[j] + 0.5 * fl1_fx * ts_yyz_yyz_0[j];

                tdz_yyyz_yyz_0[j] = pa_y[j] * tdz_yyz_yyz_0[j] + fl1_fx * tdz_yz_yyz_0[j] + fl1_fx * tdz_yyz_yz_0[j];

                tdx_yyyz_yzz_0[j] = pa_y[j] * tdx_yyz_yzz_0[j] + fl1_fx * tdx_yz_yzz_0[j] + 0.5 * fl1_fx * tdx_yyz_zz_0[j];

                tdy_yyyz_yzz_0[j] = pa_y[j] * tdy_yyz_yzz_0[j] + fl1_fx * tdy_yz_yzz_0[j] + 0.5 * fl1_fx * tdy_yyz_zz_0[j] + 0.5 * fl1_fx * ts_yyz_yzz_0[j];

                tdz_yyyz_yzz_0[j] = pa_y[j] * tdz_yyz_yzz_0[j] + fl1_fx * tdz_yz_yzz_0[j] + 0.5 * fl1_fx * tdz_yyz_zz_0[j];

                tdx_yyyz_zzz_0[j] = pa_y[j] * tdx_yyz_zzz_0[j] + fl1_fx * tdx_yz_zzz_0[j];

                tdy_yyyz_zzz_0[j] = pa_y[j] * tdy_yyz_zzz_0[j] + fl1_fx * tdy_yz_zzz_0[j] + 0.5 * fl1_fx * ts_yyz_zzz_0[j];

                tdz_yyyz_zzz_0[j] = pa_y[j] * tdz_yyz_zzz_0[j] + fl1_fx * tdz_yz_zzz_0[j];

                tdx_yyzz_xxx_0[j] = pa_y[j] * tdx_yzz_xxx_0[j] + 0.5 * fl1_fx * tdx_zz_xxx_0[j];

                tdy_yyzz_xxx_0[j] = pa_y[j] * tdy_yzz_xxx_0[j] + 0.5 * fl1_fx * tdy_zz_xxx_0[j] + 0.5 * fl1_fx * ts_yzz_xxx_0[j];

                tdz_yyzz_xxx_0[j] = pa_y[j] * tdz_yzz_xxx_0[j] + 0.5 * fl1_fx * tdz_zz_xxx_0[j];

                tdx_yyzz_xxy_0[j] = pa_y[j] * tdx_yzz_xxy_0[j] + 0.5 * fl1_fx * tdx_zz_xxy_0[j] + 0.5 * fl1_fx * tdx_yzz_xx_0[j];

                tdy_yyzz_xxy_0[j] = pa_y[j] * tdy_yzz_xxy_0[j] + 0.5 * fl1_fx * tdy_zz_xxy_0[j] + 0.5 * fl1_fx * tdy_yzz_xx_0[j] + 0.5 * fl1_fx * ts_yzz_xxy_0[j];

                tdz_yyzz_xxy_0[j] = pa_y[j] * tdz_yzz_xxy_0[j] + 0.5 * fl1_fx * tdz_zz_xxy_0[j] + 0.5 * fl1_fx * tdz_yzz_xx_0[j];

                tdx_yyzz_xxz_0[j] = pa_y[j] * tdx_yzz_xxz_0[j] + 0.5 * fl1_fx * tdx_zz_xxz_0[j];

                tdy_yyzz_xxz_0[j] = pa_y[j] * tdy_yzz_xxz_0[j] + 0.5 * fl1_fx * tdy_zz_xxz_0[j] + 0.5 * fl1_fx * ts_yzz_xxz_0[j];

                tdz_yyzz_xxz_0[j] = pa_y[j] * tdz_yzz_xxz_0[j] + 0.5 * fl1_fx * tdz_zz_xxz_0[j];

                tdx_yyzz_xyy_0[j] = pa_y[j] * tdx_yzz_xyy_0[j] + 0.5 * fl1_fx * tdx_zz_xyy_0[j] + fl1_fx * tdx_yzz_xy_0[j];

                tdy_yyzz_xyy_0[j] = pa_y[j] * tdy_yzz_xyy_0[j] + 0.5 * fl1_fx * tdy_zz_xyy_0[j] + fl1_fx * tdy_yzz_xy_0[j] + 0.5 * fl1_fx * ts_yzz_xyy_0[j];

                tdz_yyzz_xyy_0[j] = pa_y[j] * tdz_yzz_xyy_0[j] + 0.5 * fl1_fx * tdz_zz_xyy_0[j] + fl1_fx * tdz_yzz_xy_0[j];

                tdx_yyzz_xyz_0[j] = pa_y[j] * tdx_yzz_xyz_0[j] + 0.5 * fl1_fx * tdx_zz_xyz_0[j] + 0.5 * fl1_fx * tdx_yzz_xz_0[j];

                tdy_yyzz_xyz_0[j] = pa_y[j] * tdy_yzz_xyz_0[j] + 0.5 * fl1_fx * tdy_zz_xyz_0[j] + 0.5 * fl1_fx * tdy_yzz_xz_0[j] + 0.5 * fl1_fx * ts_yzz_xyz_0[j];

                tdz_yyzz_xyz_0[j] = pa_y[j] * tdz_yzz_xyz_0[j] + 0.5 * fl1_fx * tdz_zz_xyz_0[j] + 0.5 * fl1_fx * tdz_yzz_xz_0[j];

                tdx_yyzz_xzz_0[j] = pa_y[j] * tdx_yzz_xzz_0[j] + 0.5 * fl1_fx * tdx_zz_xzz_0[j];

                tdy_yyzz_xzz_0[j] = pa_y[j] * tdy_yzz_xzz_0[j] + 0.5 * fl1_fx * tdy_zz_xzz_0[j] + 0.5 * fl1_fx * ts_yzz_xzz_0[j];

                tdz_yyzz_xzz_0[j] = pa_y[j] * tdz_yzz_xzz_0[j] + 0.5 * fl1_fx * tdz_zz_xzz_0[j];

                tdx_yyzz_yyy_0[j] = pa_y[j] * tdx_yzz_yyy_0[j] + 0.5 * fl1_fx * tdx_zz_yyy_0[j] + 1.5 * fl1_fx * tdx_yzz_yy_0[j];

                tdy_yyzz_yyy_0[j] = pa_y[j] * tdy_yzz_yyy_0[j] + 0.5 * fl1_fx * tdy_zz_yyy_0[j] + 1.5 * fl1_fx * tdy_yzz_yy_0[j] + 0.5 * fl1_fx * ts_yzz_yyy_0[j];

                tdz_yyzz_yyy_0[j] = pa_y[j] * tdz_yzz_yyy_0[j] + 0.5 * fl1_fx * tdz_zz_yyy_0[j] + 1.5 * fl1_fx * tdz_yzz_yy_0[j];

                tdx_yyzz_yyz_0[j] = pa_y[j] * tdx_yzz_yyz_0[j] + 0.5 * fl1_fx * tdx_zz_yyz_0[j] + fl1_fx * tdx_yzz_yz_0[j];

                tdy_yyzz_yyz_0[j] = pa_y[j] * tdy_yzz_yyz_0[j] + 0.5 * fl1_fx * tdy_zz_yyz_0[j] + fl1_fx * tdy_yzz_yz_0[j] + 0.5 * fl1_fx * ts_yzz_yyz_0[j];

                tdz_yyzz_yyz_0[j] = pa_y[j] * tdz_yzz_yyz_0[j] + 0.5 * fl1_fx * tdz_zz_yyz_0[j] + fl1_fx * tdz_yzz_yz_0[j];

                tdx_yyzz_yzz_0[j] = pa_y[j] * tdx_yzz_yzz_0[j] + 0.5 * fl1_fx * tdx_zz_yzz_0[j] + 0.5 * fl1_fx * tdx_yzz_zz_0[j];

                tdy_yyzz_yzz_0[j] = pa_y[j] * tdy_yzz_yzz_0[j] + 0.5 * fl1_fx * tdy_zz_yzz_0[j] + 0.5 * fl1_fx * tdy_yzz_zz_0[j] + 0.5 * fl1_fx * ts_yzz_yzz_0[j];

                tdz_yyzz_yzz_0[j] = pa_y[j] * tdz_yzz_yzz_0[j] + 0.5 * fl1_fx * tdz_zz_yzz_0[j] + 0.5 * fl1_fx * tdz_yzz_zz_0[j];

                tdx_yyzz_zzz_0[j] = pa_y[j] * tdx_yzz_zzz_0[j] + 0.5 * fl1_fx * tdx_zz_zzz_0[j];

                tdy_yyzz_zzz_0[j] = pa_y[j] * tdy_yzz_zzz_0[j] + 0.5 * fl1_fx * tdy_zz_zzz_0[j] + 0.5 * fl1_fx * ts_yzz_zzz_0[j];

                tdz_yyzz_zzz_0[j] = pa_y[j] * tdz_yzz_zzz_0[j] + 0.5 * fl1_fx * tdz_zz_zzz_0[j];

                tdx_yzzz_xxx_0[j] = pa_y[j] * tdx_zzz_xxx_0[j];

                tdy_yzzz_xxx_0[j] = pa_y[j] * tdy_zzz_xxx_0[j] + 0.5 * fl1_fx * ts_zzz_xxx_0[j];

                tdz_yzzz_xxx_0[j] = pa_y[j] * tdz_zzz_xxx_0[j];

                tdx_yzzz_xxy_0[j] = pa_y[j] * tdx_zzz_xxy_0[j] + 0.5 * fl1_fx * tdx_zzz_xx_0[j];

                tdy_yzzz_xxy_0[j] = pa_y[j] * tdy_zzz_xxy_0[j] + 0.5 * fl1_fx * tdy_zzz_xx_0[j] + 0.5 * fl1_fx * ts_zzz_xxy_0[j];

                tdz_yzzz_xxy_0[j] = pa_y[j] * tdz_zzz_xxy_0[j] + 0.5 * fl1_fx * tdz_zzz_xx_0[j];

                tdx_yzzz_xxz_0[j] = pa_y[j] * tdx_zzz_xxz_0[j];

                tdy_yzzz_xxz_0[j] = pa_y[j] * tdy_zzz_xxz_0[j] + 0.5 * fl1_fx * ts_zzz_xxz_0[j];

                tdz_yzzz_xxz_0[j] = pa_y[j] * tdz_zzz_xxz_0[j];

                tdx_yzzz_xyy_0[j] = pa_y[j] * tdx_zzz_xyy_0[j] + fl1_fx * tdx_zzz_xy_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGF_400_450(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto)
    {
        // Batch of Integrals (400,450)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 90); 

            auto tdy_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 90); 

            auto tdz_zzz_xxx_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 90); 

            auto tdx_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 91); 

            auto tdy_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 91); 

            auto tdz_zzz_xxy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 91); 

            auto tdx_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 92); 

            auto tdy_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 92); 

            auto tdz_zzz_xxz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 92); 

            auto tdx_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 93); 

            auto tdy_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 93); 

            auto tdz_zzz_xyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 93); 

            auto tdx_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 94); 

            auto tdy_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 94); 

            auto tdz_zzz_xyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 94); 

            auto tdx_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 95); 

            auto tdy_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 95); 

            auto tdz_zzz_xzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 95); 

            auto tdx_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 96); 

            auto tdy_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 96); 

            auto tdz_zzz_yyy_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 96); 

            auto tdx_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 97); 

            auto tdy_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 97); 

            auto tdz_zzz_yyz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 97); 

            auto tdx_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 98); 

            auto tdy_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 98); 

            auto tdz_zzz_yzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 98); 

            auto tdx_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * idx + 99); 

            auto tdy_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 100 * bdim + 100 * idx + 99); 

            auto tdz_zzz_zzz_0 = primBuffer.data(pidx_d_3_3_m0 + 200 * bdim + 100 * idx + 99); 

            auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50); 

            auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51); 

            auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52); 

            auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 53); 

            auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53); 

            auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53); 

            auto tdx_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 54); 

            auto tdy_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdx_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 55); 

            auto tdy_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 55); 

            auto tdz_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 55); 

            auto tdx_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 56); 

            auto tdy_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 56); 

            auto tdz_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 56); 

            auto tdx_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 57); 

            auto tdy_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 57); 

            auto tdz_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 57); 

            auto tdx_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 58); 

            auto tdy_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 58); 

            auto tdz_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 58); 

            auto tdx_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 59); 

            auto tdy_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 59); 

            auto tdz_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 59); 

            auto tdx_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 54); 

            auto tdy_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdx_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 55); 

            auto tdy_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tdz_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tdx_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 56); 

            auto tdy_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tdz_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tdx_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 57); 

            auto tdy_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tdz_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tdx_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 58); 

            auto tdy_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tdz_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tdx_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 59); 

            auto tdy_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tdz_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            auto ts_zzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 90); 

            auto ts_zzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 91); 

            auto ts_zzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 92); 

            auto ts_zzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 93); 

            auto ts_zzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 94); 

            auto ts_zzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 95); 

            auto ts_zzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 96); 

            auto ts_zzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 97); 

            auto ts_zzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 98); 

            auto ts_zzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 99); 

            // set up pointers to integrals

            auto tdy_yzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 133); 

            auto tdz_yzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 133); 

            auto tdx_yzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 134); 

            auto tdy_yzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 134); 

            auto tdz_yzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 134); 

            auto tdx_yzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 135); 

            auto tdy_yzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 135); 

            auto tdz_yzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 135); 

            auto tdx_yzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 136); 

            auto tdy_yzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 136); 

            auto tdz_yzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 136); 

            auto tdx_yzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 137); 

            auto tdy_yzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 137); 

            auto tdz_yzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 137); 

            auto tdx_yzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 138); 

            auto tdy_yzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 138); 

            auto tdz_yzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 138); 

            auto tdx_yzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 139); 

            auto tdy_yzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 139); 

            auto tdz_yzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 139); 

            auto tdx_zzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 140); 

            auto tdy_zzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 140); 

            auto tdz_zzzz_xxx_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 140); 

            auto tdx_zzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 141); 

            auto tdy_zzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 141); 

            auto tdz_zzzz_xxy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 141); 

            auto tdx_zzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 142); 

            auto tdy_zzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 142); 

            auto tdz_zzzz_xxz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 142); 

            auto tdx_zzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 143); 

            auto tdy_zzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 143); 

            auto tdz_zzzz_xyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 143); 

            auto tdx_zzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 144); 

            auto tdy_zzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 144); 

            auto tdz_zzzz_xyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 144); 

            auto tdx_zzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 145); 

            auto tdy_zzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 145); 

            auto tdz_zzzz_xzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 145); 

            auto tdx_zzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 146); 

            auto tdy_zzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 146); 

            auto tdz_zzzz_yyy_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 146); 

            auto tdx_zzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 147); 

            auto tdy_zzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 147); 

            auto tdz_zzzz_yyz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 147); 

            auto tdx_zzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 148); 

            auto tdy_zzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 148); 

            auto tdz_zzzz_yzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 148); 

            auto tdx_zzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * idx + 149); 

            auto tdy_zzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 150 * bdim + 150 * idx + 149); 

            auto tdz_zzzz_zzz_0 = primBuffer.data(pidx_d_4_3_m0 + 300 * bdim + 150 * idx + 149); 

            // Batch of Integrals (400,450)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_yzzz_xyz_0, tdx_yzzz_xzz_0, tdx_yzzz_yyy_0, \
                                     tdx_yzzz_yyz_0, tdx_yzzz_yzz_0, tdx_yzzz_zzz_0, tdx_zz_xxx_0, tdx_zz_xxy_0, \
                                     tdx_zz_xxz_0, tdx_zz_xyy_0, tdx_zz_xyz_0, tdx_zz_xzz_0, tdx_zz_yyy_0, tdx_zz_yyz_0, \
                                     tdx_zz_yzz_0, tdx_zz_zzz_0, tdx_zzz_xx_0, tdx_zzz_xxx_0, tdx_zzz_xxy_0, \
                                     tdx_zzz_xxz_0, tdx_zzz_xy_0, tdx_zzz_xyy_0, tdx_zzz_xyz_0, tdx_zzz_xz_0, \
                                     tdx_zzz_xzz_0, tdx_zzz_yy_0, tdx_zzz_yyy_0, tdx_zzz_yyz_0, tdx_zzz_yz_0, \
                                     tdx_zzz_yzz_0, tdx_zzz_zz_0, tdx_zzz_zzz_0, tdx_zzzz_xxx_0, tdx_zzzz_xxy_0, \
                                     tdx_zzzz_xxz_0, tdx_zzzz_xyy_0, tdx_zzzz_xyz_0, tdx_zzzz_xzz_0, tdx_zzzz_yyy_0, \
                                     tdx_zzzz_yyz_0, tdx_zzzz_yzz_0, tdx_zzzz_zzz_0, tdy_yzzz_xyy_0, tdy_yzzz_xyz_0, \
                                     tdy_yzzz_xzz_0, tdy_yzzz_yyy_0, tdy_yzzz_yyz_0, tdy_yzzz_yzz_0, tdy_yzzz_zzz_0, \
                                     tdy_zz_xxx_0, tdy_zz_xxy_0, tdy_zz_xxz_0, tdy_zz_xyy_0, tdy_zz_xyz_0, tdy_zz_xzz_0, \
                                     tdy_zz_yyy_0, tdy_zz_yyz_0, tdy_zz_yzz_0, tdy_zz_zzz_0, tdy_zzz_xx_0, \
                                     tdy_zzz_xxx_0, tdy_zzz_xxy_0, tdy_zzz_xxz_0, tdy_zzz_xy_0, tdy_zzz_xyy_0, \
                                     tdy_zzz_xyz_0, tdy_zzz_xz_0, tdy_zzz_xzz_0, tdy_zzz_yy_0, tdy_zzz_yyy_0, \
                                     tdy_zzz_yyz_0, tdy_zzz_yz_0, tdy_zzz_yzz_0, tdy_zzz_zz_0, tdy_zzz_zzz_0, \
                                     tdy_zzzz_xxx_0, tdy_zzzz_xxy_0, tdy_zzzz_xxz_0, tdy_zzzz_xyy_0, tdy_zzzz_xyz_0, \
                                     tdy_zzzz_xzz_0, tdy_zzzz_yyy_0, tdy_zzzz_yyz_0, tdy_zzzz_yzz_0, tdy_zzzz_zzz_0, \
                                     tdz_yzzz_xyy_0, tdz_yzzz_xyz_0, tdz_yzzz_xzz_0, tdz_yzzz_yyy_0, tdz_yzzz_yyz_0, \
                                     tdz_yzzz_yzz_0, tdz_yzzz_zzz_0, tdz_zz_xxx_0, tdz_zz_xxy_0, tdz_zz_xxz_0, \
                                     tdz_zz_xyy_0, tdz_zz_xyz_0, tdz_zz_xzz_0, tdz_zz_yyy_0, tdz_zz_yyz_0, tdz_zz_yzz_0, \
                                     tdz_zz_zzz_0, tdz_zzz_xx_0, tdz_zzz_xxx_0, tdz_zzz_xxy_0, tdz_zzz_xxz_0, \
                                     tdz_zzz_xy_0, tdz_zzz_xyy_0, tdz_zzz_xyz_0, tdz_zzz_xz_0, tdz_zzz_xzz_0, \
                                     tdz_zzz_yy_0, tdz_zzz_yyy_0, tdz_zzz_yyz_0, tdz_zzz_yz_0, tdz_zzz_yzz_0, \
                                     tdz_zzz_zz_0, tdz_zzz_zzz_0, tdz_zzzz_xxx_0, tdz_zzzz_xxy_0, tdz_zzzz_xxz_0, \
                                     tdz_zzzz_xyy_0, tdz_zzzz_xyz_0, tdz_zzzz_xzz_0, tdz_zzzz_yyy_0, tdz_zzzz_yyz_0, \
                                     tdz_zzzz_yzz_0, tdz_zzzz_zzz_0, ts_zzz_xxx_0, ts_zzz_xxy_0, ts_zzz_xxz_0, \
                                     ts_zzz_xyy_0, ts_zzz_xyz_0, ts_zzz_xzz_0, ts_zzz_yyy_0, ts_zzz_yyz_0, ts_zzz_yzz_0, \
                                     ts_zzz_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdy_yzzz_xyy_0[j] = pa_y[j] * tdy_zzz_xyy_0[j] + fl1_fx * tdy_zzz_xy_0[j] + 0.5 * fl1_fx * ts_zzz_xyy_0[j];

                tdz_yzzz_xyy_0[j] = pa_y[j] * tdz_zzz_xyy_0[j] + fl1_fx * tdz_zzz_xy_0[j];

                tdx_yzzz_xyz_0[j] = pa_y[j] * tdx_zzz_xyz_0[j] + 0.5 * fl1_fx * tdx_zzz_xz_0[j];

                tdy_yzzz_xyz_0[j] = pa_y[j] * tdy_zzz_xyz_0[j] + 0.5 * fl1_fx * tdy_zzz_xz_0[j] + 0.5 * fl1_fx * ts_zzz_xyz_0[j];

                tdz_yzzz_xyz_0[j] = pa_y[j] * tdz_zzz_xyz_0[j] + 0.5 * fl1_fx * tdz_zzz_xz_0[j];

                tdx_yzzz_xzz_0[j] = pa_y[j] * tdx_zzz_xzz_0[j];

                tdy_yzzz_xzz_0[j] = pa_y[j] * tdy_zzz_xzz_0[j] + 0.5 * fl1_fx * ts_zzz_xzz_0[j];

                tdz_yzzz_xzz_0[j] = pa_y[j] * tdz_zzz_xzz_0[j];

                tdx_yzzz_yyy_0[j] = pa_y[j] * tdx_zzz_yyy_0[j] + 1.5 * fl1_fx * tdx_zzz_yy_0[j];

                tdy_yzzz_yyy_0[j] = pa_y[j] * tdy_zzz_yyy_0[j] + 1.5 * fl1_fx * tdy_zzz_yy_0[j] + 0.5 * fl1_fx * ts_zzz_yyy_0[j];

                tdz_yzzz_yyy_0[j] = pa_y[j] * tdz_zzz_yyy_0[j] + 1.5 * fl1_fx * tdz_zzz_yy_0[j];

                tdx_yzzz_yyz_0[j] = pa_y[j] * tdx_zzz_yyz_0[j] + fl1_fx * tdx_zzz_yz_0[j];

                tdy_yzzz_yyz_0[j] = pa_y[j] * tdy_zzz_yyz_0[j] + fl1_fx * tdy_zzz_yz_0[j] + 0.5 * fl1_fx * ts_zzz_yyz_0[j];

                tdz_yzzz_yyz_0[j] = pa_y[j] * tdz_zzz_yyz_0[j] + fl1_fx * tdz_zzz_yz_0[j];

                tdx_yzzz_yzz_0[j] = pa_y[j] * tdx_zzz_yzz_0[j] + 0.5 * fl1_fx * tdx_zzz_zz_0[j];

                tdy_yzzz_yzz_0[j] = pa_y[j] * tdy_zzz_yzz_0[j] + 0.5 * fl1_fx * tdy_zzz_zz_0[j] + 0.5 * fl1_fx * ts_zzz_yzz_0[j];

                tdz_yzzz_yzz_0[j] = pa_y[j] * tdz_zzz_yzz_0[j] + 0.5 * fl1_fx * tdz_zzz_zz_0[j];

                tdx_yzzz_zzz_0[j] = pa_y[j] * tdx_zzz_zzz_0[j];

                tdy_yzzz_zzz_0[j] = pa_y[j] * tdy_zzz_zzz_0[j] + 0.5 * fl1_fx * ts_zzz_zzz_0[j];

                tdz_yzzz_zzz_0[j] = pa_y[j] * tdz_zzz_zzz_0[j];

                tdx_zzzz_xxx_0[j] = pa_z[j] * tdx_zzz_xxx_0[j] + 1.5 * fl1_fx * tdx_zz_xxx_0[j];

                tdy_zzzz_xxx_0[j] = pa_z[j] * tdy_zzz_xxx_0[j] + 1.5 * fl1_fx * tdy_zz_xxx_0[j];

                tdz_zzzz_xxx_0[j] = pa_z[j] * tdz_zzz_xxx_0[j] + 1.5 * fl1_fx * tdz_zz_xxx_0[j] + 0.5 * fl1_fx * ts_zzz_xxx_0[j];

                tdx_zzzz_xxy_0[j] = pa_z[j] * tdx_zzz_xxy_0[j] + 1.5 * fl1_fx * tdx_zz_xxy_0[j];

                tdy_zzzz_xxy_0[j] = pa_z[j] * tdy_zzz_xxy_0[j] + 1.5 * fl1_fx * tdy_zz_xxy_0[j];

                tdz_zzzz_xxy_0[j] = pa_z[j] * tdz_zzz_xxy_0[j] + 1.5 * fl1_fx * tdz_zz_xxy_0[j] + 0.5 * fl1_fx * ts_zzz_xxy_0[j];

                tdx_zzzz_xxz_0[j] = pa_z[j] * tdx_zzz_xxz_0[j] + 1.5 * fl1_fx * tdx_zz_xxz_0[j] + 0.5 * fl1_fx * tdx_zzz_xx_0[j];

                tdy_zzzz_xxz_0[j] = pa_z[j] * tdy_zzz_xxz_0[j] + 1.5 * fl1_fx * tdy_zz_xxz_0[j] + 0.5 * fl1_fx * tdy_zzz_xx_0[j];

                tdz_zzzz_xxz_0[j] = pa_z[j] * tdz_zzz_xxz_0[j] + 1.5 * fl1_fx * tdz_zz_xxz_0[j] + 0.5 * fl1_fx * tdz_zzz_xx_0[j] + 0.5 * fl1_fx * ts_zzz_xxz_0[j];

                tdx_zzzz_xyy_0[j] = pa_z[j] * tdx_zzz_xyy_0[j] + 1.5 * fl1_fx * tdx_zz_xyy_0[j];

                tdy_zzzz_xyy_0[j] = pa_z[j] * tdy_zzz_xyy_0[j] + 1.5 * fl1_fx * tdy_zz_xyy_0[j];

                tdz_zzzz_xyy_0[j] = pa_z[j] * tdz_zzz_xyy_0[j] + 1.5 * fl1_fx * tdz_zz_xyy_0[j] + 0.5 * fl1_fx * ts_zzz_xyy_0[j];

                tdx_zzzz_xyz_0[j] = pa_z[j] * tdx_zzz_xyz_0[j] + 1.5 * fl1_fx * tdx_zz_xyz_0[j] + 0.5 * fl1_fx * tdx_zzz_xy_0[j];

                tdy_zzzz_xyz_0[j] = pa_z[j] * tdy_zzz_xyz_0[j] + 1.5 * fl1_fx * tdy_zz_xyz_0[j] + 0.5 * fl1_fx * tdy_zzz_xy_0[j];

                tdz_zzzz_xyz_0[j] = pa_z[j] * tdz_zzz_xyz_0[j] + 1.5 * fl1_fx * tdz_zz_xyz_0[j] + 0.5 * fl1_fx * tdz_zzz_xy_0[j] + 0.5 * fl1_fx * ts_zzz_xyz_0[j];

                tdx_zzzz_xzz_0[j] = pa_z[j] * tdx_zzz_xzz_0[j] + 1.5 * fl1_fx * tdx_zz_xzz_0[j] + fl1_fx * tdx_zzz_xz_0[j];

                tdy_zzzz_xzz_0[j] = pa_z[j] * tdy_zzz_xzz_0[j] + 1.5 * fl1_fx * tdy_zz_xzz_0[j] + fl1_fx * tdy_zzz_xz_0[j];

                tdz_zzzz_xzz_0[j] = pa_z[j] * tdz_zzz_xzz_0[j] + 1.5 * fl1_fx * tdz_zz_xzz_0[j] + fl1_fx * tdz_zzz_xz_0[j] + 0.5 * fl1_fx * ts_zzz_xzz_0[j];

                tdx_zzzz_yyy_0[j] = pa_z[j] * tdx_zzz_yyy_0[j] + 1.5 * fl1_fx * tdx_zz_yyy_0[j];

                tdy_zzzz_yyy_0[j] = pa_z[j] * tdy_zzz_yyy_0[j] + 1.5 * fl1_fx * tdy_zz_yyy_0[j];

                tdz_zzzz_yyy_0[j] = pa_z[j] * tdz_zzz_yyy_0[j] + 1.5 * fl1_fx * tdz_zz_yyy_0[j] + 0.5 * fl1_fx * ts_zzz_yyy_0[j];

                tdx_zzzz_yyz_0[j] = pa_z[j] * tdx_zzz_yyz_0[j] + 1.5 * fl1_fx * tdx_zz_yyz_0[j] + 0.5 * fl1_fx * tdx_zzz_yy_0[j];

                tdy_zzzz_yyz_0[j] = pa_z[j] * tdy_zzz_yyz_0[j] + 1.5 * fl1_fx * tdy_zz_yyz_0[j] + 0.5 * fl1_fx * tdy_zzz_yy_0[j];

                tdz_zzzz_yyz_0[j] = pa_z[j] * tdz_zzz_yyz_0[j] + 1.5 * fl1_fx * tdz_zz_yyz_0[j] + 0.5 * fl1_fx * tdz_zzz_yy_0[j] + 0.5 * fl1_fx * ts_zzz_yyz_0[j];

                tdx_zzzz_yzz_0[j] = pa_z[j] * tdx_zzz_yzz_0[j] + 1.5 * fl1_fx * tdx_zz_yzz_0[j] + fl1_fx * tdx_zzz_yz_0[j];

                tdy_zzzz_yzz_0[j] = pa_z[j] * tdy_zzz_yzz_0[j] + 1.5 * fl1_fx * tdy_zz_yzz_0[j] + fl1_fx * tdy_zzz_yz_0[j];

                tdz_zzzz_yzz_0[j] = pa_z[j] * tdz_zzz_yzz_0[j] + 1.5 * fl1_fx * tdz_zz_yzz_0[j] + fl1_fx * tdz_zzz_yz_0[j] + 0.5 * fl1_fx * ts_zzz_yzz_0[j];

                tdx_zzzz_zzz_0[j] = pa_z[j] * tdx_zzz_zzz_0[j] + 1.5 * fl1_fx * tdx_zz_zzz_0[j] + 1.5 * fl1_fx * tdx_zzz_zz_0[j];

                tdy_zzzz_zzz_0[j] = pa_z[j] * tdy_zzz_zzz_0[j] + 1.5 * fl1_fx * tdy_zz_zzz_0[j] + 1.5 * fl1_fx * tdy_zzz_zz_0[j];

                tdz_zzzz_zzz_0[j] = pa_z[j] * tdz_zzz_zzz_0[j] + 1.5 * fl1_fx * tdz_zz_zzz_0[j] + 1.5 * fl1_fx * tdz_zzz_zz_0[j] + 0.5 * fl1_fx * ts_zzz_zzz_0[j];
            }

            idx++;
        }
    }


} // ediprecfunc namespace

