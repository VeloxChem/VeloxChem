//                                                                              
//                           VELOXCHEM 1.0-RC2                                  
//         ----------------------------------------------------                 
//                     An Electronic Structure Code                             
//                                                                              
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.         
//  Contact: https://veloxchem.org/contact                                      
//                                                                              
//  SPDX-License-Identifier: LGPL-3.0-or-later                                  
//                                                                              
//  This file is part of VeloxChem.                                             
//                                                                              
//  VeloxChem is free software: you can redistribute it and/or modify it under  
//  the terms of the GNU Lesser General Public License as published by the Free 
//  Software Foundation, either version 3 of the License, or (at your option)   
//  any later version.                                                          
//                                                                              
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT    
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public         
//  License for more details.                                                   
//                                                                              
//  You should have received a copy of the GNU Lesser General Public License    
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.          

#include <cstdint>

#include "Buffer.hpp"
#include "BinnedGtoPairBlock.hpp"
#include "DiagEriRecFacts.hpp"
#include "EriVRRForSSSS.hpp"
#include "EriVRRForSSSP.hpp"
#include "EriVRRForSSSD.hpp"
#include "EriVRRForSSSF.hpp"
#include "EriVRRForSPSD.hpp"
#include "EriVRRForSPSF.hpp"
#include "EriHRRForSPPD.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostSPPD(      T*                                 intsBuffer,
             const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
             const int32_t                            bPosition,
             const int32_t                            ePosition) -> void
{
    // set up dimensions

    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();

    const auto ncpairs = ePosition - bPosition;

    // allocate Obara-Saika factors

    BufferHostMY<T, 4> osfacts(ncpairs);

    // allocate distances

    BufferHostMY<T, 3> rpq(ncpairs); 

    BufferHostMY<T, 3> rcd(ncpairs);

    BufferHostMY<T, 3> rpb(ncpairs);

    BufferHostMY<T, 3> rqd(ncpairs);

    BufferHostMY<T, 3> rwp(ncpairs);

    BufferHostMY<T, 3> rwq(ncpairs);

    // allocate coordinates

    BufferHostMY<T, 3> rw(ncpairs); 

    // allocate Boys function data

    BufferHostX<T> bargs(ncpairs);

    BufferHostXY<T> bvals(5, ncpairs);

    CBoysFunc<T, 4> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(5, ncpairs);

    BufferHostXY<T> pbufSSSP0(3, ncpairs);

0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 
    BufferHostXY<T> pbufSSSP1(3, ncpairs);

0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 
    BufferHostXY<T> pbufSSSP2(3, ncpairs);

0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 
    BufferHostXY<T> pbufSSSP3(3, ncpairs);

0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 
    BufferHostXY<T> pbufSSSD0(6, ncpairs);

0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 
    BufferHostXY<T> pbufSSSD1(6, ncpairs);

0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 
    BufferHostXY<T> pbufSSSD2(4, ncpairs);

0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xx_2 
    BufferHostXY<T> pbufSSSF0(10, ncpairs);

0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 
    BufferHostXY<T> pbufSSSF1(10, ncpairs);

0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 
    BufferHostXY<T> pbufSPSD0(18, ncpairs);

0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 
    BufferHostXY<T> pbufSPSF0(30, ncpairs);

0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSPSD(18, ncpairs);

0_z_0_zz 0_z_0_yz 0_z_0_yy 0_z_0_xz 0_z_0_xy 0_z_0_xx 0_y_0_zz 0_y_0_yz 0_y_0_yy 0_y_0_xz 0_y_0_xy 0_y_0_xx 0_x_0_zz 0_x_0_yz 0_x_0_yy 0_x_0_xz 0_x_0_xy 0_x_0_xx 
    BufferHostXY<T> cbufSPSF(30, ncpairs);

0_z_0_zzz 0_z_0_yzz 0_z_0_yyz 0_z_0_yyy 0_z_0_xzz 0_z_0_xyz 0_z_0_xyy 0_z_0_xxz 0_z_0_xxy 0_z_0_xxx 0_y_0_zzz 0_y_0_yzz 0_y_0_yyz 0_y_0_yyy 0_y_0_xzz 0_y_0_xyz 0_y_0_xyy 0_y_0_xxz 0_y_0_xxy 0_y_0_xxx 0_x_0_zzz 0_x_0_yzz 0_x_0_yyz 0_x_0_yyy 0_x_0_xzz 0_x_0_xyz 0_x_0_xyy 0_x_0_xxz 0_x_0_xxy 0_x_0_xxx 
    BufferHostXY<T> cbufSPPD(54, ncpairs);

0_z_z_zz 0_z_z_yz 0_z_z_yy 0_z_z_xz 0_z_z_xy 0_z_z_xx 0_z_y_zz 0_z_y_yz 0_z_y_yy 0_z_y_xz 0_z_y_xy 0_z_y_xx 0_z_x_zz 0_z_x_yz 0_z_x_yy 0_z_x_xz 0_z_x_xy 0_z_x_xx 0_y_z_zz 0_y_z_yz 0_y_z_yy 0_y_z_xz 0_y_z_xy 0_y_z_xx 0_y_y_zz 0_y_y_yz 0_y_y_yy 0_y_y_xz 0_y_y_xy 0_y_y_xx 0_y_x_zz 0_y_x_yz 0_y_x_yy 0_y_x_xz 0_y_x_xy 0_y_x_xx 0_x_z_zz 0_x_z_yz 0_x_z_yy 0_x_z_xz 0_x_z_xy 0_x_z_xx 0_x_y_zz 0_x_y_yz 0_x_y_yy 0_x_y_xz 0_x_y_xy 0_x_y_xx 0_x_x_zz 0_x_x_yz 0_x_x_yy 0_x_x_xz 0_x_x_xy 0_x_x_xx 
    for (int32_t i = 0; i < nppairs; i++)
    {
        derirec::compHostDistancesPT(rpb, gtoPairBlock, bPosition, ePosition, i);

        for (int j = i; j < nppairs; j++)
        {
            derirec::compHostDistancesPQ(rpq, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsRho(frho, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsNorm(fnorm, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsZeta(fze, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsPartialZeta(fe, gtoPairBlock, bPosition, ePosition, j);

            derirec::compHostDistancesPT(rqd, gtoPairBlock, bPosition, ePosition, j);

            derirec::compHostCoordinatesW(rw, gtoPairBlock, bPosition, ePosition, i, j);

            if (i == j)
            {
                rwp.setZero();

                rwq.setZero();
            }
            else
            {
                derirec::compHostDistancesWT(rwp, rw, gtoPairBlock, bPosition, ePosition, i);

                derirec::compHostDistancesWT(rwq, rw, gtoPairBlock, bPosition, ePosition, j);
            }

            derirec::compHostBoysArguments(bargs, rpq, frho, ncpairs);

            bftable.compute(bvals, bargs);

            derirec::compHostVRRForSSSS(pbufSSSS, bvals, frho, fnorm, ncpairs);

signature:
0_z_0_zz_0 0_z_0_zzz_0 0_z_0_yz_0 0_z_0_yzz_0 0_z_0_yy_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xz_0 0_z_0_xzz_0 0_z_0_xy_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xx_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_z_z_zz_0 0_z_z_yz_0 0_z_z_yy_0 0_z_z_xz_0 0_z_z_xy_0 0_z_z_xx_0 0_z_y_zz_0 0_z_y_yz_0 0_z_y_yy_0 0_z_y_xz_0 0_z_y_xy_0 0_z_y_xx_0 0_z_x_zz_0 0_z_x_yz_0 0_z_x_yy_0 0_z_x_xz_0 0_z_x_xy_0 0_z_x_xx_0 0_y_0_zz_0 0_y_0_zzz_0 0_y_0_yz_0 0_y_0_yzz_0 0_y_0_yy_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xz_0 0_y_0_xzz_0 0_y_0_xy_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xx_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_y_z_zz_0 0_y_z_yz_0 0_y_z_yy_0 0_y_z_xz_0 0_y_z_xy_0 0_y_z_xx_0 0_y_y_zz_0 0_y_y_yz_0 0_y_y_yy_0 0_y_y_xz_0 0_y_y_xy_0 0_y_y_xx_0 0_y_x_zz_0 0_y_x_yz_0 0_y_x_yy_0 0_y_x_xz_0 0_y_x_xy_0 0_y_x_xx_0 0_x_0_zz_0 0_x_0_zzz_0 0_x_0_yz_0 0_x_0_yzz_0 0_x_0_yy_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xz_0 0_x_0_xzz_0 0_x_0_xy_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xx_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 0_x_z_zz_0 0_x_z_yz_0 0_x_z_yy_0 0_x_z_xz_0 0_x_z_xy_0 0_x_z_xx_0 0_x_y_zz_0 0_x_y_yz_0 0_x_y_yy_0 0_x_y_xz_0 0_x_y_xy_0 0_x_y_xx_0 0_x_x_zz_0 0_x_x_yz_0 0_x_x_yy_0 0_x_x_xz_0 0_x_x_xy_0 0_x_x_xx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xx_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

signature:
0_0_0_0_0 

signature:
0_0_0_0_0 

signature:
0_0_0_0_0 

signature:
0_0_0_0_0 

signature:
0_0_0_0_0 

INTEGRAL:0 : 1 : 0 : SSSS_0 Y SSSS_1 Y SSSP_0 Y 

0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSS_0 SSSS_1 SSSP_0 

SSSS_0 : 0_0_0_0_0 

SSSS_1 : 0_0_0_0_1 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

INTEGRAL:0 : 1 : 1 : SSSS_1 Y SSSS_2 Y SSSP_1 Y 

0_0_0_0_1 0_0_0_0_2 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSS_1 SSSS_2 SSSP_1 

SSSS_1 : 0_0_0_0_1 

SSSS_2 : 0_0_0_0_2 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

INTEGRAL:0 : 1 : 2 : SSSS_2 Y SSSS_3 Y SSSP_2 Y 

0_0_0_0_2 0_0_0_0_3 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSS_2 SSSS_3 SSSP_2 

SSSS_2 : 0_0_0_0_2 

SSSS_3 : 0_0_0_0_3 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

INTEGRAL:0 : 1 : 3 : SSSS_3 Y SSSS_4 Y SSSP_3 Y 

0_0_0_0_3 0_0_0_0_4 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSS_3 SSSS_4 SSSP_3 

SSSS_3 : 0_0_0_0_3 

SSSS_4 : 0_0_0_0_4 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

INTEGRAL:0 : 2 : 0 : SSSS_0 Y SSSS_1 Y SSSP_0 Y SSSP_1 Y SSSD_0 Y 

0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSS_0 SSSS_1 SSSP_0 SSSP_1 SSSD_0 

SSSS_0 : 0_0_0_0_0 

SSSS_1 : 0_0_0_0_1 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

INTEGRAL:0 : 2 : 1 : SSSS_1 Y SSSS_2 Y SSSP_1 Y SSSP_2 Y SSSD_1 Y 

0_0_0_0_1 0_0_0_0_2 0_0_0_z_1 0_0_0_z_2 0_0_0_zz_1 0_0_0_y_1 0_0_0_y_2 0_0_0_yz_1 0_0_0_yy_1 0_0_0_x_1 0_0_0_x_2 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSS_1 SSSS_2 SSSP_1 SSSP_2 SSSD_1 

SSSS_1 : 0_0_0_0_1 

SSSS_2 : 0_0_0_0_2 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

INTEGRAL:0 : 2 : 2 : SSSS_2 Y SSSS_3 Y SSSP_2 Y SSSP_3 Y SSSD_2 Y 

0_0_0_0_2 0_0_0_0_3 0_0_0_z_2 0_0_0_z_3 0_0_0_zz_2 0_0_0_y_2 0_0_0_y_3 0_0_0_yz_2 0_0_0_yy_2 0_0_0_x_2 0_0_0_x_3 0_0_0_xx_2 

SSSS_2 SSSS_3 SSSP_2 SSSP_3 SSSD_2 

SSSS_2 : 0_0_0_0_2 

SSSS_3 : 0_0_0_0_3 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xx_2 

INTEGRAL:0 : 3 : 0 : SSSP_0 Y SSSP_1 Y SSSD_0 N SSSD_1 N SSSF_0 Y 

0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSP_0 SSSP_1 SSSD_0 SSSD_1 SSSF_0 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

INTEGRAL:0 : 3 : 1 : SSSP_1 Y SSSP_2 Y SSSD_1 N SSSD_2 Y SSSF_1 Y 

0_0_0_z_1 0_0_0_z_2 0_0_0_zz_1 0_0_0_zz_2 0_0_0_zzz_1 0_0_0_y_1 0_0_0_y_2 0_0_0_yz_1 0_0_0_yz_2 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yy_2 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_x_1 0_0_0_x_2 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xx_2 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSP_1 SSSP_2 SSSD_1 SSSD_2 SSSF_1 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xx_2 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

INTEGRAL:1 : 2 : 0 : SSSP_1 Y SSSD_0 Y SSSD_1 Y SPSD_0 Y 

0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SSSP_1 SSSD_0 SSSD_1 SPSD_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

INTEGRAL:1 : 3 : 0 : SSSD_1 Y SSSF_0 Y SSSF_1 Y SPSF_0 Y 

0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

SSSD_1 SSSF_0 SSSF_1 SPSF_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SPSF_0 : 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

        }
    }
}


} // derirec namespace
