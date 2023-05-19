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
#include "EriVRRForSPSS.hpp"
#include "EriVRRForSPSP.hpp"
#include "EriVRRForSDSS.hpp"
#include "EriVRRForSDSP.hpp"
#include "EriVRRForSFSP.hpp"
#include "EriHRRForPDSP.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostPDSP(      T*                                 intsBuffer,
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

    BufferHostMY<T, 3> rab(ncpairs);

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
    BufferHostXY<T> pbufSPSS1(3, ncpairs);

0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 
    BufferHostXY<T> pbufSPSS2(3, ncpairs);

0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 
    BufferHostXY<T> pbufSPSP0(9, ncpairs);

0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 
    BufferHostXY<T> pbufSPSP1(9, ncpairs);

0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 
    BufferHostXY<T> pbufSPSP2(9, ncpairs);

0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 
    BufferHostXY<T> pbufSDSS1(3, ncpairs);

0_zz_0_0_1 0_yy_0_0_1 0_xx_0_0_1 
    BufferHostXY<T> pbufSDSP0(18, ncpairs);

0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 
    BufferHostXY<T> pbufSDSP1(12, ncpairs);

0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yz_0_z_1 0_yz_0_y_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xz_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 
    BufferHostXY<T> pbufSFSP0(30, ncpairs);

0_zzz_0_z_0 0_zzz_0_y_0 0_zzz_0_x_0 0_yzz_0_z_0 0_yzz_0_y_0 0_yzz_0_x_0 0_yyz_0_z_0 0_yyz_0_y_0 0_yyz_0_x_0 0_yyy_0_z_0 0_yyy_0_y_0 0_yyy_0_x_0 0_xzz_0_z_0 0_xzz_0_y_0 0_xzz_0_x_0 0_xyz_0_z_0 0_xyz_0_y_0 0_xyz_0_x_0 0_xyy_0_z_0 0_xyy_0_y_0 0_xyy_0_x_0 0_xxz_0_z_0 0_xxz_0_y_0 0_xxz_0_x_0 0_xxy_0_z_0 0_xxy_0_y_0 0_xxy_0_x_0 0_xxx_0_z_0 0_xxx_0_y_0 0_xxx_0_x_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSDSP(18, ncpairs);

0_zz_0_z 0_zz_0_y 0_zz_0_x 0_yz_0_z 0_yz_0_y 0_yz_0_x 0_yy_0_z 0_yy_0_y 0_yy_0_x 0_xz_0_z 0_xz_0_y 0_xz_0_x 0_xy_0_z 0_xy_0_y 0_xy_0_x 0_xx_0_z 0_xx_0_y 0_xx_0_x 
    BufferHostXY<T> cbufSFSP(30, ncpairs);

0_zzz_0_z 0_zzz_0_y 0_zzz_0_x 0_yzz_0_z 0_yzz_0_y 0_yzz_0_x 0_yyz_0_z 0_yyz_0_y 0_yyz_0_x 0_yyy_0_z 0_yyy_0_y 0_yyy_0_x 0_xzz_0_z 0_xzz_0_y 0_xzz_0_x 0_xyz_0_z 0_xyz_0_y 0_xyz_0_x 0_xyy_0_z 0_xyy_0_y 0_xyy_0_x 0_xxz_0_z 0_xxz_0_y 0_xxz_0_x 0_xxy_0_z 0_xxy_0_y 0_xxy_0_x 0_xxx_0_z 0_xxx_0_y 0_xxx_0_x 
    BufferHostXY<T> cbufPDSP(54, ncpairs);

z_zz_0_z z_zz_0_y z_zz_0_x z_yz_0_z z_yz_0_y z_yz_0_x z_yy_0_z z_yy_0_y z_yy_0_x z_xz_0_z z_xz_0_y z_xz_0_x z_xy_0_z z_xy_0_y z_xy_0_x z_xx_0_z z_xx_0_y z_xx_0_x y_zz_0_z y_zz_0_y y_zz_0_x y_yz_0_z y_yz_0_y y_yz_0_x y_yy_0_z y_yy_0_y y_yy_0_x y_xz_0_z y_xz_0_y y_xz_0_x y_xy_0_z y_xy_0_y y_xy_0_x y_xx_0_z y_xx_0_y y_xx_0_x x_zz_0_z x_zz_0_y x_zz_0_x x_yz_0_z x_yz_0_y x_yz_0_x x_yy_0_z x_yy_0_y x_yy_0_x x_xz_0_z x_xz_0_y x_xz_0_x x_xy_0_z x_xy_0_y x_xy_0_x x_xx_0_z x_xx_0_y x_xx_0_x 
    for (int32_t i = 0; i < nppairs; i++)
    {
        derirec::compHostDistancesPT(rpb, gtoPairBlock, bPosition, ePosition, i);

        derirec::compHostFactorsPartialZeta(fz, gtoPairBlock, bPosition, ePosition, i);

        for (int j = i; j < nppairs; j++)
        {
            derirec::compHostDistancesPQ(rpq, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsRho(frho, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsNorm(fnorm, gtoPairBlock, bPosition, ePosition, i, j);

            derirec::compHostFactorsZeta(fze, gtoPairBlock, bPosition, ePosition, i, j);

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
0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_zzz_0_z_0 0_zzz_0_y_0 0_zzz_0_x_0 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yzz_0_z_0 0_yzz_0_y_0 0_yzz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_yyz_0_z_0 0_yyz_0_y_0 0_yyz_0_x_0 0_yyy_0_z_0 0_yyy_0_y_0 0_yyy_0_x_0 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xzz_0_z_0 0_xzz_0_y_0 0_xzz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xyz_0_z_0 0_xyz_0_y_0 0_xyz_0_x_0 0_xyy_0_z_0 0_xyy_0_y_0 0_xyy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 0_xxz_0_z_0 0_xxz_0_y_0 0_xxz_0_x_0 0_xxy_0_z_0 0_xxy_0_y_0 0_xxy_0_x_0 0_xxx_0_z_0 0_xxx_0_y_0 0_xxx_0_x_0 z_zz_0_z_0 z_zz_0_y_0 z_zz_0_x_0 z_yz_0_z_0 z_yz_0_y_0 z_yz_0_x_0 z_yy_0_z_0 z_yy_0_y_0 z_yy_0_x_0 z_xz_0_z_0 z_xz_0_y_0 z_xz_0_x_0 z_xy_0_z_0 z_xy_0_y_0 z_xy_0_x_0 z_xx_0_z_0 z_xx_0_y_0 z_xx_0_x_0 y_zz_0_z_0 y_zz_0_y_0 y_zz_0_x_0 y_yz_0_z_0 y_yz_0_y_0 y_yz_0_x_0 y_yy_0_z_0 y_yy_0_y_0 y_yy_0_x_0 y_xz_0_z_0 y_xz_0_y_0 y_xz_0_x_0 y_xy_0_z_0 y_xy_0_y_0 y_xy_0_x_0 y_xx_0_z_0 y_xx_0_y_0 y_xx_0_x_0 x_zz_0_z_0 x_zz_0_y_0 x_zz_0_x_0 x_yz_0_z_0 x_yz_0_y_0 x_yz_0_x_0 x_yy_0_z_0 x_yy_0_y_0 x_yy_0_x_0 x_xz_0_z_0 x_xz_0_y_0 x_xz_0_x_0 x_xy_0_z_0 x_xy_0_y_0 x_xy_0_x_0 x_xx_0_z_0 x_xx_0_y_0 x_xx_0_x_0 

signature:
0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_0_1 0_zz_0_z_0 0_zz_0_z_1 0_zz_0_y_0 0_zz_0_y_1 0_zz_0_x_0 0_zz_0_x_1 0_zzz_0_z_0 0_zzz_0_y_0 0_zzz_0_x_0 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yz_0_z_0 0_yz_0_z_1 0_yz_0_y_0 0_yz_0_y_1 0_yz_0_x_0 0_yzz_0_z_0 0_yzz_0_y_0 0_yzz_0_x_0 0_yy_0_0_1 0_yy_0_z_0 0_yy_0_z_1 0_yy_0_y_0 0_yy_0_y_1 0_yy_0_x_0 0_yy_0_x_1 0_yyz_0_z_0 0_yyz_0_y_0 0_yyz_0_x_0 0_yyy_0_z_0 0_yyy_0_y_0 0_yyy_0_x_0 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xz_0_x_1 0_xzz_0_z_0 0_xzz_0_y_0 0_xzz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xyz_0_z_0 0_xyz_0_y_0 0_xyz_0_x_0 0_xyy_0_z_0 0_xyy_0_y_0 0_xyy_0_x_0 0_xx_0_0_1 0_xx_0_z_0 0_xx_0_z_1 0_xx_0_y_0 0_xx_0_y_1 0_xx_0_x_0 0_xx_0_x_1 0_xxz_0_z_0 0_xxz_0_y_0 0_xxz_0_x_0 0_xxy_0_z_0 0_xxy_0_y_0 0_xxy_0_x_0 0_xxx_0_z_0 0_xxx_0_y_0 0_xxx_0_x_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yz_0_z_0 0_yz_0_y_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_z_0_0_1 0_zz_0_0_0 0_y_0_0_0 0_y_0_0_1 0_yy_0_0_0 0_x_0_0_0 0_x_0_0_1 0_xx_0_0_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

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

INTEGRAL:1 : 0 : 1 : SSSS_1 Y SSSS_2 Y SPSS_1 Y 

0_0_0_0_1 0_0_0_0_2 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SSSS_1 SSSS_2 SPSS_1 

SSSS_1 : 0_0_0_0_1 

SSSS_2 : 0_0_0_0_2 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

INTEGRAL:1 : 0 : 2 : SSSS_2 Y SSSS_3 Y SPSS_2 Y 

0_0_0_0_2 0_0_0_0_3 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SSSS_2 SSSS_3 SPSS_2 

SSSS_2 : 0_0_0_0_2 

SSSS_3 : 0_0_0_0_3 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

INTEGRAL:1 : 1 : 0 : SSSS_1 Y SSSP_0 Y SSSP_1 Y SPSP_0 Y 

0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

SSSS_1 SSSP_0 SSSP_1 SPSP_0 

SSSS_1 : 0_0_0_0_1 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SPSP_0 : 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

INTEGRAL:1 : 1 : 1 : SSSS_2 Y SSSP_1 Y SSSP_2 Y SPSP_1 Y 

0_0_0_0_2 0_0_0_z_1 0_0_0_z_2 0_0_0_y_1 0_0_0_y_2 0_0_0_x_1 0_0_0_x_2 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SSSS_2 SSSP_1 SSSP_2 SPSP_1 

SSSS_2 : 0_0_0_0_2 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

INTEGRAL:1 : 1 : 2 : SSSS_3 Y SSSP_2 Y SSSP_3 Y SPSP_2 Y 

0_0_0_0_3 0_0_0_z_2 0_0_0_z_3 0_0_0_y_2 0_0_0_y_3 0_0_0_x_2 0_0_0_x_3 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

SSSS_3 SSSP_2 SSSP_3 SPSP_2 

SSSS_3 : 0_0_0_0_3 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SPSP_2 : 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

INTEGRAL:2 : 0 : 1 : SSSS_1 Y SSSS_2 Y SPSS_1 Y SPSS_2 Y SDSS_1 Y 

0_0_0_0_1 0_0_0_0_2 0_z_0_0_1 0_z_0_0_2 0_zz_0_0_1 0_y_0_0_1 0_y_0_0_2 0_yy_0_0_1 0_x_0_0_1 0_x_0_0_2 0_xx_0_0_1 

SSSS_1 SSSS_2 SPSS_1 SPSS_2 SDSS_1 

SSSS_1 : 0_0_0_0_1 

SSSS_2 : 0_0_0_0_2 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SDSS_1 : 0_zz_0_0_1 0_yy_0_0_1 0_xx_0_0_1 

INTEGRAL:2 : 1 : 0 : SSSP_0 Y SSSP_1 Y SPSS_1 Y SPSP_0 Y SPSP_1 Y SDSP_0 Y 

0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

SSSP_0 SSSP_1 SPSS_1 SPSP_0 SPSP_1 SDSP_0 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SPSS_1 : 0_z_0_0_1 0_y_0_0_1 0_x_0_0_1 

SPSP_0 : 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SDSP_0 : 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

INTEGRAL:2 : 1 : 1 : SSSP_1 Y SSSP_2 Y SPSS_2 Y SPSP_1 Y SPSP_2 Y SDSP_1 Y 

0_0_0_z_1 0_0_0_z_2 0_0_0_y_1 0_0_0_y_2 0_0_0_x_1 0_0_0_x_2 0_z_0_0_2 0_z_0_z_1 0_z_0_z_2 0_z_0_y_1 0_z_0_y_2 0_z_0_x_1 0_z_0_x_2 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_y_0_0_2 0_y_0_z_1 0_y_0_z_2 0_y_0_y_1 0_y_0_y_2 0_y_0_x_1 0_y_0_x_2 0_yz_0_z_1 0_yz_0_y_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_x_0_0_2 0_x_0_z_1 0_x_0_z_2 0_x_0_y_1 0_x_0_y_2 0_x_0_x_1 0_x_0_x_2 0_xz_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SSSP_1 SSSP_2 SPSS_2 SPSP_1 SPSP_2 SDSP_1 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SPSP_2 : 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yz_0_z_1 0_yz_0_y_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xz_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

INTEGRAL:3 : 1 : 0 : SPSP_0 Y SPSP_1 Y SDSS_1 Y SDSP_0 N SDSP_1 Y SFSP_0 Y 

0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_0_1 0_zz_0_z_0 0_zz_0_z_1 0_zz_0_y_0 0_zz_0_y_1 0_zz_0_x_0 0_zz_0_x_1 0_zzz_0_z_0 0_zzz_0_y_0 0_zzz_0_x_0 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yz_0_z_0 0_yz_0_z_1 0_yz_0_y_0 0_yz_0_y_1 0_yzz_0_z_0 0_yzz_0_y_0 0_yzz_0_x_0 0_yy_0_0_1 0_yy_0_z_0 0_yy_0_z_1 0_yy_0_y_0 0_yy_0_y_1 0_yy_0_x_0 0_yy_0_x_1 0_yyz_0_z_0 0_yyz_0_y_0 0_yyz_0_x_0 0_yyy_0_z_0 0_yyy_0_y_0 0_yyy_0_x_0 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xz_0_x_0 0_xz_0_x_1 0_xzz_0_z_0 0_xzz_0_y_0 0_xzz_0_x_0 0_xyz_0_z_0 0_xyz_0_y_0 0_xyz_0_x_0 0_xyy_0_z_0 0_xyy_0_y_0 0_xyy_0_x_0 0_xx_0_0_1 0_xx_0_z_0 0_xx_0_z_1 0_xx_0_y_0 0_xx_0_y_1 0_xx_0_x_0 0_xx_0_x_1 0_xxz_0_z_0 0_xxz_0_y_0 0_xxz_0_x_0 0_xxy_0_z_0 0_xxy_0_y_0 0_xxy_0_x_0 0_xxx_0_z_0 0_xxx_0_y_0 0_xxx_0_x_0 

SPSP_0 SPSP_1 SDSS_1 SDSP_0 SDSP_1 SFSP_0 

SPSP_0 : 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SDSS_1 : 0_zz_0_0_1 0_yy_0_0_1 0_xx_0_0_1 

SDSP_0 : 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_yz_0_z_0 0_yz_0_y_0 0_yz_0_x_0 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_xz_0_z_0 0_xz_0_y_0 0_xz_0_x_0 0_xy_0_z_0 0_xy_0_y_0 0_xy_0_x_0 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yz_0_z_1 0_yz_0_y_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xz_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SFSP_0 : 0_zzz_0_z_0 0_zzz_0_y_0 0_zzz_0_x_0 0_yzz_0_z_0 0_yzz_0_y_0 0_yzz_0_x_0 0_yyz_0_z_0 0_yyz_0_y_0 0_yyz_0_x_0 0_yyy_0_z_0 0_yyy_0_y_0 0_yyy_0_x_0 0_xzz_0_z_0 0_xzz_0_y_0 0_xzz_0_x_0 0_xyz_0_z_0 0_xyz_0_y_0 0_xyz_0_x_0 0_xyy_0_z_0 0_xyy_0_y_0 0_xyy_0_x_0 0_xxz_0_z_0 0_xxz_0_y_0 0_xxz_0_x_0 0_xxy_0_z_0 0_xxy_0_y_0 0_xxy_0_x_0 0_xxx_0_z_0 0_xxx_0_y_0 0_xxx_0_x_0 

        }
    }
}


} // derirec namespace
