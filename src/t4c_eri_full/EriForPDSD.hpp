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
#include "EriVRRForSPSS.hpp"
#include "EriVRRForSPSP.hpp"
#include "EriVRRForSPSD.hpp"
#include "EriVRRForSDSP.hpp"
#include "EriVRRForSDSD.hpp"
#include "EriVRRForSFSD.hpp"
#include "EriHRRForPDSD.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostPDSD(      T*                                 intsBuffer,
             const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock,
             const int32_t                            bPosition,
             const int32_t                            ePosition) -> void
{
    // set up dimensions

    const auto nppairs = gtoPairBlock->getNumberOfPrimPairs();

    const auto ncpairs = ePosition - bPosition;

    // allocate Obara-Saika factors

    BufferHostMY<T, 6> osfacts(ncpairs);

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

    BufferHostXY<T> bvals(6, ncpairs);

    CBoysFunc<T, 5> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(6, ncpairs);

    BufferHostXY<T> pbufSSSP0(3, ncpairs);

0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 
    BufferHostXY<T> pbufSSSP1(3, ncpairs);

0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 
    BufferHostXY<T> pbufSSSP2(3, ncpairs);

0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 
    BufferHostXY<T> pbufSSSP3(3, ncpairs);

0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 
    BufferHostXY<T> pbufSSSP4(3, ncpairs);

0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 
    BufferHostXY<T> pbufSSSD0(6, ncpairs);

0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 
    BufferHostXY<T> pbufSSSD1(6, ncpairs);

0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 
    BufferHostXY<T> pbufSSSD2(6, ncpairs);

0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 
    BufferHostXY<T> pbufSSSD3(6, ncpairs);

0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 
    BufferHostXY<T> pbufSPSS2(3, ncpairs);

0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 
    BufferHostXY<T> pbufSPSP1(9, ncpairs);

0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 
    BufferHostXY<T> pbufSPSP2(9, ncpairs);

0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 
    BufferHostXY<T> pbufSPSD0(18, ncpairs);

0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 
    BufferHostXY<T> pbufSPSD1(18, ncpairs);

0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 
    BufferHostXY<T> pbufSPSD2(18, ncpairs);

0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 
    BufferHostXY<T> pbufSDSP1(9, ncpairs);

0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 
    BufferHostXY<T> pbufSDSD0(36, ncpairs);

0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 
    BufferHostXY<T> pbufSDSD1(24, ncpairs);

0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_yy_1 0_zz_0_xz_1 0_zz_0_xy_1 0_zz_0_xx_1 0_yz_0_zz_1 0_yz_0_yz_1 0_yz_0_yy_1 0_yy_0_zz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xx_1 0_xz_0_xz_1 0_xz_0_xx_1 0_xy_0_xy_1 0_xx_0_zz_1 0_xx_0_yz_1 0_xx_0_yy_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 
    BufferHostXY<T> pbufSFSD0(60, ncpairs);

0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSDSD(36, ncpairs);

0_zz_0_zz 0_zz_0_yz 0_zz_0_yy 0_zz_0_xz 0_zz_0_xy 0_zz_0_xx 0_yz_0_zz 0_yz_0_yz 0_yz_0_yy 0_yz_0_xz 0_yz_0_xy 0_yz_0_xx 0_yy_0_zz 0_yy_0_yz 0_yy_0_yy 0_yy_0_xz 0_yy_0_xy 0_yy_0_xx 0_xz_0_zz 0_xz_0_yz 0_xz_0_yy 0_xz_0_xz 0_xz_0_xy 0_xz_0_xx 0_xy_0_zz 0_xy_0_yz 0_xy_0_yy 0_xy_0_xz 0_xy_0_xy 0_xy_0_xx 0_xx_0_zz 0_xx_0_yz 0_xx_0_yy 0_xx_0_xz 0_xx_0_xy 0_xx_0_xx 
    BufferHostXY<T> cbufSFSD(60, ncpairs);

0_zzz_0_zz 0_zzz_0_yz 0_zzz_0_yy 0_zzz_0_xz 0_zzz_0_xy 0_zzz_0_xx 0_yzz_0_zz 0_yzz_0_yz 0_yzz_0_yy 0_yzz_0_xz 0_yzz_0_xy 0_yzz_0_xx 0_yyz_0_zz 0_yyz_0_yz 0_yyz_0_yy 0_yyz_0_xz 0_yyz_0_xy 0_yyz_0_xx 0_yyy_0_zz 0_yyy_0_yz 0_yyy_0_yy 0_yyy_0_xz 0_yyy_0_xy 0_yyy_0_xx 0_xzz_0_zz 0_xzz_0_yz 0_xzz_0_yy 0_xzz_0_xz 0_xzz_0_xy 0_xzz_0_xx 0_xyz_0_zz 0_xyz_0_yz 0_xyz_0_yy 0_xyz_0_xz 0_xyz_0_xy 0_xyz_0_xx 0_xyy_0_zz 0_xyy_0_yz 0_xyy_0_yy 0_xyy_0_xz 0_xyy_0_xy 0_xyy_0_xx 0_xxz_0_zz 0_xxz_0_yz 0_xxz_0_yy 0_xxz_0_xz 0_xxz_0_xy 0_xxz_0_xx 0_xxy_0_zz 0_xxy_0_yz 0_xxy_0_yy 0_xxy_0_xz 0_xxy_0_xy 0_xxy_0_xx 0_xxx_0_zz 0_xxx_0_yz 0_xxx_0_yy 0_xxx_0_xz 0_xxx_0_xy 0_xxx_0_xx 
    BufferHostXY<T> cbufPDSD(108, ncpairs);

z_zz_0_zz z_zz_0_yz z_zz_0_yy z_zz_0_xz z_zz_0_xy z_zz_0_xx z_yz_0_zz z_yz_0_yz z_yz_0_yy z_yz_0_xz z_yz_0_xy z_yz_0_xx z_yy_0_zz z_yy_0_yz z_yy_0_yy z_yy_0_xz z_yy_0_xy z_yy_0_xx z_xz_0_zz z_xz_0_yz z_xz_0_yy z_xz_0_xz z_xz_0_xy z_xz_0_xx z_xy_0_zz z_xy_0_yz z_xy_0_yy z_xy_0_xz z_xy_0_xy z_xy_0_xx z_xx_0_zz z_xx_0_yz z_xx_0_yy z_xx_0_xz z_xx_0_xy z_xx_0_xx y_zz_0_zz y_zz_0_yz y_zz_0_yy y_zz_0_xz y_zz_0_xy y_zz_0_xx y_yz_0_zz y_yz_0_yz y_yz_0_yy y_yz_0_xz y_yz_0_xy y_yz_0_xx y_yy_0_zz y_yy_0_yz y_yy_0_yy y_yy_0_xz y_yy_0_xy y_yy_0_xx y_xz_0_zz y_xz_0_yz y_xz_0_yy y_xz_0_xz y_xz_0_xy y_xz_0_xx y_xy_0_zz y_xy_0_yz y_xy_0_yy y_xy_0_xz y_xy_0_xy y_xy_0_xx y_xx_0_zz y_xx_0_yz y_xx_0_yy y_xx_0_xz y_xx_0_xy y_xx_0_xx x_zz_0_zz x_zz_0_yz x_zz_0_yy x_zz_0_xz x_zz_0_xy x_zz_0_xx x_yz_0_zz x_yz_0_yz x_yz_0_yy x_yz_0_xz x_yz_0_xy x_yz_0_xx x_yy_0_zz x_yy_0_yz x_yy_0_yy x_yy_0_xz x_yy_0_xy x_yy_0_xx x_xz_0_zz x_xz_0_yz x_xz_0_yy x_xz_0_xz x_xz_0_xy x_xz_0_xx x_xy_0_zz x_xy_0_yz x_xy_0_yy x_xy_0_xz x_xy_0_xy x_xy_0_xx x_xx_0_zz x_xx_0_yz x_xx_0_yy x_xx_0_xz x_xx_0_xy x_xx_0_xx 
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
0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 z_zz_0_zz_0 z_zz_0_yz_0 z_zz_0_yy_0 z_zz_0_xz_0 z_zz_0_xy_0 z_zz_0_xx_0 z_yz_0_zz_0 z_yz_0_yz_0 z_yz_0_yy_0 z_yz_0_xz_0 z_yz_0_xy_0 z_yz_0_xx_0 z_yy_0_zz_0 z_yy_0_yz_0 z_yy_0_yy_0 z_yy_0_xz_0 z_yy_0_xy_0 z_yy_0_xx_0 z_xz_0_zz_0 z_xz_0_yz_0 z_xz_0_yy_0 z_xz_0_xz_0 z_xz_0_xy_0 z_xz_0_xx_0 z_xy_0_zz_0 z_xy_0_yz_0 z_xy_0_yy_0 z_xy_0_xz_0 z_xy_0_xy_0 z_xy_0_xx_0 z_xx_0_zz_0 z_xx_0_yz_0 z_xx_0_yy_0 z_xx_0_xz_0 z_xx_0_xy_0 z_xx_0_xx_0 y_zz_0_zz_0 y_zz_0_yz_0 y_zz_0_yy_0 y_zz_0_xz_0 y_zz_0_xy_0 y_zz_0_xx_0 y_yz_0_zz_0 y_yz_0_yz_0 y_yz_0_yy_0 y_yz_0_xz_0 y_yz_0_xy_0 y_yz_0_xx_0 y_yy_0_zz_0 y_yy_0_yz_0 y_yy_0_yy_0 y_yy_0_xz_0 y_yy_0_xy_0 y_yy_0_xx_0 y_xz_0_zz_0 y_xz_0_yz_0 y_xz_0_yy_0 y_xz_0_xz_0 y_xz_0_xy_0 y_xz_0_xx_0 y_xy_0_zz_0 y_xy_0_yz_0 y_xy_0_yy_0 y_xy_0_xz_0 y_xy_0_xy_0 y_xy_0_xx_0 y_xx_0_zz_0 y_xx_0_yz_0 y_xx_0_yy_0 y_xx_0_xz_0 y_xx_0_xy_0 y_xx_0_xx_0 x_zz_0_zz_0 x_zz_0_yz_0 x_zz_0_yy_0 x_zz_0_xz_0 x_zz_0_xy_0 x_zz_0_xx_0 x_yz_0_zz_0 x_yz_0_yz_0 x_yz_0_yy_0 x_yz_0_xz_0 x_yz_0_xy_0 x_yz_0_xx_0 x_yy_0_zz_0 x_yy_0_yz_0 x_yy_0_yy_0 x_yy_0_xz_0 x_yy_0_xy_0 x_yy_0_xx_0 x_xz_0_zz_0 x_xz_0_yz_0 x_xz_0_yy_0 x_xz_0_xz_0 x_xz_0_xy_0 x_xz_0_xx_0 x_xy_0_zz_0 x_xy_0_yz_0 x_xy_0_yy_0 x_xy_0_xz_0 x_xy_0_xy_0 x_xy_0_xx_0 x_xx_0_zz_0 x_xx_0_yz_0 x_xx_0_yy_0 x_xx_0_xz_0 x_xx_0_xy_0 x_xx_0_xx_0 

signature:
0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_yy_0 0_zz_0_yy_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zz_0_xx_0 0_zz_0_xx_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_y_0_zz_0 0_y_0_zz_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_zz_1 0_yz_0_yz_0 0_yz_0_yz_1 0_yz_0_yy_0 0_yz_0_yy_1 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yy_0_z_1 0_yy_0_zz_0 0_yy_0_zz_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yy_0_xx_0 0_yy_0_xx_1 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_x_0_zz_0 0_x_0_zz_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xz_1 0_xz_0_xy_0 0_xz_0_xx_0 0_xz_0_xx_1 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xy_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xx_0_z_1 0_xx_0_zz_0 0_xx_0_zz_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_yy_0 0_xx_0_yy_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xz_0_xx_0 0_xy_0_xy_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_z_0_0_0 0_y_0_0_0 0_x_0_0_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

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

INTEGRAL:0 : 1 : 4 : SSSS_4 Y SSSS_5 Y SSSP_4 Y 

0_0_0_0_4 0_0_0_0_5 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSS_4 SSSS_5 SSSP_4 

SSSS_4 : 0_0_0_0_4 

SSSS_5 : 0_0_0_0_5 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

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

0_0_0_0_2 0_0_0_0_3 0_0_0_z_2 0_0_0_z_3 0_0_0_zz_2 0_0_0_y_2 0_0_0_y_3 0_0_0_yz_2 0_0_0_yy_2 0_0_0_x_2 0_0_0_x_3 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSS_2 SSSS_3 SSSP_2 SSSP_3 SSSD_2 

SSSS_2 : 0_0_0_0_2 

SSSS_3 : 0_0_0_0_3 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

INTEGRAL:0 : 2 : 3 : SSSS_3 Y SSSS_4 Y SSSP_3 Y SSSP_4 Y SSSD_3 Y 

0_0_0_0_3 0_0_0_0_4 0_0_0_z_3 0_0_0_z_4 0_0_0_zz_3 0_0_0_y_3 0_0_0_y_4 0_0_0_yz_3 0_0_0_yy_3 0_0_0_x_3 0_0_0_x_4 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSS_3 SSSS_4 SSSP_3 SSSP_4 SSSD_3 

SSSS_3 : 0_0_0_0_3 

SSSS_4 : 0_0_0_0_4 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

INTEGRAL:1 : 0 : 2 : SSSS_2 Y SSSS_3 Y SPSS_2 Y 

0_0_0_0_2 0_0_0_0_3 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SSSS_2 SSSS_3 SPSS_2 

SSSS_2 : 0_0_0_0_2 

SSSS_3 : 0_0_0_0_3 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

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

INTEGRAL:1 : 2 : 0 : SSSP_1 Y SSSD_0 Y SSSD_1 Y SPSD_0 Y 

0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SSSP_1 SSSD_0 SSSD_1 SPSD_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

INTEGRAL:1 : 2 : 1 : SSSP_2 Y SSSD_1 Y SSSD_2 Y SPSD_1 Y 

0_0_0_z_2 0_0_0_zz_1 0_0_0_zz_2 0_0_0_y_2 0_0_0_yz_1 0_0_0_yz_2 0_0_0_yy_1 0_0_0_yy_2 0_0_0_x_2 0_0_0_xz_1 0_0_0_xz_2 0_0_0_xy_1 0_0_0_xy_2 0_0_0_xx_1 0_0_0_xx_2 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SSSP_2 SSSD_1 SSSD_2 SPSD_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

INTEGRAL:1 : 2 : 2 : SSSP_3 Y SSSD_2 Y SSSD_3 Y SPSD_2 Y 

0_0_0_z_3 0_0_0_zz_2 0_0_0_zz_3 0_0_0_y_3 0_0_0_yz_2 0_0_0_yz_3 0_0_0_yy_2 0_0_0_yy_3 0_0_0_x_3 0_0_0_xz_2 0_0_0_xz_3 0_0_0_xy_2 0_0_0_xy_3 0_0_0_xx_2 0_0_0_xx_3 0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 

SSSP_3 SSSD_2 SSSD_3 SPSD_2 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SPSD_2 : 0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 

INTEGRAL:2 : 1 : 1 : SSSP_1 Y SSSP_2 Y SPSS_2 Y SPSP_1 Y SPSP_2 Y SDSP_1 Y 

0_0_0_z_1 0_0_0_z_2 0_0_0_y_1 0_0_0_y_2 0_0_0_x_1 0_0_0_x_2 0_z_0_0_2 0_z_0_z_1 0_z_0_z_2 0_z_0_y_1 0_z_0_y_2 0_z_0_x_1 0_z_0_x_2 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_y_0_0_2 0_y_0_z_1 0_y_0_z_2 0_y_0_y_1 0_y_0_y_2 0_y_0_x_1 0_y_0_x_2 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_x_0_0_2 0_x_0_z_1 0_x_0_z_2 0_x_0_y_1 0_x_0_y_2 0_x_0_x_1 0_x_0_x_2 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SSSP_1 SSSP_2 SPSS_2 SPSP_1 SPSP_2 SDSP_1 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SPSS_2 : 0_z_0_0_2 0_y_0_0_2 0_x_0_0_2 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SPSP_2 : 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

INTEGRAL:2 : 2 : 0 : SSSD_0 Y SSSD_1 Y SPSP_1 Y SPSD_0 Y SPSD_1 Y SDSD_0 Y 

0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

SSSD_0 SSSD_1 SPSP_1 SPSD_0 SPSD_1 SDSD_0 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SDSD_0 : 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

INTEGRAL:2 : 2 : 1 : SSSD_1 Y SSSD_2 Y SPSP_2 Y SPSD_1 Y SPSD_2 Y SDSD_1 Y 

0_0_0_zz_1 0_0_0_zz_2 0_0_0_yz_1 0_0_0_yz_2 0_0_0_yy_1 0_0_0_yy_2 0_0_0_xz_1 0_0_0_xz_2 0_0_0_xy_1 0_0_0_xy_2 0_0_0_xx_1 0_0_0_xx_2 0_z_0_z_2 0_z_0_zz_1 0_z_0_zz_2 0_z_0_y_2 0_z_0_yz_1 0_z_0_yz_2 0_z_0_yy_1 0_z_0_yy_2 0_z_0_x_2 0_z_0_xz_1 0_z_0_xz_2 0_z_0_xy_1 0_z_0_xy_2 0_z_0_xx_1 0_z_0_xx_2 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_yy_1 0_zz_0_xz_1 0_zz_0_xy_1 0_zz_0_xx_1 0_y_0_z_2 0_y_0_zz_1 0_y_0_zz_2 0_y_0_y_2 0_y_0_yz_1 0_y_0_yz_2 0_y_0_yy_1 0_y_0_yy_2 0_y_0_x_2 0_y_0_xz_1 0_y_0_xz_2 0_y_0_xy_1 0_y_0_xy_2 0_y_0_xx_1 0_y_0_xx_2 0_yz_0_zz_1 0_yz_0_yz_1 0_yz_0_yy_1 0_yy_0_zz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xx_1 0_x_0_z_2 0_x_0_zz_1 0_x_0_zz_2 0_x_0_y_2 0_x_0_yz_1 0_x_0_yz_2 0_x_0_yy_1 0_x_0_yy_2 0_x_0_x_2 0_x_0_xz_1 0_x_0_xz_2 0_x_0_xy_1 0_x_0_xy_2 0_x_0_xx_1 0_x_0_xx_2 0_xz_0_xz_1 0_xz_0_xx_1 0_xy_0_xy_1 0_xx_0_zz_1 0_xx_0_yz_1 0_xx_0_yy_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SSSD_1 SSSD_2 SPSP_2 SPSD_1 SPSD_2 SDSD_1 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SPSP_2 : 0_z_0_z_2 0_z_0_y_2 0_z_0_x_2 0_y_0_z_2 0_y_0_y_2 0_y_0_x_2 0_x_0_z_2 0_x_0_y_2 0_x_0_x_2 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SPSD_2 : 0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_yy_1 0_zz_0_xz_1 0_zz_0_xy_1 0_zz_0_xx_1 0_yz_0_zz_1 0_yz_0_yz_1 0_yz_0_yy_1 0_yy_0_zz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xx_1 0_xz_0_xz_1 0_xz_0_xx_1 0_xy_0_xy_1 0_xx_0_zz_1 0_xx_0_yz_1 0_xx_0_yy_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

INTEGRAL:3 : 2 : 0 : SPSD_0 Y SPSD_1 Y SDSP_1 Y SDSD_0 N SDSD_1 Y SFSD_0 Y 

0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_yy_0 0_zz_0_yy_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zz_0_xx_0 0_zz_0_xx_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_y_0_zz_0 0_y_0_zz_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_zz_1 0_yz_0_yz_0 0_yz_0_yz_1 0_yz_0_yy_0 0_yz_0_yy_1 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yy_0_z_1 0_yy_0_zz_0 0_yy_0_zz_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yy_0_xx_0 0_yy_0_xx_1 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_x_0_zz_0 0_x_0_zz_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xz_0_xz_1 0_xz_0_xx_0 0_xz_0_xx_1 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xx_0_z_1 0_xx_0_zz_0 0_xx_0_zz_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_yy_0 0_xx_0_yy_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

SPSD_0 SPSD_1 SDSP_1 SDSD_0 SDSD_1 SFSD_0 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SDSD_0 : 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_yy_1 0_zz_0_xz_1 0_zz_0_xy_1 0_zz_0_xx_1 0_yz_0_zz_1 0_yz_0_yz_1 0_yz_0_yy_1 0_yy_0_zz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xx_1 0_xz_0_xz_1 0_xz_0_xx_1 0_xy_0_xy_1 0_xx_0_zz_1 0_xx_0_yz_1 0_xx_0_yy_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SFSD_0 : 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

        }
    }
}


} // derirec namespace
