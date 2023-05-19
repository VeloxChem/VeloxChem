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
#include "EriVRRForSPSS.hpp"
#include "EriVRRForSPSP.hpp"
#include "EriVRRForSPSD.hpp"
#include "EriVRRForSPSF.hpp"
#include "EriVRRForSDSP.hpp"
#include "EriVRRForSDSD.hpp"
#include "EriVRRForSDSF.hpp"
#include "EriVRRForSFSD.hpp"
#include "EriVRRForSFSF.hpp"
#include "EriHRRForSDPD.hpp"
#include "EriHRRForSFPD.hpp"
#include "EriHRRForPDPD.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostPDPD(      T*                                 intsBuffer,
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

    BufferHostMY<T, 3> rcd(ncpairs);

    BufferHostMY<T, 3> rpb(ncpairs);

    BufferHostMY<T, 3> rqd(ncpairs);

    BufferHostMY<T, 3> rwp(ncpairs);

    BufferHostMY<T, 3> rwq(ncpairs);

    // allocate coordinates

    BufferHostMY<T, 3> rw(ncpairs); 

    // allocate Boys function data

    BufferHostX<T> bargs(ncpairs);

    BufferHostXY<T> bvals(7, ncpairs);

    CBoysFunc<T, 6> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(7, ncpairs);

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
    BufferHostXY<T> pbufSSSP5(3, ncpairs);

0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 
    BufferHostXY<T> pbufSSSD0(6, ncpairs);

0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 
    BufferHostXY<T> pbufSSSD1(6, ncpairs);

0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 
    BufferHostXY<T> pbufSSSD2(6, ncpairs);

0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 
    BufferHostXY<T> pbufSSSD3(6, ncpairs);

0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 
    BufferHostXY<T> pbufSSSD4(4, ncpairs);

0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xx_4 
    BufferHostXY<T> pbufSSSF0(10, ncpairs);

0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 
    BufferHostXY<T> pbufSSSF1(10, ncpairs);

0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 
    BufferHostXY<T> pbufSSSF2(10, ncpairs);

0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 
    BufferHostXY<T> pbufSSSF3(10, ncpairs);

0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 
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
    BufferHostXY<T> pbufSPSF0(30, ncpairs);

0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 
    BufferHostXY<T> pbufSPSF1(30, ncpairs);

0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 
    BufferHostXY<T> pbufSPSF2(30, ncpairs);

0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_yyy_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xyy_2 0_z_0_xxz_2 0_z_0_xxy_2 0_z_0_xxx_2 0_y_0_zzz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xzz_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxz_2 0_y_0_xxy_2 0_y_0_xxx_2 0_x_0_zzz_2 0_x_0_yzz_2 0_x_0_yyz_2 0_x_0_yyy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 
    BufferHostXY<T> pbufSDSP1(9, ncpairs);

0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 
    BufferHostXY<T> pbufSDSD0(36, ncpairs);

0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 
    BufferHostXY<T> pbufSDSD1(24, ncpairs);

0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_yy_1 0_zz_0_xz_1 0_zz_0_xy_1 0_zz_0_xx_1 0_yz_0_zz_1 0_yz_0_yz_1 0_yz_0_yy_1 0_yy_0_zz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xx_1 0_xz_0_xz_1 0_xz_0_xx_1 0_xy_0_xy_1 0_xx_0_zz_1 0_xx_0_yz_1 0_xx_0_yy_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 
    BufferHostXY<T> pbufSDSF0(60, ncpairs);

0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 
    BufferHostXY<T> pbufSDSF1(40, ncpairs);

0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_yyy_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xyy_1 0_zz_0_xxz_1 0_zz_0_xxy_1 0_zz_0_xxx_1 0_yz_0_zzz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_yyy_1 0_yz_0_xyz_1 0_yy_0_zzz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xzz_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxz_1 0_yy_0_xxy_1 0_yy_0_xxx_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xz_0_xxx_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_zzz_1 0_xx_0_yzz_1 0_xx_0_yyz_1 0_xx_0_yyy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 
    BufferHostXY<T> pbufSFSD0(60, ncpairs);

0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 
    BufferHostXY<T> pbufSFSF0(100, ncpairs);

0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_yyz_0 0_zzz_0_yyy_0 0_zzz_0_xzz_0 0_zzz_0_xyz_0 0_zzz_0_xyy_0 0_zzz_0_xxz_0 0_zzz_0_xxy_0 0_zzz_0_xxx_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_yyy_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yzz_0_xyy_0 0_yzz_0_xxz_0 0_yzz_0_xxy_0 0_yzz_0_xxx_0 0_yyz_0_zzz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xzz_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xxz_0 0_yyz_0_xxy_0 0_yyz_0_xxx_0 0_yyy_0_zzz_0 0_yyy_0_yzz_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xzz_0 0_yyy_0_xyz_0 0_yyy_0_xyy_0 0_yyy_0_xxz_0 0_yyy_0_xxy_0 0_yyy_0_xxx_0 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_yyz_0 0_xzz_0_yyy_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xyy_0 0_xzz_0_xxz_0 0_xzz_0_xxy_0 0_xzz_0_xxx_0 0_xyz_0_zzz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_yyy_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxx_0 0_xyy_0_zzz_0 0_xyy_0_yzz_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xzz_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxz_0 0_xyy_0_xxy_0 0_xyy_0_xxx_0 0_xxz_0_zzz_0 0_xxz_0_yzz_0 0_xxz_0_yyz_0 0_xxz_0_yyy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xyy_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_zzz_0 0_xxy_0_yzz_0 0_xxy_0_yyz_0 0_xxy_0_yyy_0 0_xxy_0_xzz_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_zzz_0 0_xxx_0_yzz_0 0_xxx_0_yyz_0 0_xxx_0_yyy_0 0_xxx_0_xzz_0 0_xxx_0_xyz_0 0_xxx_0_xyy_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSDSD(36, ncpairs);

0_zz_0_zz 0_zz_0_yz 0_zz_0_yy 0_zz_0_xz 0_zz_0_xy 0_zz_0_xx 0_yz_0_zz 0_yz_0_yz 0_yz_0_yy 0_yz_0_xz 0_yz_0_xy 0_yz_0_xx 0_yy_0_zz 0_yy_0_yz 0_yy_0_yy 0_yy_0_xz 0_yy_0_xy 0_yy_0_xx 0_xz_0_zz 0_xz_0_yz 0_xz_0_yy 0_xz_0_xz 0_xz_0_xy 0_xz_0_xx 0_xy_0_zz 0_xy_0_yz 0_xy_0_yy 0_xy_0_xz 0_xy_0_xy 0_xy_0_xx 0_xx_0_zz 0_xx_0_yz 0_xx_0_yy 0_xx_0_xz 0_xx_0_xy 0_xx_0_xx 
    BufferHostXY<T> cbufSDSF(60, ncpairs);

0_zz_0_zzz 0_zz_0_yzz 0_zz_0_yyz 0_zz_0_yyy 0_zz_0_xzz 0_zz_0_xyz 0_zz_0_xyy 0_zz_0_xxz 0_zz_0_xxy 0_zz_0_xxx 0_yz_0_zzz 0_yz_0_yzz 0_yz_0_yyz 0_yz_0_yyy 0_yz_0_xzz 0_yz_0_xyz 0_yz_0_xyy 0_yz_0_xxz 0_yz_0_xxy 0_yz_0_xxx 0_yy_0_zzz 0_yy_0_yzz 0_yy_0_yyz 0_yy_0_yyy 0_yy_0_xzz 0_yy_0_xyz 0_yy_0_xyy 0_yy_0_xxz 0_yy_0_xxy 0_yy_0_xxx 0_xz_0_zzz 0_xz_0_yzz 0_xz_0_yyz 0_xz_0_yyy 0_xz_0_xzz 0_xz_0_xyz 0_xz_0_xyy 0_xz_0_xxz 0_xz_0_xxy 0_xz_0_xxx 0_xy_0_zzz 0_xy_0_yzz 0_xy_0_yyz 0_xy_0_yyy 0_xy_0_xzz 0_xy_0_xyz 0_xy_0_xyy 0_xy_0_xxz 0_xy_0_xxy 0_xy_0_xxx 0_xx_0_zzz 0_xx_0_yzz 0_xx_0_yyz 0_xx_0_yyy 0_xx_0_xzz 0_xx_0_xyz 0_xx_0_xyy 0_xx_0_xxz 0_xx_0_xxy 0_xx_0_xxx 
    BufferHostXY<T> cbufSFSD(60, ncpairs);

0_zzz_0_zz 0_zzz_0_yz 0_zzz_0_yy 0_zzz_0_xz 0_zzz_0_xy 0_zzz_0_xx 0_yzz_0_zz 0_yzz_0_yz 0_yzz_0_yy 0_yzz_0_xz 0_yzz_0_xy 0_yzz_0_xx 0_yyz_0_zz 0_yyz_0_yz 0_yyz_0_yy 0_yyz_0_xz 0_yyz_0_xy 0_yyz_0_xx 0_yyy_0_zz 0_yyy_0_yz 0_yyy_0_yy 0_yyy_0_xz 0_yyy_0_xy 0_yyy_0_xx 0_xzz_0_zz 0_xzz_0_yz 0_xzz_0_yy 0_xzz_0_xz 0_xzz_0_xy 0_xzz_0_xx 0_xyz_0_zz 0_xyz_0_yz 0_xyz_0_yy 0_xyz_0_xz 0_xyz_0_xy 0_xyz_0_xx 0_xyy_0_zz 0_xyy_0_yz 0_xyy_0_yy 0_xyy_0_xz 0_xyy_0_xy 0_xyy_0_xx 0_xxz_0_zz 0_xxz_0_yz 0_xxz_0_yy 0_xxz_0_xz 0_xxz_0_xy 0_xxz_0_xx 0_xxy_0_zz 0_xxy_0_yz 0_xxy_0_yy 0_xxy_0_xz 0_xxy_0_xy 0_xxy_0_xx 0_xxx_0_zz 0_xxx_0_yz 0_xxx_0_yy 0_xxx_0_xz 0_xxx_0_xy 0_xxx_0_xx 
    BufferHostXY<T> cbufSFSF(100, ncpairs);

0_zzz_0_zzz 0_zzz_0_yzz 0_zzz_0_yyz 0_zzz_0_yyy 0_zzz_0_xzz 0_zzz_0_xyz 0_zzz_0_xyy 0_zzz_0_xxz 0_zzz_0_xxy 0_zzz_0_xxx 0_yzz_0_zzz 0_yzz_0_yzz 0_yzz_0_yyz 0_yzz_0_yyy 0_yzz_0_xzz 0_yzz_0_xyz 0_yzz_0_xyy 0_yzz_0_xxz 0_yzz_0_xxy 0_yzz_0_xxx 0_yyz_0_zzz 0_yyz_0_yzz 0_yyz_0_yyz 0_yyz_0_yyy 0_yyz_0_xzz 0_yyz_0_xyz 0_yyz_0_xyy 0_yyz_0_xxz 0_yyz_0_xxy 0_yyz_0_xxx 0_yyy_0_zzz 0_yyy_0_yzz 0_yyy_0_yyz 0_yyy_0_yyy 0_yyy_0_xzz 0_yyy_0_xyz 0_yyy_0_xyy 0_yyy_0_xxz 0_yyy_0_xxy 0_yyy_0_xxx 0_xzz_0_zzz 0_xzz_0_yzz 0_xzz_0_yyz 0_xzz_0_yyy 0_xzz_0_xzz 0_xzz_0_xyz 0_xzz_0_xyy 0_xzz_0_xxz 0_xzz_0_xxy 0_xzz_0_xxx 0_xyz_0_zzz 0_xyz_0_yzz 0_xyz_0_yyz 0_xyz_0_yyy 0_xyz_0_xzz 0_xyz_0_xyz 0_xyz_0_xyy 0_xyz_0_xxz 0_xyz_0_xxy 0_xyz_0_xxx 0_xyy_0_zzz 0_xyy_0_yzz 0_xyy_0_yyz 0_xyy_0_yyy 0_xyy_0_xzz 0_xyy_0_xyz 0_xyy_0_xyy 0_xyy_0_xxz 0_xyy_0_xxy 0_xyy_0_xxx 0_xxz_0_zzz 0_xxz_0_yzz 0_xxz_0_yyz 0_xxz_0_yyy 0_xxz_0_xzz 0_xxz_0_xyz 0_xxz_0_xyy 0_xxz_0_xxz 0_xxz_0_xxy 0_xxz_0_xxx 0_xxy_0_zzz 0_xxy_0_yzz 0_xxy_0_yyz 0_xxy_0_yyy 0_xxy_0_xzz 0_xxy_0_xyz 0_xxy_0_xyy 0_xxy_0_xxz 0_xxy_0_xxy 0_xxy_0_xxx 0_xxx_0_zzz 0_xxx_0_yzz 0_xxx_0_yyz 0_xxx_0_yyy 0_xxx_0_xzz 0_xxx_0_xyz 0_xxx_0_xyy 0_xxx_0_xxz 0_xxx_0_xxy 0_xxx_0_xxx 
    BufferHostXY<T> cbufSDPD(108, ncpairs);

0_zz_z_zz 0_zz_z_yz 0_zz_z_yy 0_zz_z_xz 0_zz_z_xy 0_zz_z_xx 0_zz_y_zz 0_zz_y_yz 0_zz_y_yy 0_zz_y_xz 0_zz_y_xy 0_zz_y_xx 0_zz_x_zz 0_zz_x_yz 0_zz_x_yy 0_zz_x_xz 0_zz_x_xy 0_zz_x_xx 0_yz_z_zz 0_yz_z_yz 0_yz_z_yy 0_yz_z_xz 0_yz_z_xy 0_yz_z_xx 0_yz_y_zz 0_yz_y_yz 0_yz_y_yy 0_yz_y_xz 0_yz_y_xy 0_yz_y_xx 0_yz_x_zz 0_yz_x_yz 0_yz_x_yy 0_yz_x_xz 0_yz_x_xy 0_yz_x_xx 0_yy_z_zz 0_yy_z_yz 0_yy_z_yy 0_yy_z_xz 0_yy_z_xy 0_yy_z_xx 0_yy_y_zz 0_yy_y_yz 0_yy_y_yy 0_yy_y_xz 0_yy_y_xy 0_yy_y_xx 0_yy_x_zz 0_yy_x_yz 0_yy_x_yy 0_yy_x_xz 0_yy_x_xy 0_yy_x_xx 0_xz_z_zz 0_xz_z_yz 0_xz_z_yy 0_xz_z_xz 0_xz_z_xy 0_xz_z_xx 0_xz_y_zz 0_xz_y_yz 0_xz_y_yy 0_xz_y_xz 0_xz_y_xy 0_xz_y_xx 0_xz_x_zz 0_xz_x_yz 0_xz_x_yy 0_xz_x_xz 0_xz_x_xy 0_xz_x_xx 0_xy_z_zz 0_xy_z_yz 0_xy_z_yy 0_xy_z_xz 0_xy_z_xy 0_xy_z_xx 0_xy_y_zz 0_xy_y_yz 0_xy_y_yy 0_xy_y_xz 0_xy_y_xy 0_xy_y_xx 0_xy_x_zz 0_xy_x_yz 0_xy_x_yy 0_xy_x_xz 0_xy_x_xy 0_xy_x_xx 0_xx_z_zz 0_xx_z_yz 0_xx_z_yy 0_xx_z_xz 0_xx_z_xy 0_xx_z_xx 0_xx_y_zz 0_xx_y_yz 0_xx_y_yy 0_xx_y_xz 0_xx_y_xy 0_xx_y_xx 0_xx_x_zz 0_xx_x_yz 0_xx_x_yy 0_xx_x_xz 0_xx_x_xy 0_xx_x_xx 
    BufferHostXY<T> cbufSFPD(180, ncpairs);

0_zzz_z_zz 0_zzz_z_yz 0_zzz_z_yy 0_zzz_z_xz 0_zzz_z_xy 0_zzz_z_xx 0_zzz_y_zz 0_zzz_y_yz 0_zzz_y_yy 0_zzz_y_xz 0_zzz_y_xy 0_zzz_y_xx 0_zzz_x_zz 0_zzz_x_yz 0_zzz_x_yy 0_zzz_x_xz 0_zzz_x_xy 0_zzz_x_xx 0_yzz_z_zz 0_yzz_z_yz 0_yzz_z_yy 0_yzz_z_xz 0_yzz_z_xy 0_yzz_z_xx 0_yzz_y_zz 0_yzz_y_yz 0_yzz_y_yy 0_yzz_y_xz 0_yzz_y_xy 0_yzz_y_xx 0_yzz_x_zz 0_yzz_x_yz 0_yzz_x_yy 0_yzz_x_xz 0_yzz_x_xy 0_yzz_x_xx 0_yyz_z_zz 0_yyz_z_yz 0_yyz_z_yy 0_yyz_z_xz 0_yyz_z_xy 0_yyz_z_xx 0_yyz_y_zz 0_yyz_y_yz 0_yyz_y_yy 0_yyz_y_xz 0_yyz_y_xy 0_yyz_y_xx 0_yyz_x_zz 0_yyz_x_yz 0_yyz_x_yy 0_yyz_x_xz 0_yyz_x_xy 0_yyz_x_xx 0_yyy_z_zz 0_yyy_z_yz 0_yyy_z_yy 0_yyy_z_xz 0_yyy_z_xy 0_yyy_z_xx 0_yyy_y_zz 0_yyy_y_yz 0_yyy_y_yy 0_yyy_y_xz 0_yyy_y_xy 0_yyy_y_xx 0_yyy_x_zz 0_yyy_x_yz 0_yyy_x_yy 0_yyy_x_xz 0_yyy_x_xy 0_yyy_x_xx 0_xzz_z_zz 0_xzz_z_yz 0_xzz_z_yy 0_xzz_z_xz 0_xzz_z_xy 0_xzz_z_xx 0_xzz_y_zz 0_xzz_y_yz 0_xzz_y_yy 0_xzz_y_xz 0_xzz_y_xy 0_xzz_y_xx 0_xzz_x_zz 0_xzz_x_yz 0_xzz_x_yy 0_xzz_x_xz 0_xzz_x_xy 0_xzz_x_xx 0_xyz_z_zz 0_xyz_z_yz 0_xyz_z_yy 0_xyz_z_xz 0_xyz_z_xy 0_xyz_z_xx 0_xyz_y_zz 0_xyz_y_yz 0_xyz_y_yy 0_xyz_y_xz 0_xyz_y_xy 0_xyz_y_xx 0_xyz_x_zz 0_xyz_x_yz 0_xyz_x_yy 0_xyz_x_xz 0_xyz_x_xy 0_xyz_x_xx 0_xyy_z_zz 0_xyy_z_yz 0_xyy_z_yy 0_xyy_z_xz 0_xyy_z_xy 0_xyy_z_xx 0_xyy_y_zz 0_xyy_y_yz 0_xyy_y_yy 0_xyy_y_xz 0_xyy_y_xy 0_xyy_y_xx 0_xyy_x_zz 0_xyy_x_yz 0_xyy_x_yy 0_xyy_x_xz 0_xyy_x_xy 0_xyy_x_xx 0_xxz_z_zz 0_xxz_z_yz 0_xxz_z_yy 0_xxz_z_xz 0_xxz_z_xy 0_xxz_z_xx 0_xxz_y_zz 0_xxz_y_yz 0_xxz_y_yy 0_xxz_y_xz 0_xxz_y_xy 0_xxz_y_xx 0_xxz_x_zz 0_xxz_x_yz 0_xxz_x_yy 0_xxz_x_xz 0_xxz_x_xy 0_xxz_x_xx 0_xxy_z_zz 0_xxy_z_yz 0_xxy_z_yy 0_xxy_z_xz 0_xxy_z_xy 0_xxy_z_xx 0_xxy_y_zz 0_xxy_y_yz 0_xxy_y_yy 0_xxy_y_xz 0_xxy_y_xy 0_xxy_y_xx 0_xxy_x_zz 0_xxy_x_yz 0_xxy_x_yy 0_xxy_x_xz 0_xxy_x_xy 0_xxy_x_xx 0_xxx_z_zz 0_xxx_z_yz 0_xxx_z_yy 0_xxx_z_xz 0_xxx_z_xy 0_xxx_z_xx 0_xxx_y_zz 0_xxx_y_yz 0_xxx_y_yy 0_xxx_y_xz 0_xxx_y_xy 0_xxx_y_xx 0_xxx_x_zz 0_xxx_x_yz 0_xxx_x_yy 0_xxx_x_xz 0_xxx_x_xy 0_xxx_x_xx 
    BufferHostXY<T> cbufPDPD(324, ncpairs);

z_zz_z_zz z_zz_z_yz z_zz_z_yy z_zz_z_xz z_zz_z_xy z_zz_z_xx z_zz_y_zz z_zz_y_yz z_zz_y_yy z_zz_y_xz z_zz_y_xy z_zz_y_xx z_zz_x_zz z_zz_x_yz z_zz_x_yy z_zz_x_xz z_zz_x_xy z_zz_x_xx z_yz_z_zz z_yz_z_yz z_yz_z_yy z_yz_z_xz z_yz_z_xy z_yz_z_xx z_yz_y_zz z_yz_y_yz z_yz_y_yy z_yz_y_xz z_yz_y_xy z_yz_y_xx z_yz_x_zz z_yz_x_yz z_yz_x_yy z_yz_x_xz z_yz_x_xy z_yz_x_xx z_yy_z_zz z_yy_z_yz z_yy_z_yy z_yy_z_xz z_yy_z_xy z_yy_z_xx z_yy_y_zz z_yy_y_yz z_yy_y_yy z_yy_y_xz z_yy_y_xy z_yy_y_xx z_yy_x_zz z_yy_x_yz z_yy_x_yy z_yy_x_xz z_yy_x_xy z_yy_x_xx z_xz_z_zz z_xz_z_yz z_xz_z_yy z_xz_z_xz z_xz_z_xy z_xz_z_xx z_xz_y_zz z_xz_y_yz z_xz_y_yy z_xz_y_xz z_xz_y_xy z_xz_y_xx z_xz_x_zz z_xz_x_yz z_xz_x_yy z_xz_x_xz z_xz_x_xy z_xz_x_xx z_xy_z_zz z_xy_z_yz z_xy_z_yy z_xy_z_xz z_xy_z_xy z_xy_z_xx z_xy_y_zz z_xy_y_yz z_xy_y_yy z_xy_y_xz z_xy_y_xy z_xy_y_xx z_xy_x_zz z_xy_x_yz z_xy_x_yy z_xy_x_xz z_xy_x_xy z_xy_x_xx z_xx_z_zz z_xx_z_yz z_xx_z_yy z_xx_z_xz z_xx_z_xy z_xx_z_xx z_xx_y_zz z_xx_y_yz z_xx_y_yy z_xx_y_xz z_xx_y_xy z_xx_y_xx z_xx_x_zz z_xx_x_yz z_xx_x_yy z_xx_x_xz z_xx_x_xy z_xx_x_xx y_zz_z_zz y_zz_z_yz y_zz_z_yy y_zz_z_xz y_zz_z_xy y_zz_z_xx y_zz_y_zz y_zz_y_yz y_zz_y_yy y_zz_y_xz y_zz_y_xy y_zz_y_xx y_zz_x_zz y_zz_x_yz y_zz_x_yy y_zz_x_xz y_zz_x_xy y_zz_x_xx y_yz_z_zz y_yz_z_yz y_yz_z_yy y_yz_z_xz y_yz_z_xy y_yz_z_xx y_yz_y_zz y_yz_y_yz y_yz_y_yy y_yz_y_xz y_yz_y_xy y_yz_y_xx y_yz_x_zz y_yz_x_yz y_yz_x_yy y_yz_x_xz y_yz_x_xy y_yz_x_xx y_yy_z_zz y_yy_z_yz y_yy_z_yy y_yy_z_xz y_yy_z_xy y_yy_z_xx y_yy_y_zz y_yy_y_yz y_yy_y_yy y_yy_y_xz y_yy_y_xy y_yy_y_xx y_yy_x_zz y_yy_x_yz y_yy_x_yy y_yy_x_xz y_yy_x_xy y_yy_x_xx y_xz_z_zz y_xz_z_yz y_xz_z_yy y_xz_z_xz y_xz_z_xy y_xz_z_xx y_xz_y_zz y_xz_y_yz y_xz_y_yy y_xz_y_xz y_xz_y_xy y_xz_y_xx y_xz_x_zz y_xz_x_yz y_xz_x_yy y_xz_x_xz y_xz_x_xy y_xz_x_xx y_xy_z_zz y_xy_z_yz y_xy_z_yy y_xy_z_xz y_xy_z_xy y_xy_z_xx y_xy_y_zz y_xy_y_yz y_xy_y_yy y_xy_y_xz y_xy_y_xy y_xy_y_xx y_xy_x_zz y_xy_x_yz y_xy_x_yy y_xy_x_xz y_xy_x_xy y_xy_x_xx y_xx_z_zz y_xx_z_yz y_xx_z_yy y_xx_z_xz y_xx_z_xy y_xx_z_xx y_xx_y_zz y_xx_y_yz y_xx_y_yy y_xx_y_xz y_xx_y_xy y_xx_y_xx y_xx_x_zz y_xx_x_yz y_xx_x_yy y_xx_x_xz y_xx_x_xy y_xx_x_xx x_zz_z_zz x_zz_z_yz x_zz_z_yy x_zz_z_xz x_zz_z_xy x_zz_z_xx x_zz_y_zz x_zz_y_yz x_zz_y_yy x_zz_y_xz x_zz_y_xy x_zz_y_xx x_zz_x_zz x_zz_x_yz x_zz_x_yy x_zz_x_xz x_zz_x_xy x_zz_x_xx x_yz_z_zz x_yz_z_yz x_yz_z_yy x_yz_z_xz x_yz_z_xy x_yz_z_xx x_yz_y_zz x_yz_y_yz x_yz_y_yy x_yz_y_xz x_yz_y_xy x_yz_y_xx x_yz_x_zz x_yz_x_yz x_yz_x_yy x_yz_x_xz x_yz_x_xy x_yz_x_xx x_yy_z_zz x_yy_z_yz x_yy_z_yy x_yy_z_xz x_yy_z_xy x_yy_z_xx x_yy_y_zz x_yy_y_yz x_yy_y_yy x_yy_y_xz x_yy_y_xy x_yy_y_xx x_yy_x_zz x_yy_x_yz x_yy_x_yy x_yy_x_xz x_yy_x_xy x_yy_x_xx x_xz_z_zz x_xz_z_yz x_xz_z_yy x_xz_z_xz x_xz_z_xy x_xz_z_xx x_xz_y_zz x_xz_y_yz x_xz_y_yy x_xz_y_xz x_xz_y_xy x_xz_y_xx x_xz_x_zz x_xz_x_yz x_xz_x_yy x_xz_x_xz x_xz_x_xy x_xz_x_xx x_xy_z_zz x_xy_z_yz x_xy_z_yy x_xy_z_xz x_xy_z_xy x_xy_z_xx x_xy_y_zz x_xy_y_yz x_xy_y_yy x_xy_y_xz x_xy_y_xy x_xy_y_xx x_xy_x_zz x_xy_x_yz x_xy_x_yy x_xy_x_xz x_xy_x_xy x_xy_x_xx x_xx_z_zz x_xx_z_yz x_xx_z_yy x_xx_z_xz x_xx_z_xy x_xx_z_xx x_xx_y_zz x_xx_y_yz x_xx_y_yy x_xx_y_xz x_xx_y_xy x_xx_y_xx x_xx_x_zz x_xx_x_yz x_xx_x_yy x_xx_x_xz x_xx_x_xy x_xx_x_xx 
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
0_zz_z_zz_0 0_zz_z_yz_0 0_zz_z_yy_0 0_zz_z_xz_0 0_zz_z_xy_0 0_zz_z_xx_0 0_zz_y_zz_0 0_zz_y_yz_0 0_zz_y_yy_0 0_zz_y_xz_0 0_zz_y_xy_0 0_zz_y_xx_0 0_zz_x_zz_0 0_zz_x_yz_0 0_zz_x_yy_0 0_zz_x_xz_0 0_zz_x_xy_0 0_zz_x_xx_0 0_zzz_z_zz_0 0_zzz_z_yz_0 0_zzz_z_yy_0 0_zzz_z_xz_0 0_zzz_z_xy_0 0_zzz_z_xx_0 0_zzz_y_zz_0 0_zzz_y_yz_0 0_zzz_y_yy_0 0_zzz_y_xz_0 0_zzz_y_xy_0 0_zzz_y_xx_0 0_zzz_x_zz_0 0_zzz_x_yz_0 0_zzz_x_yy_0 0_zzz_x_xz_0 0_zzz_x_xy_0 0_zzz_x_xx_0 0_yz_z_zz_0 0_yz_z_yz_0 0_yz_z_yy_0 0_yz_z_xz_0 0_yz_z_xy_0 0_yz_z_xx_0 0_yz_y_zz_0 0_yz_y_yz_0 0_yz_y_yy_0 0_yz_y_xz_0 0_yz_y_xy_0 0_yz_y_xx_0 0_yz_x_zz_0 0_yz_x_yz_0 0_yz_x_yy_0 0_yz_x_xz_0 0_yz_x_xy_0 0_yz_x_xx_0 0_yzz_z_zz_0 0_yzz_z_yz_0 0_yzz_z_yy_0 0_yzz_z_xz_0 0_yzz_z_xy_0 0_yzz_z_xx_0 0_yzz_y_zz_0 0_yzz_y_yz_0 0_yzz_y_yy_0 0_yzz_y_xz_0 0_yzz_y_xy_0 0_yzz_y_xx_0 0_yzz_x_zz_0 0_yzz_x_yz_0 0_yzz_x_yy_0 0_yzz_x_xz_0 0_yzz_x_xy_0 0_yzz_x_xx_0 0_yy_z_zz_0 0_yy_z_yz_0 0_yy_z_yy_0 0_yy_z_xz_0 0_yy_z_xy_0 0_yy_z_xx_0 0_yy_y_zz_0 0_yy_y_yz_0 0_yy_y_yy_0 0_yy_y_xz_0 0_yy_y_xy_0 0_yy_y_xx_0 0_yy_x_zz_0 0_yy_x_yz_0 0_yy_x_yy_0 0_yy_x_xz_0 0_yy_x_xy_0 0_yy_x_xx_0 0_yyz_z_zz_0 0_yyz_z_yz_0 0_yyz_z_yy_0 0_yyz_z_xz_0 0_yyz_z_xy_0 0_yyz_z_xx_0 0_yyz_y_zz_0 0_yyz_y_yz_0 0_yyz_y_yy_0 0_yyz_y_xz_0 0_yyz_y_xy_0 0_yyz_y_xx_0 0_yyz_x_zz_0 0_yyz_x_yz_0 0_yyz_x_yy_0 0_yyz_x_xz_0 0_yyz_x_xy_0 0_yyz_x_xx_0 0_yyy_z_zz_0 0_yyy_z_yz_0 0_yyy_z_yy_0 0_yyy_z_xz_0 0_yyy_z_xy_0 0_yyy_z_xx_0 0_yyy_y_zz_0 0_yyy_y_yz_0 0_yyy_y_yy_0 0_yyy_y_xz_0 0_yyy_y_xy_0 0_yyy_y_xx_0 0_yyy_x_zz_0 0_yyy_x_yz_0 0_yyy_x_yy_0 0_yyy_x_xz_0 0_yyy_x_xy_0 0_yyy_x_xx_0 0_xz_z_zz_0 0_xz_z_yz_0 0_xz_z_yy_0 0_xz_z_xz_0 0_xz_z_xy_0 0_xz_z_xx_0 0_xz_y_zz_0 0_xz_y_yz_0 0_xz_y_yy_0 0_xz_y_xz_0 0_xz_y_xy_0 0_xz_y_xx_0 0_xz_x_zz_0 0_xz_x_yz_0 0_xz_x_yy_0 0_xz_x_xz_0 0_xz_x_xy_0 0_xz_x_xx_0 0_xzz_z_zz_0 0_xzz_z_yz_0 0_xzz_z_yy_0 0_xzz_z_xz_0 0_xzz_z_xy_0 0_xzz_z_xx_0 0_xzz_y_zz_0 0_xzz_y_yz_0 0_xzz_y_yy_0 0_xzz_y_xz_0 0_xzz_y_xy_0 0_xzz_y_xx_0 0_xzz_x_zz_0 0_xzz_x_yz_0 0_xzz_x_yy_0 0_xzz_x_xz_0 0_xzz_x_xy_0 0_xzz_x_xx_0 0_xy_z_zz_0 0_xy_z_yz_0 0_xy_z_yy_0 0_xy_z_xz_0 0_xy_z_xy_0 0_xy_z_xx_0 0_xy_y_zz_0 0_xy_y_yz_0 0_xy_y_yy_0 0_xy_y_xz_0 0_xy_y_xy_0 0_xy_y_xx_0 0_xy_x_zz_0 0_xy_x_yz_0 0_xy_x_yy_0 0_xy_x_xz_0 0_xy_x_xy_0 0_xy_x_xx_0 0_xyz_z_zz_0 0_xyz_z_yz_0 0_xyz_z_yy_0 0_xyz_z_xz_0 0_xyz_z_xy_0 0_xyz_z_xx_0 0_xyz_y_zz_0 0_xyz_y_yz_0 0_xyz_y_yy_0 0_xyz_y_xz_0 0_xyz_y_xy_0 0_xyz_y_xx_0 0_xyz_x_zz_0 0_xyz_x_yz_0 0_xyz_x_yy_0 0_xyz_x_xz_0 0_xyz_x_xy_0 0_xyz_x_xx_0 0_xyy_z_zz_0 0_xyy_z_yz_0 0_xyy_z_yy_0 0_xyy_z_xz_0 0_xyy_z_xy_0 0_xyy_z_xx_0 0_xyy_y_zz_0 0_xyy_y_yz_0 0_xyy_y_yy_0 0_xyy_y_xz_0 0_xyy_y_xy_0 0_xyy_y_xx_0 0_xyy_x_zz_0 0_xyy_x_yz_0 0_xyy_x_yy_0 0_xyy_x_xz_0 0_xyy_x_xy_0 0_xyy_x_xx_0 0_xx_z_zz_0 0_xx_z_yz_0 0_xx_z_yy_0 0_xx_z_xz_0 0_xx_z_xy_0 0_xx_z_xx_0 0_xx_y_zz_0 0_xx_y_yz_0 0_xx_y_yy_0 0_xx_y_xz_0 0_xx_y_xy_0 0_xx_y_xx_0 0_xx_x_zz_0 0_xx_x_yz_0 0_xx_x_yy_0 0_xx_x_xz_0 0_xx_x_xy_0 0_xx_x_xx_0 0_xxz_z_zz_0 0_xxz_z_yz_0 0_xxz_z_yy_0 0_xxz_z_xz_0 0_xxz_z_xy_0 0_xxz_z_xx_0 0_xxz_y_zz_0 0_xxz_y_yz_0 0_xxz_y_yy_0 0_xxz_y_xz_0 0_xxz_y_xy_0 0_xxz_y_xx_0 0_xxz_x_zz_0 0_xxz_x_yz_0 0_xxz_x_yy_0 0_xxz_x_xz_0 0_xxz_x_xy_0 0_xxz_x_xx_0 0_xxy_z_zz_0 0_xxy_z_yz_0 0_xxy_z_yy_0 0_xxy_z_xz_0 0_xxy_z_xy_0 0_xxy_z_xx_0 0_xxy_y_zz_0 0_xxy_y_yz_0 0_xxy_y_yy_0 0_xxy_y_xz_0 0_xxy_y_xy_0 0_xxy_y_xx_0 0_xxy_x_zz_0 0_xxy_x_yz_0 0_xxy_x_yy_0 0_xxy_x_xz_0 0_xxy_x_xy_0 0_xxy_x_xx_0 0_xxx_z_zz_0 0_xxx_z_yz_0 0_xxx_z_yy_0 0_xxx_z_xz_0 0_xxx_z_xy_0 0_xxx_z_xx_0 0_xxx_y_zz_0 0_xxx_y_yz_0 0_xxx_y_yy_0 0_xxx_y_xz_0 0_xxx_y_xy_0 0_xxx_y_xx_0 0_xxx_x_zz_0 0_xxx_x_yz_0 0_xxx_x_yy_0 0_xxx_x_xz_0 0_xxx_x_xy_0 0_xxx_x_xx_0 z_zz_z_zz_0 z_zz_z_yz_0 z_zz_z_yy_0 z_zz_z_xz_0 z_zz_z_xy_0 z_zz_z_xx_0 z_zz_y_zz_0 z_zz_y_yz_0 z_zz_y_yy_0 z_zz_y_xz_0 z_zz_y_xy_0 z_zz_y_xx_0 z_zz_x_zz_0 z_zz_x_yz_0 z_zz_x_yy_0 z_zz_x_xz_0 z_zz_x_xy_0 z_zz_x_xx_0 z_yz_z_zz_0 z_yz_z_yz_0 z_yz_z_yy_0 z_yz_z_xz_0 z_yz_z_xy_0 z_yz_z_xx_0 z_yz_y_zz_0 z_yz_y_yz_0 z_yz_y_yy_0 z_yz_y_xz_0 z_yz_y_xy_0 z_yz_y_xx_0 z_yz_x_zz_0 z_yz_x_yz_0 z_yz_x_yy_0 z_yz_x_xz_0 z_yz_x_xy_0 z_yz_x_xx_0 z_yy_z_zz_0 z_yy_z_yz_0 z_yy_z_yy_0 z_yy_z_xz_0 z_yy_z_xy_0 z_yy_z_xx_0 z_yy_y_zz_0 z_yy_y_yz_0 z_yy_y_yy_0 z_yy_y_xz_0 z_yy_y_xy_0 z_yy_y_xx_0 z_yy_x_zz_0 z_yy_x_yz_0 z_yy_x_yy_0 z_yy_x_xz_0 z_yy_x_xy_0 z_yy_x_xx_0 z_xz_z_zz_0 z_xz_z_yz_0 z_xz_z_yy_0 z_xz_z_xz_0 z_xz_z_xy_0 z_xz_z_xx_0 z_xz_y_zz_0 z_xz_y_yz_0 z_xz_y_yy_0 z_xz_y_xz_0 z_xz_y_xy_0 z_xz_y_xx_0 z_xz_x_zz_0 z_xz_x_yz_0 z_xz_x_yy_0 z_xz_x_xz_0 z_xz_x_xy_0 z_xz_x_xx_0 z_xy_z_zz_0 z_xy_z_yz_0 z_xy_z_yy_0 z_xy_z_xz_0 z_xy_z_xy_0 z_xy_z_xx_0 z_xy_y_zz_0 z_xy_y_yz_0 z_xy_y_yy_0 z_xy_y_xz_0 z_xy_y_xy_0 z_xy_y_xx_0 z_xy_x_zz_0 z_xy_x_yz_0 z_xy_x_yy_0 z_xy_x_xz_0 z_xy_x_xy_0 z_xy_x_xx_0 z_xx_z_zz_0 z_xx_z_yz_0 z_xx_z_yy_0 z_xx_z_xz_0 z_xx_z_xy_0 z_xx_z_xx_0 z_xx_y_zz_0 z_xx_y_yz_0 z_xx_y_yy_0 z_xx_y_xz_0 z_xx_y_xy_0 z_xx_y_xx_0 z_xx_x_zz_0 z_xx_x_yz_0 z_xx_x_yy_0 z_xx_x_xz_0 z_xx_x_xy_0 z_xx_x_xx_0 y_zz_z_zz_0 y_zz_z_yz_0 y_zz_z_yy_0 y_zz_z_xz_0 y_zz_z_xy_0 y_zz_z_xx_0 y_zz_y_zz_0 y_zz_y_yz_0 y_zz_y_yy_0 y_zz_y_xz_0 y_zz_y_xy_0 y_zz_y_xx_0 y_zz_x_zz_0 y_zz_x_yz_0 y_zz_x_yy_0 y_zz_x_xz_0 y_zz_x_xy_0 y_zz_x_xx_0 y_yz_z_zz_0 y_yz_z_yz_0 y_yz_z_yy_0 y_yz_z_xz_0 y_yz_z_xy_0 y_yz_z_xx_0 y_yz_y_zz_0 y_yz_y_yz_0 y_yz_y_yy_0 y_yz_y_xz_0 y_yz_y_xy_0 y_yz_y_xx_0 y_yz_x_zz_0 y_yz_x_yz_0 y_yz_x_yy_0 y_yz_x_xz_0 y_yz_x_xy_0 y_yz_x_xx_0 y_yy_z_zz_0 y_yy_z_yz_0 y_yy_z_yy_0 y_yy_z_xz_0 y_yy_z_xy_0 y_yy_z_xx_0 y_yy_y_zz_0 y_yy_y_yz_0 y_yy_y_yy_0 y_yy_y_xz_0 y_yy_y_xy_0 y_yy_y_xx_0 y_yy_x_zz_0 y_yy_x_yz_0 y_yy_x_yy_0 y_yy_x_xz_0 y_yy_x_xy_0 y_yy_x_xx_0 y_xz_z_zz_0 y_xz_z_yz_0 y_xz_z_yy_0 y_xz_z_xz_0 y_xz_z_xy_0 y_xz_z_xx_0 y_xz_y_zz_0 y_xz_y_yz_0 y_xz_y_yy_0 y_xz_y_xz_0 y_xz_y_xy_0 y_xz_y_xx_0 y_xz_x_zz_0 y_xz_x_yz_0 y_xz_x_yy_0 y_xz_x_xz_0 y_xz_x_xy_0 y_xz_x_xx_0 y_xy_z_zz_0 y_xy_z_yz_0 y_xy_z_yy_0 y_xy_z_xz_0 y_xy_z_xy_0 y_xy_z_xx_0 y_xy_y_zz_0 y_xy_y_yz_0 y_xy_y_yy_0 y_xy_y_xz_0 y_xy_y_xy_0 y_xy_y_xx_0 y_xy_x_zz_0 y_xy_x_yz_0 y_xy_x_yy_0 y_xy_x_xz_0 y_xy_x_xy_0 y_xy_x_xx_0 y_xx_z_zz_0 y_xx_z_yz_0 y_xx_z_yy_0 y_xx_z_xz_0 y_xx_z_xy_0 y_xx_z_xx_0 y_xx_y_zz_0 y_xx_y_yz_0 y_xx_y_yy_0 y_xx_y_xz_0 y_xx_y_xy_0 y_xx_y_xx_0 y_xx_x_zz_0 y_xx_x_yz_0 y_xx_x_yy_0 y_xx_x_xz_0 y_xx_x_xy_0 y_xx_x_xx_0 x_zz_z_zz_0 x_zz_z_yz_0 x_zz_z_yy_0 x_zz_z_xz_0 x_zz_z_xy_0 x_zz_z_xx_0 x_zz_y_zz_0 x_zz_y_yz_0 x_zz_y_yy_0 x_zz_y_xz_0 x_zz_y_xy_0 x_zz_y_xx_0 x_zz_x_zz_0 x_zz_x_yz_0 x_zz_x_yy_0 x_zz_x_xz_0 x_zz_x_xy_0 x_zz_x_xx_0 x_yz_z_zz_0 x_yz_z_yz_0 x_yz_z_yy_0 x_yz_z_xz_0 x_yz_z_xy_0 x_yz_z_xx_0 x_yz_y_zz_0 x_yz_y_yz_0 x_yz_y_yy_0 x_yz_y_xz_0 x_yz_y_xy_0 x_yz_y_xx_0 x_yz_x_zz_0 x_yz_x_yz_0 x_yz_x_yy_0 x_yz_x_xz_0 x_yz_x_xy_0 x_yz_x_xx_0 x_yy_z_zz_0 x_yy_z_yz_0 x_yy_z_yy_0 x_yy_z_xz_0 x_yy_z_xy_0 x_yy_z_xx_0 x_yy_y_zz_0 x_yy_y_yz_0 x_yy_y_yy_0 x_yy_y_xz_0 x_yy_y_xy_0 x_yy_y_xx_0 x_yy_x_zz_0 x_yy_x_yz_0 x_yy_x_yy_0 x_yy_x_xz_0 x_yy_x_xy_0 x_yy_x_xx_0 x_xz_z_zz_0 x_xz_z_yz_0 x_xz_z_yy_0 x_xz_z_xz_0 x_xz_z_xy_0 x_xz_z_xx_0 x_xz_y_zz_0 x_xz_y_yz_0 x_xz_y_yy_0 x_xz_y_xz_0 x_xz_y_xy_0 x_xz_y_xx_0 x_xz_x_zz_0 x_xz_x_yz_0 x_xz_x_yy_0 x_xz_x_xz_0 x_xz_x_xy_0 x_xz_x_xx_0 x_xy_z_zz_0 x_xy_z_yz_0 x_xy_z_yy_0 x_xy_z_xz_0 x_xy_z_xy_0 x_xy_z_xx_0 x_xy_y_zz_0 x_xy_y_yz_0 x_xy_y_yy_0 x_xy_y_xz_0 x_xy_y_xy_0 x_xy_y_xx_0 x_xy_x_zz_0 x_xy_x_yz_0 x_xy_x_yy_0 x_xy_x_xz_0 x_xy_x_xy_0 x_xy_x_xx_0 x_xx_z_zz_0 x_xx_z_yz_0 x_xx_z_yy_0 x_xx_z_xz_0 x_xx_z_xy_0 x_xx_z_xx_0 x_xx_y_zz_0 x_xx_y_yz_0 x_xx_y_yy_0 x_xx_y_xz_0 x_xx_y_xy_0 x_xx_y_xx_0 x_xx_x_zz_0 x_xx_x_yz_0 x_xx_x_yy_0 x_xx_x_xz_0 x_xx_x_xy_0 x_xx_x_xx_0 

signature:
0_zzz_0_zz_0 0_zzz_0_zzz_0 0_zzz_0_yz_0 0_zzz_0_yzz_0 0_zzz_0_yy_0 0_zzz_0_yyz_0 0_zzz_0_yyy_0 0_zzz_0_xz_0 0_zzz_0_xzz_0 0_zzz_0_xy_0 0_zzz_0_xyz_0 0_zzz_0_xyy_0 0_zzz_0_xx_0 0_zzz_0_xxz_0 0_zzz_0_xxy_0 0_zzz_0_xxx_0 0_zzz_z_zz_0 0_zzz_z_yz_0 0_zzz_z_yy_0 0_zzz_z_xz_0 0_zzz_z_xy_0 0_zzz_z_xx_0 0_zzz_y_zz_0 0_zzz_y_yz_0 0_zzz_y_yy_0 0_zzz_y_xz_0 0_zzz_y_xy_0 0_zzz_y_xx_0 0_zzz_x_zz_0 0_zzz_x_yz_0 0_zzz_x_yy_0 0_zzz_x_xz_0 0_zzz_x_xy_0 0_zzz_x_xx_0 0_yzz_0_zz_0 0_yzz_0_zzz_0 0_yzz_0_yz_0 0_yzz_0_yzz_0 0_yzz_0_yy_0 0_yzz_0_yyz_0 0_yzz_0_yyy_0 0_yzz_0_xz_0 0_yzz_0_xzz_0 0_yzz_0_xy_0 0_yzz_0_xyz_0 0_yzz_0_xyy_0 0_yzz_0_xx_0 0_yzz_0_xxz_0 0_yzz_0_xxy_0 0_yzz_0_xxx_0 0_yzz_z_zz_0 0_yzz_z_yz_0 0_yzz_z_yy_0 0_yzz_z_xz_0 0_yzz_z_xy_0 0_yzz_z_xx_0 0_yzz_y_zz_0 0_yzz_y_yz_0 0_yzz_y_yy_0 0_yzz_y_xz_0 0_yzz_y_xy_0 0_yzz_y_xx_0 0_yzz_x_zz_0 0_yzz_x_yz_0 0_yzz_x_yy_0 0_yzz_x_xz_0 0_yzz_x_xy_0 0_yzz_x_xx_0 0_yyz_0_zz_0 0_yyz_0_zzz_0 0_yyz_0_yz_0 0_yyz_0_yzz_0 0_yyz_0_yy_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xz_0 0_yyz_0_xzz_0 0_yyz_0_xy_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xx_0 0_yyz_0_xxz_0 0_yyz_0_xxy_0 0_yyz_0_xxx_0 0_yyz_z_zz_0 0_yyz_z_yz_0 0_yyz_z_yy_0 0_yyz_z_xz_0 0_yyz_z_xy_0 0_yyz_z_xx_0 0_yyz_y_zz_0 0_yyz_y_yz_0 0_yyz_y_yy_0 0_yyz_y_xz_0 0_yyz_y_xy_0 0_yyz_y_xx_0 0_yyz_x_zz_0 0_yyz_x_yz_0 0_yyz_x_yy_0 0_yyz_x_xz_0 0_yyz_x_xy_0 0_yyz_x_xx_0 0_yyy_0_zz_0 0_yyy_0_zzz_0 0_yyy_0_yz_0 0_yyy_0_yzz_0 0_yyy_0_yy_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xz_0 0_yyy_0_xzz_0 0_yyy_0_xy_0 0_yyy_0_xyz_0 0_yyy_0_xyy_0 0_yyy_0_xx_0 0_yyy_0_xxz_0 0_yyy_0_xxy_0 0_yyy_0_xxx_0 0_yyy_z_zz_0 0_yyy_z_yz_0 0_yyy_z_yy_0 0_yyy_z_xz_0 0_yyy_z_xy_0 0_yyy_z_xx_0 0_yyy_y_zz_0 0_yyy_y_yz_0 0_yyy_y_yy_0 0_yyy_y_xz_0 0_yyy_y_xy_0 0_yyy_y_xx_0 0_yyy_x_zz_0 0_yyy_x_yz_0 0_yyy_x_yy_0 0_yyy_x_xz_0 0_yyy_x_xy_0 0_yyy_x_xx_0 0_xzz_0_zz_0 0_xzz_0_zzz_0 0_xzz_0_yz_0 0_xzz_0_yzz_0 0_xzz_0_yy_0 0_xzz_0_yyz_0 0_xzz_0_yyy_0 0_xzz_0_xz_0 0_xzz_0_xzz_0 0_xzz_0_xy_0 0_xzz_0_xyz_0 0_xzz_0_xyy_0 0_xzz_0_xx_0 0_xzz_0_xxz_0 0_xzz_0_xxy_0 0_xzz_0_xxx_0 0_xzz_z_zz_0 0_xzz_z_yz_0 0_xzz_z_yy_0 0_xzz_z_xz_0 0_xzz_z_xy_0 0_xzz_z_xx_0 0_xzz_y_zz_0 0_xzz_y_yz_0 0_xzz_y_yy_0 0_xzz_y_xz_0 0_xzz_y_xy_0 0_xzz_y_xx_0 0_xzz_x_zz_0 0_xzz_x_yz_0 0_xzz_x_yy_0 0_xzz_x_xz_0 0_xzz_x_xy_0 0_xzz_x_xx_0 0_xyz_0_zz_0 0_xyz_0_zzz_0 0_xyz_0_yz_0 0_xyz_0_yzz_0 0_xyz_0_yy_0 0_xyz_0_yyz_0 0_xyz_0_yyy_0 0_xyz_0_xz_0 0_xyz_0_xzz_0 0_xyz_0_xy_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xx_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxx_0 0_xyz_z_zz_0 0_xyz_z_yz_0 0_xyz_z_yy_0 0_xyz_z_xz_0 0_xyz_z_xy_0 0_xyz_z_xx_0 0_xyz_y_zz_0 0_xyz_y_yz_0 0_xyz_y_yy_0 0_xyz_y_xz_0 0_xyz_y_xy_0 0_xyz_y_xx_0 0_xyz_x_zz_0 0_xyz_x_yz_0 0_xyz_x_yy_0 0_xyz_x_xz_0 0_xyz_x_xy_0 0_xyz_x_xx_0 0_xyy_0_zz_0 0_xyy_0_zzz_0 0_xyy_0_yz_0 0_xyy_0_yzz_0 0_xyy_0_yy_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xz_0 0_xyy_0_xzz_0 0_xyy_0_xy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xx_0 0_xyy_0_xxz_0 0_xyy_0_xxy_0 0_xyy_0_xxx_0 0_xyy_z_zz_0 0_xyy_z_yz_0 0_xyy_z_yy_0 0_xyy_z_xz_0 0_xyy_z_xy_0 0_xyy_z_xx_0 0_xyy_y_zz_0 0_xyy_y_yz_0 0_xyy_y_yy_0 0_xyy_y_xz_0 0_xyy_y_xy_0 0_xyy_y_xx_0 0_xyy_x_zz_0 0_xyy_x_yz_0 0_xyy_x_yy_0 0_xyy_x_xz_0 0_xyy_x_xy_0 0_xyy_x_xx_0 0_xxz_0_zz_0 0_xxz_0_zzz_0 0_xxz_0_yz_0 0_xxz_0_yzz_0 0_xxz_0_yy_0 0_xxz_0_yyz_0 0_xxz_0_yyy_0 0_xxz_0_xz_0 0_xxz_0_xzz_0 0_xxz_0_xy_0 0_xxz_0_xyz_0 0_xxz_0_xyy_0 0_xxz_0_xx_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxz_z_zz_0 0_xxz_z_yz_0 0_xxz_z_yy_0 0_xxz_z_xz_0 0_xxz_z_xy_0 0_xxz_z_xx_0 0_xxz_y_zz_0 0_xxz_y_yz_0 0_xxz_y_yy_0 0_xxz_y_xz_0 0_xxz_y_xy_0 0_xxz_y_xx_0 0_xxz_x_zz_0 0_xxz_x_yz_0 0_xxz_x_yy_0 0_xxz_x_xz_0 0_xxz_x_xy_0 0_xxz_x_xx_0 0_xxy_0_zz_0 0_xxy_0_zzz_0 0_xxy_0_yz_0 0_xxy_0_yzz_0 0_xxy_0_yy_0 0_xxy_0_yyz_0 0_xxy_0_yyy_0 0_xxy_0_xz_0 0_xxy_0_xzz_0 0_xxy_0_xy_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xx_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxy_z_zz_0 0_xxy_z_yz_0 0_xxy_z_yy_0 0_xxy_z_xz_0 0_xxy_z_xy_0 0_xxy_z_xx_0 0_xxy_y_zz_0 0_xxy_y_yz_0 0_xxy_y_yy_0 0_xxy_y_xz_0 0_xxy_y_xy_0 0_xxy_y_xx_0 0_xxy_x_zz_0 0_xxy_x_yz_0 0_xxy_x_yy_0 0_xxy_x_xz_0 0_xxy_x_xy_0 0_xxy_x_xx_0 0_xxx_0_zz_0 0_xxx_0_zzz_0 0_xxx_0_yz_0 0_xxx_0_yzz_0 0_xxx_0_yy_0 0_xxx_0_yyz_0 0_xxx_0_yyy_0 0_xxx_0_xz_0 0_xxx_0_xzz_0 0_xxx_0_xy_0 0_xxx_0_xyz_0 0_xxx_0_xyy_0 0_xxx_0_xx_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 0_xxx_z_zz_0 0_xxx_z_yz_0 0_xxx_z_yy_0 0_xxx_z_xz_0 0_xxx_z_xy_0 0_xxx_z_xx_0 0_xxx_y_zz_0 0_xxx_y_yz_0 0_xxx_y_yy_0 0_xxx_y_xz_0 0_xxx_y_xy_0 0_xxx_y_xx_0 0_xxx_x_zz_0 0_xxx_x_yz_0 0_xxx_x_yy_0 0_xxx_x_xz_0 0_xxx_x_xy_0 0_xxx_x_xx_0 

signature:
0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_yyy_0 0_z_0_yyy_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xyy_0 0_z_0_xyy_1 0_z_0_xxz_0 0_z_0_xxz_1 0_z_0_xxy_0 0_z_0_xxy_1 0_z_0_xxx_0 0_z_0_xxx_1 0_zz_0_zz_1 0_zz_0_zzz_0 0_zz_0_zzz_1 0_zz_0_yz_1 0_zz_0_yzz_0 0_zz_0_yzz_1 0_zz_0_yy_1 0_zz_0_yyz_0 0_zz_0_yyz_1 0_zz_0_yyy_0 0_zz_0_yyy_1 0_zz_0_xz_1 0_zz_0_xzz_0 0_zz_0_xzz_1 0_zz_0_xy_1 0_zz_0_xyz_0 0_zz_0_xyz_1 0_zz_0_xyy_0 0_zz_0_xyy_1 0_zz_0_xx_1 0_zz_0_xxz_0 0_zz_0_xxz_1 0_zz_0_xxy_0 0_zz_0_xxy_1 0_zz_0_xxx_0 0_zz_0_xxx_1 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_yyz_0 0_zzz_0_yyy_0 0_zzz_0_xzz_0 0_zzz_0_xyz_0 0_zzz_0_xyy_0 0_zzz_0_xxz_0 0_zzz_0_xxy_0 0_zzz_0_xxx_0 0_y_0_zzz_0 0_y_0_zzz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xzz_0 0_y_0_xzz_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xxz_0 0_y_0_xxz_1 0_y_0_xxy_0 0_y_0_xxy_1 0_y_0_xxx_0 0_y_0_xxx_1 0_yz_0_zz_1 0_yz_0_zzz_0 0_yz_0_zzz_1 0_yz_0_yz_1 0_yz_0_yzz_0 0_yz_0_yzz_1 0_yz_0_yy_1 0_yz_0_yyz_0 0_yz_0_yyz_1 0_yz_0_yyy_0 0_yz_0_yyy_1 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyz_1 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_yyy_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yzz_0_xyy_0 0_yzz_0_xxz_0 0_yzz_0_xxy_0 0_yzz_0_xxx_0 0_yy_0_zz_1 0_yy_0_zzz_0 0_yy_0_zzz_1 0_yy_0_yz_1 0_yy_0_yzz_0 0_yy_0_yzz_1 0_yy_0_yy_1 0_yy_0_yyz_0 0_yy_0_yyz_1 0_yy_0_yyy_0 0_yy_0_yyy_1 0_yy_0_xz_1 0_yy_0_xzz_0 0_yy_0_xzz_1 0_yy_0_xy_1 0_yy_0_xyz_0 0_yy_0_xyz_1 0_yy_0_xyy_0 0_yy_0_xyy_1 0_yy_0_xx_1 0_yy_0_xxz_0 0_yy_0_xxz_1 0_yy_0_xxy_0 0_yy_0_xxy_1 0_yy_0_xxx_0 0_yy_0_xxx_1 0_yyz_0_zzz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xzz_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xxz_0 0_yyz_0_xxy_0 0_yyz_0_xxx_0 0_yyy_0_zzz_0 0_yyy_0_yzz_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xzz_0 0_yyy_0_xyz_0 0_yyy_0_xyy_0 0_yyy_0_xxz_0 0_yyy_0_xxy_0 0_yyy_0_xxx_0 0_x_0_zzz_0 0_x_0_zzz_1 0_x_0_yzz_0 0_x_0_yzz_1 0_x_0_yyz_0 0_x_0_yyz_1 0_x_0_yyy_0 0_x_0_yyy_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xz_1 0_xz_0_xzz_0 0_xz_0_xzz_1 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xx_1 0_xz_0_xxz_0 0_xz_0_xxz_1 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xz_0_xxx_1 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_yyz_0 0_xzz_0_yyy_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xyy_0 0_xzz_0_xxz_0 0_xzz_0_xxy_0 0_xzz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xy_1 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xyy_1 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxy_1 0_xy_0_xxx_0 0_xyz_0_zzz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_yyy_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxx_0 0_xyy_0_zzz_0 0_xyy_0_yzz_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xzz_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxz_0 0_xyy_0_xxy_0 0_xyy_0_xxx_0 0_xx_0_zz_1 0_xx_0_zzz_0 0_xx_0_zzz_1 0_xx_0_yz_1 0_xx_0_yzz_0 0_xx_0_yzz_1 0_xx_0_yy_1 0_xx_0_yyz_0 0_xx_0_yyz_1 0_xx_0_yyy_0 0_xx_0_yyy_1 0_xx_0_xz_1 0_xx_0_xzz_0 0_xx_0_xzz_1 0_xx_0_xy_1 0_xx_0_xyz_0 0_xx_0_xyz_1 0_xx_0_xyy_0 0_xx_0_xyy_1 0_xx_0_xx_1 0_xx_0_xxz_0 0_xx_0_xxz_1 0_xx_0_xxy_0 0_xx_0_xxy_1 0_xx_0_xxx_0 0_xx_0_xxx_1 0_xxz_0_zzz_0 0_xxz_0_yzz_0 0_xxz_0_yyz_0 0_xxz_0_yyy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xyy_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_zzz_0 0_xxy_0_yzz_0 0_xxy_0_yyz_0 0_xxy_0_yyy_0 0_xxy_0_xzz_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_zzz_0 0_xxx_0_yzz_0 0_xxx_0_yyz_0 0_xxx_0_yyy_0 0_xxx_0_xzz_0 0_xxx_0_xyz_0 0_xxx_0_xyy_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

signature:
0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_yy_0 0_zz_0_yy_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zz_0_xx_0 0_zz_0_xx_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_y_0_zz_0 0_y_0_zz_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_zz_1 0_yz_0_yz_0 0_yz_0_yz_1 0_yz_0_yy_0 0_yz_0_yy_1 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yy_0_z_1 0_yy_0_zz_0 0_yy_0_zz_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yy_0_xx_0 0_yy_0_xx_1 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_x_0_zz_0 0_x_0_zz_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xz_1 0_xz_0_xy_0 0_xz_0_xx_0 0_xz_0_xx_1 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xy_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xx_0_z_1 0_xx_0_zz_0 0_xx_0_zz_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_yy_0 0_xx_0_yy_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

signature:
0_zz_0_zz_0 0_zz_0_zzz_0 0_zz_0_yz_0 0_zz_0_yzz_0 0_zz_0_yy_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xz_0 0_zz_0_xzz_0 0_zz_0_xy_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xx_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_zz_z_zz_0 0_zz_z_yz_0 0_zz_z_yy_0 0_zz_z_xz_0 0_zz_z_xy_0 0_zz_z_xx_0 0_zz_y_zz_0 0_zz_y_yz_0 0_zz_y_yy_0 0_zz_y_xz_0 0_zz_y_xy_0 0_zz_y_xx_0 0_zz_x_zz_0 0_zz_x_yz_0 0_zz_x_yy_0 0_zz_x_xz_0 0_zz_x_xy_0 0_zz_x_xx_0 0_yz_0_zz_0 0_yz_0_zzz_0 0_yz_0_yz_0 0_yz_0_yzz_0 0_yz_0_yy_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xz_0 0_yz_0_xzz_0 0_yz_0_xy_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xx_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yz_z_zz_0 0_yz_z_yz_0 0_yz_z_yy_0 0_yz_z_xz_0 0_yz_z_xy_0 0_yz_z_xx_0 0_yz_y_zz_0 0_yz_y_yz_0 0_yz_y_yy_0 0_yz_y_xz_0 0_yz_y_xy_0 0_yz_y_xx_0 0_yz_x_zz_0 0_yz_x_yz_0 0_yz_x_yy_0 0_yz_x_xz_0 0_yz_x_xy_0 0_yz_x_xx_0 0_yy_0_zz_0 0_yy_0_zzz_0 0_yy_0_yz_0 0_yy_0_yzz_0 0_yy_0_yy_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xz_0 0_yy_0_xzz_0 0_yy_0_xy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xx_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_yy_z_zz_0 0_yy_z_yz_0 0_yy_z_yy_0 0_yy_z_xz_0 0_yy_z_xy_0 0_yy_z_xx_0 0_yy_y_zz_0 0_yy_y_yz_0 0_yy_y_yy_0 0_yy_y_xz_0 0_yy_y_xy_0 0_yy_y_xx_0 0_yy_x_zz_0 0_yy_x_yz_0 0_yy_x_yy_0 0_yy_x_xz_0 0_yy_x_xy_0 0_yy_x_xx_0 0_xz_0_zz_0 0_xz_0_zzz_0 0_xz_0_yz_0 0_xz_0_yzz_0 0_xz_0_yy_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xz_0 0_xz_0_xzz_0 0_xz_0_xy_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xx_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xz_z_zz_0 0_xz_z_yz_0 0_xz_z_yy_0 0_xz_z_xz_0 0_xz_z_xy_0 0_xz_z_xx_0 0_xz_y_zz_0 0_xz_y_yz_0 0_xz_y_yy_0 0_xz_y_xz_0 0_xz_y_xy_0 0_xz_y_xx_0 0_xz_x_zz_0 0_xz_x_yz_0 0_xz_x_yy_0 0_xz_x_xz_0 0_xz_x_xy_0 0_xz_x_xx_0 0_xy_0_zz_0 0_xy_0_zzz_0 0_xy_0_yz_0 0_xy_0_yzz_0 0_xy_0_yy_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xz_0 0_xy_0_xzz_0 0_xy_0_xy_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xx_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xy_z_zz_0 0_xy_z_yz_0 0_xy_z_yy_0 0_xy_z_xz_0 0_xy_z_xy_0 0_xy_z_xx_0 0_xy_y_zz_0 0_xy_y_yz_0 0_xy_y_yy_0 0_xy_y_xz_0 0_xy_y_xy_0 0_xy_y_xx_0 0_xy_x_zz_0 0_xy_x_yz_0 0_xy_x_yy_0 0_xy_x_xz_0 0_xy_x_xy_0 0_xy_x_xx_0 0_xx_0_zz_0 0_xx_0_zzz_0 0_xx_0_yz_0 0_xx_0_yzz_0 0_xx_0_yy_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xz_0 0_xx_0_xzz_0 0_xx_0_xy_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xx_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 0_xx_z_zz_0 0_xx_z_yz_0 0_xx_z_yy_0 0_xx_z_xz_0 0_xx_z_xy_0 0_xx_z_xx_0 0_xx_y_zz_0 0_xx_y_yz_0 0_xx_y_yy_0 0_xx_y_xz_0 0_xx_y_xy_0 0_xx_y_xx_0 0_xx_x_zz_0 0_xx_x_yz_0 0_xx_x_yy_0 0_xx_x_xz_0 0_xx_x_xy_0 0_xx_x_xx_0 

signature:
0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_yyy_0 0_z_0_yyy_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xyy_0 0_z_0_xyy_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_z_0_xxy_0 0_z_0_xxy_1 0_z_0_xxx_0 0_z_0_xxx_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_y_0_zz_1 0_y_0_zzz_0 0_y_0_zzz_1 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xzz_0 0_y_0_xzz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxz_0 0_y_0_xxz_1 0_y_0_xxy_0 0_y_0_xxy_1 0_y_0_xxx_0 0_y_0_xxx_1 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xyz_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_x_0_zz_1 0_x_0_zzz_0 0_x_0_zzz_1 0_x_0_yz_1 0_x_0_yzz_0 0_x_0_yzz_1 0_x_0_yy_1 0_x_0_yyz_0 0_x_0_yyz_1 0_x_0_yyy_0 0_x_0_yyy_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_xzz_0 0_xz_0_xxz_0 0_xz_0_xxx_0 0_xy_0_xyy_0 0_xy_0_xxy_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

signature:
0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_yyy_0 0_z_0_yyy_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xyy_0 0_z_0_xyy_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_z_0_xxy_0 0_z_0_xxy_1 0_z_0_xxx_0 0_z_0_xxx_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_y_0_zz_1 0_y_0_zzz_0 0_y_0_zzz_1 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xzz_0 0_y_0_xzz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxz_0 0_y_0_xxz_1 0_y_0_xxy_0 0_y_0_xxy_1 0_y_0_xxx_0 0_y_0_xxx_1 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_x_0_zz_1 0_x_0_zzz_0 0_x_0_zzz_1 0_x_0_yz_1 0_x_0_yzz_0 0_x_0_yzz_1 0_x_0_yy_1 0_x_0_yyz_0 0_x_0_yyz_1 0_x_0_yyy_0 0_x_0_yyy_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xz_0_xx_0 0_xy_0_xy_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_0_1 0_z_0_z_0 0_z_0_z_1 0_z_0_y_0 0_z_0_y_1 0_z_0_x_0 0_z_0_x_1 0_zz_0_z_0 0_zz_0_y_0 0_zz_0_x_0 0_y_0_0_1 0_y_0_z_0 0_y_0_z_1 0_y_0_y_0 0_y_0_y_1 0_y_0_x_0 0_y_0_x_1 0_yy_0_z_0 0_yy_0_y_0 0_yy_0_x_0 0_x_0_0_1 0_x_0_z_0 0_x_0_z_1 0_x_0_y_0 0_x_0_y_1 0_x_0_x_0 0_x_0_x_1 0_xx_0_z_0 0_xx_0_y_0 0_xx_0_x_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

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
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xx_0 

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

INTEGRAL:0 : 1 : 5 : SSSS_5 Y SSSS_6 Y SSSP_5 Y 

0_0_0_0_5 0_0_0_0_6 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSS_5 SSSS_6 SSSP_5 

SSSS_5 : 0_0_0_0_5 

SSSS_6 : 0_0_0_0_6 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

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

INTEGRAL:0 : 2 : 4 : SSSS_4 Y SSSS_5 Y SSSP_4 Y SSSP_5 Y SSSD_4 Y 

0_0_0_0_4 0_0_0_0_5 0_0_0_z_4 0_0_0_z_5 0_0_0_zz_4 0_0_0_y_4 0_0_0_y_5 0_0_0_yz_4 0_0_0_yy_4 0_0_0_x_4 0_0_0_x_5 0_0_0_xx_4 

SSSS_4 SSSS_5 SSSP_4 SSSP_5 SSSD_4 

SSSS_4 : 0_0_0_0_4 

SSSS_5 : 0_0_0_0_5 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xx_4 

INTEGRAL:0 : 3 : 0 : SSSP_0 Y SSSP_1 Y SSSD_0 N SSSD_1 N SSSF_0 Y 

0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSP_0 SSSP_1 SSSD_0 SSSD_1 SSSF_0 

SSSP_0 : 0_0_0_z_0 0_0_0_y_0 0_0_0_x_0 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

INTEGRAL:0 : 3 : 1 : SSSP_1 Y SSSP_2 Y SSSD_1 N SSSD_2 N SSSF_1 Y 

0_0_0_z_1 0_0_0_z_2 0_0_0_zz_1 0_0_0_zz_2 0_0_0_zzz_1 0_0_0_y_1 0_0_0_y_2 0_0_0_yz_1 0_0_0_yz_2 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yy_2 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_x_1 0_0_0_x_2 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xx_2 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSP_1 SSSP_2 SSSD_1 SSSD_2 SSSF_1 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

INTEGRAL:0 : 3 : 2 : SSSP_2 Y SSSP_3 Y SSSD_2 N SSSD_3 N SSSF_2 Y 

0_0_0_z_2 0_0_0_z_3 0_0_0_zz_2 0_0_0_zz_3 0_0_0_zzz_2 0_0_0_y_2 0_0_0_y_3 0_0_0_yz_2 0_0_0_yz_3 0_0_0_yzz_2 0_0_0_yy_2 0_0_0_yy_3 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_x_2 0_0_0_x_3 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xx_2 0_0_0_xx_3 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSP_2 SSSP_3 SSSD_2 SSSD_3 SSSF_2 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

INTEGRAL:0 : 3 : 3 : SSSP_3 Y SSSP_4 Y SSSD_3 N SSSD_4 Y SSSF_3 Y 

0_0_0_z_3 0_0_0_z_4 0_0_0_zz_3 0_0_0_zz_4 0_0_0_zzz_3 0_0_0_y_3 0_0_0_y_4 0_0_0_yz_3 0_0_0_yz_4 0_0_0_yzz_3 0_0_0_yy_3 0_0_0_yy_4 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_x_3 0_0_0_x_4 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xx_3 0_0_0_xx_4 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SSSP_3 SSSP_4 SSSD_3 SSSD_4 SSSF_3 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xx_4 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

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

INTEGRAL:1 : 3 : 0 : SSSD_1 Y SSSF_0 Y SSSF_1 Y SPSF_0 Y 

0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

SSSD_1 SSSF_0 SSSF_1 SPSF_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SPSF_0 : 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

INTEGRAL:1 : 3 : 1 : SSSD_2 Y SSSF_1 Y SSSF_2 Y SPSF_1 Y 

0_0_0_zz_2 0_0_0_zzz_1 0_0_0_zzz_2 0_0_0_yz_2 0_0_0_yzz_1 0_0_0_yzz_2 0_0_0_yy_2 0_0_0_yyz_1 0_0_0_yyz_2 0_0_0_yyy_1 0_0_0_yyy_2 0_0_0_xz_2 0_0_0_xzz_1 0_0_0_xzz_2 0_0_0_xy_2 0_0_0_xyz_1 0_0_0_xyz_2 0_0_0_xyy_1 0_0_0_xyy_2 0_0_0_xx_2 0_0_0_xxz_1 0_0_0_xxz_2 0_0_0_xxy_1 0_0_0_xxy_2 0_0_0_xxx_1 0_0_0_xxx_2 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SSSD_2 SSSF_1 SSSF_2 SPSF_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

INTEGRAL:1 : 3 : 2 : SSSD_3 Y SSSF_2 Y SSSF_3 Y SPSF_2 Y 

0_0_0_zz_3 0_0_0_zzz_2 0_0_0_zzz_3 0_0_0_yz_3 0_0_0_yzz_2 0_0_0_yzz_3 0_0_0_yy_3 0_0_0_yyz_2 0_0_0_yyz_3 0_0_0_yyy_2 0_0_0_yyy_3 0_0_0_xz_3 0_0_0_xzz_2 0_0_0_xzz_3 0_0_0_xy_3 0_0_0_xyz_2 0_0_0_xyz_3 0_0_0_xyy_2 0_0_0_xyy_3 0_0_0_xx_3 0_0_0_xxz_2 0_0_0_xxz_3 0_0_0_xxy_2 0_0_0_xxy_3 0_0_0_xxx_2 0_0_0_xxx_3 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_yyy_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xyy_2 0_z_0_xxz_2 0_z_0_xxy_2 0_z_0_xxx_2 0_y_0_zzz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xzz_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxz_2 0_y_0_xxy_2 0_y_0_xxx_2 0_x_0_zzz_2 0_x_0_yzz_2 0_x_0_yyz_2 0_x_0_yyy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

SSSD_3 SSSF_2 SSSF_3 SPSF_2 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SPSF_2 : 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_yyy_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xyy_2 0_z_0_xxz_2 0_z_0_xxy_2 0_z_0_xxx_2 0_y_0_zzz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xzz_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxz_2 0_y_0_xxy_2 0_y_0_xxx_2 0_x_0_zzz_2 0_x_0_yzz_2 0_x_0_yyz_2 0_x_0_yyy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

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

INTEGRAL:2 : 3 : 0 : SSSF_0 Y SSSF_1 Y SPSD_1 Y SPSF_0 Y SPSF_1 Y SDSF_0 Y 

0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_yyy_0 0_z_0_yyy_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xyy_0 0_z_0_xyy_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_z_0_xxy_0 0_z_0_xxy_1 0_z_0_xxx_0 0_z_0_xxx_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_y_0_zz_1 0_y_0_zzz_0 0_y_0_zzz_1 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xzz_0 0_y_0_xzz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxz_0 0_y_0_xxz_1 0_y_0_xxy_0 0_y_0_xxy_1 0_y_0_xxx_0 0_y_0_xxx_1 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_x_0_zz_1 0_x_0_zzz_0 0_x_0_zzz_1 0_x_0_yz_1 0_x_0_yzz_0 0_x_0_yzz_1 0_x_0_yy_1 0_x_0_yyz_0 0_x_0_yyz_1 0_x_0_yyy_0 0_x_0_yyy_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

SSSF_0 SSSF_1 SPSD_1 SPSF_0 SPSF_1 SDSF_0 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SPSF_0 : 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SDSF_0 : 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

INTEGRAL:2 : 3 : 1 : SSSF_1 Y SSSF_2 Y SPSD_2 Y SPSF_1 Y SPSF_2 Y SDSF_1 Y 

0_0_0_zzz_1 0_0_0_zzz_2 0_0_0_yzz_1 0_0_0_yzz_2 0_0_0_yyz_1 0_0_0_yyz_2 0_0_0_yyy_1 0_0_0_yyy_2 0_0_0_xzz_1 0_0_0_xzz_2 0_0_0_xyz_1 0_0_0_xyz_2 0_0_0_xyy_1 0_0_0_xyy_2 0_0_0_xxz_1 0_0_0_xxz_2 0_0_0_xxy_1 0_0_0_xxy_2 0_0_0_xxx_1 0_0_0_xxx_2 0_z_0_zz_2 0_z_0_zzz_1 0_z_0_zzz_2 0_z_0_yz_2 0_z_0_yzz_1 0_z_0_yzz_2 0_z_0_yy_2 0_z_0_yyz_1 0_z_0_yyz_2 0_z_0_yyy_1 0_z_0_yyy_2 0_z_0_xz_2 0_z_0_xzz_1 0_z_0_xzz_2 0_z_0_xy_2 0_z_0_xyz_1 0_z_0_xyz_2 0_z_0_xyy_1 0_z_0_xyy_2 0_z_0_xx_2 0_z_0_xxz_1 0_z_0_xxz_2 0_z_0_xxy_1 0_z_0_xxy_2 0_z_0_xxx_1 0_z_0_xxx_2 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_yyy_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xyy_1 0_zz_0_xxz_1 0_zz_0_xxy_1 0_zz_0_xxx_1 0_y_0_zz_2 0_y_0_zzz_1 0_y_0_zzz_2 0_y_0_yz_2 0_y_0_yzz_1 0_y_0_yzz_2 0_y_0_yy_2 0_y_0_yyz_1 0_y_0_yyz_2 0_y_0_yyy_1 0_y_0_yyy_2 0_y_0_xz_2 0_y_0_xzz_1 0_y_0_xzz_2 0_y_0_xy_2 0_y_0_xyz_1 0_y_0_xyz_2 0_y_0_xyy_1 0_y_0_xyy_2 0_y_0_xx_2 0_y_0_xxz_1 0_y_0_xxz_2 0_y_0_xxy_1 0_y_0_xxy_2 0_y_0_xxx_1 0_y_0_xxx_2 0_yz_0_zzz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_yyy_1 0_yz_0_xyz_1 0_yy_0_zzz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xzz_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxz_1 0_yy_0_xxy_1 0_yy_0_xxx_1 0_x_0_zz_2 0_x_0_zzz_1 0_x_0_zzz_2 0_x_0_yz_2 0_x_0_yzz_1 0_x_0_yzz_2 0_x_0_yy_2 0_x_0_yyz_1 0_x_0_yyz_2 0_x_0_yyy_1 0_x_0_yyy_2 0_x_0_xz_2 0_x_0_xzz_1 0_x_0_xzz_2 0_x_0_xy_2 0_x_0_xyz_1 0_x_0_xyz_2 0_x_0_xyy_1 0_x_0_xyy_2 0_x_0_xx_2 0_x_0_xxz_1 0_x_0_xxz_2 0_x_0_xxy_1 0_x_0_xxy_2 0_x_0_xxx_1 0_x_0_xxx_2 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xz_0_xxx_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_zzz_1 0_xx_0_yzz_1 0_xx_0_yyz_1 0_xx_0_yyy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

SSSF_1 SSSF_2 SPSD_2 SPSF_1 SPSF_2 SDSF_1 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SPSD_2 : 0_z_0_zz_2 0_z_0_yz_2 0_z_0_yy_2 0_z_0_xz_2 0_z_0_xy_2 0_z_0_xx_2 0_y_0_zz_2 0_y_0_yz_2 0_y_0_yy_2 0_y_0_xz_2 0_y_0_xy_2 0_y_0_xx_2 0_x_0_zz_2 0_x_0_yz_2 0_x_0_yy_2 0_x_0_xz_2 0_x_0_xy_2 0_x_0_xx_2 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SPSF_2 : 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_yyy_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xyy_2 0_z_0_xxz_2 0_z_0_xxy_2 0_z_0_xxx_2 0_y_0_zzz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xzz_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxz_2 0_y_0_xxy_2 0_y_0_xxx_2 0_x_0_zzz_2 0_x_0_yzz_2 0_x_0_yyz_2 0_x_0_yyy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

SDSF_1 : 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_yyy_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xyy_1 0_zz_0_xxz_1 0_zz_0_xxy_1 0_zz_0_xxx_1 0_yz_0_zzz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_yyy_1 0_yz_0_xyz_1 0_yy_0_zzz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xzz_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxz_1 0_yy_0_xxy_1 0_yy_0_xxx_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xz_0_xxx_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_zzz_1 0_xx_0_yzz_1 0_xx_0_yyz_1 0_xx_0_yyy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

INTEGRAL:3 : 2 : 0 : SPSD_0 Y SPSD_1 Y SDSP_1 Y SDSD_0 N SDSD_1 Y SFSD_0 Y 

0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_yy_0 0_zz_0_yy_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zz_0_xx_0 0_zz_0_xx_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_y_0_zz_0 0_y_0_zz_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_zz_1 0_yz_0_yz_0 0_yz_0_yz_1 0_yz_0_yy_0 0_yz_0_yy_1 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yy_0_z_1 0_yy_0_zz_0 0_yy_0_zz_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yy_0_xx_0 0_yy_0_xx_1 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_x_0_zz_0 0_x_0_zz_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xz_0_xz_1 0_xz_0_xx_0 0_xz_0_xx_1 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xx_0_z_1 0_xx_0_zz_0 0_xx_0_zz_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_yy_0 0_xx_0_yy_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

SPSD_0 SPSD_1 SDSP_1 SDSD_0 SDSD_1 SFSD_0 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SDSD_0 : 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_yy_1 0_zz_0_xz_1 0_zz_0_xy_1 0_zz_0_xx_1 0_yz_0_zz_1 0_yz_0_yz_1 0_yz_0_yy_1 0_yy_0_zz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xx_1 0_xz_0_xz_1 0_xz_0_xx_1 0_xy_0_xy_1 0_xx_0_zz_1 0_xx_0_yz_1 0_xx_0_yy_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SFSD_0 : 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

INTEGRAL:3 : 3 : 0 : SPSF_0 Y SPSF_1 Y SDSD_1 N SDSF_0 N SDSF_1 Y SFSF_0 Y 

0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_yyy_0 0_z_0_yyy_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xyy_0 0_z_0_xyy_1 0_z_0_xxz_0 0_z_0_xxz_1 0_z_0_xxy_0 0_z_0_xxy_1 0_z_0_xxx_0 0_z_0_xxx_1 0_zz_0_zz_1 0_zz_0_zzz_0 0_zz_0_zzz_1 0_zz_0_yz_1 0_zz_0_yzz_0 0_zz_0_yzz_1 0_zz_0_yy_1 0_zz_0_yyz_0 0_zz_0_yyz_1 0_zz_0_yyy_0 0_zz_0_yyy_1 0_zz_0_xz_1 0_zz_0_xzz_0 0_zz_0_xzz_1 0_zz_0_xy_1 0_zz_0_xyz_0 0_zz_0_xyz_1 0_zz_0_xyy_0 0_zz_0_xyy_1 0_zz_0_xx_1 0_zz_0_xxz_0 0_zz_0_xxz_1 0_zz_0_xxy_0 0_zz_0_xxy_1 0_zz_0_xxx_0 0_zz_0_xxx_1 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_yyz_0 0_zzz_0_yyy_0 0_zzz_0_xzz_0 0_zzz_0_xyz_0 0_zzz_0_xyy_0 0_zzz_0_xxz_0 0_zzz_0_xxy_0 0_zzz_0_xxx_0 0_y_0_zzz_0 0_y_0_zzz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xzz_0 0_y_0_xzz_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xxz_0 0_y_0_xxz_1 0_y_0_xxy_0 0_y_0_xxy_1 0_y_0_xxx_0 0_y_0_xxx_1 0_yz_0_zzz_0 0_yz_0_zzz_1 0_yz_0_yz_1 0_yz_0_yzz_0 0_yz_0_yzz_1 0_yz_0_yyz_0 0_yz_0_yyz_1 0_yz_0_yyy_0 0_yz_0_yyy_1 0_yz_0_xyz_0 0_yz_0_xyz_1 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_yyy_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yzz_0_xyy_0 0_yzz_0_xxz_0 0_yzz_0_xxy_0 0_yzz_0_xxx_0 0_yy_0_zz_1 0_yy_0_zzz_0 0_yy_0_zzz_1 0_yy_0_yz_1 0_yy_0_yzz_0 0_yy_0_yzz_1 0_yy_0_yy_1 0_yy_0_yyz_0 0_yy_0_yyz_1 0_yy_0_yyy_0 0_yy_0_yyy_1 0_yy_0_xz_1 0_yy_0_xzz_0 0_yy_0_xzz_1 0_yy_0_xy_1 0_yy_0_xyz_0 0_yy_0_xyz_1 0_yy_0_xyy_0 0_yy_0_xyy_1 0_yy_0_xx_1 0_yy_0_xxz_0 0_yy_0_xxz_1 0_yy_0_xxy_0 0_yy_0_xxy_1 0_yy_0_xxx_0 0_yy_0_xxx_1 0_yyz_0_zzz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xzz_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xxz_0 0_yyz_0_xxy_0 0_yyz_0_xxx_0 0_yyy_0_zzz_0 0_yyy_0_yzz_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xzz_0 0_yyy_0_xyz_0 0_yyy_0_xyy_0 0_yyy_0_xxz_0 0_yyy_0_xxy_0 0_yyy_0_xxx_0 0_x_0_zzz_0 0_x_0_zzz_1 0_x_0_yzz_0 0_x_0_yzz_1 0_x_0_yyz_0 0_x_0_yyz_1 0_x_0_yyy_0 0_x_0_yyy_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_xzz_0 0_xz_0_xzz_1 0_xz_0_xxz_0 0_xz_0_xxz_1 0_xz_0_xxx_0 0_xz_0_xxx_1 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_yyz_0 0_xzz_0_yyy_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xyy_0 0_xzz_0_xxz_0 0_xzz_0_xxy_0 0_xzz_0_xxx_0 0_xy_0_xyy_0 0_xy_0_xyy_1 0_xy_0_xxy_0 0_xy_0_xxy_1 0_xyz_0_zzz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_yyy_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxx_0 0_xyy_0_zzz_0 0_xyy_0_yzz_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xzz_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxz_0 0_xyy_0_xxy_0 0_xyy_0_xxx_0 0_xx_0_zz_1 0_xx_0_zzz_0 0_xx_0_zzz_1 0_xx_0_yz_1 0_xx_0_yzz_0 0_xx_0_yzz_1 0_xx_0_yy_1 0_xx_0_yyz_0 0_xx_0_yyz_1 0_xx_0_yyy_0 0_xx_0_yyy_1 0_xx_0_xz_1 0_xx_0_xzz_0 0_xx_0_xzz_1 0_xx_0_xy_1 0_xx_0_xyz_0 0_xx_0_xyz_1 0_xx_0_xyy_0 0_xx_0_xyy_1 0_xx_0_xx_1 0_xx_0_xxz_0 0_xx_0_xxz_1 0_xx_0_xxy_0 0_xx_0_xxy_1 0_xx_0_xxx_0 0_xx_0_xxx_1 0_xxz_0_zzz_0 0_xxz_0_yzz_0 0_xxz_0_yyz_0 0_xxz_0_yyy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xyy_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_zzz_0 0_xxy_0_yzz_0 0_xxy_0_yyz_0 0_xxy_0_yyy_0 0_xxy_0_xzz_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_zzz_0 0_xxx_0_yzz_0 0_xxx_0_yyz_0 0_xxx_0_yyy_0 0_xxx_0_xzz_0 0_xxx_0_xyz_0 0_xxx_0_xyy_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

SPSF_0 SPSF_1 SDSD_1 SDSF_0 SDSF_1 SFSF_0 

SPSF_0 : 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_yy_1 0_zz_0_xz_1 0_zz_0_xy_1 0_zz_0_xx_1 0_yz_0_zz_1 0_yz_0_yz_1 0_yz_0_yy_1 0_yy_0_zz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xx_1 0_xz_0_xz_1 0_xz_0_xx_1 0_xy_0_xy_1 0_xx_0_zz_1 0_xx_0_yz_1 0_xx_0_yy_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SDSF_0 : 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

SDSF_1 : 0_zz_0_zzz_1 0_zz_0_yzz_1 0_zz_0_yyz_1 0_zz_0_yyy_1 0_zz_0_xzz_1 0_zz_0_xyz_1 0_zz_0_xyy_1 0_zz_0_xxz_1 0_zz_0_xxy_1 0_zz_0_xxx_1 0_yz_0_zzz_1 0_yz_0_yzz_1 0_yz_0_yyz_1 0_yz_0_yyy_1 0_yz_0_xyz_1 0_yy_0_zzz_1 0_yy_0_yzz_1 0_yy_0_yyz_1 0_yy_0_yyy_1 0_yy_0_xzz_1 0_yy_0_xyz_1 0_yy_0_xyy_1 0_yy_0_xxz_1 0_yy_0_xxy_1 0_yy_0_xxx_1 0_xz_0_xzz_1 0_xz_0_xxz_1 0_xz_0_xxx_1 0_xy_0_xyy_1 0_xy_0_xxy_1 0_xx_0_zzz_1 0_xx_0_yzz_1 0_xx_0_yyz_1 0_xx_0_yyy_1 0_xx_0_xzz_1 0_xx_0_xyz_1 0_xx_0_xyy_1 0_xx_0_xxz_1 0_xx_0_xxy_1 0_xx_0_xxx_1 

SFSF_0 : 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_yyz_0 0_zzz_0_yyy_0 0_zzz_0_xzz_0 0_zzz_0_xyz_0 0_zzz_0_xyy_0 0_zzz_0_xxz_0 0_zzz_0_xxy_0 0_zzz_0_xxx_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_yyy_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yzz_0_xyy_0 0_yzz_0_xxz_0 0_yzz_0_xxy_0 0_yzz_0_xxx_0 0_yyz_0_zzz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xzz_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xxz_0 0_yyz_0_xxy_0 0_yyz_0_xxx_0 0_yyy_0_zzz_0 0_yyy_0_yzz_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xzz_0 0_yyy_0_xyz_0 0_yyy_0_xyy_0 0_yyy_0_xxz_0 0_yyy_0_xxy_0 0_yyy_0_xxx_0 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_yyz_0 0_xzz_0_yyy_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xyy_0 0_xzz_0_xxz_0 0_xzz_0_xxy_0 0_xzz_0_xxx_0 0_xyz_0_zzz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_yyy_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxx_0 0_xyy_0_zzz_0 0_xyy_0_yzz_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xzz_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxz_0 0_xyy_0_xxy_0 0_xyy_0_xxx_0 0_xxz_0_zzz_0 0_xxz_0_yzz_0 0_xxz_0_yyz_0 0_xxz_0_yyy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xyy_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_zzz_0 0_xxy_0_yzz_0 0_xxy_0_yyz_0 0_xxy_0_yyy_0 0_xxy_0_xzz_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_zzz_0 0_xxx_0_yzz_0 0_xxx_0_yyz_0 0_xxx_0_yyy_0 0_xxx_0_xzz_0 0_xxx_0_xyz_0 0_xxx_0_xyy_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

        }
    }
}


} // derirec namespace
