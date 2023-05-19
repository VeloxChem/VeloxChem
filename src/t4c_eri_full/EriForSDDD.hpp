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
#include "EriVRRForSSSG.hpp"
#include "EriVRRForSPSP.hpp"
#include "EriVRRForSPSD.hpp"
#include "EriVRRForSPSF.hpp"
#include "EriVRRForSPSG.hpp"
#include "EriVRRForSDSD.hpp"
#include "EriVRRForSDSF.hpp"
#include "EriVRRForSDSG.hpp"
#include "EriHRRForSDPD.hpp"
#include "EriHRRForSDPF.hpp"
#include "EriHRRForSDDD.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostSDDD(      T*                                 intsBuffer,
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
    BufferHostXY<T> pbufSSSD3(4, ncpairs);

0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xx_3 
    BufferHostXY<T> pbufSSSD4(3, ncpairs);

0_0_0_zz_4 0_0_0_yy_4 0_0_0_xx_4 
    BufferHostXY<T> pbufSSSF0(10, ncpairs);

0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 
    BufferHostXY<T> pbufSSSF1(10, ncpairs);

0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 
    BufferHostXY<T> pbufSSSF2(10, ncpairs);

0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 
    BufferHostXY<T> pbufSSSF3(8, ncpairs);

0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxx_3 
    BufferHostXY<T> pbufSSSG0(15, ncpairs);

0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 
    BufferHostXY<T> pbufSSSG1(15, ncpairs);

0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 
    BufferHostXY<T> pbufSSSG2(15, ncpairs);

0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 
    BufferHostXY<T> pbufSPSP1(9, ncpairs);

0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 
    BufferHostXY<T> pbufSPSD0(18, ncpairs);

0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 
    BufferHostXY<T> pbufSPSD1(18, ncpairs);

0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 
    BufferHostXY<T> pbufSPSF0(30, ncpairs);

0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 
    BufferHostXY<T> pbufSPSF1(30, ncpairs);

0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 
    BufferHostXY<T> pbufSPSG0(45, ncpairs);

0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 
    BufferHostXY<T> pbufSPSG1(45, ncpairs);

0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_yyyy_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xyyy_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxyy_1 0_z_0_xxxz_1 0_z_0_xxxy_1 0_z_0_xxxx_1 0_y_0_zzzz_1 0_y_0_yzzz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xzzz_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxzz_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxz_1 0_y_0_xxxy_1 0_y_0_xxxx_1 0_x_0_zzzz_1 0_x_0_yzzz_1 0_x_0_yyzz_1 0_x_0_yyyz_1 0_x_0_yyyy_1 0_x_0_xzzz_1 0_x_0_xyzz_1 0_x_0_xyyz_1 0_x_0_xyyy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 
    BufferHostXY<T> pbufSDSD0(36, ncpairs);

0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 
    BufferHostXY<T> pbufSDSF0(60, ncpairs);

0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 
    BufferHostXY<T> pbufSDSG0(90, ncpairs);

0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzzz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxzz_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSDSD(36, ncpairs);

0_zz_0_zz 0_zz_0_yz 0_zz_0_yy 0_zz_0_xz 0_zz_0_xy 0_zz_0_xx 0_yz_0_zz 0_yz_0_yz 0_yz_0_yy 0_yz_0_xz 0_yz_0_xy 0_yz_0_xx 0_yy_0_zz 0_yy_0_yz 0_yy_0_yy 0_yy_0_xz 0_yy_0_xy 0_yy_0_xx 0_xz_0_zz 0_xz_0_yz 0_xz_0_yy 0_xz_0_xz 0_xz_0_xy 0_xz_0_xx 0_xy_0_zz 0_xy_0_yz 0_xy_0_yy 0_xy_0_xz 0_xy_0_xy 0_xy_0_xx 0_xx_0_zz 0_xx_0_yz 0_xx_0_yy 0_xx_0_xz 0_xx_0_xy 0_xx_0_xx 
    BufferHostXY<T> cbufSDSF(60, ncpairs);

0_zz_0_zzz 0_zz_0_yzz 0_zz_0_yyz 0_zz_0_yyy 0_zz_0_xzz 0_zz_0_xyz 0_zz_0_xyy 0_zz_0_xxz 0_zz_0_xxy 0_zz_0_xxx 0_yz_0_zzz 0_yz_0_yzz 0_yz_0_yyz 0_yz_0_yyy 0_yz_0_xzz 0_yz_0_xyz 0_yz_0_xyy 0_yz_0_xxz 0_yz_0_xxy 0_yz_0_xxx 0_yy_0_zzz 0_yy_0_yzz 0_yy_0_yyz 0_yy_0_yyy 0_yy_0_xzz 0_yy_0_xyz 0_yy_0_xyy 0_yy_0_xxz 0_yy_0_xxy 0_yy_0_xxx 0_xz_0_zzz 0_xz_0_yzz 0_xz_0_yyz 0_xz_0_yyy 0_xz_0_xzz 0_xz_0_xyz 0_xz_0_xyy 0_xz_0_xxz 0_xz_0_xxy 0_xz_0_xxx 0_xy_0_zzz 0_xy_0_yzz 0_xy_0_yyz 0_xy_0_yyy 0_xy_0_xzz 0_xy_0_xyz 0_xy_0_xyy 0_xy_0_xxz 0_xy_0_xxy 0_xy_0_xxx 0_xx_0_zzz 0_xx_0_yzz 0_xx_0_yyz 0_xx_0_yyy 0_xx_0_xzz 0_xx_0_xyz 0_xx_0_xyy 0_xx_0_xxz 0_xx_0_xxy 0_xx_0_xxx 
    BufferHostXY<T> cbufSDSG(90, ncpairs);

0_zz_0_zzzz 0_zz_0_yzzz 0_zz_0_yyzz 0_zz_0_yyyz 0_zz_0_yyyy 0_zz_0_xzzz 0_zz_0_xyzz 0_zz_0_xyyz 0_zz_0_xyyy 0_zz_0_xxzz 0_zz_0_xxyz 0_zz_0_xxyy 0_zz_0_xxxz 0_zz_0_xxxy 0_zz_0_xxxx 0_yz_0_zzzz 0_yz_0_yzzz 0_yz_0_yyzz 0_yz_0_yyyz 0_yz_0_yyyy 0_yz_0_xzzz 0_yz_0_xyzz 0_yz_0_xyyz 0_yz_0_xyyy 0_yz_0_xxzz 0_yz_0_xxyz 0_yz_0_xxyy 0_yz_0_xxxz 0_yz_0_xxxy 0_yz_0_xxxx 0_yy_0_zzzz 0_yy_0_yzzz 0_yy_0_yyzz 0_yy_0_yyyz 0_yy_0_yyyy 0_yy_0_xzzz 0_yy_0_xyzz 0_yy_0_xyyz 0_yy_0_xyyy 0_yy_0_xxzz 0_yy_0_xxyz 0_yy_0_xxyy 0_yy_0_xxxz 0_yy_0_xxxy 0_yy_0_xxxx 0_xz_0_zzzz 0_xz_0_yzzz 0_xz_0_yyzz 0_xz_0_yyyz 0_xz_0_yyyy 0_xz_0_xzzz 0_xz_0_xyzz 0_xz_0_xyyz 0_xz_0_xyyy 0_xz_0_xxzz 0_xz_0_xxyz 0_xz_0_xxyy 0_xz_0_xxxz 0_xz_0_xxxy 0_xz_0_xxxx 0_xy_0_zzzz 0_xy_0_yzzz 0_xy_0_yyzz 0_xy_0_yyyz 0_xy_0_yyyy 0_xy_0_xzzz 0_xy_0_xyzz 0_xy_0_xyyz 0_xy_0_xyyy 0_xy_0_xxzz 0_xy_0_xxyz 0_xy_0_xxyy 0_xy_0_xxxz 0_xy_0_xxxy 0_xy_0_xxxx 0_xx_0_zzzz 0_xx_0_yzzz 0_xx_0_yyzz 0_xx_0_yyyz 0_xx_0_yyyy 0_xx_0_xzzz 0_xx_0_xyzz 0_xx_0_xyyz 0_xx_0_xyyy 0_xx_0_xxzz 0_xx_0_xxyz 0_xx_0_xxyy 0_xx_0_xxxz 0_xx_0_xxxy 0_xx_0_xxxx 
    BufferHostXY<T> cbufSDPD(108, ncpairs);

0_zz_z_zz 0_zz_z_yz 0_zz_z_yy 0_zz_z_xz 0_zz_z_xy 0_zz_z_xx 0_zz_y_zz 0_zz_y_yz 0_zz_y_yy 0_zz_y_xz 0_zz_y_xy 0_zz_y_xx 0_zz_x_zz 0_zz_x_yz 0_zz_x_yy 0_zz_x_xz 0_zz_x_xy 0_zz_x_xx 0_yz_z_zz 0_yz_z_yz 0_yz_z_yy 0_yz_z_xz 0_yz_z_xy 0_yz_z_xx 0_yz_y_zz 0_yz_y_yz 0_yz_y_yy 0_yz_y_xz 0_yz_y_xy 0_yz_y_xx 0_yz_x_zz 0_yz_x_yz 0_yz_x_yy 0_yz_x_xz 0_yz_x_xy 0_yz_x_xx 0_yy_z_zz 0_yy_z_yz 0_yy_z_yy 0_yy_z_xz 0_yy_z_xy 0_yy_z_xx 0_yy_y_zz 0_yy_y_yz 0_yy_y_yy 0_yy_y_xz 0_yy_y_xy 0_yy_y_xx 0_yy_x_zz 0_yy_x_yz 0_yy_x_yy 0_yy_x_xz 0_yy_x_xy 0_yy_x_xx 0_xz_z_zz 0_xz_z_yz 0_xz_z_yy 0_xz_z_xz 0_xz_z_xy 0_xz_z_xx 0_xz_y_zz 0_xz_y_yz 0_xz_y_yy 0_xz_y_xz 0_xz_y_xy 0_xz_y_xx 0_xz_x_zz 0_xz_x_yz 0_xz_x_yy 0_xz_x_xz 0_xz_x_xy 0_xz_x_xx 0_xy_z_zz 0_xy_z_yz 0_xy_z_yy 0_xy_z_xz 0_xy_z_xy 0_xy_z_xx 0_xy_y_zz 0_xy_y_yz 0_xy_y_yy 0_xy_y_xz 0_xy_y_xy 0_xy_y_xx 0_xy_x_zz 0_xy_x_yz 0_xy_x_yy 0_xy_x_xz 0_xy_x_xy 0_xy_x_xx 0_xx_z_zz 0_xx_z_yz 0_xx_z_yy 0_xx_z_xz 0_xx_z_xy 0_xx_z_xx 0_xx_y_zz 0_xx_y_yz 0_xx_y_yy 0_xx_y_xz 0_xx_y_xy 0_xx_y_xx 0_xx_x_zz 0_xx_x_yz 0_xx_x_yy 0_xx_x_xz 0_xx_x_xy 0_xx_x_xx 
    BufferHostXY<T> cbufSDPF(150, ncpairs);

0_zz_z_zzz 0_zz_z_yzz 0_zz_z_yyz 0_zz_z_xzz 0_zz_z_xyz 0_zz_z_xxz 0_zz_y_zzz 0_zz_y_yzz 0_zz_y_yyz 0_zz_y_yyy 0_zz_y_xzz 0_zz_y_xyz 0_zz_y_xyy 0_zz_y_xxz 0_zz_y_xxy 0_zz_x_zzz 0_zz_x_yzz 0_zz_x_yyz 0_zz_x_yyy 0_zz_x_xzz 0_zz_x_xyz 0_zz_x_xyy 0_zz_x_xxz 0_zz_x_xxy 0_zz_x_xxx 0_yz_z_zzz 0_yz_z_yzz 0_yz_z_yyz 0_yz_z_xzz 0_yz_z_xyz 0_yz_z_xxz 0_yz_y_zzz 0_yz_y_yzz 0_yz_y_yyz 0_yz_y_yyy 0_yz_y_xzz 0_yz_y_xyz 0_yz_y_xyy 0_yz_y_xxz 0_yz_y_xxy 0_yz_x_zzz 0_yz_x_yzz 0_yz_x_yyz 0_yz_x_yyy 0_yz_x_xzz 0_yz_x_xyz 0_yz_x_xyy 0_yz_x_xxz 0_yz_x_xxy 0_yz_x_xxx 0_yy_z_zzz 0_yy_z_yzz 0_yy_z_yyz 0_yy_z_xzz 0_yy_z_xyz 0_yy_z_xxz 0_yy_y_zzz 0_yy_y_yzz 0_yy_y_yyz 0_yy_y_yyy 0_yy_y_xzz 0_yy_y_xyz 0_yy_y_xyy 0_yy_y_xxz 0_yy_y_xxy 0_yy_x_zzz 0_yy_x_yzz 0_yy_x_yyz 0_yy_x_yyy 0_yy_x_xzz 0_yy_x_xyz 0_yy_x_xyy 0_yy_x_xxz 0_yy_x_xxy 0_yy_x_xxx 0_xz_z_zzz 0_xz_z_yzz 0_xz_z_yyz 0_xz_z_xzz 0_xz_z_xyz 0_xz_z_xxz 0_xz_y_zzz 0_xz_y_yzz 0_xz_y_yyz 0_xz_y_yyy 0_xz_y_xzz 0_xz_y_xyz 0_xz_y_xyy 0_xz_y_xxz 0_xz_y_xxy 0_xz_x_zzz 0_xz_x_yzz 0_xz_x_yyz 0_xz_x_yyy 0_xz_x_xzz 0_xz_x_xyz 0_xz_x_xyy 0_xz_x_xxz 0_xz_x_xxy 0_xz_x_xxx 0_xy_z_zzz 0_xy_z_yzz 0_xy_z_yyz 0_xy_z_xzz 0_xy_z_xyz 0_xy_z_xxz 0_xy_y_zzz 0_xy_y_yzz 0_xy_y_yyz 0_xy_y_yyy 0_xy_y_xzz 0_xy_y_xyz 0_xy_y_xyy 0_xy_y_xxz 0_xy_y_xxy 0_xy_x_zzz 0_xy_x_yzz 0_xy_x_yyz 0_xy_x_yyy 0_xy_x_xzz 0_xy_x_xyz 0_xy_x_xyy 0_xy_x_xxz 0_xy_x_xxy 0_xy_x_xxx 0_xx_z_zzz 0_xx_z_yzz 0_xx_z_yyz 0_xx_z_xzz 0_xx_z_xyz 0_xx_z_xxz 0_xx_y_zzz 0_xx_y_yzz 0_xx_y_yyz 0_xx_y_yyy 0_xx_y_xzz 0_xx_y_xyz 0_xx_y_xyy 0_xx_y_xxz 0_xx_y_xxy 0_xx_x_zzz 0_xx_x_yzz 0_xx_x_yyz 0_xx_x_yyy 0_xx_x_xzz 0_xx_x_xyz 0_xx_x_xyy 0_xx_x_xxz 0_xx_x_xxy 0_xx_x_xxx 
    BufferHostXY<T> cbufSDDD(216, ncpairs);

0_zz_zz_zz 0_zz_zz_yz 0_zz_zz_yy 0_zz_zz_xz 0_zz_zz_xy 0_zz_zz_xx 0_zz_yz_zz 0_zz_yz_yz 0_zz_yz_yy 0_zz_yz_xz 0_zz_yz_xy 0_zz_yz_xx 0_zz_yy_zz 0_zz_yy_yz 0_zz_yy_yy 0_zz_yy_xz 0_zz_yy_xy 0_zz_yy_xx 0_zz_xz_zz 0_zz_xz_yz 0_zz_xz_yy 0_zz_xz_xz 0_zz_xz_xy 0_zz_xz_xx 0_zz_xy_zz 0_zz_xy_yz 0_zz_xy_yy 0_zz_xy_xz 0_zz_xy_xy 0_zz_xy_xx 0_zz_xx_zz 0_zz_xx_yz 0_zz_xx_yy 0_zz_xx_xz 0_zz_xx_xy 0_zz_xx_xx 0_yz_zz_zz 0_yz_zz_yz 0_yz_zz_yy 0_yz_zz_xz 0_yz_zz_xy 0_yz_zz_xx 0_yz_yz_zz 0_yz_yz_yz 0_yz_yz_yy 0_yz_yz_xz 0_yz_yz_xy 0_yz_yz_xx 0_yz_yy_zz 0_yz_yy_yz 0_yz_yy_yy 0_yz_yy_xz 0_yz_yy_xy 0_yz_yy_xx 0_yz_xz_zz 0_yz_xz_yz 0_yz_xz_yy 0_yz_xz_xz 0_yz_xz_xy 0_yz_xz_xx 0_yz_xy_zz 0_yz_xy_yz 0_yz_xy_yy 0_yz_xy_xz 0_yz_xy_xy 0_yz_xy_xx 0_yz_xx_zz 0_yz_xx_yz 0_yz_xx_yy 0_yz_xx_xz 0_yz_xx_xy 0_yz_xx_xx 0_yy_zz_zz 0_yy_zz_yz 0_yy_zz_yy 0_yy_zz_xz 0_yy_zz_xy 0_yy_zz_xx 0_yy_yz_zz 0_yy_yz_yz 0_yy_yz_yy 0_yy_yz_xz 0_yy_yz_xy 0_yy_yz_xx 0_yy_yy_zz 0_yy_yy_yz 0_yy_yy_yy 0_yy_yy_xz 0_yy_yy_xy 0_yy_yy_xx 0_yy_xz_zz 0_yy_xz_yz 0_yy_xz_yy 0_yy_xz_xz 0_yy_xz_xy 0_yy_xz_xx 0_yy_xy_zz 0_yy_xy_yz 0_yy_xy_yy 0_yy_xy_xz 0_yy_xy_xy 0_yy_xy_xx 0_yy_xx_zz 0_yy_xx_yz 0_yy_xx_yy 0_yy_xx_xz 0_yy_xx_xy 0_yy_xx_xx 0_xz_zz_zz 0_xz_zz_yz 0_xz_zz_yy 0_xz_zz_xz 0_xz_zz_xy 0_xz_zz_xx 0_xz_yz_zz 0_xz_yz_yz 0_xz_yz_yy 0_xz_yz_xz 0_xz_yz_xy 0_xz_yz_xx 0_xz_yy_zz 0_xz_yy_yz 0_xz_yy_yy 0_xz_yy_xz 0_xz_yy_xy 0_xz_yy_xx 0_xz_xz_zz 0_xz_xz_yz 0_xz_xz_yy 0_xz_xz_xz 0_xz_xz_xy 0_xz_xz_xx 0_xz_xy_zz 0_xz_xy_yz 0_xz_xy_yy 0_xz_xy_xz 0_xz_xy_xy 0_xz_xy_xx 0_xz_xx_zz 0_xz_xx_yz 0_xz_xx_yy 0_xz_xx_xz 0_xz_xx_xy 0_xz_xx_xx 0_xy_zz_zz 0_xy_zz_yz 0_xy_zz_yy 0_xy_zz_xz 0_xy_zz_xy 0_xy_zz_xx 0_xy_yz_zz 0_xy_yz_yz 0_xy_yz_yy 0_xy_yz_xz 0_xy_yz_xy 0_xy_yz_xx 0_xy_yy_zz 0_xy_yy_yz 0_xy_yy_yy 0_xy_yy_xz 0_xy_yy_xy 0_xy_yy_xx 0_xy_xz_zz 0_xy_xz_yz 0_xy_xz_yy 0_xy_xz_xz 0_xy_xz_xy 0_xy_xz_xx 0_xy_xy_zz 0_xy_xy_yz 0_xy_xy_yy 0_xy_xy_xz 0_xy_xy_xy 0_xy_xy_xx 0_xy_xx_zz 0_xy_xx_yz 0_xy_xx_yy 0_xy_xx_xz 0_xy_xx_xy 0_xy_xx_xx 0_xx_zz_zz 0_xx_zz_yz 0_xx_zz_yy 0_xx_zz_xz 0_xx_zz_xy 0_xx_zz_xx 0_xx_yz_zz 0_xx_yz_yz 0_xx_yz_yy 0_xx_yz_xz 0_xx_yz_xy 0_xx_yz_xx 0_xx_yy_zz 0_xx_yy_yz 0_xx_yy_yy 0_xx_yy_xz 0_xx_yy_xy 0_xx_yy_xx 0_xx_xz_zz 0_xx_xz_yz 0_xx_xz_yy 0_xx_xz_xz 0_xx_xz_xy 0_xx_xz_xx 0_xx_xy_zz 0_xx_xy_yz 0_xx_xy_yy 0_xx_xy_xz 0_xx_xy_xy 0_xx_xy_xx 0_xx_xx_zz 0_xx_xx_yz 0_xx_xx_yy 0_xx_xx_xz 0_xx_xx_xy 0_xx_xx_xx 
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
0_zz_z_zz_0 0_zz_z_zzz_0 0_zz_z_yz_0 0_zz_z_yzz_0 0_zz_z_yy_0 0_zz_z_yyz_0 0_zz_z_xz_0 0_zz_z_xzz_0 0_zz_z_xy_0 0_zz_z_xyz_0 0_zz_z_xx_0 0_zz_z_xxz_0 0_zz_zz_zz_0 0_zz_zz_yz_0 0_zz_zz_yy_0 0_zz_zz_xz_0 0_zz_zz_xy_0 0_zz_zz_xx_0 0_zz_y_zz_0 0_zz_y_zzz_0 0_zz_y_yz_0 0_zz_y_yzz_0 0_zz_y_yy_0 0_zz_y_yyz_0 0_zz_y_yyy_0 0_zz_y_xz_0 0_zz_y_xzz_0 0_zz_y_xy_0 0_zz_y_xyz_0 0_zz_y_xyy_0 0_zz_y_xx_0 0_zz_y_xxz_0 0_zz_y_xxy_0 0_zz_yz_zz_0 0_zz_yz_yz_0 0_zz_yz_yy_0 0_zz_yz_xz_0 0_zz_yz_xy_0 0_zz_yz_xx_0 0_zz_yy_zz_0 0_zz_yy_yz_0 0_zz_yy_yy_0 0_zz_yy_xz_0 0_zz_yy_xy_0 0_zz_yy_xx_0 0_zz_x_zz_0 0_zz_x_zzz_0 0_zz_x_yz_0 0_zz_x_yzz_0 0_zz_x_yy_0 0_zz_x_yyz_0 0_zz_x_yyy_0 0_zz_x_xz_0 0_zz_x_xzz_0 0_zz_x_xy_0 0_zz_x_xyz_0 0_zz_x_xyy_0 0_zz_x_xx_0 0_zz_x_xxz_0 0_zz_x_xxy_0 0_zz_x_xxx_0 0_zz_xz_zz_0 0_zz_xz_yz_0 0_zz_xz_yy_0 0_zz_xz_xz_0 0_zz_xz_xy_0 0_zz_xz_xx_0 0_zz_xy_zz_0 0_zz_xy_yz_0 0_zz_xy_yy_0 0_zz_xy_xz_0 0_zz_xy_xy_0 0_zz_xy_xx_0 0_zz_xx_zz_0 0_zz_xx_yz_0 0_zz_xx_yy_0 0_zz_xx_xz_0 0_zz_xx_xy_0 0_zz_xx_xx_0 0_yz_z_zz_0 0_yz_z_zzz_0 0_yz_z_yz_0 0_yz_z_yzz_0 0_yz_z_yy_0 0_yz_z_yyz_0 0_yz_z_xz_0 0_yz_z_xzz_0 0_yz_z_xy_0 0_yz_z_xyz_0 0_yz_z_xx_0 0_yz_z_xxz_0 0_yz_zz_zz_0 0_yz_zz_yz_0 0_yz_zz_yy_0 0_yz_zz_xz_0 0_yz_zz_xy_0 0_yz_zz_xx_0 0_yz_y_zz_0 0_yz_y_zzz_0 0_yz_y_yz_0 0_yz_y_yzz_0 0_yz_y_yy_0 0_yz_y_yyz_0 0_yz_y_yyy_0 0_yz_y_xz_0 0_yz_y_xzz_0 0_yz_y_xy_0 0_yz_y_xyz_0 0_yz_y_xyy_0 0_yz_y_xx_0 0_yz_y_xxz_0 0_yz_y_xxy_0 0_yz_yz_zz_0 0_yz_yz_yz_0 0_yz_yz_yy_0 0_yz_yz_xz_0 0_yz_yz_xy_0 0_yz_yz_xx_0 0_yz_yy_zz_0 0_yz_yy_yz_0 0_yz_yy_yy_0 0_yz_yy_xz_0 0_yz_yy_xy_0 0_yz_yy_xx_0 0_yz_x_zz_0 0_yz_x_zzz_0 0_yz_x_yz_0 0_yz_x_yzz_0 0_yz_x_yy_0 0_yz_x_yyz_0 0_yz_x_yyy_0 0_yz_x_xz_0 0_yz_x_xzz_0 0_yz_x_xy_0 0_yz_x_xyz_0 0_yz_x_xyy_0 0_yz_x_xx_0 0_yz_x_xxz_0 0_yz_x_xxy_0 0_yz_x_xxx_0 0_yz_xz_zz_0 0_yz_xz_yz_0 0_yz_xz_yy_0 0_yz_xz_xz_0 0_yz_xz_xy_0 0_yz_xz_xx_0 0_yz_xy_zz_0 0_yz_xy_yz_0 0_yz_xy_yy_0 0_yz_xy_xz_0 0_yz_xy_xy_0 0_yz_xy_xx_0 0_yz_xx_zz_0 0_yz_xx_yz_0 0_yz_xx_yy_0 0_yz_xx_xz_0 0_yz_xx_xy_0 0_yz_xx_xx_0 0_yy_z_zz_0 0_yy_z_zzz_0 0_yy_z_yz_0 0_yy_z_yzz_0 0_yy_z_yy_0 0_yy_z_yyz_0 0_yy_z_xz_0 0_yy_z_xzz_0 0_yy_z_xy_0 0_yy_z_xyz_0 0_yy_z_xx_0 0_yy_z_xxz_0 0_yy_zz_zz_0 0_yy_zz_yz_0 0_yy_zz_yy_0 0_yy_zz_xz_0 0_yy_zz_xy_0 0_yy_zz_xx_0 0_yy_y_zz_0 0_yy_y_zzz_0 0_yy_y_yz_0 0_yy_y_yzz_0 0_yy_y_yy_0 0_yy_y_yyz_0 0_yy_y_yyy_0 0_yy_y_xz_0 0_yy_y_xzz_0 0_yy_y_xy_0 0_yy_y_xyz_0 0_yy_y_xyy_0 0_yy_y_xx_0 0_yy_y_xxz_0 0_yy_y_xxy_0 0_yy_yz_zz_0 0_yy_yz_yz_0 0_yy_yz_yy_0 0_yy_yz_xz_0 0_yy_yz_xy_0 0_yy_yz_xx_0 0_yy_yy_zz_0 0_yy_yy_yz_0 0_yy_yy_yy_0 0_yy_yy_xz_0 0_yy_yy_xy_0 0_yy_yy_xx_0 0_yy_x_zz_0 0_yy_x_zzz_0 0_yy_x_yz_0 0_yy_x_yzz_0 0_yy_x_yy_0 0_yy_x_yyz_0 0_yy_x_yyy_0 0_yy_x_xz_0 0_yy_x_xzz_0 0_yy_x_xy_0 0_yy_x_xyz_0 0_yy_x_xyy_0 0_yy_x_xx_0 0_yy_x_xxz_0 0_yy_x_xxy_0 0_yy_x_xxx_0 0_yy_xz_zz_0 0_yy_xz_yz_0 0_yy_xz_yy_0 0_yy_xz_xz_0 0_yy_xz_xy_0 0_yy_xz_xx_0 0_yy_xy_zz_0 0_yy_xy_yz_0 0_yy_xy_yy_0 0_yy_xy_xz_0 0_yy_xy_xy_0 0_yy_xy_xx_0 0_yy_xx_zz_0 0_yy_xx_yz_0 0_yy_xx_yy_0 0_yy_xx_xz_0 0_yy_xx_xy_0 0_yy_xx_xx_0 0_xz_z_zz_0 0_xz_z_zzz_0 0_xz_z_yz_0 0_xz_z_yzz_0 0_xz_z_yy_0 0_xz_z_yyz_0 0_xz_z_xz_0 0_xz_z_xzz_0 0_xz_z_xy_0 0_xz_z_xyz_0 0_xz_z_xx_0 0_xz_z_xxz_0 0_xz_zz_zz_0 0_xz_zz_yz_0 0_xz_zz_yy_0 0_xz_zz_xz_0 0_xz_zz_xy_0 0_xz_zz_xx_0 0_xz_y_zz_0 0_xz_y_zzz_0 0_xz_y_yz_0 0_xz_y_yzz_0 0_xz_y_yy_0 0_xz_y_yyz_0 0_xz_y_yyy_0 0_xz_y_xz_0 0_xz_y_xzz_0 0_xz_y_xy_0 0_xz_y_xyz_0 0_xz_y_xyy_0 0_xz_y_xx_0 0_xz_y_xxz_0 0_xz_y_xxy_0 0_xz_yz_zz_0 0_xz_yz_yz_0 0_xz_yz_yy_0 0_xz_yz_xz_0 0_xz_yz_xy_0 0_xz_yz_xx_0 0_xz_yy_zz_0 0_xz_yy_yz_0 0_xz_yy_yy_0 0_xz_yy_xz_0 0_xz_yy_xy_0 0_xz_yy_xx_0 0_xz_x_zz_0 0_xz_x_zzz_0 0_xz_x_yz_0 0_xz_x_yzz_0 0_xz_x_yy_0 0_xz_x_yyz_0 0_xz_x_yyy_0 0_xz_x_xz_0 0_xz_x_xzz_0 0_xz_x_xy_0 0_xz_x_xyz_0 0_xz_x_xyy_0 0_xz_x_xx_0 0_xz_x_xxz_0 0_xz_x_xxy_0 0_xz_x_xxx_0 0_xz_xz_zz_0 0_xz_xz_yz_0 0_xz_xz_yy_0 0_xz_xz_xz_0 0_xz_xz_xy_0 0_xz_xz_xx_0 0_xz_xy_zz_0 0_xz_xy_yz_0 0_xz_xy_yy_0 0_xz_xy_xz_0 0_xz_xy_xy_0 0_xz_xy_xx_0 0_xz_xx_zz_0 0_xz_xx_yz_0 0_xz_xx_yy_0 0_xz_xx_xz_0 0_xz_xx_xy_0 0_xz_xx_xx_0 0_xy_z_zz_0 0_xy_z_zzz_0 0_xy_z_yz_0 0_xy_z_yzz_0 0_xy_z_yy_0 0_xy_z_yyz_0 0_xy_z_xz_0 0_xy_z_xzz_0 0_xy_z_xy_0 0_xy_z_xyz_0 0_xy_z_xx_0 0_xy_z_xxz_0 0_xy_zz_zz_0 0_xy_zz_yz_0 0_xy_zz_yy_0 0_xy_zz_xz_0 0_xy_zz_xy_0 0_xy_zz_xx_0 0_xy_y_zz_0 0_xy_y_zzz_0 0_xy_y_yz_0 0_xy_y_yzz_0 0_xy_y_yy_0 0_xy_y_yyz_0 0_xy_y_yyy_0 0_xy_y_xz_0 0_xy_y_xzz_0 0_xy_y_xy_0 0_xy_y_xyz_0 0_xy_y_xyy_0 0_xy_y_xx_0 0_xy_y_xxz_0 0_xy_y_xxy_0 0_xy_yz_zz_0 0_xy_yz_yz_0 0_xy_yz_yy_0 0_xy_yz_xz_0 0_xy_yz_xy_0 0_xy_yz_xx_0 0_xy_yy_zz_0 0_xy_yy_yz_0 0_xy_yy_yy_0 0_xy_yy_xz_0 0_xy_yy_xy_0 0_xy_yy_xx_0 0_xy_x_zz_0 0_xy_x_zzz_0 0_xy_x_yz_0 0_xy_x_yzz_0 0_xy_x_yy_0 0_xy_x_yyz_0 0_xy_x_yyy_0 0_xy_x_xz_0 0_xy_x_xzz_0 0_xy_x_xy_0 0_xy_x_xyz_0 0_xy_x_xyy_0 0_xy_x_xx_0 0_xy_x_xxz_0 0_xy_x_xxy_0 0_xy_x_xxx_0 0_xy_xz_zz_0 0_xy_xz_yz_0 0_xy_xz_yy_0 0_xy_xz_xz_0 0_xy_xz_xy_0 0_xy_xz_xx_0 0_xy_xy_zz_0 0_xy_xy_yz_0 0_xy_xy_yy_0 0_xy_xy_xz_0 0_xy_xy_xy_0 0_xy_xy_xx_0 0_xy_xx_zz_0 0_xy_xx_yz_0 0_xy_xx_yy_0 0_xy_xx_xz_0 0_xy_xx_xy_0 0_xy_xx_xx_0 0_xx_z_zz_0 0_xx_z_zzz_0 0_xx_z_yz_0 0_xx_z_yzz_0 0_xx_z_yy_0 0_xx_z_yyz_0 0_xx_z_xz_0 0_xx_z_xzz_0 0_xx_z_xy_0 0_xx_z_xyz_0 0_xx_z_xx_0 0_xx_z_xxz_0 0_xx_zz_zz_0 0_xx_zz_yz_0 0_xx_zz_yy_0 0_xx_zz_xz_0 0_xx_zz_xy_0 0_xx_zz_xx_0 0_xx_y_zz_0 0_xx_y_zzz_0 0_xx_y_yz_0 0_xx_y_yzz_0 0_xx_y_yy_0 0_xx_y_yyz_0 0_xx_y_yyy_0 0_xx_y_xz_0 0_xx_y_xzz_0 0_xx_y_xy_0 0_xx_y_xyz_0 0_xx_y_xyy_0 0_xx_y_xx_0 0_xx_y_xxz_0 0_xx_y_xxy_0 0_xx_yz_zz_0 0_xx_yz_yz_0 0_xx_yz_yy_0 0_xx_yz_xz_0 0_xx_yz_xy_0 0_xx_yz_xx_0 0_xx_yy_zz_0 0_xx_yy_yz_0 0_xx_yy_yy_0 0_xx_yy_xz_0 0_xx_yy_xy_0 0_xx_yy_xx_0 0_xx_x_zz_0 0_xx_x_zzz_0 0_xx_x_yz_0 0_xx_x_yzz_0 0_xx_x_yy_0 0_xx_x_yyz_0 0_xx_x_yyy_0 0_xx_x_xz_0 0_xx_x_xzz_0 0_xx_x_xy_0 0_xx_x_xyz_0 0_xx_x_xyy_0 0_xx_x_xx_0 0_xx_x_xxz_0 0_xx_x_xxy_0 0_xx_x_xxx_0 0_xx_xz_zz_0 0_xx_xz_yz_0 0_xx_xz_yy_0 0_xx_xz_xz_0 0_xx_xz_xy_0 0_xx_xz_xx_0 0_xx_xy_zz_0 0_xx_xy_yz_0 0_xx_xy_yy_0 0_xx_xy_xz_0 0_xx_xy_xy_0 0_xx_xy_xx_0 0_xx_xx_zz_0 0_xx_xx_yz_0 0_xx_xx_yy_0 0_xx_xx_xz_0 0_xx_xx_xy_0 0_xx_xx_xx_0 

signature:
0_zz_0_zzz_0 0_zz_0_zzzz_0 0_zz_0_yzz_0 0_zz_0_yzzz_0 0_zz_0_yyz_0 0_zz_0_yyzz_0 0_zz_0_yyy_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzz_0 0_zz_0_xzzz_0 0_zz_0_xyz_0 0_zz_0_xyzz_0 0_zz_0_xyy_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxz_0 0_zz_0_xxzz_0 0_zz_0_xxy_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxx_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_zz_z_zzz_0 0_zz_z_yzz_0 0_zz_z_yyz_0 0_zz_z_xzz_0 0_zz_z_xyz_0 0_zz_z_xxz_0 0_zz_y_zzz_0 0_zz_y_yzz_0 0_zz_y_yyz_0 0_zz_y_yyy_0 0_zz_y_xzz_0 0_zz_y_xyz_0 0_zz_y_xyy_0 0_zz_y_xxz_0 0_zz_y_xxy_0 0_zz_x_zzz_0 0_zz_x_yzz_0 0_zz_x_yyz_0 0_zz_x_yyy_0 0_zz_x_xzz_0 0_zz_x_xyz_0 0_zz_x_xyy_0 0_zz_x_xxz_0 0_zz_x_xxy_0 0_zz_x_xxx_0 0_yz_0_zzz_0 0_yz_0_zzzz_0 0_yz_0_yzz_0 0_yz_0_yzzz_0 0_yz_0_yyz_0 0_yz_0_yyzz_0 0_yz_0_yyy_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzz_0 0_yz_0_xzzz_0 0_yz_0_xyz_0 0_yz_0_xyzz_0 0_yz_0_xyy_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxz_0 0_yz_0_xxzz_0 0_yz_0_xxy_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxx_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yz_z_zzz_0 0_yz_z_yzz_0 0_yz_z_yyz_0 0_yz_z_xzz_0 0_yz_z_xyz_0 0_yz_z_xxz_0 0_yz_y_zzz_0 0_yz_y_yzz_0 0_yz_y_yyz_0 0_yz_y_yyy_0 0_yz_y_xzz_0 0_yz_y_xyz_0 0_yz_y_xyy_0 0_yz_y_xxz_0 0_yz_y_xxy_0 0_yz_x_zzz_0 0_yz_x_yzz_0 0_yz_x_yyz_0 0_yz_x_yyy_0 0_yz_x_xzz_0 0_yz_x_xyz_0 0_yz_x_xyy_0 0_yz_x_xxz_0 0_yz_x_xxy_0 0_yz_x_xxx_0 0_yy_0_zzz_0 0_yy_0_zzzz_0 0_yy_0_yzz_0 0_yy_0_yzzz_0 0_yy_0_yyz_0 0_yy_0_yyzz_0 0_yy_0_yyy_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzz_0 0_yy_0_xzzz_0 0_yy_0_xyz_0 0_yy_0_xyzz_0 0_yy_0_xyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxz_0 0_yy_0_xxzz_0 0_yy_0_xxy_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxx_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_yy_z_zzz_0 0_yy_z_yzz_0 0_yy_z_yyz_0 0_yy_z_xzz_0 0_yy_z_xyz_0 0_yy_z_xxz_0 0_yy_y_zzz_0 0_yy_y_yzz_0 0_yy_y_yyz_0 0_yy_y_yyy_0 0_yy_y_xzz_0 0_yy_y_xyz_0 0_yy_y_xyy_0 0_yy_y_xxz_0 0_yy_y_xxy_0 0_yy_x_zzz_0 0_yy_x_yzz_0 0_yy_x_yyz_0 0_yy_x_yyy_0 0_yy_x_xzz_0 0_yy_x_xyz_0 0_yy_x_xyy_0 0_yy_x_xxz_0 0_yy_x_xxy_0 0_yy_x_xxx_0 0_xz_0_zzz_0 0_xz_0_zzzz_0 0_xz_0_yzz_0 0_xz_0_yzzz_0 0_xz_0_yyz_0 0_xz_0_yyzz_0 0_xz_0_yyy_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzz_0 0_xz_0_xzzz_0 0_xz_0_xyz_0 0_xz_0_xyzz_0 0_xz_0_xyy_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxz_0 0_xz_0_xxzz_0 0_xz_0_xxy_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxx_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xz_z_zzz_0 0_xz_z_yzz_0 0_xz_z_yyz_0 0_xz_z_xzz_0 0_xz_z_xyz_0 0_xz_z_xxz_0 0_xz_y_zzz_0 0_xz_y_yzz_0 0_xz_y_yyz_0 0_xz_y_yyy_0 0_xz_y_xzz_0 0_xz_y_xyz_0 0_xz_y_xyy_0 0_xz_y_xxz_0 0_xz_y_xxy_0 0_xz_x_zzz_0 0_xz_x_yzz_0 0_xz_x_yyz_0 0_xz_x_yyy_0 0_xz_x_xzz_0 0_xz_x_xyz_0 0_xz_x_xyy_0 0_xz_x_xxz_0 0_xz_x_xxy_0 0_xz_x_xxx_0 0_xy_0_zzz_0 0_xy_0_zzzz_0 0_xy_0_yzz_0 0_xy_0_yzzz_0 0_xy_0_yyz_0 0_xy_0_yyzz_0 0_xy_0_yyy_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzz_0 0_xy_0_xzzz_0 0_xy_0_xyz_0 0_xy_0_xyzz_0 0_xy_0_xyy_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxz_0 0_xy_0_xxzz_0 0_xy_0_xxy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxx_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xy_z_zzz_0 0_xy_z_yzz_0 0_xy_z_yyz_0 0_xy_z_xzz_0 0_xy_z_xyz_0 0_xy_z_xxz_0 0_xy_y_zzz_0 0_xy_y_yzz_0 0_xy_y_yyz_0 0_xy_y_yyy_0 0_xy_y_xzz_0 0_xy_y_xyz_0 0_xy_y_xyy_0 0_xy_y_xxz_0 0_xy_y_xxy_0 0_xy_x_zzz_0 0_xy_x_yzz_0 0_xy_x_yyz_0 0_xy_x_yyy_0 0_xy_x_xzz_0 0_xy_x_xyz_0 0_xy_x_xyy_0 0_xy_x_xxz_0 0_xy_x_xxy_0 0_xy_x_xxx_0 0_xx_0_zzz_0 0_xx_0_zzzz_0 0_xx_0_yzz_0 0_xx_0_yzzz_0 0_xx_0_yyz_0 0_xx_0_yyzz_0 0_xx_0_yyy_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzz_0 0_xx_0_xzzz_0 0_xx_0_xyz_0 0_xx_0_xyzz_0 0_xx_0_xyy_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxz_0 0_xx_0_xxzz_0 0_xx_0_xxy_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxx_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 0_xx_z_zzz_0 0_xx_z_yzz_0 0_xx_z_yyz_0 0_xx_z_xzz_0 0_xx_z_xyz_0 0_xx_z_xxz_0 0_xx_y_zzz_0 0_xx_y_yzz_0 0_xx_y_yyz_0 0_xx_y_yyy_0 0_xx_y_xzz_0 0_xx_y_xyz_0 0_xx_y_xyy_0 0_xx_y_xxz_0 0_xx_y_xxy_0 0_xx_x_zzz_0 0_xx_x_yzz_0 0_xx_x_yyz_0 0_xx_x_yyy_0 0_xx_x_xzz_0 0_xx_x_xyz_0 0_xx_x_xyy_0 0_xx_x_xxz_0 0_xx_x_xxy_0 0_xx_x_xxx_0 

signature:
0_zz_0_zz_0 0_zz_0_zzz_0 0_zz_0_yz_0 0_zz_0_yzz_0 0_zz_0_yy_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xz_0 0_zz_0_xzz_0 0_zz_0_xy_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xx_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_zz_z_zz_0 0_zz_z_yz_0 0_zz_z_yy_0 0_zz_z_xz_0 0_zz_z_xy_0 0_zz_z_xx_0 0_zz_y_zz_0 0_zz_y_yz_0 0_zz_y_yy_0 0_zz_y_xz_0 0_zz_y_xy_0 0_zz_y_xx_0 0_zz_x_zz_0 0_zz_x_yz_0 0_zz_x_yy_0 0_zz_x_xz_0 0_zz_x_xy_0 0_zz_x_xx_0 0_yz_0_zz_0 0_yz_0_zzz_0 0_yz_0_yz_0 0_yz_0_yzz_0 0_yz_0_yy_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xz_0 0_yz_0_xzz_0 0_yz_0_xy_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xx_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yz_z_zz_0 0_yz_z_yz_0 0_yz_z_yy_0 0_yz_z_xz_0 0_yz_z_xy_0 0_yz_z_xx_0 0_yz_y_zz_0 0_yz_y_yz_0 0_yz_y_yy_0 0_yz_y_xz_0 0_yz_y_xy_0 0_yz_y_xx_0 0_yz_x_zz_0 0_yz_x_yz_0 0_yz_x_yy_0 0_yz_x_xz_0 0_yz_x_xy_0 0_yz_x_xx_0 0_yy_0_zz_0 0_yy_0_zzz_0 0_yy_0_yz_0 0_yy_0_yzz_0 0_yy_0_yy_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xz_0 0_yy_0_xzz_0 0_yy_0_xy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xx_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_yy_z_zz_0 0_yy_z_yz_0 0_yy_z_yy_0 0_yy_z_xz_0 0_yy_z_xy_0 0_yy_z_xx_0 0_yy_y_zz_0 0_yy_y_yz_0 0_yy_y_yy_0 0_yy_y_xz_0 0_yy_y_xy_0 0_yy_y_xx_0 0_yy_x_zz_0 0_yy_x_yz_0 0_yy_x_yy_0 0_yy_x_xz_0 0_yy_x_xy_0 0_yy_x_xx_0 0_xz_0_zz_0 0_xz_0_zzz_0 0_xz_0_yz_0 0_xz_0_yzz_0 0_xz_0_yy_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xz_0 0_xz_0_xzz_0 0_xz_0_xy_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xx_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xz_z_zz_0 0_xz_z_yz_0 0_xz_z_yy_0 0_xz_z_xz_0 0_xz_z_xy_0 0_xz_z_xx_0 0_xz_y_zz_0 0_xz_y_yz_0 0_xz_y_yy_0 0_xz_y_xz_0 0_xz_y_xy_0 0_xz_y_xx_0 0_xz_x_zz_0 0_xz_x_yz_0 0_xz_x_yy_0 0_xz_x_xz_0 0_xz_x_xy_0 0_xz_x_xx_0 0_xy_0_zz_0 0_xy_0_zzz_0 0_xy_0_yz_0 0_xy_0_yzz_0 0_xy_0_yy_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xz_0 0_xy_0_xzz_0 0_xy_0_xy_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xx_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xy_z_zz_0 0_xy_z_yz_0 0_xy_z_yy_0 0_xy_z_xz_0 0_xy_z_xy_0 0_xy_z_xx_0 0_xy_y_zz_0 0_xy_y_yz_0 0_xy_y_yy_0 0_xy_y_xz_0 0_xy_y_xy_0 0_xy_y_xx_0 0_xy_x_zz_0 0_xy_x_yz_0 0_xy_x_yy_0 0_xy_x_xz_0 0_xy_x_xy_0 0_xy_x_xx_0 0_xx_0_zz_0 0_xx_0_zzz_0 0_xx_0_yz_0 0_xx_0_yzz_0 0_xx_0_yy_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xz_0 0_xx_0_xzz_0 0_xx_0_xy_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xx_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 0_xx_z_zz_0 0_xx_z_yz_0 0_xx_z_yy_0 0_xx_z_xz_0 0_xx_z_xy_0 0_xx_z_xx_0 0_xx_y_zz_0 0_xx_y_yz_0 0_xx_y_yy_0 0_xx_y_xz_0 0_xx_y_xy_0 0_xx_y_xx_0 0_xx_x_zz_0 0_xx_x_yz_0 0_xx_x_yy_0 0_xx_x_xz_0 0_xx_x_xy_0 0_xx_x_xx_0 

signature:
0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyy_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_yyyy_0 0_z_0_yyyy_1 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyy_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xyyy_0 0_z_0_xyyy_1 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxy_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxyy_0 0_z_0_xxyy_1 0_z_0_xxx_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_z_0_xxxy_0 0_z_0_xxxy_1 0_z_0_xxxx_0 0_z_0_xxxx_1 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_y_0_zzz_1 0_y_0_zzzz_0 0_y_0_zzzz_1 0_y_0_yzz_1 0_y_0_yzzz_0 0_y_0_yzzz_1 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xzz_1 0_y_0_xzzz_0 0_y_0_xzzz_1 0_y_0_xyz_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxz_1 0_y_0_xxzz_0 0_y_0_xxzz_1 0_y_0_xxy_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxx_1 0_y_0_xxxz_0 0_y_0_xxxz_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_y_0_xxxx_0 0_y_0_xxxx_1 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzzz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_x_0_zzz_1 0_x_0_zzzz_0 0_x_0_zzzz_1 0_x_0_yzz_1 0_x_0_yzzz_0 0_x_0_yzzz_1 0_x_0_yyz_1 0_x_0_yyzz_0 0_x_0_yyzz_1 0_x_0_yyy_1 0_x_0_yyyz_0 0_x_0_yyyz_1 0_x_0_yyyy_0 0_x_0_yyyy_1 0_x_0_xzz_1 0_x_0_xzzz_0 0_x_0_xzzz_1 0_x_0_xyz_1 0_x_0_xyzz_0 0_x_0_xyzz_1 0_x_0_xyy_1 0_x_0_xyyz_0 0_x_0_xyyz_1 0_x_0_xyyy_0 0_x_0_xyyy_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxzz_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

signature:
0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_yyy_0 0_z_0_yyy_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xyy_0 0_z_0_xyy_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_z_0_xxy_0 0_z_0_xxy_1 0_z_0_xxx_0 0_z_0_xxx_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_y_0_zz_1 0_y_0_zzz_0 0_y_0_zzz_1 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xzz_0 0_y_0_xzz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxz_0 0_y_0_xxz_1 0_y_0_xxy_0 0_y_0_xxy_1 0_y_0_xxx_0 0_y_0_xxx_1 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_x_0_zz_1 0_x_0_zzz_0 0_x_0_zzz_1 0_x_0_yz_1 0_x_0_yzz_0 0_x_0_yzz_1 0_x_0_yy_1 0_x_0_yyz_0 0_x_0_yyz_1 0_x_0_yyy_0 0_x_0_yyy_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

signature:
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

signature:
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

signature:
0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_y_0 0_0_0_y_1 0_0_0_x_0 0_0_0_x_1 0_z_0_z_0 0_z_0_y_0 0_z_0_x_0 0_y_0_z_0 0_y_0_y_0 0_y_0_x_0 0_x_0_z_0 0_x_0_y_0 0_x_0_x_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xz_0 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xy_0 0_0_0_xyz_0 0_0_0_xyzz_0 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxy_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxx_0 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxx_0 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxx_0 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xzz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xx_0 

signature:
0_0_0_0_0 0_0_0_0_1 0_0_0_z_0 0_0_0_z_1 0_0_0_zz_0 0_0_0_y_0 0_0_0_y_1 0_0_0_yz_0 0_0_0_yy_0 0_0_0_x_0 0_0_0_x_1 0_0_0_xx_0 

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

0_0_0_0_3 0_0_0_0_4 0_0_0_z_3 0_0_0_z_4 0_0_0_zz_3 0_0_0_y_3 0_0_0_y_4 0_0_0_yz_3 0_0_0_yy_3 0_0_0_x_3 0_0_0_x_4 0_0_0_xx_3 

SSSS_3 SSSS_4 SSSP_3 SSSP_4 SSSD_3 

SSSS_3 : 0_0_0_0_3 

SSSS_4 : 0_0_0_0_4 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xx_3 

INTEGRAL:0 : 2 : 4 : SSSS_4 Y SSSS_5 Y SSSP_4 Y SSSP_5 Y SSSD_4 Y 

0_0_0_0_4 0_0_0_0_5 0_0_0_z_4 0_0_0_z_5 0_0_0_zz_4 0_0_0_y_4 0_0_0_y_5 0_0_0_yy_4 0_0_0_x_4 0_0_0_x_5 0_0_0_xx_4 

SSSS_4 SSSS_5 SSSP_4 SSSP_5 SSSD_4 

SSSS_4 : 0_0_0_0_4 

SSSS_5 : 0_0_0_0_5 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSD_4 : 0_0_0_zz_4 0_0_0_yy_4 0_0_0_xx_4 

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

INTEGRAL:0 : 3 : 2 : SSSP_2 Y SSSP_3 Y SSSD_2 N SSSD_3 Y SSSF_2 Y 

0_0_0_z_2 0_0_0_z_3 0_0_0_zz_2 0_0_0_zz_3 0_0_0_zzz_2 0_0_0_y_2 0_0_0_y_3 0_0_0_yz_2 0_0_0_yz_3 0_0_0_yzz_2 0_0_0_yy_2 0_0_0_yy_3 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_x_2 0_0_0_x_3 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xx_2 0_0_0_xx_3 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSP_2 SSSP_3 SSSD_2 SSSD_3 SSSF_2 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xx_3 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

INTEGRAL:0 : 3 : 3 : SSSP_3 Y SSSP_4 Y SSSD_3 N SSSD_4 Y SSSF_3 Y 

0_0_0_z_3 0_0_0_z_4 0_0_0_zz_3 0_0_0_zz_4 0_0_0_zzz_3 0_0_0_y_3 0_0_0_y_4 0_0_0_yzz_3 0_0_0_yy_3 0_0_0_yy_4 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_x_3 0_0_0_x_4 0_0_0_xzz_3 0_0_0_xyy_3 0_0_0_xx_3 0_0_0_xx_4 0_0_0_xxz_3 0_0_0_xxx_3 

SSSP_3 SSSP_4 SSSD_3 SSSD_4 SSSF_3 

SSSP_3 : 0_0_0_z_3 0_0_0_y_3 0_0_0_x_3 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xx_3 

SSSD_4 : 0_0_0_zz_4 0_0_0_yy_4 0_0_0_xx_4 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxx_3 

INTEGRAL:0 : 4 : 0 : SSSD_0 N SSSD_1 N SSSF_0 N SSSF_1 N SSSG_0 Y 

0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxx_0 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

SSSD_0 SSSD_1 SSSF_0 SSSF_1 SSSG_0 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSG_0 : 0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

INTEGRAL:0 : 4 : 1 : SSSD_1 N SSSD_2 N SSSF_1 N SSSF_2 N SSSG_1 Y 

0_0_0_zz_1 0_0_0_zz_2 0_0_0_zzz_1 0_0_0_zzz_2 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzz_2 0_0_0_yzzz_1 0_0_0_yy_1 0_0_0_yy_2 0_0_0_yyz_1 0_0_0_yyz_2 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyy_2 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzz_2 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyy_2 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xx_1 0_0_0_xx_2 0_0_0_xxz_1 0_0_0_xxz_2 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxx_2 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SSSD_1 SSSD_2 SSSF_1 SSSF_2 SSSG_1 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

INTEGRAL:0 : 4 : 2 : SSSD_2 N SSSD_3 N SSSF_2 N SSSF_3 Y SSSG_2 Y 

0_0_0_zz_2 0_0_0_zz_3 0_0_0_zzz_2 0_0_0_zzz_3 0_0_0_zzzz_2 0_0_0_yzz_2 0_0_0_yzz_3 0_0_0_yzzz_2 0_0_0_yy_2 0_0_0_yy_3 0_0_0_yyz_2 0_0_0_yyz_3 0_0_0_yyzz_2 0_0_0_yyy_2 0_0_0_yyy_3 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzz_2 0_0_0_xzz_3 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyy_2 0_0_0_xyy_3 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xx_2 0_0_0_xx_3 0_0_0_xxz_2 0_0_0_xxz_3 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxx_2 0_0_0_xxx_3 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SSSD_2 SSSD_3 SSSF_2 SSSF_3 SSSG_2 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xx_3 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxx_3 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

INTEGRAL:1 : 1 : 1 : SSSS_2 Y SSSP_1 Y SSSP_2 Y SPSP_1 Y 

0_0_0_0_2 0_0_0_z_1 0_0_0_z_2 0_0_0_y_1 0_0_0_y_2 0_0_0_x_1 0_0_0_x_2 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SSSS_2 SSSP_1 SSSP_2 SPSP_1 

SSSS_2 : 0_0_0_0_2 

SSSP_1 : 0_0_0_z_1 0_0_0_y_1 0_0_0_x_1 

SSSP_2 : 0_0_0_z_2 0_0_0_y_2 0_0_0_x_2 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

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

INTEGRAL:1 : 4 : 0 : SSSF_1 Y SSSG_0 Y SSSG_1 Y SPSG_0 Y 

0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

SSSF_1 SSSG_0 SSSG_1 SPSG_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SSSG_0 : 0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SPSG_0 : 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

INTEGRAL:1 : 4 : 1 : SSSF_2 Y SSSG_1 Y SSSG_2 Y SPSG_1 Y 

0_0_0_zzz_2 0_0_0_zzzz_1 0_0_0_zzzz_2 0_0_0_yzz_2 0_0_0_yzzz_1 0_0_0_yzzz_2 0_0_0_yyz_2 0_0_0_yyzz_1 0_0_0_yyzz_2 0_0_0_yyy_2 0_0_0_yyyz_1 0_0_0_yyyz_2 0_0_0_yyyy_1 0_0_0_yyyy_2 0_0_0_xzz_2 0_0_0_xzzz_1 0_0_0_xzzz_2 0_0_0_xyz_2 0_0_0_xyzz_1 0_0_0_xyzz_2 0_0_0_xyy_2 0_0_0_xyyz_1 0_0_0_xyyz_2 0_0_0_xyyy_1 0_0_0_xyyy_2 0_0_0_xxz_2 0_0_0_xxzz_1 0_0_0_xxzz_2 0_0_0_xxy_2 0_0_0_xxyz_1 0_0_0_xxyz_2 0_0_0_xxyy_1 0_0_0_xxyy_2 0_0_0_xxx_2 0_0_0_xxxz_1 0_0_0_xxxz_2 0_0_0_xxxy_1 0_0_0_xxxy_2 0_0_0_xxxx_1 0_0_0_xxxx_2 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_yyyy_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xyyy_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxyy_1 0_z_0_xxxz_1 0_z_0_xxxy_1 0_z_0_xxxx_1 0_y_0_zzzz_1 0_y_0_yzzz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xzzz_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxzz_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxz_1 0_y_0_xxxy_1 0_y_0_xxxx_1 0_x_0_zzzz_1 0_x_0_yzzz_1 0_x_0_yyzz_1 0_x_0_yyyz_1 0_x_0_yyyy_1 0_x_0_xzzz_1 0_x_0_xyzz_1 0_x_0_xyyz_1 0_x_0_xyyy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SSSF_2 SSSG_1 SSSG_2 SPSG_1 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_yyyy_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xyyy_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxyy_1 0_z_0_xxxz_1 0_z_0_xxxy_1 0_z_0_xxxx_1 0_y_0_zzzz_1 0_y_0_yzzz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xzzz_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxzz_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxz_1 0_y_0_xxxy_1 0_y_0_xxxx_1 0_x_0_zzzz_1 0_x_0_yzzz_1 0_x_0_yyzz_1 0_x_0_yyyz_1 0_x_0_yyyy_1 0_x_0_xzzz_1 0_x_0_xyzz_1 0_x_0_xyyz_1 0_x_0_xyyy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

INTEGRAL:2 : 2 : 0 : SSSD_0 Y SSSD_1 Y SPSP_1 Y SPSD_0 Y SPSD_1 Y SDSD_0 Y 

0_0_0_zz_0 0_0_0_zz_1 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yy_0 0_0_0_yy_1 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xx_0 0_0_0_xx_1 0_z_0_z_1 0_z_0_zz_0 0_z_0_zz_1 0_z_0_y_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_x_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_y_0_z_1 0_y_0_zz_0 0_y_0_zz_1 0_y_0_y_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_x_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_x_0_z_1 0_x_0_zz_0 0_x_0_zz_1 0_x_0_y_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_x_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

SSSD_0 SSSD_1 SPSP_1 SPSD_0 SPSD_1 SDSD_0 

SSSD_0 : 0_0_0_zz_0 0_0_0_yz_0 0_0_0_yy_0 0_0_0_xz_0 0_0_0_xy_0 0_0_0_xx_0 

SSSD_1 : 0_0_0_zz_1 0_0_0_yz_1 0_0_0_yy_1 0_0_0_xz_1 0_0_0_xy_1 0_0_0_xx_1 

SPSP_1 : 0_z_0_z_1 0_z_0_y_1 0_z_0_x_1 0_y_0_z_1 0_y_0_y_1 0_y_0_x_1 0_x_0_z_1 0_x_0_y_1 0_x_0_x_1 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SDSD_0 : 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

INTEGRAL:2 : 3 : 0 : SSSF_0 Y SSSF_1 Y SPSD_1 Y SPSF_0 Y SPSF_1 Y SDSF_0 Y 

0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxx_0 0_0_0_xxx_1 0_z_0_zz_1 0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yy_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_yyy_0 0_z_0_yyy_1 0_z_0_xz_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xy_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xyy_0 0_z_0_xyy_1 0_z_0_xx_1 0_z_0_xxz_0 0_z_0_xxz_1 0_z_0_xxy_0 0_z_0_xxy_1 0_z_0_xxx_0 0_z_0_xxx_1 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_y_0_zz_1 0_y_0_zzz_0 0_y_0_zzz_1 0_y_0_yz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yy_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xz_1 0_y_0_xzz_0 0_y_0_xzz_1 0_y_0_xy_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xx_1 0_y_0_xxz_0 0_y_0_xxz_1 0_y_0_xxy_0 0_y_0_xxy_1 0_y_0_xxx_0 0_y_0_xxx_1 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_x_0_zz_1 0_x_0_zzz_0 0_x_0_zzz_1 0_x_0_yz_1 0_x_0_yzz_0 0_x_0_yzz_1 0_x_0_yy_1 0_x_0_yyz_0 0_x_0_yyz_1 0_x_0_yyy_0 0_x_0_yyy_1 0_x_0_xz_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xy_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xx_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

SSSF_0 SSSF_1 SPSD_1 SPSF_0 SPSF_1 SDSF_0 

SSSF_0 : 0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 

SSSF_1 : 0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SPSF_0 : 0_z_0_zzz_0 0_z_0_yzz_0 0_z_0_yyz_0 0_z_0_yyy_0 0_z_0_xzz_0 0_z_0_xyz_0 0_z_0_xyy_0 0_z_0_xxz_0 0_z_0_xxy_0 0_z_0_xxx_0 0_y_0_zzz_0 0_y_0_yzz_0 0_y_0_yyz_0 0_y_0_yyy_0 0_y_0_xzz_0 0_y_0_xyz_0 0_y_0_xyy_0 0_y_0_xxz_0 0_y_0_xxy_0 0_y_0_xxx_0 0_x_0_zzz_0 0_x_0_yzz_0 0_x_0_yyz_0 0_x_0_yyy_0 0_x_0_xzz_0 0_x_0_xyz_0 0_x_0_xyy_0 0_x_0_xxz_0 0_x_0_xxy_0 0_x_0_xxx_0 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SDSF_0 : 0_zz_0_zzz_0 0_zz_0_yzz_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xzz_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_yz_0_zzz_0 0_yz_0_yzz_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yy_0_zzz_0 0_yy_0_yzz_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xzz_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xzz_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xx_0_zzz_0 0_xx_0_yzz_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xzz_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 

INTEGRAL:2 : 4 : 0 : SSSG_0 Y SSSG_1 Y SPSF_1 Y SPSG_0 Y SPSG_1 Y SDSG_0 Y 

0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyy_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_yyyy_0 0_z_0_yyyy_1 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyy_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xyyy_0 0_z_0_xyyy_1 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxy_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxyy_0 0_z_0_xxyy_1 0_z_0_xxx_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_z_0_xxxy_0 0_z_0_xxxy_1 0_z_0_xxxx_0 0_z_0_xxxx_1 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_y_0_zzz_1 0_y_0_zzzz_0 0_y_0_zzzz_1 0_y_0_yzz_1 0_y_0_yzzz_0 0_y_0_yzzz_1 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xzz_1 0_y_0_xzzz_0 0_y_0_xzzz_1 0_y_0_xyz_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxz_1 0_y_0_xxzz_0 0_y_0_xxzz_1 0_y_0_xxy_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxx_1 0_y_0_xxxz_0 0_y_0_xxxz_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_y_0_xxxx_0 0_y_0_xxxx_1 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzzz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_x_0_zzz_1 0_x_0_zzzz_0 0_x_0_zzzz_1 0_x_0_yzz_1 0_x_0_yzzz_0 0_x_0_yzzz_1 0_x_0_yyz_1 0_x_0_yyzz_0 0_x_0_yyzz_1 0_x_0_yyy_1 0_x_0_yyyz_0 0_x_0_yyyz_1 0_x_0_yyyy_0 0_x_0_yyyy_1 0_x_0_xzz_1 0_x_0_xzzz_0 0_x_0_xzzz_1 0_x_0_xyz_1 0_x_0_xyzz_0 0_x_0_xyzz_1 0_x_0_xyy_1 0_x_0_xyyz_0 0_x_0_xyyz_1 0_x_0_xyyy_0 0_x_0_xyyy_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxzz_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

SSSG_0 SSSG_1 SPSF_1 SPSG_0 SPSG_1 SDSG_0 

SSSG_0 : 0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SPSG_0 : 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_yyyy_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xyyy_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxyy_1 0_z_0_xxxz_1 0_z_0_xxxy_1 0_z_0_xxxx_1 0_y_0_zzzz_1 0_y_0_yzzz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xzzz_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxzz_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxz_1 0_y_0_xxxy_1 0_y_0_xxxx_1 0_x_0_zzzz_1 0_x_0_yzzz_1 0_x_0_yyzz_1 0_x_0_yyyz_1 0_x_0_yyyy_1 0_x_0_xzzz_1 0_x_0_xyzz_1 0_x_0_xyyz_1 0_x_0_xyyy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SDSG_0 : 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzzz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxzz_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

        }
    }
}


} // derirec namespace
