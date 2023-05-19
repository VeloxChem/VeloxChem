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
#include "EriVRRForSPSS.hpp"
#include "EriVRRForSPSP.hpp"
#include "EriVRRForSPSD.hpp"
#include "EriVRRForSPSF.hpp"
#include "EriVRRForSPSG.hpp"
#include "EriVRRForSDSP.hpp"
#include "EriVRRForSDSD.hpp"
#include "EriVRRForSDSF.hpp"
#include "EriVRRForSDSG.hpp"
#include "EriVRRForSFSD.hpp"
#include "EriVRRForSFSF.hpp"
#include "EriVRRForSFSG.hpp"
#include "EriHRRForSDPD.hpp"
#include "EriHRRForSDPF.hpp"
#include "EriHRRForSDDD.hpp"
#include "EriHRRForSFPD.hpp"
#include "EriHRRForSFPF.hpp"
#include "EriHRRForSFDD.hpp"
#include "EriHRRForPDDD.hpp"

namespace derirec { // derirec namespace

template <typename T>
auto
compHostPDDD(      T*                                 intsBuffer,
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

    BufferHostXY<T> bvals(8, ncpairs);

    CBoysFunc<T, 7> bftable;

    // Primitive integral buffers

    BufferHostXY<T> pbufSSSS(8, ncpairs);

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
    BufferHostXY<T> pbufSSSP6(3, ncpairs);

0_0_0_z_6 0_0_0_y_6 0_0_0_x_6 
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
    BufferHostXY<T> pbufSSSD5(3, ncpairs);

0_0_0_zz_5 0_0_0_yy_5 0_0_0_xx_5 
    BufferHostXY<T> pbufSSSF0(10, ncpairs);

0_0_0_zzz_0 0_0_0_yzz_0 0_0_0_yyz_0 0_0_0_yyy_0 0_0_0_xzz_0 0_0_0_xyz_0 0_0_0_xyy_0 0_0_0_xxz_0 0_0_0_xxy_0 0_0_0_xxx_0 
    BufferHostXY<T> pbufSSSF1(10, ncpairs);

0_0_0_zzz_1 0_0_0_yzz_1 0_0_0_yyz_1 0_0_0_yyy_1 0_0_0_xzz_1 0_0_0_xyz_1 0_0_0_xyy_1 0_0_0_xxz_1 0_0_0_xxy_1 0_0_0_xxx_1 
    BufferHostXY<T> pbufSSSF2(10, ncpairs);

0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 
    BufferHostXY<T> pbufSSSF3(10, ncpairs);

0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 
    BufferHostXY<T> pbufSSSF4(8, ncpairs);

0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxx_4 
    BufferHostXY<T> pbufSSSG0(15, ncpairs);

0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 
    BufferHostXY<T> pbufSSSG1(15, ncpairs);

0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 
    BufferHostXY<T> pbufSSSG2(15, ncpairs);

0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 
    BufferHostXY<T> pbufSSSG3(15, ncpairs);

0_0_0_zzzz_3 0_0_0_yzzz_3 0_0_0_yyzz_3 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 
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
    BufferHostXY<T> pbufSPSG0(45, ncpairs);

0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 
    BufferHostXY<T> pbufSPSG1(45, ncpairs);

0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_yyyy_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xyyy_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxyy_1 0_z_0_xxxz_1 0_z_0_xxxy_1 0_z_0_xxxx_1 0_y_0_zzzz_1 0_y_0_yzzz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xzzz_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxzz_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxz_1 0_y_0_xxxy_1 0_y_0_xxxx_1 0_x_0_zzzz_1 0_x_0_yzzz_1 0_x_0_yyzz_1 0_x_0_yyyz_1 0_x_0_yyyy_1 0_x_0_xzzz_1 0_x_0_xyzz_1 0_x_0_xyyz_1 0_x_0_xyyy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 
    BufferHostXY<T> pbufSPSG2(45, ncpairs);

0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_yyyz_2 0_z_0_yyyy_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xyyy_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_z_0_xxyy_2 0_z_0_xxxz_2 0_z_0_xxxy_2 0_z_0_xxxx_2 0_y_0_zzzz_2 0_y_0_yzzz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xzzz_2 0_y_0_xyzz_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxzz_2 0_y_0_xxyz_2 0_y_0_xxyy_2 0_y_0_xxxz_2 0_y_0_xxxy_2 0_y_0_xxxx_2 0_x_0_zzzz_2 0_x_0_yzzz_2 0_x_0_yyzz_2 0_x_0_yyyz_2 0_x_0_yyyy_2 0_x_0_xzzz_2 0_x_0_xyzz_2 0_x_0_xyyz_2 0_x_0_xyyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 
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
    BufferHostXY<T> pbufSDSG0(90, ncpairs);

0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzzz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxzz_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 
    BufferHostXY<T> pbufSDSG1(60, ncpairs);

0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_yyyz_1 0_zz_0_yyyy_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xyyz_1 0_zz_0_xyyy_1 0_zz_0_xxzz_1 0_zz_0_xxyz_1 0_zz_0_xxyy_1 0_zz_0_xxxz_1 0_zz_0_xxxy_1 0_zz_0_xxxx_1 0_yz_0_zzzz_1 0_yz_0_yzzz_1 0_yz_0_yyzz_1 0_yz_0_yyyz_1 0_yz_0_yyyy_1 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_zzzz_1 0_yy_0_yzzz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xzzz_1 0_yy_0_xyzz_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxzz_1 0_yy_0_xxyz_1 0_yy_0_xxyy_1 0_yy_0_xxxz_1 0_yy_0_xxxy_1 0_yy_0_xxxx_1 0_xz_0_xzzz_1 0_xz_0_xxzz_1 0_xz_0_xxxz_1 0_xz_0_xxxx_1 0_xy_0_xyyy_1 0_xy_0_xxyy_1 0_xy_0_xxxy_1 0_xx_0_zzzz_1 0_xx_0_yzzz_1 0_xx_0_yyzz_1 0_xx_0_yyyz_1 0_xx_0_yyyy_1 0_xx_0_xzzz_1 0_xx_0_xyzz_1 0_xx_0_xyyz_1 0_xx_0_xyyy_1 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 
    BufferHostXY<T> pbufSFSD0(60, ncpairs);

0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 
    BufferHostXY<T> pbufSFSF0(100, ncpairs);

0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_yyz_0 0_zzz_0_yyy_0 0_zzz_0_xzz_0 0_zzz_0_xyz_0 0_zzz_0_xyy_0 0_zzz_0_xxz_0 0_zzz_0_xxy_0 0_zzz_0_xxx_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_yyy_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yzz_0_xyy_0 0_yzz_0_xxz_0 0_yzz_0_xxy_0 0_yzz_0_xxx_0 0_yyz_0_zzz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xzz_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xxz_0 0_yyz_0_xxy_0 0_yyz_0_xxx_0 0_yyy_0_zzz_0 0_yyy_0_yzz_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xzz_0 0_yyy_0_xyz_0 0_yyy_0_xyy_0 0_yyy_0_xxz_0 0_yyy_0_xxy_0 0_yyy_0_xxx_0 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_yyz_0 0_xzz_0_yyy_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xyy_0 0_xzz_0_xxz_0 0_xzz_0_xxy_0 0_xzz_0_xxx_0 0_xyz_0_zzz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_yyy_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxx_0 0_xyy_0_zzz_0 0_xyy_0_yzz_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xzz_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxz_0 0_xyy_0_xxy_0 0_xyy_0_xxx_0 0_xxz_0_zzz_0 0_xxz_0_yzz_0 0_xxz_0_yyz_0 0_xxz_0_yyy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xyy_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_zzz_0 0_xxy_0_yzz_0 0_xxy_0_yyz_0 0_xxy_0_yyy_0 0_xxy_0_xzz_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_zzz_0 0_xxx_0_yzz_0 0_xxx_0_yyz_0 0_xxx_0_yyy_0 0_xxx_0_xzz_0 0_xxx_0_xyz_0 0_xxx_0_xyy_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 
    BufferHostXY<T> pbufSFSG0(150, ncpairs);

0_zzz_0_zzzz_0 0_zzz_0_yzzz_0 0_zzz_0_yyzz_0 0_zzz_0_yyyz_0 0_zzz_0_yyyy_0 0_zzz_0_xzzz_0 0_zzz_0_xyzz_0 0_zzz_0_xyyz_0 0_zzz_0_xyyy_0 0_zzz_0_xxzz_0 0_zzz_0_xxyz_0 0_zzz_0_xxyy_0 0_zzz_0_xxxz_0 0_zzz_0_xxxy_0 0_zzz_0_xxxx_0 0_yzz_0_zzzz_0 0_yzz_0_yzzz_0 0_yzz_0_yyzz_0 0_yzz_0_yyyz_0 0_yzz_0_yyyy_0 0_yzz_0_xzzz_0 0_yzz_0_xyzz_0 0_yzz_0_xyyz_0 0_yzz_0_xyyy_0 0_yzz_0_xxzz_0 0_yzz_0_xxyz_0 0_yzz_0_xxyy_0 0_yzz_0_xxxz_0 0_yzz_0_xxxy_0 0_yzz_0_xxxx_0 0_yyz_0_zzzz_0 0_yyz_0_yzzz_0 0_yyz_0_yyzz_0 0_yyz_0_yyyz_0 0_yyz_0_yyyy_0 0_yyz_0_xzzz_0 0_yyz_0_xyzz_0 0_yyz_0_xyyz_0 0_yyz_0_xyyy_0 0_yyz_0_xxzz_0 0_yyz_0_xxyz_0 0_yyz_0_xxyy_0 0_yyz_0_xxxz_0 0_yyz_0_xxxy_0 0_yyz_0_xxxx_0 0_yyy_0_zzzz_0 0_yyy_0_yzzz_0 0_yyy_0_yyzz_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xzzz_0 0_yyy_0_xyzz_0 0_yyy_0_xyyz_0 0_yyy_0_xyyy_0 0_yyy_0_xxzz_0 0_yyy_0_xxyz_0 0_yyy_0_xxyy_0 0_yyy_0_xxxz_0 0_yyy_0_xxxy_0 0_yyy_0_xxxx_0 0_xzz_0_zzzz_0 0_xzz_0_yzzz_0 0_xzz_0_yyzz_0 0_xzz_0_yyyz_0 0_xzz_0_yyyy_0 0_xzz_0_xzzz_0 0_xzz_0_xyzz_0 0_xzz_0_xyyz_0 0_xzz_0_xyyy_0 0_xzz_0_xxzz_0 0_xzz_0_xxyz_0 0_xzz_0_xxyy_0 0_xzz_0_xxxz_0 0_xzz_0_xxxy_0 0_xzz_0_xxxx_0 0_xyz_0_zzzz_0 0_xyz_0_yzzz_0 0_xyz_0_yyzz_0 0_xyz_0_yyyz_0 0_xyz_0_yyyy_0 0_xyz_0_xzzz_0 0_xyz_0_xyzz_0 0_xyz_0_xyyz_0 0_xyz_0_xyyy_0 0_xyz_0_xxzz_0 0_xyz_0_xxyz_0 0_xyz_0_xxyy_0 0_xyz_0_xxxz_0 0_xyz_0_xxxy_0 0_xyz_0_xxxx_0 0_xyy_0_zzzz_0 0_xyy_0_yzzz_0 0_xyy_0_yyzz_0 0_xyy_0_yyyz_0 0_xyy_0_yyyy_0 0_xyy_0_xzzz_0 0_xyy_0_xyzz_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxzz_0 0_xyy_0_xxyz_0 0_xyy_0_xxyy_0 0_xyy_0_xxxz_0 0_xyy_0_xxxy_0 0_xyy_0_xxxx_0 0_xxz_0_zzzz_0 0_xxz_0_yzzz_0 0_xxz_0_yyzz_0 0_xxz_0_yyyz_0 0_xxz_0_yyyy_0 0_xxz_0_xzzz_0 0_xxz_0_xyzz_0 0_xxz_0_xyyz_0 0_xxz_0_xyyy_0 0_xxz_0_xxzz_0 0_xxz_0_xxyz_0 0_xxz_0_xxyy_0 0_xxz_0_xxxz_0 0_xxz_0_xxxy_0 0_xxz_0_xxxx_0 0_xxy_0_zzzz_0 0_xxy_0_yzzz_0 0_xxy_0_yyzz_0 0_xxy_0_yyyz_0 0_xxy_0_yyyy_0 0_xxy_0_xzzz_0 0_xxy_0_xyzz_0 0_xxy_0_xyyz_0 0_xxy_0_xyyy_0 0_xxy_0_xxzz_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxxz_0 0_xxy_0_xxxy_0 0_xxy_0_xxxx_0 0_xxx_0_zzzz_0 0_xxx_0_yzzz_0 0_xxx_0_yyzz_0 0_xxx_0_yyyz_0 0_xxx_0_yyyy_0 0_xxx_0_xzzz_0 0_xxx_0_xyzz_0 0_xxx_0_xyyz_0 0_xxx_0_xyyy_0 0_xxx_0_xxzz_0 0_xxx_0_xxyz_0 0_xxx_0_xxyy_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 
    // Contracted integral buffers

    BufferHostXY<T> cbufSDSD(36, ncpairs);

0_zz_0_zz 0_zz_0_yz 0_zz_0_yy 0_zz_0_xz 0_zz_0_xy 0_zz_0_xx 0_yz_0_zz 0_yz_0_yz 0_yz_0_yy 0_yz_0_xz 0_yz_0_xy 0_yz_0_xx 0_yy_0_zz 0_yy_0_yz 0_yy_0_yy 0_yy_0_xz 0_yy_0_xy 0_yy_0_xx 0_xz_0_zz 0_xz_0_yz 0_xz_0_yy 0_xz_0_xz 0_xz_0_xy 0_xz_0_xx 0_xy_0_zz 0_xy_0_yz 0_xy_0_yy 0_xy_0_xz 0_xy_0_xy 0_xy_0_xx 0_xx_0_zz 0_xx_0_yz 0_xx_0_yy 0_xx_0_xz 0_xx_0_xy 0_xx_0_xx 
    BufferHostXY<T> cbufSDSF(60, ncpairs);

0_zz_0_zzz 0_zz_0_yzz 0_zz_0_yyz 0_zz_0_yyy 0_zz_0_xzz 0_zz_0_xyz 0_zz_0_xyy 0_zz_0_xxz 0_zz_0_xxy 0_zz_0_xxx 0_yz_0_zzz 0_yz_0_yzz 0_yz_0_yyz 0_yz_0_yyy 0_yz_0_xzz 0_yz_0_xyz 0_yz_0_xyy 0_yz_0_xxz 0_yz_0_xxy 0_yz_0_xxx 0_yy_0_zzz 0_yy_0_yzz 0_yy_0_yyz 0_yy_0_yyy 0_yy_0_xzz 0_yy_0_xyz 0_yy_0_xyy 0_yy_0_xxz 0_yy_0_xxy 0_yy_0_xxx 0_xz_0_zzz 0_xz_0_yzz 0_xz_0_yyz 0_xz_0_yyy 0_xz_0_xzz 0_xz_0_xyz 0_xz_0_xyy 0_xz_0_xxz 0_xz_0_xxy 0_xz_0_xxx 0_xy_0_zzz 0_xy_0_yzz 0_xy_0_yyz 0_xy_0_yyy 0_xy_0_xzz 0_xy_0_xyz 0_xy_0_xyy 0_xy_0_xxz 0_xy_0_xxy 0_xy_0_xxx 0_xx_0_zzz 0_xx_0_yzz 0_xx_0_yyz 0_xx_0_yyy 0_xx_0_xzz 0_xx_0_xyz 0_xx_0_xyy 0_xx_0_xxz 0_xx_0_xxy 0_xx_0_xxx 
    BufferHostXY<T> cbufSDSG(90, ncpairs);

0_zz_0_zzzz 0_zz_0_yzzz 0_zz_0_yyzz 0_zz_0_yyyz 0_zz_0_yyyy 0_zz_0_xzzz 0_zz_0_xyzz 0_zz_0_xyyz 0_zz_0_xyyy 0_zz_0_xxzz 0_zz_0_xxyz 0_zz_0_xxyy 0_zz_0_xxxz 0_zz_0_xxxy 0_zz_0_xxxx 0_yz_0_zzzz 0_yz_0_yzzz 0_yz_0_yyzz 0_yz_0_yyyz 0_yz_0_yyyy 0_yz_0_xzzz 0_yz_0_xyzz 0_yz_0_xyyz 0_yz_0_xyyy 0_yz_0_xxzz 0_yz_0_xxyz 0_yz_0_xxyy 0_yz_0_xxxz 0_yz_0_xxxy 0_yz_0_xxxx 0_yy_0_zzzz 0_yy_0_yzzz 0_yy_0_yyzz 0_yy_0_yyyz 0_yy_0_yyyy 0_yy_0_xzzz 0_yy_0_xyzz 0_yy_0_xyyz 0_yy_0_xyyy 0_yy_0_xxzz 0_yy_0_xxyz 0_yy_0_xxyy 0_yy_0_xxxz 0_yy_0_xxxy 0_yy_0_xxxx 0_xz_0_zzzz 0_xz_0_yzzz 0_xz_0_yyzz 0_xz_0_yyyz 0_xz_0_yyyy 0_xz_0_xzzz 0_xz_0_xyzz 0_xz_0_xyyz 0_xz_0_xyyy 0_xz_0_xxzz 0_xz_0_xxyz 0_xz_0_xxyy 0_xz_0_xxxz 0_xz_0_xxxy 0_xz_0_xxxx 0_xy_0_zzzz 0_xy_0_yzzz 0_xy_0_yyzz 0_xy_0_yyyz 0_xy_0_yyyy 0_xy_0_xzzz 0_xy_0_xyzz 0_xy_0_xyyz 0_xy_0_xyyy 0_xy_0_xxzz 0_xy_0_xxyz 0_xy_0_xxyy 0_xy_0_xxxz 0_xy_0_xxxy 0_xy_0_xxxx 0_xx_0_zzzz 0_xx_0_yzzz 0_xx_0_yyzz 0_xx_0_yyyz 0_xx_0_yyyy 0_xx_0_xzzz 0_xx_0_xyzz 0_xx_0_xyyz 0_xx_0_xyyy 0_xx_0_xxzz 0_xx_0_xxyz 0_xx_0_xxyy 0_xx_0_xxxz 0_xx_0_xxxy 0_xx_0_xxxx 
    BufferHostXY<T> cbufSFSD(60, ncpairs);

0_zzz_0_zz 0_zzz_0_yz 0_zzz_0_yy 0_zzz_0_xz 0_zzz_0_xy 0_zzz_0_xx 0_yzz_0_zz 0_yzz_0_yz 0_yzz_0_yy 0_yzz_0_xz 0_yzz_0_xy 0_yzz_0_xx 0_yyz_0_zz 0_yyz_0_yz 0_yyz_0_yy 0_yyz_0_xz 0_yyz_0_xy 0_yyz_0_xx 0_yyy_0_zz 0_yyy_0_yz 0_yyy_0_yy 0_yyy_0_xz 0_yyy_0_xy 0_yyy_0_xx 0_xzz_0_zz 0_xzz_0_yz 0_xzz_0_yy 0_xzz_0_xz 0_xzz_0_xy 0_xzz_0_xx 0_xyz_0_zz 0_xyz_0_yz 0_xyz_0_yy 0_xyz_0_xz 0_xyz_0_xy 0_xyz_0_xx 0_xyy_0_zz 0_xyy_0_yz 0_xyy_0_yy 0_xyy_0_xz 0_xyy_0_xy 0_xyy_0_xx 0_xxz_0_zz 0_xxz_0_yz 0_xxz_0_yy 0_xxz_0_xz 0_xxz_0_xy 0_xxz_0_xx 0_xxy_0_zz 0_xxy_0_yz 0_xxy_0_yy 0_xxy_0_xz 0_xxy_0_xy 0_xxy_0_xx 0_xxx_0_zz 0_xxx_0_yz 0_xxx_0_yy 0_xxx_0_xz 0_xxx_0_xy 0_xxx_0_xx 
    BufferHostXY<T> cbufSFSF(100, ncpairs);

0_zzz_0_zzz 0_zzz_0_yzz 0_zzz_0_yyz 0_zzz_0_yyy 0_zzz_0_xzz 0_zzz_0_xyz 0_zzz_0_xyy 0_zzz_0_xxz 0_zzz_0_xxy 0_zzz_0_xxx 0_yzz_0_zzz 0_yzz_0_yzz 0_yzz_0_yyz 0_yzz_0_yyy 0_yzz_0_xzz 0_yzz_0_xyz 0_yzz_0_xyy 0_yzz_0_xxz 0_yzz_0_xxy 0_yzz_0_xxx 0_yyz_0_zzz 0_yyz_0_yzz 0_yyz_0_yyz 0_yyz_0_yyy 0_yyz_0_xzz 0_yyz_0_xyz 0_yyz_0_xyy 0_yyz_0_xxz 0_yyz_0_xxy 0_yyz_0_xxx 0_yyy_0_zzz 0_yyy_0_yzz 0_yyy_0_yyz 0_yyy_0_yyy 0_yyy_0_xzz 0_yyy_0_xyz 0_yyy_0_xyy 0_yyy_0_xxz 0_yyy_0_xxy 0_yyy_0_xxx 0_xzz_0_zzz 0_xzz_0_yzz 0_xzz_0_yyz 0_xzz_0_yyy 0_xzz_0_xzz 0_xzz_0_xyz 0_xzz_0_xyy 0_xzz_0_xxz 0_xzz_0_xxy 0_xzz_0_xxx 0_xyz_0_zzz 0_xyz_0_yzz 0_xyz_0_yyz 0_xyz_0_yyy 0_xyz_0_xzz 0_xyz_0_xyz 0_xyz_0_xyy 0_xyz_0_xxz 0_xyz_0_xxy 0_xyz_0_xxx 0_xyy_0_zzz 0_xyy_0_yzz 0_xyy_0_yyz 0_xyy_0_yyy 0_xyy_0_xzz 0_xyy_0_xyz 0_xyy_0_xyy 0_xyy_0_xxz 0_xyy_0_xxy 0_xyy_0_xxx 0_xxz_0_zzz 0_xxz_0_yzz 0_xxz_0_yyz 0_xxz_0_yyy 0_xxz_0_xzz 0_xxz_0_xyz 0_xxz_0_xyy 0_xxz_0_xxz 0_xxz_0_xxy 0_xxz_0_xxx 0_xxy_0_zzz 0_xxy_0_yzz 0_xxy_0_yyz 0_xxy_0_yyy 0_xxy_0_xzz 0_xxy_0_xyz 0_xxy_0_xyy 0_xxy_0_xxz 0_xxy_0_xxy 0_xxy_0_xxx 0_xxx_0_zzz 0_xxx_0_yzz 0_xxx_0_yyz 0_xxx_0_yyy 0_xxx_0_xzz 0_xxx_0_xyz 0_xxx_0_xyy 0_xxx_0_xxz 0_xxx_0_xxy 0_xxx_0_xxx 
    BufferHostXY<T> cbufSFSG(150, ncpairs);

0_zzz_0_zzzz 0_zzz_0_yzzz 0_zzz_0_yyzz 0_zzz_0_yyyz 0_zzz_0_yyyy 0_zzz_0_xzzz 0_zzz_0_xyzz 0_zzz_0_xyyz 0_zzz_0_xyyy 0_zzz_0_xxzz 0_zzz_0_xxyz 0_zzz_0_xxyy 0_zzz_0_xxxz 0_zzz_0_xxxy 0_zzz_0_xxxx 0_yzz_0_zzzz 0_yzz_0_yzzz 0_yzz_0_yyzz 0_yzz_0_yyyz 0_yzz_0_yyyy 0_yzz_0_xzzz 0_yzz_0_xyzz 0_yzz_0_xyyz 0_yzz_0_xyyy 0_yzz_0_xxzz 0_yzz_0_xxyz 0_yzz_0_xxyy 0_yzz_0_xxxz 0_yzz_0_xxxy 0_yzz_0_xxxx 0_yyz_0_zzzz 0_yyz_0_yzzz 0_yyz_0_yyzz 0_yyz_0_yyyz 0_yyz_0_yyyy 0_yyz_0_xzzz 0_yyz_0_xyzz 0_yyz_0_xyyz 0_yyz_0_xyyy 0_yyz_0_xxzz 0_yyz_0_xxyz 0_yyz_0_xxyy 0_yyz_0_xxxz 0_yyz_0_xxxy 0_yyz_0_xxxx 0_yyy_0_zzzz 0_yyy_0_yzzz 0_yyy_0_yyzz 0_yyy_0_yyyz 0_yyy_0_yyyy 0_yyy_0_xzzz 0_yyy_0_xyzz 0_yyy_0_xyyz 0_yyy_0_xyyy 0_yyy_0_xxzz 0_yyy_0_xxyz 0_yyy_0_xxyy 0_yyy_0_xxxz 0_yyy_0_xxxy 0_yyy_0_xxxx 0_xzz_0_zzzz 0_xzz_0_yzzz 0_xzz_0_yyzz 0_xzz_0_yyyz 0_xzz_0_yyyy 0_xzz_0_xzzz 0_xzz_0_xyzz 0_xzz_0_xyyz 0_xzz_0_xyyy 0_xzz_0_xxzz 0_xzz_0_xxyz 0_xzz_0_xxyy 0_xzz_0_xxxz 0_xzz_0_xxxy 0_xzz_0_xxxx 0_xyz_0_zzzz 0_xyz_0_yzzz 0_xyz_0_yyzz 0_xyz_0_yyyz 0_xyz_0_yyyy 0_xyz_0_xzzz 0_xyz_0_xyzz 0_xyz_0_xyyz 0_xyz_0_xyyy 0_xyz_0_xxzz 0_xyz_0_xxyz 0_xyz_0_xxyy 0_xyz_0_xxxz 0_xyz_0_xxxy 0_xyz_0_xxxx 0_xyy_0_zzzz 0_xyy_0_yzzz 0_xyy_0_yyzz 0_xyy_0_yyyz 0_xyy_0_yyyy 0_xyy_0_xzzz 0_xyy_0_xyzz 0_xyy_0_xyyz 0_xyy_0_xyyy 0_xyy_0_xxzz 0_xyy_0_xxyz 0_xyy_0_xxyy 0_xyy_0_xxxz 0_xyy_0_xxxy 0_xyy_0_xxxx 0_xxz_0_zzzz 0_xxz_0_yzzz 0_xxz_0_yyzz 0_xxz_0_yyyz 0_xxz_0_yyyy 0_xxz_0_xzzz 0_xxz_0_xyzz 0_xxz_0_xyyz 0_xxz_0_xyyy 0_xxz_0_xxzz 0_xxz_0_xxyz 0_xxz_0_xxyy 0_xxz_0_xxxz 0_xxz_0_xxxy 0_xxz_0_xxxx 0_xxy_0_zzzz 0_xxy_0_yzzz 0_xxy_0_yyzz 0_xxy_0_yyyz 0_xxy_0_yyyy 0_xxy_0_xzzz 0_xxy_0_xyzz 0_xxy_0_xyyz 0_xxy_0_xyyy 0_xxy_0_xxzz 0_xxy_0_xxyz 0_xxy_0_xxyy 0_xxy_0_xxxz 0_xxy_0_xxxy 0_xxy_0_xxxx 0_xxx_0_zzzz 0_xxx_0_yzzz 0_xxx_0_yyzz 0_xxx_0_yyyz 0_xxx_0_yyyy 0_xxx_0_xzzz 0_xxx_0_xyzz 0_xxx_0_xyyz 0_xxx_0_xyyy 0_xxx_0_xxzz 0_xxx_0_xxyz 0_xxx_0_xxyy 0_xxx_0_xxxz 0_xxx_0_xxxy 0_xxx_0_xxxx 
    BufferHostXY<T> cbufSDPD(108, ncpairs);

0_zz_z_zz 0_zz_z_yz 0_zz_z_yy 0_zz_z_xz 0_zz_z_xy 0_zz_z_xx 0_zz_y_zz 0_zz_y_yz 0_zz_y_yy 0_zz_y_xz 0_zz_y_xy 0_zz_y_xx 0_zz_x_zz 0_zz_x_yz 0_zz_x_yy 0_zz_x_xz 0_zz_x_xy 0_zz_x_xx 0_yz_z_zz 0_yz_z_yz 0_yz_z_yy 0_yz_z_xz 0_yz_z_xy 0_yz_z_xx 0_yz_y_zz 0_yz_y_yz 0_yz_y_yy 0_yz_y_xz 0_yz_y_xy 0_yz_y_xx 0_yz_x_zz 0_yz_x_yz 0_yz_x_yy 0_yz_x_xz 0_yz_x_xy 0_yz_x_xx 0_yy_z_zz 0_yy_z_yz 0_yy_z_yy 0_yy_z_xz 0_yy_z_xy 0_yy_z_xx 0_yy_y_zz 0_yy_y_yz 0_yy_y_yy 0_yy_y_xz 0_yy_y_xy 0_yy_y_xx 0_yy_x_zz 0_yy_x_yz 0_yy_x_yy 0_yy_x_xz 0_yy_x_xy 0_yy_x_xx 0_xz_z_zz 0_xz_z_yz 0_xz_z_yy 0_xz_z_xz 0_xz_z_xy 0_xz_z_xx 0_xz_y_zz 0_xz_y_yz 0_xz_y_yy 0_xz_y_xz 0_xz_y_xy 0_xz_y_xx 0_xz_x_zz 0_xz_x_yz 0_xz_x_yy 0_xz_x_xz 0_xz_x_xy 0_xz_x_xx 0_xy_z_zz 0_xy_z_yz 0_xy_z_yy 0_xy_z_xz 0_xy_z_xy 0_xy_z_xx 0_xy_y_zz 0_xy_y_yz 0_xy_y_yy 0_xy_y_xz 0_xy_y_xy 0_xy_y_xx 0_xy_x_zz 0_xy_x_yz 0_xy_x_yy 0_xy_x_xz 0_xy_x_xy 0_xy_x_xx 0_xx_z_zz 0_xx_z_yz 0_xx_z_yy 0_xx_z_xz 0_xx_z_xy 0_xx_z_xx 0_xx_y_zz 0_xx_y_yz 0_xx_y_yy 0_xx_y_xz 0_xx_y_xy 0_xx_y_xx 0_xx_x_zz 0_xx_x_yz 0_xx_x_yy 0_xx_x_xz 0_xx_x_xy 0_xx_x_xx 
    BufferHostXY<T> cbufSDPF(150, ncpairs);

0_zz_z_zzz 0_zz_z_yzz 0_zz_z_yyz 0_zz_z_xzz 0_zz_z_xyz 0_zz_z_xxz 0_zz_y_zzz 0_zz_y_yzz 0_zz_y_yyz 0_zz_y_yyy 0_zz_y_xzz 0_zz_y_xyz 0_zz_y_xyy 0_zz_y_xxz 0_zz_y_xxy 0_zz_x_zzz 0_zz_x_yzz 0_zz_x_yyz 0_zz_x_yyy 0_zz_x_xzz 0_zz_x_xyz 0_zz_x_xyy 0_zz_x_xxz 0_zz_x_xxy 0_zz_x_xxx 0_yz_z_zzz 0_yz_z_yzz 0_yz_z_yyz 0_yz_z_xzz 0_yz_z_xyz 0_yz_z_xxz 0_yz_y_zzz 0_yz_y_yzz 0_yz_y_yyz 0_yz_y_yyy 0_yz_y_xzz 0_yz_y_xyz 0_yz_y_xyy 0_yz_y_xxz 0_yz_y_xxy 0_yz_x_zzz 0_yz_x_yzz 0_yz_x_yyz 0_yz_x_yyy 0_yz_x_xzz 0_yz_x_xyz 0_yz_x_xyy 0_yz_x_xxz 0_yz_x_xxy 0_yz_x_xxx 0_yy_z_zzz 0_yy_z_yzz 0_yy_z_yyz 0_yy_z_xzz 0_yy_z_xyz 0_yy_z_xxz 0_yy_y_zzz 0_yy_y_yzz 0_yy_y_yyz 0_yy_y_yyy 0_yy_y_xzz 0_yy_y_xyz 0_yy_y_xyy 0_yy_y_xxz 0_yy_y_xxy 0_yy_x_zzz 0_yy_x_yzz 0_yy_x_yyz 0_yy_x_yyy 0_yy_x_xzz 0_yy_x_xyz 0_yy_x_xyy 0_yy_x_xxz 0_yy_x_xxy 0_yy_x_xxx 0_xz_z_zzz 0_xz_z_yzz 0_xz_z_yyz 0_xz_z_xzz 0_xz_z_xyz 0_xz_z_xxz 0_xz_y_zzz 0_xz_y_yzz 0_xz_y_yyz 0_xz_y_yyy 0_xz_y_xzz 0_xz_y_xyz 0_xz_y_xyy 0_xz_y_xxz 0_xz_y_xxy 0_xz_x_zzz 0_xz_x_yzz 0_xz_x_yyz 0_xz_x_yyy 0_xz_x_xzz 0_xz_x_xyz 0_xz_x_xyy 0_xz_x_xxz 0_xz_x_xxy 0_xz_x_xxx 0_xy_z_zzz 0_xy_z_yzz 0_xy_z_yyz 0_xy_z_xzz 0_xy_z_xyz 0_xy_z_xxz 0_xy_y_zzz 0_xy_y_yzz 0_xy_y_yyz 0_xy_y_yyy 0_xy_y_xzz 0_xy_y_xyz 0_xy_y_xyy 0_xy_y_xxz 0_xy_y_xxy 0_xy_x_zzz 0_xy_x_yzz 0_xy_x_yyz 0_xy_x_yyy 0_xy_x_xzz 0_xy_x_xyz 0_xy_x_xyy 0_xy_x_xxz 0_xy_x_xxy 0_xy_x_xxx 0_xx_z_zzz 0_xx_z_yzz 0_xx_z_yyz 0_xx_z_xzz 0_xx_z_xyz 0_xx_z_xxz 0_xx_y_zzz 0_xx_y_yzz 0_xx_y_yyz 0_xx_y_yyy 0_xx_y_xzz 0_xx_y_xyz 0_xx_y_xyy 0_xx_y_xxz 0_xx_y_xxy 0_xx_x_zzz 0_xx_x_yzz 0_xx_x_yyz 0_xx_x_yyy 0_xx_x_xzz 0_xx_x_xyz 0_xx_x_xyy 0_xx_x_xxz 0_xx_x_xxy 0_xx_x_xxx 
    BufferHostXY<T> cbufSDDD(216, ncpairs);

0_zz_zz_zz 0_zz_zz_yz 0_zz_zz_yy 0_zz_zz_xz 0_zz_zz_xy 0_zz_zz_xx 0_zz_yz_zz 0_zz_yz_yz 0_zz_yz_yy 0_zz_yz_xz 0_zz_yz_xy 0_zz_yz_xx 0_zz_yy_zz 0_zz_yy_yz 0_zz_yy_yy 0_zz_yy_xz 0_zz_yy_xy 0_zz_yy_xx 0_zz_xz_zz 0_zz_xz_yz 0_zz_xz_yy 0_zz_xz_xz 0_zz_xz_xy 0_zz_xz_xx 0_zz_xy_zz 0_zz_xy_yz 0_zz_xy_yy 0_zz_xy_xz 0_zz_xy_xy 0_zz_xy_xx 0_zz_xx_zz 0_zz_xx_yz 0_zz_xx_yy 0_zz_xx_xz 0_zz_xx_xy 0_zz_xx_xx 0_yz_zz_zz 0_yz_zz_yz 0_yz_zz_yy 0_yz_zz_xz 0_yz_zz_xy 0_yz_zz_xx 0_yz_yz_zz 0_yz_yz_yz 0_yz_yz_yy 0_yz_yz_xz 0_yz_yz_xy 0_yz_yz_xx 0_yz_yy_zz 0_yz_yy_yz 0_yz_yy_yy 0_yz_yy_xz 0_yz_yy_xy 0_yz_yy_xx 0_yz_xz_zz 0_yz_xz_yz 0_yz_xz_yy 0_yz_xz_xz 0_yz_xz_xy 0_yz_xz_xx 0_yz_xy_zz 0_yz_xy_yz 0_yz_xy_yy 0_yz_xy_xz 0_yz_xy_xy 0_yz_xy_xx 0_yz_xx_zz 0_yz_xx_yz 0_yz_xx_yy 0_yz_xx_xz 0_yz_xx_xy 0_yz_xx_xx 0_yy_zz_zz 0_yy_zz_yz 0_yy_zz_yy 0_yy_zz_xz 0_yy_zz_xy 0_yy_zz_xx 0_yy_yz_zz 0_yy_yz_yz 0_yy_yz_yy 0_yy_yz_xz 0_yy_yz_xy 0_yy_yz_xx 0_yy_yy_zz 0_yy_yy_yz 0_yy_yy_yy 0_yy_yy_xz 0_yy_yy_xy 0_yy_yy_xx 0_yy_xz_zz 0_yy_xz_yz 0_yy_xz_yy 0_yy_xz_xz 0_yy_xz_xy 0_yy_xz_xx 0_yy_xy_zz 0_yy_xy_yz 0_yy_xy_yy 0_yy_xy_xz 0_yy_xy_xy 0_yy_xy_xx 0_yy_xx_zz 0_yy_xx_yz 0_yy_xx_yy 0_yy_xx_xz 0_yy_xx_xy 0_yy_xx_xx 0_xz_zz_zz 0_xz_zz_yz 0_xz_zz_yy 0_xz_zz_xz 0_xz_zz_xy 0_xz_zz_xx 0_xz_yz_zz 0_xz_yz_yz 0_xz_yz_yy 0_xz_yz_xz 0_xz_yz_xy 0_xz_yz_xx 0_xz_yy_zz 0_xz_yy_yz 0_xz_yy_yy 0_xz_yy_xz 0_xz_yy_xy 0_xz_yy_xx 0_xz_xz_zz 0_xz_xz_yz 0_xz_xz_yy 0_xz_xz_xz 0_xz_xz_xy 0_xz_xz_xx 0_xz_xy_zz 0_xz_xy_yz 0_xz_xy_yy 0_xz_xy_xz 0_xz_xy_xy 0_xz_xy_xx 0_xz_xx_zz 0_xz_xx_yz 0_xz_xx_yy 0_xz_xx_xz 0_xz_xx_xy 0_xz_xx_xx 0_xy_zz_zz 0_xy_zz_yz 0_xy_zz_yy 0_xy_zz_xz 0_xy_zz_xy 0_xy_zz_xx 0_xy_yz_zz 0_xy_yz_yz 0_xy_yz_yy 0_xy_yz_xz 0_xy_yz_xy 0_xy_yz_xx 0_xy_yy_zz 0_xy_yy_yz 0_xy_yy_yy 0_xy_yy_xz 0_xy_yy_xy 0_xy_yy_xx 0_xy_xz_zz 0_xy_xz_yz 0_xy_xz_yy 0_xy_xz_xz 0_xy_xz_xy 0_xy_xz_xx 0_xy_xy_zz 0_xy_xy_yz 0_xy_xy_yy 0_xy_xy_xz 0_xy_xy_xy 0_xy_xy_xx 0_xy_xx_zz 0_xy_xx_yz 0_xy_xx_yy 0_xy_xx_xz 0_xy_xx_xy 0_xy_xx_xx 0_xx_zz_zz 0_xx_zz_yz 0_xx_zz_yy 0_xx_zz_xz 0_xx_zz_xy 0_xx_zz_xx 0_xx_yz_zz 0_xx_yz_yz 0_xx_yz_yy 0_xx_yz_xz 0_xx_yz_xy 0_xx_yz_xx 0_xx_yy_zz 0_xx_yy_yz 0_xx_yy_yy 0_xx_yy_xz 0_xx_yy_xy 0_xx_yy_xx 0_xx_xz_zz 0_xx_xz_yz 0_xx_xz_yy 0_xx_xz_xz 0_xx_xz_xy 0_xx_xz_xx 0_xx_xy_zz 0_xx_xy_yz 0_xx_xy_yy 0_xx_xy_xz 0_xx_xy_xy 0_xx_xy_xx 0_xx_xx_zz 0_xx_xx_yz 0_xx_xx_yy 0_xx_xx_xz 0_xx_xx_xy 0_xx_xx_xx 
    BufferHostXY<T> cbufSFPD(180, ncpairs);

0_zzz_z_zz 0_zzz_z_yz 0_zzz_z_yy 0_zzz_z_xz 0_zzz_z_xy 0_zzz_z_xx 0_zzz_y_zz 0_zzz_y_yz 0_zzz_y_yy 0_zzz_y_xz 0_zzz_y_xy 0_zzz_y_xx 0_zzz_x_zz 0_zzz_x_yz 0_zzz_x_yy 0_zzz_x_xz 0_zzz_x_xy 0_zzz_x_xx 0_yzz_z_zz 0_yzz_z_yz 0_yzz_z_yy 0_yzz_z_xz 0_yzz_z_xy 0_yzz_z_xx 0_yzz_y_zz 0_yzz_y_yz 0_yzz_y_yy 0_yzz_y_xz 0_yzz_y_xy 0_yzz_y_xx 0_yzz_x_zz 0_yzz_x_yz 0_yzz_x_yy 0_yzz_x_xz 0_yzz_x_xy 0_yzz_x_xx 0_yyz_z_zz 0_yyz_z_yz 0_yyz_z_yy 0_yyz_z_xz 0_yyz_z_xy 0_yyz_z_xx 0_yyz_y_zz 0_yyz_y_yz 0_yyz_y_yy 0_yyz_y_xz 0_yyz_y_xy 0_yyz_y_xx 0_yyz_x_zz 0_yyz_x_yz 0_yyz_x_yy 0_yyz_x_xz 0_yyz_x_xy 0_yyz_x_xx 0_yyy_z_zz 0_yyy_z_yz 0_yyy_z_yy 0_yyy_z_xz 0_yyy_z_xy 0_yyy_z_xx 0_yyy_y_zz 0_yyy_y_yz 0_yyy_y_yy 0_yyy_y_xz 0_yyy_y_xy 0_yyy_y_xx 0_yyy_x_zz 0_yyy_x_yz 0_yyy_x_yy 0_yyy_x_xz 0_yyy_x_xy 0_yyy_x_xx 0_xzz_z_zz 0_xzz_z_yz 0_xzz_z_yy 0_xzz_z_xz 0_xzz_z_xy 0_xzz_z_xx 0_xzz_y_zz 0_xzz_y_yz 0_xzz_y_yy 0_xzz_y_xz 0_xzz_y_xy 0_xzz_y_xx 0_xzz_x_zz 0_xzz_x_yz 0_xzz_x_yy 0_xzz_x_xz 0_xzz_x_xy 0_xzz_x_xx 0_xyz_z_zz 0_xyz_z_yz 0_xyz_z_yy 0_xyz_z_xz 0_xyz_z_xy 0_xyz_z_xx 0_xyz_y_zz 0_xyz_y_yz 0_xyz_y_yy 0_xyz_y_xz 0_xyz_y_xy 0_xyz_y_xx 0_xyz_x_zz 0_xyz_x_yz 0_xyz_x_yy 0_xyz_x_xz 0_xyz_x_xy 0_xyz_x_xx 0_xyy_z_zz 0_xyy_z_yz 0_xyy_z_yy 0_xyy_z_xz 0_xyy_z_xy 0_xyy_z_xx 0_xyy_y_zz 0_xyy_y_yz 0_xyy_y_yy 0_xyy_y_xz 0_xyy_y_xy 0_xyy_y_xx 0_xyy_x_zz 0_xyy_x_yz 0_xyy_x_yy 0_xyy_x_xz 0_xyy_x_xy 0_xyy_x_xx 0_xxz_z_zz 0_xxz_z_yz 0_xxz_z_yy 0_xxz_z_xz 0_xxz_z_xy 0_xxz_z_xx 0_xxz_y_zz 0_xxz_y_yz 0_xxz_y_yy 0_xxz_y_xz 0_xxz_y_xy 0_xxz_y_xx 0_xxz_x_zz 0_xxz_x_yz 0_xxz_x_yy 0_xxz_x_xz 0_xxz_x_xy 0_xxz_x_xx 0_xxy_z_zz 0_xxy_z_yz 0_xxy_z_yy 0_xxy_z_xz 0_xxy_z_xy 0_xxy_z_xx 0_xxy_y_zz 0_xxy_y_yz 0_xxy_y_yy 0_xxy_y_xz 0_xxy_y_xy 0_xxy_y_xx 0_xxy_x_zz 0_xxy_x_yz 0_xxy_x_yy 0_xxy_x_xz 0_xxy_x_xy 0_xxy_x_xx 0_xxx_z_zz 0_xxx_z_yz 0_xxx_z_yy 0_xxx_z_xz 0_xxx_z_xy 0_xxx_z_xx 0_xxx_y_zz 0_xxx_y_yz 0_xxx_y_yy 0_xxx_y_xz 0_xxx_y_xy 0_xxx_y_xx 0_xxx_x_zz 0_xxx_x_yz 0_xxx_x_yy 0_xxx_x_xz 0_xxx_x_xy 0_xxx_x_xx 
    BufferHostXY<T> cbufSFPF(250, ncpairs);

0_zzz_z_zzz 0_zzz_z_yzz 0_zzz_z_yyz 0_zzz_z_xzz 0_zzz_z_xyz 0_zzz_z_xxz 0_zzz_y_zzz 0_zzz_y_yzz 0_zzz_y_yyz 0_zzz_y_yyy 0_zzz_y_xzz 0_zzz_y_xyz 0_zzz_y_xyy 0_zzz_y_xxz 0_zzz_y_xxy 0_zzz_x_zzz 0_zzz_x_yzz 0_zzz_x_yyz 0_zzz_x_yyy 0_zzz_x_xzz 0_zzz_x_xyz 0_zzz_x_xyy 0_zzz_x_xxz 0_zzz_x_xxy 0_zzz_x_xxx 0_yzz_z_zzz 0_yzz_z_yzz 0_yzz_z_yyz 0_yzz_z_xzz 0_yzz_z_xyz 0_yzz_z_xxz 0_yzz_y_zzz 0_yzz_y_yzz 0_yzz_y_yyz 0_yzz_y_yyy 0_yzz_y_xzz 0_yzz_y_xyz 0_yzz_y_xyy 0_yzz_y_xxz 0_yzz_y_xxy 0_yzz_x_zzz 0_yzz_x_yzz 0_yzz_x_yyz 0_yzz_x_yyy 0_yzz_x_xzz 0_yzz_x_xyz 0_yzz_x_xyy 0_yzz_x_xxz 0_yzz_x_xxy 0_yzz_x_xxx 0_yyz_z_zzz 0_yyz_z_yzz 0_yyz_z_yyz 0_yyz_z_xzz 0_yyz_z_xyz 0_yyz_z_xxz 0_yyz_y_zzz 0_yyz_y_yzz 0_yyz_y_yyz 0_yyz_y_yyy 0_yyz_y_xzz 0_yyz_y_xyz 0_yyz_y_xyy 0_yyz_y_xxz 0_yyz_y_xxy 0_yyz_x_zzz 0_yyz_x_yzz 0_yyz_x_yyz 0_yyz_x_yyy 0_yyz_x_xzz 0_yyz_x_xyz 0_yyz_x_xyy 0_yyz_x_xxz 0_yyz_x_xxy 0_yyz_x_xxx 0_yyy_z_zzz 0_yyy_z_yzz 0_yyy_z_yyz 0_yyy_z_xzz 0_yyy_z_xyz 0_yyy_z_xxz 0_yyy_y_zzz 0_yyy_y_yzz 0_yyy_y_yyz 0_yyy_y_yyy 0_yyy_y_xzz 0_yyy_y_xyz 0_yyy_y_xyy 0_yyy_y_xxz 0_yyy_y_xxy 0_yyy_x_zzz 0_yyy_x_yzz 0_yyy_x_yyz 0_yyy_x_yyy 0_yyy_x_xzz 0_yyy_x_xyz 0_yyy_x_xyy 0_yyy_x_xxz 0_yyy_x_xxy 0_yyy_x_xxx 0_xzz_z_zzz 0_xzz_z_yzz 0_xzz_z_yyz 0_xzz_z_xzz 0_xzz_z_xyz 0_xzz_z_xxz 0_xzz_y_zzz 0_xzz_y_yzz 0_xzz_y_yyz 0_xzz_y_yyy 0_xzz_y_xzz 0_xzz_y_xyz 0_xzz_y_xyy 0_xzz_y_xxz 0_xzz_y_xxy 0_xzz_x_zzz 0_xzz_x_yzz 0_xzz_x_yyz 0_xzz_x_yyy 0_xzz_x_xzz 0_xzz_x_xyz 0_xzz_x_xyy 0_xzz_x_xxz 0_xzz_x_xxy 0_xzz_x_xxx 0_xyz_z_zzz 0_xyz_z_yzz 0_xyz_z_yyz 0_xyz_z_xzz 0_xyz_z_xyz 0_xyz_z_xxz 0_xyz_y_zzz 0_xyz_y_yzz 0_xyz_y_yyz 0_xyz_y_yyy 0_xyz_y_xzz 0_xyz_y_xyz 0_xyz_y_xyy 0_xyz_y_xxz 0_xyz_y_xxy 0_xyz_x_zzz 0_xyz_x_yzz 0_xyz_x_yyz 0_xyz_x_yyy 0_xyz_x_xzz 0_xyz_x_xyz 0_xyz_x_xyy 0_xyz_x_xxz 0_xyz_x_xxy 0_xyz_x_xxx 0_xyy_z_zzz 0_xyy_z_yzz 0_xyy_z_yyz 0_xyy_z_xzz 0_xyy_z_xyz 0_xyy_z_xxz 0_xyy_y_zzz 0_xyy_y_yzz 0_xyy_y_yyz 0_xyy_y_yyy 0_xyy_y_xzz 0_xyy_y_xyz 0_xyy_y_xyy 0_xyy_y_xxz 0_xyy_y_xxy 0_xyy_x_zzz 0_xyy_x_yzz 0_xyy_x_yyz 0_xyy_x_yyy 0_xyy_x_xzz 0_xyy_x_xyz 0_xyy_x_xyy 0_xyy_x_xxz 0_xyy_x_xxy 0_xyy_x_xxx 0_xxz_z_zzz 0_xxz_z_yzz 0_xxz_z_yyz 0_xxz_z_xzz 0_xxz_z_xyz 0_xxz_z_xxz 0_xxz_y_zzz 0_xxz_y_yzz 0_xxz_y_yyz 0_xxz_y_yyy 0_xxz_y_xzz 0_xxz_y_xyz 0_xxz_y_xyy 0_xxz_y_xxz 0_xxz_y_xxy 0_xxz_x_zzz 0_xxz_x_yzz 0_xxz_x_yyz 0_xxz_x_yyy 0_xxz_x_xzz 0_xxz_x_xyz 0_xxz_x_xyy 0_xxz_x_xxz 0_xxz_x_xxy 0_xxz_x_xxx 0_xxy_z_zzz 0_xxy_z_yzz 0_xxy_z_yyz 0_xxy_z_xzz 0_xxy_z_xyz 0_xxy_z_xxz 0_xxy_y_zzz 0_xxy_y_yzz 0_xxy_y_yyz 0_xxy_y_yyy 0_xxy_y_xzz 0_xxy_y_xyz 0_xxy_y_xyy 0_xxy_y_xxz 0_xxy_y_xxy 0_xxy_x_zzz 0_xxy_x_yzz 0_xxy_x_yyz 0_xxy_x_yyy 0_xxy_x_xzz 0_xxy_x_xyz 0_xxy_x_xyy 0_xxy_x_xxz 0_xxy_x_xxy 0_xxy_x_xxx 0_xxx_z_zzz 0_xxx_z_yzz 0_xxx_z_yyz 0_xxx_z_xzz 0_xxx_z_xyz 0_xxx_z_xxz 0_xxx_y_zzz 0_xxx_y_yzz 0_xxx_y_yyz 0_xxx_y_yyy 0_xxx_y_xzz 0_xxx_y_xyz 0_xxx_y_xyy 0_xxx_y_xxz 0_xxx_y_xxy 0_xxx_x_zzz 0_xxx_x_yzz 0_xxx_x_yyz 0_xxx_x_yyy 0_xxx_x_xzz 0_xxx_x_xyz 0_xxx_x_xyy 0_xxx_x_xxz 0_xxx_x_xxy 0_xxx_x_xxx 
    BufferHostXY<T> cbufSFDD(360, ncpairs);

0_zzz_zz_zz 0_zzz_zz_yz 0_zzz_zz_yy 0_zzz_zz_xz 0_zzz_zz_xy 0_zzz_zz_xx 0_zzz_yz_zz 0_zzz_yz_yz 0_zzz_yz_yy 0_zzz_yz_xz 0_zzz_yz_xy 0_zzz_yz_xx 0_zzz_yy_zz 0_zzz_yy_yz 0_zzz_yy_yy 0_zzz_yy_xz 0_zzz_yy_xy 0_zzz_yy_xx 0_zzz_xz_zz 0_zzz_xz_yz 0_zzz_xz_yy 0_zzz_xz_xz 0_zzz_xz_xy 0_zzz_xz_xx 0_zzz_xy_zz 0_zzz_xy_yz 0_zzz_xy_yy 0_zzz_xy_xz 0_zzz_xy_xy 0_zzz_xy_xx 0_zzz_xx_zz 0_zzz_xx_yz 0_zzz_xx_yy 0_zzz_xx_xz 0_zzz_xx_xy 0_zzz_xx_xx 0_yzz_zz_zz 0_yzz_zz_yz 0_yzz_zz_yy 0_yzz_zz_xz 0_yzz_zz_xy 0_yzz_zz_xx 0_yzz_yz_zz 0_yzz_yz_yz 0_yzz_yz_yy 0_yzz_yz_xz 0_yzz_yz_xy 0_yzz_yz_xx 0_yzz_yy_zz 0_yzz_yy_yz 0_yzz_yy_yy 0_yzz_yy_xz 0_yzz_yy_xy 0_yzz_yy_xx 0_yzz_xz_zz 0_yzz_xz_yz 0_yzz_xz_yy 0_yzz_xz_xz 0_yzz_xz_xy 0_yzz_xz_xx 0_yzz_xy_zz 0_yzz_xy_yz 0_yzz_xy_yy 0_yzz_xy_xz 0_yzz_xy_xy 0_yzz_xy_xx 0_yzz_xx_zz 0_yzz_xx_yz 0_yzz_xx_yy 0_yzz_xx_xz 0_yzz_xx_xy 0_yzz_xx_xx 0_yyz_zz_zz 0_yyz_zz_yz 0_yyz_zz_yy 0_yyz_zz_xz 0_yyz_zz_xy 0_yyz_zz_xx 0_yyz_yz_zz 0_yyz_yz_yz 0_yyz_yz_yy 0_yyz_yz_xz 0_yyz_yz_xy 0_yyz_yz_xx 0_yyz_yy_zz 0_yyz_yy_yz 0_yyz_yy_yy 0_yyz_yy_xz 0_yyz_yy_xy 0_yyz_yy_xx 0_yyz_xz_zz 0_yyz_xz_yz 0_yyz_xz_yy 0_yyz_xz_xz 0_yyz_xz_xy 0_yyz_xz_xx 0_yyz_xy_zz 0_yyz_xy_yz 0_yyz_xy_yy 0_yyz_xy_xz 0_yyz_xy_xy 0_yyz_xy_xx 0_yyz_xx_zz 0_yyz_xx_yz 0_yyz_xx_yy 0_yyz_xx_xz 0_yyz_xx_xy 0_yyz_xx_xx 0_yyy_zz_zz 0_yyy_zz_yz 0_yyy_zz_yy 0_yyy_zz_xz 0_yyy_zz_xy 0_yyy_zz_xx 0_yyy_yz_zz 0_yyy_yz_yz 0_yyy_yz_yy 0_yyy_yz_xz 0_yyy_yz_xy 0_yyy_yz_xx 0_yyy_yy_zz 0_yyy_yy_yz 0_yyy_yy_yy 0_yyy_yy_xz 0_yyy_yy_xy 0_yyy_yy_xx 0_yyy_xz_zz 0_yyy_xz_yz 0_yyy_xz_yy 0_yyy_xz_xz 0_yyy_xz_xy 0_yyy_xz_xx 0_yyy_xy_zz 0_yyy_xy_yz 0_yyy_xy_yy 0_yyy_xy_xz 0_yyy_xy_xy 0_yyy_xy_xx 0_yyy_xx_zz 0_yyy_xx_yz 0_yyy_xx_yy 0_yyy_xx_xz 0_yyy_xx_xy 0_yyy_xx_xx 0_xzz_zz_zz 0_xzz_zz_yz 0_xzz_zz_yy 0_xzz_zz_xz 0_xzz_zz_xy 0_xzz_zz_xx 0_xzz_yz_zz 0_xzz_yz_yz 0_xzz_yz_yy 0_xzz_yz_xz 0_xzz_yz_xy 0_xzz_yz_xx 0_xzz_yy_zz 0_xzz_yy_yz 0_xzz_yy_yy 0_xzz_yy_xz 0_xzz_yy_xy 0_xzz_yy_xx 0_xzz_xz_zz 0_xzz_xz_yz 0_xzz_xz_yy 0_xzz_xz_xz 0_xzz_xz_xy 0_xzz_xz_xx 0_xzz_xy_zz 0_xzz_xy_yz 0_xzz_xy_yy 0_xzz_xy_xz 0_xzz_xy_xy 0_xzz_xy_xx 0_xzz_xx_zz 0_xzz_xx_yz 0_xzz_xx_yy 0_xzz_xx_xz 0_xzz_xx_xy 0_xzz_xx_xx 0_xyz_zz_zz 0_xyz_zz_yz 0_xyz_zz_yy 0_xyz_zz_xz 0_xyz_zz_xy 0_xyz_zz_xx 0_xyz_yz_zz 0_xyz_yz_yz 0_xyz_yz_yy 0_xyz_yz_xz 0_xyz_yz_xy 0_xyz_yz_xx 0_xyz_yy_zz 0_xyz_yy_yz 0_xyz_yy_yy 0_xyz_yy_xz 0_xyz_yy_xy 0_xyz_yy_xx 0_xyz_xz_zz 0_xyz_xz_yz 0_xyz_xz_yy 0_xyz_xz_xz 0_xyz_xz_xy 0_xyz_xz_xx 0_xyz_xy_zz 0_xyz_xy_yz 0_xyz_xy_yy 0_xyz_xy_xz 0_xyz_xy_xy 0_xyz_xy_xx 0_xyz_xx_zz 0_xyz_xx_yz 0_xyz_xx_yy 0_xyz_xx_xz 0_xyz_xx_xy 0_xyz_xx_xx 0_xyy_zz_zz 0_xyy_zz_yz 0_xyy_zz_yy 0_xyy_zz_xz 0_xyy_zz_xy 0_xyy_zz_xx 0_xyy_yz_zz 0_xyy_yz_yz 0_xyy_yz_yy 0_xyy_yz_xz 0_xyy_yz_xy 0_xyy_yz_xx 0_xyy_yy_zz 0_xyy_yy_yz 0_xyy_yy_yy 0_xyy_yy_xz 0_xyy_yy_xy 0_xyy_yy_xx 0_xyy_xz_zz 0_xyy_xz_yz 0_xyy_xz_yy 0_xyy_xz_xz 0_xyy_xz_xy 0_xyy_xz_xx 0_xyy_xy_zz 0_xyy_xy_yz 0_xyy_xy_yy 0_xyy_xy_xz 0_xyy_xy_xy 0_xyy_xy_xx 0_xyy_xx_zz 0_xyy_xx_yz 0_xyy_xx_yy 0_xyy_xx_xz 0_xyy_xx_xy 0_xyy_xx_xx 0_xxz_zz_zz 0_xxz_zz_yz 0_xxz_zz_yy 0_xxz_zz_xz 0_xxz_zz_xy 0_xxz_zz_xx 0_xxz_yz_zz 0_xxz_yz_yz 0_xxz_yz_yy 0_xxz_yz_xz 0_xxz_yz_xy 0_xxz_yz_xx 0_xxz_yy_zz 0_xxz_yy_yz 0_xxz_yy_yy 0_xxz_yy_xz 0_xxz_yy_xy 0_xxz_yy_xx 0_xxz_xz_zz 0_xxz_xz_yz 0_xxz_xz_yy 0_xxz_xz_xz 0_xxz_xz_xy 0_xxz_xz_xx 0_xxz_xy_zz 0_xxz_xy_yz 0_xxz_xy_yy 0_xxz_xy_xz 0_xxz_xy_xy 0_xxz_xy_xx 0_xxz_xx_zz 0_xxz_xx_yz 0_xxz_xx_yy 0_xxz_xx_xz 0_xxz_xx_xy 0_xxz_xx_xx 0_xxy_zz_zz 0_xxy_zz_yz 0_xxy_zz_yy 0_xxy_zz_xz 0_xxy_zz_xy 0_xxy_zz_xx 0_xxy_yz_zz 0_xxy_yz_yz 0_xxy_yz_yy 0_xxy_yz_xz 0_xxy_yz_xy 0_xxy_yz_xx 0_xxy_yy_zz 0_xxy_yy_yz 0_xxy_yy_yy 0_xxy_yy_xz 0_xxy_yy_xy 0_xxy_yy_xx 0_xxy_xz_zz 0_xxy_xz_yz 0_xxy_xz_yy 0_xxy_xz_xz 0_xxy_xz_xy 0_xxy_xz_xx 0_xxy_xy_zz 0_xxy_xy_yz 0_xxy_xy_yy 0_xxy_xy_xz 0_xxy_xy_xy 0_xxy_xy_xx 0_xxy_xx_zz 0_xxy_xx_yz 0_xxy_xx_yy 0_xxy_xx_xz 0_xxy_xx_xy 0_xxy_xx_xx 0_xxx_zz_zz 0_xxx_zz_yz 0_xxx_zz_yy 0_xxx_zz_xz 0_xxx_zz_xy 0_xxx_zz_xx 0_xxx_yz_zz 0_xxx_yz_yz 0_xxx_yz_yy 0_xxx_yz_xz 0_xxx_yz_xy 0_xxx_yz_xx 0_xxx_yy_zz 0_xxx_yy_yz 0_xxx_yy_yy 0_xxx_yy_xz 0_xxx_yy_xy 0_xxx_yy_xx 0_xxx_xz_zz 0_xxx_xz_yz 0_xxx_xz_yy 0_xxx_xz_xz 0_xxx_xz_xy 0_xxx_xz_xx 0_xxx_xy_zz 0_xxx_xy_yz 0_xxx_xy_yy 0_xxx_xy_xz 0_xxx_xy_xy 0_xxx_xy_xx 0_xxx_xx_zz 0_xxx_xx_yz 0_xxx_xx_yy 0_xxx_xx_xz 0_xxx_xx_xy 0_xxx_xx_xx 
    BufferHostXY<T> cbufPDDD(648, ncpairs);

z_zz_zz_zz z_zz_zz_yz z_zz_zz_yy z_zz_zz_xz z_zz_zz_xy z_zz_zz_xx z_zz_yz_zz z_zz_yz_yz z_zz_yz_yy z_zz_yz_xz z_zz_yz_xy z_zz_yz_xx z_zz_yy_zz z_zz_yy_yz z_zz_yy_yy z_zz_yy_xz z_zz_yy_xy z_zz_yy_xx z_zz_xz_zz z_zz_xz_yz z_zz_xz_yy z_zz_xz_xz z_zz_xz_xy z_zz_xz_xx z_zz_xy_zz z_zz_xy_yz z_zz_xy_yy z_zz_xy_xz z_zz_xy_xy z_zz_xy_xx z_zz_xx_zz z_zz_xx_yz z_zz_xx_yy z_zz_xx_xz z_zz_xx_xy z_zz_xx_xx z_yz_zz_zz z_yz_zz_yz z_yz_zz_yy z_yz_zz_xz z_yz_zz_xy z_yz_zz_xx z_yz_yz_zz z_yz_yz_yz z_yz_yz_yy z_yz_yz_xz z_yz_yz_xy z_yz_yz_xx z_yz_yy_zz z_yz_yy_yz z_yz_yy_yy z_yz_yy_xz z_yz_yy_xy z_yz_yy_xx z_yz_xz_zz z_yz_xz_yz z_yz_xz_yy z_yz_xz_xz z_yz_xz_xy z_yz_xz_xx z_yz_xy_zz z_yz_xy_yz z_yz_xy_yy z_yz_xy_xz z_yz_xy_xy z_yz_xy_xx z_yz_xx_zz z_yz_xx_yz z_yz_xx_yy z_yz_xx_xz z_yz_xx_xy z_yz_xx_xx z_yy_zz_zz z_yy_zz_yz z_yy_zz_yy z_yy_zz_xz z_yy_zz_xy z_yy_zz_xx z_yy_yz_zz z_yy_yz_yz z_yy_yz_yy z_yy_yz_xz z_yy_yz_xy z_yy_yz_xx z_yy_yy_zz z_yy_yy_yz z_yy_yy_yy z_yy_yy_xz z_yy_yy_xy z_yy_yy_xx z_yy_xz_zz z_yy_xz_yz z_yy_xz_yy z_yy_xz_xz z_yy_xz_xy z_yy_xz_xx z_yy_xy_zz z_yy_xy_yz z_yy_xy_yy z_yy_xy_xz z_yy_xy_xy z_yy_xy_xx z_yy_xx_zz z_yy_xx_yz z_yy_xx_yy z_yy_xx_xz z_yy_xx_xy z_yy_xx_xx z_xz_zz_zz z_xz_zz_yz z_xz_zz_yy z_xz_zz_xz z_xz_zz_xy z_xz_zz_xx z_xz_yz_zz z_xz_yz_yz z_xz_yz_yy z_xz_yz_xz z_xz_yz_xy z_xz_yz_xx z_xz_yy_zz z_xz_yy_yz z_xz_yy_yy z_xz_yy_xz z_xz_yy_xy z_xz_yy_xx z_xz_xz_zz z_xz_xz_yz z_xz_xz_yy z_xz_xz_xz z_xz_xz_xy z_xz_xz_xx z_xz_xy_zz z_xz_xy_yz z_xz_xy_yy z_xz_xy_xz z_xz_xy_xy z_xz_xy_xx z_xz_xx_zz z_xz_xx_yz z_xz_xx_yy z_xz_xx_xz z_xz_xx_xy z_xz_xx_xx z_xy_zz_zz z_xy_zz_yz z_xy_zz_yy z_xy_zz_xz z_xy_zz_xy z_xy_zz_xx z_xy_yz_zz z_xy_yz_yz z_xy_yz_yy z_xy_yz_xz z_xy_yz_xy z_xy_yz_xx z_xy_yy_zz z_xy_yy_yz z_xy_yy_yy z_xy_yy_xz z_xy_yy_xy z_xy_yy_xx z_xy_xz_zz z_xy_xz_yz z_xy_xz_yy z_xy_xz_xz z_xy_xz_xy z_xy_xz_xx z_xy_xy_zz z_xy_xy_yz z_xy_xy_yy z_xy_xy_xz z_xy_xy_xy z_xy_xy_xx z_xy_xx_zz z_xy_xx_yz z_xy_xx_yy z_xy_xx_xz z_xy_xx_xy z_xy_xx_xx z_xx_zz_zz z_xx_zz_yz z_xx_zz_yy z_xx_zz_xz z_xx_zz_xy z_xx_zz_xx z_xx_yz_zz z_xx_yz_yz z_xx_yz_yy z_xx_yz_xz z_xx_yz_xy z_xx_yz_xx z_xx_yy_zz z_xx_yy_yz z_xx_yy_yy z_xx_yy_xz z_xx_yy_xy z_xx_yy_xx z_xx_xz_zz z_xx_xz_yz z_xx_xz_yy z_xx_xz_xz z_xx_xz_xy z_xx_xz_xx z_xx_xy_zz z_xx_xy_yz z_xx_xy_yy z_xx_xy_xz z_xx_xy_xy z_xx_xy_xx z_xx_xx_zz z_xx_xx_yz z_xx_xx_yy z_xx_xx_xz z_xx_xx_xy z_xx_xx_xx y_zz_zz_zz y_zz_zz_yz y_zz_zz_yy y_zz_zz_xz y_zz_zz_xy y_zz_zz_xx y_zz_yz_zz y_zz_yz_yz y_zz_yz_yy y_zz_yz_xz y_zz_yz_xy y_zz_yz_xx y_zz_yy_zz y_zz_yy_yz y_zz_yy_yy y_zz_yy_xz y_zz_yy_xy y_zz_yy_xx y_zz_xz_zz y_zz_xz_yz y_zz_xz_yy y_zz_xz_xz y_zz_xz_xy y_zz_xz_xx y_zz_xy_zz y_zz_xy_yz y_zz_xy_yy y_zz_xy_xz y_zz_xy_xy y_zz_xy_xx y_zz_xx_zz y_zz_xx_yz y_zz_xx_yy y_zz_xx_xz y_zz_xx_xy y_zz_xx_xx y_yz_zz_zz y_yz_zz_yz y_yz_zz_yy y_yz_zz_xz y_yz_zz_xy y_yz_zz_xx y_yz_yz_zz y_yz_yz_yz y_yz_yz_yy y_yz_yz_xz y_yz_yz_xy y_yz_yz_xx y_yz_yy_zz y_yz_yy_yz y_yz_yy_yy y_yz_yy_xz y_yz_yy_xy y_yz_yy_xx y_yz_xz_zz y_yz_xz_yz y_yz_xz_yy y_yz_xz_xz y_yz_xz_xy y_yz_xz_xx y_yz_xy_zz y_yz_xy_yz y_yz_xy_yy y_yz_xy_xz y_yz_xy_xy y_yz_xy_xx y_yz_xx_zz y_yz_xx_yz y_yz_xx_yy y_yz_xx_xz y_yz_xx_xy y_yz_xx_xx y_yy_zz_zz y_yy_zz_yz y_yy_zz_yy y_yy_zz_xz y_yy_zz_xy y_yy_zz_xx y_yy_yz_zz y_yy_yz_yz y_yy_yz_yy y_yy_yz_xz y_yy_yz_xy y_yy_yz_xx y_yy_yy_zz y_yy_yy_yz y_yy_yy_yy y_yy_yy_xz y_yy_yy_xy y_yy_yy_xx y_yy_xz_zz y_yy_xz_yz y_yy_xz_yy y_yy_xz_xz y_yy_xz_xy y_yy_xz_xx y_yy_xy_zz y_yy_xy_yz y_yy_xy_yy y_yy_xy_xz y_yy_xy_xy y_yy_xy_xx y_yy_xx_zz y_yy_xx_yz y_yy_xx_yy y_yy_xx_xz y_yy_xx_xy y_yy_xx_xx y_xz_zz_zz y_xz_zz_yz y_xz_zz_yy y_xz_zz_xz y_xz_zz_xy y_xz_zz_xx y_xz_yz_zz y_xz_yz_yz y_xz_yz_yy y_xz_yz_xz y_xz_yz_xy y_xz_yz_xx y_xz_yy_zz y_xz_yy_yz y_xz_yy_yy y_xz_yy_xz y_xz_yy_xy y_xz_yy_xx y_xz_xz_zz y_xz_xz_yz y_xz_xz_yy y_xz_xz_xz y_xz_xz_xy y_xz_xz_xx y_xz_xy_zz y_xz_xy_yz y_xz_xy_yy y_xz_xy_xz y_xz_xy_xy y_xz_xy_xx y_xz_xx_zz y_xz_xx_yz y_xz_xx_yy y_xz_xx_xz y_xz_xx_xy y_xz_xx_xx y_xy_zz_zz y_xy_zz_yz y_xy_zz_yy y_xy_zz_xz y_xy_zz_xy y_xy_zz_xx y_xy_yz_zz y_xy_yz_yz y_xy_yz_yy y_xy_yz_xz y_xy_yz_xy y_xy_yz_xx y_xy_yy_zz y_xy_yy_yz y_xy_yy_yy y_xy_yy_xz y_xy_yy_xy y_xy_yy_xx y_xy_xz_zz y_xy_xz_yz y_xy_xz_yy y_xy_xz_xz y_xy_xz_xy y_xy_xz_xx y_xy_xy_zz y_xy_xy_yz y_xy_xy_yy y_xy_xy_xz y_xy_xy_xy y_xy_xy_xx y_xy_xx_zz y_xy_xx_yz y_xy_xx_yy y_xy_xx_xz y_xy_xx_xy y_xy_xx_xx y_xx_zz_zz y_xx_zz_yz y_xx_zz_yy y_xx_zz_xz y_xx_zz_xy y_xx_zz_xx y_xx_yz_zz y_xx_yz_yz y_xx_yz_yy y_xx_yz_xz y_xx_yz_xy y_xx_yz_xx y_xx_yy_zz y_xx_yy_yz y_xx_yy_yy y_xx_yy_xz y_xx_yy_xy y_xx_yy_xx y_xx_xz_zz y_xx_xz_yz y_xx_xz_yy y_xx_xz_xz y_xx_xz_xy y_xx_xz_xx y_xx_xy_zz y_xx_xy_yz y_xx_xy_yy y_xx_xy_xz y_xx_xy_xy y_xx_xy_xx y_xx_xx_zz y_xx_xx_yz y_xx_xx_yy y_xx_xx_xz y_xx_xx_xy y_xx_xx_xx x_zz_zz_zz x_zz_zz_yz x_zz_zz_yy x_zz_zz_xz x_zz_zz_xy x_zz_zz_xx x_zz_yz_zz x_zz_yz_yz x_zz_yz_yy x_zz_yz_xz x_zz_yz_xy x_zz_yz_xx x_zz_yy_zz x_zz_yy_yz x_zz_yy_yy x_zz_yy_xz x_zz_yy_xy x_zz_yy_xx x_zz_xz_zz x_zz_xz_yz x_zz_xz_yy x_zz_xz_xz x_zz_xz_xy x_zz_xz_xx x_zz_xy_zz x_zz_xy_yz x_zz_xy_yy x_zz_xy_xz x_zz_xy_xy x_zz_xy_xx x_zz_xx_zz x_zz_xx_yz x_zz_xx_yy x_zz_xx_xz x_zz_xx_xy x_zz_xx_xx x_yz_zz_zz x_yz_zz_yz x_yz_zz_yy x_yz_zz_xz x_yz_zz_xy x_yz_zz_xx x_yz_yz_zz x_yz_yz_yz x_yz_yz_yy x_yz_yz_xz x_yz_yz_xy x_yz_yz_xx x_yz_yy_zz x_yz_yy_yz x_yz_yy_yy x_yz_yy_xz x_yz_yy_xy x_yz_yy_xx x_yz_xz_zz x_yz_xz_yz x_yz_xz_yy x_yz_xz_xz x_yz_xz_xy x_yz_xz_xx x_yz_xy_zz x_yz_xy_yz x_yz_xy_yy x_yz_xy_xz x_yz_xy_xy x_yz_xy_xx x_yz_xx_zz x_yz_xx_yz x_yz_xx_yy x_yz_xx_xz x_yz_xx_xy x_yz_xx_xx x_yy_zz_zz x_yy_zz_yz x_yy_zz_yy x_yy_zz_xz x_yy_zz_xy x_yy_zz_xx x_yy_yz_zz x_yy_yz_yz x_yy_yz_yy x_yy_yz_xz x_yy_yz_xy x_yy_yz_xx x_yy_yy_zz x_yy_yy_yz x_yy_yy_yy x_yy_yy_xz x_yy_yy_xy x_yy_yy_xx x_yy_xz_zz x_yy_xz_yz x_yy_xz_yy x_yy_xz_xz x_yy_xz_xy x_yy_xz_xx x_yy_xy_zz x_yy_xy_yz x_yy_xy_yy x_yy_xy_xz x_yy_xy_xy x_yy_xy_xx x_yy_xx_zz x_yy_xx_yz x_yy_xx_yy x_yy_xx_xz x_yy_xx_xy x_yy_xx_xx x_xz_zz_zz x_xz_zz_yz x_xz_zz_yy x_xz_zz_xz x_xz_zz_xy x_xz_zz_xx x_xz_yz_zz x_xz_yz_yz x_xz_yz_yy x_xz_yz_xz x_xz_yz_xy x_xz_yz_xx x_xz_yy_zz x_xz_yy_yz x_xz_yy_yy x_xz_yy_xz x_xz_yy_xy x_xz_yy_xx x_xz_xz_zz x_xz_xz_yz x_xz_xz_yy x_xz_xz_xz x_xz_xz_xy x_xz_xz_xx x_xz_xy_zz x_xz_xy_yz x_xz_xy_yy x_xz_xy_xz x_xz_xy_xy x_xz_xy_xx x_xz_xx_zz x_xz_xx_yz x_xz_xx_yy x_xz_xx_xz x_xz_xx_xy x_xz_xx_xx x_xy_zz_zz x_xy_zz_yz x_xy_zz_yy x_xy_zz_xz x_xy_zz_xy x_xy_zz_xx x_xy_yz_zz x_xy_yz_yz x_xy_yz_yy x_xy_yz_xz x_xy_yz_xy x_xy_yz_xx x_xy_yy_zz x_xy_yy_yz x_xy_yy_yy x_xy_yy_xz x_xy_yy_xy x_xy_yy_xx x_xy_xz_zz x_xy_xz_yz x_xy_xz_yy x_xy_xz_xz x_xy_xz_xy x_xy_xz_xx x_xy_xy_zz x_xy_xy_yz x_xy_xy_yy x_xy_xy_xz x_xy_xy_xy x_xy_xy_xx x_xy_xx_zz x_xy_xx_yz x_xy_xx_yy x_xy_xx_xz x_xy_xx_xy x_xy_xx_xx x_xx_zz_zz x_xx_zz_yz x_xx_zz_yy x_xx_zz_xz x_xx_zz_xy x_xx_zz_xx x_xx_yz_zz x_xx_yz_yz x_xx_yz_yy x_xx_yz_xz x_xx_yz_xy x_xx_yz_xx x_xx_yy_zz x_xx_yy_yz x_xx_yy_yy x_xx_yy_xz x_xx_yy_xy x_xx_yy_xx x_xx_xz_zz x_xx_xz_yz x_xx_xz_yy x_xx_xz_xz x_xx_xz_xy x_xx_xz_xx x_xx_xy_zz x_xx_xy_yz x_xx_xy_yy x_xx_xy_xz x_xx_xy_xy x_xx_xy_xx x_xx_xx_zz x_xx_xx_yz x_xx_xx_yy x_xx_xx_xz x_xx_xx_xy x_xx_xx_xx 
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
0_zz_zz_zz_0 0_zz_zz_yz_0 0_zz_zz_yy_0 0_zz_zz_xz_0 0_zz_zz_xy_0 0_zz_zz_xx_0 0_zz_yz_zz_0 0_zz_yz_yz_0 0_zz_yz_yy_0 0_zz_yz_xz_0 0_zz_yz_xy_0 0_zz_yz_xx_0 0_zz_yy_zz_0 0_zz_yy_yz_0 0_zz_yy_yy_0 0_zz_yy_xz_0 0_zz_yy_xy_0 0_zz_yy_xx_0 0_zz_xz_zz_0 0_zz_xz_yz_0 0_zz_xz_yy_0 0_zz_xz_xz_0 0_zz_xz_xy_0 0_zz_xz_xx_0 0_zz_xy_zz_0 0_zz_xy_yz_0 0_zz_xy_yy_0 0_zz_xy_xz_0 0_zz_xy_xy_0 0_zz_xy_xx_0 0_zz_xx_zz_0 0_zz_xx_yz_0 0_zz_xx_yy_0 0_zz_xx_xz_0 0_zz_xx_xy_0 0_zz_xx_xx_0 0_zzz_zz_zz_0 0_zzz_zz_yz_0 0_zzz_zz_yy_0 0_zzz_zz_xz_0 0_zzz_zz_xy_0 0_zzz_zz_xx_0 0_zzz_yz_zz_0 0_zzz_yz_yz_0 0_zzz_yz_yy_0 0_zzz_yz_xz_0 0_zzz_yz_xy_0 0_zzz_yz_xx_0 0_zzz_yy_zz_0 0_zzz_yy_yz_0 0_zzz_yy_yy_0 0_zzz_yy_xz_0 0_zzz_yy_xy_0 0_zzz_yy_xx_0 0_zzz_xz_zz_0 0_zzz_xz_yz_0 0_zzz_xz_yy_0 0_zzz_xz_xz_0 0_zzz_xz_xy_0 0_zzz_xz_xx_0 0_zzz_xy_zz_0 0_zzz_xy_yz_0 0_zzz_xy_yy_0 0_zzz_xy_xz_0 0_zzz_xy_xy_0 0_zzz_xy_xx_0 0_zzz_xx_zz_0 0_zzz_xx_yz_0 0_zzz_xx_yy_0 0_zzz_xx_xz_0 0_zzz_xx_xy_0 0_zzz_xx_xx_0 0_yz_zz_zz_0 0_yz_zz_yz_0 0_yz_zz_yy_0 0_yz_zz_xz_0 0_yz_zz_xy_0 0_yz_zz_xx_0 0_yz_yz_zz_0 0_yz_yz_yz_0 0_yz_yz_yy_0 0_yz_yz_xz_0 0_yz_yz_xy_0 0_yz_yz_xx_0 0_yz_yy_zz_0 0_yz_yy_yz_0 0_yz_yy_yy_0 0_yz_yy_xz_0 0_yz_yy_xy_0 0_yz_yy_xx_0 0_yz_xz_zz_0 0_yz_xz_yz_0 0_yz_xz_yy_0 0_yz_xz_xz_0 0_yz_xz_xy_0 0_yz_xz_xx_0 0_yz_xy_zz_0 0_yz_xy_yz_0 0_yz_xy_yy_0 0_yz_xy_xz_0 0_yz_xy_xy_0 0_yz_xy_xx_0 0_yz_xx_zz_0 0_yz_xx_yz_0 0_yz_xx_yy_0 0_yz_xx_xz_0 0_yz_xx_xy_0 0_yz_xx_xx_0 0_yzz_zz_zz_0 0_yzz_zz_yz_0 0_yzz_zz_yy_0 0_yzz_zz_xz_0 0_yzz_zz_xy_0 0_yzz_zz_xx_0 0_yzz_yz_zz_0 0_yzz_yz_yz_0 0_yzz_yz_yy_0 0_yzz_yz_xz_0 0_yzz_yz_xy_0 0_yzz_yz_xx_0 0_yzz_yy_zz_0 0_yzz_yy_yz_0 0_yzz_yy_yy_0 0_yzz_yy_xz_0 0_yzz_yy_xy_0 0_yzz_yy_xx_0 0_yzz_xz_zz_0 0_yzz_xz_yz_0 0_yzz_xz_yy_0 0_yzz_xz_xz_0 0_yzz_xz_xy_0 0_yzz_xz_xx_0 0_yzz_xy_zz_0 0_yzz_xy_yz_0 0_yzz_xy_yy_0 0_yzz_xy_xz_0 0_yzz_xy_xy_0 0_yzz_xy_xx_0 0_yzz_xx_zz_0 0_yzz_xx_yz_0 0_yzz_xx_yy_0 0_yzz_xx_xz_0 0_yzz_xx_xy_0 0_yzz_xx_xx_0 0_yy_zz_zz_0 0_yy_zz_yz_0 0_yy_zz_yy_0 0_yy_zz_xz_0 0_yy_zz_xy_0 0_yy_zz_xx_0 0_yy_yz_zz_0 0_yy_yz_yz_0 0_yy_yz_yy_0 0_yy_yz_xz_0 0_yy_yz_xy_0 0_yy_yz_xx_0 0_yy_yy_zz_0 0_yy_yy_yz_0 0_yy_yy_yy_0 0_yy_yy_xz_0 0_yy_yy_xy_0 0_yy_yy_xx_0 0_yy_xz_zz_0 0_yy_xz_yz_0 0_yy_xz_yy_0 0_yy_xz_xz_0 0_yy_xz_xy_0 0_yy_xz_xx_0 0_yy_xy_zz_0 0_yy_xy_yz_0 0_yy_xy_yy_0 0_yy_xy_xz_0 0_yy_xy_xy_0 0_yy_xy_xx_0 0_yy_xx_zz_0 0_yy_xx_yz_0 0_yy_xx_yy_0 0_yy_xx_xz_0 0_yy_xx_xy_0 0_yy_xx_xx_0 0_yyz_zz_zz_0 0_yyz_zz_yz_0 0_yyz_zz_yy_0 0_yyz_zz_xz_0 0_yyz_zz_xy_0 0_yyz_zz_xx_0 0_yyz_yz_zz_0 0_yyz_yz_yz_0 0_yyz_yz_yy_0 0_yyz_yz_xz_0 0_yyz_yz_xy_0 0_yyz_yz_xx_0 0_yyz_yy_zz_0 0_yyz_yy_yz_0 0_yyz_yy_yy_0 0_yyz_yy_xz_0 0_yyz_yy_xy_0 0_yyz_yy_xx_0 0_yyz_xz_zz_0 0_yyz_xz_yz_0 0_yyz_xz_yy_0 0_yyz_xz_xz_0 0_yyz_xz_xy_0 0_yyz_xz_xx_0 0_yyz_xy_zz_0 0_yyz_xy_yz_0 0_yyz_xy_yy_0 0_yyz_xy_xz_0 0_yyz_xy_xy_0 0_yyz_xy_xx_0 0_yyz_xx_zz_0 0_yyz_xx_yz_0 0_yyz_xx_yy_0 0_yyz_xx_xz_0 0_yyz_xx_xy_0 0_yyz_xx_xx_0 0_yyy_zz_zz_0 0_yyy_zz_yz_0 0_yyy_zz_yy_0 0_yyy_zz_xz_0 0_yyy_zz_xy_0 0_yyy_zz_xx_0 0_yyy_yz_zz_0 0_yyy_yz_yz_0 0_yyy_yz_yy_0 0_yyy_yz_xz_0 0_yyy_yz_xy_0 0_yyy_yz_xx_0 0_yyy_yy_zz_0 0_yyy_yy_yz_0 0_yyy_yy_yy_0 0_yyy_yy_xz_0 0_yyy_yy_xy_0 0_yyy_yy_xx_0 0_yyy_xz_zz_0 0_yyy_xz_yz_0 0_yyy_xz_yy_0 0_yyy_xz_xz_0 0_yyy_xz_xy_0 0_yyy_xz_xx_0 0_yyy_xy_zz_0 0_yyy_xy_yz_0 0_yyy_xy_yy_0 0_yyy_xy_xz_0 0_yyy_xy_xy_0 0_yyy_xy_xx_0 0_yyy_xx_zz_0 0_yyy_xx_yz_0 0_yyy_xx_yy_0 0_yyy_xx_xz_0 0_yyy_xx_xy_0 0_yyy_xx_xx_0 0_xz_zz_zz_0 0_xz_zz_yz_0 0_xz_zz_yy_0 0_xz_zz_xz_0 0_xz_zz_xy_0 0_xz_zz_xx_0 0_xz_yz_zz_0 0_xz_yz_yz_0 0_xz_yz_yy_0 0_xz_yz_xz_0 0_xz_yz_xy_0 0_xz_yz_xx_0 0_xz_yy_zz_0 0_xz_yy_yz_0 0_xz_yy_yy_0 0_xz_yy_xz_0 0_xz_yy_xy_0 0_xz_yy_xx_0 0_xz_xz_zz_0 0_xz_xz_yz_0 0_xz_xz_yy_0 0_xz_xz_xz_0 0_xz_xz_xy_0 0_xz_xz_xx_0 0_xz_xy_zz_0 0_xz_xy_yz_0 0_xz_xy_yy_0 0_xz_xy_xz_0 0_xz_xy_xy_0 0_xz_xy_xx_0 0_xz_xx_zz_0 0_xz_xx_yz_0 0_xz_xx_yy_0 0_xz_xx_xz_0 0_xz_xx_xy_0 0_xz_xx_xx_0 0_xzz_zz_zz_0 0_xzz_zz_yz_0 0_xzz_zz_yy_0 0_xzz_zz_xz_0 0_xzz_zz_xy_0 0_xzz_zz_xx_0 0_xzz_yz_zz_0 0_xzz_yz_yz_0 0_xzz_yz_yy_0 0_xzz_yz_xz_0 0_xzz_yz_xy_0 0_xzz_yz_xx_0 0_xzz_yy_zz_0 0_xzz_yy_yz_0 0_xzz_yy_yy_0 0_xzz_yy_xz_0 0_xzz_yy_xy_0 0_xzz_yy_xx_0 0_xzz_xz_zz_0 0_xzz_xz_yz_0 0_xzz_xz_yy_0 0_xzz_xz_xz_0 0_xzz_xz_xy_0 0_xzz_xz_xx_0 0_xzz_xy_zz_0 0_xzz_xy_yz_0 0_xzz_xy_yy_0 0_xzz_xy_xz_0 0_xzz_xy_xy_0 0_xzz_xy_xx_0 0_xzz_xx_zz_0 0_xzz_xx_yz_0 0_xzz_xx_yy_0 0_xzz_xx_xz_0 0_xzz_xx_xy_0 0_xzz_xx_xx_0 0_xy_zz_zz_0 0_xy_zz_yz_0 0_xy_zz_yy_0 0_xy_zz_xz_0 0_xy_zz_xy_0 0_xy_zz_xx_0 0_xy_yz_zz_0 0_xy_yz_yz_0 0_xy_yz_yy_0 0_xy_yz_xz_0 0_xy_yz_xy_0 0_xy_yz_xx_0 0_xy_yy_zz_0 0_xy_yy_yz_0 0_xy_yy_yy_0 0_xy_yy_xz_0 0_xy_yy_xy_0 0_xy_yy_xx_0 0_xy_xz_zz_0 0_xy_xz_yz_0 0_xy_xz_yy_0 0_xy_xz_xz_0 0_xy_xz_xy_0 0_xy_xz_xx_0 0_xy_xy_zz_0 0_xy_xy_yz_0 0_xy_xy_yy_0 0_xy_xy_xz_0 0_xy_xy_xy_0 0_xy_xy_xx_0 0_xy_xx_zz_0 0_xy_xx_yz_0 0_xy_xx_yy_0 0_xy_xx_xz_0 0_xy_xx_xy_0 0_xy_xx_xx_0 0_xyz_zz_zz_0 0_xyz_zz_yz_0 0_xyz_zz_yy_0 0_xyz_zz_xz_0 0_xyz_zz_xy_0 0_xyz_zz_xx_0 0_xyz_yz_zz_0 0_xyz_yz_yz_0 0_xyz_yz_yy_0 0_xyz_yz_xz_0 0_xyz_yz_xy_0 0_xyz_yz_xx_0 0_xyz_yy_zz_0 0_xyz_yy_yz_0 0_xyz_yy_yy_0 0_xyz_yy_xz_0 0_xyz_yy_xy_0 0_xyz_yy_xx_0 0_xyz_xz_zz_0 0_xyz_xz_yz_0 0_xyz_xz_yy_0 0_xyz_xz_xz_0 0_xyz_xz_xy_0 0_xyz_xz_xx_0 0_xyz_xy_zz_0 0_xyz_xy_yz_0 0_xyz_xy_yy_0 0_xyz_xy_xz_0 0_xyz_xy_xy_0 0_xyz_xy_xx_0 0_xyz_xx_zz_0 0_xyz_xx_yz_0 0_xyz_xx_yy_0 0_xyz_xx_xz_0 0_xyz_xx_xy_0 0_xyz_xx_xx_0 0_xyy_zz_zz_0 0_xyy_zz_yz_0 0_xyy_zz_yy_0 0_xyy_zz_xz_0 0_xyy_zz_xy_0 0_xyy_zz_xx_0 0_xyy_yz_zz_0 0_xyy_yz_yz_0 0_xyy_yz_yy_0 0_xyy_yz_xz_0 0_xyy_yz_xy_0 0_xyy_yz_xx_0 0_xyy_yy_zz_0 0_xyy_yy_yz_0 0_xyy_yy_yy_0 0_xyy_yy_xz_0 0_xyy_yy_xy_0 0_xyy_yy_xx_0 0_xyy_xz_zz_0 0_xyy_xz_yz_0 0_xyy_xz_yy_0 0_xyy_xz_xz_0 0_xyy_xz_xy_0 0_xyy_xz_xx_0 0_xyy_xy_zz_0 0_xyy_xy_yz_0 0_xyy_xy_yy_0 0_xyy_xy_xz_0 0_xyy_xy_xy_0 0_xyy_xy_xx_0 0_xyy_xx_zz_0 0_xyy_xx_yz_0 0_xyy_xx_yy_0 0_xyy_xx_xz_0 0_xyy_xx_xy_0 0_xyy_xx_xx_0 0_xx_zz_zz_0 0_xx_zz_yz_0 0_xx_zz_yy_0 0_xx_zz_xz_0 0_xx_zz_xy_0 0_xx_zz_xx_0 0_xx_yz_zz_0 0_xx_yz_yz_0 0_xx_yz_yy_0 0_xx_yz_xz_0 0_xx_yz_xy_0 0_xx_yz_xx_0 0_xx_yy_zz_0 0_xx_yy_yz_0 0_xx_yy_yy_0 0_xx_yy_xz_0 0_xx_yy_xy_0 0_xx_yy_xx_0 0_xx_xz_zz_0 0_xx_xz_yz_0 0_xx_xz_yy_0 0_xx_xz_xz_0 0_xx_xz_xy_0 0_xx_xz_xx_0 0_xx_xy_zz_0 0_xx_xy_yz_0 0_xx_xy_yy_0 0_xx_xy_xz_0 0_xx_xy_xy_0 0_xx_xy_xx_0 0_xx_xx_zz_0 0_xx_xx_yz_0 0_xx_xx_yy_0 0_xx_xx_xz_0 0_xx_xx_xy_0 0_xx_xx_xx_0 0_xxz_zz_zz_0 0_xxz_zz_yz_0 0_xxz_zz_yy_0 0_xxz_zz_xz_0 0_xxz_zz_xy_0 0_xxz_zz_xx_0 0_xxz_yz_zz_0 0_xxz_yz_yz_0 0_xxz_yz_yy_0 0_xxz_yz_xz_0 0_xxz_yz_xy_0 0_xxz_yz_xx_0 0_xxz_yy_zz_0 0_xxz_yy_yz_0 0_xxz_yy_yy_0 0_xxz_yy_xz_0 0_xxz_yy_xy_0 0_xxz_yy_xx_0 0_xxz_xz_zz_0 0_xxz_xz_yz_0 0_xxz_xz_yy_0 0_xxz_xz_xz_0 0_xxz_xz_xy_0 0_xxz_xz_xx_0 0_xxz_xy_zz_0 0_xxz_xy_yz_0 0_xxz_xy_yy_0 0_xxz_xy_xz_0 0_xxz_xy_xy_0 0_xxz_xy_xx_0 0_xxz_xx_zz_0 0_xxz_xx_yz_0 0_xxz_xx_yy_0 0_xxz_xx_xz_0 0_xxz_xx_xy_0 0_xxz_xx_xx_0 0_xxy_zz_zz_0 0_xxy_zz_yz_0 0_xxy_zz_yy_0 0_xxy_zz_xz_0 0_xxy_zz_xy_0 0_xxy_zz_xx_0 0_xxy_yz_zz_0 0_xxy_yz_yz_0 0_xxy_yz_yy_0 0_xxy_yz_xz_0 0_xxy_yz_xy_0 0_xxy_yz_xx_0 0_xxy_yy_zz_0 0_xxy_yy_yz_0 0_xxy_yy_yy_0 0_xxy_yy_xz_0 0_xxy_yy_xy_0 0_xxy_yy_xx_0 0_xxy_xz_zz_0 0_xxy_xz_yz_0 0_xxy_xz_yy_0 0_xxy_xz_xz_0 0_xxy_xz_xy_0 0_xxy_xz_xx_0 0_xxy_xy_zz_0 0_xxy_xy_yz_0 0_xxy_xy_yy_0 0_xxy_xy_xz_0 0_xxy_xy_xy_0 0_xxy_xy_xx_0 0_xxy_xx_zz_0 0_xxy_xx_yz_0 0_xxy_xx_yy_0 0_xxy_xx_xz_0 0_xxy_xx_xy_0 0_xxy_xx_xx_0 0_xxx_zz_zz_0 0_xxx_zz_yz_0 0_xxx_zz_yy_0 0_xxx_zz_xz_0 0_xxx_zz_xy_0 0_xxx_zz_xx_0 0_xxx_yz_zz_0 0_xxx_yz_yz_0 0_xxx_yz_yy_0 0_xxx_yz_xz_0 0_xxx_yz_xy_0 0_xxx_yz_xx_0 0_xxx_yy_zz_0 0_xxx_yy_yz_0 0_xxx_yy_yy_0 0_xxx_yy_xz_0 0_xxx_yy_xy_0 0_xxx_yy_xx_0 0_xxx_xz_zz_0 0_xxx_xz_yz_0 0_xxx_xz_yy_0 0_xxx_xz_xz_0 0_xxx_xz_xy_0 0_xxx_xz_xx_0 0_xxx_xy_zz_0 0_xxx_xy_yz_0 0_xxx_xy_yy_0 0_xxx_xy_xz_0 0_xxx_xy_xy_0 0_xxx_xy_xx_0 0_xxx_xx_zz_0 0_xxx_xx_yz_0 0_xxx_xx_yy_0 0_xxx_xx_xz_0 0_xxx_xx_xy_0 0_xxx_xx_xx_0 z_zz_zz_zz_0 z_zz_zz_yz_0 z_zz_zz_yy_0 z_zz_zz_xz_0 z_zz_zz_xy_0 z_zz_zz_xx_0 z_zz_yz_zz_0 z_zz_yz_yz_0 z_zz_yz_yy_0 z_zz_yz_xz_0 z_zz_yz_xy_0 z_zz_yz_xx_0 z_zz_yy_zz_0 z_zz_yy_yz_0 z_zz_yy_yy_0 z_zz_yy_xz_0 z_zz_yy_xy_0 z_zz_yy_xx_0 z_zz_xz_zz_0 z_zz_xz_yz_0 z_zz_xz_yy_0 z_zz_xz_xz_0 z_zz_xz_xy_0 z_zz_xz_xx_0 z_zz_xy_zz_0 z_zz_xy_yz_0 z_zz_xy_yy_0 z_zz_xy_xz_0 z_zz_xy_xy_0 z_zz_xy_xx_0 z_zz_xx_zz_0 z_zz_xx_yz_0 z_zz_xx_yy_0 z_zz_xx_xz_0 z_zz_xx_xy_0 z_zz_xx_xx_0 z_yz_zz_zz_0 z_yz_zz_yz_0 z_yz_zz_yy_0 z_yz_zz_xz_0 z_yz_zz_xy_0 z_yz_zz_xx_0 z_yz_yz_zz_0 z_yz_yz_yz_0 z_yz_yz_yy_0 z_yz_yz_xz_0 z_yz_yz_xy_0 z_yz_yz_xx_0 z_yz_yy_zz_0 z_yz_yy_yz_0 z_yz_yy_yy_0 z_yz_yy_xz_0 z_yz_yy_xy_0 z_yz_yy_xx_0 z_yz_xz_zz_0 z_yz_xz_yz_0 z_yz_xz_yy_0 z_yz_xz_xz_0 z_yz_xz_xy_0 z_yz_xz_xx_0 z_yz_xy_zz_0 z_yz_xy_yz_0 z_yz_xy_yy_0 z_yz_xy_xz_0 z_yz_xy_xy_0 z_yz_xy_xx_0 z_yz_xx_zz_0 z_yz_xx_yz_0 z_yz_xx_yy_0 z_yz_xx_xz_0 z_yz_xx_xy_0 z_yz_xx_xx_0 z_yy_zz_zz_0 z_yy_zz_yz_0 z_yy_zz_yy_0 z_yy_zz_xz_0 z_yy_zz_xy_0 z_yy_zz_xx_0 z_yy_yz_zz_0 z_yy_yz_yz_0 z_yy_yz_yy_0 z_yy_yz_xz_0 z_yy_yz_xy_0 z_yy_yz_xx_0 z_yy_yy_zz_0 z_yy_yy_yz_0 z_yy_yy_yy_0 z_yy_yy_xz_0 z_yy_yy_xy_0 z_yy_yy_xx_0 z_yy_xz_zz_0 z_yy_xz_yz_0 z_yy_xz_yy_0 z_yy_xz_xz_0 z_yy_xz_xy_0 z_yy_xz_xx_0 z_yy_xy_zz_0 z_yy_xy_yz_0 z_yy_xy_yy_0 z_yy_xy_xz_0 z_yy_xy_xy_0 z_yy_xy_xx_0 z_yy_xx_zz_0 z_yy_xx_yz_0 z_yy_xx_yy_0 z_yy_xx_xz_0 z_yy_xx_xy_0 z_yy_xx_xx_0 z_xz_zz_zz_0 z_xz_zz_yz_0 z_xz_zz_yy_0 z_xz_zz_xz_0 z_xz_zz_xy_0 z_xz_zz_xx_0 z_xz_yz_zz_0 z_xz_yz_yz_0 z_xz_yz_yy_0 z_xz_yz_xz_0 z_xz_yz_xy_0 z_xz_yz_xx_0 z_xz_yy_zz_0 z_xz_yy_yz_0 z_xz_yy_yy_0 z_xz_yy_xz_0 z_xz_yy_xy_0 z_xz_yy_xx_0 z_xz_xz_zz_0 z_xz_xz_yz_0 z_xz_xz_yy_0 z_xz_xz_xz_0 z_xz_xz_xy_0 z_xz_xz_xx_0 z_xz_xy_zz_0 z_xz_xy_yz_0 z_xz_xy_yy_0 z_xz_xy_xz_0 z_xz_xy_xy_0 z_xz_xy_xx_0 z_xz_xx_zz_0 z_xz_xx_yz_0 z_xz_xx_yy_0 z_xz_xx_xz_0 z_xz_xx_xy_0 z_xz_xx_xx_0 z_xy_zz_zz_0 z_xy_zz_yz_0 z_xy_zz_yy_0 z_xy_zz_xz_0 z_xy_zz_xy_0 z_xy_zz_xx_0 z_xy_yz_zz_0 z_xy_yz_yz_0 z_xy_yz_yy_0 z_xy_yz_xz_0 z_xy_yz_xy_0 z_xy_yz_xx_0 z_xy_yy_zz_0 z_xy_yy_yz_0 z_xy_yy_yy_0 z_xy_yy_xz_0 z_xy_yy_xy_0 z_xy_yy_xx_0 z_xy_xz_zz_0 z_xy_xz_yz_0 z_xy_xz_yy_0 z_xy_xz_xz_0 z_xy_xz_xy_0 z_xy_xz_xx_0 z_xy_xy_zz_0 z_xy_xy_yz_0 z_xy_xy_yy_0 z_xy_xy_xz_0 z_xy_xy_xy_0 z_xy_xy_xx_0 z_xy_xx_zz_0 z_xy_xx_yz_0 z_xy_xx_yy_0 z_xy_xx_xz_0 z_xy_xx_xy_0 z_xy_xx_xx_0 z_xx_zz_zz_0 z_xx_zz_yz_0 z_xx_zz_yy_0 z_xx_zz_xz_0 z_xx_zz_xy_0 z_xx_zz_xx_0 z_xx_yz_zz_0 z_xx_yz_yz_0 z_xx_yz_yy_0 z_xx_yz_xz_0 z_xx_yz_xy_0 z_xx_yz_xx_0 z_xx_yy_zz_0 z_xx_yy_yz_0 z_xx_yy_yy_0 z_xx_yy_xz_0 z_xx_yy_xy_0 z_xx_yy_xx_0 z_xx_xz_zz_0 z_xx_xz_yz_0 z_xx_xz_yy_0 z_xx_xz_xz_0 z_xx_xz_xy_0 z_xx_xz_xx_0 z_xx_xy_zz_0 z_xx_xy_yz_0 z_xx_xy_yy_0 z_xx_xy_xz_0 z_xx_xy_xy_0 z_xx_xy_xx_0 z_xx_xx_zz_0 z_xx_xx_yz_0 z_xx_xx_yy_0 z_xx_xx_xz_0 z_xx_xx_xy_0 z_xx_xx_xx_0 y_zz_zz_zz_0 y_zz_zz_yz_0 y_zz_zz_yy_0 y_zz_zz_xz_0 y_zz_zz_xy_0 y_zz_zz_xx_0 y_zz_yz_zz_0 y_zz_yz_yz_0 y_zz_yz_yy_0 y_zz_yz_xz_0 y_zz_yz_xy_0 y_zz_yz_xx_0 y_zz_yy_zz_0 y_zz_yy_yz_0 y_zz_yy_yy_0 y_zz_yy_xz_0 y_zz_yy_xy_0 y_zz_yy_xx_0 y_zz_xz_zz_0 y_zz_xz_yz_0 y_zz_xz_yy_0 y_zz_xz_xz_0 y_zz_xz_xy_0 y_zz_xz_xx_0 y_zz_xy_zz_0 y_zz_xy_yz_0 y_zz_xy_yy_0 y_zz_xy_xz_0 y_zz_xy_xy_0 y_zz_xy_xx_0 y_zz_xx_zz_0 y_zz_xx_yz_0 y_zz_xx_yy_0 y_zz_xx_xz_0 y_zz_xx_xy_0 y_zz_xx_xx_0 y_yz_zz_zz_0 y_yz_zz_yz_0 y_yz_zz_yy_0 y_yz_zz_xz_0 y_yz_zz_xy_0 y_yz_zz_xx_0 y_yz_yz_zz_0 y_yz_yz_yz_0 y_yz_yz_yy_0 y_yz_yz_xz_0 y_yz_yz_xy_0 y_yz_yz_xx_0 y_yz_yy_zz_0 y_yz_yy_yz_0 y_yz_yy_yy_0 y_yz_yy_xz_0 y_yz_yy_xy_0 y_yz_yy_xx_0 y_yz_xz_zz_0 y_yz_xz_yz_0 y_yz_xz_yy_0 y_yz_xz_xz_0 y_yz_xz_xy_0 y_yz_xz_xx_0 y_yz_xy_zz_0 y_yz_xy_yz_0 y_yz_xy_yy_0 y_yz_xy_xz_0 y_yz_xy_xy_0 y_yz_xy_xx_0 y_yz_xx_zz_0 y_yz_xx_yz_0 y_yz_xx_yy_0 y_yz_xx_xz_0 y_yz_xx_xy_0 y_yz_xx_xx_0 y_yy_zz_zz_0 y_yy_zz_yz_0 y_yy_zz_yy_0 y_yy_zz_xz_0 y_yy_zz_xy_0 y_yy_zz_xx_0 y_yy_yz_zz_0 y_yy_yz_yz_0 y_yy_yz_yy_0 y_yy_yz_xz_0 y_yy_yz_xy_0 y_yy_yz_xx_0 y_yy_yy_zz_0 y_yy_yy_yz_0 y_yy_yy_yy_0 y_yy_yy_xz_0 y_yy_yy_xy_0 y_yy_yy_xx_0 y_yy_xz_zz_0 y_yy_xz_yz_0 y_yy_xz_yy_0 y_yy_xz_xz_0 y_yy_xz_xy_0 y_yy_xz_xx_0 y_yy_xy_zz_0 y_yy_xy_yz_0 y_yy_xy_yy_0 y_yy_xy_xz_0 y_yy_xy_xy_0 y_yy_xy_xx_0 y_yy_xx_zz_0 y_yy_xx_yz_0 y_yy_xx_yy_0 y_yy_xx_xz_0 y_yy_xx_xy_0 y_yy_xx_xx_0 y_xz_zz_zz_0 y_xz_zz_yz_0 y_xz_zz_yy_0 y_xz_zz_xz_0 y_xz_zz_xy_0 y_xz_zz_xx_0 y_xz_yz_zz_0 y_xz_yz_yz_0 y_xz_yz_yy_0 y_xz_yz_xz_0 y_xz_yz_xy_0 y_xz_yz_xx_0 y_xz_yy_zz_0 y_xz_yy_yz_0 y_xz_yy_yy_0 y_xz_yy_xz_0 y_xz_yy_xy_0 y_xz_yy_xx_0 y_xz_xz_zz_0 y_xz_xz_yz_0 y_xz_xz_yy_0 y_xz_xz_xz_0 y_xz_xz_xy_0 y_xz_xz_xx_0 y_xz_xy_zz_0 y_xz_xy_yz_0 y_xz_xy_yy_0 y_xz_xy_xz_0 y_xz_xy_xy_0 y_xz_xy_xx_0 y_xz_xx_zz_0 y_xz_xx_yz_0 y_xz_xx_yy_0 y_xz_xx_xz_0 y_xz_xx_xy_0 y_xz_xx_xx_0 y_xy_zz_zz_0 y_xy_zz_yz_0 y_xy_zz_yy_0 y_xy_zz_xz_0 y_xy_zz_xy_0 y_xy_zz_xx_0 y_xy_yz_zz_0 y_xy_yz_yz_0 y_xy_yz_yy_0 y_xy_yz_xz_0 y_xy_yz_xy_0 y_xy_yz_xx_0 y_xy_yy_zz_0 y_xy_yy_yz_0 y_xy_yy_yy_0 y_xy_yy_xz_0 y_xy_yy_xy_0 y_xy_yy_xx_0 y_xy_xz_zz_0 y_xy_xz_yz_0 y_xy_xz_yy_0 y_xy_xz_xz_0 y_xy_xz_xy_0 y_xy_xz_xx_0 y_xy_xy_zz_0 y_xy_xy_yz_0 y_xy_xy_yy_0 y_xy_xy_xz_0 y_xy_xy_xy_0 y_xy_xy_xx_0 y_xy_xx_zz_0 y_xy_xx_yz_0 y_xy_xx_yy_0 y_xy_xx_xz_0 y_xy_xx_xy_0 y_xy_xx_xx_0 y_xx_zz_zz_0 y_xx_zz_yz_0 y_xx_zz_yy_0 y_xx_zz_xz_0 y_xx_zz_xy_0 y_xx_zz_xx_0 y_xx_yz_zz_0 y_xx_yz_yz_0 y_xx_yz_yy_0 y_xx_yz_xz_0 y_xx_yz_xy_0 y_xx_yz_xx_0 y_xx_yy_zz_0 y_xx_yy_yz_0 y_xx_yy_yy_0 y_xx_yy_xz_0 y_xx_yy_xy_0 y_xx_yy_xx_0 y_xx_xz_zz_0 y_xx_xz_yz_0 y_xx_xz_yy_0 y_xx_xz_xz_0 y_xx_xz_xy_0 y_xx_xz_xx_0 y_xx_xy_zz_0 y_xx_xy_yz_0 y_xx_xy_yy_0 y_xx_xy_xz_0 y_xx_xy_xy_0 y_xx_xy_xx_0 y_xx_xx_zz_0 y_xx_xx_yz_0 y_xx_xx_yy_0 y_xx_xx_xz_0 y_xx_xx_xy_0 y_xx_xx_xx_0 x_zz_zz_zz_0 x_zz_zz_yz_0 x_zz_zz_yy_0 x_zz_zz_xz_0 x_zz_zz_xy_0 x_zz_zz_xx_0 x_zz_yz_zz_0 x_zz_yz_yz_0 x_zz_yz_yy_0 x_zz_yz_xz_0 x_zz_yz_xy_0 x_zz_yz_xx_0 x_zz_yy_zz_0 x_zz_yy_yz_0 x_zz_yy_yy_0 x_zz_yy_xz_0 x_zz_yy_xy_0 x_zz_yy_xx_0 x_zz_xz_zz_0 x_zz_xz_yz_0 x_zz_xz_yy_0 x_zz_xz_xz_0 x_zz_xz_xy_0 x_zz_xz_xx_0 x_zz_xy_zz_0 x_zz_xy_yz_0 x_zz_xy_yy_0 x_zz_xy_xz_0 x_zz_xy_xy_0 x_zz_xy_xx_0 x_zz_xx_zz_0 x_zz_xx_yz_0 x_zz_xx_yy_0 x_zz_xx_xz_0 x_zz_xx_xy_0 x_zz_xx_xx_0 x_yz_zz_zz_0 x_yz_zz_yz_0 x_yz_zz_yy_0 x_yz_zz_xz_0 x_yz_zz_xy_0 x_yz_zz_xx_0 x_yz_yz_zz_0 x_yz_yz_yz_0 x_yz_yz_yy_0 x_yz_yz_xz_0 x_yz_yz_xy_0 x_yz_yz_xx_0 x_yz_yy_zz_0 x_yz_yy_yz_0 x_yz_yy_yy_0 x_yz_yy_xz_0 x_yz_yy_xy_0 x_yz_yy_xx_0 x_yz_xz_zz_0 x_yz_xz_yz_0 x_yz_xz_yy_0 x_yz_xz_xz_0 x_yz_xz_xy_0 x_yz_xz_xx_0 x_yz_xy_zz_0 x_yz_xy_yz_0 x_yz_xy_yy_0 x_yz_xy_xz_0 x_yz_xy_xy_0 x_yz_xy_xx_0 x_yz_xx_zz_0 x_yz_xx_yz_0 x_yz_xx_yy_0 x_yz_xx_xz_0 x_yz_xx_xy_0 x_yz_xx_xx_0 x_yy_zz_zz_0 x_yy_zz_yz_0 x_yy_zz_yy_0 x_yy_zz_xz_0 x_yy_zz_xy_0 x_yy_zz_xx_0 x_yy_yz_zz_0 x_yy_yz_yz_0 x_yy_yz_yy_0 x_yy_yz_xz_0 x_yy_yz_xy_0 x_yy_yz_xx_0 x_yy_yy_zz_0 x_yy_yy_yz_0 x_yy_yy_yy_0 x_yy_yy_xz_0 x_yy_yy_xy_0 x_yy_yy_xx_0 x_yy_xz_zz_0 x_yy_xz_yz_0 x_yy_xz_yy_0 x_yy_xz_xz_0 x_yy_xz_xy_0 x_yy_xz_xx_0 x_yy_xy_zz_0 x_yy_xy_yz_0 x_yy_xy_yy_0 x_yy_xy_xz_0 x_yy_xy_xy_0 x_yy_xy_xx_0 x_yy_xx_zz_0 x_yy_xx_yz_0 x_yy_xx_yy_0 x_yy_xx_xz_0 x_yy_xx_xy_0 x_yy_xx_xx_0 x_xz_zz_zz_0 x_xz_zz_yz_0 x_xz_zz_yy_0 x_xz_zz_xz_0 x_xz_zz_xy_0 x_xz_zz_xx_0 x_xz_yz_zz_0 x_xz_yz_yz_0 x_xz_yz_yy_0 x_xz_yz_xz_0 x_xz_yz_xy_0 x_xz_yz_xx_0 x_xz_yy_zz_0 x_xz_yy_yz_0 x_xz_yy_yy_0 x_xz_yy_xz_0 x_xz_yy_xy_0 x_xz_yy_xx_0 x_xz_xz_zz_0 x_xz_xz_yz_0 x_xz_xz_yy_0 x_xz_xz_xz_0 x_xz_xz_xy_0 x_xz_xz_xx_0 x_xz_xy_zz_0 x_xz_xy_yz_0 x_xz_xy_yy_0 x_xz_xy_xz_0 x_xz_xy_xy_0 x_xz_xy_xx_0 x_xz_xx_zz_0 x_xz_xx_yz_0 x_xz_xx_yy_0 x_xz_xx_xz_0 x_xz_xx_xy_0 x_xz_xx_xx_0 x_xy_zz_zz_0 x_xy_zz_yz_0 x_xy_zz_yy_0 x_xy_zz_xz_0 x_xy_zz_xy_0 x_xy_zz_xx_0 x_xy_yz_zz_0 x_xy_yz_yz_0 x_xy_yz_yy_0 x_xy_yz_xz_0 x_xy_yz_xy_0 x_xy_yz_xx_0 x_xy_yy_zz_0 x_xy_yy_yz_0 x_xy_yy_yy_0 x_xy_yy_xz_0 x_xy_yy_xy_0 x_xy_yy_xx_0 x_xy_xz_zz_0 x_xy_xz_yz_0 x_xy_xz_yy_0 x_xy_xz_xz_0 x_xy_xz_xy_0 x_xy_xz_xx_0 x_xy_xy_zz_0 x_xy_xy_yz_0 x_xy_xy_yy_0 x_xy_xy_xz_0 x_xy_xy_xy_0 x_xy_xy_xx_0 x_xy_xx_zz_0 x_xy_xx_yz_0 x_xy_xx_yy_0 x_xy_xx_xz_0 x_xy_xx_xy_0 x_xy_xx_xx_0 x_xx_zz_zz_0 x_xx_zz_yz_0 x_xx_zz_yy_0 x_xx_zz_xz_0 x_xx_zz_xy_0 x_xx_zz_xx_0 x_xx_yz_zz_0 x_xx_yz_yz_0 x_xx_yz_yy_0 x_xx_yz_xz_0 x_xx_yz_xy_0 x_xx_yz_xx_0 x_xx_yy_zz_0 x_xx_yy_yz_0 x_xx_yy_yy_0 x_xx_yy_xz_0 x_xx_yy_xy_0 x_xx_yy_xx_0 x_xx_xz_zz_0 x_xx_xz_yz_0 x_xx_xz_yy_0 x_xx_xz_xz_0 x_xx_xz_xy_0 x_xx_xz_xx_0 x_xx_xy_zz_0 x_xx_xy_yz_0 x_xx_xy_yy_0 x_xx_xy_xz_0 x_xx_xy_xy_0 x_xx_xy_xx_0 x_xx_xx_zz_0 x_xx_xx_yz_0 x_xx_xx_yy_0 x_xx_xx_xz_0 x_xx_xx_xy_0 x_xx_xx_xx_0 

signature:
0_zzz_z_zz_0 0_zzz_z_zzz_0 0_zzz_z_yz_0 0_zzz_z_yzz_0 0_zzz_z_yy_0 0_zzz_z_yyz_0 0_zzz_z_xz_0 0_zzz_z_xzz_0 0_zzz_z_xy_0 0_zzz_z_xyz_0 0_zzz_z_xx_0 0_zzz_z_xxz_0 0_zzz_zz_zz_0 0_zzz_zz_yz_0 0_zzz_zz_yy_0 0_zzz_zz_xz_0 0_zzz_zz_xy_0 0_zzz_zz_xx_0 0_zzz_y_zz_0 0_zzz_y_zzz_0 0_zzz_y_yz_0 0_zzz_y_yzz_0 0_zzz_y_yy_0 0_zzz_y_yyz_0 0_zzz_y_yyy_0 0_zzz_y_xz_0 0_zzz_y_xzz_0 0_zzz_y_xy_0 0_zzz_y_xyz_0 0_zzz_y_xyy_0 0_zzz_y_xx_0 0_zzz_y_xxz_0 0_zzz_y_xxy_0 0_zzz_yz_zz_0 0_zzz_yz_yz_0 0_zzz_yz_yy_0 0_zzz_yz_xz_0 0_zzz_yz_xy_0 0_zzz_yz_xx_0 0_zzz_yy_zz_0 0_zzz_yy_yz_0 0_zzz_yy_yy_0 0_zzz_yy_xz_0 0_zzz_yy_xy_0 0_zzz_yy_xx_0 0_zzz_x_zz_0 0_zzz_x_zzz_0 0_zzz_x_yz_0 0_zzz_x_yzz_0 0_zzz_x_yy_0 0_zzz_x_yyz_0 0_zzz_x_yyy_0 0_zzz_x_xz_0 0_zzz_x_xzz_0 0_zzz_x_xy_0 0_zzz_x_xyz_0 0_zzz_x_xyy_0 0_zzz_x_xx_0 0_zzz_x_xxz_0 0_zzz_x_xxy_0 0_zzz_x_xxx_0 0_zzz_xz_zz_0 0_zzz_xz_yz_0 0_zzz_xz_yy_0 0_zzz_xz_xz_0 0_zzz_xz_xy_0 0_zzz_xz_xx_0 0_zzz_xy_zz_0 0_zzz_xy_yz_0 0_zzz_xy_yy_0 0_zzz_xy_xz_0 0_zzz_xy_xy_0 0_zzz_xy_xx_0 0_zzz_xx_zz_0 0_zzz_xx_yz_0 0_zzz_xx_yy_0 0_zzz_xx_xz_0 0_zzz_xx_xy_0 0_zzz_xx_xx_0 0_yzz_z_zz_0 0_yzz_z_zzz_0 0_yzz_z_yz_0 0_yzz_z_yzz_0 0_yzz_z_yy_0 0_yzz_z_yyz_0 0_yzz_z_xz_0 0_yzz_z_xzz_0 0_yzz_z_xy_0 0_yzz_z_xyz_0 0_yzz_z_xx_0 0_yzz_z_xxz_0 0_yzz_zz_zz_0 0_yzz_zz_yz_0 0_yzz_zz_yy_0 0_yzz_zz_xz_0 0_yzz_zz_xy_0 0_yzz_zz_xx_0 0_yzz_y_zz_0 0_yzz_y_zzz_0 0_yzz_y_yz_0 0_yzz_y_yzz_0 0_yzz_y_yy_0 0_yzz_y_yyz_0 0_yzz_y_yyy_0 0_yzz_y_xz_0 0_yzz_y_xzz_0 0_yzz_y_xy_0 0_yzz_y_xyz_0 0_yzz_y_xyy_0 0_yzz_y_xx_0 0_yzz_y_xxz_0 0_yzz_y_xxy_0 0_yzz_yz_zz_0 0_yzz_yz_yz_0 0_yzz_yz_yy_0 0_yzz_yz_xz_0 0_yzz_yz_xy_0 0_yzz_yz_xx_0 0_yzz_yy_zz_0 0_yzz_yy_yz_0 0_yzz_yy_yy_0 0_yzz_yy_xz_0 0_yzz_yy_xy_0 0_yzz_yy_xx_0 0_yzz_x_zz_0 0_yzz_x_zzz_0 0_yzz_x_yz_0 0_yzz_x_yzz_0 0_yzz_x_yy_0 0_yzz_x_yyz_0 0_yzz_x_yyy_0 0_yzz_x_xz_0 0_yzz_x_xzz_0 0_yzz_x_xy_0 0_yzz_x_xyz_0 0_yzz_x_xyy_0 0_yzz_x_xx_0 0_yzz_x_xxz_0 0_yzz_x_xxy_0 0_yzz_x_xxx_0 0_yzz_xz_zz_0 0_yzz_xz_yz_0 0_yzz_xz_yy_0 0_yzz_xz_xz_0 0_yzz_xz_xy_0 0_yzz_xz_xx_0 0_yzz_xy_zz_0 0_yzz_xy_yz_0 0_yzz_xy_yy_0 0_yzz_xy_xz_0 0_yzz_xy_xy_0 0_yzz_xy_xx_0 0_yzz_xx_zz_0 0_yzz_xx_yz_0 0_yzz_xx_yy_0 0_yzz_xx_xz_0 0_yzz_xx_xy_0 0_yzz_xx_xx_0 0_yyz_z_zz_0 0_yyz_z_zzz_0 0_yyz_z_yz_0 0_yyz_z_yzz_0 0_yyz_z_yy_0 0_yyz_z_yyz_0 0_yyz_z_xz_0 0_yyz_z_xzz_0 0_yyz_z_xy_0 0_yyz_z_xyz_0 0_yyz_z_xx_0 0_yyz_z_xxz_0 0_yyz_zz_zz_0 0_yyz_zz_yz_0 0_yyz_zz_yy_0 0_yyz_zz_xz_0 0_yyz_zz_xy_0 0_yyz_zz_xx_0 0_yyz_y_zz_0 0_yyz_y_zzz_0 0_yyz_y_yz_0 0_yyz_y_yzz_0 0_yyz_y_yy_0 0_yyz_y_yyz_0 0_yyz_y_yyy_0 0_yyz_y_xz_0 0_yyz_y_xzz_0 0_yyz_y_xy_0 0_yyz_y_xyz_0 0_yyz_y_xyy_0 0_yyz_y_xx_0 0_yyz_y_xxz_0 0_yyz_y_xxy_0 0_yyz_yz_zz_0 0_yyz_yz_yz_0 0_yyz_yz_yy_0 0_yyz_yz_xz_0 0_yyz_yz_xy_0 0_yyz_yz_xx_0 0_yyz_yy_zz_0 0_yyz_yy_yz_0 0_yyz_yy_yy_0 0_yyz_yy_xz_0 0_yyz_yy_xy_0 0_yyz_yy_xx_0 0_yyz_x_zz_0 0_yyz_x_zzz_0 0_yyz_x_yz_0 0_yyz_x_yzz_0 0_yyz_x_yy_0 0_yyz_x_yyz_0 0_yyz_x_yyy_0 0_yyz_x_xz_0 0_yyz_x_xzz_0 0_yyz_x_xy_0 0_yyz_x_xyz_0 0_yyz_x_xyy_0 0_yyz_x_xx_0 0_yyz_x_xxz_0 0_yyz_x_xxy_0 0_yyz_x_xxx_0 0_yyz_xz_zz_0 0_yyz_xz_yz_0 0_yyz_xz_yy_0 0_yyz_xz_xz_0 0_yyz_xz_xy_0 0_yyz_xz_xx_0 0_yyz_xy_zz_0 0_yyz_xy_yz_0 0_yyz_xy_yy_0 0_yyz_xy_xz_0 0_yyz_xy_xy_0 0_yyz_xy_xx_0 0_yyz_xx_zz_0 0_yyz_xx_yz_0 0_yyz_xx_yy_0 0_yyz_xx_xz_0 0_yyz_xx_xy_0 0_yyz_xx_xx_0 0_yyy_z_zz_0 0_yyy_z_zzz_0 0_yyy_z_yz_0 0_yyy_z_yzz_0 0_yyy_z_yy_0 0_yyy_z_yyz_0 0_yyy_z_xz_0 0_yyy_z_xzz_0 0_yyy_z_xy_0 0_yyy_z_xyz_0 0_yyy_z_xx_0 0_yyy_z_xxz_0 0_yyy_zz_zz_0 0_yyy_zz_yz_0 0_yyy_zz_yy_0 0_yyy_zz_xz_0 0_yyy_zz_xy_0 0_yyy_zz_xx_0 0_yyy_y_zz_0 0_yyy_y_zzz_0 0_yyy_y_yz_0 0_yyy_y_yzz_0 0_yyy_y_yy_0 0_yyy_y_yyz_0 0_yyy_y_yyy_0 0_yyy_y_xz_0 0_yyy_y_xzz_0 0_yyy_y_xy_0 0_yyy_y_xyz_0 0_yyy_y_xyy_0 0_yyy_y_xx_0 0_yyy_y_xxz_0 0_yyy_y_xxy_0 0_yyy_yz_zz_0 0_yyy_yz_yz_0 0_yyy_yz_yy_0 0_yyy_yz_xz_0 0_yyy_yz_xy_0 0_yyy_yz_xx_0 0_yyy_yy_zz_0 0_yyy_yy_yz_0 0_yyy_yy_yy_0 0_yyy_yy_xz_0 0_yyy_yy_xy_0 0_yyy_yy_xx_0 0_yyy_x_zz_0 0_yyy_x_zzz_0 0_yyy_x_yz_0 0_yyy_x_yzz_0 0_yyy_x_yy_0 0_yyy_x_yyz_0 0_yyy_x_yyy_0 0_yyy_x_xz_0 0_yyy_x_xzz_0 0_yyy_x_xy_0 0_yyy_x_xyz_0 0_yyy_x_xyy_0 0_yyy_x_xx_0 0_yyy_x_xxz_0 0_yyy_x_xxy_0 0_yyy_x_xxx_0 0_yyy_xz_zz_0 0_yyy_xz_yz_0 0_yyy_xz_yy_0 0_yyy_xz_xz_0 0_yyy_xz_xy_0 0_yyy_xz_xx_0 0_yyy_xy_zz_0 0_yyy_xy_yz_0 0_yyy_xy_yy_0 0_yyy_xy_xz_0 0_yyy_xy_xy_0 0_yyy_xy_xx_0 0_yyy_xx_zz_0 0_yyy_xx_yz_0 0_yyy_xx_yy_0 0_yyy_xx_xz_0 0_yyy_xx_xy_0 0_yyy_xx_xx_0 0_xzz_z_zz_0 0_xzz_z_zzz_0 0_xzz_z_yz_0 0_xzz_z_yzz_0 0_xzz_z_yy_0 0_xzz_z_yyz_0 0_xzz_z_xz_0 0_xzz_z_xzz_0 0_xzz_z_xy_0 0_xzz_z_xyz_0 0_xzz_z_xx_0 0_xzz_z_xxz_0 0_xzz_zz_zz_0 0_xzz_zz_yz_0 0_xzz_zz_yy_0 0_xzz_zz_xz_0 0_xzz_zz_xy_0 0_xzz_zz_xx_0 0_xzz_y_zz_0 0_xzz_y_zzz_0 0_xzz_y_yz_0 0_xzz_y_yzz_0 0_xzz_y_yy_0 0_xzz_y_yyz_0 0_xzz_y_yyy_0 0_xzz_y_xz_0 0_xzz_y_xzz_0 0_xzz_y_xy_0 0_xzz_y_xyz_0 0_xzz_y_xyy_0 0_xzz_y_xx_0 0_xzz_y_xxz_0 0_xzz_y_xxy_0 0_xzz_yz_zz_0 0_xzz_yz_yz_0 0_xzz_yz_yy_0 0_xzz_yz_xz_0 0_xzz_yz_xy_0 0_xzz_yz_xx_0 0_xzz_yy_zz_0 0_xzz_yy_yz_0 0_xzz_yy_yy_0 0_xzz_yy_xz_0 0_xzz_yy_xy_0 0_xzz_yy_xx_0 0_xzz_x_zz_0 0_xzz_x_zzz_0 0_xzz_x_yz_0 0_xzz_x_yzz_0 0_xzz_x_yy_0 0_xzz_x_yyz_0 0_xzz_x_yyy_0 0_xzz_x_xz_0 0_xzz_x_xzz_0 0_xzz_x_xy_0 0_xzz_x_xyz_0 0_xzz_x_xyy_0 0_xzz_x_xx_0 0_xzz_x_xxz_0 0_xzz_x_xxy_0 0_xzz_x_xxx_0 0_xzz_xz_zz_0 0_xzz_xz_yz_0 0_xzz_xz_yy_0 0_xzz_xz_xz_0 0_xzz_xz_xy_0 0_xzz_xz_xx_0 0_xzz_xy_zz_0 0_xzz_xy_yz_0 0_xzz_xy_yy_0 0_xzz_xy_xz_0 0_xzz_xy_xy_0 0_xzz_xy_xx_0 0_xzz_xx_zz_0 0_xzz_xx_yz_0 0_xzz_xx_yy_0 0_xzz_xx_xz_0 0_xzz_xx_xy_0 0_xzz_xx_xx_0 0_xyz_z_zz_0 0_xyz_z_zzz_0 0_xyz_z_yz_0 0_xyz_z_yzz_0 0_xyz_z_yy_0 0_xyz_z_yyz_0 0_xyz_z_xz_0 0_xyz_z_xzz_0 0_xyz_z_xy_0 0_xyz_z_xyz_0 0_xyz_z_xx_0 0_xyz_z_xxz_0 0_xyz_zz_zz_0 0_xyz_zz_yz_0 0_xyz_zz_yy_0 0_xyz_zz_xz_0 0_xyz_zz_xy_0 0_xyz_zz_xx_0 0_xyz_y_zz_0 0_xyz_y_zzz_0 0_xyz_y_yz_0 0_xyz_y_yzz_0 0_xyz_y_yy_0 0_xyz_y_yyz_0 0_xyz_y_yyy_0 0_xyz_y_xz_0 0_xyz_y_xzz_0 0_xyz_y_xy_0 0_xyz_y_xyz_0 0_xyz_y_xyy_0 0_xyz_y_xx_0 0_xyz_y_xxz_0 0_xyz_y_xxy_0 0_xyz_yz_zz_0 0_xyz_yz_yz_0 0_xyz_yz_yy_0 0_xyz_yz_xz_0 0_xyz_yz_xy_0 0_xyz_yz_xx_0 0_xyz_yy_zz_0 0_xyz_yy_yz_0 0_xyz_yy_yy_0 0_xyz_yy_xz_0 0_xyz_yy_xy_0 0_xyz_yy_xx_0 0_xyz_x_zz_0 0_xyz_x_zzz_0 0_xyz_x_yz_0 0_xyz_x_yzz_0 0_xyz_x_yy_0 0_xyz_x_yyz_0 0_xyz_x_yyy_0 0_xyz_x_xz_0 0_xyz_x_xzz_0 0_xyz_x_xy_0 0_xyz_x_xyz_0 0_xyz_x_xyy_0 0_xyz_x_xx_0 0_xyz_x_xxz_0 0_xyz_x_xxy_0 0_xyz_x_xxx_0 0_xyz_xz_zz_0 0_xyz_xz_yz_0 0_xyz_xz_yy_0 0_xyz_xz_xz_0 0_xyz_xz_xy_0 0_xyz_xz_xx_0 0_xyz_xy_zz_0 0_xyz_xy_yz_0 0_xyz_xy_yy_0 0_xyz_xy_xz_0 0_xyz_xy_xy_0 0_xyz_xy_xx_0 0_xyz_xx_zz_0 0_xyz_xx_yz_0 0_xyz_xx_yy_0 0_xyz_xx_xz_0 0_xyz_xx_xy_0 0_xyz_xx_xx_0 0_xyy_z_zz_0 0_xyy_z_zzz_0 0_xyy_z_yz_0 0_xyy_z_yzz_0 0_xyy_z_yy_0 0_xyy_z_yyz_0 0_xyy_z_xz_0 0_xyy_z_xzz_0 0_xyy_z_xy_0 0_xyy_z_xyz_0 0_xyy_z_xx_0 0_xyy_z_xxz_0 0_xyy_zz_zz_0 0_xyy_zz_yz_0 0_xyy_zz_yy_0 0_xyy_zz_xz_0 0_xyy_zz_xy_0 0_xyy_zz_xx_0 0_xyy_y_zz_0 0_xyy_y_zzz_0 0_xyy_y_yz_0 0_xyy_y_yzz_0 0_xyy_y_yy_0 0_xyy_y_yyz_0 0_xyy_y_yyy_0 0_xyy_y_xz_0 0_xyy_y_xzz_0 0_xyy_y_xy_0 0_xyy_y_xyz_0 0_xyy_y_xyy_0 0_xyy_y_xx_0 0_xyy_y_xxz_0 0_xyy_y_xxy_0 0_xyy_yz_zz_0 0_xyy_yz_yz_0 0_xyy_yz_yy_0 0_xyy_yz_xz_0 0_xyy_yz_xy_0 0_xyy_yz_xx_0 0_xyy_yy_zz_0 0_xyy_yy_yz_0 0_xyy_yy_yy_0 0_xyy_yy_xz_0 0_xyy_yy_xy_0 0_xyy_yy_xx_0 0_xyy_x_zz_0 0_xyy_x_zzz_0 0_xyy_x_yz_0 0_xyy_x_yzz_0 0_xyy_x_yy_0 0_xyy_x_yyz_0 0_xyy_x_yyy_0 0_xyy_x_xz_0 0_xyy_x_xzz_0 0_xyy_x_xy_0 0_xyy_x_xyz_0 0_xyy_x_xyy_0 0_xyy_x_xx_0 0_xyy_x_xxz_0 0_xyy_x_xxy_0 0_xyy_x_xxx_0 0_xyy_xz_zz_0 0_xyy_xz_yz_0 0_xyy_xz_yy_0 0_xyy_xz_xz_0 0_xyy_xz_xy_0 0_xyy_xz_xx_0 0_xyy_xy_zz_0 0_xyy_xy_yz_0 0_xyy_xy_yy_0 0_xyy_xy_xz_0 0_xyy_xy_xy_0 0_xyy_xy_xx_0 0_xyy_xx_zz_0 0_xyy_xx_yz_0 0_xyy_xx_yy_0 0_xyy_xx_xz_0 0_xyy_xx_xy_0 0_xyy_xx_xx_0 0_xxz_z_zz_0 0_xxz_z_zzz_0 0_xxz_z_yz_0 0_xxz_z_yzz_0 0_xxz_z_yy_0 0_xxz_z_yyz_0 0_xxz_z_xz_0 0_xxz_z_xzz_0 0_xxz_z_xy_0 0_xxz_z_xyz_0 0_xxz_z_xx_0 0_xxz_z_xxz_0 0_xxz_zz_zz_0 0_xxz_zz_yz_0 0_xxz_zz_yy_0 0_xxz_zz_xz_0 0_xxz_zz_xy_0 0_xxz_zz_xx_0 0_xxz_y_zz_0 0_xxz_y_zzz_0 0_xxz_y_yz_0 0_xxz_y_yzz_0 0_xxz_y_yy_0 0_xxz_y_yyz_0 0_xxz_y_yyy_0 0_xxz_y_xz_0 0_xxz_y_xzz_0 0_xxz_y_xy_0 0_xxz_y_xyz_0 0_xxz_y_xyy_0 0_xxz_y_xx_0 0_xxz_y_xxz_0 0_xxz_y_xxy_0 0_xxz_yz_zz_0 0_xxz_yz_yz_0 0_xxz_yz_yy_0 0_xxz_yz_xz_0 0_xxz_yz_xy_0 0_xxz_yz_xx_0 0_xxz_yy_zz_0 0_xxz_yy_yz_0 0_xxz_yy_yy_0 0_xxz_yy_xz_0 0_xxz_yy_xy_0 0_xxz_yy_xx_0 0_xxz_x_zz_0 0_xxz_x_zzz_0 0_xxz_x_yz_0 0_xxz_x_yzz_0 0_xxz_x_yy_0 0_xxz_x_yyz_0 0_xxz_x_yyy_0 0_xxz_x_xz_0 0_xxz_x_xzz_0 0_xxz_x_xy_0 0_xxz_x_xyz_0 0_xxz_x_xyy_0 0_xxz_x_xx_0 0_xxz_x_xxz_0 0_xxz_x_xxy_0 0_xxz_x_xxx_0 0_xxz_xz_zz_0 0_xxz_xz_yz_0 0_xxz_xz_yy_0 0_xxz_xz_xz_0 0_xxz_xz_xy_0 0_xxz_xz_xx_0 0_xxz_xy_zz_0 0_xxz_xy_yz_0 0_xxz_xy_yy_0 0_xxz_xy_xz_0 0_xxz_xy_xy_0 0_xxz_xy_xx_0 0_xxz_xx_zz_0 0_xxz_xx_yz_0 0_xxz_xx_yy_0 0_xxz_xx_xz_0 0_xxz_xx_xy_0 0_xxz_xx_xx_0 0_xxy_z_zz_0 0_xxy_z_zzz_0 0_xxy_z_yz_0 0_xxy_z_yzz_0 0_xxy_z_yy_0 0_xxy_z_yyz_0 0_xxy_z_xz_0 0_xxy_z_xzz_0 0_xxy_z_xy_0 0_xxy_z_xyz_0 0_xxy_z_xx_0 0_xxy_z_xxz_0 0_xxy_zz_zz_0 0_xxy_zz_yz_0 0_xxy_zz_yy_0 0_xxy_zz_xz_0 0_xxy_zz_xy_0 0_xxy_zz_xx_0 0_xxy_y_zz_0 0_xxy_y_zzz_0 0_xxy_y_yz_0 0_xxy_y_yzz_0 0_xxy_y_yy_0 0_xxy_y_yyz_0 0_xxy_y_yyy_0 0_xxy_y_xz_0 0_xxy_y_xzz_0 0_xxy_y_xy_0 0_xxy_y_xyz_0 0_xxy_y_xyy_0 0_xxy_y_xx_0 0_xxy_y_xxz_0 0_xxy_y_xxy_0 0_xxy_yz_zz_0 0_xxy_yz_yz_0 0_xxy_yz_yy_0 0_xxy_yz_xz_0 0_xxy_yz_xy_0 0_xxy_yz_xx_0 0_xxy_yy_zz_0 0_xxy_yy_yz_0 0_xxy_yy_yy_0 0_xxy_yy_xz_0 0_xxy_yy_xy_0 0_xxy_yy_xx_0 0_xxy_x_zz_0 0_xxy_x_zzz_0 0_xxy_x_yz_0 0_xxy_x_yzz_0 0_xxy_x_yy_0 0_xxy_x_yyz_0 0_xxy_x_yyy_0 0_xxy_x_xz_0 0_xxy_x_xzz_0 0_xxy_x_xy_0 0_xxy_x_xyz_0 0_xxy_x_xyy_0 0_xxy_x_xx_0 0_xxy_x_xxz_0 0_xxy_x_xxy_0 0_xxy_x_xxx_0 0_xxy_xz_zz_0 0_xxy_xz_yz_0 0_xxy_xz_yy_0 0_xxy_xz_xz_0 0_xxy_xz_xy_0 0_xxy_xz_xx_0 0_xxy_xy_zz_0 0_xxy_xy_yz_0 0_xxy_xy_yy_0 0_xxy_xy_xz_0 0_xxy_xy_xy_0 0_xxy_xy_xx_0 0_xxy_xx_zz_0 0_xxy_xx_yz_0 0_xxy_xx_yy_0 0_xxy_xx_xz_0 0_xxy_xx_xy_0 0_xxy_xx_xx_0 0_xxx_z_zz_0 0_xxx_z_zzz_0 0_xxx_z_yz_0 0_xxx_z_yzz_0 0_xxx_z_yy_0 0_xxx_z_yyz_0 0_xxx_z_xz_0 0_xxx_z_xzz_0 0_xxx_z_xy_0 0_xxx_z_xyz_0 0_xxx_z_xx_0 0_xxx_z_xxz_0 0_xxx_zz_zz_0 0_xxx_zz_yz_0 0_xxx_zz_yy_0 0_xxx_zz_xz_0 0_xxx_zz_xy_0 0_xxx_zz_xx_0 0_xxx_y_zz_0 0_xxx_y_zzz_0 0_xxx_y_yz_0 0_xxx_y_yzz_0 0_xxx_y_yy_0 0_xxx_y_yyz_0 0_xxx_y_yyy_0 0_xxx_y_xz_0 0_xxx_y_xzz_0 0_xxx_y_xy_0 0_xxx_y_xyz_0 0_xxx_y_xyy_0 0_xxx_y_xx_0 0_xxx_y_xxz_0 0_xxx_y_xxy_0 0_xxx_yz_zz_0 0_xxx_yz_yz_0 0_xxx_yz_yy_0 0_xxx_yz_xz_0 0_xxx_yz_xy_0 0_xxx_yz_xx_0 0_xxx_yy_zz_0 0_xxx_yy_yz_0 0_xxx_yy_yy_0 0_xxx_yy_xz_0 0_xxx_yy_xy_0 0_xxx_yy_xx_0 0_xxx_x_zz_0 0_xxx_x_zzz_0 0_xxx_x_yz_0 0_xxx_x_yzz_0 0_xxx_x_yy_0 0_xxx_x_yyz_0 0_xxx_x_yyy_0 0_xxx_x_xz_0 0_xxx_x_xzz_0 0_xxx_x_xy_0 0_xxx_x_xyz_0 0_xxx_x_xyy_0 0_xxx_x_xx_0 0_xxx_x_xxz_0 0_xxx_x_xxy_0 0_xxx_x_xxx_0 0_xxx_xz_zz_0 0_xxx_xz_yz_0 0_xxx_xz_yy_0 0_xxx_xz_xz_0 0_xxx_xz_xy_0 0_xxx_xz_xx_0 0_xxx_xy_zz_0 0_xxx_xy_yz_0 0_xxx_xy_yy_0 0_xxx_xy_xz_0 0_xxx_xy_xy_0 0_xxx_xy_xx_0 0_xxx_xx_zz_0 0_xxx_xx_yz_0 0_xxx_xx_yy_0 0_xxx_xx_xz_0 0_xxx_xx_xy_0 0_xxx_xx_xx_0 

signature:
0_zzz_0_zzz_0 0_zzz_0_zzzz_0 0_zzz_0_yzz_0 0_zzz_0_yzzz_0 0_zzz_0_yyz_0 0_zzz_0_yyzz_0 0_zzz_0_yyy_0 0_zzz_0_yyyz_0 0_zzz_0_yyyy_0 0_zzz_0_xzz_0 0_zzz_0_xzzz_0 0_zzz_0_xyz_0 0_zzz_0_xyzz_0 0_zzz_0_xyy_0 0_zzz_0_xyyz_0 0_zzz_0_xyyy_0 0_zzz_0_xxz_0 0_zzz_0_xxzz_0 0_zzz_0_xxy_0 0_zzz_0_xxyz_0 0_zzz_0_xxyy_0 0_zzz_0_xxx_0 0_zzz_0_xxxz_0 0_zzz_0_xxxy_0 0_zzz_0_xxxx_0 0_zzz_z_zzz_0 0_zzz_z_yzz_0 0_zzz_z_yyz_0 0_zzz_z_xzz_0 0_zzz_z_xyz_0 0_zzz_z_xxz_0 0_zzz_y_zzz_0 0_zzz_y_yzz_0 0_zzz_y_yyz_0 0_zzz_y_yyy_0 0_zzz_y_xzz_0 0_zzz_y_xyz_0 0_zzz_y_xyy_0 0_zzz_y_xxz_0 0_zzz_y_xxy_0 0_zzz_x_zzz_0 0_zzz_x_yzz_0 0_zzz_x_yyz_0 0_zzz_x_yyy_0 0_zzz_x_xzz_0 0_zzz_x_xyz_0 0_zzz_x_xyy_0 0_zzz_x_xxz_0 0_zzz_x_xxy_0 0_zzz_x_xxx_0 0_yzz_0_zzz_0 0_yzz_0_zzzz_0 0_yzz_0_yzz_0 0_yzz_0_yzzz_0 0_yzz_0_yyz_0 0_yzz_0_yyzz_0 0_yzz_0_yyy_0 0_yzz_0_yyyz_0 0_yzz_0_yyyy_0 0_yzz_0_xzz_0 0_yzz_0_xzzz_0 0_yzz_0_xyz_0 0_yzz_0_xyzz_0 0_yzz_0_xyy_0 0_yzz_0_xyyz_0 0_yzz_0_xyyy_0 0_yzz_0_xxz_0 0_yzz_0_xxzz_0 0_yzz_0_xxy_0 0_yzz_0_xxyz_0 0_yzz_0_xxyy_0 0_yzz_0_xxx_0 0_yzz_0_xxxz_0 0_yzz_0_xxxy_0 0_yzz_0_xxxx_0 0_yzz_z_zzz_0 0_yzz_z_yzz_0 0_yzz_z_yyz_0 0_yzz_z_xzz_0 0_yzz_z_xyz_0 0_yzz_z_xxz_0 0_yzz_y_zzz_0 0_yzz_y_yzz_0 0_yzz_y_yyz_0 0_yzz_y_yyy_0 0_yzz_y_xzz_0 0_yzz_y_xyz_0 0_yzz_y_xyy_0 0_yzz_y_xxz_0 0_yzz_y_xxy_0 0_yzz_x_zzz_0 0_yzz_x_yzz_0 0_yzz_x_yyz_0 0_yzz_x_yyy_0 0_yzz_x_xzz_0 0_yzz_x_xyz_0 0_yzz_x_xyy_0 0_yzz_x_xxz_0 0_yzz_x_xxy_0 0_yzz_x_xxx_0 0_yyz_0_zzz_0 0_yyz_0_zzzz_0 0_yyz_0_yzz_0 0_yyz_0_yzzz_0 0_yyz_0_yyz_0 0_yyz_0_yyzz_0 0_yyz_0_yyy_0 0_yyz_0_yyyz_0 0_yyz_0_yyyy_0 0_yyz_0_xzz_0 0_yyz_0_xzzz_0 0_yyz_0_xyz_0 0_yyz_0_xyzz_0 0_yyz_0_xyy_0 0_yyz_0_xyyz_0 0_yyz_0_xyyy_0 0_yyz_0_xxz_0 0_yyz_0_xxzz_0 0_yyz_0_xxy_0 0_yyz_0_xxyz_0 0_yyz_0_xxyy_0 0_yyz_0_xxx_0 0_yyz_0_xxxz_0 0_yyz_0_xxxy_0 0_yyz_0_xxxx_0 0_yyz_z_zzz_0 0_yyz_z_yzz_0 0_yyz_z_yyz_0 0_yyz_z_xzz_0 0_yyz_z_xyz_0 0_yyz_z_xxz_0 0_yyz_y_zzz_0 0_yyz_y_yzz_0 0_yyz_y_yyz_0 0_yyz_y_yyy_0 0_yyz_y_xzz_0 0_yyz_y_xyz_0 0_yyz_y_xyy_0 0_yyz_y_xxz_0 0_yyz_y_xxy_0 0_yyz_x_zzz_0 0_yyz_x_yzz_0 0_yyz_x_yyz_0 0_yyz_x_yyy_0 0_yyz_x_xzz_0 0_yyz_x_xyz_0 0_yyz_x_xyy_0 0_yyz_x_xxz_0 0_yyz_x_xxy_0 0_yyz_x_xxx_0 0_yyy_0_zzz_0 0_yyy_0_zzzz_0 0_yyy_0_yzz_0 0_yyy_0_yzzz_0 0_yyy_0_yyz_0 0_yyy_0_yyzz_0 0_yyy_0_yyy_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xzz_0 0_yyy_0_xzzz_0 0_yyy_0_xyz_0 0_yyy_0_xyzz_0 0_yyy_0_xyy_0 0_yyy_0_xyyz_0 0_yyy_0_xyyy_0 0_yyy_0_xxz_0 0_yyy_0_xxzz_0 0_yyy_0_xxy_0 0_yyy_0_xxyz_0 0_yyy_0_xxyy_0 0_yyy_0_xxx_0 0_yyy_0_xxxz_0 0_yyy_0_xxxy_0 0_yyy_0_xxxx_0 0_yyy_z_zzz_0 0_yyy_z_yzz_0 0_yyy_z_yyz_0 0_yyy_z_xzz_0 0_yyy_z_xyz_0 0_yyy_z_xxz_0 0_yyy_y_zzz_0 0_yyy_y_yzz_0 0_yyy_y_yyz_0 0_yyy_y_yyy_0 0_yyy_y_xzz_0 0_yyy_y_xyz_0 0_yyy_y_xyy_0 0_yyy_y_xxz_0 0_yyy_y_xxy_0 0_yyy_x_zzz_0 0_yyy_x_yzz_0 0_yyy_x_yyz_0 0_yyy_x_yyy_0 0_yyy_x_xzz_0 0_yyy_x_xyz_0 0_yyy_x_xyy_0 0_yyy_x_xxz_0 0_yyy_x_xxy_0 0_yyy_x_xxx_0 0_xzz_0_zzz_0 0_xzz_0_zzzz_0 0_xzz_0_yzz_0 0_xzz_0_yzzz_0 0_xzz_0_yyz_0 0_xzz_0_yyzz_0 0_xzz_0_yyy_0 0_xzz_0_yyyz_0 0_xzz_0_yyyy_0 0_xzz_0_xzz_0 0_xzz_0_xzzz_0 0_xzz_0_xyz_0 0_xzz_0_xyzz_0 0_xzz_0_xyy_0 0_xzz_0_xyyz_0 0_xzz_0_xyyy_0 0_xzz_0_xxz_0 0_xzz_0_xxzz_0 0_xzz_0_xxy_0 0_xzz_0_xxyz_0 0_xzz_0_xxyy_0 0_xzz_0_xxx_0 0_xzz_0_xxxz_0 0_xzz_0_xxxy_0 0_xzz_0_xxxx_0 0_xzz_z_zzz_0 0_xzz_z_yzz_0 0_xzz_z_yyz_0 0_xzz_z_xzz_0 0_xzz_z_xyz_0 0_xzz_z_xxz_0 0_xzz_y_zzz_0 0_xzz_y_yzz_0 0_xzz_y_yyz_0 0_xzz_y_yyy_0 0_xzz_y_xzz_0 0_xzz_y_xyz_0 0_xzz_y_xyy_0 0_xzz_y_xxz_0 0_xzz_y_xxy_0 0_xzz_x_zzz_0 0_xzz_x_yzz_0 0_xzz_x_yyz_0 0_xzz_x_yyy_0 0_xzz_x_xzz_0 0_xzz_x_xyz_0 0_xzz_x_xyy_0 0_xzz_x_xxz_0 0_xzz_x_xxy_0 0_xzz_x_xxx_0 0_xyz_0_zzz_0 0_xyz_0_zzzz_0 0_xyz_0_yzz_0 0_xyz_0_yzzz_0 0_xyz_0_yyz_0 0_xyz_0_yyzz_0 0_xyz_0_yyy_0 0_xyz_0_yyyz_0 0_xyz_0_yyyy_0 0_xyz_0_xzz_0 0_xyz_0_xzzz_0 0_xyz_0_xyz_0 0_xyz_0_xyzz_0 0_xyz_0_xyy_0 0_xyz_0_xyyz_0 0_xyz_0_xyyy_0 0_xyz_0_xxz_0 0_xyz_0_xxzz_0 0_xyz_0_xxy_0 0_xyz_0_xxyz_0 0_xyz_0_xxyy_0 0_xyz_0_xxx_0 0_xyz_0_xxxz_0 0_xyz_0_xxxy_0 0_xyz_0_xxxx_0 0_xyz_z_zzz_0 0_xyz_z_yzz_0 0_xyz_z_yyz_0 0_xyz_z_xzz_0 0_xyz_z_xyz_0 0_xyz_z_xxz_0 0_xyz_y_zzz_0 0_xyz_y_yzz_0 0_xyz_y_yyz_0 0_xyz_y_yyy_0 0_xyz_y_xzz_0 0_xyz_y_xyz_0 0_xyz_y_xyy_0 0_xyz_y_xxz_0 0_xyz_y_xxy_0 0_xyz_x_zzz_0 0_xyz_x_yzz_0 0_xyz_x_yyz_0 0_xyz_x_yyy_0 0_xyz_x_xzz_0 0_xyz_x_xyz_0 0_xyz_x_xyy_0 0_xyz_x_xxz_0 0_xyz_x_xxy_0 0_xyz_x_xxx_0 0_xyy_0_zzz_0 0_xyy_0_zzzz_0 0_xyy_0_yzz_0 0_xyy_0_yzzz_0 0_xyy_0_yyz_0 0_xyy_0_yyzz_0 0_xyy_0_yyy_0 0_xyy_0_yyyz_0 0_xyy_0_yyyy_0 0_xyy_0_xzz_0 0_xyy_0_xzzz_0 0_xyy_0_xyz_0 0_xyy_0_xyzz_0 0_xyy_0_xyy_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxz_0 0_xyy_0_xxzz_0 0_xyy_0_xxy_0 0_xyy_0_xxyz_0 0_xyy_0_xxyy_0 0_xyy_0_xxx_0 0_xyy_0_xxxz_0 0_xyy_0_xxxy_0 0_xyy_0_xxxx_0 0_xyy_z_zzz_0 0_xyy_z_yzz_0 0_xyy_z_yyz_0 0_xyy_z_xzz_0 0_xyy_z_xyz_0 0_xyy_z_xxz_0 0_xyy_y_zzz_0 0_xyy_y_yzz_0 0_xyy_y_yyz_0 0_xyy_y_yyy_0 0_xyy_y_xzz_0 0_xyy_y_xyz_0 0_xyy_y_xyy_0 0_xyy_y_xxz_0 0_xyy_y_xxy_0 0_xyy_x_zzz_0 0_xyy_x_yzz_0 0_xyy_x_yyz_0 0_xyy_x_yyy_0 0_xyy_x_xzz_0 0_xyy_x_xyz_0 0_xyy_x_xyy_0 0_xyy_x_xxz_0 0_xyy_x_xxy_0 0_xyy_x_xxx_0 0_xxz_0_zzz_0 0_xxz_0_zzzz_0 0_xxz_0_yzz_0 0_xxz_0_yzzz_0 0_xxz_0_yyz_0 0_xxz_0_yyzz_0 0_xxz_0_yyy_0 0_xxz_0_yyyz_0 0_xxz_0_yyyy_0 0_xxz_0_xzz_0 0_xxz_0_xzzz_0 0_xxz_0_xyz_0 0_xxz_0_xyzz_0 0_xxz_0_xyy_0 0_xxz_0_xyyz_0 0_xxz_0_xyyy_0 0_xxz_0_xxz_0 0_xxz_0_xxzz_0 0_xxz_0_xxy_0 0_xxz_0_xxyz_0 0_xxz_0_xxyy_0 0_xxz_0_xxx_0 0_xxz_0_xxxz_0 0_xxz_0_xxxy_0 0_xxz_0_xxxx_0 0_xxz_z_zzz_0 0_xxz_z_yzz_0 0_xxz_z_yyz_0 0_xxz_z_xzz_0 0_xxz_z_xyz_0 0_xxz_z_xxz_0 0_xxz_y_zzz_0 0_xxz_y_yzz_0 0_xxz_y_yyz_0 0_xxz_y_yyy_0 0_xxz_y_xzz_0 0_xxz_y_xyz_0 0_xxz_y_xyy_0 0_xxz_y_xxz_0 0_xxz_y_xxy_0 0_xxz_x_zzz_0 0_xxz_x_yzz_0 0_xxz_x_yyz_0 0_xxz_x_yyy_0 0_xxz_x_xzz_0 0_xxz_x_xyz_0 0_xxz_x_xyy_0 0_xxz_x_xxz_0 0_xxz_x_xxy_0 0_xxz_x_xxx_0 0_xxy_0_zzz_0 0_xxy_0_zzzz_0 0_xxy_0_yzz_0 0_xxy_0_yzzz_0 0_xxy_0_yyz_0 0_xxy_0_yyzz_0 0_xxy_0_yyy_0 0_xxy_0_yyyz_0 0_xxy_0_yyyy_0 0_xxy_0_xzz_0 0_xxy_0_xzzz_0 0_xxy_0_xyz_0 0_xxy_0_xyzz_0 0_xxy_0_xyy_0 0_xxy_0_xyyz_0 0_xxy_0_xyyy_0 0_xxy_0_xxz_0 0_xxy_0_xxzz_0 0_xxy_0_xxy_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxx_0 0_xxy_0_xxxz_0 0_xxy_0_xxxy_0 0_xxy_0_xxxx_0 0_xxy_z_zzz_0 0_xxy_z_yzz_0 0_xxy_z_yyz_0 0_xxy_z_xzz_0 0_xxy_z_xyz_0 0_xxy_z_xxz_0 0_xxy_y_zzz_0 0_xxy_y_yzz_0 0_xxy_y_yyz_0 0_xxy_y_yyy_0 0_xxy_y_xzz_0 0_xxy_y_xyz_0 0_xxy_y_xyy_0 0_xxy_y_xxz_0 0_xxy_y_xxy_0 0_xxy_x_zzz_0 0_xxy_x_yzz_0 0_xxy_x_yyz_0 0_xxy_x_yyy_0 0_xxy_x_xzz_0 0_xxy_x_xyz_0 0_xxy_x_xyy_0 0_xxy_x_xxz_0 0_xxy_x_xxy_0 0_xxy_x_xxx_0 0_xxx_0_zzz_0 0_xxx_0_zzzz_0 0_xxx_0_yzz_0 0_xxx_0_yzzz_0 0_xxx_0_yyz_0 0_xxx_0_yyzz_0 0_xxx_0_yyy_0 0_xxx_0_yyyz_0 0_xxx_0_yyyy_0 0_xxx_0_xzz_0 0_xxx_0_xzzz_0 0_xxx_0_xyz_0 0_xxx_0_xyzz_0 0_xxx_0_xyy_0 0_xxx_0_xyyz_0 0_xxx_0_xyyy_0 0_xxx_0_xxz_0 0_xxx_0_xxzz_0 0_xxx_0_xxy_0 0_xxx_0_xxyz_0 0_xxx_0_xxyy_0 0_xxx_0_xxx_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 0_xxx_z_zzz_0 0_xxx_z_yzz_0 0_xxx_z_yyz_0 0_xxx_z_xzz_0 0_xxx_z_xyz_0 0_xxx_z_xxz_0 0_xxx_y_zzz_0 0_xxx_y_yzz_0 0_xxx_y_yyz_0 0_xxx_y_yyy_0 0_xxx_y_xzz_0 0_xxx_y_xyz_0 0_xxx_y_xyy_0 0_xxx_y_xxz_0 0_xxx_y_xxy_0 0_xxx_x_zzz_0 0_xxx_x_yzz_0 0_xxx_x_yyz_0 0_xxx_x_yyy_0 0_xxx_x_xzz_0 0_xxx_x_xyz_0 0_xxx_x_xyy_0 0_xxx_x_xxz_0 0_xxx_x_xxy_0 0_xxx_x_xxx_0 

signature:
0_zzz_0_zz_0 0_zzz_0_zzz_0 0_zzz_0_yz_0 0_zzz_0_yzz_0 0_zzz_0_yy_0 0_zzz_0_yyz_0 0_zzz_0_yyy_0 0_zzz_0_xz_0 0_zzz_0_xzz_0 0_zzz_0_xy_0 0_zzz_0_xyz_0 0_zzz_0_xyy_0 0_zzz_0_xx_0 0_zzz_0_xxz_0 0_zzz_0_xxy_0 0_zzz_0_xxx_0 0_zzz_z_zz_0 0_zzz_z_yz_0 0_zzz_z_yy_0 0_zzz_z_xz_0 0_zzz_z_xy_0 0_zzz_z_xx_0 0_zzz_y_zz_0 0_zzz_y_yz_0 0_zzz_y_yy_0 0_zzz_y_xz_0 0_zzz_y_xy_0 0_zzz_y_xx_0 0_zzz_x_zz_0 0_zzz_x_yz_0 0_zzz_x_yy_0 0_zzz_x_xz_0 0_zzz_x_xy_0 0_zzz_x_xx_0 0_yzz_0_zz_0 0_yzz_0_zzz_0 0_yzz_0_yz_0 0_yzz_0_yzz_0 0_yzz_0_yy_0 0_yzz_0_yyz_0 0_yzz_0_yyy_0 0_yzz_0_xz_0 0_yzz_0_xzz_0 0_yzz_0_xy_0 0_yzz_0_xyz_0 0_yzz_0_xyy_0 0_yzz_0_xx_0 0_yzz_0_xxz_0 0_yzz_0_xxy_0 0_yzz_0_xxx_0 0_yzz_z_zz_0 0_yzz_z_yz_0 0_yzz_z_yy_0 0_yzz_z_xz_0 0_yzz_z_xy_0 0_yzz_z_xx_0 0_yzz_y_zz_0 0_yzz_y_yz_0 0_yzz_y_yy_0 0_yzz_y_xz_0 0_yzz_y_xy_0 0_yzz_y_xx_0 0_yzz_x_zz_0 0_yzz_x_yz_0 0_yzz_x_yy_0 0_yzz_x_xz_0 0_yzz_x_xy_0 0_yzz_x_xx_0 0_yyz_0_zz_0 0_yyz_0_zzz_0 0_yyz_0_yz_0 0_yyz_0_yzz_0 0_yyz_0_yy_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xz_0 0_yyz_0_xzz_0 0_yyz_0_xy_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xx_0 0_yyz_0_xxz_0 0_yyz_0_xxy_0 0_yyz_0_xxx_0 0_yyz_z_zz_0 0_yyz_z_yz_0 0_yyz_z_yy_0 0_yyz_z_xz_0 0_yyz_z_xy_0 0_yyz_z_xx_0 0_yyz_y_zz_0 0_yyz_y_yz_0 0_yyz_y_yy_0 0_yyz_y_xz_0 0_yyz_y_xy_0 0_yyz_y_xx_0 0_yyz_x_zz_0 0_yyz_x_yz_0 0_yyz_x_yy_0 0_yyz_x_xz_0 0_yyz_x_xy_0 0_yyz_x_xx_0 0_yyy_0_zz_0 0_yyy_0_zzz_0 0_yyy_0_yz_0 0_yyy_0_yzz_0 0_yyy_0_yy_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xz_0 0_yyy_0_xzz_0 0_yyy_0_xy_0 0_yyy_0_xyz_0 0_yyy_0_xyy_0 0_yyy_0_xx_0 0_yyy_0_xxz_0 0_yyy_0_xxy_0 0_yyy_0_xxx_0 0_yyy_z_zz_0 0_yyy_z_yz_0 0_yyy_z_yy_0 0_yyy_z_xz_0 0_yyy_z_xy_0 0_yyy_z_xx_0 0_yyy_y_zz_0 0_yyy_y_yz_0 0_yyy_y_yy_0 0_yyy_y_xz_0 0_yyy_y_xy_0 0_yyy_y_xx_0 0_yyy_x_zz_0 0_yyy_x_yz_0 0_yyy_x_yy_0 0_yyy_x_xz_0 0_yyy_x_xy_0 0_yyy_x_xx_0 0_xzz_0_zz_0 0_xzz_0_zzz_0 0_xzz_0_yz_0 0_xzz_0_yzz_0 0_xzz_0_yy_0 0_xzz_0_yyz_0 0_xzz_0_yyy_0 0_xzz_0_xz_0 0_xzz_0_xzz_0 0_xzz_0_xy_0 0_xzz_0_xyz_0 0_xzz_0_xyy_0 0_xzz_0_xx_0 0_xzz_0_xxz_0 0_xzz_0_xxy_0 0_xzz_0_xxx_0 0_xzz_z_zz_0 0_xzz_z_yz_0 0_xzz_z_yy_0 0_xzz_z_xz_0 0_xzz_z_xy_0 0_xzz_z_xx_0 0_xzz_y_zz_0 0_xzz_y_yz_0 0_xzz_y_yy_0 0_xzz_y_xz_0 0_xzz_y_xy_0 0_xzz_y_xx_0 0_xzz_x_zz_0 0_xzz_x_yz_0 0_xzz_x_yy_0 0_xzz_x_xz_0 0_xzz_x_xy_0 0_xzz_x_xx_0 0_xyz_0_zz_0 0_xyz_0_zzz_0 0_xyz_0_yz_0 0_xyz_0_yzz_0 0_xyz_0_yy_0 0_xyz_0_yyz_0 0_xyz_0_yyy_0 0_xyz_0_xz_0 0_xyz_0_xzz_0 0_xyz_0_xy_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xx_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxx_0 0_xyz_z_zz_0 0_xyz_z_yz_0 0_xyz_z_yy_0 0_xyz_z_xz_0 0_xyz_z_xy_0 0_xyz_z_xx_0 0_xyz_y_zz_0 0_xyz_y_yz_0 0_xyz_y_yy_0 0_xyz_y_xz_0 0_xyz_y_xy_0 0_xyz_y_xx_0 0_xyz_x_zz_0 0_xyz_x_yz_0 0_xyz_x_yy_0 0_xyz_x_xz_0 0_xyz_x_xy_0 0_xyz_x_xx_0 0_xyy_0_zz_0 0_xyy_0_zzz_0 0_xyy_0_yz_0 0_xyy_0_yzz_0 0_xyy_0_yy_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xz_0 0_xyy_0_xzz_0 0_xyy_0_xy_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xx_0 0_xyy_0_xxz_0 0_xyy_0_xxy_0 0_xyy_0_xxx_0 0_xyy_z_zz_0 0_xyy_z_yz_0 0_xyy_z_yy_0 0_xyy_z_xz_0 0_xyy_z_xy_0 0_xyy_z_xx_0 0_xyy_y_zz_0 0_xyy_y_yz_0 0_xyy_y_yy_0 0_xyy_y_xz_0 0_xyy_y_xy_0 0_xyy_y_xx_0 0_xyy_x_zz_0 0_xyy_x_yz_0 0_xyy_x_yy_0 0_xyy_x_xz_0 0_xyy_x_xy_0 0_xyy_x_xx_0 0_xxz_0_zz_0 0_xxz_0_zzz_0 0_xxz_0_yz_0 0_xxz_0_yzz_0 0_xxz_0_yy_0 0_xxz_0_yyz_0 0_xxz_0_yyy_0 0_xxz_0_xz_0 0_xxz_0_xzz_0 0_xxz_0_xy_0 0_xxz_0_xyz_0 0_xxz_0_xyy_0 0_xxz_0_xx_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxz_z_zz_0 0_xxz_z_yz_0 0_xxz_z_yy_0 0_xxz_z_xz_0 0_xxz_z_xy_0 0_xxz_z_xx_0 0_xxz_y_zz_0 0_xxz_y_yz_0 0_xxz_y_yy_0 0_xxz_y_xz_0 0_xxz_y_xy_0 0_xxz_y_xx_0 0_xxz_x_zz_0 0_xxz_x_yz_0 0_xxz_x_yy_0 0_xxz_x_xz_0 0_xxz_x_xy_0 0_xxz_x_xx_0 0_xxy_0_zz_0 0_xxy_0_zzz_0 0_xxy_0_yz_0 0_xxy_0_yzz_0 0_xxy_0_yy_0 0_xxy_0_yyz_0 0_xxy_0_yyy_0 0_xxy_0_xz_0 0_xxy_0_xzz_0 0_xxy_0_xy_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xx_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxy_z_zz_0 0_xxy_z_yz_0 0_xxy_z_yy_0 0_xxy_z_xz_0 0_xxy_z_xy_0 0_xxy_z_xx_0 0_xxy_y_zz_0 0_xxy_y_yz_0 0_xxy_y_yy_0 0_xxy_y_xz_0 0_xxy_y_xy_0 0_xxy_y_xx_0 0_xxy_x_zz_0 0_xxy_x_yz_0 0_xxy_x_yy_0 0_xxy_x_xz_0 0_xxy_x_xy_0 0_xxy_x_xx_0 0_xxx_0_zz_0 0_xxx_0_zzz_0 0_xxx_0_yz_0 0_xxx_0_yzz_0 0_xxx_0_yy_0 0_xxx_0_yyz_0 0_xxx_0_yyy_0 0_xxx_0_xz_0 0_xxx_0_xzz_0 0_xxx_0_xy_0 0_xxx_0_xyz_0 0_xxx_0_xyy_0 0_xxx_0_xx_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 0_xxx_z_zz_0 0_xxx_z_yz_0 0_xxx_z_yy_0 0_xxx_z_xz_0 0_xxx_z_xy_0 0_xxx_z_xx_0 0_xxx_y_zz_0 0_xxx_y_yz_0 0_xxx_y_yy_0 0_xxx_y_xz_0 0_xxx_y_xy_0 0_xxx_y_xx_0 0_xxx_x_zz_0 0_xxx_x_yz_0 0_xxx_x_yy_0 0_xxx_x_xz_0 0_xxx_x_xy_0 0_xxx_x_xx_0 

signature:
0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_yyyy_0 0_z_0_yyyy_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xyyy_0 0_z_0_xyyy_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxyy_0 0_z_0_xxyy_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_z_0_xxxy_0 0_z_0_xxxy_1 0_z_0_xxxx_0 0_z_0_xxxx_1 0_zz_0_zzz_1 0_zz_0_zzzz_0 0_zz_0_zzzz_1 0_zz_0_yzz_1 0_zz_0_yzzz_0 0_zz_0_yzzz_1 0_zz_0_yyz_1 0_zz_0_yyzz_0 0_zz_0_yyzz_1 0_zz_0_yyy_1 0_zz_0_yyyz_0 0_zz_0_yyyz_1 0_zz_0_yyyy_0 0_zz_0_yyyy_1 0_zz_0_xzz_1 0_zz_0_xzzz_0 0_zz_0_xzzz_1 0_zz_0_xyz_1 0_zz_0_xyzz_0 0_zz_0_xyzz_1 0_zz_0_xyy_1 0_zz_0_xyyz_0 0_zz_0_xyyz_1 0_zz_0_xyyy_0 0_zz_0_xyyy_1 0_zz_0_xxz_1 0_zz_0_xxzz_0 0_zz_0_xxzz_1 0_zz_0_xxy_1 0_zz_0_xxyz_0 0_zz_0_xxyz_1 0_zz_0_xxyy_0 0_zz_0_xxyy_1 0_zz_0_xxx_1 0_zz_0_xxxz_0 0_zz_0_xxxz_1 0_zz_0_xxxy_0 0_zz_0_xxxy_1 0_zz_0_xxxx_0 0_zz_0_xxxx_1 0_zzz_0_zzzz_0 0_zzz_0_yzzz_0 0_zzz_0_yyzz_0 0_zzz_0_yyyz_0 0_zzz_0_yyyy_0 0_zzz_0_xzzz_0 0_zzz_0_xyzz_0 0_zzz_0_xyyz_0 0_zzz_0_xyyy_0 0_zzz_0_xxzz_0 0_zzz_0_xxyz_0 0_zzz_0_xxyy_0 0_zzz_0_xxxz_0 0_zzz_0_xxxy_0 0_zzz_0_xxxx_0 0_y_0_zzzz_0 0_y_0_zzzz_1 0_y_0_yzzz_0 0_y_0_yzzz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xzzz_0 0_y_0_xzzz_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxzz_0 0_y_0_xxzz_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxxz_0 0_y_0_xxxz_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_y_0_xxxx_0 0_y_0_xxxx_1 0_yz_0_zzz_1 0_yz_0_zzzz_0 0_yz_0_zzzz_1 0_yz_0_yzz_1 0_yz_0_yzzz_0 0_yz_0_yzzz_1 0_yz_0_yyz_1 0_yz_0_yyzz_0 0_yz_0_yyzz_1 0_yz_0_yyy_1 0_yz_0_yyyz_0 0_yz_0_yyyz_1 0_yz_0_yyyy_0 0_yz_0_yyyy_1 0_yz_0_xzzz_0 0_yz_0_xyz_1 0_yz_0_xyzz_0 0_yz_0_xyzz_1 0_yz_0_xyyz_0 0_yz_0_xyyz_1 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyz_1 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yzz_0_zzzz_0 0_yzz_0_yzzz_0 0_yzz_0_yyzz_0 0_yzz_0_yyyz_0 0_yzz_0_yyyy_0 0_yzz_0_xzzz_0 0_yzz_0_xyzz_0 0_yzz_0_xyyz_0 0_yzz_0_xyyy_0 0_yzz_0_xxzz_0 0_yzz_0_xxyz_0 0_yzz_0_xxyy_0 0_yzz_0_xxxz_0 0_yzz_0_xxxy_0 0_yzz_0_xxxx_0 0_yy_0_zzz_1 0_yy_0_zzzz_0 0_yy_0_zzzz_1 0_yy_0_yzz_1 0_yy_0_yzzz_0 0_yy_0_yzzz_1 0_yy_0_yyz_1 0_yy_0_yyzz_0 0_yy_0_yyzz_1 0_yy_0_yyy_1 0_yy_0_yyyz_0 0_yy_0_yyyz_1 0_yy_0_yyyy_0 0_yy_0_yyyy_1 0_yy_0_xzz_1 0_yy_0_xzzz_0 0_yy_0_xzzz_1 0_yy_0_xyz_1 0_yy_0_xyzz_0 0_yy_0_xyzz_1 0_yy_0_xyy_1 0_yy_0_xyyz_0 0_yy_0_xyyz_1 0_yy_0_xyyy_0 0_yy_0_xyyy_1 0_yy_0_xxz_1 0_yy_0_xxzz_0 0_yy_0_xxzz_1 0_yy_0_xxy_1 0_yy_0_xxyz_0 0_yy_0_xxyz_1 0_yy_0_xxyy_0 0_yy_0_xxyy_1 0_yy_0_xxx_1 0_yy_0_xxxz_0 0_yy_0_xxxz_1 0_yy_0_xxxy_0 0_yy_0_xxxy_1 0_yy_0_xxxx_0 0_yy_0_xxxx_1 0_yyz_0_zzzz_0 0_yyz_0_yzzz_0 0_yyz_0_yyzz_0 0_yyz_0_yyyz_0 0_yyz_0_yyyy_0 0_yyz_0_xzzz_0 0_yyz_0_xyzz_0 0_yyz_0_xyyz_0 0_yyz_0_xyyy_0 0_yyz_0_xxzz_0 0_yyz_0_xxyz_0 0_yyz_0_xxyy_0 0_yyz_0_xxxz_0 0_yyz_0_xxxy_0 0_yyz_0_xxxx_0 0_yyy_0_zzzz_0 0_yyy_0_yzzz_0 0_yyy_0_yyzz_0 0_yyy_0_yyyz_0 0_yyy_0_yyyy_0 0_yyy_0_xzzz_0 0_yyy_0_xyzz_0 0_yyy_0_xyyz_0 0_yyy_0_xyyy_0 0_yyy_0_xxzz_0 0_yyy_0_xxyz_0 0_yyy_0_xxyy_0 0_yyy_0_xxxz_0 0_yyy_0_xxxy_0 0_yyy_0_xxxx_0 0_x_0_zzzz_0 0_x_0_zzzz_1 0_x_0_yzzz_0 0_x_0_yzzz_1 0_x_0_yyzz_0 0_x_0_yyzz_1 0_x_0_yyyz_0 0_x_0_yyyz_1 0_x_0_yyyy_0 0_x_0_yyyy_1 0_x_0_xzzz_0 0_x_0_xzzz_1 0_x_0_xyzz_0 0_x_0_xyzz_1 0_x_0_xyyz_0 0_x_0_xyyz_1 0_x_0_xyyy_0 0_x_0_xyyy_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzz_1 0_xz_0_xzzz_0 0_xz_0_xzzz_1 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxz_1 0_xz_0_xxzz_0 0_xz_0_xxzz_1 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxx_1 0_xz_0_xxxz_0 0_xz_0_xxxz_1 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xz_0_xxxx_1 0_xzz_0_zzzz_0 0_xzz_0_yzzz_0 0_xzz_0_yyzz_0 0_xzz_0_yyyz_0 0_xzz_0_yyyy_0 0_xzz_0_xzzz_0 0_xzz_0_xyzz_0 0_xzz_0_xyyz_0 0_xzz_0_xyyy_0 0_xzz_0_xxzz_0 0_xzz_0_xxyz_0 0_xzz_0_xxyy_0 0_xzz_0_xxxz_0 0_xzz_0_xxxy_0 0_xzz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyy_1 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xyyy_1 0_xy_0_xxzz_0 0_xy_0_xxy_1 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxyy_1 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxy_1 0_xy_0_xxxx_0 0_xyz_0_zzzz_0 0_xyz_0_yzzz_0 0_xyz_0_yyzz_0 0_xyz_0_yyyz_0 0_xyz_0_yyyy_0 0_xyz_0_xzzz_0 0_xyz_0_xyzz_0 0_xyz_0_xyyz_0 0_xyz_0_xyyy_0 0_xyz_0_xxzz_0 0_xyz_0_xxyz_0 0_xyz_0_xxyy_0 0_xyz_0_xxxz_0 0_xyz_0_xxxy_0 0_xyz_0_xxxx_0 0_xyy_0_zzzz_0 0_xyy_0_yzzz_0 0_xyy_0_yyzz_0 0_xyy_0_yyyz_0 0_xyy_0_yyyy_0 0_xyy_0_xzzz_0 0_xyy_0_xyzz_0 0_xyy_0_xyyz_0 0_xyy_0_xyyy_0 0_xyy_0_xxzz_0 0_xyy_0_xxyz_0 0_xyy_0_xxyy_0 0_xyy_0_xxxz_0 0_xyy_0_xxxy_0 0_xyy_0_xxxx_0 0_xx_0_zzz_1 0_xx_0_zzzz_0 0_xx_0_zzzz_1 0_xx_0_yzz_1 0_xx_0_yzzz_0 0_xx_0_yzzz_1 0_xx_0_yyz_1 0_xx_0_yyzz_0 0_xx_0_yyzz_1 0_xx_0_yyy_1 0_xx_0_yyyz_0 0_xx_0_yyyz_1 0_xx_0_yyyy_0 0_xx_0_yyyy_1 0_xx_0_xzz_1 0_xx_0_xzzz_0 0_xx_0_xzzz_1 0_xx_0_xyz_1 0_xx_0_xyzz_0 0_xx_0_xyzz_1 0_xx_0_xyy_1 0_xx_0_xyyz_0 0_xx_0_xyyz_1 0_xx_0_xyyy_0 0_xx_0_xyyy_1 0_xx_0_xxz_1 0_xx_0_xxzz_0 0_xx_0_xxzz_1 0_xx_0_xxy_1 0_xx_0_xxyz_0 0_xx_0_xxyz_1 0_xx_0_xxyy_0 0_xx_0_xxyy_1 0_xx_0_xxx_1 0_xx_0_xxxz_0 0_xx_0_xxxz_1 0_xx_0_xxxy_0 0_xx_0_xxxy_1 0_xx_0_xxxx_0 0_xx_0_xxxx_1 0_xxz_0_zzzz_0 0_xxz_0_yzzz_0 0_xxz_0_yyzz_0 0_xxz_0_yyyz_0 0_xxz_0_yyyy_0 0_xxz_0_xzzz_0 0_xxz_0_xyzz_0 0_xxz_0_xyyz_0 0_xxz_0_xyyy_0 0_xxz_0_xxzz_0 0_xxz_0_xxyz_0 0_xxz_0_xxyy_0 0_xxz_0_xxxz_0 0_xxz_0_xxxy_0 0_xxz_0_xxxx_0 0_xxy_0_zzzz_0 0_xxy_0_yzzz_0 0_xxy_0_yyzz_0 0_xxy_0_yyyz_0 0_xxy_0_yyyy_0 0_xxy_0_xzzz_0 0_xxy_0_xyzz_0 0_xxy_0_xyyz_0 0_xxy_0_xyyy_0 0_xxy_0_xxzz_0 0_xxy_0_xxyz_0 0_xxy_0_xxyy_0 0_xxy_0_xxxz_0 0_xxy_0_xxxy_0 0_xxy_0_xxxx_0 0_xxx_0_zzzz_0 0_xxx_0_yzzz_0 0_xxx_0_yyzz_0 0_xxx_0_yyyz_0 0_xxx_0_yyyy_0 0_xxx_0_xzzz_0 0_xxx_0_xyzz_0 0_xxx_0_xyyz_0 0_xxx_0_xyyy_0 0_xxx_0_xxzz_0 0_xxx_0_xxyz_0 0_xxx_0_xxyy_0 0_xxx_0_xxxz_0 0_xxx_0_xxxy_0 0_xxx_0_xxxx_0 

signature:
0_z_0_zzz_0 0_z_0_zzz_1 0_z_0_yzz_0 0_z_0_yzz_1 0_z_0_yyz_0 0_z_0_yyz_1 0_z_0_yyy_0 0_z_0_yyy_1 0_z_0_xzz_0 0_z_0_xzz_1 0_z_0_xyz_0 0_z_0_xyz_1 0_z_0_xyy_0 0_z_0_xyy_1 0_z_0_xxz_0 0_z_0_xxz_1 0_z_0_xxy_0 0_z_0_xxy_1 0_z_0_xxx_0 0_z_0_xxx_1 0_zz_0_zz_1 0_zz_0_zzz_0 0_zz_0_zzz_1 0_zz_0_yz_1 0_zz_0_yzz_0 0_zz_0_yzz_1 0_zz_0_yy_1 0_zz_0_yyz_0 0_zz_0_yyz_1 0_zz_0_yyy_0 0_zz_0_yyy_1 0_zz_0_xz_1 0_zz_0_xzz_0 0_zz_0_xzz_1 0_zz_0_xy_1 0_zz_0_xyz_0 0_zz_0_xyz_1 0_zz_0_xyy_0 0_zz_0_xyy_1 0_zz_0_xx_1 0_zz_0_xxz_0 0_zz_0_xxz_1 0_zz_0_xxy_0 0_zz_0_xxy_1 0_zz_0_xxx_0 0_zz_0_xxx_1 0_zzz_0_zzz_0 0_zzz_0_yzz_0 0_zzz_0_yyz_0 0_zzz_0_yyy_0 0_zzz_0_xzz_0 0_zzz_0_xyz_0 0_zzz_0_xyy_0 0_zzz_0_xxz_0 0_zzz_0_xxy_0 0_zzz_0_xxx_0 0_y_0_zzz_0 0_y_0_zzz_1 0_y_0_yzz_0 0_y_0_yzz_1 0_y_0_yyz_0 0_y_0_yyz_1 0_y_0_yyy_0 0_y_0_yyy_1 0_y_0_xzz_0 0_y_0_xzz_1 0_y_0_xyz_0 0_y_0_xyz_1 0_y_0_xyy_0 0_y_0_xyy_1 0_y_0_xxz_0 0_y_0_xxz_1 0_y_0_xxy_0 0_y_0_xxy_1 0_y_0_xxx_0 0_y_0_xxx_1 0_yz_0_zz_1 0_yz_0_zzz_0 0_yz_0_zzz_1 0_yz_0_yz_1 0_yz_0_yzz_0 0_yz_0_yzz_1 0_yz_0_yy_1 0_yz_0_yyz_0 0_yz_0_yyz_1 0_yz_0_yyy_0 0_yz_0_yyy_1 0_yz_0_xzz_0 0_yz_0_xyz_0 0_yz_0_xyz_1 0_yz_0_xyy_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yzz_0_zzz_0 0_yzz_0_yzz_0 0_yzz_0_yyz_0 0_yzz_0_yyy_0 0_yzz_0_xzz_0 0_yzz_0_xyz_0 0_yzz_0_xyy_0 0_yzz_0_xxz_0 0_yzz_0_xxy_0 0_yzz_0_xxx_0 0_yy_0_zz_1 0_yy_0_zzz_0 0_yy_0_zzz_1 0_yy_0_yz_1 0_yy_0_yzz_0 0_yy_0_yzz_1 0_yy_0_yy_1 0_yy_0_yyz_0 0_yy_0_yyz_1 0_yy_0_yyy_0 0_yy_0_yyy_1 0_yy_0_xz_1 0_yy_0_xzz_0 0_yy_0_xzz_1 0_yy_0_xy_1 0_yy_0_xyz_0 0_yy_0_xyz_1 0_yy_0_xyy_0 0_yy_0_xyy_1 0_yy_0_xx_1 0_yy_0_xxz_0 0_yy_0_xxz_1 0_yy_0_xxy_0 0_yy_0_xxy_1 0_yy_0_xxx_0 0_yy_0_xxx_1 0_yyz_0_zzz_0 0_yyz_0_yzz_0 0_yyz_0_yyz_0 0_yyz_0_yyy_0 0_yyz_0_xzz_0 0_yyz_0_xyz_0 0_yyz_0_xyy_0 0_yyz_0_xxz_0 0_yyz_0_xxy_0 0_yyz_0_xxx_0 0_yyy_0_zzz_0 0_yyy_0_yzz_0 0_yyy_0_yyz_0 0_yyy_0_yyy_0 0_yyy_0_xzz_0 0_yyy_0_xyz_0 0_yyy_0_xyy_0 0_yyy_0_xxz_0 0_yyy_0_xxy_0 0_yyy_0_xxx_0 0_x_0_zzz_0 0_x_0_zzz_1 0_x_0_yzz_0 0_x_0_yzz_1 0_x_0_yyz_0 0_x_0_yyz_1 0_x_0_yyy_0 0_x_0_yyy_1 0_x_0_xzz_0 0_x_0_xzz_1 0_x_0_xyz_0 0_x_0_xyz_1 0_x_0_xyy_0 0_x_0_xyy_1 0_x_0_xxz_0 0_x_0_xxz_1 0_x_0_xxy_0 0_x_0_xxy_1 0_x_0_xxx_0 0_x_0_xxx_1 0_xz_0_zzz_0 0_xz_0_yzz_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xz_1 0_xz_0_xzz_0 0_xz_0_xzz_1 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xx_1 0_xz_0_xxz_0 0_xz_0_xxz_1 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xz_0_xxx_1 0_xzz_0_zzz_0 0_xzz_0_yzz_0 0_xzz_0_yyz_0 0_xzz_0_yyy_0 0_xzz_0_xzz_0 0_xzz_0_xyz_0 0_xzz_0_xyy_0 0_xzz_0_xxz_0 0_xzz_0_xxy_0 0_xzz_0_xxx_0 0_xy_0_zzz_0 0_xy_0_yzz_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xzz_0 0_xy_0_xy_1 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xyy_1 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxy_1 0_xy_0_xxx_0 0_xyz_0_zzz_0 0_xyz_0_yzz_0 0_xyz_0_yyz_0 0_xyz_0_yyy_0 0_xyz_0_xzz_0 0_xyz_0_xyz_0 0_xyz_0_xyy_0 0_xyz_0_xxz_0 0_xyz_0_xxy_0 0_xyz_0_xxx_0 0_xyy_0_zzz_0 0_xyy_0_yzz_0 0_xyy_0_yyz_0 0_xyy_0_yyy_0 0_xyy_0_xzz_0 0_xyy_0_xyz_0 0_xyy_0_xyy_0 0_xyy_0_xxz_0 0_xyy_0_xxy_0 0_xyy_0_xxx_0 0_xx_0_zz_1 0_xx_0_zzz_0 0_xx_0_zzz_1 0_xx_0_yz_1 0_xx_0_yzz_0 0_xx_0_yzz_1 0_xx_0_yy_1 0_xx_0_yyz_0 0_xx_0_yyz_1 0_xx_0_yyy_0 0_xx_0_yyy_1 0_xx_0_xz_1 0_xx_0_xzz_0 0_xx_0_xzz_1 0_xx_0_xy_1 0_xx_0_xyz_0 0_xx_0_xyz_1 0_xx_0_xyy_0 0_xx_0_xyy_1 0_xx_0_xx_1 0_xx_0_xxz_0 0_xx_0_xxz_1 0_xx_0_xxy_0 0_xx_0_xxy_1 0_xx_0_xxx_0 0_xx_0_xxx_1 0_xxz_0_zzz_0 0_xxz_0_yzz_0 0_xxz_0_yyz_0 0_xxz_0_yyy_0 0_xxz_0_xzz_0 0_xxz_0_xyz_0 0_xxz_0_xyy_0 0_xxz_0_xxz_0 0_xxz_0_xxy_0 0_xxz_0_xxx_0 0_xxy_0_zzz_0 0_xxy_0_yzz_0 0_xxy_0_yyz_0 0_xxy_0_yyy_0 0_xxy_0_xzz_0 0_xxy_0_xyz_0 0_xxy_0_xyy_0 0_xxy_0_xxz_0 0_xxy_0_xxy_0 0_xxy_0_xxx_0 0_xxx_0_zzz_0 0_xxx_0_yzz_0 0_xxx_0_yyz_0 0_xxx_0_yyy_0 0_xxx_0_xzz_0 0_xxx_0_xyz_0 0_xxx_0_xyy_0 0_xxx_0_xxz_0 0_xxx_0_xxy_0 0_xxx_0_xxx_0 

signature:
0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_yy_0 0_zz_0_yy_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zz_0_xx_0 0_zz_0_xx_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_y_0_zz_0 0_y_0_zz_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_zz_1 0_yz_0_yz_0 0_yz_0_yz_1 0_yz_0_yy_0 0_yz_0_yy_1 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yy_0_z_1 0_yy_0_zz_0 0_yy_0_zz_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yy_0_xx_0 0_yy_0_xx_1 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_x_0_zz_0 0_x_0_zz_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xz_1 0_xz_0_xy_0 0_xz_0_xx_0 0_xz_0_xx_1 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xy_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xx_0_z_1 0_xx_0_zz_0 0_xx_0_zz_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_yy_0 0_xx_0_yy_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

signature:
0_zz_z_zz_0 0_zz_z_zzz_0 0_zz_z_yz_0 0_zz_z_yzz_0 0_zz_z_yy_0 0_zz_z_yyz_0 0_zz_z_xz_0 0_zz_z_xzz_0 0_zz_z_xy_0 0_zz_z_xyz_0 0_zz_z_xx_0 0_zz_z_xxz_0 0_zz_zz_zz_0 0_zz_zz_yz_0 0_zz_zz_yy_0 0_zz_zz_xz_0 0_zz_zz_xy_0 0_zz_zz_xx_0 0_zz_y_zz_0 0_zz_y_zzz_0 0_zz_y_yz_0 0_zz_y_yzz_0 0_zz_y_yy_0 0_zz_y_yyz_0 0_zz_y_yyy_0 0_zz_y_xz_0 0_zz_y_xzz_0 0_zz_y_xy_0 0_zz_y_xyz_0 0_zz_y_xyy_0 0_zz_y_xx_0 0_zz_y_xxz_0 0_zz_y_xxy_0 0_zz_yz_zz_0 0_zz_yz_yz_0 0_zz_yz_yy_0 0_zz_yz_xz_0 0_zz_yz_xy_0 0_zz_yz_xx_0 0_zz_yy_zz_0 0_zz_yy_yz_0 0_zz_yy_yy_0 0_zz_yy_xz_0 0_zz_yy_xy_0 0_zz_yy_xx_0 0_zz_x_zz_0 0_zz_x_zzz_0 0_zz_x_yz_0 0_zz_x_yzz_0 0_zz_x_yy_0 0_zz_x_yyz_0 0_zz_x_yyy_0 0_zz_x_xz_0 0_zz_x_xzz_0 0_zz_x_xy_0 0_zz_x_xyz_0 0_zz_x_xyy_0 0_zz_x_xx_0 0_zz_x_xxz_0 0_zz_x_xxy_0 0_zz_x_xxx_0 0_zz_xz_zz_0 0_zz_xz_yz_0 0_zz_xz_yy_0 0_zz_xz_xz_0 0_zz_xz_xy_0 0_zz_xz_xx_0 0_zz_xy_zz_0 0_zz_xy_yz_0 0_zz_xy_yy_0 0_zz_xy_xz_0 0_zz_xy_xy_0 0_zz_xy_xx_0 0_zz_xx_zz_0 0_zz_xx_yz_0 0_zz_xx_yy_0 0_zz_xx_xz_0 0_zz_xx_xy_0 0_zz_xx_xx_0 0_yz_z_zz_0 0_yz_z_zzz_0 0_yz_z_yz_0 0_yz_z_yzz_0 0_yz_z_yy_0 0_yz_z_yyz_0 0_yz_z_xz_0 0_yz_z_xzz_0 0_yz_z_xy_0 0_yz_z_xyz_0 0_yz_z_xx_0 0_yz_z_xxz_0 0_yz_zz_zz_0 0_yz_zz_yz_0 0_yz_zz_yy_0 0_yz_zz_xz_0 0_yz_zz_xy_0 0_yz_zz_xx_0 0_yz_y_zz_0 0_yz_y_zzz_0 0_yz_y_yz_0 0_yz_y_yzz_0 0_yz_y_yy_0 0_yz_y_yyz_0 0_yz_y_yyy_0 0_yz_y_xz_0 0_yz_y_xzz_0 0_yz_y_xy_0 0_yz_y_xyz_0 0_yz_y_xyy_0 0_yz_y_xx_0 0_yz_y_xxz_0 0_yz_y_xxy_0 0_yz_yz_zz_0 0_yz_yz_yz_0 0_yz_yz_yy_0 0_yz_yz_xz_0 0_yz_yz_xy_0 0_yz_yz_xx_0 0_yz_yy_zz_0 0_yz_yy_yz_0 0_yz_yy_yy_0 0_yz_yy_xz_0 0_yz_yy_xy_0 0_yz_yy_xx_0 0_yz_x_zz_0 0_yz_x_zzz_0 0_yz_x_yz_0 0_yz_x_yzz_0 0_yz_x_yy_0 0_yz_x_yyz_0 0_yz_x_yyy_0 0_yz_x_xz_0 0_yz_x_xzz_0 0_yz_x_xy_0 0_yz_x_xyz_0 0_yz_x_xyy_0 0_yz_x_xx_0 0_yz_x_xxz_0 0_yz_x_xxy_0 0_yz_x_xxx_0 0_yz_xz_zz_0 0_yz_xz_yz_0 0_yz_xz_yy_0 0_yz_xz_xz_0 0_yz_xz_xy_0 0_yz_xz_xx_0 0_yz_xy_zz_0 0_yz_xy_yz_0 0_yz_xy_yy_0 0_yz_xy_xz_0 0_yz_xy_xy_0 0_yz_xy_xx_0 0_yz_xx_zz_0 0_yz_xx_yz_0 0_yz_xx_yy_0 0_yz_xx_xz_0 0_yz_xx_xy_0 0_yz_xx_xx_0 0_yy_z_zz_0 0_yy_z_zzz_0 0_yy_z_yz_0 0_yy_z_yzz_0 0_yy_z_yy_0 0_yy_z_yyz_0 0_yy_z_xz_0 0_yy_z_xzz_0 0_yy_z_xy_0 0_yy_z_xyz_0 0_yy_z_xx_0 0_yy_z_xxz_0 0_yy_zz_zz_0 0_yy_zz_yz_0 0_yy_zz_yy_0 0_yy_zz_xz_0 0_yy_zz_xy_0 0_yy_zz_xx_0 0_yy_y_zz_0 0_yy_y_zzz_0 0_yy_y_yz_0 0_yy_y_yzz_0 0_yy_y_yy_0 0_yy_y_yyz_0 0_yy_y_yyy_0 0_yy_y_xz_0 0_yy_y_xzz_0 0_yy_y_xy_0 0_yy_y_xyz_0 0_yy_y_xyy_0 0_yy_y_xx_0 0_yy_y_xxz_0 0_yy_y_xxy_0 0_yy_yz_zz_0 0_yy_yz_yz_0 0_yy_yz_yy_0 0_yy_yz_xz_0 0_yy_yz_xy_0 0_yy_yz_xx_0 0_yy_yy_zz_0 0_yy_yy_yz_0 0_yy_yy_yy_0 0_yy_yy_xz_0 0_yy_yy_xy_0 0_yy_yy_xx_0 0_yy_x_zz_0 0_yy_x_zzz_0 0_yy_x_yz_0 0_yy_x_yzz_0 0_yy_x_yy_0 0_yy_x_yyz_0 0_yy_x_yyy_0 0_yy_x_xz_0 0_yy_x_xzz_0 0_yy_x_xy_0 0_yy_x_xyz_0 0_yy_x_xyy_0 0_yy_x_xx_0 0_yy_x_xxz_0 0_yy_x_xxy_0 0_yy_x_xxx_0 0_yy_xz_zz_0 0_yy_xz_yz_0 0_yy_xz_yy_0 0_yy_xz_xz_0 0_yy_xz_xy_0 0_yy_xz_xx_0 0_yy_xy_zz_0 0_yy_xy_yz_0 0_yy_xy_yy_0 0_yy_xy_xz_0 0_yy_xy_xy_0 0_yy_xy_xx_0 0_yy_xx_zz_0 0_yy_xx_yz_0 0_yy_xx_yy_0 0_yy_xx_xz_0 0_yy_xx_xy_0 0_yy_xx_xx_0 0_xz_z_zz_0 0_xz_z_zzz_0 0_xz_z_yz_0 0_xz_z_yzz_0 0_xz_z_yy_0 0_xz_z_yyz_0 0_xz_z_xz_0 0_xz_z_xzz_0 0_xz_z_xy_0 0_xz_z_xyz_0 0_xz_z_xx_0 0_xz_z_xxz_0 0_xz_zz_zz_0 0_xz_zz_yz_0 0_xz_zz_yy_0 0_xz_zz_xz_0 0_xz_zz_xy_0 0_xz_zz_xx_0 0_xz_y_zz_0 0_xz_y_zzz_0 0_xz_y_yz_0 0_xz_y_yzz_0 0_xz_y_yy_0 0_xz_y_yyz_0 0_xz_y_yyy_0 0_xz_y_xz_0 0_xz_y_xzz_0 0_xz_y_xy_0 0_xz_y_xyz_0 0_xz_y_xyy_0 0_xz_y_xx_0 0_xz_y_xxz_0 0_xz_y_xxy_0 0_xz_yz_zz_0 0_xz_yz_yz_0 0_xz_yz_yy_0 0_xz_yz_xz_0 0_xz_yz_xy_0 0_xz_yz_xx_0 0_xz_yy_zz_0 0_xz_yy_yz_0 0_xz_yy_yy_0 0_xz_yy_xz_0 0_xz_yy_xy_0 0_xz_yy_xx_0 0_xz_x_zz_0 0_xz_x_zzz_0 0_xz_x_yz_0 0_xz_x_yzz_0 0_xz_x_yy_0 0_xz_x_yyz_0 0_xz_x_yyy_0 0_xz_x_xz_0 0_xz_x_xzz_0 0_xz_x_xy_0 0_xz_x_xyz_0 0_xz_x_xyy_0 0_xz_x_xx_0 0_xz_x_xxz_0 0_xz_x_xxy_0 0_xz_x_xxx_0 0_xz_xz_zz_0 0_xz_xz_yz_0 0_xz_xz_yy_0 0_xz_xz_xz_0 0_xz_xz_xy_0 0_xz_xz_xx_0 0_xz_xy_zz_0 0_xz_xy_yz_0 0_xz_xy_yy_0 0_xz_xy_xz_0 0_xz_xy_xy_0 0_xz_xy_xx_0 0_xz_xx_zz_0 0_xz_xx_yz_0 0_xz_xx_yy_0 0_xz_xx_xz_0 0_xz_xx_xy_0 0_xz_xx_xx_0 0_xy_z_zz_0 0_xy_z_zzz_0 0_xy_z_yz_0 0_xy_z_yzz_0 0_xy_z_yy_0 0_xy_z_yyz_0 0_xy_z_xz_0 0_xy_z_xzz_0 0_xy_z_xy_0 0_xy_z_xyz_0 0_xy_z_xx_0 0_xy_z_xxz_0 0_xy_zz_zz_0 0_xy_zz_yz_0 0_xy_zz_yy_0 0_xy_zz_xz_0 0_xy_zz_xy_0 0_xy_zz_xx_0 0_xy_y_zz_0 0_xy_y_zzz_0 0_xy_y_yz_0 0_xy_y_yzz_0 0_xy_y_yy_0 0_xy_y_yyz_0 0_xy_y_yyy_0 0_xy_y_xz_0 0_xy_y_xzz_0 0_xy_y_xy_0 0_xy_y_xyz_0 0_xy_y_xyy_0 0_xy_y_xx_0 0_xy_y_xxz_0 0_xy_y_xxy_0 0_xy_yz_zz_0 0_xy_yz_yz_0 0_xy_yz_yy_0 0_xy_yz_xz_0 0_xy_yz_xy_0 0_xy_yz_xx_0 0_xy_yy_zz_0 0_xy_yy_yz_0 0_xy_yy_yy_0 0_xy_yy_xz_0 0_xy_yy_xy_0 0_xy_yy_xx_0 0_xy_x_zz_0 0_xy_x_zzz_0 0_xy_x_yz_0 0_xy_x_yzz_0 0_xy_x_yy_0 0_xy_x_yyz_0 0_xy_x_yyy_0 0_xy_x_xz_0 0_xy_x_xzz_0 0_xy_x_xy_0 0_xy_x_xyz_0 0_xy_x_xyy_0 0_xy_x_xx_0 0_xy_x_xxz_0 0_xy_x_xxy_0 0_xy_x_xxx_0 0_xy_xz_zz_0 0_xy_xz_yz_0 0_xy_xz_yy_0 0_xy_xz_xz_0 0_xy_xz_xy_0 0_xy_xz_xx_0 0_xy_xy_zz_0 0_xy_xy_yz_0 0_xy_xy_yy_0 0_xy_xy_xz_0 0_xy_xy_xy_0 0_xy_xy_xx_0 0_xy_xx_zz_0 0_xy_xx_yz_0 0_xy_xx_yy_0 0_xy_xx_xz_0 0_xy_xx_xy_0 0_xy_xx_xx_0 0_xx_z_zz_0 0_xx_z_zzz_0 0_xx_z_yz_0 0_xx_z_yzz_0 0_xx_z_yy_0 0_xx_z_yyz_0 0_xx_z_xz_0 0_xx_z_xzz_0 0_xx_z_xy_0 0_xx_z_xyz_0 0_xx_z_xx_0 0_xx_z_xxz_0 0_xx_zz_zz_0 0_xx_zz_yz_0 0_xx_zz_yy_0 0_xx_zz_xz_0 0_xx_zz_xy_0 0_xx_zz_xx_0 0_xx_y_zz_0 0_xx_y_zzz_0 0_xx_y_yz_0 0_xx_y_yzz_0 0_xx_y_yy_0 0_xx_y_yyz_0 0_xx_y_yyy_0 0_xx_y_xz_0 0_xx_y_xzz_0 0_xx_y_xy_0 0_xx_y_xyz_0 0_xx_y_xyy_0 0_xx_y_xx_0 0_xx_y_xxz_0 0_xx_y_xxy_0 0_xx_yz_zz_0 0_xx_yz_yz_0 0_xx_yz_yy_0 0_xx_yz_xz_0 0_xx_yz_xy_0 0_xx_yz_xx_0 0_xx_yy_zz_0 0_xx_yy_yz_0 0_xx_yy_yy_0 0_xx_yy_xz_0 0_xx_yy_xy_0 0_xx_yy_xx_0 0_xx_x_zz_0 0_xx_x_zzz_0 0_xx_x_yz_0 0_xx_x_yzz_0 0_xx_x_yy_0 0_xx_x_yyz_0 0_xx_x_yyy_0 0_xx_x_xz_0 0_xx_x_xzz_0 0_xx_x_xy_0 0_xx_x_xyz_0 0_xx_x_xyy_0 0_xx_x_xx_0 0_xx_x_xxz_0 0_xx_x_xxy_0 0_xx_x_xxx_0 0_xx_xz_zz_0 0_xx_xz_yz_0 0_xx_xz_yy_0 0_xx_xz_xz_0 0_xx_xz_xy_0 0_xx_xz_xx_0 0_xx_xy_zz_0 0_xx_xy_yz_0 0_xx_xy_yy_0 0_xx_xy_xz_0 0_xx_xy_xy_0 0_xx_xy_xx_0 0_xx_xx_zz_0 0_xx_xx_yz_0 0_xx_xx_yy_0 0_xx_xx_xz_0 0_xx_xx_xy_0 0_xx_xx_xx_0 

signature:
0_zz_0_zzz_0 0_zz_0_zzzz_0 0_zz_0_yzz_0 0_zz_0_yzzz_0 0_zz_0_yyz_0 0_zz_0_yyzz_0 0_zz_0_yyy_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzz_0 0_zz_0_xzzz_0 0_zz_0_xyz_0 0_zz_0_xyzz_0 0_zz_0_xyy_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxz_0 0_zz_0_xxzz_0 0_zz_0_xxy_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxx_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_zz_z_zzz_0 0_zz_z_yzz_0 0_zz_z_yyz_0 0_zz_z_xzz_0 0_zz_z_xyz_0 0_zz_z_xxz_0 0_zz_y_zzz_0 0_zz_y_yzz_0 0_zz_y_yyz_0 0_zz_y_yyy_0 0_zz_y_xzz_0 0_zz_y_xyz_0 0_zz_y_xyy_0 0_zz_y_xxz_0 0_zz_y_xxy_0 0_zz_x_zzz_0 0_zz_x_yzz_0 0_zz_x_yyz_0 0_zz_x_yyy_0 0_zz_x_xzz_0 0_zz_x_xyz_0 0_zz_x_xyy_0 0_zz_x_xxz_0 0_zz_x_xxy_0 0_zz_x_xxx_0 0_yz_0_zzz_0 0_yz_0_zzzz_0 0_yz_0_yzz_0 0_yz_0_yzzz_0 0_yz_0_yyz_0 0_yz_0_yyzz_0 0_yz_0_yyy_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzz_0 0_yz_0_xzzz_0 0_yz_0_xyz_0 0_yz_0_xyzz_0 0_yz_0_xyy_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxz_0 0_yz_0_xxzz_0 0_yz_0_xxy_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxx_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yz_z_zzz_0 0_yz_z_yzz_0 0_yz_z_yyz_0 0_yz_z_xzz_0 0_yz_z_xyz_0 0_yz_z_xxz_0 0_yz_y_zzz_0 0_yz_y_yzz_0 0_yz_y_yyz_0 0_yz_y_yyy_0 0_yz_y_xzz_0 0_yz_y_xyz_0 0_yz_y_xyy_0 0_yz_y_xxz_0 0_yz_y_xxy_0 0_yz_x_zzz_0 0_yz_x_yzz_0 0_yz_x_yyz_0 0_yz_x_yyy_0 0_yz_x_xzz_0 0_yz_x_xyz_0 0_yz_x_xyy_0 0_yz_x_xxz_0 0_yz_x_xxy_0 0_yz_x_xxx_0 0_yy_0_zzz_0 0_yy_0_zzzz_0 0_yy_0_yzz_0 0_yy_0_yzzz_0 0_yy_0_yyz_0 0_yy_0_yyzz_0 0_yy_0_yyy_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzz_0 0_yy_0_xzzz_0 0_yy_0_xyz_0 0_yy_0_xyzz_0 0_yy_0_xyy_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxz_0 0_yy_0_xxzz_0 0_yy_0_xxy_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxx_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_yy_z_zzz_0 0_yy_z_yzz_0 0_yy_z_yyz_0 0_yy_z_xzz_0 0_yy_z_xyz_0 0_yy_z_xxz_0 0_yy_y_zzz_0 0_yy_y_yzz_0 0_yy_y_yyz_0 0_yy_y_yyy_0 0_yy_y_xzz_0 0_yy_y_xyz_0 0_yy_y_xyy_0 0_yy_y_xxz_0 0_yy_y_xxy_0 0_yy_x_zzz_0 0_yy_x_yzz_0 0_yy_x_yyz_0 0_yy_x_yyy_0 0_yy_x_xzz_0 0_yy_x_xyz_0 0_yy_x_xyy_0 0_yy_x_xxz_0 0_yy_x_xxy_0 0_yy_x_xxx_0 0_xz_0_zzz_0 0_xz_0_zzzz_0 0_xz_0_yzz_0 0_xz_0_yzzz_0 0_xz_0_yyz_0 0_xz_0_yyzz_0 0_xz_0_yyy_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzz_0 0_xz_0_xzzz_0 0_xz_0_xyz_0 0_xz_0_xyzz_0 0_xz_0_xyy_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxz_0 0_xz_0_xxzz_0 0_xz_0_xxy_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxx_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xz_z_zzz_0 0_xz_z_yzz_0 0_xz_z_yyz_0 0_xz_z_xzz_0 0_xz_z_xyz_0 0_xz_z_xxz_0 0_xz_y_zzz_0 0_xz_y_yzz_0 0_xz_y_yyz_0 0_xz_y_yyy_0 0_xz_y_xzz_0 0_xz_y_xyz_0 0_xz_y_xyy_0 0_xz_y_xxz_0 0_xz_y_xxy_0 0_xz_x_zzz_0 0_xz_x_yzz_0 0_xz_x_yyz_0 0_xz_x_yyy_0 0_xz_x_xzz_0 0_xz_x_xyz_0 0_xz_x_xyy_0 0_xz_x_xxz_0 0_xz_x_xxy_0 0_xz_x_xxx_0 0_xy_0_zzz_0 0_xy_0_zzzz_0 0_xy_0_yzz_0 0_xy_0_yzzz_0 0_xy_0_yyz_0 0_xy_0_yyzz_0 0_xy_0_yyy_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzz_0 0_xy_0_xzzz_0 0_xy_0_xyz_0 0_xy_0_xyzz_0 0_xy_0_xyy_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxz_0 0_xy_0_xxzz_0 0_xy_0_xxy_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxx_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xy_z_zzz_0 0_xy_z_yzz_0 0_xy_z_yyz_0 0_xy_z_xzz_0 0_xy_z_xyz_0 0_xy_z_xxz_0 0_xy_y_zzz_0 0_xy_y_yzz_0 0_xy_y_yyz_0 0_xy_y_yyy_0 0_xy_y_xzz_0 0_xy_y_xyz_0 0_xy_y_xyy_0 0_xy_y_xxz_0 0_xy_y_xxy_0 0_xy_x_zzz_0 0_xy_x_yzz_0 0_xy_x_yyz_0 0_xy_x_yyy_0 0_xy_x_xzz_0 0_xy_x_xyz_0 0_xy_x_xyy_0 0_xy_x_xxz_0 0_xy_x_xxy_0 0_xy_x_xxx_0 0_xx_0_zzz_0 0_xx_0_zzzz_0 0_xx_0_yzz_0 0_xx_0_yzzz_0 0_xx_0_yyz_0 0_xx_0_yyzz_0 0_xx_0_yyy_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzz_0 0_xx_0_xzzz_0 0_xx_0_xyz_0 0_xx_0_xyzz_0 0_xx_0_xyy_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxz_0 0_xx_0_xxzz_0 0_xx_0_xxy_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxx_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 0_xx_z_zzz_0 0_xx_z_yzz_0 0_xx_z_yyz_0 0_xx_z_xzz_0 0_xx_z_xyz_0 0_xx_z_xxz_0 0_xx_y_zzz_0 0_xx_y_yzz_0 0_xx_y_yyz_0 0_xx_y_yyy_0 0_xx_y_xzz_0 0_xx_y_xyz_0 0_xx_y_xyy_0 0_xx_y_xxz_0 0_xx_y_xxy_0 0_xx_x_zzz_0 0_xx_x_yzz_0 0_xx_x_yyz_0 0_xx_x_yyy_0 0_xx_x_xzz_0 0_xx_x_xyz_0 0_xx_x_xyy_0 0_xx_x_xxz_0 0_xx_x_xxy_0 0_xx_x_xxx_0 

signature:
0_zz_0_zz_0 0_zz_0_zzz_0 0_zz_0_yz_0 0_zz_0_yzz_0 0_zz_0_yy_0 0_zz_0_yyz_0 0_zz_0_yyy_0 0_zz_0_xz_0 0_zz_0_xzz_0 0_zz_0_xy_0 0_zz_0_xyz_0 0_zz_0_xyy_0 0_zz_0_xx_0 0_zz_0_xxz_0 0_zz_0_xxy_0 0_zz_0_xxx_0 0_zz_z_zz_0 0_zz_z_yz_0 0_zz_z_yy_0 0_zz_z_xz_0 0_zz_z_xy_0 0_zz_z_xx_0 0_zz_y_zz_0 0_zz_y_yz_0 0_zz_y_yy_0 0_zz_y_xz_0 0_zz_y_xy_0 0_zz_y_xx_0 0_zz_x_zz_0 0_zz_x_yz_0 0_zz_x_yy_0 0_zz_x_xz_0 0_zz_x_xy_0 0_zz_x_xx_0 0_yz_0_zz_0 0_yz_0_zzz_0 0_yz_0_yz_0 0_yz_0_yzz_0 0_yz_0_yy_0 0_yz_0_yyz_0 0_yz_0_yyy_0 0_yz_0_xz_0 0_yz_0_xzz_0 0_yz_0_xy_0 0_yz_0_xyz_0 0_yz_0_xyy_0 0_yz_0_xx_0 0_yz_0_xxz_0 0_yz_0_xxy_0 0_yz_0_xxx_0 0_yz_z_zz_0 0_yz_z_yz_0 0_yz_z_yy_0 0_yz_z_xz_0 0_yz_z_xy_0 0_yz_z_xx_0 0_yz_y_zz_0 0_yz_y_yz_0 0_yz_y_yy_0 0_yz_y_xz_0 0_yz_y_xy_0 0_yz_y_xx_0 0_yz_x_zz_0 0_yz_x_yz_0 0_yz_x_yy_0 0_yz_x_xz_0 0_yz_x_xy_0 0_yz_x_xx_0 0_yy_0_zz_0 0_yy_0_zzz_0 0_yy_0_yz_0 0_yy_0_yzz_0 0_yy_0_yy_0 0_yy_0_yyz_0 0_yy_0_yyy_0 0_yy_0_xz_0 0_yy_0_xzz_0 0_yy_0_xy_0 0_yy_0_xyz_0 0_yy_0_xyy_0 0_yy_0_xx_0 0_yy_0_xxz_0 0_yy_0_xxy_0 0_yy_0_xxx_0 0_yy_z_zz_0 0_yy_z_yz_0 0_yy_z_yy_0 0_yy_z_xz_0 0_yy_z_xy_0 0_yy_z_xx_0 0_yy_y_zz_0 0_yy_y_yz_0 0_yy_y_yy_0 0_yy_y_xz_0 0_yy_y_xy_0 0_yy_y_xx_0 0_yy_x_zz_0 0_yy_x_yz_0 0_yy_x_yy_0 0_yy_x_xz_0 0_yy_x_xy_0 0_yy_x_xx_0 0_xz_0_zz_0 0_xz_0_zzz_0 0_xz_0_yz_0 0_xz_0_yzz_0 0_xz_0_yy_0 0_xz_0_yyz_0 0_xz_0_yyy_0 0_xz_0_xz_0 0_xz_0_xzz_0 0_xz_0_xy_0 0_xz_0_xyz_0 0_xz_0_xyy_0 0_xz_0_xx_0 0_xz_0_xxz_0 0_xz_0_xxy_0 0_xz_0_xxx_0 0_xz_z_zz_0 0_xz_z_yz_0 0_xz_z_yy_0 0_xz_z_xz_0 0_xz_z_xy_0 0_xz_z_xx_0 0_xz_y_zz_0 0_xz_y_yz_0 0_xz_y_yy_0 0_xz_y_xz_0 0_xz_y_xy_0 0_xz_y_xx_0 0_xz_x_zz_0 0_xz_x_yz_0 0_xz_x_yy_0 0_xz_x_xz_0 0_xz_x_xy_0 0_xz_x_xx_0 0_xy_0_zz_0 0_xy_0_zzz_0 0_xy_0_yz_0 0_xy_0_yzz_0 0_xy_0_yy_0 0_xy_0_yyz_0 0_xy_0_yyy_0 0_xy_0_xz_0 0_xy_0_xzz_0 0_xy_0_xy_0 0_xy_0_xyz_0 0_xy_0_xyy_0 0_xy_0_xx_0 0_xy_0_xxz_0 0_xy_0_xxy_0 0_xy_0_xxx_0 0_xy_z_zz_0 0_xy_z_yz_0 0_xy_z_yy_0 0_xy_z_xz_0 0_xy_z_xy_0 0_xy_z_xx_0 0_xy_y_zz_0 0_xy_y_yz_0 0_xy_y_yy_0 0_xy_y_xz_0 0_xy_y_xy_0 0_xy_y_xx_0 0_xy_x_zz_0 0_xy_x_yz_0 0_xy_x_yy_0 0_xy_x_xz_0 0_xy_x_xy_0 0_xy_x_xx_0 0_xx_0_zz_0 0_xx_0_zzz_0 0_xx_0_yz_0 0_xx_0_yzz_0 0_xx_0_yy_0 0_xx_0_yyz_0 0_xx_0_yyy_0 0_xx_0_xz_0 0_xx_0_xzz_0 0_xx_0_xy_0 0_xx_0_xyz_0 0_xx_0_xyy_0 0_xx_0_xx_0 0_xx_0_xxz_0 0_xx_0_xxy_0 0_xx_0_xxx_0 0_xx_z_zz_0 0_xx_z_yz_0 0_xx_z_yy_0 0_xx_z_xz_0 0_xx_z_xy_0 0_xx_z_xx_0 0_xx_y_zz_0 0_xx_y_yz_0 0_xx_y_yy_0 0_xx_y_xz_0 0_xx_y_xy_0 0_xx_y_xx_0 0_xx_x_zz_0 0_xx_x_yz_0 0_xx_x_yy_0 0_xx_x_xz_0 0_xx_x_xy_0 0_xx_x_xx_0 

signature:
0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyy_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_yyyy_0 0_z_0_yyyy_1 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyy_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xyyy_0 0_z_0_xyyy_1 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxy_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxyy_0 0_z_0_xxyy_1 0_z_0_xxx_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_z_0_xxxy_0 0_z_0_xxxy_1 0_z_0_xxxx_0 0_z_0_xxxx_1 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_y_0_zzz_1 0_y_0_zzzz_0 0_y_0_zzzz_1 0_y_0_yzz_1 0_y_0_yzzz_0 0_y_0_yzzz_1 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xzz_1 0_y_0_xzzz_0 0_y_0_xzzz_1 0_y_0_xyz_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxz_1 0_y_0_xxzz_0 0_y_0_xxzz_1 0_y_0_xxy_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxx_1 0_y_0_xxxz_0 0_y_0_xxxz_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_y_0_xxxx_0 0_y_0_xxxx_1 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xxyz_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_x_0_zzz_1 0_x_0_zzzz_0 0_x_0_zzzz_1 0_x_0_yzz_1 0_x_0_yzzz_0 0_x_0_yzzz_1 0_x_0_yyz_1 0_x_0_yyzz_0 0_x_0_yyzz_1 0_x_0_yyy_1 0_x_0_yyyz_0 0_x_0_yyyz_1 0_x_0_yyyy_0 0_x_0_yyyy_1 0_x_0_xzz_1 0_x_0_xzzz_0 0_x_0_xzzz_1 0_x_0_xyz_1 0_x_0_xyzz_0 0_x_0_xyzz_1 0_x_0_xyy_1 0_x_0_xyyz_0 0_x_0_xyyz_1 0_x_0_xyyy_0 0_x_0_xyyy_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_xzzz_0 0_xz_0_xxzz_0 0_xz_0_xxxz_0 0_xz_0_xxxx_0 0_xy_0_xyyy_0 0_xy_0_xxyy_0 0_xy_0_xxxy_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

signature:
0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyy_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_yyyy_0 0_z_0_yyyy_1 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyy_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xyyy_0 0_z_0_xyyy_1 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxy_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxyy_0 0_z_0_xxyy_1 0_z_0_xxx_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_z_0_xxxy_0 0_z_0_xxxy_1 0_z_0_xxxx_0 0_z_0_xxxx_1 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_y_0_zzz_1 0_y_0_zzzz_0 0_y_0_zzzz_1 0_y_0_yzz_1 0_y_0_yzzz_0 0_y_0_yzzz_1 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xzz_1 0_y_0_xzzz_0 0_y_0_xzzz_1 0_y_0_xyz_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxz_1 0_y_0_xxzz_0 0_y_0_xxzz_1 0_y_0_xxy_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxx_1 0_y_0_xxxz_0 0_y_0_xxxz_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_y_0_xxxx_0 0_y_0_xxxx_1 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzzz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_x_0_zzz_1 0_x_0_zzzz_0 0_x_0_zzzz_1 0_x_0_yzz_1 0_x_0_yzzz_0 0_x_0_yzzz_1 0_x_0_yyz_1 0_x_0_yyzz_0 0_x_0_yyzz_1 0_x_0_yyy_1 0_x_0_yyyz_0 0_x_0_yyyz_1 0_x_0_yyyy_0 0_x_0_yyyy_1 0_x_0_xzz_1 0_x_0_xzzz_0 0_x_0_xzzz_1 0_x_0_xyz_1 0_x_0_xyzz_0 0_x_0_xyzz_1 0_x_0_xyy_1 0_x_0_xyyz_0 0_x_0_xyyz_1 0_x_0_xyyy_0 0_x_0_xyyy_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxzz_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

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
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

signature:
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

signature:
0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

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
0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xz_0 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xy_0 0_0_0_xyz_0 0_0_0_xyzz_0 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxy_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxx_0 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

signature:
0_0_0_zz_0 0_0_0_zz_1 0_0_0_zzz_0 0_0_0_zzz_1 0_0_0_zzzz_0 0_0_0_yz_0 0_0_0_yz_1 0_0_0_yzz_0 0_0_0_yzz_1 0_0_0_yzzz_0 0_0_0_yy_0 0_0_0_yy_1 0_0_0_yyz_0 0_0_0_yyz_1 0_0_0_yyzz_0 0_0_0_yyy_0 0_0_0_yyy_1 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xz_0 0_0_0_xz_1 0_0_0_xzz_0 0_0_0_xzz_1 0_0_0_xzzz_0 0_0_0_xy_0 0_0_0_xy_1 0_0_0_xyz_0 0_0_0_xyz_1 0_0_0_xyzz_0 0_0_0_xyy_0 0_0_0_xyy_1 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xx_0 0_0_0_xx_1 0_0_0_xxz_0 0_0_0_xxz_1 0_0_0_xxzz_0 0_0_0_xxy_0 0_0_0_xxy_1 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxx_0 0_0_0_xxx_1 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

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

INTEGRAL:0 : 1 : 6 : SSSS_6 Y SSSS_7 Y SSSP_6 Y 

0_0_0_0_6 0_0_0_0_7 0_0_0_z_6 0_0_0_y_6 0_0_0_x_6 

SSSS_6 SSSS_7 SSSP_6 

SSSS_6 : 0_0_0_0_6 

SSSS_7 : 0_0_0_0_7 

SSSP_6 : 0_0_0_z_6 0_0_0_y_6 0_0_0_x_6 

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

INTEGRAL:0 : 2 : 5 : SSSS_5 Y SSSS_6 Y SSSP_5 Y SSSP_6 Y SSSD_5 Y 

0_0_0_0_5 0_0_0_0_6 0_0_0_z_5 0_0_0_z_6 0_0_0_zz_5 0_0_0_y_5 0_0_0_y_6 0_0_0_yy_5 0_0_0_x_5 0_0_0_x_6 0_0_0_xx_5 

SSSS_5 SSSS_6 SSSP_5 SSSP_6 SSSD_5 

SSSS_5 : 0_0_0_0_5 

SSSS_6 : 0_0_0_0_6 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSP_6 : 0_0_0_z_6 0_0_0_y_6 0_0_0_x_6 

SSSD_5 : 0_0_0_zz_5 0_0_0_yy_5 0_0_0_xx_5 

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

INTEGRAL:0 : 3 : 4 : SSSP_4 Y SSSP_5 Y SSSD_4 N SSSD_5 Y SSSF_4 Y 

0_0_0_z_4 0_0_0_z_5 0_0_0_zz_4 0_0_0_zz_5 0_0_0_zzz_4 0_0_0_y_4 0_0_0_y_5 0_0_0_yzz_4 0_0_0_yy_4 0_0_0_yy_5 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_x_4 0_0_0_x_5 0_0_0_xzz_4 0_0_0_xyy_4 0_0_0_xx_4 0_0_0_xx_5 0_0_0_xxz_4 0_0_0_xxx_4 

SSSP_4 SSSP_5 SSSD_4 SSSD_5 SSSF_4 

SSSP_4 : 0_0_0_z_4 0_0_0_y_4 0_0_0_x_4 

SSSP_5 : 0_0_0_z_5 0_0_0_y_5 0_0_0_x_5 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xx_4 

SSSD_5 : 0_0_0_zz_5 0_0_0_yy_5 0_0_0_xx_5 

SSSF_4 : 0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxx_4 

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

INTEGRAL:0 : 4 : 2 : SSSD_2 N SSSD_3 N SSSF_2 N SSSF_3 N SSSG_2 Y 

0_0_0_zz_2 0_0_0_zz_3 0_0_0_zzz_2 0_0_0_zzz_3 0_0_0_zzzz_2 0_0_0_yzz_2 0_0_0_yzz_3 0_0_0_yzzz_2 0_0_0_yy_2 0_0_0_yy_3 0_0_0_yyz_2 0_0_0_yyz_3 0_0_0_yyzz_2 0_0_0_yyy_2 0_0_0_yyy_3 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzz_2 0_0_0_xzz_3 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyy_2 0_0_0_xyy_3 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xx_2 0_0_0_xx_3 0_0_0_xxz_2 0_0_0_xxz_3 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxx_2 0_0_0_xxx_3 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SSSD_2 SSSD_3 SSSF_2 SSSF_3 SSSG_2 

SSSD_2 : 0_0_0_zz_2 0_0_0_yz_2 0_0_0_yy_2 0_0_0_xz_2 0_0_0_xy_2 0_0_0_xx_2 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSF_2 : 0_0_0_zzz_2 0_0_0_yzz_2 0_0_0_yyz_2 0_0_0_yyy_2 0_0_0_xzz_2 0_0_0_xyz_2 0_0_0_xyy_2 0_0_0_xxz_2 0_0_0_xxy_2 0_0_0_xxx_2 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

INTEGRAL:0 : 4 : 3 : SSSD_3 N SSSD_4 N SSSF_3 N SSSF_4 Y SSSG_3 Y 

0_0_0_zz_3 0_0_0_zz_4 0_0_0_zzz_3 0_0_0_zzz_4 0_0_0_zzzz_3 0_0_0_yzz_3 0_0_0_yzz_4 0_0_0_yzzz_3 0_0_0_yy_3 0_0_0_yy_4 0_0_0_yyz_3 0_0_0_yyz_4 0_0_0_yyzz_3 0_0_0_yyy_3 0_0_0_yyy_4 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzz_3 0_0_0_xzz_4 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyy_3 0_0_0_xyy_4 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xx_3 0_0_0_xx_4 0_0_0_xxz_3 0_0_0_xxz_4 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxx_3 0_0_0_xxx_4 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 

SSSD_3 SSSD_4 SSSF_3 SSSF_4 SSSG_3 

SSSD_3 : 0_0_0_zz_3 0_0_0_yz_3 0_0_0_yy_3 0_0_0_xz_3 0_0_0_xy_3 0_0_0_xx_3 

SSSD_4 : 0_0_0_zz_4 0_0_0_yz_4 0_0_0_yy_4 0_0_0_xx_4 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SSSF_4 : 0_0_0_zzz_4 0_0_0_yzz_4 0_0_0_yyz_4 0_0_0_yyy_4 0_0_0_xzz_4 0_0_0_xyy_4 0_0_0_xxz_4 0_0_0_xxx_4 

SSSG_3 : 0_0_0_zzzz_3 0_0_0_yzzz_3 0_0_0_yyzz_3 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 

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

INTEGRAL:1 : 4 : 2 : SSSF_3 Y SSSG_2 Y SSSG_3 Y SPSG_2 Y 

0_0_0_zzz_3 0_0_0_zzzz_2 0_0_0_zzzz_3 0_0_0_yzz_3 0_0_0_yzzz_2 0_0_0_yzzz_3 0_0_0_yyz_3 0_0_0_yyzz_2 0_0_0_yyzz_3 0_0_0_yyy_3 0_0_0_yyyz_2 0_0_0_yyyz_3 0_0_0_yyyy_2 0_0_0_yyyy_3 0_0_0_xzz_3 0_0_0_xzzz_2 0_0_0_xzzz_3 0_0_0_xyz_3 0_0_0_xyzz_2 0_0_0_xyzz_3 0_0_0_xyy_3 0_0_0_xyyz_2 0_0_0_xyyz_3 0_0_0_xyyy_2 0_0_0_xyyy_3 0_0_0_xxz_3 0_0_0_xxzz_2 0_0_0_xxzz_3 0_0_0_xxy_3 0_0_0_xxyz_2 0_0_0_xxyz_3 0_0_0_xxyy_2 0_0_0_xxyy_3 0_0_0_xxx_3 0_0_0_xxxz_2 0_0_0_xxxz_3 0_0_0_xxxy_2 0_0_0_xxxy_3 0_0_0_xxxx_2 0_0_0_xxxx_3 0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_yyyz_2 0_z_0_yyyy_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xyyy_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_z_0_xxyy_2 0_z_0_xxxz_2 0_z_0_xxxy_2 0_z_0_xxxx_2 0_y_0_zzzz_2 0_y_0_yzzz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xzzz_2 0_y_0_xyzz_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxzz_2 0_y_0_xxyz_2 0_y_0_xxyy_2 0_y_0_xxxz_2 0_y_0_xxxy_2 0_y_0_xxxx_2 0_x_0_zzzz_2 0_x_0_yzzz_2 0_x_0_yyzz_2 0_x_0_yyyz_2 0_x_0_yyyy_2 0_x_0_xzzz_2 0_x_0_xyzz_2 0_x_0_xyyz_2 0_x_0_xyyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 

SSSF_3 SSSG_2 SSSG_3 SPSG_2 

SSSF_3 : 0_0_0_zzz_3 0_0_0_yzz_3 0_0_0_yyz_3 0_0_0_yyy_3 0_0_0_xzz_3 0_0_0_xyz_3 0_0_0_xyy_3 0_0_0_xxz_3 0_0_0_xxy_3 0_0_0_xxx_3 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SSSG_3 : 0_0_0_zzzz_3 0_0_0_yzzz_3 0_0_0_yyzz_3 0_0_0_yyyz_3 0_0_0_yyyy_3 0_0_0_xzzz_3 0_0_0_xyzz_3 0_0_0_xyyz_3 0_0_0_xyyy_3 0_0_0_xxzz_3 0_0_0_xxyz_3 0_0_0_xxyy_3 0_0_0_xxxz_3 0_0_0_xxxy_3 0_0_0_xxxx_3 

SPSG_2 : 0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_yyyz_2 0_z_0_yyyy_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xyyy_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_z_0_xxyy_2 0_z_0_xxxz_2 0_z_0_xxxy_2 0_z_0_xxxx_2 0_y_0_zzzz_2 0_y_0_yzzz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xzzz_2 0_y_0_xyzz_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxzz_2 0_y_0_xxyz_2 0_y_0_xxyy_2 0_y_0_xxxz_2 0_y_0_xxxy_2 0_y_0_xxxx_2 0_x_0_zzzz_2 0_x_0_yzzz_2 0_x_0_yyzz_2 0_x_0_yyyz_2 0_x_0_yyyy_2 0_x_0_xzzz_2 0_x_0_xyzz_2 0_x_0_xyyz_2 0_x_0_xyyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 

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

INTEGRAL:2 : 4 : 0 : SSSG_0 Y SSSG_1 Y SPSF_1 Y SPSG_0 Y SPSG_1 Y SDSG_0 Y 

0_0_0_zzzz_0 0_0_0_zzzz_1 0_0_0_yzzz_0 0_0_0_yzzz_1 0_0_0_yyzz_0 0_0_0_yyzz_1 0_0_0_yyyz_0 0_0_0_yyyz_1 0_0_0_yyyy_0 0_0_0_yyyy_1 0_0_0_xzzz_0 0_0_0_xzzz_1 0_0_0_xyzz_0 0_0_0_xyzz_1 0_0_0_xyyz_0 0_0_0_xyyz_1 0_0_0_xyyy_0 0_0_0_xyyy_1 0_0_0_xxzz_0 0_0_0_xxzz_1 0_0_0_xxyz_0 0_0_0_xxyz_1 0_0_0_xxyy_0 0_0_0_xxyy_1 0_0_0_xxxz_0 0_0_0_xxxz_1 0_0_0_xxxy_0 0_0_0_xxxy_1 0_0_0_xxxx_0 0_0_0_xxxx_1 0_z_0_zzz_1 0_z_0_zzzz_0 0_z_0_zzzz_1 0_z_0_yzz_1 0_z_0_yzzz_0 0_z_0_yzzz_1 0_z_0_yyz_1 0_z_0_yyzz_0 0_z_0_yyzz_1 0_z_0_yyy_1 0_z_0_yyyz_0 0_z_0_yyyz_1 0_z_0_yyyy_0 0_z_0_yyyy_1 0_z_0_xzz_1 0_z_0_xzzz_0 0_z_0_xzzz_1 0_z_0_xyz_1 0_z_0_xyzz_0 0_z_0_xyzz_1 0_z_0_xyy_1 0_z_0_xyyz_0 0_z_0_xyyz_1 0_z_0_xyyy_0 0_z_0_xyyy_1 0_z_0_xxz_1 0_z_0_xxzz_0 0_z_0_xxzz_1 0_z_0_xxy_1 0_z_0_xxyz_0 0_z_0_xxyz_1 0_z_0_xxyy_0 0_z_0_xxyy_1 0_z_0_xxx_1 0_z_0_xxxz_0 0_z_0_xxxz_1 0_z_0_xxxy_0 0_z_0_xxxy_1 0_z_0_xxxx_0 0_z_0_xxxx_1 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_y_0_zzz_1 0_y_0_zzzz_0 0_y_0_zzzz_1 0_y_0_yzz_1 0_y_0_yzzz_0 0_y_0_yzzz_1 0_y_0_yyz_1 0_y_0_yyzz_0 0_y_0_yyzz_1 0_y_0_yyy_1 0_y_0_yyyz_0 0_y_0_yyyz_1 0_y_0_yyyy_0 0_y_0_yyyy_1 0_y_0_xzz_1 0_y_0_xzzz_0 0_y_0_xzzz_1 0_y_0_xyz_1 0_y_0_xyzz_0 0_y_0_xyzz_1 0_y_0_xyy_1 0_y_0_xyyz_0 0_y_0_xyyz_1 0_y_0_xyyy_0 0_y_0_xyyy_1 0_y_0_xxz_1 0_y_0_xxzz_0 0_y_0_xxzz_1 0_y_0_xxy_1 0_y_0_xxyz_0 0_y_0_xxyz_1 0_y_0_xxyy_0 0_y_0_xxyy_1 0_y_0_xxx_1 0_y_0_xxxz_0 0_y_0_xxxz_1 0_y_0_xxxy_0 0_y_0_xxxy_1 0_y_0_xxxx_0 0_y_0_xxxx_1 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzzz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_x_0_zzz_1 0_x_0_zzzz_0 0_x_0_zzzz_1 0_x_0_yzz_1 0_x_0_yzzz_0 0_x_0_yzzz_1 0_x_0_yyz_1 0_x_0_yyzz_0 0_x_0_yyzz_1 0_x_0_yyy_1 0_x_0_yyyz_0 0_x_0_yyyz_1 0_x_0_yyyy_0 0_x_0_yyyy_1 0_x_0_xzz_1 0_x_0_xzzz_0 0_x_0_xzzz_1 0_x_0_xyz_1 0_x_0_xyzz_0 0_x_0_xyzz_1 0_x_0_xyy_1 0_x_0_xyyz_0 0_x_0_xyyz_1 0_x_0_xyyy_0 0_x_0_xyyy_1 0_x_0_xxz_1 0_x_0_xxzz_0 0_x_0_xxzz_1 0_x_0_xxy_1 0_x_0_xxyz_0 0_x_0_xxyz_1 0_x_0_xxyy_0 0_x_0_xxyy_1 0_x_0_xxx_1 0_x_0_xxxz_0 0_x_0_xxxz_1 0_x_0_xxxy_0 0_x_0_xxxy_1 0_x_0_xxxx_0 0_x_0_xxxx_1 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxzz_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

SSSG_0 SSSG_1 SPSF_1 SPSG_0 SPSG_1 SDSG_0 

SSSG_0 : 0_0_0_zzzz_0 0_0_0_yzzz_0 0_0_0_yyzz_0 0_0_0_yyyz_0 0_0_0_yyyy_0 0_0_0_xzzz_0 0_0_0_xyzz_0 0_0_0_xyyz_0 0_0_0_xyyy_0 0_0_0_xxzz_0 0_0_0_xxyz_0 0_0_0_xxyy_0 0_0_0_xxxz_0 0_0_0_xxxy_0 0_0_0_xxxx_0 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SPSF_1 : 0_z_0_zzz_1 0_z_0_yzz_1 0_z_0_yyz_1 0_z_0_yyy_1 0_z_0_xzz_1 0_z_0_xyz_1 0_z_0_xyy_1 0_z_0_xxz_1 0_z_0_xxy_1 0_z_0_xxx_1 0_y_0_zzz_1 0_y_0_yzz_1 0_y_0_yyz_1 0_y_0_yyy_1 0_y_0_xzz_1 0_y_0_xyz_1 0_y_0_xyy_1 0_y_0_xxz_1 0_y_0_xxy_1 0_y_0_xxx_1 0_x_0_zzz_1 0_x_0_yzz_1 0_x_0_yyz_1 0_x_0_yyy_1 0_x_0_xzz_1 0_x_0_xyz_1 0_x_0_xyy_1 0_x_0_xxz_1 0_x_0_xxy_1 0_x_0_xxx_1 

SPSG_0 : 0_z_0_zzzz_0 0_z_0_yzzz_0 0_z_0_yyzz_0 0_z_0_yyyz_0 0_z_0_yyyy_0 0_z_0_xzzz_0 0_z_0_xyzz_0 0_z_0_xyyz_0 0_z_0_xyyy_0 0_z_0_xxzz_0 0_z_0_xxyz_0 0_z_0_xxyy_0 0_z_0_xxxz_0 0_z_0_xxxy_0 0_z_0_xxxx_0 0_y_0_zzzz_0 0_y_0_yzzz_0 0_y_0_yyzz_0 0_y_0_yyyz_0 0_y_0_yyyy_0 0_y_0_xzzz_0 0_y_0_xyzz_0 0_y_0_xyyz_0 0_y_0_xyyy_0 0_y_0_xxzz_0 0_y_0_xxyz_0 0_y_0_xxyy_0 0_y_0_xxxz_0 0_y_0_xxxy_0 0_y_0_xxxx_0 0_x_0_zzzz_0 0_x_0_yzzz_0 0_x_0_yyzz_0 0_x_0_yyyz_0 0_x_0_yyyy_0 0_x_0_xzzz_0 0_x_0_xyzz_0 0_x_0_xyyz_0 0_x_0_xyyy_0 0_x_0_xxzz_0 0_x_0_xxyz_0 0_x_0_xxyy_0 0_x_0_xxxz_0 0_x_0_xxxy_0 0_x_0_xxxx_0 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_yyyy_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xyyy_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxyy_1 0_z_0_xxxz_1 0_z_0_xxxy_1 0_z_0_xxxx_1 0_y_0_zzzz_1 0_y_0_yzzz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xzzz_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxzz_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxz_1 0_y_0_xxxy_1 0_y_0_xxxx_1 0_x_0_zzzz_1 0_x_0_yzzz_1 0_x_0_yyzz_1 0_x_0_yyyz_1 0_x_0_yyyy_1 0_x_0_xzzz_1 0_x_0_xyzz_1 0_x_0_xyyz_1 0_x_0_xyyy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SDSG_0 : 0_zz_0_zzzz_0 0_zz_0_yzzz_0 0_zz_0_yyzz_0 0_zz_0_yyyz_0 0_zz_0_yyyy_0 0_zz_0_xzzz_0 0_zz_0_xyzz_0 0_zz_0_xyyz_0 0_zz_0_xyyy_0 0_zz_0_xxzz_0 0_zz_0_xxyz_0 0_zz_0_xxyy_0 0_zz_0_xxxz_0 0_zz_0_xxxy_0 0_zz_0_xxxx_0 0_yz_0_zzzz_0 0_yz_0_yzzz_0 0_yz_0_yyzz_0 0_yz_0_yyyz_0 0_yz_0_yyyy_0 0_yz_0_xzzz_0 0_yz_0_xyzz_0 0_yz_0_xyyz_0 0_yz_0_xyyy_0 0_yz_0_xxzz_0 0_yz_0_xxyz_0 0_yz_0_xxyy_0 0_yz_0_xxxz_0 0_yz_0_xxxy_0 0_yz_0_xxxx_0 0_yy_0_zzzz_0 0_yy_0_yzzz_0 0_yy_0_yyzz_0 0_yy_0_yyyz_0 0_yy_0_yyyy_0 0_yy_0_xzzz_0 0_yy_0_xyzz_0 0_yy_0_xyyz_0 0_yy_0_xyyy_0 0_yy_0_xxzz_0 0_yy_0_xxyz_0 0_yy_0_xxyy_0 0_yy_0_xxxz_0 0_yy_0_xxxy_0 0_yy_0_xxxx_0 0_xz_0_zzzz_0 0_xz_0_yzzz_0 0_xz_0_yyzz_0 0_xz_0_yyyz_0 0_xz_0_yyyy_0 0_xz_0_xzzz_0 0_xz_0_xyzz_0 0_xz_0_xyyz_0 0_xz_0_xyyy_0 0_xz_0_xxzz_0 0_xz_0_xxyz_0 0_xz_0_xxyy_0 0_xz_0_xxxz_0 0_xz_0_xxxy_0 0_xz_0_xxxx_0 0_xy_0_zzzz_0 0_xy_0_yzzz_0 0_xy_0_yyzz_0 0_xy_0_yyyz_0 0_xy_0_yyyy_0 0_xy_0_xzzz_0 0_xy_0_xyzz_0 0_xy_0_xyyz_0 0_xy_0_xyyy_0 0_xy_0_xxzz_0 0_xy_0_xxyz_0 0_xy_0_xxyy_0 0_xy_0_xxxz_0 0_xy_0_xxxy_0 0_xy_0_xxxx_0 0_xx_0_zzzz_0 0_xx_0_yzzz_0 0_xx_0_yyzz_0 0_xx_0_yyyz_0 0_xx_0_yyyy_0 0_xx_0_xzzz_0 0_xx_0_xyzz_0 0_xx_0_xyyz_0 0_xx_0_xyyy_0 0_xx_0_xxzz_0 0_xx_0_xxyz_0 0_xx_0_xxyy_0 0_xx_0_xxxz_0 0_xx_0_xxxy_0 0_xx_0_xxxx_0 

INTEGRAL:2 : 4 : 1 : SSSG_1 Y SSSG_2 Y SPSF_2 Y SPSG_1 Y SPSG_2 Y SDSG_1 Y 

0_0_0_zzzz_1 0_0_0_zzzz_2 0_0_0_yzzz_1 0_0_0_yzzz_2 0_0_0_yyzz_1 0_0_0_yyzz_2 0_0_0_yyyz_1 0_0_0_yyyz_2 0_0_0_yyyy_1 0_0_0_yyyy_2 0_0_0_xzzz_1 0_0_0_xzzz_2 0_0_0_xyzz_1 0_0_0_xyzz_2 0_0_0_xyyz_1 0_0_0_xyyz_2 0_0_0_xyyy_1 0_0_0_xyyy_2 0_0_0_xxzz_1 0_0_0_xxzz_2 0_0_0_xxyz_1 0_0_0_xxyz_2 0_0_0_xxyy_1 0_0_0_xxyy_2 0_0_0_xxxz_1 0_0_0_xxxz_2 0_0_0_xxxy_1 0_0_0_xxxy_2 0_0_0_xxxx_1 0_0_0_xxxx_2 0_z_0_zzz_2 0_z_0_zzzz_1 0_z_0_zzzz_2 0_z_0_yzz_2 0_z_0_yzzz_1 0_z_0_yzzz_2 0_z_0_yyz_2 0_z_0_yyzz_1 0_z_0_yyzz_2 0_z_0_yyy_2 0_z_0_yyyz_1 0_z_0_yyyz_2 0_z_0_yyyy_1 0_z_0_yyyy_2 0_z_0_xzz_2 0_z_0_xzzz_1 0_z_0_xzzz_2 0_z_0_xyz_2 0_z_0_xyzz_1 0_z_0_xyzz_2 0_z_0_xyy_2 0_z_0_xyyz_1 0_z_0_xyyz_2 0_z_0_xyyy_1 0_z_0_xyyy_2 0_z_0_xxz_2 0_z_0_xxzz_1 0_z_0_xxzz_2 0_z_0_xxy_2 0_z_0_xxyz_1 0_z_0_xxyz_2 0_z_0_xxyy_1 0_z_0_xxyy_2 0_z_0_xxx_2 0_z_0_xxxz_1 0_z_0_xxxz_2 0_z_0_xxxy_1 0_z_0_xxxy_2 0_z_0_xxxx_1 0_z_0_xxxx_2 0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_yyyz_1 0_zz_0_yyyy_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xyyz_1 0_zz_0_xyyy_1 0_zz_0_xxzz_1 0_zz_0_xxyz_1 0_zz_0_xxyy_1 0_zz_0_xxxz_1 0_zz_0_xxxy_1 0_zz_0_xxxx_1 0_y_0_zzz_2 0_y_0_zzzz_1 0_y_0_zzzz_2 0_y_0_yzz_2 0_y_0_yzzz_1 0_y_0_yzzz_2 0_y_0_yyz_2 0_y_0_yyzz_1 0_y_0_yyzz_2 0_y_0_yyy_2 0_y_0_yyyz_1 0_y_0_yyyz_2 0_y_0_yyyy_1 0_y_0_yyyy_2 0_y_0_xzz_2 0_y_0_xzzz_1 0_y_0_xzzz_2 0_y_0_xyz_2 0_y_0_xyzz_1 0_y_0_xyzz_2 0_y_0_xyy_2 0_y_0_xyyz_1 0_y_0_xyyz_2 0_y_0_xyyy_1 0_y_0_xyyy_2 0_y_0_xxz_2 0_y_0_xxzz_1 0_y_0_xxzz_2 0_y_0_xxy_2 0_y_0_xxyz_1 0_y_0_xxyz_2 0_y_0_xxyy_1 0_y_0_xxyy_2 0_y_0_xxx_2 0_y_0_xxxz_1 0_y_0_xxxz_2 0_y_0_xxxy_1 0_y_0_xxxy_2 0_y_0_xxxx_1 0_y_0_xxxx_2 0_yz_0_zzzz_1 0_yz_0_yzzz_1 0_yz_0_yyzz_1 0_yz_0_yyyz_1 0_yz_0_yyyy_1 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_zzzz_1 0_yy_0_yzzz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xzzz_1 0_yy_0_xyzz_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxzz_1 0_yy_0_xxyz_1 0_yy_0_xxyy_1 0_yy_0_xxxz_1 0_yy_0_xxxy_1 0_yy_0_xxxx_1 0_x_0_zzz_2 0_x_0_zzzz_1 0_x_0_zzzz_2 0_x_0_yzz_2 0_x_0_yzzz_1 0_x_0_yzzz_2 0_x_0_yyz_2 0_x_0_yyzz_1 0_x_0_yyzz_2 0_x_0_yyy_2 0_x_0_yyyz_1 0_x_0_yyyz_2 0_x_0_yyyy_1 0_x_0_yyyy_2 0_x_0_xzz_2 0_x_0_xzzz_1 0_x_0_xzzz_2 0_x_0_xyz_2 0_x_0_xyzz_1 0_x_0_xyzz_2 0_x_0_xyy_2 0_x_0_xyyz_1 0_x_0_xyyz_2 0_x_0_xyyy_1 0_x_0_xyyy_2 0_x_0_xxz_2 0_x_0_xxzz_1 0_x_0_xxzz_2 0_x_0_xxy_2 0_x_0_xxyz_1 0_x_0_xxyz_2 0_x_0_xxyy_1 0_x_0_xxyy_2 0_x_0_xxx_2 0_x_0_xxxz_1 0_x_0_xxxz_2 0_x_0_xxxy_1 0_x_0_xxxy_2 0_x_0_xxxx_1 0_x_0_xxxx_2 0_xz_0_xzzz_1 0_xz_0_xxzz_1 0_xz_0_xxxz_1 0_xz_0_xxxx_1 0_xy_0_xyyy_1 0_xy_0_xxyy_1 0_xy_0_xxxy_1 0_xx_0_zzzz_1 0_xx_0_yzzz_1 0_xx_0_yyzz_1 0_xx_0_yyyz_1 0_xx_0_yyyy_1 0_xx_0_xzzz_1 0_xx_0_xyzz_1 0_xx_0_xyyz_1 0_xx_0_xyyy_1 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 

SSSG_1 SSSG_2 SPSF_2 SPSG_1 SPSG_2 SDSG_1 

SSSG_1 : 0_0_0_zzzz_1 0_0_0_yzzz_1 0_0_0_yyzz_1 0_0_0_yyyz_1 0_0_0_yyyy_1 0_0_0_xzzz_1 0_0_0_xyzz_1 0_0_0_xyyz_1 0_0_0_xyyy_1 0_0_0_xxzz_1 0_0_0_xxyz_1 0_0_0_xxyy_1 0_0_0_xxxz_1 0_0_0_xxxy_1 0_0_0_xxxx_1 

SSSG_2 : 0_0_0_zzzz_2 0_0_0_yzzz_2 0_0_0_yyzz_2 0_0_0_yyyz_2 0_0_0_yyyy_2 0_0_0_xzzz_2 0_0_0_xyzz_2 0_0_0_xyyz_2 0_0_0_xyyy_2 0_0_0_xxzz_2 0_0_0_xxyz_2 0_0_0_xxyy_2 0_0_0_xxxz_2 0_0_0_xxxy_2 0_0_0_xxxx_2 

SPSF_2 : 0_z_0_zzz_2 0_z_0_yzz_2 0_z_0_yyz_2 0_z_0_yyy_2 0_z_0_xzz_2 0_z_0_xyz_2 0_z_0_xyy_2 0_z_0_xxz_2 0_z_0_xxy_2 0_z_0_xxx_2 0_y_0_zzz_2 0_y_0_yzz_2 0_y_0_yyz_2 0_y_0_yyy_2 0_y_0_xzz_2 0_y_0_xyz_2 0_y_0_xyy_2 0_y_0_xxz_2 0_y_0_xxy_2 0_y_0_xxx_2 0_x_0_zzz_2 0_x_0_yzz_2 0_x_0_yyz_2 0_x_0_yyy_2 0_x_0_xzz_2 0_x_0_xyz_2 0_x_0_xyy_2 0_x_0_xxz_2 0_x_0_xxy_2 0_x_0_xxx_2 

SPSG_1 : 0_z_0_zzzz_1 0_z_0_yzzz_1 0_z_0_yyzz_1 0_z_0_yyyz_1 0_z_0_yyyy_1 0_z_0_xzzz_1 0_z_0_xyzz_1 0_z_0_xyyz_1 0_z_0_xyyy_1 0_z_0_xxzz_1 0_z_0_xxyz_1 0_z_0_xxyy_1 0_z_0_xxxz_1 0_z_0_xxxy_1 0_z_0_xxxx_1 0_y_0_zzzz_1 0_y_0_yzzz_1 0_y_0_yyzz_1 0_y_0_yyyz_1 0_y_0_yyyy_1 0_y_0_xzzz_1 0_y_0_xyzz_1 0_y_0_xyyz_1 0_y_0_xyyy_1 0_y_0_xxzz_1 0_y_0_xxyz_1 0_y_0_xxyy_1 0_y_0_xxxz_1 0_y_0_xxxy_1 0_y_0_xxxx_1 0_x_0_zzzz_1 0_x_0_yzzz_1 0_x_0_yyzz_1 0_x_0_yyyz_1 0_x_0_yyyy_1 0_x_0_xzzz_1 0_x_0_xyzz_1 0_x_0_xyyz_1 0_x_0_xyyy_1 0_x_0_xxzz_1 0_x_0_xxyz_1 0_x_0_xxyy_1 0_x_0_xxxz_1 0_x_0_xxxy_1 0_x_0_xxxx_1 

SPSG_2 : 0_z_0_zzzz_2 0_z_0_yzzz_2 0_z_0_yyzz_2 0_z_0_yyyz_2 0_z_0_yyyy_2 0_z_0_xzzz_2 0_z_0_xyzz_2 0_z_0_xyyz_2 0_z_0_xyyy_2 0_z_0_xxzz_2 0_z_0_xxyz_2 0_z_0_xxyy_2 0_z_0_xxxz_2 0_z_0_xxxy_2 0_z_0_xxxx_2 0_y_0_zzzz_2 0_y_0_yzzz_2 0_y_0_yyzz_2 0_y_0_yyyz_2 0_y_0_yyyy_2 0_y_0_xzzz_2 0_y_0_xyzz_2 0_y_0_xyyz_2 0_y_0_xyyy_2 0_y_0_xxzz_2 0_y_0_xxyz_2 0_y_0_xxyy_2 0_y_0_xxxz_2 0_y_0_xxxy_2 0_y_0_xxxx_2 0_x_0_zzzz_2 0_x_0_yzzz_2 0_x_0_yyzz_2 0_x_0_yyyz_2 0_x_0_yyyy_2 0_x_0_xzzz_2 0_x_0_xyzz_2 0_x_0_xyyz_2 0_x_0_xyyy_2 0_x_0_xxzz_2 0_x_0_xxyz_2 0_x_0_xxyy_2 0_x_0_xxxz_2 0_x_0_xxxy_2 0_x_0_xxxx_2 

SDSG_1 : 0_zz_0_zzzz_1 0_zz_0_yzzz_1 0_zz_0_yyzz_1 0_zz_0_yyyz_1 0_zz_0_yyyy_1 0_zz_0_xzzz_1 0_zz_0_xyzz_1 0_zz_0_xyyz_1 0_zz_0_xyyy_1 0_zz_0_xxzz_1 0_zz_0_xxyz_1 0_zz_0_xxyy_1 0_zz_0_xxxz_1 0_zz_0_xxxy_1 0_zz_0_xxxx_1 0_yz_0_zzzz_1 0_yz_0_yzzz_1 0_yz_0_yyzz_1 0_yz_0_yyyz_1 0_yz_0_yyyy_1 0_yz_0_xyzz_1 0_yz_0_xyyz_1 0_yz_0_xxyz_1 0_yy_0_zzzz_1 0_yy_0_yzzz_1 0_yy_0_yyzz_1 0_yy_0_yyyz_1 0_yy_0_yyyy_1 0_yy_0_xzzz_1 0_yy_0_xyzz_1 0_yy_0_xyyz_1 0_yy_0_xyyy_1 0_yy_0_xxzz_1 0_yy_0_xxyz_1 0_yy_0_xxyy_1 0_yy_0_xxxz_1 0_yy_0_xxxy_1 0_yy_0_xxxx_1 0_xz_0_xzzz_1 0_xz_0_xxzz_1 0_xz_0_xxxz_1 0_xz_0_xxxx_1 0_xy_0_xyyy_1 0_xy_0_xxyy_1 0_xy_0_xxxy_1 0_xx_0_zzzz_1 0_xx_0_yzzz_1 0_xx_0_yyzz_1 0_xx_0_yyyz_1 0_xx_0_yyyy_1 0_xx_0_xzzz_1 0_xx_0_xyzz_1 0_xx_0_xyyz_1 0_xx_0_xyyy_1 0_xx_0_xxzz_1 0_xx_0_xxyz_1 0_xx_0_xxyy_1 0_xx_0_xxxz_1 0_xx_0_xxxy_1 0_xx_0_xxxx_1 

INTEGRAL:3 : 2 : 0 : SPSD_0 Y SPSD_1 Y SDSP_1 Y SDSD_0 N SDSD_1 Y SFSD_0 Y 

0_z_0_zz_0 0_z_0_zz_1 0_z_0_yz_0 0_z_0_yz_1 0_z_0_yy_0 0_z_0_yy_1 0_z_0_xz_0 0_z_0_xz_1 0_z_0_xy_0 0_z_0_xy_1 0_z_0_xx_0 0_z_0_xx_1 0_zz_0_z_1 0_zz_0_zz_0 0_zz_0_zz_1 0_zz_0_y_1 0_zz_0_yz_0 0_zz_0_yz_1 0_zz_0_yy_0 0_zz_0_yy_1 0_zz_0_x_1 0_zz_0_xz_0 0_zz_0_xz_1 0_zz_0_xy_0 0_zz_0_xy_1 0_zz_0_xx_0 0_zz_0_xx_1 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_y_0_zz_0 0_y_0_zz_1 0_y_0_yz_0 0_y_0_yz_1 0_y_0_yy_0 0_y_0_yy_1 0_y_0_xz_0 0_y_0_xz_1 0_y_0_xy_0 0_y_0_xy_1 0_y_0_xx_0 0_y_0_xx_1 0_yz_0_zz_0 0_yz_0_zz_1 0_yz_0_yz_0 0_yz_0_yz_1 0_yz_0_yy_0 0_yz_0_yy_1 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yy_0_z_1 0_yy_0_zz_0 0_yy_0_zz_1 0_yy_0_y_1 0_yy_0_yz_0 0_yy_0_yz_1 0_yy_0_yy_0 0_yy_0_yy_1 0_yy_0_x_1 0_yy_0_xz_0 0_yy_0_xz_1 0_yy_0_xy_0 0_yy_0_xy_1 0_yy_0_xx_0 0_yy_0_xx_1 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_x_0_zz_0 0_x_0_zz_1 0_x_0_yz_0 0_x_0_yz_1 0_x_0_yy_0 0_x_0_yy_1 0_x_0_xz_0 0_x_0_xz_1 0_x_0_xy_0 0_x_0_xy_1 0_x_0_xx_0 0_x_0_xx_1 0_xz_0_xz_0 0_xz_0_xz_1 0_xz_0_xx_0 0_xz_0_xx_1 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xy_0_xy_0 0_xy_0_xy_1 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xx_0_z_1 0_xx_0_zz_0 0_xx_0_zz_1 0_xx_0_y_1 0_xx_0_yz_0 0_xx_0_yz_1 0_xx_0_yy_0 0_xx_0_yy_1 0_xx_0_x_1 0_xx_0_xz_0 0_xx_0_xz_1 0_xx_0_xy_0 0_xx_0_xy_1 0_xx_0_xx_0 0_xx_0_xx_1 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

SPSD_0 SPSD_1 SDSP_1 SDSD_0 SDSD_1 SFSD_0 

SPSD_0 : 0_z_0_zz_0 0_z_0_yz_0 0_z_0_yy_0 0_z_0_xz_0 0_z_0_xy_0 0_z_0_xx_0 0_y_0_zz_0 0_y_0_yz_0 0_y_0_yy_0 0_y_0_xz_0 0_y_0_xy_0 0_y_0_xx_0 0_x_0_zz_0 0_x_0_yz_0 0_x_0_yy_0 0_x_0_xz_0 0_x_0_xy_0 0_x_0_xx_0 

SPSD_1 : 0_z_0_zz_1 0_z_0_yz_1 0_z_0_yy_1 0_z_0_xz_1 0_z_0_xy_1 0_z_0_xx_1 0_y_0_zz_1 0_y_0_yz_1 0_y_0_yy_1 0_y_0_xz_1 0_y_0_xy_1 0_y_0_xx_1 0_x_0_zz_1 0_x_0_yz_1 0_x_0_yy_1 0_x_0_xz_1 0_x_0_xy_1 0_x_0_xx_1 

SDSP_1 : 0_zz_0_z_1 0_zz_0_y_1 0_zz_0_x_1 0_yy_0_z_1 0_yy_0_y_1 0_yy_0_x_1 0_xx_0_z_1 0_xx_0_y_1 0_xx_0_x_1 

SDSD_0 : 0_zz_0_zz_0 0_zz_0_yz_0 0_zz_0_yy_0 0_zz_0_xz_0 0_zz_0_xy_0 0_zz_0_xx_0 0_yz_0_zz_0 0_yz_0_yz_0 0_yz_0_yy_0 0_yz_0_xz_0 0_yz_0_xy_0 0_yz_0_xx_0 0_yy_0_zz_0 0_yy_0_yz_0 0_yy_0_yy_0 0_yy_0_xz_0 0_yy_0_xy_0 0_yy_0_xx_0 0_xz_0_zz_0 0_xz_0_yz_0 0_xz_0_yy_0 0_xz_0_xz_0 0_xz_0_xy_0 0_xz_0_xx_0 0_xy_0_zz_0 0_xy_0_yz_0 0_xy_0_yy_0 0_xy_0_xz_0 0_xy_0_xy_0 0_xy_0_xx_0 0_xx_0_zz_0 0_xx_0_yz_0 0_xx_0_yy_0 0_xx_0_xz_0 0_xx_0_xy_0 0_xx_0_xx_0 

SDSD_1 : 0_zz_0_zz_1 0_zz_0_yz_1 0_zz_0_yy_1 0_zz_0_xz_1 0_zz_0_xy_1 0_zz_0_xx_1 0_yz_0_zz_1 0_yz_0_yz_1 0_yz_0_yy_1 0_yy_0_zz_1 0_yy_0_yz_1 0_yy_0_yy_1 0_yy_0_xz_1 0_yy_0_xy_1 0_yy_0_xx_1 0_xz_0_xz_1 0_xz_0_xx_1 0_xy_0_xy_1 0_xx_0_zz_1 0_xx_0_yz_1 0_xx_0_yy_1 0_xx_0_xz_1 0_xx_0_xy_1 0_xx_0_xx_1 

SFSD_0 : 0_zzz_0_zz_0 0_zzz_0_yz_0 0_zzz_0_yy_0 0_zzz_0_xz_0 0_zzz_0_xy_0 0_zzz_0_xx_0 0_yzz_0_zz_0 0_yzz_0_yz_0 0_yzz_0_yy_0 0_yzz_0_xz_0 0_yzz_0_xy_0 0_yzz_0_xx_0 0_yyz_0_zz_0 0_yyz_0_yz_0 0_yyz_0_yy_0 0_yyz_0_xz_0 0_yyz_0_xy_0 0_yyz_0_xx_0 0_yyy_0_zz_0 0_yyy_0_yz_0 0_yyy_0_yy_0 0_yyy_0_xz_0 0_yyy_0_xy_0 0_yyy_0_xx_0 0_xzz_0_zz_0 0_xzz_0_yz_0 0_xzz_0_yy_0 0_xzz_0_xz_0 0_xzz_0_xy_0 0_xzz_0_xx_0 0_xyz_0_zz_0 0_xyz_0_yz_0 0_xyz_0_yy_0 0_xyz_0_xz_0 0_xyz_0_xy_0 0_xyz_0_xx_0 0_xyy_0_zz_0 0_xyy_0_yz_0 0_xyy_0_yy_0 0_xyy_0_xz_0 0_xyy_0_xy_0 0_xyy_0_xx_0 0_xxz_0_zz_0 0_xxz_0_yz_0 0_xxz_0_yy_0 0_xxz_0_xz_0 0_xxz_0_xy_0 0_xxz_0_xx_0 0_xxy_0_zz_0 0_xxy_0_yz_0 0_xxy_0_yy_0 0_xxy_0_xz_0 0_xxy_0_xy_0 0_xxy_0_xx_0 0_xxx_0_zz_0 0_xxx_0_yz_0 0_xxx_0_yy_0 0_xxx_0_xz_0 0_xxx_0_xy_0 0_xxx_0_xx_0 

